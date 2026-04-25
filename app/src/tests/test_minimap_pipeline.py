"""
test_minimap_pipeline.py

Full Case 2 pipeline: minimap2-based coordinate lifting from reference to query.

Pipeline:
    1. Load ref + query GenBank
    2. Auto-detect virus → load alias config
    3. Parse ref features (mat_peptide or CDS) → normalize names via alias
    4. Write ref + query as FASTA → run minimap2 → build ref_to_query map
    5. Lift ref features onto query → validate start/stop codons
    6. (Optional) Compare against actual query annotation if present
    7. Output: FASTA per gene + TSV summary

Usage:
    python -m app.src.tests.test_minimap_pipeline \
        --reference app/data/FMD_ref_test.gb \
        --query app/data/FMD_FJ175661_Anno.gb \
        --alias-registry app/config/virus_alias_registry.json \
        --output-dir output/minimap_test
"""

import argparse
import csv
import tempfile
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from app.src.alignment.minimap_runner import run_minimap2
from app.src.alignment.sam_lifter import build_ref_to_query_map, get_primary_alignment
from app.src.annotation.alias_registry import detect_alias_config_for_record, get_detected_virus_name
from app.src.annotation.annotation_strategy import choose_strategy
from app.src.annotation.extractor import extract_all_lifted
from app.src.annotation.gene_alias import (
    apply_alias_to_features,
    build_canonical_to_ref_map,
    apply_ref_naming,
    load_alias_lookup,
)
from app.src.io.genbank_parser import (
    load_single_genbank,
    parse_cds_features,
    parse_mat_peptides,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _write_fasta(record: SeqRecord, path: Path) -> None:
    SeqIO.write(record, str(path), "fasta")


def _parse_ref_features(ref_record: SeqRecord):
    """Parse ref features using annotation strategy."""
    _, feature_type = choose_strategy(ref_record)
    if feature_type == "mat_peptide":
        return parse_mat_peptides(ref_record), "mat_peptide"
    elif feature_type == "CDS":
        return parse_cds_features(ref_record), "CDS"
    else:
        return [], None


def _normalize_ref_features(ref_features, alias_lookup, use_ref_naming: bool):
    """Apply alias normalization + optional ref naming to ref features."""
    if not alias_lookup:
        return ref_features

    normalized = apply_alias_to_features(ref_features, alias_lookup)

    if use_ref_naming:
        canonical_to_ref = build_canonical_to_ref_map(ref_features, alias_lookup)
        normalized = apply_ref_naming(normalized, canonical_to_ref)

    return normalized


def _parse_actual_annotation(query_record: SeqRecord, alias_lookup: dict):
    """
    Parse actual annotation from query for comparison.
    Tries mat_peptide first, then CDS.
    Returns list of dicts with name (canonical), start, end.
    """
    _, feature_type = choose_strategy(query_record)
    if feature_type == "mat_peptide":
        features = parse_mat_peptides(query_record)
    elif feature_type == "CDS":
        features = parse_cds_features(query_record)
    else:
        return []

    if alias_lookup:
        # Keep features in canonical form for matching against lifted results
        features = apply_alias_to_features(features, alias_lookup)

    return features


def _run_alignment(ref_record: SeqRecord, query_record: SeqRecord, tmp_dir: Path):
    """Write FASTA files, run minimap2, return ref_to_query map."""
    ref_fasta = tmp_dir / "ref.fa"
    query_fasta = tmp_dir / "query.fa"
    sam_path = tmp_dir / "alignment.sam"

    _write_fasta(ref_record, ref_fasta)
    _write_fasta(query_record, query_fasta)

    run_minimap2(ref_fasta, query_fasta, sam_path, quiet=True)

    aln = get_primary_alignment(sam_path)
    return build_ref_to_query_map(aln)


# ---------------------------------------------------------------------------
# Output writers
# ---------------------------------------------------------------------------

def write_fasta_output(lifted_results, query_record: SeqRecord, output_dir: Path) -> Path:
    """Write one FASTA file with all successfully lifted sequences."""
    fasta_path = output_dir / f"{query_record.id}_lifted.fa"
    records = []
    for r in lifted_results:
        if r["status"] not in ("ok", "invalid_boundaries") or not r.get("sequence"):
            continue
        rec = SeqRecord(
            Seq(r["sequence"]),
            id=f"{query_record.id}|{r['name']}",
            description=f"status={r['status']} coverage={r['coverage']}",
        )
        records.append(rec)
    SeqIO.write(records, str(fasta_path), "fasta")
    return fasta_path


def write_tsv_output(lifted_results, output_dir: Path, query_id: str) -> Path:
    """Write TSV summary of lifted features."""
    tsv_path = output_dir / f"{query_id}_lifted.tsv"
    fieldnames = [
        "name", "status", "coverage",
        "ref_start", "ref_end",
        "query_start", "query_end",
        "has_start_codon", "has_stop_codon",
        "rescue_offset", "seq_len",
    ]
    with open(tsv_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for r in lifted_results:
            writer.writerow({
                "name": r["name"],
                "status": r["status"],
                "coverage": r["coverage"],
                "ref_start": r.get("ref_start"),
                "ref_end": r.get("ref_end"),
                "query_start": r.get("start"),
                "query_end": r.get("end"),
                "has_start_codon": r.get("has_start_codon", ""),
                "has_stop_codon": r.get("has_stop_codon", ""),
                "rescue_offset": r.get("rescue_offset", ""),
                "seq_len": len(r["sequence"]) if r.get("sequence") else "",
            })
    return tsv_path


def write_comparison_tsv(comparison_rows, output_dir: Path, query_id: str) -> Path:
    """Write side-by-side comparison of lifted vs actual annotation."""
    tsv_path = output_dir / f"{query_id}_comparison.tsv"
    fieldnames = [
        "name",
        "lifted_start", "lifted_end", "lifted_len",
        "actual_start", "actual_end", "actual_len",
        "start_diff", "end_diff", "len_diff",
        "lifted_status", "rescue_offset",
    ]
    with open(tsv_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(comparison_rows)
    return tsv_path


# ---------------------------------------------------------------------------
# Comparison logic
# ---------------------------------------------------------------------------

def build_comparison(lifted_results, actual_features) -> list:
    """
    Match lifted results to actual annotation by canonical name, produce comparison rows.

    Lifted features may have ref raw names in "name" and canonical name in "canonical_name".
    Actual features are normalized to canonical names.
    Matching uses canonical_name if available, else name.
    """
    actual_by_canonical = {f["name"]: f for f in actual_features}
    rows = []

    for r in lifted_results:
        # Prefer canonical_name for lookup (set by apply_ref_naming); fallback to name
        match_key = r.get("canonical_name") or r["name"]
        actual = actual_by_canonical.get(match_key)
        name = r["name"]  # display ref name in output

        lifted_start = r.get("start")
        lifted_end = r.get("end")
        lifted_len = (lifted_end - lifted_start + 1) if (lifted_start and lifted_end) else None

        if actual:
            actual_start = actual["start"]
            actual_end = actual["end"]
            actual_len = actual_end - actual_start + 1
            start_diff = (lifted_start - actual_start) if lifted_start else None
            end_diff = (lifted_end - actual_end) if lifted_end else None
            len_diff = (lifted_len - actual_len) if lifted_len else None
        else:
            actual_start = actual_end = actual_len = None
            start_diff = end_diff = len_diff = None

        rows.append({
            "name": name,
            "lifted_start": lifted_start,
            "lifted_end": lifted_end,
            "lifted_len": lifted_len,
            "actual_start": actual_start,
            "actual_end": actual_end,
            "actual_len": actual_len,
            "start_diff": start_diff,
            "end_diff": end_diff,
            "len_diff": len_diff,
            "lifted_status": r["status"],
            "rescue_offset": r.get("rescue_offset"),
        })

    return rows


def print_comparison(comparison_rows, virus_name: str, query_id: str):
    """Print comparison table to stdout."""
    print(f"\n{'='*80}")
    print(f"COMPARISON: {virus_name} | query={query_id}")
    print(f"{'='*80}")
    print(f"{'Name':<20} {'Lifted':>12} {'Actual':>12} {'dStart':>7} {'dEnd':>7} {'dLen':>7}  {'Status':<20} {'Rescue'}")
    print(f"{'-'*90}")

    for r in comparison_rows:
        lifted = f"{r['lifted_start']}-{r['lifted_end']}" if r["lifted_start"] else "unmapped"
        actual = f"{r['actual_start']}-{r['actual_end']}" if r["actual_start"] else "n/a"
        ds = f"{r['start_diff']:+d}" if r["start_diff"] is not None else "n/a"
        de = f"{r['end_diff']:+d}" if r["end_diff"] is not None else "n/a"
        dl = f"{r['len_diff']:+d}" if r["len_diff"] is not None else "n/a"
        rescue = f"offset={r['rescue_offset']:+d}" if r.get("rescue_offset") is not None else ""
        print(f"{r['name']:<20} {lifted:>12} {actual:>12} {ds:>7} {de:>7} {dl:>7}  {r['lifted_status']:<20} {rescue}")


def print_lifted_summary(lifted_results, virus_name: str, query_id: str):
    """Print lifted results when no actual annotation to compare against."""
    print(f"\n{'='*70}")
    print(f"LIFTED FEATURES: {virus_name} | query={query_id}")
    print(f"{'='*70}")
    print(f"{'Name':<20} {'Start':>8} {'End':>8} {'Len':>7} {'Cov':>6}  Status")
    print(f"{'-'*70}")
    for r in lifted_results:
        start = r.get("start") or ""
        end = r.get("end") or ""
        length = (r["end"] - r["start"] + 1) if (r.get("start") and r.get("end")) else ""
        cov = f"{r['coverage']:.3f}"
        print(f"{r['name']:<20} {str(start):>8} {str(end):>8} {str(length):>7} {cov:>6}  {r['status']}")


# ---------------------------------------------------------------------------
# Pipeline
# ---------------------------------------------------------------------------

def run_minimap_pipeline(
    ref_record: SeqRecord,
    query_record: SeqRecord,
    alias_lookup: dict,
    virus_name: str,
    output_dir: Path,
    compare: bool = True,
):
    """Run full Case 2 minimap pipeline for one query record."""
    print(f"\n--- Processing: {query_record.id} ---")

    # 1. Parse + normalize ref features
    ref_features, ref_feature_type = _parse_ref_features(ref_record)
    if not ref_features:
        print(f"  [SKIP] Reference has no parseable features.")
        return

    print(f"  Ref features ({ref_feature_type}): {len(ref_features)}")

    normalized_ref = _normalize_ref_features(ref_features, alias_lookup, use_ref_naming=True)

    # 2. Align
    with tempfile.TemporaryDirectory() as tmp:
        tmp_dir = Path(tmp)
        try:
            ref_to_query = _run_alignment(ref_record, query_record, tmp_dir)
        except Exception as e:
            print(f"  [ERROR] Alignment failed: {e}")
            return

    print(f"  Mapped positions: {len(ref_to_query)}")

    # 3. Lift + validate (skip codon check for mat_peptide — cleavage products have no ATG/stop)
    validate_codons = (ref_feature_type == "CDS")
    lifted = extract_all_lifted(normalized_ref, query_record, ref_to_query,
                                validate_codons=validate_codons)

    ok_count = sum(1 for r in lifted if r["status"] == "ok")
    inv_count = sum(1 for r in lifted if r["status"] == "invalid_boundaries")
    low_count = sum(1 for r in lifted if r["status"] == "low_coverage")
    unm_count = sum(1 for r in lifted if r["status"] == "unmapped")

    print(f"  Lift results: ok={ok_count}, invalid_boundaries={inv_count}, "
          f"low_coverage={low_count}, unmapped={unm_count}")

    # 4. Compare vs actual annotation (if present + requested)
    if compare:
        actual_features = _parse_actual_annotation(query_record, alias_lookup)
        if actual_features:
            comparison = build_comparison(lifted, actual_features)
            print_comparison(comparison, virus_name, query_record.id)
            cmp_path = write_comparison_tsv(comparison, output_dir, query_record.id)
            print(f"\n  Comparison TSV: {cmp_path}")
        else:
            print_lifted_summary(lifted, virus_name, query_record.id)
    else:
        print_lifted_summary(lifted, virus_name, query_record.id)

    # 5. Write outputs
    fasta_path = write_fasta_output(lifted, query_record, output_dir)
    tsv_path = write_tsv_output(lifted, output_dir, query_record.id)

    print(f"  FASTA output: {fasta_path}")
    print(f"  TSV output  : {tsv_path}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description="Case 2 minimap2 pipeline: lift ref features onto query genome."
    )
    parser.add_argument("--reference", required=True, help="Reference GenBank file.")
    parser.add_argument("--query", required=True, help="Query GenBank file.")
    parser.add_argument("--alias-registry", default="app/config/virus_alias_registry.json")
    parser.add_argument("--output-dir", default="output/minimap_pipeline")
    parser.add_argument(
        "--no-compare",
        action="store_true",
        help="Skip comparison with actual annotation (for truly unannotated queries).",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    ref_path = Path(args.reference)
    query_path = Path(args.query)
    registry_path = Path(args.alias_registry)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    ref_record = load_single_genbank(ref_path)
    query_record = load_single_genbank(query_path)

    print(f"Reference : {ref_record.id} ({len(ref_record.seq)} bp)")
    print(f"Query     : {query_record.id} ({len(query_record.seq)} bp)")

    # Detect virus from ref
    virus_name = get_detected_virus_name(ref_record, registry_path) or "Unknown"
    alias_config_path = detect_alias_config_for_record(ref_record, registry_path)

    if alias_config_path:
        alias_lookup = load_alias_lookup(alias_config_path)
        print(f"Virus     : {virus_name} (alias config: {alias_config_path})")
    else:
        alias_lookup = {}
        print(f"Virus     : {virus_name} (no alias config found — raw names will be used)")

    run_minimap_pipeline(
        ref_record=ref_record,
        query_record=query_record,
        alias_lookup=alias_lookup,
        virus_name=virus_name,
        output_dir=output_dir,
        compare=not args.no_compare,
    )


if __name__ == "__main__":
    main()
