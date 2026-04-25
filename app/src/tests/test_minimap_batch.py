"""
test_minimap_batch.py

Batch accuracy test for Case 2 minimap pipeline.

For each annotated record in a multi-record GenBank file:
    - Run minimap2 lift from reference (with asm5 -> asm10 -> asm20 fallback)
    - Filter out sub-features from actual annotation before comparison
    - Compare lifted coordinates vs actual annotation (ground truth)
    - Skip match if lifted/actual length ratio is too different (size_mismatch)
    - Report per-record and aggregate accuracy

Only processes records that have actual annotation (mat_peptide or CDS > 1)
so comparison is meaningful.

Usage:
    python -m app.src.tests.test_minimap_batch \
        --reference app/data/FMD_ref_test.gb \
        --query app/data/FMD_test.gb \
        --alias-registry app/config/virus_alias_registry.json \
        --output-dir output/minimap_batch/fmd
"""

import argparse
import csv
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from app.src.alignment.minimap_runner import run_minimap2_with_fallback
from app.src.alignment.sam_lifter import build_ref_to_query_map, get_primary_alignment
from app.src.annotation.alias_registry import detect_alias_config_for_record, get_detected_virus_name
from app.src.annotation.annotation_strategy import choose_strategy
from app.src.annotation.extractor import extract_all_lifted
from app.src.annotation.gene_alias import (
    apply_alias_to_features,
    apply_ref_naming,
    build_canonical_to_ref_map,
    load_alias_lookup,
)
from app.src.io.genbank_parser import (
    load_genbank_records,
    load_single_genbank,
    parse_cds_features,
    parse_mat_peptides,
)

# A lifted feature is considered a size mismatch if its length ratio vs actual
# is outside this range. e.g. 0.5 means lifted must be between 50%-200% of actual.
MIN_LEN_RATIO = 0.5
MAX_LEN_RATIO = 2.0

# A feature is considered a sub-feature of another if it is fully contained
# within it AND its length is less than this fraction of the parent.
SUBFEATURE_MAX_RATIO = 0.8


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _write_fasta(record: SeqRecord, path: Path) -> None:
    SeqIO.write(record, str(path), "fasta")


def _has_full_annotation(record: SeqRecord) -> bool:
    """Record is considered annotated if it has mat_peptide OR more than 1 CDS."""
    if parse_mat_peptides(record):
        return True
    cds = parse_cds_features(record)
    return len(cds) > 1


def _parse_ref_features(ref_record: SeqRecord):
    _, feature_type = choose_strategy(ref_record)
    if feature_type == "mat_peptide":
        return parse_mat_peptides(ref_record), "mat_peptide"
    elif feature_type == "CDS":
        return parse_cds_features(ref_record), "CDS"
    return [], None


def _normalize_ref_features(ref_features, alias_lookup):
    if not alias_lookup:
        return ref_features
    normalized = apply_alias_to_features(ref_features, alias_lookup)
    canonical_to_ref = build_canonical_to_ref_map(ref_features, alias_lookup)
    return apply_ref_naming(normalized, canonical_to_ref)


def _filter_subfeatures(features: List[Dict]) -> Tuple[List[Dict], List[Dict]]:
    """
    Detect and separate sub-features from the feature list.

    A feature is a sub-feature if:
      - it is fully contained within another feature (start >= parent.start AND end <= parent.end)
      - its length < SUBFEATURE_MAX_RATIO * parent length

    Returns: (main_features, sub_features)
    """
    main, subs = [], []
    for i, f in enumerate(features):
        f_start, f_end = f["start"], f["end"]
        f_len = f_end - f_start + 1
        is_sub = False
        for j, parent in enumerate(features):
            if i == j:
                continue
            p_start, p_end = parent["start"], parent["end"]
            p_len = p_end - p_start + 1
            if f_start >= p_start and f_end <= p_end and f_len < SUBFEATURE_MAX_RATIO * p_len:
                is_sub = True
                break
        if is_sub:
            subs.append(f)
        else:
            main.append(f)
    return main, subs


def _parse_actual_annotation(
    query_record: SeqRecord,
    alias_lookup: dict,
    ref_feature_type: str,
) -> Tuple[List[Dict], List[Dict]]:
    """
    Parse actual annotation from query, normalized to canonical names.
    Returns (main_features, sub_features).
    """
    if ref_feature_type == "mat_peptide":
        features = parse_mat_peptides(query_record)
        if not features:
            features = parse_cds_features(query_record)
    else:
        features = parse_cds_features(query_record)

    if alias_lookup and features:
        features = apply_alias_to_features(features, alias_lookup)

    return _filter_subfeatures(features)


def _run_alignment(ref_record: SeqRecord, query_record: SeqRecord, tmp_dir: Path) -> Tuple[Dict, str]:
    """Run minimap2 with fallback presets. Returns (ref_to_query, preset_used)."""
    ref_fasta = tmp_dir / "ref.fa"
    query_fasta = tmp_dir / "query.fa"
    sam_path = tmp_dir / "alignment.sam"
    _write_fasta(ref_record, ref_fasta)
    _write_fasta(query_record, query_fasta)
    preset_used = run_minimap2_with_fallback(ref_fasta, query_fasta, sam_path, quiet=True)
    aln = get_primary_alignment(sam_path)
    return build_ref_to_query_map(aln), preset_used


# ---------------------------------------------------------------------------
# Comparison
# ---------------------------------------------------------------------------

def _compare(lifted_results: List[Dict], actual_features: List[Dict]) -> List[Dict]:
    """
    Match lifted results to actual annotation by canonical name.

    Matching rules:
    - Match by canonical_name (set by apply_ref_naming) or name
    - If length ratio of lifted vs actual is outside [MIN_LEN_RATIO, MAX_LEN_RATIO]
      → mark as size_mismatch, do not compute diffs
    """
    # Build lookup: canonical -> actual feature
    # If multiple features share same canonical (shouldn't after sub-feature filter),
    # keep the largest one (most likely the main feature).
    actual_by_canonical: Dict[str, Dict] = {}
    for f in actual_features:
        key = f["name"]
        existing = actual_by_canonical.get(key)
        if existing is None:
            actual_by_canonical[key] = f
        else:
            # Keep larger one
            if (f["end"] - f["start"]) > (existing["end"] - existing["start"]):
                actual_by_canonical[key] = f

    rows = []
    for r in lifted_results:
        match_key = r.get("canonical_name") or r["name"]
        actual = actual_by_canonical.get(match_key)

        lifted_start = r.get("start")
        lifted_end = r.get("end")
        lifted_len = (lifted_end - lifted_start + 1) if (lifted_start and lifted_end) else None

        if actual:
            actual_start = actual["start"]
            actual_end = actual["end"]
            actual_len = actual_end - actual_start + 1

            # Size mismatch check
            if lifted_len is not None:
                ratio = lifted_len / actual_len
                size_ok = MIN_LEN_RATIO <= ratio <= MAX_LEN_RATIO
            else:
                size_ok = False

            if size_ok:
                start_diff = (lifted_start - actual_start) if lifted_start is not None else None
                end_diff = (lifted_end - actual_end) if lifted_end is not None else None
                len_diff = (lifted_len - actual_len) if lifted_len is not None else None
                matched = True
                match_note = ""
            else:
                start_diff = end_diff = len_diff = None
                matched = False
                match_note = "size_mismatch"
        else:
            actual_start = actual_end = actual_len = None
            start_diff = end_diff = len_diff = None
            matched = False
            match_note = "not_in_actual"

        rows.append({
            "name": r["name"],
            "canonical_name": r.get("canonical_name", ""),
            "lifted_start": lifted_start,
            "lifted_end": lifted_end,
            "lifted_len": lifted_len,
            "actual_start": actual_start,
            "actual_end": actual_end,
            "actual_len": actual_len,
            "start_diff": start_diff,
            "end_diff": end_diff,
            "len_diff": len_diff,
            "status": r["status"],
            "rescue_offset": r.get("rescue_offset"),
            "matched": matched,
            "match_note": match_note,
        })

    return rows


def _accuracy_stats(comparison_rows: List[Dict]) -> Dict:
    total = len(comparison_rows)
    matched = sum(1 for r in comparison_rows if r["matched"])
    perfect = sum(
        1 for r in comparison_rows
        if r["matched"] and r["start_diff"] == 0 and r["end_diff"] == 0
    )
    rescued = sum(1 for r in comparison_rows if r["status"] == "ok_rescued")
    invalid = sum(1 for r in comparison_rows if r["status"] == "invalid_boundaries")
    unmapped = sum(1 for r in comparison_rows if r["status"] == "unmapped")
    size_mismatch = sum(1 for r in comparison_rows if r["match_note"] == "size_mismatch")
    not_in_actual = sum(1 for r in comparison_rows if r["match_note"] == "not_in_actual")

    return {
        "total": total,
        "matched": matched,
        "perfect": perfect,
        "rescued": rescued,
        "invalid_boundaries": invalid,
        "unmapped": unmapped,
        "size_mismatch": size_mismatch,
        "not_in_actual": not_in_actual,
    }


# ---------------------------------------------------------------------------
# Per-record processing
# ---------------------------------------------------------------------------

def process_record(ref_record, query_record, normalized_ref, ref_feature_type,
                   alias_lookup, output_dir):
    with tempfile.TemporaryDirectory() as tmp:
        try:
            ref_to_query, preset_used = _run_alignment(ref_record, query_record, Path(tmp))
        except Exception as e:
            print(f"  {query_record.id:<20}  [ALIGNMENT FAILED] {e}")
            return None

    validate_codons = (ref_feature_type == "CDS")
    lifted = extract_all_lifted(normalized_ref, query_record, ref_to_query,
                                validate_codons=validate_codons)

    actual_main, actual_subs = _parse_actual_annotation(query_record, alias_lookup, ref_feature_type)
    comparison = _compare(lifted, actual_main)
    stats = _accuracy_stats(comparison)

    total = stats["total"]
    perfect_pct = 100 * stats["perfect"] / total if total else 0

    extras = []
    if preset_used != "asm5":
        extras.append(f"preset={preset_used}")
    if stats["rescued"]:
        extras.append(f"rescued={stats['rescued']}")
    if stats["invalid_boundaries"]:
        extras.append(f"invalid={stats['invalid_boundaries']}")
    if stats["unmapped"]:
        extras.append(f"unmapped={stats['unmapped']}")
    if stats["size_mismatch"]:
        extras.append(f"size_mismatch={stats['size_mismatch']}")
    if stats["not_in_actual"]:
        extras.append(f"not_in_actual={stats['not_in_actual']}")
    if actual_subs:
        extras.append(f"subfeatures_skipped={len(actual_subs)}")

    extra_str = f"  [{', '.join(extras)}]" if extras else ""

    print(f"  {query_record.id:<20}  "
          f"perfect={stats['perfect']}/{total} ({perfect_pct:.0f}%)"
          f"{extra_str}")

    # Print non-perfect details
    for r in comparison:
        if r["matched"] and r["start_diff"] == 0 and r["end_diff"] == 0:
            continue
        ds = f"{r['start_diff']:+d}" if r["start_diff"] is not None else "n/a"
        de = f"{r['end_diff']:+d}" if r["end_diff"] is not None else "n/a"
        rescue = f" rescue={r['rescue_offset']:+d}" if r.get("rescue_offset") is not None else ""
        note = f" [{r['match_note']}]" if r["match_note"] else ""
        print(f"    -> {r['name']:<35} dStart={ds:>5} dEnd={de:>5}  {r['status']}{rescue}{note}")

    return comparison, stats


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def write_batch_tsv(all_rows, output_dir: Path, stem: str):
    tsv_path = output_dir / f"{stem}_batch_comparison.tsv"
    fieldnames = [
        "record_id", "name", "canonical_name",
        "lifted_start", "lifted_end", "lifted_len",
        "actual_start", "actual_end", "actual_len",
        "start_diff", "end_diff", "len_diff",
        "status", "rescue_offset", "matched", "match_note",
    ]
    with open(tsv_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(all_rows)
    return tsv_path


def write_summary_tsv(summary_rows, output_dir: Path, stem: str):
    tsv_path = output_dir / f"{stem}_batch_summary.tsv"
    fieldnames = [
        "record_id", "total", "matched", "perfect",
        "rescued", "invalid_boundaries", "unmapped",
        "size_mismatch", "not_in_actual", "perfect_pct",
    ]
    with open(tsv_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(summary_rows)
    return tsv_path


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description="Batch minimap accuracy test over annotated records in a GenBank file."
    )
    parser.add_argument("--reference", required=True)
    parser.add_argument("--query", required=True, help="Multi-record GenBank file.")
    parser.add_argument("--alias-registry", default="app/config/virus_alias_registry.json")
    parser.add_argument("--output-dir", default="output/minimap_batch")
    return parser.parse_args()


def main():
    args = parse_args()

    ref_path = Path(args.reference)
    query_path = Path(args.query)
    registry_path = Path(args.alias_registry)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    ref_record = load_single_genbank(ref_path)
    all_records = load_genbank_records(query_path)

    annotated = [r for r in all_records if _has_full_annotation(r)]
    skipped = len(all_records) - len(annotated)

    virus_name = get_detected_virus_name(ref_record, registry_path) or "Unknown"
    alias_config_path = detect_alias_config_for_record(ref_record, registry_path)
    alias_lookup = load_alias_lookup(alias_config_path) if alias_config_path else {}

    ref_features, ref_feature_type = _parse_ref_features(ref_record)
    normalized_ref = _normalize_ref_features(ref_features, alias_lookup)

    print(f"Virus     : {virus_name}")
    print(f"Reference : {ref_record.id}  ({ref_feature_type}, {len(ref_features)} features)")
    print(f"Query file: {query_path.name}  ({len(all_records)} records, "
          f"{len(annotated)} annotated, {skipped} skipped)\n")

    all_comparison_rows = []
    summary_rows = []

    for query_record in annotated:
        result = process_record(
            ref_record, query_record,
            normalized_ref, ref_feature_type,
            alias_lookup, output_dir,
        )
        if result is None:
            continue

        comparison, stats = result
        for row in comparison:
            all_comparison_rows.append({"record_id": query_record.id, **row})

        perfect_pct = round(100 * stats["perfect"] / stats["total"], 1) if stats["total"] else 0
        summary_rows.append({"record_id": query_record.id, **stats, "perfect_pct": perfect_pct})

    if not summary_rows:
        print("\nNo records processed.")
        return

    total_features = sum(r["total"] for r in summary_rows)
    total_perfect = sum(r["perfect"] for r in summary_rows)
    total_rescued = sum(r["rescued"] for r in summary_rows)
    total_invalid = sum(r["invalid_boundaries"] for r in summary_rows)
    total_size_mm = sum(r["size_mismatch"] for r in summary_rows)
    total_not_actual = sum(r["not_in_actual"] for r in summary_rows)

    print(f"\n{'='*60}")
    print(f"AGGREGATE ({len(summary_rows)} records, {total_features} feature-lifts)")
    print(f"  Perfect (dStart=0, dEnd=0) : {total_perfect}/{total_features} ({100*total_perfect/total_features:.1f}%)")
    if total_rescued:
        print(f"  Rescued (ok_rescued)       : {total_rescued}")
    if total_invalid:
        print(f"  Invalid boundaries         : {total_invalid}")
    if total_size_mm:
        print(f"  Size mismatch (skipped)    : {total_size_mm}")
    if total_not_actual:
        print(f"  Not found in actual        : {total_not_actual}")

    batch_tsv = write_batch_tsv(all_comparison_rows, output_dir, query_path.stem)
    summary_tsv = write_summary_tsv(summary_rows, output_dir, query_path.stem)
    print(f"\nBatch TSV   : {batch_tsv}")
    print(f"Summary TSV : {summary_tsv}")


if __name__ == "__main__":
    main()
