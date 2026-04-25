from pathlib import Path
from typing import Dict, List, Optional, Tuple
import argparse
import csv

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm

from app.src.annotation.alias_registry import (
    detect_alias_config_for_record,
    get_detected_virus_name,
)
from app.src.annotation.annotation_strategy import choose_strategy
from app.src.annotation.gene_alias import (
    apply_alias_to_features,
    load_alias_lookup,
)
from app.src.io.genbank_parser import (
    load_single_genbank,
    load_genbank_records,
    parse_cds_features,
    parse_mat_peptides,
)
from app.src.lifting.tblastn_lifter import lift_all_tblastn


"""
Module: main.py

Purpose:
    CLI entry point for reference-guided viral CDS transfer using tblastn.

Workflow:
    1. Load reference and query GenBank records
    2. Parse reference CDS features
    3. Optionally normalize feature names using alias config
    4. For each query genome: translate ref proteins and search via tblastn
    5. Lift CDS coordinates from reference to query using protein homology
    6. Extract transferred sequences with codon validation and rescue
    7. Write TSV annotation output

Alias behavior:
    - If --alias-config is provided, that config is used directly
    - Otherwise, the program tries to auto-detect the correct alias config
      using config/virus_alias_registry.json
    - If no alias config is found, raw feature names are preserved
"""


# ---------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    """
    Parse command-line arguments.

    Returns:
        Parsed CLI arguments
    """
    parser = argparse.ArgumentParser(
        prog="viralift",
        description="Reference-guided viral CDS transfer using tblastn.",
        epilog=(
            "Examples:\n"
            "  python -m app.src.main "
            "--reference data/PRRS_ref_test.gb "
            "--query data/PRRSV_test.gb "
            "--output output/prrsv_multi\n\n"
            "  python -m app.src.main "
            "--reference data/PRRS_ref_test.gb "
            "--query data/PRRSV_test.gb "
            "--output output/prrsv_multi "
            "--alias-config config/prrsv_alias.json\n"
        ),
        formatter_class=argparse.RawTextHelpFormatter,
    )

    parser.add_argument(
        "--reference",
        required=True,
        help="Reference GenBank file (single record).",
    )
    parser.add_argument(
        "--query",
        required=True,
        help="Query GenBank file (one or more records).",
    )
    parser.add_argument(
        "--output",
        default="output/run",
        help="Output directory. Default: output/run",
    )
    parser.add_argument(
        "--feature-type",
        default="CDS",
        choices=["CDS"],
        help="Feature type to transfer. Default: CDS",
    )
    parser.add_argument(
        "--min-coverage",
        type=float,
        default=0.5,
        help="Minimum protein coverage threshold for lifted features. Default: 0.5",
    )
    parser.add_argument(
        "--min-identity",
        type=float,
        default=0.3,
        help="Minimum protein identity threshold for lifted features. Default: 0.3",
    )
    parser.add_argument(
        "--evalue",
        type=float,
        default=1e-5,
        help="E-value threshold for tblastn search. Default: 1e-5",
    )
    parser.add_argument(
        "--rescue-window",
        type=int,
        default=50,
        help="Window size (bp) for start codon rescue. Default: 50",
    )
    parser.add_argument(
        "--alias-config",
        required=False,
        help="Optional path to alias JSON config file.",
    )
    parser.add_argument(
        "--alias-registry",
        default="config/virus_alias_registry.json",
        help="Path to virus alias registry JSON file. Default: config/virus_alias_registry.json",
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="Reduce console output but still show progress and final summary.",
    )

    return parser.parse_args()


# ---------------------------------------------------------------------
# Alias utilities
# ---------------------------------------------------------------------

def prepare_reference_features(
    ref_record: SeqRecord,
    alias_config_arg: Optional[str],
    alias_registry_arg: str,
) -> Tuple[List[Dict], Optional[Path], Optional[str]]:
    """
    Parse reference features and optionally normalize their names using alias config.

    Uses annotation strategy to determine whether to parse CDS or mat_peptide features.

    Priority:
        1. Use user-provided alias config if available
        2. Otherwise auto-detect config using alias registry
        3. If detection fails, keep raw names

    Args:
        ref_record: Reference SeqRecord
        alias_config_arg: Optional CLI value for --alias-config
        alias_registry_arg: CLI value for --alias-registry

    Returns:
        Tuple of:
            - ref_features: Parsed (and optionally alias-normalized) reference features
            - alias_config_path: Path to alias config actually used, or None
            - detected_virus_name: Virus name detected from registry, or None
    """
    # Determine feature type based on annotation strategy
    _, feature_type = choose_strategy(ref_record)

    if feature_type == "mat_peptide":
        ref_features = parse_mat_peptides(ref_record)
    elif feature_type == "CDS":
        ref_features = parse_cds_features(ref_record)
    else:
        raise ValueError("Reference record has no CDS or mat_peptide features.")

    if not ref_features:
        raise ValueError(f"Reference record has no {feature_type} features.")

    alias_config_path: Optional[Path] = None
    detected_virus_name: Optional[str] = None

    # Case 1: user explicitly provides alias config
    if alias_config_arg:
        alias_config_path = Path(alias_config_arg)

    # Case 2: auto-detect from registry
    else:
        registry_path = Path(alias_registry_arg)

        try:
            alias_config_path = detect_alias_config_for_record(ref_record, registry_path)

            if alias_config_path is not None:
                detected_virus_name = get_detected_virus_name(ref_record, registry_path)

        except FileNotFoundError:
            alias_config_path = None
        except ValueError:
            alias_config_path = None

    # Apply alias normalization if config exists
    if alias_config_path is not None:
        alias_lookup = load_alias_lookup(alias_config_path)
        ref_features = apply_alias_to_features(ref_features, alias_lookup)

    return ref_features, alias_config_path, detected_virus_name


# ---------------------------------------------------------------------
# Output writers
# ---------------------------------------------------------------------

def write_results_fasta(all_results: List[Tuple[str, List]], out_path: Path) -> None:
    """
    Write successfully extracted CDS sequences from all query records to a FASTA file.

    Args:
        all_results: List of (query_id, lifted_features)
        out_path: Output FASTA path
    """
    records: List[SeqRecord] = []

    for query_id, results in all_results:
        for lifted in results:
            if lifted.status not in ("ok", "ok_rescued"):
                continue

            if not lifted.sequence:
                continue

            record_id = f"{query_id}|{lifted.name}|{lifted.method}"
            records.append(
                SeqRecord(
                    Seq(lifted.sequence),
                    id=record_id,
                    description="",
                )
            )

    SeqIO.write(records, str(out_path), "fasta")


def write_results_tsv(all_results: List[Tuple[str, List]], out_path: Path) -> None:
    """
    Write extracted feature results from all query records to a TSV file.

    Args:
        all_results: List of (query_id, lifted_features)
        out_path: Output TSV path
    """
    rows: List[Dict] = []

    for query_id, results in all_results:
        for lifted in results:
            rows.append(
                {
                    "query_id": query_id,
                    "name": lifted.name,
                    "canonical_name": lifted.canonical_name or "",
                    "ref_start": lifted.ref_start,
                    "ref_end": lifted.ref_end,
                    "start": lifted.query_start,
                    "end": lifted.query_end,
                    "strand": lifted.strand,
                    "method": lifted.method,
                    "status": lifted.status,
                    "coverage": lifted.coverage,
                    "identity": lifted.identity,
                    "score": lifted.score,
                    "has_start_codon": lifted.has_start_codon,
                    "has_stop_codon": lifted.has_stop_codon,
                    "rescue_offset": lifted.rescue_offset,
                    "length": len(lifted.sequence) if lifted.sequence else None,
                }
            )

    if not rows:
        return

    with open(out_path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


# ---------------------------------------------------------------------
# Summary utilities
# ---------------------------------------------------------------------

def summarize_counts(all_results: List[Tuple[str, List]]) -> Dict[str, int]:
    """
    Summarize result statuses across all processed records.

    Args:
        all_results: List of (query_id, lifted_features)

    Returns:
        Dictionary with counts by status
    """
    summary = {
        "ok": 0,
        "ok_rescued": 0,
        "invalid_boundaries": 0,
        "low_coverage": 0,
        "no_hit": 0,
        "translation_fail": 0,
    }

    for _, results in all_results:
        for lifted in results:
            status = lifted.status
            if status in summary:
                summary[status] += 1

    return summary


# ---------------------------------------------------------------------
# Core processing
# ---------------------------------------------------------------------



def process_one_query_record(
    ref_record: SeqRecord,
    query_record: SeqRecord,
    ref_cds: List[Dict],
    ref_feature_type: str,
    min_coverage: float,
    min_identity: float = 0.3,
    evalue: float = 1e-5,
    rescue_window: int = 50,
    quiet: bool = False,
) -> List:
    """
    Process one query genome using reference-guided tblastn transfer.

    Steps:
        1. Determine if codon validation is needed based on feature type
        2. Use tblastn to lift all reference features to query genome
        3. Return LiftedFeature objects with validation results

    Args:
        ref_record: Reference genome record
        query_record: Query genome record
        ref_cds: Parsed reference CDS features
        ref_feature_type: Type of reference features (CDS or mat_peptide)
        min_coverage: Minimum accepted protein coverage
        min_identity: Minimum accepted protein identity
        evalue: E-value threshold for tblastn
        rescue_window: Window size for start codon rescue
        quiet: Whether to reduce console output

    Returns:
        List of LiftedFeature objects
    """
    # mat_peptide features don't need codon validation
    validate_codons = (ref_feature_type == "CDS")

    results = lift_all_tblastn(
        ref_features=ref_cds,
        ref_record=ref_record,
        query_record=query_record,
        min_coverage=min_coverage,
        min_identity=min_identity,
        evalue=evalue,
        rescue_window=rescue_window,
        validate_codons=validate_codons,
    )

    return results


# ---------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------

def main() -> None:
    """Run the viral CDS transfer pipeline."""
    args = parse_args()

    outdir = Path(args.output)
    outdir.mkdir(parents=True, exist_ok=True)

    ref_record = load_single_genbank(Path(args.reference))
    query_records = load_genbank_records(Path(args.query))

    if not query_records:
        raise ValueError("No query records found.")

    ref_features, alias_config_path, detected_virus_name = prepare_reference_features(
        ref_record=ref_record,
        alias_config_arg=args.alias_config,
        alias_registry_arg=args.alias_registry,
    )

    # Detect reference feature type for validation
    _, ref_feature_type = choose_strategy(ref_record)

    print("ViraLift")
    print(f"  Reference record : {ref_record.id}")
    print(f"  Feature type     : {ref_feature_type}")
    print(f"  Reference features : {len(ref_features)}")
    print(f"  Query records    : {len(query_records)}")
    print(f"  Min coverage     : {args.min_coverage}")
    print(f"  Min identity     : {args.min_identity}")
    print(f"  E-value          : {args.evalue}")
    print(f"  Output folder    : {outdir}")

    if args.alias_config:
        print(f"  Alias config     : {alias_config_path} (user-provided)")
    elif alias_config_path:
        print(f"  Alias config     : {alias_config_path} (auto-detected)")
        if detected_virus_name:
            print(f"  Detected virus   : {detected_virus_name}")
    else:
        print("  Alias config     : none (using raw names)")

    all_results: List[Tuple[str, List[Dict]]] = []

    iterator = tqdm(
        query_records,
        desc="Processing records",
        unit="record",
        ncols=90,
    )

    for query_record in iterator:
        if not args.quiet:
            iterator.set_postfix_str(query_record.id)

        results = process_one_query_record(
            ref_record=ref_record,
            query_record=query_record,
            ref_cds=ref_features,
            ref_feature_type=ref_feature_type,
            min_coverage=args.min_coverage,
            min_identity=args.min_identity,
            evalue=args.evalue,
            rescue_window=args.rescue_window,
            quiet=args.quiet,
        )
        all_results.append((query_record.id, results))

    # fasta_out = outdir / "extracted_cds.fasta"
    tsv_out = outdir / "extracted_cds.tsv"

    # write_results_fasta(all_results, fasta_out)  # FASTA output not essential for gene normalization
    write_results_tsv(all_results, tsv_out)

    summary = summarize_counts(all_results)

    print("\nRun summary")
    print(f"  Query records processed : {len(query_records)}")
    print(f"  OK                      : {summary['ok']}")
    print(f"  OK (rescued)            : {summary['ok_rescued']}")
    print(f"  Invalid boundaries      : {summary['invalid_boundaries']}")
    print(f"  Low coverage            : {summary['low_coverage']}")
    print(f"  No hit                  : {summary['no_hit']}")
    print(f"  Translation fail        : {summary['translation_fail']}")

    print("\nOutput files")
    print(f"  TSV   : {tsv_out}")


if __name__ == "__main__":
    main()