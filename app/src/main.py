import argparse
from pathlib import Path
from typing import List, Tuple

from tqdm import tqdm

from app.src.features.annotation_strategy import get_best_feature_type, get_feature_type, get_strategy
from app.src.features.direct_extractor import direct_extract_with_alias
from app.src.features.ref_loader import prepare_reference_features
from app.src.io.genbank_parser import load_single_genbank, load_genbank_records
from app.src.io.result_writer import summarize_counts, write_results_tsv
from app.src.lifting.tblastn_lifter import process_one_query_record


"""
Module: main.py

Purpose:
    CLI entry point for ViraLift — reference-guided viral gene annotation transfer.

Workflow:
    1. Parse CLI arguments
    2. Load reference + query GenBank records
    3. Prepare reference features (detect type, load alias config)
    4. For each query record: decide strategy (direct / tblastn) and process
    5. Write TSV output + print run summary
"""


# ---------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
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
    parser.add_argument("--reference", required=True,
                        help="Reference GenBank file (single record).")
    parser.add_argument("--query",     required=True,
                        help="Query GenBank file (one or more records).")
    parser.add_argument("--output",    default="output/run",
                        help="Output directory. Default: output/run")
    parser.add_argument("--min-coverage", type=float, default=0.5,
                        help="Minimum protein coverage threshold. Default: 0.5")
    parser.add_argument("--min-identity", type=float, default=0.3,
                        help="Minimum protein identity threshold. Default: 0.3")
    parser.add_argument("--evalue",       type=float, default=1e-5,
                        help="E-value threshold for tblastn. Default: 1e-5")
    parser.add_argument("--rescue-window", type=int, default=50,
                        help="Window size (bp) for start codon rescue. Default: 50")
    parser.add_argument("--alias-config",
                        help="Optional path to alias JSON config file.")
    parser.add_argument("--alias-registry",
                        default="app/config/virus_alias_registry.json",
                        help="Path to virus alias registry. Default: app/config/virus_alias_registry.json")
    parser.add_argument("--quiet", action="store_true",
                        help="Reduce console output.")
    return parser.parse_args()


# ---------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------

def main() -> None:
    """Run the ViraLift annotation transfer pipeline."""
    args = parse_args()

    outdir = Path(args.output)
    outdir.mkdir(parents=True, exist_ok=True)

    ref_record     = load_single_genbank(Path(args.reference))
    query_records  = load_genbank_records(Path(args.query))

    if not query_records:
        raise ValueError("No query records found.")

    ref_features, alias_config_path, detected_virus_name, alias_lookup = (
        prepare_reference_features(
            ref_record=ref_record,
            alias_config_arg=args.alias_config,
            alias_registry_arg=args.alias_registry,
        )
    )
    ref_feature_type = get_best_feature_type(ref_record, alias_lookup) or get_feature_type(ref_record)

    print("ViraLift")
    print(f"  Reference record   : {ref_record.id}")
    print(f"  Feature type       : {ref_feature_type}")
    print(f"  Reference features : {len(ref_features)}")
    print(f"  Query records      : {len(query_records)}")
    print(f"  Min coverage       : {args.min_coverage}")
    print(f"  Min identity       : {args.min_identity}")
    print(f"  E-value            : {args.evalue}")
    print(f"  Output folder      : {outdir}")

    if args.alias_config:
        print(f"  Alias config       : {alias_config_path} (user-provided)")
    elif alias_config_path:
        print(f"  Alias config       : {alias_config_path} (auto-detected)")
        if detected_virus_name:
            print(f"  Detected virus     : {detected_virus_name}")
    else:
        print("  Alias config       : none (using raw names)")

    all_results: List[Tuple] = []
    direct_count = lifted_count = 0

    iterator = tqdm(query_records, desc="Processing records", unit="record", dynamic_ncols=True)

    for query_record in iterator:
        if not args.quiet:
            iterator.set_postfix_str(query_record.id)

        strategy = get_strategy(query_record, ref_feature_type, alias_lookup)

        if strategy == "direct":
            results = direct_extract_with_alias(
                query_record=query_record,
                query_feature_type=get_best_feature_type(query_record, alias_lookup),
                ref_features=ref_features,
                alias_lookup=alias_lookup,
            )
            direct_count += 1
        else:
            results = process_one_query_record(
                ref_record=ref_record,
                query_record=query_record,
                ref_cds=ref_features,
                ref_feature_type=ref_feature_type,
                min_coverage=args.min_coverage,
                min_identity=args.min_identity,
                evalue=args.evalue,
                rescue_window=args.rescue_window,
            )
            lifted_count += 1

        all_results.append((query_record.id, results))

    tsv_out = outdir / "extracted_cds.tsv"
    write_results_tsv(all_results, tsv_out)

    summary = summarize_counts(all_results)

    print("\nRun summary")
    print(f"  Query records processed : {len(query_records)}")
    print(f"    Direct (annotated)    : {direct_count}")
    print(f"    Lifted (tblastn)      : {lifted_count}")
    print(f"  OK                      : {summary['ok']}")
    print(f"  OK (rescued)            : {summary['ok_rescued']}")
    print(f"  Invalid boundaries      : {summary['invalid_boundaries']}")
    print(f"  Low coverage            : {summary['low_coverage']}")
    print(f"  No hit                  : {summary['no_hit']}")
    print(f"  Translation fail        : {summary['translation_fail']}")
    print(f"  Unresolved names        : {summary['unresolved_name']}")
    print(f"  Ambiguous names         : {summary['ambiguous_name']}")
    print(f"  Not in reference        : {summary['not_in_reference']}")
    print(f"\n  TSV : {tsv_out}")


if __name__ == "__main__":
    main()
