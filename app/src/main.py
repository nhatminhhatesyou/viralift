import argparse
from pathlib import Path

from tqdm import tqdm

from app.src.alias.alias_registry import DEFAULT_REGISTRY_PATH
from app.src.features.ref_loader import prepare_reference_features
from app.src.io.genbank_parser import load_single_genbank, load_genbank_records
from app.src.io.result_writer import write_results_tsv
from app.src.lifting.base import STATUS_LABELS
from app.src.pipeline import PipelineConfig, run_pipeline


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
    parser.add_argument("--rescue-window", type=int, default=200,
                        help="Window size (bp) for start codon rescue. Default: 200")
    parser.add_argument("--alias-config",
                        help="Optional path to alias JSON config file.")
    parser.add_argument("--alias-registry",
                        default=str(DEFAULT_REGISTRY_PATH),
                        help="Path to virus alias registry. "
                             "Default: the registry shipped with the package.")
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

    ref_features, ref_feature_type, alias_config_path, detected_virus_name, alias_lookup = (
        prepare_reference_features(
            ref_record=ref_record,
            alias_config_arg=args.alias_config,
            alias_registry_arg=args.alias_registry,
        )
    )

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

    iterator = tqdm(query_records, desc="Processing records", unit="record", dynamic_ncols=True)

    def update_progress(index: int, total: int, record_id: str) -> None:
        if index >= total:
            iterator.update(total - iterator.n)
            return
        if not args.quiet:
            iterator.set_postfix_str(record_id)
        iterator.update(index - iterator.n)

    run_result = run_pipeline(
        ref_record=ref_record,
        query_records=query_records,
        ref_features=ref_features,
        ref_feature_type=ref_feature_type,
        alias_lookup=alias_lookup,
        config=PipelineConfig(
            min_coverage=args.min_coverage,
            min_identity=args.min_identity,
            evalue=args.evalue,
            rescue_window=args.rescue_window,
        ),
        progress_callback=update_progress,
    )

    tsv_out = outdir / "extracted_cds.tsv"
    write_results_tsv(run_result.all_results, tsv_out)

    print("\nRun summary")
    print(f"  Query records processed : {len(query_records)}")
    print(f"    Direct (annotated)    : {run_result.direct_count}")
    print(f"    Lifted (tblastn)      : {run_result.lifted_count}")
    for status, count in run_result.summary.items():
        label = STATUS_LABELS.get(status, status)
        print(f"  {label:<24}: {count}")
    print(f"\n  TSV : {tsv_out}")


if __name__ == "__main__":
    main()
