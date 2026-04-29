import argparse
from pathlib import Path

from app.src.alias.alias_payload import run_alias_pipeline, save_payload
from app.src.io.genbank_parser import load_single_genbank, load_genbank_records


def parse_args():
    parser = argparse.ArgumentParser(
        description="Test Case 1 alias normalization pipeline: reference + query -> normalized names + LLM payload."
    )
    parser.add_argument("--reference", required=True, help="Reference GenBank file (single record).")
    parser.add_argument("--query", required=True, help="Query GenBank file (single or multi-record).")
    parser.add_argument("--alias-registry", default="app/config/virus_alias_registry.json")
    parser.add_argument("--output-dir", default="output/alias_payloads")
    parser.add_argument(
        "--canonical",
        action="store_true",
        help="Output canonical names instead of reference names. Default: use reference naming.",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    ref_record = load_single_genbank(Path(args.reference))
    query_records = load_genbank_records(Path(args.query))

    print(f"Reference : {ref_record.id}")
    print(f"Query     : {Path(args.query).name} ({len(query_records)} records)\n")

    results, payload = run_alias_pipeline(
        ref_record=ref_record,
        query_records=query_records,
        registry_path=Path(args.alias_registry),
        use_ref_naming=not args.canonical,
    )

    # Per-record output
    for r in results:
        status = r["status"]
        rid = r["record_id"]

        if status == "no_annotation":
            print(f"  [{rid}] no annotation -> minimap case")
        elif status == "new_virus":
            print(f"  [{rid}] virus not in registry -> build_alias_map ({r['feature_type']}, {r['total']} features)")
        elif status == "all_resolved":
            print(f"  [{rid}] {r['feature_type']} | {r['resolved']}/{r['total']} resolved | all good")
        elif status == "has_unresolved":
            unresolved = r["total"] - r["resolved"]
            print(f"  [{rid}] {r['feature_type']} | resolved: {r['resolved']}/{r['total']} | unresolved: {unresolved}")

    # Summary
    counts = {}
    for r in results:
        counts[r["status"]] = counts.get(r["status"], 0) + 1

    print(f"\n--- Summary ---")
    print(f"All resolved    : {counts.get('all_resolved', 0)}")
    print(f"Has unresolved  : {counts.get('has_unresolved', 0)}")
    print(f"New virus       : {counts.get('new_virus', 0)}")
    print(f"No annotation   : {counts.get('no_annotation', 0)}")

    if payload is None:
        print("\nAll features resolved. No payload needed.")
        return

    output_dir = Path(args.output_dir)
    query_stem = Path(args.query).stem
    output_path = output_dir / f"{query_stem}_alias_payload.json"
    save_payload(payload, output_path)
    print(f"\nPayload saved : {output_path}")


if __name__ == "__main__":
    main()
