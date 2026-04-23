import argparse
import json
from pathlib import Path

from app.src.annotation.alias_payload import (
    build_map_aliases_payload,
    build_new_virus_payload,
    get_unresolved_features,
)
from app.src.annotation.alias_registry import (
    detect_alias_config_for_record,
    get_detected_virus_name,
)
from app.src.annotation.annotation_strategy import choose_strategy
from app.src.annotation.gene_alias import apply_alias_to_features, load_alias_lookup
from app.src.io.genbank_parser import (
    load_genbank_records,
    parse_cds_features,
    parse_mat_peptides,
)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Check alias coverage for all records in a GenBank file and save a single payload for LLM."
    )
    parser.add_argument("--input", required=True, help="Input GenBank file (single or multi-record).")
    parser.add_argument("--alias-registry", default="app/config/virus_alias_registry.json")
    parser.add_argument("--output-dir", default="output/alias_payloads")
    return parser.parse_args()


def process_record(record, registry_path):
    """
    Process one record. Returns a payload dict if unresolved features exist, else None.
    Also returns a status string for summary.
    """
    strategy, feature_type = choose_strategy(record)

    if feature_type == "mat_peptide":
        raw_features = parse_mat_peptides(record)
    elif feature_type == "CDS":
        raw_features = parse_cds_features(record)
    else:
        print(f"  [{record.id}] no annotation -> skip (minimap case)")
        return None, "no_annotation"

    alias_config_path = detect_alias_config_for_record(record, registry_path)
    virus_name = get_detected_virus_name(record, registry_path)

    if alias_config_path is None:
        print(f"  [{record.id}] virus not in registry -> build_alias_map ({feature_type}, {len(raw_features)} features)")
        payload = build_new_virus_payload(
            record=record,
            all_features=raw_features,
            feature_type=feature_type,
        )
        return payload, "new_virus"

    alias_lookup = load_alias_lookup(alias_config_path)
    normalized = apply_alias_to_features(raw_features, alias_lookup)
    unresolved = get_unresolved_features(normalized)
    resolved_count = len(normalized) - len(unresolved)

    if not unresolved:
        print(f"  [{record.id}] {virus_name} | {feature_type} | {resolved_count}/{len(normalized)} resolved | all good")
        return None, "all_resolved"

    print(f"  [{record.id}] {virus_name} | {feature_type} | resolved: {resolved_count}/{len(normalized)} | unresolved: {len(unresolved)}")

    seen = set()
    canonical_names = []
    for name in alias_lookup.values():
        if name not in seen:
            seen.add(name)
            canonical_names.append(name)

    payload = build_map_aliases_payload(
        record=record,
        unresolved_features=unresolved,
        existing_canonical_names=canonical_names,
        feature_type=feature_type,
    )
    return payload, "has_unresolved"


def main():
    args = parse_args()

    input_path = Path(args.input)
    registry_path = Path(args.alias_registry)
    output_dir = Path(args.output_dir)

    records = load_genbank_records(input_path)
    print(f"File    : {input_path}")
    print(f"Records : {len(records)}\n")

    counts = {"all_resolved": 0, "has_unresolved": 0, "new_virus": 0, "no_annotation": 0}
    payloads = []

    for record in records:
        payload, status = process_record(record, registry_path)
        counts[status] += 1
        if payload is not None:
            payloads.append(payload)

    print(f"\n--- Summary ---")
    print(f"All resolved    : {counts['all_resolved']}")
    print(f"Has unresolved  : {counts['has_unresolved']}")
    print(f"New virus       : {counts['new_virus']}")
    print(f"No annotation   : {counts['no_annotation']}")

    if not payloads:
        print("\nAll records resolved. No payload needed.")
        return

    output = {
        "source_file": input_path.name,
        "total_records": len(records),
        "records_with_unresolved": len(payloads),
        "records": payloads,
    }

    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / f"{input_path.stem}_alias_payload.json"

    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(output, f, indent=2, ensure_ascii=False)

    print(f"\nPayload saved : {output_path}")


if __name__ == "__main__":
    main()
