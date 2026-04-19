import argparse
from pathlib import Path

from app.src.annotation.alias_payload import (
    build_map_aliases_payload,
    build_new_virus_payload,
    get_unresolved_features,
    save_payload,
)
from app.src.annotation.alias_registry import (
    detect_alias_config_for_record,
    get_detected_virus_name,
)
from app.src.annotation.annotation_strategy import choose_strategy
from app.src.annotation.gene_alias import apply_alias_to_features, load_alias_lookup
from app.src.io.genbank_parser import (
    load_single_genbank,
    parse_cds_features,
    parse_mat_peptides,
)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract and save alias payload for LLM processing."
    )
    parser.add_argument("--input", required=True, help="Input GenBank file (single record).")
    parser.add_argument("--alias-registry", default="app/config/virus_alias_registry.json")
    parser.add_argument("--output-dir", default="output/alias_payloads")
    return parser.parse_args()


def main():
    args = parse_args()

    input_path = Path(args.input)
    registry_path = Path(args.alias_registry)
    output_dir = Path(args.output_dir)

    record = load_single_genbank(input_path)

    strategy, feature_type = choose_strategy(record)

    if feature_type == "mat_peptide":
        raw_features = parse_mat_peptides(record)
    elif feature_type == "CDS":
        raw_features = parse_cds_features(record)
    else:
        print("No CDS or mat_peptide features found. Nothing to extract.")
        return

    print(f"Record       : {record.id}")
    print(f"Strategy     : {strategy} ({feature_type})")
    print(f"Features     : {len(raw_features)}")

    alias_config_path = detect_alias_config_for_record(record, registry_path)
    virus_name = get_detected_virus_name(record, registry_path)

    if alias_config_path is None:
        # Virus not in registry -> build_alias_map payload
        print(f"Virus        : not in registry -> task: build_alias_map")

        payload = build_new_virus_payload(
            record=record,
            all_features=raw_features,
            feature_type=feature_type,
        )
    else:
        # Virus known -> apply alias, find unresolved -> map_aliases payload
        print(f"Virus        : {virus_name}")
        print(f"Alias config : {alias_config_path}")

        alias_lookup = load_alias_lookup(alias_config_path)
        normalized = apply_alias_to_features(raw_features, alias_lookup)
        unresolved = get_unresolved_features(normalized)

        resolved_count = len(normalized) - len(unresolved)
        print(f"Resolved     : {resolved_count}/{len(normalized)}")
        print(f"Unresolved   : {len(unresolved)}")

        if not unresolved:
            print("\nAll features resolved. No payload needed.")
            return

        existing_canonical_names = list(
            load_alias_lookup(alias_config_path).values()
        )
        # deduplicate while preserving order
        seen = set()
        canonical_names = []
        for name in existing_canonical_names:
            if name not in seen:
                seen.add(name)
                canonical_names.append(name)

        payload = build_map_aliases_payload(
            record=record,
            unresolved_features=unresolved,
            existing_canonical_names=canonical_names,
            feature_type=feature_type,
        )

    output_path = output_dir / f"{record.id}_alias_payload.json"
    save_payload(payload, output_path)
    print(f"\nPayload saved: {output_path}")


if __name__ == "__main__":
    main()
