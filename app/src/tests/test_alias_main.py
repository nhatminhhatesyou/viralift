from pathlib import Path
import argparse

from app.src.annotation.alias_registry import (
    detect_alias_config_for_record,
    get_detected_virus_name,
)
from app.src.annotation.gene_alias import (
    load_alias_lookup,
    apply_alias_to_features,
)
from app.src.annotation.annotation_strategy import choose_strategy
from app.src.io.genbank_parser import (
    load_single_genbank,
    parse_cds_features,
    parse_mat_peptides,
)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Test alias normalization on a real GenBank reference record."
    )

    parser.add_argument(
        "--input",
        required=True,
        help="Input GenBank file (single record). Example: app/data/PRRS_ref_test.gb",
    )

    parser.add_argument(
        "--alias-config",
        required=False,
        help="Optional alias JSON config file. Example: app/config/prrsv_alias.json",
    )

    parser.add_argument(
        "--alias-registry",
        default="app/config/virus_alias_registry.json",
        help="Path to virus alias registry JSON file. Default: app/config/virus_alias_registry.json",
    )

    return parser.parse_args()


def summarize_alias_results(features):
    """
    Summarize alias normalization results.

    Returns:
        Dictionary with:
        - total
        - alias_count
        - raw_count
        - renamed_items
    """
    alias_count = 0
    raw_count = 0
    renamed_items = []

    for feature in features:
        raw_name = feature.get("raw_name")
        final_name = feature.get("name")
        name_source = feature.get("name_source", "raw")

        if name_source == "alias":
            alias_count += 1
            renamed_items.append((raw_name, final_name))
        else:
            raw_count += 1

    return {
        "total": len(features),
        "alias_count": alias_count,
        "raw_count": raw_count,
        "renamed_items": renamed_items,
    }


def resolve_alias_config(record, args):
    """
    Resolve which alias config to use.

    Priority:
        1. user-provided --alias-config
        2. auto-detect from registry
        3. no alias config
    """
    if args.alias_config:
        return Path(args.alias_config), None, "manual"

    registry_path = Path(args.alias_registry)

    alias_config_path = detect_alias_config_for_record(record, registry_path)
    detected_virus_name = get_detected_virus_name(record, registry_path)

    if alias_config_path is not None:
        return alias_config_path, detected_virus_name, "auto"

    return None, None, "none"


def main():
    args = parse_args()

    input_path = Path(args.input)
    record = load_single_genbank(input_path)

    strategy, feature_type = choose_strategy(record)

    if feature_type == "mat_peptide":
        raw_features = parse_mat_peptides(record)
    elif feature_type == "CDS":
        raw_features = parse_cds_features(record)
    else:
        raw_features = []

    if not raw_features:
        raise ValueError(f"No {feature_type or 'CDS/mat_peptide'} features found in input file.")

    alias_config_path, detected_virus_name, mode = resolve_alias_config(record, args)

    if alias_config_path is not None:
        alias_lookup = load_alias_lookup(alias_config_path)
        normalized_features = apply_alias_to_features(raw_features, alias_lookup)
    else:
        normalized_features = []
        for feature in raw_features:
            new_feature = feature.copy()
            new_feature["raw_name"] = feature.get("name")
            new_feature["name"] = feature.get("name")
            new_feature["name_source"] = "raw"
            normalized_features.append(new_feature)

    summary = summarize_alias_results(normalized_features)

    print("\n=== ALIAS NORMALIZATION TEST ===\n")
    print(f"Input file        : {input_path}")
    print(f"Record ID         : {record.id}")
    print(f"Record description: {record.description}")

    if mode == "manual":
        print("Alias mode        : manual")
        print(f"Alias config      : {alias_config_path}")
    elif mode == "auto":
        print("Alias mode        : auto-detect")
        print(f"Detected virus    : {detected_virus_name}")
        print(f"Alias config      : {alias_config_path}")
    else:
        print("Alias mode        : none")
        print("Alias config      : not found")
        print(f"Registry path     : {args.alias_registry}")

    print(f"Strategy          : {strategy} ({feature_type})")
    print(f"Total features    : {summary['total']}")
    print(f"Alias matched     : {summary['alias_count']}")
    print(f"Kept as raw       : {summary['raw_count']}")

    print("\n=== FEATURE NAME MAPPING ===\n")
    for i, feature in enumerate(normalized_features, start=1):
        raw_name = feature.get("raw_name")
        final_name = feature.get("name")
        name_source = feature.get("name_source")
        start = feature.get("start")
        end = feature.get("end")

        print(f"[{i}] {raw_name} -> {final_name} ({name_source}) [{start}-{end}]")

    print("\n=== ALIAS-HIT FEATURES ===\n")
    if summary["renamed_items"]:
        for i, (raw_name, final_name) in enumerate(summary["renamed_items"], start=1):
            print(f"[{i}] {raw_name} -> {final_name}")
    else:
        print("No alias matches found.")


if __name__ == "__main__":
    main()