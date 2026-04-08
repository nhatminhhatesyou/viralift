from pathlib import Path
import argparse

from src.annotation.gene_alias import (
    load_alias_lookup,
    apply_alias_to_features,
)
from src.io.genbank_parser import load_single_genbank, parse_cds_features


def parse_args():
    parser = argparse.ArgumentParser(
        description="Check alias normalization on GenBank CDS features."
    )
    parser.add_argument(
        "--input",
        required=True,
        help="Input GenBank file (single record). Example: data/PRRS_ref_test.gb",
    )
    parser.add_argument(
        "--alias-config",
        required=True,
        help="Alias JSON config file. Example: config/prrsv_alias.json",
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
        name_source = feature.get("name_source")

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


def main():
    args = parse_args()

    input_path = Path(args.input)
    alias_config_path = Path(args.alias_config)

    record = load_single_genbank(input_path)
    raw_features = parse_cds_features(record)

    if not raw_features:
        raise ValueError("No CDS features found in input file.")

    alias_lookup = load_alias_lookup(alias_config_path)
    normalized_features = apply_alias_to_features(raw_features, alias_lookup)

    summary = summarize_alias_results(normalized_features)

    print("\n=== ALIAS NORMALIZATION SUMMARY ===\n")
    print(f"Input file        : {input_path}")
    print(f"Record ID         : {record.id}")
    print(f"Alias config      : {alias_config_path}")
    print(f"Total CDS         : {summary['total']}")
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

    if summary["renamed_items"]:
        print("\n=== ALIAS-HIT FEATURES ===\n")
        for i, (raw_name, final_name) in enumerate(summary["renamed_items"], start=1):
            print(f"[{i}] {raw_name} -> {final_name}")
    else:
        print("\n=== ALIAS-HIT FEATURES ===\n")
        print("No alias matches found.")


if __name__ == "__main__":
    main()