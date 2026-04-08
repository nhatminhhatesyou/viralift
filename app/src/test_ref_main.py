from pathlib import Path
import argparse

from Bio import SeqIO


def parse_args():
    parser = argparse.ArgumentParser(
        description="Print summary of reference GenBank features (CDS and mat_peptide)."
    )
    parser.add_argument(
        "--reference",
        required=True,
        help="Reference GenBank file (single record)",
    )
    return parser.parse_args()


def load_single_genbank(path: Path):
    records = list(SeqIO.parse(str(path), "genbank"))
    if not records:
        raise ValueError(f"No records found in {path}")
    if len(records) > 1:
        raise ValueError(f"Expected 1 record in reference file, found {len(records)}")
    return records[0]


def get_feature_name(feature):
    qualifiers = feature.qualifiers

    for key in ["gene", "product", "label", "standard_name", "note"]:
        values = qualifiers.get(key)
        if values:
            return values[0]

    return "unknown"


def get_feature_start_end_strand(feature):
    start = int(feature.location.start) + 1
    end = int(feature.location.end)
    strand = feature.location.strand
    return start, end, strand


def summarize_features(record, feature_type):
    items = []

    for feature in record.features:
        if feature.type != feature_type:
            continue

        name = get_feature_name(feature)
        start, end, strand = get_feature_start_end_strand(feature)

        items.append(
            {
                "type": feature.type,
                "name": name,
                "start": start,
                "end": end,
                "strand": strand,
            }
        )

    return items


def print_feature_table(title, items):
    print(f"\n{title}: {len(items)}")

    if not items:
        print("(none)")
        return

    idx_width = max(3, len(str(len(items))))
    name_width = max(10, max(len(str(item["name"])) for item in items))
    start_width = max(7, max(len(str(item["start"])) for item in items))
    end_width = max(5, max(len(str(item["end"])) for item in items))
    strand_width = max(6, max(len(str(item["strand"])) for item in items))

    header = (
        f"{'No.':<{idx_width}} | "
        f"{'Name':<{name_width}} | "
        f"{'Start':>{start_width}} | "
        f"{'End':>{end_width}} | "
        f"{'Strand':>{strand_width}}"
    )

    sep = (
        f"{'-' * idx_width}-+-"
        f"{'-' * name_width}-+-"
        f"{'-' * start_width}-+-"
        f"{'-' * end_width}-+-"
        f"{'-' * strand_width}"
    )

    print(sep)
    print(header)
    print(sep)

    for i, item in enumerate(items, 1):
        print(
            f"{i:<{idx_width}} | "
            f"{item['name']:<{name_width}} | "
            f"{item['start']:>{start_width}} | "
            f"{item['end']:>{end_width}} | "
            f"{item['strand']:>{strand_width}}"
        )

    print(sep)


def main():
    args = parse_args()
    ref_path = Path(args.reference)

    record = load_single_genbank(ref_path)

    cds_list = summarize_features(record, "CDS")
    mat_peptides = summarize_features(record, "mat_peptide")

    print("\n=== REFERENCE SUMMARY ===\n")
    print(f"Record ID         : {record.id}")
    print(f"Description       : {record.description}")
    print(f"Genome length     : {len(record.seq)}")
    print(f"Total CDS         : {len(cds_list)}")
    print(f"Total mat_peptide : {len(mat_peptides)}")

    print_feature_table("CDS FEATURES", cds_list)
    print_feature_table("MAT_PEPTIDE FEATURES", mat_peptides)


if __name__ == "__main__":
    main()