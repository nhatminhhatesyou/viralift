from pathlib import Path
import argparse

from Bio import SeqIO


def parse_args():
    parser = argparse.ArgumentParser(
        description="Check input GenBank records for CDS and mat_peptide annotations."
    )
    parser.add_argument(
        "--input",
        required=True,
        help="Input GenBank file (single or multi-record)",
    )
    return parser.parse_args()


def load_genbank_records(path: Path):
    return list(SeqIO.parse(str(path), "genbank"))


def get_feature_name(feature):
    qualifiers = feature.qualifiers

    for key in ["gene", "product", "label", "standard_name", "note"]:
        values = qualifiers.get(key)
        if values:
            return values[0]

    return "unknown"


def summarize_record(record):
    cds_features = []
    mat_peptides = []

    for feature in record.features:
        if feature.type == "CDS":
            cds_features.append(feature)
        elif feature.type == "mat_peptide":
            mat_peptides.append(feature)

    mat_names = [get_feature_name(f) for f in mat_peptides]

    return {
        "record_id": record.id,
        "has_cds": len(cds_features) > 0,
        "cds_count": len(cds_features),
        "has_mat_peptide": len(mat_peptides) > 0,
        "mat_peptide_count": len(mat_peptides),
        "mat_peptide_names": mat_names,
    }


def main():
    args = parse_args()
    records = load_genbank_records(Path(args.input))

    if not records:
        raise ValueError("No records found in input file.")

    print("\n=== INPUT ANNOTATION SUMMARY ===\n")

    full_annot_count = 0
    cds_only_count = 0
    no_cds_count = 0

    for record in records:
        info = summarize_record(record)

        print(f"Record: {info['record_id']}")
        print(f"  CDS count           : {info['cds_count']}")
        print(f"  Has mat_peptide     : {'yes' if info['has_mat_peptide'] else 'no'}")

        if info["has_cds"] and info["has_mat_peptide"]:
            full_annot_count += 1
            print("  Annotation status   : FULL (CDS has been split into mat_peptide)")
            print(f"  mat_peptide count   : {info['mat_peptide_count']}")
            print("  mat_peptide names   :")
            for i, name in enumerate(info["mat_peptide_names"], 1):
                print(f"    [{i}] {name}")

        elif info["has_cds"] and not info["has_mat_peptide"]:
            cds_only_count += 1
            print("  Annotation status   : CDS only (no mat_peptide annotation)")

        else:
            no_cds_count += 1
            print("  Annotation status   : No CDS")

        print()

    print("=== SUMMARY ===")
    print(f"Total records                    : {len(records)}")
    print(f"Full annotation (CDS + mat_peptide): {full_annot_count}")
    print(f"CDS only                         : {cds_only_count}")
    print(f"No CDS                           : {no_cds_count}")


if __name__ == "__main__":
    main()