from pathlib import Path
import argparse

from src.genbank_parser import (
    load_single_genbank,
    load_genbank_records,
    parse_cds_features,
    parse_mat_peptides,
)
from src.feature_renamer import match_features


def parse_args():
    parser = argparse.ArgumentParser(
        description="Debug feature matching between query and reference."
    )
    parser.add_argument("--reference", required=True, help="Reference GenBank file")
    parser.add_argument("--query", required=True, help="Query GenBank file")
    return parser.parse_args()


def print_match_results(title, matches):
    print(f"\n=== {title} ===")

    if not matches:
        print("(none)")
        return

    for i, item in enumerate(matches, 1):
        q = item["query"]
        r = item["best_ref"]
        d = item["score_details"]

        print(f"\n[{i}] Query: {q['name']}  ({q['start']}-{q['end']})")
        print(f"    Status        : {item['status']}")

        if r is None:
            print("    Best ref      : None")
            continue

        print(f"    Best ref      : {r['name']}  ({r['start']}-{r['end']})")
        print(f"    Total score   : {d['score']}")
        print(f"    Length score  : {d['length_score']}")
        print(f"    Position score: {d['position_score']}")
        print(f"    Order score   : {d['order_score']}")


def main():
    args = parse_args()

    ref_record = load_single_genbank(Path(args.reference))
    query_records = load_genbank_records(Path(args.query))

    ref_cds = parse_cds_features(ref_record)
    ref_mat = parse_mat_peptides(ref_record)

    for record in query_records:
        print("\n" + "=" * 70)
        print(f"RECORD: {record.id}")
        print("=" * 70)

        query_cds = parse_cds_features(record)
        query_mat = parse_mat_peptides(record)

        cds_matches = match_features(query_cds, ref_cds)
        mat_matches = match_features(query_mat, ref_mat)

        print_match_results("CDS MATCHES", cds_matches)
        print_match_results("MAT_PEPTIDE MATCHES", mat_matches)


if __name__ == "__main__":
    main()