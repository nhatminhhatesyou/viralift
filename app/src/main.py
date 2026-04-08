from pathlib import Path
import argparse
import csv

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm

from src.extractor import extract_all_lifted
from src.fasta_writer import write_record_to_fasta
from src.genbank_parser import (
    load_single_genbank,
    load_genbank_records,
    parse_cds_features,
)
from src.minimap_runner import run_minimap2
from src.sam_lifter import get_primary_alignment, build_ref_to_query_map


def parse_args():
    parser = argparse.ArgumentParser(
        prog="viralift",
        description="Reference-guided viral CDS transfer using minimap2.",
        epilog=(
            "Example:\n"
            "  python -m src.main "
            "--reference data/PRRS_ref_test.gb "
            "--query data/PRRSV_test.gb "
            "--output output/prrsv_multi\n"
        ),
        formatter_class=argparse.RawTextHelpFormatter,
    )

    parser.add_argument(
        "--reference",
        required=True,
        help="Reference GenBank file (single record).",
    )
    parser.add_argument(
        "--query",
        required=True,
        help="Query GenBank file (multi-record).",
    )
    parser.add_argument(
        "--output",
        default="output/run",
        help="Output directory. Default: output/run",
    )
    parser.add_argument(
        "--feature-type",
        default="CDS",
        choices=["CDS"],
        help="Feature type to transfer. Default: CDS",
    )
    parser.add_argument(
        "--min-coverage",
        type=float,
        default=0.8,
        help="Minimum coverage for lifted features. Default: 0.8",
    )
    parser.add_argument(
        "--keep-temp",
        action="store_true",
        help="Keep temporary FASTA and SAM files for debugging.",
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="Reduce console output but still show progress and final summary.",
    )

    return parser.parse_args()


def write_results_fasta(all_results, out_path: Path) -> None:
    """Write extracted results from all query records to one FASTA file."""
    records = []

    for query_id, results in all_results:
        for item in results:
            if item["status"] != "ok":
                continue
            if not item.get("sequence"):
                continue

            record_id = f"{query_id}|{item['name']}|{item['method']}"
            record = SeqRecord(Seq(item["sequence"]), id=record_id, description="")
            records.append(record)

    SeqIO.write(records, str(out_path), "fasta")


def write_results_tsv(all_results, out_path: Path) -> None:
    """Write extracted results from all query records to one TSV file."""
    rows = []

    for query_id, results in all_results:
        for item in results:
            rows.append(
                {
                    "query_id": query_id,
                    "name": item.get("name"),
                    "gene": item.get("gene"),
                    "product": item.get("product"),
                    "start": item.get("start"),
                    "end": item.get("end"),
                    "strand": item.get("strand"),
                    "method": item.get("method"),
                    "status": item.get("status"),
                    "coverage": item.get("coverage"),
                    "length": len(item["sequence"]) if item.get("sequence") else None,
                }
            )

    if not rows:
        return

    with open(out_path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def summarize_counts(all_results):
    """Summarize output status counts."""
    ok = 0
    no_alignment = 0
    low_coverage = 0
    unmapped = 0

    for _, results in all_results:
        for item in results:
            status = item.get("status")
            if status == "ok":
                ok += 1
            elif status == "no_alignment":
                no_alignment += 1
            elif status == "low_coverage":
                low_coverage += 1
            elif status == "unmapped":
                unmapped += 1

    return ok, no_alignment, low_coverage, unmapped


def process_one_query_record(
    ref_record,
    query_record,
    ref_cds,
    outdir: Path,
    min_coverage: float,
    keep_temp: bool = False,
    quiet: bool = False,
):
    """Process one query record using minimap2 only."""
    tmp_dir = outdir / "tmp"
    tmp_dir.mkdir(parents=True, exist_ok=True)

    ref_fa = tmp_dir / f"{query_record.id}_ref.fa"
    query_fa = tmp_dir / f"{query_record.id}_query.fa"
    sam_path = tmp_dir / f"{query_record.id}_alignment.sam"

    write_record_to_fasta(ref_record, ref_fa)
    write_record_to_fasta(query_record, query_fa)

    run_minimap2(ref_fa, query_fa, sam_path, quiet=quiet)

    try:
        aln = get_primary_alignment(sam_path)
        ref_to_query = build_ref_to_query_map(aln)

        results = extract_all_lifted(
            ref_cds=ref_cds,
            query_record=query_record,
            ref_to_query=ref_to_query,
            min_coverage=min_coverage,
        )

    except ValueError:
        results = [
            {
                "name": feature["name"],
                "gene": feature.get("gene"),
                "product": feature.get("product"),
                "start": None,
                "end": None,
                "strand": feature["strand"],
                "sequence": None,
                "method": "minimap_transfer",
                "status": "no_alignment",
                "coverage": 0.0,
            }
            for feature in ref_cds
        ]

    if not keep_temp:
        ref_fa.unlink(missing_ok=True)
        query_fa.unlink(missing_ok=True)
        sam_path.unlink(missing_ok=True)

    return results


def main():
    args = parse_args()

    outdir = Path(args.output)
    outdir.mkdir(parents=True, exist_ok=True)

    ref_record = load_single_genbank(Path(args.reference))
    query_records = load_genbank_records(Path(args.query))

    if not query_records:
        raise ValueError("No query records found.")

    ref_cds = parse_cds_features(ref_record)

    if not ref_cds:
        raise ValueError("Reference record has no CDS features.")

    print("ViraLift")
    print(f"  Reference record : {ref_record.id}")
    print(f"  Feature type     : {args.feature_type}")
    print(f"  Reference CDS    : {len(ref_cds)}")
    print(f"  Query records    : {len(query_records)}")
    print(f"  Output folder    : {outdir}")

    all_results = []
    total = len(query_records)

    iterator = tqdm(
        query_records,
        desc="Processing records",
        unit="record",
        ncols=90,
    )

    for query_record in iterator:
        if not args.quiet:
            iterator.set_postfix_str(query_record.id)

        results = process_one_query_record(
            ref_record=ref_record,
            query_record=query_record,
            ref_cds=ref_cds,
            outdir=outdir,
            min_coverage=args.min_coverage,
            keep_temp=args.keep_temp,
            quiet=args.quiet,
        )
        all_results.append((query_record.id, results))

    fasta_out = outdir / "extracted_cds.fasta"
    tsv_out = outdir / "extracted_cds.tsv"

    write_results_fasta(all_results, fasta_out)
    write_results_tsv(all_results, tsv_out)

    ok, no_alignment, low_coverage, unmapped = summarize_counts(all_results)

    print("\nRun summary")
    print(f"  Query records processed : {total}")
    print(f"  OK                      : {ok}")
    print(f"  No alignment            : {no_alignment}")
    print(f"  Low coverage            : {low_coverage}")
    print(f"  Unmapped                : {unmapped}")

    print("\nOutput files")
    print(f"  FASTA : {fasta_out}")
    print(f"  TSV   : {tsv_out}")


if __name__ == "__main__":
    main()