from pathlib import Path
import argparse
import csv

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

from src.annotation_strategy import choose_strategy
from src.direct_extractor import extract_all_direct
from src.extractor import extract_all_lifted
from src.fasta_writer import write_record_to_fasta
from src.feature_renamer import rename_query_cds_by_reference_order
from src.genbank_parser import load_single_genbank, load_genbank_records, parse_cds_features
from src.minimap_runner import run_minimap2
from src.sam_lifter import get_primary_alignment, build_ref_to_query_map


def parse_args():
    parser = argparse.ArgumentParser(description="Reference-guided viral CDS extraction tool.")
    parser.add_argument("--ref-gb", required=True, help="Reference GenBank file (single record)")
    parser.add_argument("--query-gb", required=True, help="Query GenBank file (multi-record)")
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("--min-coverage", type=float, default=0.8, help="Minimum lifted CDS coverage")
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


def process_one_query_record(ref_record, query_record, outdir: Path, min_coverage: float):
    """Process one query record using either direct or minimap strategy."""
    ref_cds = parse_cds_features(ref_record)
    query_cds = parse_cds_features(query_record)

    strategy = choose_strategy(ref_record, query_record)
    print(f"[{query_record.id}] strategy = {strategy}")

    if strategy == "direct":
        renamed_query_cds = rename_query_cds_by_reference_order(ref_cds, query_cds)
        return extract_all_direct(query_record, renamed_query_cds)

    # minimap fallback
    tmp_dir = outdir / "tmp"
    tmp_dir.mkdir(parents=True, exist_ok=True)

    ref_fa = tmp_dir / f"{query_record.id}_ref.fa"
    query_fa = tmp_dir / f"{query_record.id}_query.fa"
    sam_path = tmp_dir / f"{query_record.id}_alignment.sam"

    write_record_to_fasta(ref_record, ref_fa)
    write_record_to_fasta(query_record, query_fa)

    run_minimap2(ref_fa, query_fa, sam_path)

    try:
        aln = get_primary_alignment(sam_path)
        ref_to_query = build_ref_to_query_map(aln)

        return extract_all_lifted(
            ref_cds=ref_cds,
            query_record=query_record,
            ref_to_query=ref_to_query,
            min_coverage=min_coverage,
        )
    except ValueError:
        return [
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

def main():
    args = parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    ref_record = load_single_genbank(Path(args.ref_gb))
    query_records = load_genbank_records(Path(args.query_gb))

    if not query_records:
        raise ValueError("No query records found.")

    print(f"Using reference record: {ref_record.id}")
    print(f"Total query records: {len(query_records)}")

    all_results = []

    for query_record in query_records:
        results = process_one_query_record(
            ref_record=ref_record,
            query_record=query_record,
            outdir=outdir,
            min_coverage=args.min_coverage,
        )
        all_results.append((query_record.id, results))

    fasta_out = outdir / "extracted_cds.fasta"
    tsv_out = outdir / "extracted_cds.tsv"

    write_results_fasta(all_results, fasta_out)
    write_results_tsv(all_results, tsv_out)

    print(f"Done. FASTA: {fasta_out}")
    print(f"Done. TSV: {tsv_out}")


if __name__ == "__main__":
    main()