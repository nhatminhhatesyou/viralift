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
from src.genbank_parser import load_single_genbank, parse_cds_features
from src.minimap_runner import run_minimap2
from src.sam_lifter import get_primary_alignment, build_ref_to_query_map


def parse_args():
    parser = argparse.ArgumentParser(description="Reference-guided viral CDS extraction tool.")
    parser.add_argument("--ref-gb", required=True, help="Reference GenBank file")
    parser.add_argument("--query-gb", required=True, help="Query GenBank file")
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("--min-coverage", type=float, default=0.8, help="Minimum lifted CDS coverage")
    return parser.parse_args()


def write_results_fasta(results, query_id: str, out_path: Path) -> None:
    """Write extracted results to FASTA."""
    records = []

    for item in results:
        if item["status"] != "ok":
            continue
        if not item["sequence"]:
            continue

        record_id = f"{query_id}|{item['name']}|{item['method']}"
        record = SeqRecord(Seq(item["sequence"]), id=record_id, description="")
        records.append(record)

    SeqIO.write(records, str(out_path), "fasta")


def write_results_tsv(results, query_id: str, out_path: Path) -> None:
    """Write extracted results to TSV."""
    rows = []

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

    with open(out_path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main():
    args = parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    ref_record = load_single_genbank(Path(args.ref_gb))
    query_record = load_single_genbank(Path(args.query_gb))

    ref_cds = parse_cds_features(ref_record)
    query_cds = parse_cds_features(query_record)

    strategy = choose_strategy(ref_record, query_record)
    print(f"Selected strategy: {strategy}")

    if strategy == "direct":
        renamed_query_cds = rename_query_cds_by_reference_order(ref_cds, query_cds)
        results = extract_all_direct(query_record, renamed_query_cds)

    else:
        ref_fa = outdir / "tmp_ref.fa"
        query_fa = outdir / "tmp_query.fa"
        sam_path = outdir / "alignment.sam"

        write_record_to_fasta(ref_record, ref_fa)
        write_record_to_fasta(query_record, query_fa)

        run_minimap2(ref_fa, query_fa, sam_path)

        aln = get_primary_alignment(sam_path)
        ref_to_query = build_ref_to_query_map(aln)

        results = extract_all_lifted(
            ref_cds=ref_cds,
            query_record=query_record,
            ref_to_query=ref_to_query,
            min_coverage=args.min_coverage,
        )

    fasta_out = outdir / "extracted_cds.fasta"
    tsv_out = outdir / "extracted_cds.tsv"

    write_results_fasta(results, query_record.id, fasta_out)
    write_results_tsv(results, query_record.id, tsv_out)

    print(f"Done. FASTA: {fasta_out}")
    print(f"Done. TSV: {tsv_out}")


if __name__ == "__main__":
    main()