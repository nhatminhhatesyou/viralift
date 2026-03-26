from pathlib import Path

from src.genbank_parser import load_genbank_records
from src.fasta_writer import write_record_to_fasta
from src.minimap_runner import run_minimap2


if __name__ == "__main__":
    app_dir = Path(__file__).resolve().parents[1]
    data_dir = app_dir / "data"

    # Load one reference record
    ref_records = load_genbank_records(data_dir / "FMD_ref_test.gb")
    ref_record = ref_records[0]

    # Load one query record
    query_records = load_genbank_records(data_dir / "FMD_test.gb")
    query_record = query_records[0]

    # Write temporary FASTA files
    ref_fa = data_dir / "tmp_ref.fa"
    query_fa = data_dir / "tmp_query.fa"
    out_sam = data_dir / "tmp_alignment.sam"

    write_record_to_fasta(ref_record, ref_fa)
    write_record_to_fasta(query_record, query_fa)

    # Run minimap2
    run_minimap2(ref_fa, query_fa, out_sam)

    print("minimap2 finished successfully.")
    print(f"SAM file written to: {out_sam}")