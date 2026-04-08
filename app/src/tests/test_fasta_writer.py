from pathlib import Path

from app.src.io.genbank_parser import load_genbank_records
from app.src.io.fasta_writer import write_record_to_fasta


if __name__ == "__main__":
    app_dir = Path(__file__).resolve().parents[1]
    data_dir = app_dir / "data"
    out_dir = app_dir / "data"

    # Test with the first FMDV record
    fmd_records = load_genbank_records(data_dir / "FMD_test.gb")
    first_fmd = fmd_records[0]
    write_record_to_fasta(first_fmd, out_dir / "FMD_test_record1.fa")

    # Test with the first PRRSV record
    prrsv_records = load_genbank_records(data_dir / "PRRSV_test.gb")
    first_prrsv = prrsv_records[0]
    write_record_to_fasta(first_prrsv, out_dir / "PRRSV_test_record1.fa")

    print("FASTA files written successfully.")