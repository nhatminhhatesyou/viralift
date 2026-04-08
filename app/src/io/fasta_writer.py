from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def write_record_to_fasta(record: SeqRecord, out_path: Path) -> None:
    """Write one SeqRecord to a FASTA file."""
    SeqIO.write(record, str(out_path), "fasta")