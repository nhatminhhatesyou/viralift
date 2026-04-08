from Bio.SeqRecord import SeqRecord
from src.genbank_parser import parse_cds_features


def has_matching_cds_structure(ref_record: SeqRecord, query_record: SeqRecord) -> bool:
    """Check if query has the same number of CDS features as reference."""
    ref_cds = parse_cds_features(ref_record)
    query_cds = parse_cds_features(query_record)

    if len(ref_cds) == 0:
        return False

    return len(ref_cds) == len(query_cds)


def choose_strategy(ref_record: SeqRecord, query_record: SeqRecord) -> str:
    """Choose extraction strategy."""
    if has_matching_cds_structure(ref_record, query_record):
        return "direct"

    return "minimap"