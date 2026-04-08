from Bio.SeqRecord import SeqRecord
from app.src.io.genbank_parser import parse_cds_features


"""
Module: annotation_strategy.py

Purpose:
    Decide which extraction strategy to use for a query genome:
    - "direct": when query and reference have matching CDS structure
    - "minimap": fallback using alignment-based coordinate transfer

Notes:
    - This module only decides the strategy, it does NOT perform extraction.
    - Alias-based naming and feature normalization should be handled separately.
"""


def has_matching_cds_structure(ref_record: SeqRecord, query_record: SeqRecord) -> bool:
    """
    Check whether the query genome has the same CDS structure as the reference.

    Criteria:
        - Same number of CDS features

    Args:
        ref_record: Reference genome (GenBank record)
        query_record: Query genome (GenBank record)

    Returns:
        True if structures match, False otherwise
    """
    ref_cds = parse_cds_features(ref_record)
    query_cds = parse_cds_features(query_record)

    # If reference has no CDS, cannot compare structure
    if not ref_cds:
        return False

    return len(ref_cds) == len(query_cds)


def choose_strategy(ref_record: SeqRecord, query_record: SeqRecord) -> str:
    """
    Select extraction strategy based on CDS structure similarity.

    Strategy:
        - "direct": if CDS counts match (assumes similar genome structure)
        - "minimap": otherwise, use alignment-based transfer

    Args:
        ref_record: Reference genome
        query_record: Query genome

    Returns:
        Strategy name ("direct" or "minimap")
    """
    if has_matching_cds_structure(ref_record, query_record):
        return "direct"

    return "minimap"