from typing import Dict, List

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


"""
Module: direct_extractor.py

Purpose:
    Extract feature sequences directly from query annotation coordinates,
    without alignment-based coordinate lifting.

Use case:
    This module is intended for cases where query features are already
    annotated and their coordinates can be used directly.

Notes:
    - This module assumes feature coordinates are already on the query genome.
    - It does NOT perform alias resolution or feature matching.
    - Upstream code may optionally provide:
        - raw_name
        - ref_name
        - ref_gene
        - ref_product
"""


def extract_feature_sequence(record: SeqRecord, feature: Dict) -> Seq:
    """
    Extract a single feature sequence from a query record.

    Args:
        record: Query genome record
        feature: Feature dictionary containing at least:
            - start (1-based)
            - end (1-based, inclusive)
            - strand ("+" or "-")

    Returns:
        Extracted nucleotide sequence.
    """
    seq = record.seq[feature["start"] - 1: feature["end"]]

    if feature["strand"] == "-":
        seq = seq.reverse_complement()

    return seq


def extract_all_direct(record: SeqRecord, features: List[Dict]) -> List[Dict]:
    """
    Extract all feature sequences directly from query coordinates.

    Args:
        record: Query genome record
        features: List of feature dictionaries already defined on the query genome

    Returns:
        List of extracted feature result dictionaries
    """
    results: List[Dict] = []

    for feature in features:
        seq = extract_feature_sequence(record, feature)

        results.append(
            {
                "raw_name": feature.get("raw_name", feature.get("name")),
                "name": feature["name"],
                "ref_name": feature.get("ref_name", feature["name"]),
                "gene": feature.get("ref_gene", feature.get("gene")),
                "product": feature.get("ref_product", feature.get("product")),
                "start": feature["start"],
                "end": feature["end"],
                "strand": feature["strand"],
                "sequence": str(seq),
                "method": "direct_annotation",
                "status": "ok",
            }
        )

    return results