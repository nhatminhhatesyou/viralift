from typing import Dict, List, Optional, Tuple

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from app.src.annotation.validator import validate_cds_boundaries, rescue_start_codon


"""
Module: extractor.py

Purpose:
    Lift CDS (coding sequence) features from a reference genome
    onto a query genome using a position mapping (ref_to_query),
    and extract corresponding nucleotide sequences.

Notes:
    - This module handles coordinate transformation and sequence extraction only.
    - It does NOT perform gene name normalization or alias resolution.
    - Naming consistency should be handled upstream (e.g., alias mapping layer).
"""


def lift_feature_coordinates(
    feature: Dict,
    ref_to_query: Dict[int, int]
) -> Tuple[Optional[int], Optional[int], float]:
    """
    Map a feature from reference coordinates to query coordinates.

    Args:
        feature: Dictionary containing feature metadata with "start" and "end".
        ref_to_query: Mapping from reference positions to query positions.

    Returns:
        query_start: Minimum mapped position on query (None if unmapped).
        query_end: Maximum mapped position on query (None if unmapped).
        coverage: Fraction of feature positions successfully mapped.
    """
    mapped_positions = []

    # Iterate through each nucleotide position in the reference feature
    for ref_pos in range(feature["start"], feature["end"] + 1):
        if ref_pos in ref_to_query:
            mapped_positions.append(ref_to_query[ref_pos])

    feature_length = feature["end"] - feature["start"] + 1
    coverage = len(mapped_positions) / feature_length if feature_length > 0 else 0.0

    if not mapped_positions:
        return None, None, coverage

    return min(mapped_positions), max(mapped_positions), coverage


def extract_lifted_sequence(
    query_record: SeqRecord,
    query_start: int,
    query_end: int,
    strand: str
) -> Seq:
    """
    Extract a sequence from the query genome based on lifted coordinates.

    Args:
        query_record: Biopython SeqRecord containing query genome.
        query_start: Start position (1-based).
        query_end: End position (1-based, inclusive).
        strand: "+" or "-" indicating strand orientation.

    Returns:
        Extracted sequence (reverse complemented if strand is "-").
    """
    # Convert 1-based inclusive coordinates to Python slicing (0-based, end-exclusive)
    seq = query_record.seq[query_start - 1: query_end]

    if strand == "-":
        seq = seq.reverse_complement()

    return seq


def extract_all_lifted(
    ref_cds: List[Dict],
    query_record: SeqRecord,
    ref_to_query: Dict[int, int],
    min_coverage: float = 0.8,
    validate_codons: bool = True,
    rescue_window: int = 50,
) -> List[Dict]:
    """
    Lift all CDS features from reference to query genome and extract sequences.

    Workflow:
        1. Lift coordinates using ref_to_query mapping
        2. Evaluate coverage
        3. Extract sequence if valid
        4. Return structured result per feature

    Args:
        ref_cds: List of reference CDS features.
        query_record: Query genome sequence.
        ref_to_query: Mapping from reference positions to query positions.
        min_coverage: Minimum coverage threshold to accept a lifted feature.

    Returns:
        List of dictionaries with lifted feature information.
    """
    results = []

    for feature in ref_cds:
        query_start, query_end, coverage = lift_feature_coordinates(feature, ref_to_query)

        base_result = {
            "name": feature["name"],
            "canonical_name": feature.get("canonical_name"),  # set by apply_ref_naming
            "gene": feature.get("gene"),
            "product": feature.get("product"),
            "strand": feature["strand"],
            "method": "minimap_transfer",
            "coverage": round(coverage, 4),

            # Debug fields (useful for tracing issues)
            "ref_start": feature["start"],
            "ref_end": feature["end"],
        }

        # Case 1: completely unmapped
        if query_start is None or query_end is None:
            results.append({
                **base_result,
                "start": None,
                "end": None,
                "sequence": None,
                "status": "unmapped",
            })
            continue

        # Case 2: low coverage
        if coverage < min_coverage:
            results.append({
                **base_result,
                "start": query_start,
                "end": query_end,
                "sequence": None,
                "status": "low_coverage",
            })
            continue

        # Case 3: valid feature → extract sequence + optionally validate boundaries
        seq = extract_lifted_sequence(
            query_record,
            query_start,
            query_end,
            feature["strand"]
        )
        seq_str = str(seq)

        if validate_codons:
            validation = validate_cds_boundaries(seq_str)

            if validation["valid"]:
                status = "ok"
                extra = {
                    "has_start_codon": True,
                    "has_stop_codon": True,
                    "rescue_offset": None,
                }
            elif not validation["has_start_codon"]:
                # Try to find nearest ATG within expanding window
                rescued = rescue_start_codon(
                    query_record, query_start, query_end,
                    feature["strand"], max_window=rescue_window,
                )
                if rescued:
                    new_start, seq_str, offset = rescued
                    query_start = new_start
                    revalidation = validate_cds_boundaries(seq_str)
                    status = "ok_rescued" if revalidation["has_stop_codon"] else "invalid_boundaries"
                    extra = {
                        "has_start_codon": True,
                        "has_stop_codon": revalidation["has_stop_codon"],
                        "rescue_offset": offset,
                    }
                else:
                    status = "invalid_boundaries"
                    extra = {
                        "has_start_codon": False,
                        "has_stop_codon": validation["has_stop_codon"],
                        "rescue_offset": None,
                    }
            else:
                # Has start but no stop — rescue not attempted (different problem)
                status = "invalid_boundaries"
                extra = {
                    "has_start_codon": True,
                    "has_stop_codon": False,
                    "rescue_offset": None,
                }
        else:
            status = "ok"
            extra = {}

        results.append({
            **base_result,
            "start": query_start,
            "end": query_end,
            "sequence": seq_str,
            "status": status,
            **extra,
        })

    return results