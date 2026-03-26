from typing import Dict, List, Optional, Tuple

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def lift_feature_coordinates(feature: Dict, ref_to_query: Dict[int, int]) -> Tuple[Optional[int], Optional[int], float]:
    """
    Lift one feature from reference coordinates to query coordinates.
    Returns query_start, query_end, coverage.
    """
    mapped_positions = []

    for ref_pos in range(feature["start"], feature["end"] + 1):
        if ref_pos in ref_to_query:
            mapped_positions.append(ref_to_query[ref_pos])

    feature_length = feature["end"] - feature["start"] + 1
    coverage = len(mapped_positions) / feature_length if feature_length > 0 else 0.0

    if not mapped_positions:
        return None, None, coverage

    return min(mapped_positions), max(mapped_positions), coverage


def extract_lifted_sequence(query_record: SeqRecord, query_start: int, query_end: int, strand: str) -> Seq:
    """Extract a lifted feature sequence from the query record."""
    seq = query_record.seq[query_start - 1: query_end]

    if strand == "-":
        seq = seq.reverse_complement()

    return seq


def extract_all_lifted(
    ref_cds: List[Dict],
    query_record: SeqRecord,
    ref_to_query: Dict[int, int],
    min_coverage: float = 0.8,
) -> List[Dict]:
    """Extract all lifted CDS features from query using reference coordinates."""
    results = []

    for feature in ref_cds:
        query_start, query_end, coverage = lift_feature_coordinates(feature, ref_to_query)

        if query_start is None or query_end is None:
            results.append(
                {
                    "name": feature["name"],
                    "gene": feature.get("gene"),
                    "product": feature.get("product"),
                    "start": None,
                    "end": None,
                    "strand": feature["strand"],
                    "sequence": None,
                    "method": "minimap_transfer",
                    "status": "unmapped",
                    "coverage": round(coverage, 4),
                }
            )
            continue

        if coverage < min_coverage:
            results.append(
                {
                    "name": feature["name"],
                    "gene": feature.get("gene"),
                    "product": feature.get("product"),
                    "start": query_start,
                    "end": query_end,
                    "strand": feature["strand"],
                    "sequence": None,
                    "method": "minimap_transfer",
                    "status": "low_coverage",
                    "coverage": round(coverage, 4),
                }
            )
            continue

        seq = extract_lifted_sequence(query_record, query_start, query_end, feature["strand"])

        results.append(
            {
                "name": feature["name"],
                "gene": feature.get("gene"),
                "product": feature.get("product"),
                "start": query_start,
                "end": query_end,
                "strand": feature["strand"],
                "sequence": str(seq),
                "method": "minimap_transfer",
                "status": "ok",
                "coverage": round(coverage, 4),
            }
        )

    return results