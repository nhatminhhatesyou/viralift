from pathlib import Path
from typing import Dict

import pysam


def get_primary_alignment(sam_path: Path) -> pysam.AlignedSegment:
    """Return the first primary alignment."""
    with pysam.AlignmentFile(str(sam_path), "r") as sam_file:
        for aln in sam_file:
            if aln.is_unmapped:
                continue
            if aln.is_secondary or aln.is_supplementary:
                continue
            return aln

    raise ValueError("No primary alignment found.")


def build_ref_to_query_map(aln: pysam.AlignedSegment) -> Dict[int, int]:
    """Build mapping from 1-based reference positions to 1-based query positions."""
    ref_to_query = {}

    for query_pos, ref_pos in aln.get_aligned_pairs(matches_only=False):
        if query_pos is None or ref_pos is None:
            continue

        ref_to_query[ref_pos + 1] = query_pos + 1

    return ref_to_query