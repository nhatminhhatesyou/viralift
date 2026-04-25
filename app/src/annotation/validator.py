from typing import Dict, Optional, Tuple

from Bio.SeqRecord import SeqRecord

STOP_CODONS = {"TAA", "TAG", "TGA"}


def validate_cds_boundaries(sequence: Optional[str]) -> Dict:
    """
    Validate that a CDS sequence has a proper start and stop codon.

    Returns a dict with:
        valid           -- True if both start and stop codon are present
        has_start_codon -- True if sequence starts with ATG
        has_stop_codon  -- True if sequence ends with TAA/TAG/TGA
    """
    if not sequence or len(sequence) < 6:
        return {"valid": False, "has_start_codon": False, "has_stop_codon": False}

    seq = sequence.upper()
    has_start = seq[:3] == "ATG"
    has_stop = seq[-3:] in STOP_CODONS

    return {
        "valid": has_start and has_stop,
        "has_start_codon": has_start,
        "has_stop_codon": has_stop,
    }


def rescue_start_codon(
    query_record: SeqRecord,
    query_start: int,
    query_end: int,
    strand: str,
    max_window: int = 50,
) -> Optional[Tuple[int, str, int]]:
    """
    Try to find the nearest ATG to the lifted start position by expanding search.

    Scans positions offset 1, 2, ... max_window in both directions from
    query_start. Returns the closest ATG found.

    Args:
        query_record: Query genome SeqRecord
        query_start:  Lifted start (1-based)
        query_end:    Lifted end (1-based, inclusive)
        strand:       "+" or "-"
        max_window:   Max distance to search (bp)

    Returns:
        (new_start, new_sequence, offset_used) if ATG found, else None
        offset_used is negative = upstream, positive = downstream
    """
    genome_len = len(query_record.seq)

    for offset in range(1, max_window + 1):
        for direction in (-1, +1):  # upstream first (more common for frameshift)
            new_start = query_start + direction * offset

            # Bounds check
            if new_start < 1 or new_start > genome_len:
                continue
            if new_start > query_end:
                continue

            # Extract candidate sequence
            if strand == "+":
                candidate = str(query_record.seq[new_start - 1: query_end]).upper()
            else:
                candidate = str(
                    query_record.seq[query_start - 1: query_end + offset].reverse_complement()
                    if direction == +1
                    else query_record.seq[new_start - 1: query_end].reverse_complement()
                ).upper()

            if candidate[:3] == "ATG":
                return new_start, candidate, direction * offset

    return None


def rescue_stop_codon(
    query_record: SeqRecord,
    query_start: int,
    query_end: int,
    strand: str,
    max_codons: int = 30,
) -> Optional[Tuple[int, str, int]]:
    """
    Scan forward from query_end in-frame to find the nearest stop codon.

    Used when tblastn HSP is truncated before the stop codon — either because
    the C-terminus is divergent (HSP ends early) or the +3 fix wasn't enough.

    Args:
        query_record: Query genome SeqRecord
        query_start:  Lifted start (1-based)
        query_end:    Current end (1-based, inclusive) — expected to lack stop codon
        strand:       "+" or "-"
        max_codons:   Max codons to scan forward (default 30 = 90bp)

    Returns:
        (new_end, new_sequence, codons_extended) if stop found, else None
    """
    genome_len = len(query_record.seq)

    for n in range(1, max_codons + 1):
        extension = n * 3
        new_end = query_end + extension

        if new_end > genome_len:
            break

        if strand == "+":
            candidate = str(query_record.seq[query_start - 1: new_end]).upper()
        else:
            candidate = str(
                query_record.seq[query_start - 1: new_end].reverse_complement()
            ).upper()

        if candidate[-3:] in STOP_CODONS:
            return new_end, candidate, n

    return None
