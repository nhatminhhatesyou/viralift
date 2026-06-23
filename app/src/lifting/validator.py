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
        in_frame        -- True if sequence length is divisible by 3
    """
    if not sequence or len(sequence) < 6:
        return {
            "valid": False,
            "has_start_codon": False,
            "has_stop_codon": False,
            "in_frame": False,
        }

    seq = sequence.upper()
    has_start = seq[:3] == "ATG"
    has_stop = seq[-3:] in STOP_CODONS
    in_frame = len(seq) % 3 == 0

    return {
        "valid": has_start and has_stop and in_frame,
        "has_start_codon": has_start,
        "has_stop_codon": has_stop,
        "in_frame": in_frame,
    }


def rescue_start_codon(
    query_record: SeqRecord,
    query_start: int,
    query_end: int,
    strand: str,
    max_window: int = 50,
    expected_length: Optional[int] = None,
) -> Optional[Tuple[int, str, int]]:
    """
    Try to find a plausible ATG around the lifted start position.

    Earlier versions returned the nearest ATG. That can pick an internal ATG
    when the HSP starts late, especially for short CDS features. This function
    now scores all ATG candidates in the rescue window and prefers candidates
    that preserve CDS frame and, when available, reference CDS length.

    Args:
        query_record: Query genome SeqRecord
        query_start:  Lifted start (1-based)
        query_end:    Lifted end (1-based, inclusive)
        strand:       "+" or "-"
        max_window:   Max distance to search (bp)
        expected_length: Expected CDS length in bp, including stop codon.

    Returns:
        (new_start, new_sequence, offset_used) if ATG found, else None
        offset_used is negative = upstream, positive = downstream
    """
    genome_len = len(query_record.seq)
    candidates = []

    def add_candidate(new_start: int, offset_used: int) -> None:
        if new_start < 1 or new_start > genome_len:
            return
        if new_start > query_end:
            return
        if strand == "+":
            candidate = str(query_record.seq[new_start - 1: query_end]).upper()
        else:
            candidate = str(
                query_record.seq[new_start - 1: query_end].reverse_complement()
            ).upper()
        if candidate[:3] != "ATG":
            return
        length = len(candidate)
        length_delta = (
            abs(length - expected_length)
            if expected_length is not None
            else 0
        )
        candidates.append((
            length % 3 != 0,
            length_delta,
            abs(offset_used),
            0 if offset_used < 0 else 1,
            new_start,
            candidate,
            offset_used,
        ))

    if expected_length is not None and strand == "+":
        expected_start = query_end - expected_length + 1
        add_candidate(expected_start, expected_start - query_start)

    for offset in range(1, max_window + 1):
        for direction in (-1, +1):  # upstream first (more common for frameshift)
            new_start = query_start + direction * offset

            # Bounds check
            if new_start < 1 or new_start > genome_len:
                continue
            if new_start > query_end:
                continue

            add_candidate(new_start, direction * offset)

    if not candidates:
        return None

    _, _, _, _, new_start, candidate, offset = min(candidates)
    return new_start, candidate, offset


def rescue_stop_codon(
    query_record: SeqRecord,
    query_start: int,
    query_end: int,
    strand: str,
    max_codons: int = 30,
    expected_length: Optional[int] = None,
) -> Optional[Tuple[int, str, int]]:
    """
    Scan forward from query_end in-frame to find a plausible stop codon.

    Used when tblastn HSP is truncated before the stop codon — either because
    the C-terminus is divergent (HSP ends early) or the +3 fix wasn't enough.
    Earlier versions returned the first stop codon found. That can over-trust
    a premature nearby stop. This function now scores all stop candidates in
    the scan window and, when reference length is available, prefers the stop
    whose resulting CDS length is closest to the reference CDS length.

    Args:
        query_record: Query genome SeqRecord
        query_start:  Lifted start (1-based)
        query_end:    Current end (1-based, inclusive) — expected to lack stop codon
        strand:       "+" or "-"
        max_codons:   Max codons to scan forward (default 30 = 90bp)
        expected_length: Expected CDS length in bp, including stop codon.

    Returns:
        (new_end, new_sequence, codons_extended) if stop found, else None
    """
    genome_len = len(query_record.seq)
    candidates = []

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
            length = len(candidate)
            length_delta = (
                abs(length - expected_length)
                if expected_length is not None
                else 0
            )
            candidates.append((
                length % 3 != 0,
                length_delta,
                n,
                new_end,
                candidate,
            ))

    if not candidates:
        return None

    _, _, codons_extended, new_end, candidate = min(candidates)
    return new_end, candidate, codons_extended
