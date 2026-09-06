"""Regression tests for tblastn subject-strand detection.

Biopython's NCBIXML reports sbjct_start <= sbjct_end even for minus-strand hits;
the real strand is the sign of the subject frame (hsp.frame[1]). The old logic
inferred strand from coordinate order and so mislabelled every minus-strand hit
as "+", which broke reverse-complement-deposited genomes (e.g. PEDV PX984909.1).
"""
from app.src.lifting.tblastn_lifter import _hsp_strand, _hsp_to_genome_coords


class _FakeHSP:
    def __init__(self, sbjct_start, sbjct_end, frame):
        self.sbjct_start = sbjct_start
        self.sbjct_end = sbjct_end
        self.frame = frame


def test_minus_strand_from_negative_frame_despite_ascending_coords():
    # The exact bug case: coords low->high but subject frame negative.
    h = _FakeHSP(453, 1775, (0, -3))
    assert _hsp_strand(h) == "-"
    assert _hsp_to_genome_coords(h) == (453, 1775, "-")


def test_plus_strand_from_positive_frame():
    h = _FakeHSP(453, 1775, (0, 2))
    assert _hsp_strand(h) == "+"
    assert _hsp_to_genome_coords(h) == (453, 1775, "+")


def test_coords_are_low_to_high_regardless_of_input_order():
    h = _FakeHSP(1775, 453, (0, -1))
    start, end, strand = _hsp_to_genome_coords(h)
    assert (start, end) == (453, 1775)
    assert strand == "-"


def test_fallback_to_coord_order_when_frame_missing():
    # Older/partial HSP objects without a usable frame fall back to coord order.
    assert _hsp_strand(_FakeHSP(453, 1775, None)) == "+"
    assert _hsp_strand(_FakeHSP(1775, 453, None)) == "-"
