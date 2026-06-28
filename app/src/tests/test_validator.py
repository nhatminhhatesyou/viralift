from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from app.src.lifting.validator import rescue_start_codon, rescue_stop_codon


def test_rescue_stop_codon_prefers_reference_length_when_available():
    record = SeqRecord(Seq("ATGAAACCCTAAAAATGA"), id="query")

    rescued = rescue_stop_codon(
        query_record=record,
        query_start=1,
        query_end=9,
        strand="+",
        max_codons=3,
        expected_length=18,
    )

    assert rescued is not None
    new_start, new_end, sequence, codons_extended = rescued
    assert new_start == 1
    assert new_end == 18
    assert sequence == "ATGAAACCCTAAAAATGA"
    assert codons_extended == 3


def test_rescue_stop_codon_keeps_nearest_stop_without_reference_length():
    record = SeqRecord(Seq("ATGAAACCCTAAAAATGA"), id="query")

    rescued = rescue_stop_codon(
        query_record=record,
        query_start=1,
        query_end=9,
        strand="+",
        max_codons=3,
    )

    assert rescued is not None
    new_start, new_end, sequence, codons_extended = rescued
    assert new_start == 1
    assert new_end == 12
    assert sequence == "ATGAAACCCTAA"
    assert codons_extended == 1


def test_rescue_stop_codon_extends_negative_strand_toward_lower_coordinate():
    # Reverse complement of the full biological CDS ATGAAACCCTAA.
    # The lifted span 4..12 translates to ATGAAACCC and lacks the stop codon,
    # so stop rescue on the negative strand must extend query_start down to 1.
    record = SeqRecord(Seq("TTAGGGTTTCAT"), id="query")

    rescued = rescue_stop_codon(
        query_record=record,
        query_start=4,
        query_end=12,
        strand="-",
        max_codons=3,
    )

    assert rescued is not None
    new_start, new_end, sequence, codons_extended = rescued
    assert new_start == 1
    assert new_end == 12
    assert sequence == "ATGAAACCCTAA"
    assert codons_extended == 1


def test_rescue_start_codon_extends_negative_strand_toward_higher_coordinate():
    # Reverse complement of the full biological CDS ATGAAACCCTAA.
    # The lifted span 1..9 translates to AAACCCTAA and lacks the start codon,
    # so start rescue on the negative strand must extend query_end up to 12.
    record = SeqRecord(Seq("TTAGGGTTTCAT"), id="query")

    rescued = rescue_start_codon(
        query_record=record,
        query_start=1,
        query_end=9,
        strand="-",
        max_window=3,
    )

    assert rescued is not None
    new_start, new_end, sequence, offset = rescued
    assert new_start == 1
    assert new_end == 12
    assert sequence == "ATGAAACCCTAA"
    assert offset == 3
