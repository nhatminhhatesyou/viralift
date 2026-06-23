from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from app.src.lifting.validator import rescue_stop_codon


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
    new_end, sequence, codons_extended = rescued
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
    new_end, sequence, codons_extended = rescued
    assert new_end == 12
    assert sequence == "ATGAAACCCTAA"
    assert codons_extended == 1
