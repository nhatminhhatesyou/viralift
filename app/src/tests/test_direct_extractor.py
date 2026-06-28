from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from app.src.features.direct_extractor import direct_extract_with_alias


def _cds(start_zero_based: int, end_one_based: int, gene: str) -> SeqFeature:
    return SeqFeature(
        FeatureLocation(start_zero_based, end_one_based, strand=1),
        type="CDS",
        qualifiers={"gene": [gene]},
    )


def test_direct_extract_collapses_duplicate_canonical_features_without_ref_match():
    record = SeqRecord(Seq("A" * 300), id="query")
    record.features = [
        _cds(9, 210, "ORF1ab"),
        _cds(9, 120, "ORF1ab"),
        _cds(219, 260, "S"),
    ]

    ref_features = [{"name": "S", "start": 220, "end": 260, "length": 41}]
    alias_lookup = {"orf1ab": "ORF1ab", "s": "S"}

    results = direct_extract_with_alias(
        query_record=record,
        query_feature_type="CDS",
        ref_features=ref_features,
        alias_lookup=alias_lookup,
    )

    names = [feature.name for feature in results]
    assert names == ["ORF1ab", "S"]

    orf1ab = results[0]
    assert orf1ab.query_start == 10
    assert orf1ab.query_end == 210
    assert orf1ab.status == "not_in_reference"


def test_direct_extract_reports_alias_source_value_instead_of_generic_display_name():
    record = SeqRecord(Seq("A" * 900), id="query")
    record.features = [
        SeqFeature(
            FeatureLocation(99, 300, strand=1),
            type="CDS",
            qualifiers={"product": ["unknown"], "note": ["ORF3"]},
        )
    ]

    ref_features = [{"name": "ORF3", "start": 100, "end": 300, "length": 201}]
    alias_lookup = {
        "unknown": "__ignored__",
        "orf3": "ORF3",
    }

    results = direct_extract_with_alias(
        query_record=record,
        query_feature_type="CDS",
        ref_features=ref_features,
        alias_lookup=alias_lookup,
    )

    assert len(results) == 1
    assert results[0].name == "ORF3"
    assert results[0].source_name == "ORF3"
    assert results[0].status == "ok"
