from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from app.src.alias.gene_alias import normalize_text
from ui.services import _scan_unknown_names


def _record_with_cds(*, product: str, note: str) -> SeqRecord:
    record = SeqRecord(Seq("ATG" + "AAA" * 20 + "TAA"), id="query1")
    record.features.append(
        SeqFeature(
            FeatureLocation(0, len(record.seq), strand=1),
            type="CDS",
            qualifiers={
                "product": [product],
                "note": [note],
            },
        )
    )
    return record


def test_scan_unknown_names_skips_feature_resolved_by_semicolon_candidate():
    record = _record_with_cds(
        product="unglycosylated membrane protein",
        note="ORF6; M",
    )
    alias_lookup = {
        normalize_text("ORF6"): "ORF6",
        normalize_text("M"): "ORF6",
    }

    unknown = _scan_unknown_names(
        query_records=[record],
        alias_lookup=alias_lookup,
        ignored_names=set(),
    )

    assert unknown == {}


def test_scan_unknown_names_keeps_feature_when_semicolon_candidate_is_unresolved():
    record = _record_with_cds(
        product="unglycosylated membrane protein",
        note="ORF6; M",
    )

    unknown = _scan_unknown_names(
        query_records=[record],
        alias_lookup={},
        ignored_names=set(),
    )

    assert unknown == {
        "unglycosylated membrane protein": {
            "records": ["query1"],
            "candidates": ["unglycosylated membrane protein", "ORF6; M"],
            "ambiguous": False,
        }
    }
