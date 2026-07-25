import json

from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from app.src.alias.gene_alias import normalize_text
from ui.services import _save_to_alias_config, _scan_unknown_names


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
        excluded_names=set(),
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
        excluded_names=set(),
    )

    assert unknown == {
        "unglycosylated membrane protein": {
            "records": ["query1"],
            "candidates": ["unglycosylated membrane protein", "ORF6; M"],
            "ambiguous": False,
        }
    }


def test_save_to_alias_config_uses_backup_aware_save(tmp_path):
    config_path = tmp_path / "virus_alias.json"
    config_path.write_text(
        json.dumps({"virus": "Test", "canonical_names": {"ORF6": []}}),
        encoding="utf-8",
    )

    written = _save_to_alias_config(
        config_path,
        {"unglycosylated membrane protein": "ORF6"},
    )

    assert written == 1
    saved = json.loads(config_path.read_text(encoding="utf-8"))
    assert saved["canonical_names"]["ORF6"] == ["unglycosylated membrane protein"]
    backups = list((tmp_path / "backups").glob("virus_alias.*.json"))
    assert len(backups) == 1
