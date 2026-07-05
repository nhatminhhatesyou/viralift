from app.src.alias.alias_bootstrap import deduplicate_suggestions


def _suggestion(raw_value: str, canonical_name: str, score: int = 9):
    return {
        "record_id": f"{canonical_name}_record",
        "raw_value": raw_value,
        "field": "gene",
        "canonical_name": canonical_name,
        "suggested_action": "save_alias",
        "confidence": "high",
        "score": score,
        "reason": "coordinate IoU 0.96; strong field: gene",
        "default_save": True,
    }


def test_deduplicate_demotes_raw_name_seen_in_multiple_canonicals():
    rows = deduplicate_suggestions([
        _suggestion("HNZK1", "ORF3"),
        _suggestion("HNZK1", "S"),
    ])

    assert len(rows) == 2
    assert {row["canonical_name"] for row in rows} == {"ORF3", "S"}
    for row in rows:
        assert row["suggested_action"] == "manual_review"
        assert row["confidence"] == "medium"
        assert row["default_save"] is False
        assert "multiple canonical targets" in row["reason"]


def test_deduplicate_keeps_single_target_alias_saveable():
    rows = deduplicate_suggestions([
        _suggestion("S", "S"),
        _suggestion("S", "S"),
    ])

    assert len(rows) == 1
    assert rows[0]["canonical_name"] == "S"
    assert rows[0]["suggested_action"] == "save_alias"
    assert rows[0]["confidence"] == "high"
    assert rows[0]["support_count"] == 1
