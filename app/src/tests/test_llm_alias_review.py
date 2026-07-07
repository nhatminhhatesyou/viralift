from app.src.llm.alias_review import (
    needs_llm_review,
    review_unresolved_names,
    review_uncertain_alias_suggestions,
)
from app.src.alias.alias_payload import build_uncertain_suggestion_review_payload
from app.src.alias.alias_payload import build_unresolved_name_review_payload
from app.src.llm.config import LLMConfig
from app.src.llm.provider import MockLLMProvider


def _row(raw_value, action, confidence, canonical="M", field="product"):
    return {
        "record_id": "Q1",
        "raw_value": raw_value,
        "field": field,
        "canonical_name": canonical,
        "suggested_action": action,
        "confidence": confidence,
        "reason": "test reason",
        "score": 3,
        "default_save": action == "save_alias",
        "support_count": 1,
        "support_records": "Q1",
        "iou": 0.96,
        "coverage": 1.0,
        "identity": 99.0,
    }


def test_needs_llm_review_only_for_uncertain_or_risky_rows():
    assert not needs_llm_review(_row("M", "save_alias", "high", field="gene"))
    assert needs_llm_review(_row("membrane protein", "manual_review", "medium"))
    assert needs_llm_review(_row("ORF6; M", "save_alias", "high"))
    assert not needs_llm_review(_row("ORF1a protein", "save_alias", "high", canonical="ORF1a"))
    assert not needs_llm_review(_row("ORF1ab", "save_alias", "high", canonical="ORF1ab"))
    assert not needs_llm_review(_row("hypothetical protein", "ignore", "low"))
    assert needs_llm_review(_row("polyprotein", "ignore", "low", canonical="ORF1a"))
    assert needs_llm_review(_row("glycoprotein", "ignore", "low", canonical="S"))
    assert needs_llm_review(_row("structural protein", "ignore", "low", canonical="S"))
    assert needs_llm_review(_row("ORF2a", "ignore", "low", canonical="ORF2"))


def test_uncertain_payload_marks_matching_available_canonical():
    payload = build_uncertain_suggestion_review_payload(
        virus_name="test virus",
        canonical_names=["ORF1a", "ORF1ab", "N"],
        suggestions=[
            _row("ORF1a protein", "manual_review", "medium", canonical="ORF1a"),
            _row("ORF1ab", "manual_review", "medium", canonical="ORF1a"),
        ],
    )

    assert payload["available_canonicals"] == ["ORF1a", "ORF1ab", "N"]
    assert payload["suggestions"][0]["matching_available_canonical"] == "ORF1a"
    assert payload["suggestions"][1]["matching_available_canonical"] == "ORF1ab"


def test_uncertain_payload_marks_specific_polyprotein_as_orf_match():
    payload = build_uncertain_suggestion_review_payload(
        virus_name="test virus",
        canonical_names=["ORF1a", "ORF1b", "ORF1ab"],
        suggestions=[
            _row("polyprotein 1a", "ignore", "low", canonical="ORF1a"),
            _row("replicase polyprotein 1b", "ignore", "low", canonical="ORF1b"),
            _row("polyprotein 1ab", "ignore", "low", canonical="ORF1ab"),
        ],
    )

    assert payload["suggestions"][0]["matching_available_canonical"] == "ORF1a"
    assert payload["suggestions"][1]["matching_available_canonical"] == "ORF1b"
    assert payload["suggestions"][2]["matching_available_canonical"] == "ORF1ab"


def test_uncertain_payload_includes_cross_canonical_context():
    row = _row("small membrane protein", "manual_review", "medium", canonical="E")
    row["shared_across_canonicals"] = True
    row["cross_canonical_targets"] = ["E", "M"]
    row["cross_canonical_target_count"] = 2

    payload = build_uncertain_suggestion_review_payload(
        virus_name="test virus",
        canonical_names=["E", "M"],
        suggestions=[row],
    )

    suggestion = payload["suggestions"][0]
    assert suggestion["shared_across_canonicals"] is True
    assert suggestion["cross_canonical_targets"] == ["E", "M"]
    assert suggestion["cross_canonical_target_count"] == 2


def test_unresolved_payload_marks_combined_orf_when_orf1ab_available():
    payload = build_unresolved_name_review_payload(
        virus_name="test virus",
        canonical_names=["ORF1a", "ORF1b", "ORF1ab"],
        unknown_items={
            "ORF1a/1b polyprotein": {
                "records": ["Q1"],
                "candidates": ["ORF1a/1b polyprotein"],
            },
            "contains ORF1a and ORF1b": {
                "records": ["Q2"],
                "candidates": ["contains ORF1a and ORF1b"],
            },
        },
    )

    matches = {
        row["raw_value"]: row["matching_available_canonical"]
        for row in payload["suggestions"]
    }
    assert matches["ORF1a/1b polyprotein"] == "ORF1ab"
    assert matches["contains ORF1a and ORF1b"] == "ORF1ab"
    first = payload["suggestions"][0]
    assert first["candidate_values_count"] == 1
    assert first["high_support_unresolved"] is False


def test_unresolved_payload_marks_high_support_mixed_candidates():
    payload = build_unresolved_name_review_payload(
        virus_name="test virus",
        canonical_names=["S", "M"],
        unknown_items={
            "membrane-associated protein": {
                "records": ["Q1", "Q2", "Q3"],
                "candidates": ["membrane-associated protein", "small membrane protein"],
            },
        },
    )

    row = payload["suggestions"][0]
    assert row["support_count"] == 3
    assert row["candidate_values_count"] == 2
    assert row["has_multiple_candidate_values"] is True
    assert row["high_support_unresolved"] is True


def test_review_uncertain_alias_suggestions_merges_mock_reviews():
    suggestions = [
        _row("M", "save_alias", "high", field="gene"),
        _row("unglycosylated membrane protein", "manual_review", "medium"),
    ]
    provider = MockLLMProvider([
        {
            "review_id": "placeholder",
            "recommendation": "save_alias",
            "canonical_name": "M",
            "confidence": "medium",
            "reason": "Membrane protein wording supports M in this context.",
        }
    ])

    # Prime the generated review id without duplicating its hash logic.
    primed, _ = review_uncertain_alias_suggestions(
        suggestions,
        virus_name="test virus",
        canonical_names=["M", "N"],
        config=LLMConfig(enabled=False),
        provider=MockLLMProvider([]),
    )
    review_id = primed[1]["llm_review_id"]
    provider = MockLLMProvider([
        {
            "review_id": review_id,
            "recommendation": "save_alias",
            "canonical_name": "M",
            "confidence": "medium",
            "reason": "Membrane protein wording supports M in this context.",
        }
    ])

    reviewed, diagnostics = review_uncertain_alias_suggestions(
        suggestions,
        virus_name="test virus",
        canonical_names=["M", "N"],
        config=LLMConfig(enabled=True, api_key="test"),
        provider=provider,
    )

    assert diagnostics["status"] == "reviewed"
    assert diagnostics["submitted_rows"] == 1
    assert len(provider.calls) == 1
    assert len(provider.calls[0]["suggestions"]) == 1
    assert provider.calls[0]["suggestions"][0]["raw_value"] == "unglycosylated membrane protein"
    assert reviewed[0]["llm_reviewed"] is False
    assert reviewed[1]["llm_reviewed"] is True
    assert reviewed[1]["llm_action"] == "save_alias"
    assert reviewed[1]["llm_canonical_name"] == "M"


def test_review_uncertain_alias_suggestions_uses_cache():
    suggestions = [_row("unglycosylated membrane protein", "manual_review", "medium")]
    cache = {}
    provider = MockLLMProvider([])

    first, first_diag = review_uncertain_alias_suggestions(
        suggestions,
        virus_name="test virus",
        canonical_names=["M"],
        config=LLMConfig(enabled=True, api_key="test"),
        provider=provider,
        cache=cache,
    )
    second, second_diag = review_uncertain_alias_suggestions(
        suggestions,
        virus_name="test virus",
        canonical_names=["M"],
        config=LLMConfig(enabled=True, api_key="test"),
        provider=provider,
        cache=cache,
    )

    assert first_diag["cache_hit"] is False
    assert second_diag["cache_hit"] is True
    assert len(provider.calls) == 1
    assert first[0]["llm_review_id"] == second[0]["llm_review_id"]


def test_review_unresolved_names_returns_recommendations_by_raw_name():
    unknown_items = {
        "ORF1a/1b": {
            "records": ["Q1"],
            "candidates": ["ORF1a/1b"],
        }
    }

    primed, _ = review_unresolved_names(
        unknown_items=unknown_items,
        ambiguous_items={},
        virus_name="test virus",
        canonical_names=["ORF1a", "ORF1b", "ORF1ab"],
        config=LLMConfig(enabled=False),
        provider=MockLLMProvider([]),
    )
    assert primed == {}

    payload = build_unresolved_name_review_payload(
        virus_name="test virus",
        canonical_names=["ORF1a", "ORF1b", "ORF1ab"],
        unknown_items=unknown_items,
    )
    review_id = payload["suggestions"][0]["review_id"]
    provider = MockLLMProvider([
        {
            "review_id": review_id,
            "recommendation": "save_alias",
            "canonical_name": "ORF1ab",
            "confidence": "medium",
            "reason": "Combined ORF1a/1b wording maps to ORF1ab.",
        }
    ])

    reviews, diagnostics = review_unresolved_names(
        unknown_items=unknown_items,
        ambiguous_items={},
        virus_name="test virus",
        canonical_names=["ORF1a", "ORF1b", "ORF1ab"],
        config=LLMConfig(enabled=True, api_key="test"),
        provider=provider,
    )

    assert diagnostics["status"] == "reviewed"
    assert diagnostics["submitted_rows"] == 1
    assert reviews["ORF1a/1b"]["action"] == "save_alias"
    assert reviews["ORF1a/1b"]["canonical_name"] == "ORF1ab"
    assert reviews["orf1a/1b"]["canonical_name"] == "ORF1ab"


def test_review_unresolved_names_can_be_found_by_candidate_value():
    unknown_items = {
        "combined_orf_label": {
            "records": ["Q1"],
            "candidates": ["ORF1a/1b"],
        }
    }
    payload = build_unresolved_name_review_payload(
        virus_name="test virus",
        canonical_names=["ORF1a", "ORF1b", "ORF1ab"],
        unknown_items=unknown_items,
    )
    review_id = payload["suggestions"][0]["review_id"]
    provider = MockLLMProvider([
        {
            "review_id": review_id,
            "recommendation": "save_alias",
            "canonical_name": "ORF1ab",
            "confidence": "high",
            "reason": "Combined ORF1a/1b wording maps to ORF1ab.",
        }
    ])

    reviews, diagnostics = review_unresolved_names(
        unknown_items=unknown_items,
        ambiguous_items={},
        virus_name="test virus",
        canonical_names=["ORF1a", "ORF1b", "ORF1ab"],
        config=LLMConfig(enabled=True, api_key="test"),
        provider=provider,
    )

    assert diagnostics["reviewed_rows"] == 1
    assert reviews["combined_orf_label"]["canonical_name"] == "ORF1ab"
    assert reviews["ORF1a/1b"]["canonical_name"] == "ORF1ab"
    assert reviews["orf1a/1b"]["canonical_name"] == "ORF1ab"


def test_review_uncertain_alias_suggestions_skips_when_disabled():
    suggestions = [_row("unglycosylated membrane protein", "manual_review", "medium")]
    reviewed, diagnostics = review_uncertain_alias_suggestions(
        suggestions,
        virus_name="test virus",
        canonical_names=["M"],
        config=LLMConfig(enabled=False, api_key=None),
    )

    assert diagnostics["status"] == "disabled"
    assert diagnostics["submitted_rows"] == 1
    assert reviewed[0]["llm_reviewed"] is False
