import pandas as pd

from ui.stages.bootstrap_alias import (
    ACTION_LABELS,
    _default_suggestion_action,
    _resolve_suggestion_default,
)


CANONICALS = ["ORF1a", "ORF1b", "ORF1ab", "E", "M", "N"]


def _row(**overrides):
    base = {
        "raw_value": "raw",
        "canonical_name": "M",
        "suggested_action": "manual_review",
        "default_save": False,
        "llm_reviewed": False,
        "llm_action": None,
        "llm_confidence": None,
        "llm_canonical_name": None,
    }
    base.update(overrides)
    return pd.Series(base)


def test_default_action_falls_back_to_deterministic_save():
    row = _row(suggested_action="save_alias", default_save=True)
    assert _default_suggestion_action(row, CANONICALS) == "save"


def test_default_action_falls_back_to_deterministic_ignore():
    row = _row(suggested_action="ignore", default_save=False)
    assert _default_suggestion_action(row, CANONICALS) == "ignore"


def test_default_action_falls_back_to_skip_when_no_strong_signal():
    row = _row(suggested_action="manual_review", default_save=False)
    assert _default_suggestion_action(row, CANONICALS) == "skip"


def test_llm_save_wins_only_at_medium_or_high_confidence():
    high = _row(
        suggested_action="ignore", default_save=False,
        llm_reviewed=True, llm_action="save_alias", llm_confidence="high",
        llm_canonical_name="M",
    )
    low = _row(
        suggested_action="ignore", default_save=False,
        llm_reviewed=True, llm_action="save_alias", llm_confidence="low",
        llm_canonical_name="M",
    )
    assert _default_suggestion_action(high, CANONICALS) == "save"
    # Low-confidence LLM save must not override the deterministic default.
    assert _default_suggestion_action(low, CANONICALS) == "ignore"


def test_llm_save_ignored_when_canonical_not_available():
    row = _row(
        suggested_action="ignore", default_save=False,
        llm_reviewed=True, llm_action="save_alias", llm_confidence="high",
        llm_canonical_name="NOT_A_REAL_CANONICAL",
    )
    assert _default_suggestion_action(row, CANONICALS) == "ignore"


def test_legacy_llm_move_to_ambiguous_maps_to_exclude():
    row = _row(
        suggested_action="manual_review", default_save=False,
        llm_reviewed=True, llm_action="move_to_ambiguous", llm_confidence="medium",
    )
    assert _default_suggestion_action(row, CANONICALS) == "ignore"


def test_llm_skip_wins_regardless_of_confidence():
    row = _row(
        suggested_action="save_alias", default_save=True,
        llm_reviewed=True, llm_action="skip", llm_confidence="low",
    )
    # Even though the deterministic scorer said save, an LLM skip (any
    # confidence) always wins — matches the original checkbox-era behavior
    # where "skip" was applied without a confidence gate.
    assert _default_suggestion_action(row, CANONICALS) == "skip"


def test_resolve_suggestion_default_swaps_canonical_on_llm_save():
    row = _row(
        canonical_name="ORF1a", suggested_action="save_alias", default_save=True,
        llm_reviewed=True, llm_action="save_alias", llm_confidence="high",
        llm_canonical_name="ORF1ab",
    )
    resolved = _resolve_suggestion_default(row, CANONICALS)
    assert resolved["action"] == ACTION_LABELS["save"]
    assert resolved["canonical_name"] == "ORF1ab"


def test_resolve_suggestion_default_keeps_original_canonical_when_not_llm_save():
    row = _row(canonical_name="M", suggested_action="ignore", default_save=False)
    resolved = _resolve_suggestion_default(row, CANONICALS)
    assert resolved["action"] == ACTION_LABELS["ignore"]
    assert resolved["canonical_name"] == "M"
