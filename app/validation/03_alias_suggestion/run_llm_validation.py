"""
Step 2: call the real LLM alias-review feature on the prepared datasets and
score it against ground truth.

USAGE (run from the viralift/ project root, with a working OPENAI_API_KEY
in .env and network access to api.openai.com):

    python app/validation/07_llm_alias_validation/build_dataset.py   # if not already run
    python app/validation/07_llm_alias_validation/run_llm_validation.py

Add --mock to do a dry run with a trivial mock provider (no network / no API
key needed) to sanity check the script mechanics without spending API calls.
"""
import argparse
import sys
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from app.src.llm.alias_review import (  # noqa: E402
    review_uncertain_alias_suggestions,
    review_unresolved_names,
)
from app.src.llm.config import LLMConfig  # noqa: E402
from app.src.llm.provider import OpenAILLMProvider  # noqa: E402

from track_a_ground_truth import (  # noqa: E402
    TRACK_A_UNKNOWN_ITEMS,
    TRACK_A_AMBIGUOUS_ITEMS,
    TRACK_A_POSITIONAL_ITEMS,
    TRACK_A_GROUND_TRUTH,
)

from app.src.alias.alias_payload import build_position_context  # noqa: E402
from app.src.alias.gene_alias import lookup_field_value  # noqa: E402
from app.src.io.genbank_parser import load_genbank_records  # noqa: E402

PED_QUERY_GB = REPO_ROOT / "app/data/PED/PED_100seqs.gb"
PED_REF_GB = REPO_ROOT / "app/data/PED/PED_ref_1.gb"
# 5' to 3' gene order of the PEDV reference, used only to turn resolved
# neighbours into a slot. Derived from the reference record, not from any
# ground-truth mapping of the names under review.
PED_GENE_ORDER = ["ORF1a", "ORF1b", "S", "ORF3", "E", "M", "N"]


def _enrich_with_position(items: dict) -> dict:
    """Attach coordinate context to unresolved-name rows.

    Without this the model sees only the annotation string, which is often
    undecidable ("sM" reads as either spike or membrane). Position is not.
    """
    from app.validation._shared.validation_utils import load_reference_bundle

    records = load_genbank_records(PED_QUERY_GB)
    alias_lookup = load_reference_bundle(PED_REF_GB)["alias_lookup"]

    def resolve(value):
        resolved = lookup_field_value(value, alias_lookup)
        if not resolved or str(resolved).startswith("__"):
            return None
        return resolved

    enriched = {}
    for name, info in (items or {}).items():
        context = build_position_context(
            records, name, resolve, ref_gene_order=PED_GENE_ORDER
        )
        enriched[name] = {**info, "position": context} if context else dict(info)
    return enriched

OUT_DIR = Path(__file__).resolve().parent / "outputs"
PED_CANONICALS = ["ORF1a", "ORF1b", "ORF1ab", "S", "ORF3", "E", "M", "N"]


def _load_dotenv(env_path: Path):
    import os
    if not env_path.exists():
        return
    for line in env_path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, _, value = line.partition("=")
        os.environ.setdefault(key.strip(), value.strip())


class _MockProvider:
    """Dry-run provider: guesses save_alias when matching_available_canonical
    is present, otherwise ignore. Only for sanity-checking the script, NOT a
    substitute for the real evaluation."""

    def review_alias_suggestions(self, payload):
        reviews = []
        for row in payload.get("suggestions", []):
            match = row.get("matching_available_canonical")
            reviews.append({
                "review_id": row.get("review_id"),
                "recommendation": "save_alias" if match else "ignore",
                "canonical_name": match,
                "confidence": "medium",
                "reason": "mock provider heuristic",
            })
        return {"reviews": reviews}


def _operative_bucket(row) -> str:
    """Replicate the exact tick logic from ui/stages/bootstrap_alias.py so
    metrics reflect what actually happens in the app, not just the raw LLM
    label."""
    if not row.get("llm_reviewed"):
        return "deterministic:" + str(row.get("suggested_action"))
    action = row.get("llm_action")
    conf = row.get("llm_confidence")
    if action == "save_alias" and conf in ("medium", "high") and row.get("llm_canonical_name") in PED_CANONICALS:
        return "save"
    if action == "ignore" and conf in ("medium", "high"):
        return "ignore"
    if action in ("skip", "move_to_ambiguous"):
        return "skip"  # NOTE: both map to the same inert bucket in the current UI
    return "deterministic:" + str(row.get("suggested_action"))


def _expected_action(status: str) -> str:
    return {
        "canonical": "save_alias",
        "ignored": "ignore",
        "ambiguous": "move_to_ambiguous",
        "unresolved": "ignore",
    }.get(status, "ignore")


def run_track_b(provider, cache):
    df = pd.read_csv(OUT_DIR / "llm_suggestion_candidates_uncertain.tsv", sep="\t")
    suggestions = df.to_dict("records")
    reviewed, diagnostics = review_uncertain_alias_suggestions(
        suggestions,
        virus_name="PED",
        canonical_names=PED_CANONICALS,
        provider=provider,
        cache=cache,
    )
    out = pd.DataFrame(reviewed)
    out["expected_action"] = out["ground_truth_status"].map(_expected_action)
    out["operative_bucket"] = out.apply(_operative_bucket, axis=1)
    out["action_correct"] = out["llm_action"] == out["expected_action"]
    canonical_check = (out["llm_action"] == "save_alias") & (out["expected_action"] == "save_alias")
    out["canonical_correct"] = None
    out.loc[canonical_check, "canonical_correct"] = (
        out.loc[canonical_check, "llm_canonical_name"] == out.loc[canonical_check, "ground_truth_resolved_to"]
    )
    out["dangerous_false_positive"] = (out["llm_action"] == "save_alias") & (out["expected_action"] != "save_alias")
    out.to_csv(OUT_DIR / "track_b_results.tsv", sep="\t", index=False)
    print("=== Track B (uncertain coordinate-supported suggestions) ===")
    print(f"Diagnostics: {diagnostics}")
    print(out[[
        "field", "raw_value", "expected_action", "llm_action", "llm_confidence",
        "operative_bucket", "action_correct", "canonical_correct", "dangerous_false_positive",
    ]].to_string())
    return out, diagnostics


def run_track_a(provider, cache, with_position: bool = True):
    unknown = {**TRACK_A_UNKNOWN_ITEMS, **TRACK_A_POSITIONAL_ITEMS}
    ambiguous = dict(TRACK_A_AMBIGUOUS_ITEMS)
    if with_position:
        unknown, ambiguous = _enrich_with_position(unknown), _enrich_with_position(ambiguous)
    reviews, diagnostics = review_unresolved_names(
        unknown_items=unknown,
        ambiguous_items=ambiguous,
        virus_name="PED",
        canonical_names=PED_CANONICALS,
        provider=provider,
        cache=cache,
    )
    rows = []
    all_items = {**unknown, **ambiguous}
    for rep in all_items:
        review = reviews.get(rep, {})
        gt = TRACK_A_GROUND_TRUTH.get(rep, {})
        rows.append({
            "raw_value": rep,
            "expected_action": gt.get("action"),
            "expected_canonical": gt.get("canonical"),
            "llm_action": review.get("action"),
            "llm_canonical_name": review.get("canonical_name"),
            "llm_confidence": review.get("confidence"),
            "llm_reason": review.get("reason"),
            "action_correct": review.get("action") == gt.get("action"),
            "canonical_correct": (
                review.get("canonical_name") == gt.get("canonical")
                if gt.get("action") == "save_alias" else None
            ),
            "has_position_evidence": bool((all_items.get(rep) or {}).get("position")),
            "dangerous_false_positive": (
                review.get("action") == "save_alias" and gt.get("action") != "save_alias"
            ),
        })
    out = pd.DataFrame(rows)
    suffix = "with_position" if with_position else "name_only"
    out.to_csv(OUT_DIR / f"track_a_results_{suffix}.tsv", sep="\t", index=False)
    print(f"\n=== Track A ({suffix}) ===")
    print(f"Diagnostics: {diagnostics}")
    print(out.to_string())
    return out, diagnostics


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mock", action="store_true", help="dry run with a mock provider, no API calls")
    args = parser.parse_args()

    _load_dotenv(REPO_ROOT / ".env")

    cache = {}
    if args.mock:
        provider = _MockProvider()
        print("Running in --mock mode (no real API calls; sanity check only).\n")
    else:
        config = LLMConfig.from_env()
        if not config.api_key:
            print("ERROR: OPENAI_API_KEY not found in environment or .env. Aborting.")
            sys.exit(1)
        provider = OpenAILLMProvider(config)
        print(f"Using real OpenAILLMProvider, model={config.model}, fallback={config.fallback_model}\n")

    b_out, b_diag = run_track_b(provider, cache)
    a_out, a_diag = run_track_a(provider, cache)

    summary = pd.DataFrame([
        {"track": "B_uncertain_suggestions", "n": len(b_out), "action_accuracy": b_out["action_correct"].mean(),
         "canonical_accuracy_given_save": b_out["canonical_correct"].dropna().mean() if b_out["canonical_correct"].notna().any() else None,
         "dangerous_false_positives": int(b_out["dangerous_false_positive"].sum())},
        {"track": "A_unresolved_names", "n": len(a_out), "action_accuracy": a_out["action_correct"].mean(),
         "canonical_accuracy_given_save": None, "dangerous_false_positives": None},
    ])
    summary.to_csv(OUT_DIR / "summary.tsv", sep="\t", index=False)
    print("\n=== Summary ===")
    print(summary.to_string())


if __name__ == "__main__":
    main()
