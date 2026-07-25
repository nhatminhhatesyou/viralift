"""
Alias-config reconstruction validation.

Question this answers: if the tool is handed a brand-new virus (only the
reference gene names, no alias map), how close does its coordinate-driven
suggestion feature come to reproducing the curated "gold" alias config?

Design (three steps, mirrors the agreed plan):

  1. SEED   — strip the gold config down to canonical names only, exactly the
              state a user faces on an unseen virus. `build_seed_alias_config_from_ref`.
  2. SUGGEST— run the real `build_coordinate_supported_alias_suggestions` from
              that seed. This is the tool doing its job: query annotations +
              tblastn coordinate evidence, no gold names leaked in.
  3. DIFF   — for every raw annotation string in the corpus, compare where the
              suggestion maps it vs where the gold config maps it, bucketed by
              error type (the buckets have very different real-world costs).

Why not circular: the gold config is curated independently (expert + literature
+ coordinates). The suggestion only ever sees the seed. The diff measures the
tool's naming/decision layer, not the coordinate layer (already validated at
~99.7% in 02_lifting_accuracy).

Requires a working `tblastn` on PATH (step 2 runs real lifting). Run from the
viralift/ project root.
"""
import argparse
import sys
from collections import Counter
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT))

from app.src.alias.alias_bootstrap import (  # noqa: E402
    apply_approved_alias_suggestions,
    build_coordinate_supported_alias_suggestions,
    build_seed_alias_config_from_ref,
    write_new_alias_config,
)
from app.src.alias.gene_alias import (  # noqa: E402
    build_alias_lookup,
    normalize_text,
)
from app.src.features.ref_loader import prepare_reference_features  # noqa: E402
from app.src.io.genbank_parser import load_genbank_records  # noqa: E402
from app.validation._shared.validation_utils import collect_alias_rows  # noqa: E402
from app.src.llm.alias_review import review_uncertain_alias_suggestions  # noqa: E402
from app.src.llm.config import LLMConfig  # noqa: E402
import json
import tempfile

OUT_DIR = Path(__file__).resolve().parent / "outputs"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# One entry per virus: gold config + reference + query corpus.
DATASETS = {
    "PEDV": {
        "gold": "app/config/porcine_epidemic_diarrhea_virus_alias.json",
        "ref": "app/data/PED/PED_ref_1.gb",
        "query": "app/data/PED/PED_100seqs.gb",
    },
    "PRRSV": {
        "gold": "app/config/prrsv_alias.json",
        "ref": "app/data/PRRS/PRRS_ref_test.gb",
        "query": "app/data/PRRS/PRRS_100seq_anno.gb",
    },
    "FMDV": {
        "gold": "app/config/fmdv_alias.json",
        "ref": "app/data/FMD/FMD_ref_test.gb",
        "query": "app/data/FMD/FMD_100seq_anno.gb",
    },
}


def gold_mapping(gold_config: dict) -> dict:
    """{normalised raw name -> gold verdict}. Verdict is a canonical name or 'excluded'."""
    lookup = build_alias_lookup(gold_config)
    mapping = {}
    for norm, canonical in lookup.items():
        # Sentinels (__excluded__/__ambiguous__) collapse to a single 'excluded'
        # bucket: for reconstruction we only care "is this a usable gene alias".
        mapping[norm] = canonical if not str(canonical).startswith("__") else "excluded"
    for name in gold_config.get("excluded_names", []) or []:
        mapping.setdefault(_norm(name), "excluded")
    return mapping


def _norm(value: str) -> str:
    return normalize_text(value)


# --- Auto-approval policy -------------------------------------------------
# The real flow has a human step (the user ticks save/exclude/skip in the UI)
# that cannot be automated. The harness stops at the tool's *recommendation*
# and applies one explicit, declared policy in the human's place:
#
#   accept every row the tool would auto-save; exclude nothing on the tool's
#   behalf; leave everything else unapproved.
#
# The resulting number is therefore an UPPER BOUND: "how well would the config
# turn out if the user accepted exactly the tool's save recommendations". It is
# not the outcome for a real, selective user. This assumption is stated in the
# notebook alongside the metric.
def _operative_action(row: dict) -> str:
    """
    The action that actually takes effect, mirroring ui/stages/bootstrap_alias.py.

    Deterministic rows use `suggested_action`; rows the LLM reviewed use
    `llm_action` (medium/high confidence) — the same precedence the UI applies
    when it pre-ticks the suggestion table.
    """
    if row.get("llm_reviewed"):
        action = row.get("llm_action")
        conf = row.get("llm_confidence")
        if action == "save_alias" and conf in ("medium", "high"):
            return "save_alias"
        if action == "ignore" and conf in ("medium", "high"):
            return "ignore"
        if action in ("skip", "move_to_ambiguous"):
            return "skip"
    return row.get("suggested_action") or "skip"


def build_reconstructed_config(
    seed_config_path: Path,
    reviewed_suggestions: list,
) -> dict:
    """
    Turn the tool's reviewed suggestions into an actual config, via the SAME
    function the UI calls (`apply_approved_alias_suggestions`), under the
    declared auto-approval policy above. No reimplementation of tool logic.
    """
    approved = [r for r in reviewed_suggestions if _operative_action(r) == "save_alias"]
    excluded = [r for r in reviewed_suggestions if _operative_action(r) == "ignore"]
    return apply_approved_alias_suggestions(
        seed_config_path,
        approved_rows=approved,
        excluded_rows=excluded,
    )


def build_config_temp(name: str, cfg: dict, min_iou: float, llm_provider=None):
    """
    Run the real suggestion pipeline for one virus from a bare seed and return
    (config_temp dict, corpus_norms set).

    config_temp is what the tool would produce for an unseen virus. corpus_norms
    is every normalised annotation string that actually occurs in the query
    records — the recall denominator gate (a truth alias absent from the corpus
    cannot possibly be suggested, so it must not count as a miss).
    """
    ref_record = load_genbank_records(REPO_ROOT / cfg["ref"])[0]
    query_records = load_genbank_records(REPO_ROOT / cfg["query"])

    ref_features, ref_ft, _, _, _ = prepare_reference_features(
        ref_record,
        alias_config_arg=REPO_ROOT / cfg["gold"],
        alias_registry_arg=REPO_ROOT / "app/config/virus_alias_registry.json",
    )
    # STEP 1 — seed: canonical names only, no aliases (unseen-virus state).
    seed = build_seed_alias_config_from_ref(ref_record, ref_features, virus_name=name)
    seed_canonicals = list(seed["canonical_names"].keys())
    tmp_dir = Path(tempfile.mkdtemp())
    seed_path = tmp_dir / "seed.json"
    write_new_alias_config(seed, seed_path)

    # STEP 2 — real suggestion feature (coordinate layer) + real LLM-review layer.
    suggestions = build_coordinate_supported_alias_suggestions(
        ref_record=ref_record,
        query_records=query_records,
        ref_features=ref_features,
        ref_feature_type=ref_ft,
        seed_canonical_names=seed_canonicals,
        min_iou=min_iou,
    )
    # Production caps LLM review at config.max_rows (default 20) per pass — a
    # first-time bootstrap of a virus with many uncertain names needs several
    # passes. For validation we want the tool's full capability in one shot, so
    # raise the cap to cover every uncertain row; otherwise rows beyond the cap
    # stay at 'manual_review' and, under the accept-recommendation policy, fall
    # through to neither saved nor excluded (a pagination artifact, not a tool
    # error). The 20-row cap itself is reported as a product finding.
    import dataclasses
    review_config = dataclasses.replace(LLMConfig.from_env(), max_rows=100000)
    reviewed, _ = review_uncertain_alias_suggestions(
        suggestions,
        virus_name=name,
        canonical_names=seed_canonicals,
        excluded_names=seed.get("excluded_names", []),
        config=review_config,
        provider=llm_provider,
    )

    # Diagnostic dump: what happened to every surfaced name — was it reviewed by
    # the LLM, and what did it decide? Lets us see per-name why something ended
    # up saved / excluded / None, instead of guessing.
    diag = pd.DataFrame([{
        "raw_value": r.get("raw_value"),
        "norm": _norm(r.get("raw_value", "")),
        "coord_canonical": r.get("canonical_name"),
        "iou": r.get("iou"),
        "cross_canon_n": r.get("cross_canonical_target_count"),
        "det_action": r.get("suggested_action"),
        "llm_reviewed": r.get("llm_reviewed"),
        "llm_action": r.get("llm_action"),
        "llm_confidence": r.get("llm_confidence"),
        "operative": _operative_action(r),
        "llm_reason": r.get("llm_reason"),
    } for r in reviewed])
    diag.to_csv(OUT_DIR / f"reviewed_{name}.tsv", sep="\t", index=False)

    # STEP 3 — build config_temp via the real config-builder, save it.
    config_temp = build_reconstructed_config(seed_path, reviewed)
    out_path = OUT_DIR / f"config_temp_{name}.json"
    out_path.write_text(json.dumps(config_temp, indent=2, ensure_ascii=False) + "\n")

    corpus_norms = {
        _norm(r.raw_value)
        for r in collect_alias_rows(query_records).itertuples()
    }
    return config_temp, corpus_norms


def compare_per_canonical(config_temp: dict, config_truth: dict, corpus_norms: set):
    """
    Two-directional per-canonical comparison (config_temp vs config_truth).

    PRECISION (iterate config_temp): each name the tool put under canonical X —
    does truth agree it belongs to X?
        correct       tool→X, truth→X
        wrong_gene    tool→X, truth→Y≠X          (mis-assigned — costly)
        false_save    tool→X, truth excluded it  (saved a non-alias — costly)
        not_in_truth  tool→X, truth silent       (truth may be incomplete — check)

    RECALL (iterate config_truth, gated by corpus): each real alias truth has
    under X that ACTUALLY occurs in the 100 records — did config_temp keep it at X?
        found         truth→X, temp→X
        missed        truth→X, temp did not      (sM-type failure surfaces here)
      A truth alias absent from the corpus is parked, never scored.

    Both sides read verdicts through build_alias_lookup on real configs.
    """
    temp_lu = build_alias_lookup(config_temp)
    truth_lu = build_alias_lookup(config_truth)

    def verdict(lu, norm):
        v = lu.get(norm)
        if v is None:
            return None
        return v if not str(v).startswith("__") else "excluded"

    # invert temp lookup to {canonical: [names]} for precision
    temp_by_canon = {}
    for norm, canon in temp_lu.items():
        if canon and not str(canon).startswith("__"):
            temp_by_canon.setdefault(canon, []).append(norm)
    truth_by_canon = {}
    for norm, canon in truth_lu.items():
        if canon and not str(canon).startswith("__"):
            truth_by_canon.setdefault(canon, []).append(norm)

    # Canonicals the tool could actually work with: those present in the seed,
    # i.e. the reference gene set (config_temp always carries the seed canonicals).
    # A truth canonical finer than the reference (e.g. FMDV 3B1/3B2/3B3 when the
    # ref only annotates 3B) can never be reconstructed from the ref — its aliases
    # are parked, not counted as misses.
    seed_canonicals = set(config_temp.get("canonical_names", {}))

    canonicals = sorted(set(temp_by_canon) | set(truth_by_canon))
    rows, detail = [], []
    for c in canonicals:
        # precision
        p_correct = p_wrong = p_false = p_notintruth = 0
        for n in temp_by_canon.get(c, []):
            tv = verdict(truth_lu, n)
            if tv == c:
                p_correct += 1; tag = "correct"
            elif tv == "excluded":
                p_false += 1; tag = "false_save"
            elif tv is None:
                p_notintruth += 1; tag = "not_in_truth"
            else:
                p_wrong += 1; tag = f"wrong_gene({tv})"
            if tag != "correct":
                detail.append({"canonical": c, "side": "precision", "name": n, "issue": tag})
        # recall (gated by: alias occurs in corpus AND canonical exists in ref)
        r_found = r_missed = r_parked = 0
        canonical_reconstructable = c in seed_canonicals
        for n in truth_by_canon.get(c, []):
            if not canonical_reconstructable:
                r_parked += 1
                detail.append({"canonical": c, "side": "recall", "name": n,
                               "issue": "parked_canonical_not_in_ref"})
                continue
            if n not in corpus_norms:
                r_parked += 1
                detail.append({"canonical": c, "side": "recall", "name": n, "issue": "parked_not_in_corpus"})
                continue
            if verdict(temp_lu, n) == c:
                r_found += 1
            else:
                r_missed += 1
                detail.append({"canonical": c, "side": "recall", "name": n,
                               "issue": f"missed(temp→{verdict(temp_lu, n)})"})
        temp_total = len(temp_by_canon.get(c, []))
        recall_denom = r_found + r_missed
        rows.append({
            "canonical": c,
            "temp_saved": temp_total,
            "precision_correct": p_correct,
            "wrong_gene": p_wrong,
            "false_save": p_false,
            "not_in_truth": p_notintruth,
            "precision_pct": round(p_correct / temp_total * 100, 1) if temp_total else None,
            "recall_denom": recall_denom,
            "recall_found": r_found,
            "missed": r_missed,
            "recall_pct": round(r_found / recall_denom * 100, 1) if recall_denom else None,
            "parked": r_parked,
        })
    return pd.DataFrame(rows), pd.DataFrame(detail)


def run_virus(name: str, cfg: dict, min_iou: float, llm_provider=None):
    """Full run for one virus: returns (per_canonical_df, detail_df)."""
    config_temp, corpus = build_config_temp(name, cfg, min_iou, llm_provider)
    truth = json.loads((REPO_ROOT / cfg["gold"]).read_text())
    per_canon, detail = compare_per_canonical(config_temp, truth, corpus)
    per_canon.to_csv(OUT_DIR / f"reconstruction_{name}.tsv", sep="\t", index=False)
    detail.to_csv(OUT_DIR / f"reconstruction_{name}_detail.tsv", sep="\t", index=False)
    return per_canon, detail


def summarize(name: str, per_canon: pd.DataFrame) -> None:
    tot = per_canon[["temp_saved", "precision_correct", "wrong_gene", "false_save",
                     "not_in_truth", "recall_denom", "recall_found", "missed", "parked"]].sum()
    prec = tot["precision_correct"] / tot["temp_saved"] * 100 if tot["temp_saved"] else 0
    rec = tot["recall_found"] / tot["recall_denom"] * 100 if tot["recall_denom"] else 0
    print(f"\n===== {name} =====")
    print(f"  PRECISION (config_temp → truth): {int(tot['precision_correct'])}/{int(tot['temp_saved'])} "
          f"= {prec:.1f}%")
    print(f"     wrong_gene {int(tot['wrong_gene'])}  false_save {int(tot['false_save'])}  "
          f"not_in_truth {int(tot['not_in_truth'])}")
    print(f"  RECALL (truth ∩ corpus → config_temp): {int(tot['recall_found'])}/{int(tot['recall_denom'])} "
          f"= {rec:.1f}%")
    print(f"     missed {int(tot['missed'])}  (parked, not scored: {int(tot['parked'])})")


class _MockLLMProvider:
    """Dry-run LLM stand-in: mirrors run_llm_validation._MockProvider.

    Approves save_alias when the row already carries a matching canonical from
    the coordinate layer, else ignore. NOT a substitute for the real model —
    only lets the harness exercise the full call path without an API key. The
    notebook labels any run using this as 'deterministic + mock LLM'.
    """

    def review_alias_suggestions(self, payload):
        reviews = []
        for row in payload.get("suggestions", []):
            match = row.get("matching_available_canonical") or row.get("canonical_candidate")
            reviews.append({
                "review_id": row.get("review_id"),
                "recommendation": "save_alias" if match else "skip",
                "canonical_name": match,
                "confidence": "medium",
                "reason": "mock provider heuristic",
            })
        return {"reviews": reviews}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--virus", choices=list(DATASETS), help="one virus (default: all)")
    ap.add_argument("--min-iou", type=float, default=0.90)
    ap.add_argument("--mock", action="store_true",
                    help="use a mock LLM (no API key); result is deterministic+mock, not full")
    args = ap.parse_args()

    provider = None
    if args.mock:
        provider = _MockLLMProvider()
        print("Running with MOCK LLM (no API calls). Result is an approximation.\n")
    elif not LLMConfig.from_env().available:
        provider = _MockLLMProvider()
        print("No OPENAI_API_KEY found — falling back to MOCK LLM. Set a key for the real run.\n")

    targets = [args.virus] if args.virus else list(DATASETS)
    for name in targets:
        try:
            per_canon, _ = run_virus(name, DATASETS[name], args.min_iou, llm_provider=provider)
            summarize(name, per_canon)
        except Exception as exc:  # noqa: BLE001
            import traceback
            print(f"\n===== {name} ===== FAILED: {type(exc).__name__}: {exc}")
            traceback.print_exc()
    print(f"\nPer-canonical TSVs written to {OUT_DIR}/reconstruction_<virus>.tsv")


if __name__ == "__main__":
    main()
