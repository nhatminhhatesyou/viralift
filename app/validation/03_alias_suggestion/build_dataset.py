"""
Step 1 (offline, no API calls): build a coordinate-supported suggestion
dataset for the LLM alias-review validation, reusing already-validated PED
data from app/validation/06_ped_validation instead of re-running tblastn.

Rationale:
- Coordinate accuracy (tblastn) is already validated separately at ~99%
  (see PED_VALIDATION_REPORT.md section 4). This script isolates the
  *naming/classification* layer that the LLM review feature actually sits on
  top of, so we don't need to re-run tblastn to test it meaningfully.
- ped_alias_qualifier_rows.tsv already has, per query record and per
  qualifier field, the raw annotation text and its already-validated
  ground-truth resolution (resolved_to / status).
- We derive a real coordinate_canonical + iou per row by overlapping the
  query feature's coordinates against the ref_1 gene model
  (ped_ref_model.tsv), i.e. simulating what tblastn+IoU matching would have
  produced, using real coordinates rather than fabricated numbers.
"""
import sys
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT))

from app.src.alias.alias_classifier import classify_alias_candidate  # noqa: E402
from app.src.llm.alias_review import needs_llm_review  # noqa: E402
from app.src.alias.gene_alias import normalize_text  # noqa: E402

PED_DIR = REPO_ROOT / "app" / "validation" / "06_ped_validation" / "outputs" / "alias"
OUT_DIR = Path(__file__).resolve().parent / "outputs"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def compute_iou(a_start, a_end, b_start, b_end) -> float:
    a_start, a_end = sorted([a_start, a_end])
    b_start, b_end = sorted([b_start, b_end])
    inter = max(0, min(a_end, b_end) - max(a_start, b_start))
    union = (a_end - a_start) + (b_end - b_start) - inter
    if union <= 0:
        return 0.0
    return inter / union


def main():
    qual_rows = pd.read_csv(PED_DIR / "ped_alias_qualifier_rows.tsv", sep="\t")
    ref_model = pd.read_csv(PED_DIR / "ped_ref_model.tsv", sep="\t")
    ref1 = ref_model[ref_model["ref"] == "ref_1"]

    query = qual_rows[qual_rows["dataset"] == "query_100"].copy()

    # Per-row coordinate canonical + iou vs ref_1 gene model (best-IoU gene).
    def best_ref_match(row):
        best_name, best_iou = None, 0.0
        for _, gene in ref1.iterrows():
            iou = compute_iou(row["feature_start"], row["feature_end"], gene["start"], gene["end"])
            if iou > best_iou:
                best_iou, best_name = iou, gene["name"]
        return pd.Series({"coordinate_canonical": best_name, "iou": round(best_iou, 4)})

    query = query.join(query.apply(best_ref_match, axis=1))

    # Only keep rows that would have cleared the real min_iou=0.90 bootstrap
    # threshold (build_coordinate_supported_alias_suggestions default).
    query = query[query["iou"] >= 0.90].copy()

    # Aggregate to one suggestion row per (field, raw_value), like
    # deduplicate_suggestions() does in the real pipeline.
    def agg(group):
        return pd.Series({
            "canonical_name": group["coordinate_canonical"].mode().iat[0],
            "iou": group["iou"].max(),
            "ground_truth_resolved_to": group["resolved_to"].mode().iat[0],
            "ground_truth_status": group["status"].mode().iat[0],
            "support_count": group["record_id"].nunique(),
            "support_records": ", ".join(sorted(group["record_id"].unique())[:10]),
        })

    grouped = (
        query.groupby(["field", "raw_value"], as_index=False)
        .apply(agg, include_groups=False)
        .reset_index(drop=True)
    )

    rows = []
    for _, r in grouped.iterrows():
        classified = classify_alias_candidate(
            raw_value=r["raw_value"],
            field=r["field"],
            canonical_name=r["canonical_name"],
            evidence={"iou": r["iou"], "strand_match": True},
        )
        rows.append({
            "field": r["field"],
            "raw_value": r["raw_value"],
            "canonical_name": r["canonical_name"],
            "iou": r["iou"],
            "support_count": r["support_count"],
            "support_records": r["support_records"],
            "ground_truth_resolved_to": r["ground_truth_resolved_to"],
            "ground_truth_status": r["ground_truth_status"],
            **classified,
        })

    df = pd.DataFrame(rows)
    df["needs_llm_review"] = df.apply(lambda r: needs_llm_review(r.to_dict()), axis=1)

    df.to_csv(OUT_DIR / "llm_suggestion_candidates_full.tsv", sep="\t", index=False)
    uncertain = df[df["needs_llm_review"]].copy()
    uncertain.to_csv(OUT_DIR / "llm_suggestion_candidates_uncertain.tsv", sep="\t", index=False)

    print(f"Total distinct (field, raw_value) suggestion rows (iou>=0.90): {len(df)}")
    print(f"Deterministic action distribution:\n{df['suggested_action'].value_counts()}")
    print(f"\nFlagged as needs_llm_review: {len(uncertain)} / {len(df)}")
    print(uncertain[["field", "raw_value", "canonical_name", "suggested_action", "confidence", "ground_truth_resolved_to", "ground_truth_status"]].to_string())


if __name__ == "__main__":
    main()
