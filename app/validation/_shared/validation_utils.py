from __future__ import annotations

import sys
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import pandas as pd
from Bio.SeqRecord import SeqRecord


ROOT = Path(__file__).resolve().parents[3]
APP = ROOT / "app"
DATA = APP / "data"
CROSS_CHECK = DATA / "cross_check"
CONFIG = APP / "config"

if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from app.src.alias.alias_registry import detect_alias_config_for_record, get_detected_virus_name
from app.src.alias.gene_alias import (
    EXCLUSION_SENTINELS,
    apply_alias_to_feature,
    apply_alias_to_features,
    load_alias_lookup,
    lookup_field_value,
    normalize_text,
)
from app.src.features.annotation_strategy import get_strategy, select_feature_type
from app.src.features.direct_extractor import direct_extract_with_alias
from app.src.features.ref_loader import prepare_reference_features
from app.src.io.genbank_parser import (
    _LOOKUP_QUALIFIER_KEYS,
    load_genbank_records,
    load_single_genbank,
    parse_cds_features,
    parse_mat_peptides,
)
from app.src.lifting.tblastn_lifter import lift_all_tblastn, process_one_query_record
from app.src.pipeline import run_pipeline, PipelineConfig


def iou(a_start: Optional[int], a_end: Optional[int], b_start: Optional[int], b_end: Optional[int]) -> float:
    if None in (a_start, a_end, b_start, b_end):
        return 0.0
    left = max(int(a_start), int(b_start))
    right = min(int(a_end), int(b_end))
    inter = max(0, right - left + 1)
    union = max(int(a_end), int(b_end)) - min(int(a_start), int(b_start)) + 1
    return inter / union if union else 0.0


def filter_subfeatures(features: Sequence[Dict], max_ratio: float = 0.8) -> Tuple[List[Dict], List[Dict]]:
    main, subs = [], []
    for i, feature in enumerate(features):
        feature_len = feature["end"] - feature["start"] + 1
        is_subfeature = any(
            j != i
            and feature["start"] >= parent["start"]
            and feature["end"] <= parent["end"]
            and feature_len < max_ratio * (parent["end"] - parent["start"] + 1)
            for j, parent in enumerate(features)
        )
        (subs if is_subfeature else main).append(feature)
    return main, subs


def parse_features_for_type(record: SeqRecord, feature_type: Optional[str]) -> List[Dict]:
    if feature_type == "mat_peptide":
        return parse_mat_peptides(record)
    if feature_type == "CDS":
        return parse_cds_features(record)
    return []


def load_reference_bundle(ref_path: Path, registry_path: Path = CONFIG / "virus_alias_registry.json") -> Dict:
    ref_record = load_single_genbank(ref_path)
    alias_config_path = detect_alias_config_for_record(ref_record, registry_path)
    if alias_config_path is not None and not alias_config_path.is_absolute():
        alias_config_path = ROOT / alias_config_path
    ref_features, ref_feature_type, alias_config_path, virus_name, alias_lookup = prepare_reference_features(
        ref_record=ref_record,
        alias_config_arg=str(alias_config_path) if alias_config_path else None,
        alias_registry_arg=str(registry_path),
    )
    return {
        "record": ref_record,
        "features": ref_features,
        "feature_type": ref_feature_type,
        "alias_config_path": alias_config_path,
        "alias_lookup": alias_lookup,
        "virus_name": virus_name or get_detected_virus_name(ref_record, registry_path),
    }


def parse_truth_features(
    record: SeqRecord,
    alias_lookup: Dict[str, str],
    feature_type: Optional[str],
    filter_nested: bool = True,
    target_names: Optional[Sequence[str]] = None,
    keep_extra_names: Optional[Sequence[str]] = None,
) -> Tuple[List[Dict], List[Dict]]:
    features = parse_features_for_type(record, feature_type)
    if alias_lookup and features:
        features = apply_alias_to_features(features, alias_lookup)
    if target_names is not None:
        allowed = set(target_names)
        if keep_extra_names:
            allowed.update(keep_extra_names)
        kept = [feature for feature in features if feature.get("name") in allowed]
        dropped = [feature for feature in features if feature.get("name") not in allowed]
        return kept, dropped
    if not filter_nested:
        return features, []
    return filter_subfeatures(features)


def truth_target_names(ref_features: Sequence[Dict], virus_label: Optional[str] = None) -> Tuple[List[str], List[str]]:
    # Truth targets = the gene names the reference actually carries (set R).
    # Accuracy denominator is R ∩ query-truth: a gene is only evaluable on a
    # record whose truth contains it. Granularity/merge cases (e.g. a query
    # ORF1ab where the ref has ORF1a/ORF1b) fall outside R and are reported
    # separately by the failure-attribution notebook via a generic overlap
    # rule — no virus- or gene-name is hard-coded here.
    target_names = sorted({feature["name"] for feature in ref_features if feature.get("name")})
    return target_names, []


def should_use_target_truth_filter(feature_type: Optional[str], virus_label: Optional[str] = None) -> bool:
    # Apply the ref-gene truth filter for coding-sequence evaluation, generically
    # for every virus (virus_label kept only for signature compatibility).
    return feature_type == "CDS"


def parse_validation_truth_features(
    record: SeqRecord,
    alias_lookup: Dict[str, str],
    feature_type: Optional[str],
    ref_features: Sequence[Dict],
    virus_label: Optional[str] = None,
) -> Tuple[List[Dict], List[Dict]]:
    if should_use_target_truth_filter(feature_type, virus_label):
        targets, extras = truth_target_names(ref_features, virus_label)
        return parse_truth_features(
            record,
            alias_lookup,
            feature_type,
            filter_nested=False,
            target_names=targets,
            keep_extra_names=extras,
        )
    return parse_truth_features(record, alias_lookup, feature_type)


def lifted_to_rows(record_id: str, lifted: Iterable, method: str) -> List[Dict]:
    rows = []
    for item in lifted:
        if hasattr(item, "to_dict"):
            data = item.to_dict()
        else:
            data = dict(item)
        rows.append({
            "record_id": record_id,
            "method": method,
            "ref_name": data.get("ref_name") or data.get("name"),
            "pred_name": data.get("name"),
            "source_name": data.get("source_name"),
            "pred_start": data.get("start") or data.get("query_start"),
            "pred_end": data.get("end") or data.get("query_end"),
            "strand": data.get("strand"),
            "status": data.get("status"),
            "coverage": data.get("coverage"),
            "identity": data.get("identity"),
            "score": data.get("score"),
            # biological validity the tool already computed — carried through so the
            # coord-correct metric can require a valid CDS, not just a good IoU.
            "has_start_codon": data.get("has_start_codon"),
            "has_stop_codon": data.get("has_stop_codon"),
            "in_frame": data.get("in_frame"),
            # diagnostics: where the raw HSP landed vs the final coords, and rescue actions
            "raw_start": data.get("raw_start"),
            "raw_end": data.get("raw_end"),
            "n_term_missing_aa": data.get("n_term_missing_aa"),
            "rescue_target": data.get("rescue_target"),
            "rescue_offset": data.get("rescue_offset"),
            "rescue_action": data.get("rescue_action"),
        })
    return rows


def best_overlap(pred_row: Dict, truth_features: Sequence[Dict]) -> Tuple[Optional[Dict], float]:
    best_feature, best_iou = None, 0.0
    for truth in truth_features:
        score = iou(pred_row.get("pred_start"), pred_row.get("pred_end"), truth["start"], truth["end"])
        if score > best_iou:
            best_feature, best_iou = truth, score
    return best_feature, best_iou


def dedupe_truth_by_name(truth_features: Sequence[Dict]) -> List[Dict]:
    """Collapse truth features that share a name, keeping the longest.

    `compare_predictions_to_truth` already resolves truth by name this way (see its
    `truth_by_name` loop): at most ONE truth feature per name can ever be matched. A
    truth-anchored harness that counts every raw truth feature therefore inflates its own
    denominator whenever a record annotates the same gene twice -- the duplicate is an
    automatic miss for every tool, by construction rather than by performance.

    Real case: PRRSV `AF331831.1` labels BOTH ORF1a and ORF1b "RNA polymerase", so ORF1b
    appears twice in truth after alias resolution. Counting it twice made the comparison
    denominator 797 where the accuracy harness (which counts (record, gene) presence via
    `build_truth_presence`) counts 796.

    Use this in truth-anchored harnesses that build their denominator from the truth list, so
    the denominator matches what the comparator can actually match.
    """
    best: Dict[str, Dict] = {}
    for feature in truth_features:
        existing = best.get(feature["name"])
        if existing is None or (feature["end"] - feature["start"]) > (existing["end"] - existing["start"]):
            best[feature["name"]] = feature
    return [feature for feature in truth_features if best.get(feature["name"]) is feature]


def compare_predictions_to_truth(
    predictions: Sequence[Dict],
    truth_features: Sequence[Dict],
    iou_threshold: float = 0.90,
    codon_required_names: Optional[set] = None,
    bp_tolerance: int = 6,
) -> pd.DataFrame:
    # codon_required_names: genes whose REFERENCE feature is itself a clean CDS
    # (valid start/stop/frame). Only these are held to the codon check. Genes whose
    # reference is partial / frameshift / overlapping (e.g. PRRSV ORF1b, which has no
    # independent ATG start) are exempt -- an "invalid codon" there is biology, not a
    # tool error. None => hold every codon-computable prediction to the check.
    truth_by_name: Dict[str, Dict] = {}
    for truth in truth_features:
        existing = truth_by_name.get(truth["name"])
        if existing is None or (truth["end"] - truth["start"]) > (existing["end"] - existing["start"]):
            truth_by_name[truth["name"]] = truth

    rows = []
    for pred in predictions:
        same_name_truth = truth_by_name.get(pred.get("pred_name"))
        overlap_truth, best_iou = best_overlap(pred, truth_features)
        truth = same_name_truth
        match_iou = (
            iou(pred.get("pred_start"), pred.get("pred_end"), truth["start"], truth["end"])
            if truth else 0.0
        )
        # Codon/frame validity: the lifter reports has_start_codon / has_stop_codon /
        # in_frame (None for non-CDS such as mat_peptide, where it does not apply).
        # A prediction is codon-valid unless one of these is explicitly False.
        codon_valid = not (
            pred.get("has_start_codon") is False
            or pred.get("has_stop_codon") is False
            or pred.get("in_frame") is False
        )
        # Only hold this gene to the codon check if its reference is a clean CDS.
        codon_required = codon_required_names is None or pred.get("pred_name") in codon_required_names
        codon_ok = codon_valid or not codon_required
        exact = (
            truth is not None
            and pred.get("pred_start") == truth["start"]
            and pred.get("pred_end") == truth["end"]
        )
        # Length / boundary diagnostics (for coord-only & failure inspection).
        ps, pe = pred.get("pred_start"), pred.get("pred_end")
        ts, te = (truth["start"], truth["end"]) if truth else (None, None)
        pred_len = (pe - ps + 1) if ps is not None and pe is not None else None
        truth_len = (te - ts + 1) if truth else None
        d_start = (ps - ts) if (truth and ps is not None) else None
        d_end = (pe - te) if (truth and pe is not None) else None

        # Absolute boundary tolerance. IoU is scale-dependent: the SAME few-bp
        # annotation-convention difference leaves a long feature at IoU~0.99 (coord-only)
        # but pushes a very short one below the cutoff (failed) -- e.g. a 2-codon
        # cleavage-boundary convention on an 18 aa peptide. A prediction whose BOTH
        # boundaries sit within `bp_tolerance` of truth is coordinate-correct regardless
        # of feature length. Generic: no gene or virus name, no truth used beyond the
        # comparison itself.
        within_bp_tol = (
            d_start is not None and d_end is not None
            and abs(d_start) <= bp_tolerance and abs(d_end) <= bp_tolerance
        )
        # Coordinate-correct = right region (IoU >= threshold OR within bp tolerance)
        # AND a valid CDS (start + stop + frame) WHERE the reference gene warrants it.
        coord_correct = truth is not None and (
            exact or ((match_iou >= iou_threshold or within_bp_tol) and codon_ok)
        )
        rows.append({
            **pred,
            "truth_name": truth["name"] if truth else None,
            "truth_start": truth["start"] if truth else None,
            "truth_end": truth["end"] if truth else None,
            "best_overlap_name": overlap_truth["name"] if overlap_truth else None,
            "best_iou": round(best_iou, 4),
            "iou": round(match_iou, 4),
            "name_match": truth is not None,
            "codon_valid": codon_valid,
            "codon_required": codon_required,
            "within_bp_tolerance": within_bp_tol,
            "coord_correct": coord_correct,
            "exact_match": exact,
            "pred_len": pred_len,
            "truth_len": truth_len,
            "pred_minus_truth_len": (pred_len - truth_len) if pred_len is not None and truth_len is not None else None,
            "delta_start": d_start,
            "delta_end": d_end,
            "failure_mode": classify_prediction(
                pred, truth, match_iou, overlap_truth, best_iou, iou_threshold, codon_ok, within_bp_tol
            ),
        })
    return pd.DataFrame(rows)


def classify_prediction(
    pred: Dict,
    truth: Optional[Dict],
    match_iou: float,
    overlap_truth: Optional[Dict],
    best_iou: float,
    iou_threshold: float,
    codon_ok: bool = True,
    within_bp_tol: bool = False,
) -> str:
    if pred.get("pred_start") is None or pred.get("pred_end") is None:
        return "Not lifted"
    # within_bp_tol: both boundaries within the absolute bp tolerance -> treated as the
    # right region even when IoU is low (short features, where IoU is over-sensitive).
    if truth and (match_iou >= iou_threshold or within_bp_tol):
        if pred.get("pred_start") == truth["start"] and pred.get("pred_end") == truth["end"]:
            return "Correct"
        # codon_ok already folds in the reference-driven exemption: it is False only when
        # this gene is codon-required AND the lifted CDS has a bad start/stop/frame. A
        # frameshift gene (e.g. ORF1b, exempt) is codon_ok=True -> "Boundary offset", not
        # "Invalid codon", even though it has no ATG start.
        if not codon_ok:
            # right region but broken CDS on a codon-required gene — not coord-correct
            return "Invalid codon"
        return "Boundary offset"
    if truth:
        return "Wrong coords"
    if overlap_truth and best_iou >= iou_threshold:
        return "Alias/name gap"
    if overlap_truth and best_iou > 0:
        return "Possible overlap"
    return "Not in truth"


def summarize_comparison(df: pd.DataFrame, group_cols: Sequence[str]) -> pd.DataFrame:
    return (
        df.groupby(list(group_cols), dropna=False)
        .agg(
            total=("pred_name", "size"),
            exact=("exact_match", "sum"),
            coord_correct=("coord_correct", "sum"),
            mean_iou=("iou", "mean"),
        )
        .reset_index()
        .assign(
            exact_pct=lambda d: (d["exact"] / d["total"] * 100).round(2),
            coord_pct=lambda d: (d["coord_correct"] / d["total"] * 100).round(2),
            mean_iou=lambda d: d["mean_iou"].round(4),
        )
    )


def summarize_evaluable(df: pd.DataFrame, group_cols: Sequence[str]) -> pd.DataFrame:
    """Accuracy over the EVALUABLE set only (the paper denominator).

    Keeps predictions whose gene name is present in the query truth
    (``name_match == True``) — i.e. genes in R ∩ query-truth. Predictions with no
    same-name truth (granularity / name-gap, e.g. an ORF1a lifted onto a record
    whose truth only has the merged ORF1ab) are dropped from the denominator and
    analysed separately by the failure-attribution notebook, instead of being
    counted as errors. Generic: no virus or gene name is referenced.
    """
    evaluable = df[df["name_match"].fillna(False)] if "name_match" in df.columns else df
    return summarize_comparison(evaluable, group_cols)


def ref_codon_checkable_genes(bundle: Dict) -> set:
    """Reference genes whose OWN CDS is codon-clean (valid start + stop + frame).

    Only these genes are held to the codon check when scoring predictions. Reference
    features that are partial / frameshift / overlapping (e.g. PRRSV ORF1b, annotated
    `<7687..12072`, which has no independent ATG start) fail this and are therefore
    exempt — an "invalid codon" on such a gene reflects its biology, not a tool error.
    Empty set for non-CDS reference feature types (mat_peptide), where codon validity
    does not apply. Fully data-driven from the reference; no gene name is hard-coded.
    """
    if bundle.get("feature_type") != "CDS":
        return set()
    from Bio.Seq import Seq
    from app.src.lifting.validator import validate_cds_boundaries
    seq = str(bundle["record"].seq)
    genes = set()
    for f in bundle["features"]:
        s, e, strand, name = f.get("start"), f.get("end"), f.get("strand"), f.get("name")
        if s is None or e is None or not name:
            continue
        sub = seq[s - 1:e]
        if strand in ("-", -1):
            sub = str(Seq(sub).reverse_complement())
        if validate_cds_boundaries(sub).get("valid"):
            genes.add(name)
    return genes


# --------------------------------------------------------------------------- #
# Truth-based per-gene accuracy (ported from the archived accuracy notebooks). #
# Denominator = records where the reference gene is present in the query truth #
# (R ∩ query-truth), aligned by ref_name. Catches "no prediction at all" cases #
# and surfaces name-gap predictions separately. Fully generic (no gene names). #
# --------------------------------------------------------------------------- #

def build_truth_presence(bundle: Dict, query_records: Sequence[SeqRecord],
                         virus_label: Optional[str] = None) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Return (ref_model, truth_presence).

    ref_model: one row per reference gene with ref_start/ref_end/ref_len.
    truth_presence: one row per (record, ref-gene) with has_truth_gene flag.
    """
    ref_features = bundle["features"]
    ref_model = pd.DataFrame(ref_features).copy()
    ref_model["ref_len"] = ref_model["end"] - ref_model["start"] + 1
    ref_model = (
        ref_model.rename(columns={"name": "ref_name", "start": "ref_start", "end": "ref_end"})
        [["ref_name", "ref_start", "ref_end", "ref_len", "strand"]]
        .sort_values(["ref_start", "ref_end"])
        .reset_index(drop=True)
    )
    ref_gene_names = ref_model["ref_name"].dropna().astype(str).tolist()
    rows = []
    for record in query_records:
        truth_features, _ = parse_validation_truth_features(
            record, bundle["alias_lookup"], bundle["feature_type"], ref_features, virus_label,
        )
        truth_names = {f["name"] for f in truth_features}
        for gene in ref_gene_names:
            rows.append({"record_id": record.id, "gene": gene, "has_truth_gene": gene in truth_names})
    return ref_model, pd.DataFrame(rows)


def attach_ref_lengths(pred_df: pd.DataFrame, ref_model: pd.DataFrame) -> pd.DataFrame:
    """Add ref_len + pred/truth-vs-ref length deltas by joining on ref_name."""
    if pred_df.empty:
        return pred_df
    ref_len = ref_model.set_index("ref_name")["ref_len"]
    out = pred_df.copy()
    out["ref_len"] = out["ref_name"].map(ref_len)
    out["pred_minus_ref_len"] = out["pred_len"] - out["ref_len"]
    out["truth_minus_ref_len"] = out["truth_len"] - out["ref_len"]
    return out


def summarize_per_gene(pred_df: pd.DataFrame, truth_presence: pd.DataFrame,
                       ref_model: pd.DataFrame, virus_name: str) -> pd.DataFrame:
    """Per-gene buckets on the truth-available denominator.

    total = records where the gene is in truth AND a prediction exists;
    exact / coord_only / failed are mutually exclusive; no_hit and
    extra_predictions_without_truth (gene predicted where truth lacks it) tracked separately.
    """
    rows = []
    for gene in ref_model["ref_name"].tolist():
        truth_records = set(
            truth_presence[(truth_presence["gene"] == gene) & truth_presence["has_truth_gene"]]["record_id"]
        )
        gene_rows = pred_df[(pred_df["ref_name"] == gene) & pred_df["record_id"].isin(truth_records)]
        total = len(gene_rows)
        exact = int(gene_rows["exact_match"].fillna(False).sum()) if total else 0
        coord_correct = int(gene_rows["coord_correct"].fillna(False).sum()) if total else 0
        coord_only = coord_correct - exact
        failed = total - coord_correct
        no_hit = int(gene_rows["status"].fillna("").eq("no_hit").sum()) if total and "status" in gene_rows else 0
        extra = int(len(pred_df[(pred_df["ref_name"] == gene) & ~pred_df["record_id"].isin(truth_records)]))
        rows.append({
            "virus": virus_name, "gene": gene, "total": total,
            "exact": exact, "coord_only": coord_only, "coord_correct": coord_correct,
            "failed": failed, "no_hit": no_hit, "extra_predictions_without_truth": extra,
            "exact_pct": round(exact / total * 100, 2) if total else None,
            "coord_only_pct": round(coord_only / total * 100, 2) if total else None,
            "failed_pct": round(failed / total * 100, 2) if total else None,
            "accuracy_pct": round((exact + coord_only) / total * 100, 2) if total else None,
        })
    return pd.DataFrame(rows)


def summarize_overall(per_gene_df: pd.DataFrame, group_col: str = "virus") -> pd.DataFrame:
    """Roll per-gene buckets up to per-virus + an ALL row."""
    def agg(df, label):
        total = int(df["total"].sum()); exact = int(df["exact"].sum())
        coord_only = int(df["coord_only"].sum()); coord_correct = int(df["coord_correct"].sum())
        failed = int(df["failed"].sum())
        return {
            group_col: label, "total": total, "exact": exact, "coord_only": coord_only,
            "coord_correct": coord_correct, "failed": failed,
            "no_hit": int(df["no_hit"].sum()) if "no_hit" in df else 0,
            "exact_pct": round(exact / total * 100, 2) if total else None,
            "coord_only_pct": round(coord_only / total * 100, 2) if total else None,
            "failed_pct": round(failed / total * 100, 2) if total else None,
            "accuracy_pct": round((exact + coord_only) / total * 100, 2) if total else None,
        }
    rows = [agg(df, label) for label, df in per_gene_df.groupby(group_col, sort=False)]
    rows.append(agg(per_gene_df, "ALL"))
    return pd.DataFrame(rows)


def adjudicate_failures(df: pd.DataFrame, short_peptide_bp: int = 150) -> pd.DataFrame:
    """First-pass, GENERIC adjudication of non-coord-correct cases.

    Splits each failure into blame = ref_truth / tool / review using only ref/pred/truth
    length deltas (no virus- or gene-specific rule). Rules (ported from the archived
    breakdown notebooks, de-hardcoded):

      - no hit / not lifted                          -> tool  (no_hit)
      - truth feature absent                         -> ref_truth (truth_feature_absent)
      - pred length == ref, truth differs            -> ref_truth (ref_query_boundary_convention)
            · if the reference gene is short (<= short_peptide_bp) the same convention gap
              is flagged short_peptide_iou_artifact (a few-bp offset drops IoU below cutoff)
      - truth length == ref, prediction differs      -> tool  (prediction_boundary_error)
      - pred shorter than BOTH ref & truth, low cov  -> tool  (truncation)
      - differs from both, otherwise                 -> review

    Needs columns: pred_len, truth_len, ref_len, pred_minus_ref_len, truth_minus_ref_len,
    pred_minus_truth_len, coverage, status, failure_mode. Returns df + blame/cause/explanation.
    """
    def _one(row):
        status = str(row.get("status") or "")
        mode = row.get("failure_mode")
        pr, tr, pt = row.get("pred_minus_ref_len"), row.get("truth_minus_ref_len"), row.get("pred_minus_truth_len")
        ref_len, cov = row.get("ref_len"), row.get("coverage")
        if status == "no_hit" or mode == "Not lifted":
            return ("tool", "no_hit", "Feature in truth but tblastn produced no usable hit.")
        if pd.isna(row.get("truth_len")):
            return ("ref_truth", "truth_feature_absent", "Query truth lacks the same-name feature.")
        if pd.notna(pr) and pd.notna(tr):
            if pr == 0 and tr != 0:
                if pd.notna(ref_len) and ref_len <= short_peptide_bp:
                    return ("ref_truth", "short_peptide_iou_artifact",
                            "Prediction follows the reference; query truth differs by a few bp, but the gene "
                            "is short so IoU drops below the cutoff and it registers as failed.")
                return ("ref_truth", "ref_query_boundary_convention",
                        "Prediction length follows the reference; query truth uses a different boundary.")
            if tr == 0 and pr != 0:
                return ("tool", "prediction_boundary_error",
                        "Truth length follows the reference; the prediction boundary differs -> tool-side.")
            if pr != 0 and tr != 0:
                if pd.notna(pt) and pt < 0 and (cov if pd.notna(cov) else 1) < 0.9:
                    return ("tool", "truncation",
                            "Prediction shorter than BOTH ref and truth with low coverage -> tool truncation.")
                return ("review", "differs_from_ref_and_truth",
                        "Prediction differs from both ref and truth; manual review.")
        return ("review", "unclassified", "Not enough length evidence; manual review.")

    if df.empty:
        return df.assign(blame=[], cause=[], explanation=[])
    df = df.reset_index(drop=True)   # align index before concat, else NaN/duplicate rows
    adj = df.apply(lambda r: pd.Series(_one(r), index=["blame", "cause", "explanation"]), axis=1)
    return pd.concat([df, adj], axis=1)


def run_lifting_accuracy(virus_label: str, ref_path: Path, query_path: Path, output_dir: Path,
                         registry_path: Path = CONFIG / "virus_alias_registry.json",
                         progress: bool = True):
    """One-call lifting-accuracy: run tblastn, build truth presence, per-gene + overall buckets.

    Returns (per_pred, per_gene, overall, truth_presence, ref_model). per_pred carries the
    diagnostic columns (pred_len, delta_start/end, pred_minus_ref_len, ...) for inspection.
    """
    per_pred, _ = run_tblastn_against_truth(virus_label, ref_path, query_path, output_dir,
                                            registry_path=registry_path, progress=progress)
    bundle = load_reference_bundle(ref_path, registry_path)
    records = load_genbank_records(query_path)
    ref_model, truth_presence = build_truth_presence(bundle, records, virus_label)
    per_pred = attach_ref_lengths(per_pred, ref_model)
    per_gene = summarize_per_gene(per_pred, truth_presence, ref_model, virus_label)
    overall = summarize_overall(per_gene)
    return per_pred, per_gene, overall, truth_presence, ref_model


def run_tblastn_against_truth(
    virus_label: str,
    ref_path: Path,
    query_path: Path,
    output_dir: Path,
    registry_path: Path = CONFIG / "virus_alias_registry.json",
    progress: bool = True,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    from tqdm.auto import tqdm

    bundle = load_reference_bundle(ref_path, registry_path)
    query_records = load_genbank_records(query_path)
    codon_required = ref_codon_checkable_genes(bundle)
    rows = []

    iterator = tqdm(
        query_records,
        desc=f"tblastn  {virus_label}",
        unit="rec",
        dynamic_ncols=True,
        leave=True,
    ) if progress else query_records

    for record in iterator:
        truth, _ = parse_validation_truth_features(
            record,
            bundle["alias_lookup"],
            bundle["feature_type"],
            bundle["features"],
            virus_label,
        )
        if not truth:
            continue
        if progress:
            iterator.set_postfix_str(record.id, refresh=True)
        # Force the tblastn engine on every record (this measures lifting quality,
        # not routing) but with the tool's default config, not incidental defaults.
        cfg = PipelineConfig()
        lifted = lift_all_tblastn(
            ref_features=bundle["features"],
            ref_record=bundle["record"],
            query_record=record,
            min_coverage=cfg.min_coverage,
            min_identity=cfg.min_identity,
            evalue=cfg.evalue,
            rescue_window=cfg.rescue_window,
            validate_codons=(bundle["feature_type"] == "CDS"),
        )
        preds = lifted_to_rows(record.id, lifted, "tblastn")
        compared = compare_predictions_to_truth(preds, truth, codon_required_names=codon_required)
        compared.insert(0, "virus", virus_label)
        rows.append(compared)

    result = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
    summary = summarize_comparison(result, ["virus", "method"]) if not result.empty else pd.DataFrame()
    output_dir.mkdir(parents=True, exist_ok=True)
    result.to_csv(output_dir / "tblastn_vs_truth.tsv", sep="\t", index=False)
    summary.to_csv(output_dir / "tblastn_summary.tsv", sep="\t", index=False)
    return result, summary


def run_production_pipeline_against_truth(
    virus_label: str,
    ref_path: Path,
    query_path: Path,
    output_dir: Path,
    registry_path: Path = CONFIG / "virus_alias_registry.json",
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    bundle = load_reference_bundle(ref_path, registry_path)
    query_records = load_genbank_records(query_path)

    # Run the EXACT shipping pipeline (same routing + same default config the CLI
    # and UI use) so validation measures what users actually run — no re-implemented
    # routing, no magic parameters. run_pipeline handles direct vs tblastn internally.
    run_result = run_pipeline(
        ref_record=bundle["record"],
        query_records=query_records,
        ref_features=bundle["features"],
        ref_feature_type=bundle["feature_type"],
        alias_lookup=bundle["alias_lookup"],
        config=PipelineConfig(),
    )
    lifted_by_record = dict(run_result.all_results)
    codon_required = ref_codon_checkable_genes(bundle)

    rows = []
    for record in query_records:
        truth_type = select_feature_type(record, bundle["alias_lookup"]) or bundle["feature_type"]
        truth, _ = parse_validation_truth_features(
            record,
            bundle["alias_lookup"],
            truth_type,
            bundle["features"],
            virus_label,
        )
        if not truth:
            continue
        # method label only — same deterministic call the pipeline made internally
        strategy, _ = get_strategy(record, bundle["alias_lookup"])
        lifted = lifted_by_record.get(record.id, [])
        preds = lifted_to_rows(record.id, lifted, strategy)
        compared = compare_predictions_to_truth(preds, truth, codon_required_names=codon_required)
        compared.insert(0, "virus", virus_label)
        rows.append(compared)
    result = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
    summary = summarize_comparison(result, ["virus", "method"]) if not result.empty else pd.DataFrame()
    output_dir.mkdir(parents=True, exist_ok=True)
    result.to_csv(output_dir / "e2e_per_prediction.tsv", sep="\t", index=False)
    summary.to_csv(output_dir / "e2e_summary.tsv", sep="\t", index=False)
    return result, summary


def collect_alias_rows(records: Sequence[SeqRecord], feature_types: Sequence[str] = ("CDS", "mat_peptide")) -> pd.DataFrame:
    rows = []
    for record in records:
        for feature_type in feature_types:
            for feature in parse_features_for_type(record, feature_type):
                for field in _LOOKUP_QUALIFIER_KEYS:
                    value = feature.get(field)
                    if value:
                        rows.append({
                            "record_id": record.id,
                            "feature_type": feature_type,
                            "field": field,
                            "raw_value": value,
                        })
    return pd.DataFrame(rows)


def evaluate_alias_rows(rows: pd.DataFrame, alias_lookup: Dict[str, str], virus_label: str) -> pd.DataFrame:
    if rows.empty:
        return rows

    def resolve(value: str) -> str:
        hit = lookup_field_value(value, alias_lookup)
        if hit in EXCLUSION_SENTINELS:
            return "excluded"
        if hit is None:
            return "unresolved"
        return hit

    result = rows.copy()
    result.insert(0, "virus", virus_label)
    result["resolved_to"] = result["raw_value"].map(resolve)
    result["is_canonical"] = ~result["resolved_to"].isin(["excluded", "unresolved"])
    return result
