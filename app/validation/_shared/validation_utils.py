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
    AMBIGUOUS_SENTINEL,
    IGNORED_SENTINEL,
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
    target_names = sorted({feature["name"] for feature in ref_features if feature.get("name")})
    extra_names: List[str] = []
    if virus_label and virus_label.upper().startswith("PRRS"):
        # PRRSV records may annotate ORF1a/ORF1b as one ORF1ab feature.
        # Keep it in truth so validation can identify granularity mismatches
        # instead of treating the region as absent.
        extra_names.append("ORF1ab")
    return target_names, extra_names


def should_use_target_truth_filter(feature_type: Optional[str], virus_label: Optional[str] = None) -> bool:
    if virus_label and virus_label.upper().startswith("PRRS"):
        return True
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
        })
    return rows


def best_overlap(pred_row: Dict, truth_features: Sequence[Dict]) -> Tuple[Optional[Dict], float]:
    best_feature, best_iou = None, 0.0
    for truth in truth_features:
        score = iou(pred_row.get("pred_start"), pred_row.get("pred_end"), truth["start"], truth["end"])
        if score > best_iou:
            best_feature, best_iou = truth, score
    return best_feature, best_iou


def compare_predictions_to_truth(
    predictions: Sequence[Dict],
    truth_features: Sequence[Dict],
    iou_threshold: float = 0.90,
) -> pd.DataFrame:
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
        coord_correct = truth is not None and match_iou >= iou_threshold
        exact = (
            truth is not None
            and pred.get("pred_start") == truth["start"]
            and pred.get("pred_end") == truth["end"]
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
            "coord_correct": coord_correct,
            "exact_match": exact,
            "failure_mode": classify_prediction(pred, truth, match_iou, overlap_truth, best_iou, iou_threshold),
        })
    return pd.DataFrame(rows)


def classify_prediction(
    pred: Dict,
    truth: Optional[Dict],
    match_iou: float,
    overlap_truth: Optional[Dict],
    best_iou: float,
    iou_threshold: float,
) -> str:
    if pred.get("pred_start") is None or pred.get("pred_end") is None:
        return "Not lifted"
    if truth and match_iou >= iou_threshold:
        if pred.get("pred_start") == truth["start"] and pred.get("pred_end") == truth["end"]:
            return "Correct"
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
        lifted = lift_all_tblastn(
            ref_features=bundle["features"],
            ref_record=bundle["record"],
            query_record=record,
            validate_codons=(bundle["feature_type"] == "CDS"),
        )
        preds = lifted_to_rows(record.id, lifted, "tblastn")
        compared = compare_predictions_to_truth(preds, truth)
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
        strategy, query_feature_type = get_strategy(record, bundle["alias_lookup"])
        if strategy == "direct":
            lifted = direct_extract_with_alias(record, query_feature_type, bundle["features"], bundle["alias_lookup"])
        else:
            lifted = process_one_query_record(
                ref_record=bundle["record"],
                query_record=record,
                ref_cds=bundle["features"],
                ref_feature_type=bundle["feature_type"],
                min_coverage=0.5,
                min_identity=0.3,
                evalue=1e-5,
                rescue_window=50,
            )
        preds = lifted_to_rows(record.id, lifted, strategy)
        compared = compare_predictions_to_truth(preds, truth)
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
        if hit == IGNORED_SENTINEL:
            return "ignored"
        if hit == AMBIGUOUS_SENTINEL:
            return "ambiguous"
        if hit is None:
            return "unresolved"
        return hit

    result = rows.copy()
    result.insert(0, "virus", virus_label)
    result["resolved_to"] = result["raw_value"].map(resolve)
    result["is_canonical"] = ~result["resolved_to"].isin(["ignored", "ambiguous", "unresolved"])
    return result
