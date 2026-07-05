import hashlib
import json
import re
from typing import Dict, Iterable, List, Optional, Tuple

from app.src.alias.alias_classifier import DESCRIPTIVE_REVIEW_TERMS, GENERIC_NAME_BLACKLIST
from app.src.alias.alias_payload import build_uncertain_suggestion_review_payload
from app.src.alias.gene_alias import normalize_text
from app.src.llm.config import LLMConfig
from app.src.llm.provider import LLMProviderError, NoopLLMProvider, OpenAILLMProvider


LLM_REVIEW_ACTIONS = {"save_alias", "ignore", "skip", "move_to_ambiguous"}
LLM_CONFIDENCES = {"low", "medium", "high"}


def needs_llm_review(row: Dict) -> bool:
    """
    Return True only for alias suggestions where deterministic scoring is not
    enough or the name has a known risky shape.
    """
    action = row.get("suggested_action")
    confidence = row.get("confidence")
    raw_value = str(row.get("raw_value") or "")

    if action == "manual_review" or confidence == "medium":
        return True

    if action == "ignore" and _looks_specific_but_low_confidence(raw_value):
        return True

    if _has_direct_canonical_signal(raw_value, row.get("canonical_name")):
        return False

    if action == "ignore" and confidence == "low":
        return _should_review_low_confidence_ignore(row)

    if confidence == "high" and _has_risky_shape(raw_value, row.get("canonical_name")):
        return True

    return False


def review_uncertain_alias_suggestions(
    suggestions: List[Dict],
    virus_name: str,
    canonical_names: Iterable[str],
    *,
    ignored_names: Optional[Iterable[str]] = None,
    ambiguous_names: Optional[Iterable[str]] = None,
    config: Optional[LLMConfig] = None,
    provider=None,
    cache: Optional[Dict[str, Dict]] = None,
) -> Tuple[List[Dict], Dict]:
    """
    Run optional LLM review for uncertain suggestion rows and merge advisory
    results back into the suggestions list.
    """
    config = config or LLMConfig.from_env()
    enriched = [_with_default_llm_fields(row) for row in suggestions]
    uncertain = [
        row for row in enriched
        if needs_llm_review(row)
    ][:config.max_rows]

    diagnostics = {
        "enabled": config.enabled,
        "available": config.available if provider is None else True,
        "model": config.model,
        "fallback_model": config.fallback_model,
        "uncertain_rows": len([row for row in enriched if needs_llm_review(row)]),
        "submitted_rows": len(uncertain),
        "reviewed_rows": 0,
        "status": "skipped",
        "error": None,
        "cache_hit": False,
    }

    if not uncertain:
        diagnostics["status"] = "no_uncertain_rows"
        return enriched, diagnostics

    if provider is None and not config.available:
        diagnostics["status"] = "missing_api_key" if config.enabled else "disabled"
        return enriched, diagnostics

    payload = build_uncertain_suggestion_review_payload(
        virus_name=virus_name,
        canonical_names=list(canonical_names),
        suggestions=uncertain,
        ignored_names=list(ignored_names or []),
        ambiguous_names=list(ambiguous_names or []),
    )
    cache_key = alias_review_cache_key(payload)
    if cache is not None and cache_key in cache:
        response = cache[cache_key]
        diagnostics["cache_hit"] = True
    else:
        provider = provider or OpenAILLMProvider(config)
        try:
            response = provider.review_alias_suggestions(payload)
        except LLMProviderError as exc:
            diagnostics["status"] = "error"
            diagnostics["error"] = str(exc)
            return enriched, diagnostics
        if cache is not None:
            cache[cache_key] = response

    reviews = _valid_reviews(response.get("reviews", []), set(canonical_names))
    merged = merge_alias_reviews(enriched, reviews)
    diagnostics["reviewed_rows"] = sum(1 for row in merged if row.get("llm_reviewed"))
    diagnostics["status"] = "reviewed"
    return merged, diagnostics


def merge_alias_reviews(suggestions: List[Dict], reviews: List[Dict]) -> List[Dict]:
    by_id = {review.get("review_id"): review for review in reviews}
    merged = []
    for row in suggestions:
        review = by_id.get(row.get("llm_review_id"))
        if not review:
            merged.append(row)
            continue
        merged.append({
            **row,
            "llm_reviewed": True,
            "llm_action": review.get("recommendation"),
            "llm_canonical_name": review.get("canonical_name") or "",
            "llm_confidence": review.get("confidence"),
            "llm_reason": review.get("reason"),
        })
    return merged


def alias_review_cache_key(payload: Dict) -> str:
    raw = json.dumps(payload, sort_keys=True, ensure_ascii=False)
    return hashlib.sha256(raw.encode("utf-8")).hexdigest()


def _with_default_llm_fields(row: Dict) -> Dict:
    raw = row.get("raw_value") or ""
    field = row.get("field") or ""
    canonical = row.get("canonical_name") or ""
    record = row.get("record_id") or ""
    digest = hashlib.sha1(
        "|".join([normalize_text(raw), field, canonical, str(record)]).encode("utf-8")
    ).hexdigest()[:12]
    return {
        **row,
        "llm_review_id": row.get("llm_review_id") or f"alias_{digest}",
        "llm_reviewed": bool(row.get("llm_reviewed", False)),
        "llm_action": row.get("llm_action", ""),
        "llm_canonical_name": row.get("llm_canonical_name", ""),
        "llm_confidence": row.get("llm_confidence", ""),
        "llm_reason": row.get("llm_reason", ""),
    }


def _valid_reviews(reviews: Iterable[Dict], canonical_names: set) -> List[Dict]:
    valid = []
    for review in reviews or []:
        if review.get("recommendation") not in LLM_REVIEW_ACTIONS:
            continue
        if review.get("confidence") not in LLM_CONFIDENCES:
            continue
        canonical = review.get("canonical_name")
        if review.get("recommendation") == "save_alias" and canonical not in canonical_names:
            continue
        valid.append(review)
    return valid


def _has_risky_shape(raw_value: str, canonical_name: Optional[str]) -> bool:
    raw = raw_value or ""
    if any(sep in raw for sep in [";", "/", ","]):
        return True
    if re.search(r"\bor\b", raw, flags=re.IGNORECASE):
        return True

    raw_norm = normalize_text(raw)
    canonical_norm = normalize_text(canonical_name or "")
    if _orf_family_root(raw_norm) and _orf_family_root(raw_norm) == _orf_family_root(canonical_norm):
        return raw_norm != canonical_norm

    words = set(re.findall(r"[a-zA-Z]+", raw.lower()))
    return bool(words & DESCRIPTIVE_REVIEW_TERMS and raw_norm != canonical_norm)


def _has_direct_canonical_signal(raw_value: str, canonical_name: Optional[str]) -> bool:
    """
    Detect raw aliases that clearly name the current canonical.

    ORF-family names are inherently risky, but rows like "ORF1a protein" for
    canonical ORF1a should not be escalated just because ORF1ab also exists.
    """
    raw_norm = normalize_text(raw_value or "")
    canonical_norm = normalize_text(canonical_name or "")
    if not raw_norm or not canonical_norm:
        return False
    if raw_norm == canonical_norm:
        return True
    suffixes = (
        "protein",
        "polyprotein",
        "gene",
        "cds",
        "openreadingframe",
    )
    return any(raw_norm == f"{canonical_norm}{suffix}" for suffix in suffixes)


def _looks_specific_but_low_confidence(raw_value: str) -> bool:
    raw = raw_value or ""
    raw_norm = normalize_text(raw)
    if not raw_norm or raw_norm in {normalize_text(v) for v in GENERIC_NAME_BLACKLIST}:
        return False
    if re.search(r"\borf\s*\d+[a-zA-Z]?\b", raw, flags=re.IGNORECASE):
        return True
    if re.fullmatch(r"[A-Za-z]{1,4}\d+[A-Za-z]?", raw.strip()):
        return True
    if _has_risky_shape(raw, None):
        return True
    return False


def _should_review_low_confidence_ignore(row: Dict) -> bool:
    raw = str(row.get("raw_value") or "")
    raw_norm = normalize_text(raw)
    hard_ignored = {
        normalize_text(value)
        for value in {
            "unknown",
            "unknown protein",
            "hypothetical protein",
            "putative protein",
            "protein",
            "viral protein",
            "gene",
            "orf",
            "open reading frame",
        }
    }
    if raw_norm in hard_ignored:
        return False

    support = int(row.get("support_count") or 0)
    iou = float(row.get("iou") or 0.0)
    has_coordinate_support = support >= 2 or iou >= 0.90
    words = set(re.findall(r"[a-zA-Z]+", raw.lower()))
    return bool(has_coordinate_support and (words & DESCRIPTIVE_REVIEW_TERMS))


def _orf_family_root(value: str) -> Optional[str]:
    match = re.search(r"\borf\s*(\d+)", value or "", flags=re.IGNORECASE)
    return f"orf{match.group(1)}" if match else None
