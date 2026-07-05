import re
from typing import Dict

from app.src.alias.gene_alias import normalize_text


STRONG_FIELDS = {"gene", "label", "standard_name"}
WEAK_FIELDS = {"product", "note", "locus_tag"}

GENERIC_NAME_BLACKLIST = {
    "gene",
    "protein",
    "viral protein",
    "unknown",
    "unknown protein",
    "hypothetical protein",
    "putative protein",
    "open reading frame",
    "orf",
}

DESCRIPTIVE_REVIEW_TERMS = {
    "polymerase",
    "protease",
    "helicase",
    "nucleocapsid",
    "capsid",
    "matrix",
    "envelope",
    "glycoprotein",
    "membrane",
    "replicase",
    "attachment",
    "polyprotein",
    "structural",
    "nonstructural",
}

DESCRIPTIVE_CANONICAL_ALIASES = {
    "envelope protein": {"E"},
    "envelope glycoprotein": {"E"},
    "membrane protein": {"M"},
    "nucleocapsid protein": {"N"},
    "spike protein": {"S"},
    "spike glycoprotein": {"S"},
}


def classify_alias_candidate(
    raw_value: str,
    field: str,
    canonical_name: str,
    evidence: Dict,
) -> Dict:
    """
    Classify one raw qualifier string as save_alias/manual_review/ignore.

    Coordinate evidence supports the feature-level canonical assignment, but the
    raw string is scored independently to avoid saving generic names globally.
    """
    raw_norm = normalize_text(raw_value)
    canonical_norm = normalize_text(canonical_name)
    raw_words = _word_tokens(raw_value)

    score = 0
    reasons = []
    has_name_specific_evidence = False

    iou = float(evidence.get("iou") or 0.0)
    if iou >= 0.95:
        score += 5
        reasons.append(f"coordinate IoU {iou:.2f}")
    elif iou >= 0.90:
        score += 4
        reasons.append(f"coordinate IoU {iou:.2f}")

    if evidence.get("strand_match"):
        score += 1
        reasons.append("same strand")

    if field in STRONG_FIELDS:
        score += 3
        reasons.append(f"strong field: {field}")
    elif field in WEAK_FIELDS:
        score -= 1
        reasons.append(f"weak/descriptive field: {field}")

    if raw_norm == canonical_norm:
        score += 5
        has_name_specific_evidence = True
        reasons.append("exact canonical text")
    elif canonical_norm and _raw_contains_canonical(raw_value, raw_norm, canonical_norm):
        score += 4
        has_name_specific_evidence = True
        reasons.append("contains canonical name")
    elif _looks_like_short_gene_symbol(raw_value, canonical_name):
        score += 3
        has_name_specific_evidence = True
        reasons.append("short gene symbol consistent with canonical")

    descriptive_match = _descriptive_alias_matches_canonical(raw_value, canonical_name)
    if descriptive_match:
        score += 5
        has_name_specific_evidence = True
        reasons.append("descriptive synonym matches canonical")

    if _is_generic_name(raw_value) and not descriptive_match:
        score -= 8
        reasons.append("generic name")
    elif raw_words & DESCRIPTIVE_REVIEW_TERMS:
        if has_name_specific_evidence:
            score += 1
            reasons.append("descriptive term with specific gene name")
        else:
            score -= 4
            reasons.append("descriptive biological term")

    if _looks_like_locus_tag(raw_value):
        score -= 6
        reasons.append("locus tag-like value")

    if score >= 8 and not has_name_specific_evidence:
        score = 5
        reasons.append("no name-specific alias evidence")

    if score >= 8:
        action = "save_alias"
        confidence = "high"
        default_save = True
    elif score >= 3:
        action = "manual_review"
        confidence = "medium"
        default_save = False
    else:
        action = "ignore"
        confidence = "low"
        default_save = False

    return {
        "suggested_action": action,
        "confidence": confidence,
        "score": score,
        "reason": "; ".join(reasons) if reasons else "no strong alias evidence",
        "default_save": default_save,
    }


def _word_tokens(value: str) -> set:
    return set(re.findall(r"[a-zA-Z]+", (value or "").lower()))


def _is_generic_name(value: str) -> bool:
    raw = value or ""
    normalized = normalize_text(raw)
    normalized_blacklist = {normalize_text(v) for v in GENERIC_NAME_BLACKLIST}
    if normalized in normalized_blacklist:
        return True

    # Do not demote specific ORF labels like ORF3/ORF 3/ORF1a to generic "orf".
    if re.search(r"\d", raw):
        return False

    normalized_words = " ".join(re.findall(r"[a-zA-Z]+", raw.lower())).strip()
    return normalized_words in GENERIC_NAME_BLACKLIST


def _raw_contains_canonical(raw_value: str, raw_norm: str, canonical_norm: str) -> bool:
    if len(canonical_norm) > 1:
        return canonical_norm in raw_norm
    tokens = {normalize_text(token) for token in re.findall(r"[A-Za-z0-9]+", raw_value or "")}
    return canonical_norm in tokens


def _looks_like_short_gene_symbol(raw_value: str, canonical_name: str) -> bool:
    raw = (raw_value or "").strip()
    canonical = (canonical_name or "").strip()
    if not re.search(r"\d", raw):
        return False
    canonical_digits = re.findall(r"\d+[a-zA-Z]*", canonical.lower())
    if not canonical_digits:
        return False

    raw_tokens = re.findall(r"[a-zA-Z]+[0-9]+[a-zA-Z]*", raw.lower())
    if len(raw) <= 12:
        raw_tokens.append(raw.lower())

    for token in raw_tokens:
        token_digits = re.findall(r"\d+[a-zA-Z]*", token)
        if token_digits and token_digits[-1] == canonical_digits[-1]:
            return True
    return False


def _descriptive_alias_matches_canonical(raw_value: str, canonical_name: str) -> bool:
    normalized = normalize_text(raw_value or "")
    canonical_norm = normalize_text(canonical_name or "")
    for raw_alias, canonicals in DESCRIPTIVE_CANONICAL_ALIASES.items():
        if normalized == normalize_text(raw_alias):
            return canonical_norm in {normalize_text(canonical) for canonical in canonicals}
    return False


def _looks_like_locus_tag(value: str) -> bool:
    raw = (value or "").strip()
    return bool(re.fullmatch(r"[A-Z]{2,}[_-]?[0-9]{3,}[A-Z0-9_.-]*", raw))
