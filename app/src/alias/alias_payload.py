import hashlib
import json
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from Bio.SeqRecord import SeqRecord

from app.src.alias.gene_alias import normalize_text


"""
Module: alias_payload.py

Purpose:
    Extract structured feature information from GenBank records
    and build JSON payloads intended for LLM-assisted alias building.

Two payload types:
    - "map_aliases":     virus is known, some feature names are unresolved
    - "build_alias_map": virus is new, need to build alias map from scratch

Main entry point:
    run_alias_pipeline() — processes ref + query records, returns per-record
    results and an aggregated LLM payload if anything is unresolved.

Notes:
    - Payloads are compact and intended for optional review-only LLM calls.
    - Nucleotide sequences are never included.
    - Only naming-relevant fields are extracted per feature.
"""


# ---------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------

def _extract_feature_info(feature: Dict) -> Dict:
    return {
        "name": feature.get("name"),
        "gene": feature.get("gene"),
        "product": feature.get("product"),
        "length": feature.get("length"),
        "order": feature.get("order"),
    }


def _get_organism(record: SeqRecord) -> str:
    return record.annotations.get("organism", "") or record.description or record.id


def _unique_ordered(values: List[str]) -> List[str]:
    seen = set()
    result = []
    for v in values:
        if v not in seen:
            seen.add(v)
            result.append(v)
    return result


# ---------------------------------------------------------------------
# Payload builders
# ---------------------------------------------------------------------

def build_map_aliases_payload(
    record: SeqRecord,
    unresolved_features: List[Dict],
    existing_canonical_names: List[str],
    feature_type: str,
) -> Dict:
    """
    Build payload for virus known, some feature names not yet in alias map.
    """
    return {
        "task": "map_aliases",
        "record_id": record.id,
        "organism": _get_organism(record),
        "feature_type": feature_type,
        "existing_canonical_names": existing_canonical_names,
        "unresolved_features": [_extract_feature_info(f) for f in unresolved_features],
    }


def build_new_virus_payload(
    record: SeqRecord,
    all_features: List[Dict],
    feature_type: str,
) -> Dict:
    """
    Build payload for virus not yet in registry. LLM builds alias map from scratch.
    """
    return {
        "task": "build_alias_map",
        "record_id": record.id,
        "organism": _get_organism(record),
        "feature_type": feature_type,
        "all_features": [_extract_feature_info(f) for f in all_features],
    }


def build_uncertain_suggestion_review_payload(
    virus_name: str,
    canonical_names: List[str],
    suggestions: List[Dict],
    excluded_names: Optional[List[str]] = None,
) -> Dict:
    """
    Build a compact payload for LLM review of uncertain alias suggestions.

    This payload is intentionally row-level and small: it contains only the
    deterministic suggestion context needed to advise the user. It never
    includes nucleotide sequence or full GenBank records.
    """
    return {
        "task": "review_uncertain_alias_suggestions",
        "virus": virus_name,
        "available_canonicals": canonical_names,
        "excluded_names": excluded_names or [],
        "instructions": (
            "Each row is a raw string taken from the gene/product qualifier of an "
            "annotated gene feature in a real record: it is the label some lab used "
            "for a gene at those coordinates. Your job is to normalize these "
            "inconsistent labels onto available_canonicals. "
            # --- the central point: descriptive labels ARE the aliases --------
            "CRUCIAL. Labs routinely name genes with plain descriptive phrases — a "
            "membrane protein, a nucleocapsid protein, a glycoprotein, a matrix "
            "protein, an envelope protein, a bare gene letter. These descriptive "
            "labels ARE exactly the aliases this task exists to capture; they are NOT "
            "'descriptive context' to be discarded. Sounding descriptive, functional, "
            "or plain is NEVER by itself a reason to ignore or skip a row. If you drop "
            "'membrane protein' or 'nucleocapsid protein' or a lone gene letter, the "
            "records that use that label can never be resolved — which defeats the "
            "whole purpose. "
            # --- scope + decision rule, by evidence not by surface form -------
            "SCOPE. This alias map applies ONLY to the virus in the 'virus' field; "
            "every alias is looked up solely within this virus's available_canonicals. "
            "So the test is not whether a label is unambiguous in general biology, but "
            "whether it is unambiguous WITHIN this virus, decided by the evidence "
            "fields. DECISION RULE: recommend save_alias to canonical_candidate "
            "whenever canonical_candidate is set, iou is reasonably high, and "
            "cross_canonical_target_count is 1 (no other available_canonical competes "
            "for it) — regardless of how short, plain, abbreviated, or descriptive the "
            "string is. If matching_available_canonical is present, prefer it. "
            "query_name is weak context and may describe a parent or combined ORF; do "
            "not reject a raw label solely because query_name differs. "
            # --- what genuinely warrants ignore/skip --------------------------
            "Use ignore ONLY for: (a) a label that resolves to more than one "
            "available_canonical by coordinate — shared_across_canonicals true or "
            "cross_canonical_target_count greater than 1 — which is genuinely "
            "ambiguous within this virus; or (b) a string that is not a gene label at "
            "all, such as a free-text comment, a mutation or assembly remark, or a "
            "purely generic word with no coordinate resolution (for example "
            "'polyprotein' or 'unknown protein' alone with no canonical_candidate). "
            "Use skip only when there is no canonical_candidate and no available "
            "canonical fits. Do NOT ignore or skip a label merely for being "
            "descriptive, generic-sounding, short, or lacking an ORF number when the "
            "coordinates already tie it to exactly one canonical."
        ),
        "suggestions": [
            _extract_suggestion_review_info(row, canonical_names)
            for row in suggestions
        ],
    }


def build_unresolved_name_review_payload(
    virus_name: str,
    canonical_names: List[str],
    unknown_items: Dict[str, Dict],
    ambiguous_items: Optional[Dict[str, Dict]] = None,
    excluded_names: Optional[List[str]] = None,
) -> Dict:
    """
    Build a payload for LLM review of names that did not get coordinate-backed
    alias suggestions.

    This is used on the manual Name review page. It is intentionally advisory:
    the model can suggest a mapping, but the UI still requires user approval.
    """
    unresolved_rows = []
    for representative, info in (unknown_items or {}).items():
        unresolved_rows.append(
            _extract_unresolved_review_info(
                representative,
                info,
                canonical_names,
                is_ambiguous=False,
            )
        )
    for representative, info in (ambiguous_items or {}).items():
        unresolved_rows.append(
            _extract_unresolved_review_info(
                representative,
                info,
                canonical_names,
                is_ambiguous=True,
            )
        )

    return {
        "task": "review_unresolved_names",
        "virus": virus_name,
        "available_canonicals": canonical_names,
        "excluded_names": excluded_names or [],
        "instructions": (
            "Review unresolved query annotation names that were not resolved by the "
            "alias config and may not have coordinate-backed tblastn suggestions. "
            "Recommend save_alias only when the representative name or one of its "
            "candidate qualifiers clearly belongs to exactly one available canonical. "
            "Recommend ignore for broad descriptions, comments, locus-like values, "
            "or parent annotations that should not become a reusable alias. Recommend "
            "skip when the evidence is insufficient. Treat "
            "available_canonicals as authoritative and do not invent new canonicals. "
            "Use support_count and candidate_values_count as context: high support "
            "with one clear target can justify save_alias; high support with broad "
            "or mixed wording is stronger evidence for ignore/exclude "
            "than for skip. "
            "If a name explicitly denotes a span covering two or more adjacent "
            "canonicals (for example by joining their identifiers with a slash, a "
            "dash, or wording like 'contains X and Y'), prefer a single available "
            "canonical that represents that combined span when one exists; otherwise "
            "skip and let the user add a new canonical first. "
            # --- positional evidence ---------------------------------------
            "POSITIONAL EVIDENCE. Some rows carry coordinate context: feature_type, "
            "start, end, length_bp, neighbor_before, neighbor_after, is_subfeature_of, "
            "ref_slot_by_position and spans_reference_genes. When these are present they "
            "outrank the wording of the name itself, because annotation strings are "
            "submitted inconsistently while genomic position is not. Apply them as follows. "
            "(1) If is_subfeature_of is set, the feature lies wholly inside another coding "
            "feature, so it is a cleavage product or domain, not a gene: recommend ignore "
            "even when the name looks specific and biologically real, for example a "
            "polymerase or protease name sitting inside a replicase ORF. "
            "(2) If ref_slot_by_position names exactly one available canonical, or if "
            "neighbor_before and neighbor_after bracket exactly one missing slot in the "
            "reference gene order, prefer save_alias to that canonical even when the raw "
            "string is abbreviated, unfamiliar, or superficially suggests a different gene. "
            "(3) If spans_reference_genes lists several adjacent reference genes and a "
            "combined canonical covering them is available, prefer that combined canonical. "
            "(4) Treat length_bp as a sanity check against the expected size of the "
            "candidate canonical; a large mismatch is evidence against save_alias. "
            "(5) When positional evidence and the name disagree, say so in reason and "
            "follow the position. When positional evidence is absent, fall back to the "
            "name-based rules above and stay conservative."
        ),
        "suggestions": unresolved_rows,
    }


# ---------------------------------------------------------------------
# Feature helpers
# ---------------------------------------------------------------------

def get_unresolved_features(normalized_features: List[Dict]) -> List[Dict]:
    """Return features that were not resolved by alias lookup (name_source == 'raw')."""
    return [f for f in normalized_features if f.get("name_source") == "raw"]


_POSITION_QUALIFIER_KEYS = ("gene", "product", "note", "standard_name", "label")


def build_position_context(
    records,
    raw_value: str,
    resolve_name,
    ref_gene_order: Optional[List[str]] = None,
    coding_types: Tuple[str, ...] = ("CDS",),
) -> Dict:
    """
    Derive coordinate context for an unresolved annotation name.

    Args:
        records: iterable of SeqRecord to search for the raw name.
        raw_value: the unresolved annotation string.
        resolve_name: callable(str) -> Optional[str], maps a raw annotation
            value to a canonical name (or None). Used only to name the
            *neighbouring* features, never the row under review.
        ref_gene_order: canonical gene order of the reference, 5' to 3'.
            Used to turn a neighbour pair into a single missing slot.
        coding_types: feature types treated as genes for neighbour/containment
            logic. Anything outside this set can be a sub-feature.

    Returns a dict suitable for `info["position"]`. Empty when the name is not
    found. No ground-truth mapping of `raw_value` is consulted anywhere.
    """
    target = normalize_text(raw_value)
    for record in records or []:
        features = []
        for feature in record.features:
            if feature.type not in set(coding_types) | {"mat_peptide"}:
                continue
            values = [
                str(value)
                for key in _POSITION_QUALIFIER_KEYS
                for value in feature.qualifiers.get(key, [])
            ]
            features.append({
                "type": feature.type,
                "start": int(feature.location.start) + 1,
                "end": int(feature.location.end),
                "values": values,
                "canonical": next(
                    (
                        c for c in (resolve_name(v) for v in values)
                        if c and not (str(c).startswith("__") and str(c).endswith("__"))
                    ),
                    None,
                ),
                "is_match": any(normalize_text(v) == target for v in values),
            })

        match = next((f for f in features if f["is_match"]), None)
        if not match:
            continue

        coding = sorted(
            (f for f in features if f["type"] in coding_types and f is not match),
            key=lambda f: f["start"],
        )
        enclosing = next(
            (
                f for f in coding
                if f["start"] <= match["start"] and f["end"] >= match["end"]
                and (f["end"] - f["start"]) > (match["end"] - match["start"])
            ),
            None,
        )
        # Order by start rather than requiring disjoint ranges: overlapping
        # reading frames are normal (coronavirus ORF3/E/M routinely overlap by
        # a few dozen bp), and an end-based test silently drops the true
        # neighbour in exactly those cases.
        before = [f for f in coding if f["start"] < match["start"] and f["canonical"]]
        after = [f for f in coding if f["start"] > match["start"] and f["canonical"]]

        context = {
            "feature_type": match["type"],
            "start": match["start"],
            "end": match["end"],
            "genome_length": len(record.seq) if record.seq is not None else None,
            "neighbor_before": before[-1]["canonical"] if before else None,
            "neighbor_after": after[0]["canonical"] if after else None,
            "is_subfeature_of": (
                enclosing["canonical"] or "an unnamed coding feature"
            ) if enclosing else None,
        }
        candidates = _slot_candidates(
            context["neighbor_before"], context["neighbor_after"], ref_gene_order
        )
        context["slot_candidates"] = candidates
        # Only commit to a single slot when the neighbours leave exactly one
        # gap. A wider gap stays a candidate list so the model can weigh it
        # against length_bp instead of being handed a guess.
        context["ref_slot_by_position"] = candidates[0] if len(candidates) == 1 else None
        return context
    return {}


def _slot_candidates(
    before: Optional[str],
    after: Optional[str],
    ref_gene_order: Optional[List[str]],
) -> List[str]:
    """Reference genes lying strictly between two resolved neighbours."""
    if not ref_gene_order:
        return []
    try:
        i = ref_gene_order.index(before) if before else -1
        j = ref_gene_order.index(after) if after else len(ref_gene_order)
    except ValueError:
        return []
    return ref_gene_order[i + 1:j]


def _extract_suggestion_review_info(row: Dict, canonical_names: Optional[List[str]] = None) -> Dict:
    return {
        "review_id": row.get("llm_review_id"),
        "raw_value": row.get("raw_value"),
        "field": row.get("field"),
        "matching_available_canonical": _matching_available_canonical(
            row.get("raw_value"),
            canonical_names or [],
        ),
        "deterministic_action": row.get("suggested_action"),
        "deterministic_confidence": row.get("confidence"),
        "deterministic_reason": row.get("reason"),
        "deterministic_score": row.get("score"),
        "canonical_candidate": row.get("canonical_name"),
        "shared_across_canonicals": bool(row.get("shared_across_canonicals")),
        "cross_canonical_targets": row.get("cross_canonical_targets") or [],
        "cross_canonical_target_count": row.get("cross_canonical_target_count") or 0,
        "query_feature_type": row.get("query_feature_type"),
        "query_name": row.get("query_name"),
        "support_count": row.get("support_count"),
        "support_records": row.get("support_records"),
        "iou": row.get("iou"),
        "coverage": row.get("coverage"),
        "identity": row.get("identity"),
    }


def _extract_unresolved_review_info(
    representative: str,
    info: Dict,
    canonical_names: List[str],
    is_ambiguous: bool,
) -> Dict:
    candidates = [str(value) for value in info.get("candidates", []) if value]
    matching = (
        _matching_available_canonical(representative, canonical_names)
        or _first_matching_candidate(candidates, canonical_names)
    )
    return {
        "review_id": _unresolved_review_id(representative),
        "raw_value": representative,
        "candidate_values": candidates,
        "matching_available_canonical": matching,
        "deterministic_action": "ambiguous" if is_ambiguous else "unresolved",
        "deterministic_confidence": "low",
        "deterministic_reason": (
            "raw name maps to multiple canonicals"
            if is_ambiguous
            else "raw name was not resolved by alias lookup"
        ),
        "canonical_candidate": matching,
        "support_count": len(info.get("records", [])),
        "support_records": list(info.get("records", []))[:12],
        "candidate_values_count": len(candidates),
        "has_multiple_candidate_values": len(candidates) > 1,
        "high_support_unresolved": len(info.get("records", [])) >= 3,
        "is_ambiguous": is_ambiguous,
        # --- positional evidence -------------------------------------------
        # A bare annotation string is often undecidable ("sM" could read as
        # spike or membrane). Where the feature actually sits in the genome
        # usually is decidable. These fields are optional: they only appear
        # when the caller supplied coordinate context, so older callers keep
        # working unchanged.
        **_positional_evidence(info),
    }


def _positional_evidence(info: Dict) -> Dict:
    """
    Emit coordinate-derived context for an unresolved name, when available.

    Rationale: name-only review cannot distinguish a real gene alias from a
    mature-peptide sub-part or resolve names whose wording collides with a
    different canonical. Position can. Every field here is derived from the
    query record's own annotation layout plus the reference gene order --
    never from any ground-truth mapping of the name itself.
    """
    context = info.get("position") or {}
    if not context:
        return {}

    start = context.get("start")
    end = context.get("end")
    evidence = {
        "feature_type": context.get("feature_type"),
        "start": start,
        "end": end,
        "length_bp": (end - start + 1) if (start is not None and end is not None) else None,
        "genome_length_bp": context.get("genome_length"),
        # Canonical names of the nearest already-resolved features on each
        # side. With a known reference gene order this is often enough to
        # pin the slot exactly.
        "neighbor_before": context.get("neighbor_before"),
        "neighbor_after": context.get("neighbor_after"),
        # Set when this feature lies wholly inside another coding feature,
        # i.e. it is a cleavage product / sub-part rather than a gene.
        "is_subfeature_of": context.get("is_subfeature_of"),
        # Which reference gene slot the coordinates fall into, by position
        # only (no name matching involved).
        "ref_slot_by_position": context.get("ref_slot_by_position"),
        "slot_candidates": context.get("slot_candidates") or None,
        "spans_reference_genes": context.get("spans_reference_genes") or None,
    }
    return {key: value for key, value in evidence.items() if value is not None}


def _unresolved_review_id(representative: str) -> str:
    digest = hashlib.sha1(normalize_text(representative).encode("utf-8")).hexdigest()[:12]
    return f"unresolved_{digest}"


def _first_matching_candidate(candidates: List[str], canonical_names: List[str]) -> Optional[str]:
    for candidate in candidates:
        matched = _matching_available_canonical(candidate, canonical_names)
        if matched:
            return matched
    return None


def _matching_available_canonical(raw_value: str, canonical_names: List[str]) -> Optional[str]:
    raw_norm = normalize_text(raw_value or "")
    if not raw_norm:
        return None
    matches = []
    for canonical in canonical_names:
        canonical_norm = normalize_text(canonical)
        if not canonical_norm:
            continue
        if raw_norm == canonical_norm:
            matches.append((0, canonical))
        elif raw_norm in {
            f"{canonical_norm}protein",
            f"{canonical_norm}polyprotein",
            f"{canonical_norm}gene",
            f"{canonical_norm}cds",
        }:
            matches.append((1, canonical))
        elif _descriptive_orf_polyprotein_matches(raw_norm, canonical_norm):
            matches.append((2, canonical))
        elif _combined_orf_matches(raw_value or "", canonical_norm):
            matches.append((3, canonical))
    if not matches:
        return None
    return sorted(matches, key=lambda item: (item[0], len(item[1])))[0][1]


def _descriptive_orf_polyprotein_matches(raw_norm: str, canonical_norm: str) -> bool:
    if not raw_norm.startswith("polyprotein") and "polyprotein" not in raw_norm:
        return False

    canonical_match = re.fullmatch(r"orf(\d+[a-z]*)", canonical_norm)
    if not canonical_match:
        return False

    raw_match = re.search(r"polyprotein(?:orf)?(\d+[a-z]*)$", raw_norm)
    return bool(raw_match and raw_match.group(1) == canonical_match.group(1))


def _combined_orf_matches(raw_value: str, canonical_norm: str) -> bool:
    if canonical_norm != "orf1ab":
        return False
    raw = (raw_value or "").lower()
    raw_compact = normalize_text(raw)
    if raw_compact in {"orf1a/1b", "orf1a/b"}:
        return True
    if re.search(r"orf\s*1a\s*/\s*(?:orf\s*)?1?b", raw):
        return True
    return bool(
        re.search(r"contains\s+orf\s*1a\s+and\s+orf\s*1b", raw)
        or re.search(r"orf\s*1a\s+and\s+orf\s*1b", raw)
    )


# ---------------------------------------------------------------------
# Pipeline
# ---------------------------------------------------------------------

def run_alias_pipeline(
    ref_record: SeqRecord,
    query_records: List[SeqRecord],
    registry_path: Path,
    use_ref_naming: bool = True,
) -> Tuple[List[Dict], Optional[Dict]]:
    """
    Run alias normalization pipeline over all query records.

    Uses the reference record to detect the virus and load the alias config.
    Each query record is then normalized using that config.

    Args:
        ref_record:     Reference GenBank record (used for virus detection)
        query_records:  List of query GenBank records to normalize
        registry_path:  Path to virus_alias_registry.json

    Returns:
        Tuple of:
            - results: list of per-record dicts with keys:
                record_id, strategy, feature_type, status,
                resolved, total, normalized_features
            - payload: aggregated LLM payload dict, or None if all resolved
    """
    from app.src.alias.alias_registry import (
        detect_alias_config_for_record,
        get_detected_virus_name,
    )
    from app.src.features.annotation_strategy import select_feature_type
    from app.src.alias.gene_alias import (
        apply_alias_to_features,
        apply_ref_naming,
        build_canonical_to_ref_map,
        load_alias_lookup,
    )
    from app.src.io.genbank_parser import parse_cds_features, parse_mat_peptides

    # Detect virus and load alias config from reference
    alias_config_path = detect_alias_config_for_record(ref_record, registry_path)
    virus_name = get_detected_virus_name(ref_record, registry_path)

    alias_lookup = load_alias_lookup(alias_config_path) if alias_config_path else {}
    canonical_names = _unique_ordered(list(alias_lookup.values())) if alias_lookup else []

    # Determine expected feature type from reference
    ref_feature_type = select_feature_type(ref_record, alias_lookup or None)

    # Build canonical -> ref name map for output naming
    if use_ref_naming and alias_lookup:
        if ref_feature_type == "mat_peptide":
            ref_features = parse_mat_peptides(ref_record)
        elif ref_feature_type == "CDS":
            ref_features = parse_cds_features(ref_record)
        else:
            ref_features = []
        canonical_to_ref = build_canonical_to_ref_map(ref_features, alias_lookup)
    else:
        canonical_to_ref = {}

    results = []
    payload_records = []

    for record in query_records:
        feature_type = select_feature_type(record, alias_lookup or None)

        needs_lifting = feature_type is None

        if feature_type == "mat_peptide":
            raw_features = parse_mat_peptides(record)
        elif feature_type == "CDS":
            raw_features = parse_cds_features(record)
        else:
            results.append({
                "record_id": record.id,
                "strategy": "tblastn",
                "feature_type": None,
                "status": "no_annotation",
                "resolved": 0,
                "total": 0,
                "normalized_features": [],
            })
            continue

        if alias_config_path is None:
            # Virus not in registry — build alias map from scratch using ref
            results.append({
                "record_id": record.id,
                "strategy": "direct" if not needs_lifting else "tblastn",
                "feature_type": feature_type,
                "status": "new_virus",
                "resolved": 0,
                "total": len(raw_features),
                "normalized_features": raw_features,
            })
            payload_records.append(build_new_virus_payload(record, raw_features, feature_type))
            continue

        normalized = apply_alias_to_features(raw_features, alias_lookup)
        if canonical_to_ref:
            normalized = apply_ref_naming(normalized, canonical_to_ref)
        unresolved = get_unresolved_features(normalized)
        resolved_count = len(normalized) - len(unresolved)

        if not unresolved:
            results.append({
                "record_id": record.id,
                "strategy": "direct" if not needs_lifting else "tblastn",
                "feature_type": feature_type,
                "status": "all_resolved",
                "resolved": resolved_count,
                "total": len(normalized),
                "normalized_features": normalized,
            })
            continue

        results.append({
            "record_id": record.id,
            "strategy": "direct" if not needs_lifting else "tblastn",
            "feature_type": feature_type,
            "status": "has_unresolved",
            "resolved": resolved_count,
            "total": len(normalized),
            "normalized_features": normalized,
        })
        payload_records.append(build_map_aliases_payload(
            record=record,
            unresolved_features=unresolved,
            existing_canonical_names=canonical_names,
            feature_type=feature_type,
        ))

    if not payload_records:
        return results, None

    aggregated_payload = {
        "virus": virus_name,
        "alias_config": str(alias_config_path) if alias_config_path else None,
        "records_with_unresolved": len(payload_records),
        "records": payload_records,
    }

    return results, aggregated_payload


# ---------------------------------------------------------------------
# IO
# ---------------------------------------------------------------------

def save_payload(payload: Dict, output_path: Path) -> None:
    """Save a payload dictionary to a JSON file."""
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, ensure_ascii=False)
