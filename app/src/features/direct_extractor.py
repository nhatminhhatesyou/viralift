from typing import Dict, List

from Bio.SeqRecord import SeqRecord

from app.src.alias.gene_alias import apply_alias_to_features
from app.src.io.genbank_parser import parse_cds_features, parse_mat_peptides
from app.src.lifting.base import LiftedFeature


"""
Module: direct_extractor.py

Purpose:
    Extract feature sequences directly from query annotation coordinates
    and normalize names using the alias lookup — no alignment needed.

Use case:
    Query records that already carry gene-level annotation (CDS or mat_peptide).
    Much faster than tblastn lifting; just parse, rename, and slice sequence.
"""


def direct_extract_with_alias(
    query_record: SeqRecord,
    query_feature_type: str,
    ref_features: List[Dict],
    alias_lookup: Dict[str, str],
) -> List[LiftedFeature]:
    """
    Extract annotated features from a query record and normalize their names.

    Args:
        query_record:       Query genome with existing annotation.
        query_feature_type: "CDS" or "mat_peptide".
        ref_features:       Reference features (alias-normalized).
                            Used to populate ref_start / ref_end by name match.
        alias_lookup:       Alias normalization lookup {raw_name: canonical}.

    Returns:
        List of LiftedFeature objects with method="direct". Status values:
        - ok: name resolved and exists in the reference feature set
        - unresolved_name: name was not resolved by the alias lookup
        - ambiguous_name: name is known ambiguous and needs user resolution
        - not_in_reference: name resolved, but the selected reference lacks it
    """
    if query_feature_type == "mat_peptide":
        query_features = parse_mat_peptides(query_record)
    else:
        query_features = parse_cds_features(query_record)

    if alias_lookup:
        query_features = apply_alias_to_features(query_features, alias_lookup)

    # canonical name → ref feature (for ref_start / ref_end)
    ref_by_name: Dict[str, Dict] = {f["name"]: f for f in ref_features}

    ignored: List[str] = []
    results: List[LiftedFeature] = []
    for qf in query_features:
        if qf.get("name_source") == "ignored":
            ignored.append(qf.get("raw_name") or qf["name"])
            continue

        name = qf["name"]
        name_source = qf.get("name_source")
        ref_match = ref_by_name.get(name)

        if name_source == "ambiguous":
            status = "ambiguous_name"
        elif alias_lookup and name_source == "raw":
            status = "unresolved_name"
        elif ref_match is None:
            status = "not_in_reference"
        else:
            status = "ok"

        start  = qf["start"]
        end    = qf["end"]
        strand = qf.get("strand", "+")

        seq = query_record.seq[start - 1: end]
        if strand == "-":
            seq = seq.reverse_complement()

        results.append(LiftedFeature(
            name=name,
            source_name=(
                qf.get("raw_name")
                if name_source in ("alias", "product_alias", "alias_conflict_resolved")
                else None
            ),
            ref_start=ref_match["start"] if ref_match else None,
            ref_end=ref_match["end"]   if ref_match else None,
            strand=strand,
            method="direct",
            query_start=start,
            query_end=end,
            sequence=str(seq),
            coverage=1.0,
            status=status,
        ))

    if ignored:
        print(
            f"  [WARN] {query_record.id}: {len(ignored)} feature(s) skipped "
            f"(ambiguous name, cannot resolve to canonical): {ignored}"
        )

    return results
