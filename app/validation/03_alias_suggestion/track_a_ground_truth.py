"""
Track A ground truth: the 10 known ambiguous/ignored/unresolved raw names
from the existing PED alias validation
(06_ped_validation/outputs/alias/ped_alias_names_to_review.tsv), with the
expected correct LLM recommendation for each, reasoned from
PED_VALIDATION_REPORT.md and the alias config's own combined-ORF matching
rule (ORF1a + ORF1b wording -> ORF1ab, since ORF1ab is a real PED canonical).

Only 10 items -> too small alone for statistically meaningful precision/
recall, but each item's ground truth is already independently justified in
writing, so this is a clean, low-ambiguity sanity check track.
"""

TRACK_A_UNKNOWN_ITEMS = {
    "101-bp deletion results in frameshift and premature termination": {
        "records": ["PV533621.1"],
        "candidates": ["101-bp deletion results in frameshift and premature termination"],
        "ambiguous": False,
    },
    "1a polyprotein and 1b polyprotein": {
        "records": ["JX112709.1"],
        "candidates": ["1a polyprotein and 1b polyprotein"],
        "ambiguous": False,
    },
    "contains ORF1a and ORF1b": {
        "records": ["KY928065.1"],
        "candidates": ["contains ORF1a and ORF1b"],
        "ambiguous": False,
    },
    "truncated hypothetical protein": {
        "records": ["PV533621.1"],
        "candidates": ["truncated hypothetical protein"],
        "ambiguous": False,
    },
}

TRACK_A_AMBIGUOUS_ITEMS = {
    "mp": {
        "records": ["LT898413.1", "MF577027.1"],
        "candidates": ["mp"],
        "ambiguous": True,
    },
}

# Cases added after the coordinate-evidence audit of the PED alias config.
# Each was previously missed or mishandled by name-only review; each is
# decidable once positional context is supplied.
TRACK_A_POSITIONAL_ITEMS = {
    "sM": {"records": ["EF185992.1"], "candidates": ["sM"], "ambiguous": False},
    "POL1": {"records": ["EF185992.1"], "candidates": ["POL1"], "ambiguous": False},
    "nucleoprotein": {"records": ["KM189367.2"], "candidates": ["nucleoprotein"], "ambiguous": False},
    "RNA-dependent RNA polymerase": {
        "records": ["LR812926.1", "LR812930.1", "LR812932.1"],
        "candidates": ["RNA-dependent RNA polymerase"],
        "ambiguous": False,
    },
    "putative coronavirus nsp12": {
        "records": ["LR812926.1"], "candidates": ["putative coronavirus nsp12"], "ambiguous": False,
    },
}

# ground_truth_action: one of save_alias / ignore / skip / move_to_ambiguous
# (matches LLM_REVIEW_ACTIONS in app/src/llm/alias_review.py)
TRACK_A_GROUND_TRUTH = {
    "101-bp deletion results in frameshift and premature termination": {
        "action": "ignore", "canonical": None,
        "note": "Descriptive mutation note, not a gene name.",
    },
    "1a polyprotein and 1b polyprotein": {
        "action": "save_alias", "canonical": "ORF1ab",
        "note": "Combined ORF1a+ORF1b wording; ORF1ab is a real PED canonical (granularity mismatch case).",
    },
    "contains ORF1a and ORF1b": {
        "action": "save_alias", "canonical": "ORF1ab",
        "note": "Same combined-ORF pattern as above.",
    },
    "truncated hypothetical protein": {
        "action": "ignore", "canonical": None,
        "note": "Generic + truncation caveat, not a specific gene.",
    },
    "mp": {
        "action": "save_alias", "canonical": "ORF3",
        "note": (
            "REVISED after coordinate audit. Previously scored move_to_ambiguous on the "
            "assumption that 'mp' clashes between ORF3 accessory and membrane protein. "
            "Position settles it: in MF577027.1 'mp' is 24751..25425 (675bp), flanked by S "
            "and E, and the PEDV gene order is 5'UTR-ORF1a-ORF1b-S-ORF3-E-M-N-3'UTR, so the "
            "only slot between S and E is ORF3. M sits later at 25644..26324. The old label "
            "was a limitation of name-only reasoning, not a genuine ambiguity."
        ),
    },
    "sM": {
        "action": "save_alias", "canonical": "E",
        "note": (
            "sM (small membrane) is a documented synonym of the coronavirus envelope "
            "protein E. Position confirms: 25444..25674 (231bp) in the slot between ORF3 "
            "and M. Name-only review reads 'M' in the string and drifts toward membrane "
            "protein; the config had it in excluded_names, which silently dropped gene E."
        ),
    },
    "POL1": {
        "action": "save_alias", "canonical": "ORF1ab",
        "note": (
            "Spans 297..20641 = 20345bp, the full replicase, matching records that split "
            "it into ORF1a + ORF1b. A single feature covering both is ORF1ab (pp1ab)."
        ),
    },
    "nucleoprotein": {
        "action": "save_alias", "canonical": "N",
        "note": "Nucleoprotein = nucleocapsid = N. Position confirms: last CDS, directly after M.",
    },
    "RNA-dependent RNA polymerase": {
        "action": "ignore", "canonical": None,
        "note": (
            "Real protein name but NOT a gene-level feature: it is a mat_peptide at "
            "12595..15374 lying wholly inside ORF1ab, i.e. nsp12 cleaved from pp1ab. "
            "Mapping it to ORF1b would be a false positive. This is the case name-only "
            "review is most likely to get actively WRONG rather than merely miss."
        ),
    },
    "putative coronavirus nsp12": {
        "action": "ignore", "canonical": None,
        "note": "mat_peptide inside ORF1ab; cleavage product, not a gene.",
    },
    # also include the 5 "ignored" rows already resolved deterministically
    # (not LLM-uncertain, but useful as an easy-case sanity check if the
    # unresolved-name LLM assist is ever run on them)
    "replicase polyprotein": {"action": "ignore", "canonical": None, "note": "Generic replicase description."},
    "hypothetical protein": {"action": "ignore", "canonical": None, "note": "Generic, no gene-specific info."},
    "polyprotein": {"action": "ignore", "canonical": None, "note": "Generic, no gene-specific info."},
    "HNZK1": {"action": "ignore", "canonical": None, "note": "Strain/isolate prefix, not a gene."},
    "unknown": {"action": "ignore", "canonical": None, "note": "Placeholder text, no gene-specific info."},
}
