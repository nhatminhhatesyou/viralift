from dataclasses import dataclass
from typing import Optional


# Single source of truth for every status a LiftedFeature can carry.
# summarize_counts(), main.py and the UI all iterate over this list so a new
# status can never be silently dropped from the run summary again.
TBLASTN_STATUSES = (
    "ok",
    "ok_rescued",
    "ok_extrapolated",
    "ok_no_start_codon",
    "ok_gap_filled",
    "invalid_boundaries",
    "low_coverage",
    "low_identity",
    "low_coverage_and_identity",
    "no_hit",
    "translation_fail",
)

DIRECT_STATUSES = (
    "unresolved_name",
    "not_in_reference",
)

# Ordered list of all statuses, deduplicated while preserving order.
ALL_STATUSES = tuple(
    dict.fromkeys(TBLASTN_STATUSES + DIRECT_STATUSES)
)

STATUS_META = {
    "ok": {"label": "OK", "category": "pass", "is_pass": True},
    "ok_rescued": {"label": "OK (rescued)", "category": "pass", "is_pass": True},
    "ok_extrapolated": {"label": "OK (extrapolated)", "category": "pass", "is_pass": True},
    "ok_no_start_codon": {"label": "OK (frameshift / no start codon)", "category": "pass", "is_pass": True},
    "ok_gap_filled": {"label": "OK (gap-filled from neighbors)", "category": "pass", "is_pass": True},
    "invalid_boundaries": {"label": "Invalid boundaries", "category": "review", "is_pass": False},
    "low_coverage": {"label": "Low coverage", "category": "review", "is_pass": False},
    "low_identity": {"label": "Low identity", "category": "review", "is_pass": False},
    "low_coverage_and_identity": {"label": "Low coverage + identity", "category": "review", "is_pass": False},
    "no_hit": {"label": "No hit", "category": "review", "is_pass": False},
    "translation_fail": {"label": "Translation fail", "category": "review", "is_pass": False},
    "unresolved_name": {"label": "Unresolved names", "category": "review", "is_pass": False},
    "not_in_reference": {"label": "Not in reference", "category": "pass", "is_pass": True},
}

STATUS_LABELS = {
    status: meta["label"]
    for status, meta in STATUS_META.items()
}

PASS_STATUSES = {
    status
    for status, meta in STATUS_META.items()
    if meta.get("is_pass")
}


@dataclass
class LiftedFeature:
    """
    Unified output of the lifting/extraction pipeline (tblastn or direct).

    Coordinates are 1-based inclusive, same convention as GenBank.
    """
    name: str                        # canonical / standardized gene name
    source_name: Optional[str]       # raw original name from the annotation (before alias)

    ref_start: Optional[int]
    ref_end: Optional[int]
    strand: str                      # "+" or "-"

    query_start: Optional[int]       # None if unmapped
    query_end: Optional[int]
    sequence: Optional[str]          # extracted nucleotide sequence

    coverage: float                  # fraction of ref feature mapped
    status: str                      # ok | ok_rescued | invalid_boundaries | low_coverage | no_hit | ...
    method: str                      # "tblastn" | "direct"

    has_start_codon: Optional[bool] = None
    has_stop_codon: Optional[bool] = None
    in_frame: Optional[bool] = None
    rescue_offset: Optional[int] = None
    rescue_target: Optional[str] = None
    rescue_action: Optional[str] = None

    # Diagnostics: raw merged-HSP coordinates (before any rescue/extrapolation) and
    # how many N-terminal reference residues the HSP failed to align.
    raw_start: Optional[int] = None
    raw_end: Optional[int] = None
    n_term_missing_aa: Optional[int] = None

    # Extra engine-specific info
    identity: Optional[float] = None     # tblastn: % identity of best HSP
    score: Optional[float] = None        # tblastn: bit score

    def to_dict(self) -> dict:
        return {
            "name": self.name,
            "source_name": self.source_name,
            "ref_start": self.ref_start,
            "ref_end": self.ref_end,
            "strand": self.strand,
            "start": self.query_start,
            "end": self.query_end,
            "sequence": self.sequence,
            "coverage": self.coverage,
            "status": self.status,
            "method": self.method,
            "has_start_codon": self.has_start_codon,
            "has_stop_codon": self.has_stop_codon,
            "in_frame": self.in_frame,
            "rescue_offset": self.rescue_offset,
            "rescue_target": self.rescue_target,
            "rescue_action": self.rescue_action,
            "raw_start": self.raw_start,
            "raw_end": self.raw_end,
            "n_term_missing_aa": self.n_term_missing_aa,
            "identity": self.identity,
            "score": self.score,
        }
