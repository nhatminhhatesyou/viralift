from dataclasses import dataclass, field
from typing import Optional


@dataclass
class LiftedFeature:
    """
    Unified output from any lifting engine (tblastn, minimap, etc.).

    Coordinates are 1-based inclusive, same convention as GenBank.
    """
    name: str                        # canonical / standardized gene name
    source_name: Optional[str]       # raw original name from the annotation (before alias)

    ref_start: int
    ref_end: int
    strand: str                      # "+" or "-"

    query_start: Optional[int]       # None if unmapped
    query_end: Optional[int]
    sequence: Optional[str]          # extracted nucleotide sequence

    coverage: float                  # fraction of ref feature mapped
    status: str                      # ok | ok_rescued | invalid_boundaries | low_coverage | no_hit | ...
    method: str                      # "tblastn" | "minimap_transfer"

    has_start_codon: Optional[bool] = None
    has_stop_codon: Optional[bool] = None
    in_frame: Optional[bool] = None
    rescue_offset: Optional[int] = None

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
            "identity": self.identity,
            "score": self.score,
        }
