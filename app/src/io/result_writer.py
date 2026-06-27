import csv
from pathlib import Path
from typing import Dict, List, Tuple

from app.src.lifting.base import ALL_STATUSES


"""
Module: result_writer.py

Purpose:
    Write pipeline results to disk and compute run-level summaries.

Functions:
    summarize_counts()  — aggregate status counts across all records.
    write_results_tsv() — full annotation table as TSV.
"""


def summarize_counts(all_results: List[Tuple[str, List]]) -> Dict[str, int]:
    """
    Aggregate feature status counts across all processed query records.

    The summary is seeded from ALL_STATUSES (the single source of truth in
    base.py), so any status a LiftedFeature can carry is always represented —
    no status can be silently dropped. Any unexpected status still gets counted
    under its own key rather than being discarded.

    Args:
        all_results: List of (query_id, lifted_features).

    Returns:
        Dict with counts per status key.
    """
    summary = {status: 0 for status in ALL_STATUSES}
    for _, results in all_results:
        for lifted in results:
            status = lifted.status
            summary[status] = summary.get(status, 0) + 1
    return summary


def write_results_tsv(all_results: List[Tuple[str, List]], out_path: Path) -> None:
    """
    Write extracted feature results from all query records to a TSV file.

    Args:
        all_results: List of (query_id, lifted_features).
        out_path:    Output TSV path.
    """
    fieldnames = [
        "query_id", "name", "source_name", "ref_start", "ref_end",
        "start", "end", "strand", "method", "status", "coverage",
        "identity", "score", "has_start_codon", "has_stop_codon",
        "in_frame", "rescue_offset", "length",
    ]

    rows: List[Dict] = []
    for query_id, results in all_results:
        for lifted in results:
            rows.append({
                "query_id":       query_id,
                "name":           lifted.name,
                "source_name":    lifted.source_name or "",
                "ref_start":      lifted.ref_start,
                "ref_end":        lifted.ref_end,
                "start":          lifted.query_start,
                "end":            lifted.query_end,
                "strand":         lifted.strand,
                "method":         lifted.method,
                "status":         lifted.status,
                "coverage":       lifted.coverage,
                "identity":       lifted.identity,
                "score":          lifted.score,
                "has_start_codon": lifted.has_start_codon,
                "has_stop_codon":  lifted.has_stop_codon,
                "in_frame":        lifted.in_frame,
                "rescue_offset":   lifted.rescue_offset,
                "length":         len(lifted.sequence) if lifted.sequence else None,
            })

    # Always write the file — even with no rows, emit a header-only TSV so
    # downstream consumers can rely on the output path existing.
    with open(out_path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

