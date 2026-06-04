import csv
from pathlib import Path
from typing import Dict, List, Tuple


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

    Args:
        all_results: List of (query_id, lifted_features).

    Returns:
        Dict with counts per status key.
    """
    summary = {
        "ok": 0,
        "ok_rescued": 0,
        "invalid_boundaries": 0,
        "low_coverage": 0,
        "no_hit": 0,
        "translation_fail": 0,
        "unresolved_name": 0,
        "ambiguous_name": 0,
        "not_in_reference": 0,
    }
    for _, results in all_results:
        for lifted in results:
            status = lifted.status
            if status in summary:
                summary[status] += 1
    return summary


def write_results_tsv(all_results: List[Tuple[str, List]], out_path: Path) -> None:
    """
    Write extracted feature results from all query records to a TSV file.

    Args:
        all_results: List of (query_id, lifted_features).
        out_path:    Output TSV path.
    """
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

    if not rows:
        return

    with open(out_path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

