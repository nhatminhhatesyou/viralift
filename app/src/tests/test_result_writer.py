"""
Unit tests for result_writer — regression coverage for the status-summary bug
where ok_extrapolated features were silently dropped, and for the empty-output
behavior.
"""
from pathlib import Path

from app.src.io.result_writer import summarize_counts, write_results_tsv
from app.src.lifting.base import LiftedFeature, ALL_STATUSES


def _lf(status: str, sequence: str = "ATGAAATAA") -> LiftedFeature:
    return LiftedFeature(
        name="GENE",
        source_name=None,
        ref_start=1,
        ref_end=len(sequence),
        strand="+",
        query_start=1,
        query_end=len(sequence),
        sequence=sequence,
        coverage=1.0,
        status=status,
        method="tblastn",
    )


def test_summarize_counts_includes_ok_extrapolated():
    # Regression: ok_extrapolated was produced by the mat_peptide path but
    # missing from the summary dict, so it never showed up in the run report.
    summary = summarize_counts([("q1", [_lf("ok_extrapolated")])])
    assert summary["ok_extrapolated"] == 1


def test_summarize_counts_covers_every_known_status():
    # Every status the pipeline can emit must be representable in the summary.
    results = [(f"q{i}", [_lf(s)]) for i, s in enumerate(ALL_STATUSES)]
    summary = summarize_counts(results)
    for status in ALL_STATUSES:
        assert summary[status] == 1


def test_summarize_counts_handles_unexpected_status_gracefully():
    summary = summarize_counts([("q1", [_lf("some_future_status")])])
    assert summary["some_future_status"] == 1


def test_write_results_tsv_writes_header_when_empty(tmp_path: Path):
    # Regression: previously returned early and wrote no file at all.
    out = tmp_path / "out.tsv"
    write_results_tsv([], out)
    assert out.exists()
    header = out.read_text(encoding="utf-8").splitlines()[0]
    assert header.startswith("query_id\t")
    assert "status" in header


def test_write_results_tsv_writes_rows(tmp_path: Path):
    out = tmp_path / "out.tsv"
    write_results_tsv([("q1", [_lf("ok")])], out)
    lines = out.read_text(encoding="utf-8").splitlines()
    assert len(lines) == 2  # header + one row
    assert lines[1].startswith("q1\t")
