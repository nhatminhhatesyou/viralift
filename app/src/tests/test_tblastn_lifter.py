"""
Tests for the tblastn lifter.

Split into two groups:
  - Unit tests that need no BLAST binary (environment guard).
  - A small golden integration test, skipped automatically when tblastn is not
    on PATH, that exercises the real lift path end-to-end on bundled data.
"""
import shutil
from pathlib import Path

import pytest

from app.src.lifting import tblastn_lifter as tl
from app.src.lifting.tblastn_lifter import (
    BlastNotInstalledError,
    ensure_tblastn_available,
)
from app.src.alias.alias_registry import DEFAULT_REGISTRY_PATH

HAS_TBLASTN = shutil.which("tblastn") is not None
DATA = Path(__file__).resolve().parents[2] / "data"

requires_blast = pytest.mark.skipif(
    not HAS_TBLASTN, reason="tblastn (BLAST+) not installed on PATH"
)


# --------------------------------------------------------------------------
# Environment guard — no BLAST binary needed
# --------------------------------------------------------------------------

def test_ensure_tblastn_available_raises_when_missing(monkeypatch):
    monkeypatch.setattr(tl.shutil, "which", lambda _name: None)
    with pytest.raises(BlastNotInstalledError):
        ensure_tblastn_available()


def test_ensure_tblastn_available_returns_path_when_present(monkeypatch):
    monkeypatch.setattr(tl.shutil, "which", lambda _name: "/usr/bin/tblastn")
    assert ensure_tblastn_available() == "/usr/bin/tblastn"


def test_lift_all_tblastn_fails_fast_without_blast(monkeypatch):
    # When BLAST is absent, the pipeline must raise a clear error instead of
    # silently turning every gene into no_hit.
    monkeypatch.setattr(tl.shutil, "which", lambda _name: None)
    with pytest.raises(BlastNotInstalledError):
        tl.lift_all_tblastn(
            ref_features=[{"name": "X", "start": 1, "end": 30, "strand": "+"}],
            ref_record=None,
            query_record=None,
        )


# --------------------------------------------------------------------------
# Golden integration test — requires BLAST+
# --------------------------------------------------------------------------

@requires_blast
@pytest.mark.parametrize(
    "ref_file,query_file",
    [
        ("PRRS/PRRS_ref_test.gb", "PRRS/PRRSV_test.gb"),
        ("FMD/FMD_ref_test.gb", "FMD/FMD_test.gb"),
    ],
)
def test_lift_pipeline_golden(ref_file, query_file):
    from app.src.io.genbank_parser import load_single_genbank, load_genbank_records
    from app.src.features.ref_loader import prepare_reference_features
    from app.src.features.annotation_strategy import get_strategy
    from app.src.lifting.tblastn_lifter import process_one_query_record

    ref_path = DATA / ref_file
    query_path = DATA / query_file
    if not ref_path.exists() or not query_path.exists():
        pytest.skip(f"bundled data missing: {ref_file} / {query_file}")

    ref_record = load_single_genbank(ref_path)
    query_records = load_genbank_records(query_path)
    assert query_records

    ref_features, ref_feature_type, _cfg, _virus, alias_lookup = (
        prepare_reference_features(
            ref_record=ref_record,
            alias_config_arg=None,
            alias_registry_arg=str(DEFAULT_REGISTRY_PATH),
        )
    )
    assert ref_features

    # Force the tblastn path on at least the first record and assert it produces
    # the expected number of results with valid status codes — a stable
    # regression anchor that any refactor of the post-HSP logic must preserve.
    from app.src.lifting.base import ALL_STATUSES

    results = process_one_query_record(
        ref_record=ref_record,
        query_record=query_records[0],
        ref_cds=ref_features,
        ref_feature_type=ref_feature_type,
        min_coverage=0.5,
    )
    assert len(results) == len(ref_features)
    assert all(r.status in ALL_STATUSES for r in results)
    # At least one gene should lift cleanly on a same-virus reference.
    assert any(r.status in ("ok", "ok_rescued", "ok_extrapolated") for r in results)
