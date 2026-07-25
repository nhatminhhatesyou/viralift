from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from app.src.lifting.base import LiftedFeature
from app.src.pipeline import PipelineConfig, run_pipeline


def _record(record_id: str) -> SeqRecord:
    return SeqRecord(Seq("ATGAAATAA"), id=record_id)


def _feature(status: str, method: str) -> LiftedFeature:
    return LiftedFeature(
        name="GENE",
        source_name=None,
        ref_start=1,
        ref_end=9,
        strand="+",
        query_start=1,
        query_end=9,
        sequence="ATGAAATAA",
        coverage=1.0,
        status=status,
        method=method,
    )


def test_run_pipeline_routes_direct_and_tblastn(monkeypatch):
    calls = []

    def fake_strategy(record, alias_lookup, allowed_types=None):
        return ("direct", "CDS") if record.id == "direct" else ("tblastn", None)

    def fake_direct_extract(**kwargs):
        calls.append(("direct", kwargs["query_record"].id))
        return [_feature("ok", "direct")]

    def fake_tblastn(**kwargs):
        calls.append(("tblastn", kwargs["query_record"].id))
        return [_feature("ok_rescued", "tblastn")]

    monkeypatch.setattr("app.src.pipeline.get_strategy", fake_strategy)
    monkeypatch.setattr("app.src.pipeline.direct_extract_with_alias", fake_direct_extract)
    monkeypatch.setattr("app.src.pipeline.process_one_query_record", fake_tblastn)

    result = run_pipeline(
        ref_record=_record("ref"),
        query_records=[_record("direct"), _record("lift")],
        ref_features=[{"name": "GENE", "start": 1, "end": 9}],
        ref_feature_type="CDS",
        alias_lookup={},
        config=PipelineConfig(),
    )

    assert calls == [("direct", "direct"), ("tblastn", "lift")]
    assert result.direct_count == 1
    assert result.lifted_count == 1
    assert result.summary["ok"] == 1
    assert result.summary["ok_rescued"] == 1


def test_run_pipeline_can_capture_record_errors(monkeypatch):
    monkeypatch.setattr("app.src.pipeline.get_strategy", lambda *_, **__: ("tblastn", None))

    def broken_tblastn(**_kwargs):
        raise RuntimeError("boom")

    monkeypatch.setattr("app.src.pipeline.process_one_query_record", broken_tblastn)

    result = run_pipeline(
        ref_record=_record("ref"),
        query_records=[_record("q1")],
        ref_features=[{"name": "GENE", "start": 1, "end": 9}],
        ref_feature_type="CDS",
        alias_lookup={},
        config=PipelineConfig(catch_record_errors=True),
    )

    assert result.all_results == [("q1", [])]
    assert result.errors == [{"record_id": "q1", "error": "boom"}]
