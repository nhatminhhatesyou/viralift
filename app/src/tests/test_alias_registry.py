"""
Regression tests for packaged-config path resolution.

Guards against the bug where a pip-installed CLI run from outside the repo
silently fell back to raw names because the registry / alias-config paths were
resolved relative to the current working directory.
"""
from pathlib import Path

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from app.src.alias import alias_registry as ar
from app.src.alias.alias_registry import (
    DEFAULT_REGISTRY_PATH,
    PACKAGE_CONFIG_DIR,
    detect_alias_config_for_record,
    get_detected_virus_name,
    _resolve_alias_config,
)


def test_default_registry_path_is_absolute_and_exists():
    assert DEFAULT_REGISTRY_PATH.is_absolute()
    assert DEFAULT_REGISTRY_PATH.exists()
    assert PACKAGE_CONFIG_DIR.is_dir()


def _ped_record() -> SeqRecord:
    rec = SeqRecord(Seq("ATG"), id="X", name="X",
                    description="Porcine epidemic diarrhea virus strain test")
    rec.annotations["organism"] = "Porcine epidemic diarrhea virus"
    return rec


def test_detect_resolves_config_regardless_of_cwd(monkeypatch, tmp_path):
    # Simulate an installed CLI invoked from an unrelated directory.
    monkeypatch.chdir(tmp_path)

    rec = _ped_record()
    assert get_detected_virus_name(rec, DEFAULT_REGISTRY_PATH) == \
        "Porcine epidemic diarrhea virus"

    cfg = detect_alias_config_for_record(rec, DEFAULT_REGISTRY_PATH)
    assert cfg is not None
    assert cfg.is_absolute()
    assert cfg.exists()


def test_resolve_alias_config_falls_back_to_registry_dir(monkeypatch, tmp_path):
    # Entry written as a cwd-relative path that does not exist in cwd, but the
    # basename sits next to the registry file → must resolve there.
    # chdir into a scratch dir so the cwd-relative candidate genuinely misses.
    work = tmp_path / "elsewhere"
    work.mkdir()
    monkeypatch.chdir(work)

    registry = tmp_path / "registry.json"
    registry.write_text("{}", encoding="utf-8")
    cfg_file = tmp_path / "prrsv_alias.json"
    cfg_file.write_text("{}", encoding="utf-8")

    resolved = _resolve_alias_config("app/config/prrsv_alias.json", registry)
    assert resolved == cfg_file
    assert resolved.exists()


def test_resolve_alias_config_keeps_absolute(tmp_path):
    abs_path = tmp_path / "x.json"
    assert _resolve_alias_config(str(abs_path), DEFAULT_REGISTRY_PATH) == abs_path
