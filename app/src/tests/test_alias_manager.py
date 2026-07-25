import json

import pytest

from app.src.alias.alias_manager import (
    remove_registry_entry,
    save_validated_alias_config,
)
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from app.src.alias.alias_bootstrap import (
    apply_approved_alias_suggestions,
    build_seed_alias_config_from_ref,
)


def test_remove_registry_entry_can_archive_alias_config(tmp_path):
    registry_path = tmp_path / "registry.json"
    config_path = tmp_path / "virus_alias.json"
    config_path.write_text('{"virus": "Test"}\n', encoding="utf-8")
    registry_path.write_text(
        json.dumps(
            {
                "viruses": [
                    {
                        "virus_name": "Test virus",
                        "keywords": ["test virus"],
                        "alias_config": "virus_alias.json",
                    },
                    {
                        "virus_name": "Other virus",
                        "keywords": ["other virus"],
                        "alias_config": "other_alias.json",
                    },
                ]
            }
        ),
        encoding="utf-8",
    )

    registry_backup, config_backup = remove_registry_entry(
        registry_path,
        alias_config="virus_alias.json",
        root=tmp_path,
        archive_alias_config=True,
    )

    registry = json.loads(registry_path.read_text(encoding="utf-8"))
    assert [entry["virus_name"] for entry in registry["viruses"]] == ["Other virus"]
    assert registry_backup.exists()
    assert config_backup.exists()
    assert not config_path.exists()


def test_remove_registry_entry_can_keep_alias_config(tmp_path):
    registry_path = tmp_path / "registry.json"
    config_path = tmp_path / "virus_alias.json"
    config_path.write_text('{"virus": "Test"}\n', encoding="utf-8")
    registry_path.write_text(
        json.dumps(
            {
                "viruses": [
                    {
                        "virus_name": "Test virus",
                        "keywords": ["test virus"],
                        "alias_config": "virus_alias.json",
                    },
                ]
            }
        ),
        encoding="utf-8",
    )

    _, config_backup = remove_registry_entry(
        registry_path,
        alias_config="virus_alias.json",
        root=tmp_path,
        archive_alias_config=False,
    )

    registry = json.loads(registry_path.read_text(encoding="utf-8"))
    assert registry["viruses"] == []
    assert config_backup is None
    assert config_path.exists()


def test_save_validated_alias_config_creates_backup(tmp_path):
    config_path = tmp_path / "virus_alias.json"
    config_path.write_text(
        json.dumps({"virus": "Test", "canonical_names": {"A": []}}),
        encoding="utf-8",
    )

    backup = save_validated_alias_config(
        config_path,
        {"virus": "Test", "canonical_names": {"A": ["alpha"]}},
    )

    assert backup.exists()
    assert backup.parent.name == "backups"
    saved = json.loads(config_path.read_text(encoding="utf-8"))
    assert saved["canonical_names"]["A"] == ["alpha"]


def test_save_validated_alias_config_rejects_conflicts_without_overwrite(tmp_path):
    config_path = tmp_path / "virus_alias.json"
    original = {"virus": "Test", "canonical_names": {"A": []}}
    config_path.write_text(json.dumps(original), encoding="utf-8")

    with pytest.raises(ValueError, match="validation failed"):
        save_validated_alias_config(
            config_path,
            {"virus": "Test", "canonical_names": {"A": ["x"], "B": ["x"]}},
        )

    assert json.loads(config_path.read_text(encoding="utf-8")) == original
    assert not (tmp_path / "backups").exists()


def test_apply_approved_alias_suggestions_creates_backup(tmp_path):
    config_path = tmp_path / "virus_alias.json"
    config_path.write_text(
        json.dumps({
            "virus": "Test",
            "canonical_names": {"A": []},
            "ignored_names": ["alpha"],
            "ambiguous_names": ["alpha"],
        }),
        encoding="utf-8",
    )

    apply_approved_alias_suggestions(
        config_path,
        approved_rows=[{"raw_value": "alpha", "canonical_name": "A"}],
        ignored_rows=[{"raw_value": "ignored label"}],
        ambiguous_rows=[{"raw_value": "shared label"}],
    )

    saved = json.loads(config_path.read_text(encoding="utf-8"))
    assert saved["canonical_names"]["A"] == ["alpha"]
    assert saved["excluded_names"] == ["ignored label", "shared label"]
    assert "ignored_names" not in saved
    assert "ambiguous_names" not in saved
    backups = list((tmp_path / "backups").glob("virus_alias.*.json"))
    assert len(backups) == 1


def test_seed_alias_config_does_not_ignore_contextual_descriptions():
    ref_record = SeqRecord(Seq("ATG"), id="ref")
    config = build_seed_alias_config_from_ref(
        ref_record,
        ref_features=[{"name": "S"}, {"name": "M"}, {"name": "E"}],
        virus_name="Test virus",
    )

    excluded = set(config["excluded_names"])
    assert "unknown" in excluded
    assert "hypothetical protein" in excluded
    assert "envelope protein" not in excluded
    assert "membrane protein" not in excluded
    assert "glycoprotein" not in excluded
    assert "polyprotein" not in excluded
