import json

from app.src.alias.alias_manager import remove_registry_entry


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
