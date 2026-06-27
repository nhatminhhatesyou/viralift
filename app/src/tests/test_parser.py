from pathlib import Path

import pytest

from app.src.io.genbank_parser import (
    get_record_metadata,
    load_genbank_records,
    parse_cds_features,
)

DATA = Path(__file__).resolve().parents[2] / "data"


def _summarize_file(gb_path: Path) -> int:
    """Parse one GenBank file and return the number of records. Helper, not a test."""
    records = load_genbank_records(gb_path)
    for record in records:
        meta = get_record_metadata(record)
        assert meta["id"]
        # parse_cds_features must not raise on real-world records
        parse_cds_features(record)
    return len(records)


@pytest.mark.parametrize("filename", ["FMD/FMD_test.gb", "PRRS/PRRSV_test.gb"])
def test_parse_bundled_genbank(filename: str):
    gb_path = DATA / filename
    if not gb_path.exists():
        pytest.skip(f"bundled data missing: {filename}")
    assert _summarize_file(gb_path) >= 1
