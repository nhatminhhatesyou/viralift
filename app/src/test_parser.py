from pathlib import Path
from src.genbank_parser import load_genbank_records, parse_cds_features, get_record_metadata


def test_file(gb_path: Path) -> None:
    records = load_genbank_records(gb_path)
    print(f"\nTesting file: {gb_path.name}")
    print(f"Total records: {len(records)}")

    for i, record in enumerate(records, start=1):
        meta = get_record_metadata(record)
        cds_list = parse_cds_features(record)

        print(f"\nRecord {i}")
        print(f"ID: {meta['id']}")
        print(f"Organism: {meta['organism']}")
        print(f"Length: {meta['length']}")
        print(f"CDS count: {len(cds_list)}")

        for cds in cds_list[:3]:
            print(f"  - {cds['name']} | {cds['start']}-{cds['end']} | {cds['strand']}")


if __name__ == "__main__":
    data_dir = Path("data")

    test_file(data_dir / "FMD_test.gb")
    test_file(data_dir / "PRRSV_test.gb")