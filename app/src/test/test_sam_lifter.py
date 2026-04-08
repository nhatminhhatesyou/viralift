from pathlib import Path

from app.src.alignment.sam_lifter import (
    get_primary_alignment,
    build_ref_to_query_map,
    get_alignment_summary,
)


if __name__ == "__main__":
    app_dir = Path(__file__).resolve().parents[1]
    data_dir = app_dir / "data"
    sam_path = data_dir / "tmp_alignment.sam"

    aln = get_primary_alignment(sam_path)
    summary = get_alignment_summary(aln)
    ref_to_query = build_ref_to_query_map(aln)

    print("Alignment summary:")
    for key, value in summary.items():
        print(f"{key}: {value}")

    print(f"\nMapped positions: {len(ref_to_query)}")

    print("\nFirst 10 coordinate pairs:")
    for i, (ref_pos, query_pos) in enumerate(ref_to_query.items()):
        if i >= 10:
            break
        print(f"ref {ref_pos} -> query {query_pos}")