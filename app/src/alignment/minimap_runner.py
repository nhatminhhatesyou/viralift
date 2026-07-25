from pathlib import Path
import subprocess


PRESETS_FALLBACK = ["asm5", "asm10", "asm20"]


def run_minimap2(
    ref_fasta: Path,
    query_fasta: Path,
    out_sam: Path,
    preset: str = "asm5",
    quiet: bool = False,
) -> None:
    """Run minimap2 with given preset and write SAM output."""
    cmd = [
        "minimap2",
        "-a",
        "-x", preset,
        str(ref_fasta),
        str(query_fasta),
    ]

    with open(out_sam, "w", encoding="utf-8") as sam_file:
        if quiet:
            subprocess.run(cmd, stdout=sam_file, stderr=subprocess.DEVNULL, check=True)
        else:
            subprocess.run(cmd, stdout=sam_file, check=True)


def run_minimap2_with_fallback(
    ref_fasta: Path,
    query_fasta: Path,
    out_sam: Path,
    quiet: bool = False,
) -> str:
    """
    Try minimap2 presets in order: asm5 -> asm10 -> asm20.
    Returns the preset that produced a primary alignment.
    Raises ValueError if all presets fail.
    """
    from app.src.alignment.sam_lifter import get_primary_alignment

    for preset in PRESETS_FALLBACK:
        run_minimap2(ref_fasta, query_fasta, out_sam, preset=preset, quiet=quiet)
        try:
            get_primary_alignment(out_sam)
            return preset
        except ValueError:
            continue

    raise ValueError(f"No primary alignment found with any preset: {PRESETS_FALLBACK}")
