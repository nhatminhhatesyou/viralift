from pathlib import Path
import subprocess


def run_minimap2(ref_fasta: Path, query_fasta: Path, out_sam: Path, quiet: bool = False) -> None:
    """Run minimap2 and write SAM output."""
    cmd = [
        "minimap2",
        "-a",
        "-x", "asm5",
        str(ref_fasta),
        str(query_fasta),
    ]

    with open(out_sam, "w", encoding="utf-8") as sam_file:
        if quiet:
            subprocess.run(cmd, stdout=sam_file, stderr=subprocess.DEVNULL, check=True)
        else:
            subprocess.run(cmd, stdout=sam_file, check=True)