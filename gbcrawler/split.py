"""Split a combined GenBank file into one file per virus species.

Species are matched against ViraLift's ``virus_alias_registry.json`` so each
output file lines up with a known reference + alias config. Records whose
organism doesn't match any registered virus go to ``_unmatched.gb`` for manual
review.
"""

from __future__ import annotations

import csv
import json
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

UNMATCHED = "_unmatched"


def load_ledger(path: Path) -> set:
    """Read a ledger of already-fetched accessions (one per line). Missing = empty."""
    path = Path(path)
    if not path.exists():
        return set()
    return {
        line.strip()
        for line in path.read_text().splitlines()
        if line.strip() and not line.startswith("#")
    }


def append_ledger(path: Path, accessions: List[str]) -> None:
    """Append newly-fetched accessions to the ledger, creating it if needed."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "a") as fh:
        for acc in accessions:
            fh.write(acc + "\n")


def slugify(name: str) -> str:
    """Filesystem-safe slug for a virus name."""
    slug = re.sub(r"[^A-Za-z0-9]+", "_", name.strip().lower())
    return slug.strip("_") or "virus"


def load_registry(registry_path: Path) -> List[Tuple[str, List[str]]]:
    """Return [(virus_name, [lowercased keywords]), ...] from the registry."""
    data = json.loads(Path(registry_path).read_text())
    out: List[Tuple[str, List[str]]] = []
    for v in data.get("viruses", []):
        name = v["virus_name"]
        kws = [k.lower() for k in v.get("keywords", [])]
        # Always allow the virus_name itself as a keyword.
        if name.lower() not in kws:
            kws.append(name.lower())
        # Longest keywords first so the most specific match wins.
        kws.sort(key=len, reverse=True)
        out.append((name, kws))
    return out


def _organism(record: SeqRecord) -> str:
    org = record.annotations.get("organism")
    if org:
        return org
    for feat in record.features:
        if feat.type == "source":
            vals = feat.qualifiers.get("organism")
            if vals:
                return vals[0]
    return ""


def match_virus(
    organism: str, registry: List[Tuple[str, List[str]]]
) -> Optional[str]:
    """Return the registered virus_name whose keyword is found in ``organism``."""
    org = organism.lower()
    for virus_name, keywords in registry:
        for kw in keywords:
            if kw and kw in org:
                return virus_name
    return None


def split_genbank(
    raw_path: Path,
    registry: List[Tuple[str, List[str]]],
    out_dir: Path,
    seen: Optional[set] = None,
) -> Dict[str, object]:
    """Split ``raw_path`` by species. Returns a summary dict.

    Writes ``<slug>.gb`` per matched virus, ``_unmatched.gb`` for the rest,
    and ``manifest.csv`` describing every record.

    Deduplication is by versioned accession (``record.id``, e.g. ``PP209408.1``):
    * records repeated **within this batch** are dropped (kept once);
    * records whose accession is already in ``seen`` (a ledger of previously
      fetched accessions) are dropped as cross-run duplicates.
    The manifest gains a ``status`` column = ``new`` | ``dup_in_batch`` |
    ``dup_in_ledger`` so every fetched record is accounted for.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    seen = set(seen or [])

    buckets: Dict[str, List[SeqRecord]] = {}
    manifest_rows: List[Dict[str, str]] = []
    counts: Dict[str, int] = {}
    new_accessions: List[str] = []
    batch_ids: set = set()
    dup_ledger = 0
    dup_batch = 0

    for record in SeqIO.parse(str(raw_path), "genbank"):
        acc = record.id
        organism = _organism(record)
        virus = match_virus(organism, registry)
        key = virus if virus else UNMATCHED
        out_name = f"{slugify(key)}.gb" if virus else f"{UNMATCHED}.gb"

        if acc in seen:
            status = "dup_in_ledger"
            dup_ledger += 1
        elif acc in batch_ids:
            status = "dup_in_batch"
            dup_batch += 1
        else:
            status = "new"
            batch_ids.add(acc)
            new_accessions.append(acc)
            buckets.setdefault(key, []).append(record)
            counts[key] = counts.get(key, 0) + 1

        manifest_rows.append(
            {
                "accession": acc,
                "organism": organism,
                "length": str(len(record.seq)),
                "matched_virus": virus or "",
                "output_file": out_name if status == "new" else "",
                "status": status,
            }
        )

    written_files: Dict[str, str] = {}
    for key, records in buckets.items():
        fname = f"{slugify(key)}.gb" if key != UNMATCHED else f"{UNMATCHED}.gb"
        fpath = out_dir / fname
        SeqIO.write(records, str(fpath), "genbank")
        written_files[key] = str(fpath)

    manifest_path = out_dir / "manifest.csv"
    with open(manifest_path, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["accession", "organism", "length",
                        "matched_virus", "output_file", "status"],
        )
        writer.writeheader()
        writer.writerows(manifest_rows)

    return {
        "total": len(manifest_rows),
        "new": len(new_accessions),
        "dup_in_batch": dup_batch,
        "dup_in_ledger": dup_ledger,
        "counts": counts,
        "files": written_files,
        "manifest": str(manifest_path),
        "new_accessions": new_accessions,
    }
