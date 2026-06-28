"""NCBI Entrez fetch layer for gbcrawler.

Two entry points cover both supported input modes:

* :func:`fetch_from_query`      — an NCBI search expression (taxon + filters)
* :func:`fetch_from_accessions` — an explicit list of accessions

Both stream raw GenBank text (``rettype=gbwithparts``) into a single file so
that feature tables are preserved for ViraLift's alias / tblastn routing.

NCBI etiquette is enforced here:
* ``Entrez.email`` is mandatory; an API key lifts the rate limit 3/s -> 10/s.
* requests are batched and throttled with a per-request delay.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path
from typing import List, Optional, TextIO

from Bio import Entrez


# NCBI rate limits: 3 req/s without an API key, 10 req/s with one.
# Add a small safety margin so we never trip the throttle.
_DELAY_NO_KEY = 1.0 / 3.0 + 0.05
_DELAY_KEY = 1.0 / 10.0 + 0.02


def configure(email: str, api_key: Optional[str] = None) -> None:
    """Set the global Entrez identity. ``email`` is required by NCBI."""
    if not email:
        raise ValueError(
            "NCBI requires an email address. Pass --email or set NCBI_EMAIL."
        )
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key


def _delay(api_key: Optional[str]) -> float:
    return _DELAY_KEY if api_key else _DELAY_NO_KEY


def _log(msg: str, quiet: bool) -> None:
    if not quiet:
        print(msg, file=sys.stderr, flush=True)


def count_query(
    query: str,
    *,
    db: str = "nuccore",
    api_key: Optional[str] = None,
) -> int:
    """Return how many records match ``query`` without downloading any.

    Use this to sanity-check a query before a real crawl: a count of 0 means the
    query is wrong/too strict; a huge count means it's too broad.
    """
    handle = Entrez.esearch(db=db, term=query, retmax=0)
    result = Entrez.read(handle)
    handle.close()
    return int(result["Count"])


def fetch_from_query(
    query: str,
    out_handle: TextIO,
    *,
    db: str = "nuccore",
    retmax: int = 10000,
    batch: int = 200,
    api_key: Optional[str] = None,
    quiet: bool = False,
) -> int:
    """Search NCBI and stream matching GenBank records into ``out_handle``.

    Uses the Entrez history server so large result sets are fetched without
    re-sending a giant id list. Returns the number of records written.
    """
    handle = Entrez.esearch(db=db, term=query, retmax=0, usehistory="y")
    result = Entrez.read(handle)
    handle.close()

    total = int(result["Count"])
    webenv = result["WebEnv"]
    query_key = result["QueryKey"]
    n = min(total, retmax)
    _log(f"[gbcrawler] query matched {total} records; fetching {n}.", quiet)

    written = 0
    for start in range(0, n, batch):
        chunk = min(batch, n - start)
        time.sleep(_delay(api_key))
        fetch = Entrez.efetch(
            db=db,
            rettype="gbwithparts",
            retmode="text",
            retstart=start,
            retmax=chunk,
            webenv=webenv,
            query_key=query_key,
        )
        out_handle.write(fetch.read())
        fetch.close()
        written += chunk
        _log(f"[gbcrawler]   fetched {written}/{n}", quiet)
    return written


def fetch_from_accessions(
    accessions: List[str],
    out_handle: TextIO,
    *,
    db: str = "nuccore",
    batch: int = 200,
    api_key: Optional[str] = None,
    quiet: bool = False,
) -> int:
    """Stream GenBank records for an explicit accession list into ``out_handle``."""
    accessions = [a.strip() for a in accessions if a.strip()]
    n = len(accessions)
    _log(f"[gbcrawler] fetching {n} accessions.", quiet)

    written = 0
    for start in range(0, n, batch):
        ids = accessions[start : start + batch]
        time.sleep(_delay(api_key))
        fetch = Entrez.efetch(
            db=db,
            id=",".join(ids),
            rettype="gbwithparts",
            retmode="text",
        )
        out_handle.write(fetch.read())
        fetch.close()
        written += len(ids)
        _log(f"[gbcrawler]   fetched {written}/{n}", quiet)
    return written


def read_accession_file(path: Path) -> List[str]:
    """Read accessions from a .txt/.csv (one per line, commas allowed)."""
    text = Path(path).read_text()
    out: List[str] = []
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        out.extend(p.strip() for p in line.replace(",", " ").split())
    return out
