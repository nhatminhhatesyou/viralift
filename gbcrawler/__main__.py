"""gbcrawler CLI.

Examples
--------
Crawl by NCBI search query, split by species:

    python -m gbcrawler \
        --query "txid28344[Organism:exp] AND complete genome" \
        --email you@lab.org --api-key $NCBI_API_KEY \
        --out crawl_out/

Crawl an explicit accession list:

    python -m gbcrawler --accessions accessions.txt \
        --email you@lab.org --out crawl_out/

Re-split an already-downloaded combined GenBank file (offline, no network):

    python -m gbcrawler --from-raw combined.gb --out crawl_out/

Output (per run, in --out):
    <virus_slug>.gb   one file per matched species  -> feed to ViraLift
    _unmatched.gb     records not matching any registered virus
    raw_combined.gb   everything downloaded (unless --from-raw)
    manifest.csv      accession, organism, length, matched_virus, output_file
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

# Work both as a package (`python -m gbcrawler` from viralift/) and as a direct
# script (`python __main__.py` from inside gbcrawler/).
try:
    from gbcrawler import fetch as fetch_mod
    from gbcrawler import split as split_mod
except ImportError:
    import fetch as fetch_mod
    import split as split_mod

# gbcrawler is a submodule of ViraLift (viralift/gbcrawler/), so the registry
# sits at ../app/config/ relative to this file. Override with --registry.
_DEFAULT_REGISTRY = (
    Path(__file__).resolve().parent.parent
    / "app" / "config" / "virus_alias_registry.json"
)


def parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="gbcrawler",
        description="Crawl NCBI GenBank records and split them by species for ViraLift.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    src = p.add_argument_group("input (at least one of --query / --accessions / --from-raw)")
    src.add_argument("--query", help="NCBI search expression (taxon + filters).")
    src.add_argument("--accessions", help="File of accessions (one per line, commas ok).")
    src.add_argument("--from-raw", help="Skip fetching; split an existing combined .gb.")

    net = p.add_argument_group("NCBI / fetch options")
    net.add_argument("--email", default=os.environ.get("NCBI_EMAIL"),
                     help="Contact email required by NCBI (or set NCBI_EMAIL).")
    net.add_argument("--api-key", default=os.environ.get("NCBI_API_KEY"),
                     help="NCBI API key to raise the rate limit (or set NCBI_API_KEY).")
    net.add_argument("--db", default="nuccore", help="Entrez database. Default: nuccore")
    net.add_argument("--retmax", type=int, default=10000,
                     help="Max records to fetch for a query. Default: 10000")
    net.add_argument("--batch", type=int, default=200,
                     help="Records per efetch request. Default: 200")
    net.add_argument("--count", action="store_true",
                     help="Only print how many records the --query matches, then "
                          "exit. Nothing is downloaded. Use it to check a query first.")

    out = p.add_argument_group("output")
    out.add_argument("--out", help="Output directory (required unless --count).")
    out.add_argument("--registry", default=str(_DEFAULT_REGISTRY),
                     help="ViraLift virus alias registry JSON.")
    out.add_argument("--ledger",
                     help="Accession ledger file. Records already listed here are "
                          "skipped as duplicates; new ones are appended after the run. "
                          "Use the same ledger across runs to avoid re-fetching.")
    out.add_argument("--quiet", action="store_true", help="Reduce console output.")
    return p.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)

    if not (args.query or args.accessions or args.from_raw):
        print("error: provide --query, --accessions, or --from-raw.", file=sys.stderr)
        return 2

    # Preview mode: just count how many records the query matches, then stop.
    if args.count:
        if not args.query:
            print("error: --count needs --query.", file=sys.stderr)
            return 2
        try:
            fetch_mod.configure(args.email, args.api_key)
        except ValueError as e:
            print(f"error: {e}", file=sys.stderr)
            return 2
        n = fetch_mod.count_query(args.query, db=args.db, api_key=args.api_key)
        print(f"[gbcrawler] query matches {n} records in {args.db} (nothing downloaded).")
        return 0

    if not args.out:
        print("error: --out is required (unless using --count).", file=sys.stderr)
        return 2
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    # 1. Obtain a combined GenBank file (download, or reuse an existing one).
    if args.from_raw:
        raw_path = Path(args.from_raw)
        if not raw_path.exists():
            print(f"error: --from-raw file not found: {raw_path}", file=sys.stderr)
            return 2
    else:
        try:
            fetch_mod.configure(args.email, args.api_key)
        except ValueError as e:
            print(f"error: {e}", file=sys.stderr)
            return 2

        raw_path = out_dir / "raw_combined.gb"
        total = 0
        with open(raw_path, "w") as raw:
            if args.query:
                total += fetch_mod.fetch_from_query(
                    args.query, raw, db=args.db, retmax=args.retmax,
                    batch=args.batch, api_key=args.api_key, quiet=args.quiet,
                )
            if args.accessions:
                accs = fetch_mod.read_accession_file(Path(args.accessions))
                total += fetch_mod.fetch_from_accessions(
                    accs, raw, db=args.db, batch=args.batch,
                    api_key=args.api_key, quiet=args.quiet,
                )
        if total == 0:
            print("warning: no records fetched.", file=sys.stderr)

    # 2. Split by species against the ViraLift registry (dedup vs. ledger).
    registry = split_mod.load_registry(Path(args.registry))
    seen = split_mod.load_ledger(Path(args.ledger)) if args.ledger else set()
    summary = split_mod.split_genbank(raw_path, registry, out_dir, seen=seen)

    if args.ledger:
        split_mod.append_ledger(Path(args.ledger), summary["new_accessions"])

    # 3. Report.
    print(f"\n[gbcrawler] {summary['total']} records fetched -> {out_dir}")
    print(f"  new: {summary['new']}   "
          f"dup_in_batch: {summary['dup_in_batch']}   "
          f"dup_in_ledger: {summary['dup_in_ledger']}")
    for key, count in sorted(summary["counts"].items()):
        label = key if key != split_mod.UNMATCHED else "(unmatched)"
        print(f"  {count:>5}  {label}  ->  {Path(summary['files'][key]).name}")
    print(f"  manifest: {summary['manifest']}")
    if args.ledger:
        print(f"  ledger:   {args.ledger}  (+{summary['new']} accessions)")
    print("\nNext: feed each species file to ViraLift, e.g.")
    print("  python -m app.src.main --reference <ref.gb> "
          f"--query {out_dir}/<virus>.gb --output output/run")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
