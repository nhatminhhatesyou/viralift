# 🌐 gbcrawler

> NCBI GenBank crawler — the **input layer** for ViraLift.
> Crawl sequences from NCBI nuccore, download them **with features**
> (`gbwithparts`), and split them into one `.gb` file per virus species,
> ready to drop into ViraLift.

```
NCBI (search query and/or accession list)
        │   efetch rettype=gbwithparts  (features preserved)
        ▼
   raw_combined.gb
        │   split by organism vs. ViraLift registry
        ▼
   prrsv.gb   fmdv.gb   porcine_epidemic_diarrhea_virus.gb   _unmatched.gb
        │
        ▼   each species + its reference
     ViraLift  →  extracted_cds.tsv / .fasta
```

## Why it exists

ViraLift needs GenBank input **with feature tables** (it reads `/gene`,
`/product`, `/note` for the alias path, and the sequence for the tblastn path).
A plain FASTA dump breaks both. gbcrawler always fetches `gbwithparts` and
splits by species, because ViraLift takes **one reference per species**.

## Install

```bash
pip install biopython
```

> **Layout** — gbcrawler is a submodule of ViraLift: `viralift/gbcrawler/`.
> Run all commands below from the `viralift/` folder.

## Web UI

The crawler is built into the ViraLift Streamlit app as the **"Data crawler"**
mode (see `viralift/ui/data_crawler_page.py`) — launch the app and switch to
that mode. This module (`gbcrawler/`) is the engine behind both the UI and the
CLI below.

## Usage (CLI)

Two input modes (use either or both in one run):

**By NCBI search query** (same idea as filtering on NCBI Virus by hand):

```bash
python -m gbcrawler \
    --query "txid28344[Organism:exp] AND complete genome" \
    --email you@lab.org --api-key $NCBI_API_KEY \
    --out crawl_out/
```

**By accession list** (one accession per line; commas allowed):

```bash
python -m gbcrawler --accessions accessions.txt \
    --email you@lab.org --out crawl_out/
```

**Re-split an already-downloaded file** (offline, no network):

```bash
python -m gbcrawler --from-raw raw_combined.gb --out crawl_out/
```

> NCBI requires a contact email (`--email` or `NCBI_EMAIL`). An API key
> (`--api-key` / `NCBI_API_KEY`) raises the rate limit from 3 to 10 req/s —
> recommended for large pulls. Get one in your NCBI account settings.

## Output (in `--out`)

| File | Contents |
|------|----------|
| `<virus_slug>.gb` | one file per **matched** species → feed to ViraLift |
| `_unmatched.gb` | records whose organism matched no registered virus |
| `raw_combined.gb` | everything downloaded (omitted when `--from-raw`) |
| `manifest.csv` | `accession, organism, length, matched_virus, output_file, status` |

## Controlling which sequences you get (dedup)

Every record is identified by its **versioned accession** (`PP209408.1`). The
`manifest.csv` `status` column tells you exactly what happened to each one:

| status | meaning |
|--------|---------|
| `new` | first time seen, written to a species file |
| `dup_in_batch` | the same accession appeared twice in this run (kept once) |
| `dup_in_ledger` | already fetched in a previous run (skipped) |

To make runs aware of each other, pass the **same** `--ledger` file each time.
It's a plain text list of accessions already fetched; gbcrawler skips anything
in it and appends the new ones afterwards:

```bash
# first pull
python -m gbcrawler --query '"...virus"[Organism]' --retmax 100 \
    --email you@lab.org --out run1/ --ledger prrsv_ledger.txt

# later / broader pull — overlapping records are skipped automatically
python -m gbcrawler --query '"...virus"[Organism]' --retmax 300 \
    --email you@lab.org --out run2/ --ledger prrsv_ledger.txt
```

The run summary prints `new / dup_in_batch / dup_in_ledger` counts. Without
`--ledger`, dedup still removes duplicates *within* a single run, but separate
runs won't know about each other — compare their `manifest.csv` accessions, or
adopt a ledger.

Species matching uses ViraLift's `app/config/virus_alias_registry.json`
(override with `--registry`). To support a new virus, add it to that registry
(name + keywords + alias config) — gbcrawler and ViraLift then stay in sync.

## Hand-off to ViraLift

Both run from the `viralift/` folder:

```bash
# crawl
python -m gbcrawler --from-raw combined.gb --out crawl_out/

# then standardise / extract with ViraLift
python -m app.src.main \
    --reference app/data/PRRS/PRRS_ref_test.gb \
    --query     crawl_out/prrsv.gb \
    --output    output/prrsv_run \
    --alias-config app/config/prrsv_alias.json
```

## Options

| Flag | Default | Meaning |
|------|---------|---------|
| `--query` | – | NCBI search expression |
| `--accessions` | – | file of accessions |
| `--from-raw` | – | skip fetching, split an existing `.gb` |
| `--db` | `nuccore` | Entrez database |
| `--retmax` | `10000` | max records per query |
| `--batch` | `200` | records per efetch request |
| `--email` / `--api-key` | env | NCBI identity / rate-limit key |
| `--registry` | ViraLift registry | species → file mapping |
| `--out` | *required* | output directory |
| `--quiet` | off | reduce console output |
