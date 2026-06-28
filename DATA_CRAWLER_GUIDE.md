# ViraLift Data Crawler Guide

> Navigation: [README](README.md) · [Pipeline Runner](PIPELINE_RUNNER_GUIDE.md) · **Data Crawler** · [Alias Manager](ALIAS_MANAGER_GUIDE.md)

This document explains ViraLift's data-collection stage: from the moment the user describes their criteria (virus, length, date...) until they get back GenBank files already split by species, ready to feed into the pipeline.

In the current UI this part is called **Data crawler**. The underlying engine is the `gbcrawler/` module (shared by both the Web UI and the CLI).

## Table of contents

- [What is the Data Crawler?](#what-is-the-data-crawler)
- [Why do you need the Data Crawler?](#why-do-you-need-the-data-crawler)
- [Input and output](#input-and-output)
- [Overall flow](#overall-flow)
- [Steps in the Web UI](#steps-in-the-web-ui)
- [NCBI query syntax](#ncbi-query-syntax)
- [How to tell if a query is right](#how-to-tell-if-a-query-is-right)
- [Splitting by species (registry)](#splitting-by-species-registry)
- [Deduplication: ledger and status](#deduplication-ledger-and-status)
- [Filtering by region and country](#filtering-by-region-and-country)
- [Handing output to the Pipeline Runner](#handing-output-to-the-pipeline-runner)
- [Running from the CLI](#running-from-the-cli)
- [Common cases](#common-cases)
- [Best practices](#best-practices)

## What is the Data Crawler?

The Data Crawler is ViraLift's **data preparation** stage. It takes:

```text
NCBI search criteria (virus + filters), OR a list of accessions
```

and produces output consisting of:

```text
GenBank files pre-split by species + a manifest + a dedup ledger
```

In short: the Data Crawler automates what a lab researcher does by hand on **NCBI Virus / Nucleotide** — pick a virus, filter, download — but adds: automatic species splitting, deduplication across runs, and a direct hand-off into the pipeline.

## Why do you need the Data Crawler?

Raw GenBank records have two problems that the Pipeline Runner needs to handle:

```text
1. Already annotated but gene names are inconsistent (each lab names them differently)
2. Not annotated at all
```

The Pipeline Runner solves both, but **its input must be GenBank with features** (the `gbwithparts` format), and **each species needs its own reference**. The Data Crawler is the layer in front that takes care of exactly these two things:

- Always downloads the full feature-bearing GenBank format (not bare FASTA) so the Pipeline Runner can read `/gene`, `/product`, `/note` and the sequence.
- Automatically splits records by species, one file per species, so they can be paired with the right reference + alias config.

That is why the correct flow is not `NCBI -> Pipeline`, but `NCBI -> Data Crawler -> Pipeline`.

## Input and output

### Input

| Input | Meaning |
|---|---|
| NCBI query | Criteria built from fields: virus, dataset type, length, date, region, extra condition |
| Accession list | A ready-made list of accessions (`.txt/.csv` file or pasted) — instead of a query |
| Email | NCBI requires a contact identity to use Entrez |
| API key | Optional; raises the rate limit 3 -> 10 requests/second |
| retmax | Cap on the number of records to fetch for a query |
| Registry | JSON file deciding which species are recognised and split |
| Ledger | A record of already-fetched accessions, used to dedup across runs |

### Output

By default placed in `output/gbcrawler/`:

| Output | Meaning |
|---|---|
| `raw_combined.gb` | The raw bundle of all downloaded records, before splitting / dedup |
| `<species>.gb` | One file per recognised species (e.g. `prrsv.gb`) — this is what goes into the pipeline |
| `_unmatched.gb` | Records matching no species in the registry, for manual review |
| `manifest.csv` | A table of every record: accession, organism, length, matched_virus, output_file, status |

## Overall flow

```text
User describes criteria (virus + filters)  OR  provides an accession list
        |
        v
Build the NCBI query (assemble fields into an Entrez string)
        |
        v
(optional) Count: ask NCBI how many records match, download nothing
        |
        v
Fetch: download full feature-bearing GenBank (gbwithparts) -> raw_combined.gb
        |
        v
Split: read each record's organism, match the registry -> group by species
        |
        +-- duplicate accession  -> drop (mark in status)
        |
        v
Write <species>.gb + _unmatched.gb + manifest.csv, update the ledger
        |
        v
Pick a species file -> "Use as query" -> switch to the Pipeline Runner
```

## Steps in the Web UI

Open the ViraLift app and click the **Data crawler** mode in the sidebar. In this mode, the **Email** and **API key** fields appear right in the sidebar (fill once, used for both Count and Crawl).

### 1. Build input dataset

Choose the input source: **By NCBI query** or **By accession list**.

With **By NCBI query**, fill in the fields — the app assembles the query and shows it in the **NCBI query preview** box:

| Field | Effect | Generates |
|---|---|---|
| Virus | Pick a species from the registry, or "Other / custom" to type your own | `"<organism>"[Organism]` |
| Dataset type | Complete genomes / Gene records / fragments / Both | `complete genome[Title]` or `NOT complete genome[Title]` |
| Length range | bp length range (optional) | `("lo"[SLEN] : "hi"[SLEN])` |
| Publication date | Publication date range (optional) | `("d0"[PDAT] : "d1"[PDAT])` |
| Country | Any (global) / Vietnam / Custom country | `"<country>"[All Fields]` |
| Extra NCBI condition | Free-form Entrez syntax | appended directly to the query |

There is an **Advanced: edit query manually** expander to override the query when needed.

With **By accession list**: upload a `.txt/.csv` file or paste accessions (one per line, commas allowed). The app reports how many accessions it detected.

Finally set **Max records to fetch** (retmax) — the cap on how many records to download.

The **Advanced NCBI settings** expander holds: Database (default `nuccore`), Output folder, the dedup-ledger toggle and ledger path.

### 2. Preview and crawl

- **Count records** — asks NCBI how many records match the query, **downloading nothing**. Use it to check the query first. (Disabled in accession mode.)
- **Crawl GenBank and split** — actually downloads, splits by species, updates the ledger. Has a progress bar.

### 3. Use crawler output

After crawling, this section shows:

- 4 metrics: **Total fetched / New records / Duplicate in batch / Duplicate in ledger**.
- A list of **Species files**: each species with a **Download .gb** button and a **Use as query** button (sends that file straight into the Pipeline Runner as the query input).
- A **manifest.csv** expander: the detailed table, with a download button.

## NCBI query syntax

"Filtering" in the Data Crawler really means writing the query in the **NCBI Entrez search-field syntax**. Each condition has the form `value[FIELD_NAME]`, joined by `AND` / `OR` / `NOT`. The commonly used fields (database `nuccore`):

| Field | Meaning | Example |
|---|---|---|
| `[Organism]` | Pick the species (most important) | `"Porcine reproductive and respiratory syndrome virus"[Organism]` |
| `[Title]` | Words in the record title | `complete genome[Title]`, `(ORF5[Title] OR GP5[Title])` |
| `[SLEN]` | Sequence length | `("14000"[SLEN] : "16000"[SLEN])` |
| `[PDAT]` | Publication date | `("2020/01/01"[PDAT] : "2025/12/31"[PDAT])` |
| `NOT` | Exclude | `NOT partial[Title]` |

Example of a full query (PRRSV, complete genomes, 14–16 kb, from 2020):

```text
"Porcine reproductive and respiratory syndrome virus"[Organism]
  AND complete genome[Title]
  AND ("14000"[SLEN] : "16000"[SLEN])
  AND ("2020/01/01"[PDAT] : "2025/12/31"[PDAT])
```

## How to tell if a query is right

There are three ways; use them in combination:

1. **Count records (in the UI)** — shows how many records match without downloading. Reading the number:
   - `0` -> the query is wrong or too strict (misspelled organism, wrong field).
   - very large (tens of thousands or more) -> the query is too broad (usually a missing `[Organism]`).
   - a few hundred to a few thousand -> reasonable.
2. **Cross-check on the NCBI website** — paste the exact query string into the search box at [NCBI Nucleotide](https://www.ncbi.nlm.nih.gov/nuccore); the syntax is identical. If the result count matches Count, the query is sound.
3. **Inspect `manifest.csv` after downloading** — the `organism` column shows whether you downloaded the right species; off records land in `_unmatched.gb`. This is the safety net: even if a loose query lets in another species, the tool separates it rather than mixing it into the main species file.

> Note: if you type only `prrs` **without** attaching `[Organism]`, NCBI searches the text "prrs" across all fields (`[All Fields]`) — letting in junk and missing records that use the full name. Always pick the virus from the dropdown (which attaches `[Organism]` for you).

## Splitting by species (registry)

After downloading, the tool reads each record's `organism` field and matches it against `app/config/virus_alias_registry.json` — the same file the Pipeline Runner and Alias Manager use. A record matching a species' keyword goes into that species file; non-matches go into `_unmatched.gb`.

Because they share one registry, to make the Data Crawler recognise a new virus you only need to add it to the registry (see the [Alias Manager](ALIAS_MANAGER_GUIDE.md)). The Data Crawler and the rest of ViraLift always stay in sync.

## Deduplication: ledger and status

Each record is identified by its **versioned accession** (e.g. `PP209408.1`). The tool dedups at two layers:

- **Within a single run**: records repeating an accession are dropped, keeping only one.
- **Across runs**: if you share one **ledger** file, accessions fetched in a previous run are skipped, and new accessions are appended to the ledger.

The `status` column in `manifest.csv` records the fate of each record:

| status | Meaning |
|---|---|
| `new` | Seen for the first time; written to a species file |
| `dup_in_batch` | Duplicate accession within this run (kept once) |
| `dup_in_ledger` | Already fetched in a previous run (skipped) |

So if you run `retmax 100` today and tomorrow expand to `retmax 300` with the same ledger, the 100 old records are skipped automatically and only the new ones are downloaded.

> Identification is by version: if GenBank updates a record (`.1` -> `.2`), the new version counts as a new record.

## Filtering by region and country

This is a caveat. The Data Crawler runs on the `nuccore` database, which **has no dedicated index field for country/region**. The Country field in the UI assembles to `"<country>"[All Fields]` — a text match only, so it is **approximate**:

- For a specific country (e.g. `Vietnam`) it is acceptable.
- You cannot reliably filter by "continent" (e.g. "Asia"), because records usually record a specific country name, not "Asia".

To build accurate region subsets (global / Asia / Vietnam) for analysis, two more reliable approaches:

1. Use **NCBI Virus** (web) to filter by region — it parses geo natively — then export an accession list and feed it into the Data Crawler via accession mode.
2. Download broadly (global), then **filter by the `/country` qualifier** from the downloaded records themselves.

## Handing output to the Pipeline Runner

In the results section, each species file has a **Use as query** button. Clicking it:

```text
Sets that species file as the query input
        |
        v
Switches the app to Run pipeline mode (Upload stage)
        |
        v
You just add that species' reference and run
```

This is how you connect Data Crawler -> Pipeline Runner inside one app, without exporting a file and re-uploading it. For the next stage, see the [Pipeline Runner](PIPELINE_RUNNER_GUIDE.md).

## Running from the CLI

Besides the Web UI, the `gbcrawler` engine runs from the command line (run from the `viralift/` folder):

Count only (no download):

```bash
python -m gbcrawler --count \
  --query '"Porcine reproductive and respiratory syndrome virus"[Organism] AND complete genome[Title]' \
  --email you@lab.org
```

Download + split by species, with a ledger:

```bash
python -m gbcrawler \
  --query '"Porcine reproductive and respiratory syndrome virus"[Organism] AND complete genome[Title]' \
  --retmax 1000 \
  --email you@lab.org --api-key $NCBI_API_KEY \
  --out output/gbcrawler/ --ledger output/gbcrawler_ledger.txt
```

By accession list:

```bash
python -m gbcrawler --accessions accessions.txt \
  --email you@lab.org --out output/gbcrawler/
```

Then send a species file to the Pipeline Runner:

```bash
python -m app.src.main \
  --reference app/data/PRRS/PRRS_ref_test.gb \
  --query     output/gbcrawler/prrsv.gb \
  --output    output/prrsv_run \
  --alias-config app/config/prrsv_alias.json
```

For the full flag list see `gbcrawler/README.md`.

## Common cases

### Want the complete genomes of one species

Pick the virus + Dataset type = **Complete genomes**. Click **Count** to see the number, then **Crawl**.

### Want individual gene records (e.g. ORF5), not whole genomes

Pick Dataset type = **Gene records / fragments**, or add the condition `(ORF5[Title] OR GP5[Title])` in the Extra box. Note that individual gene records often lack annotation -> the Pipeline Runner will route them through tblastn to add annotation (needs BLAST+).

### Fetch more without duplicating the previous batch

Use **the same ledger** for every run of the same project. The next run downloads only the new part; the `status` column shows how many are `dup_in_ledger`.

### Some records go to `_unmatched.gb`

It means their organism matched no species keyword in the registry. Either the query let in an unrelated species, or the registry does not have that species yet. Add the species to the registry (see the Alias Manager) and re-run if needed.

### Count returns 0

The query is wrong or too strict. Re-check the organism name, drop some filters, or cross-check directly on the NCBI website.

## Best practices

- **Always Count before you Crawl.** Confirm the number is reasonable before downloading, to avoid fetching the wrong set or too much.
- **Pick the virus from the dropdown** instead of typing, so the query always has a correct `[Organism]`.
- **Use one shared ledger per project** so batches don't overlap.
- **Register an NCBI API key** for large downloads (raises the rate from 3 -> 10 req/s).
- **Don't fully trust region filtering** on nuccore; for precise geographic subsets, filter `/country` afterward or use NCBI Virus.
- **Keep `raw_combined.gb`** for cross-checking / re-downloading; feed only the `<species>.gb` files into the pipeline.
- **Inspect `manifest.csv`** after each run to confirm you fetched the right species and to gauge the duplicate rate.
