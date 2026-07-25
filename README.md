# 🧬 ViraLift

> Standardise viral gene names and assign annotation from a reference — and handle the NCBI data-collection step too.

ViraLift is a tool in a primer-design pipeline. It solves two classic problems with GenBank data: records that are **already annotated but have inconsistent gene names** (each lab names them differently), and records that are **not annotated at all**. The output is a set of genes with standardised names + coordinates, ready for downstream steps (e.g. extracting the ORF5 CDS for the whole input).

## Three modes

ViraLift is a Streamlit app with three modes, each with its own guide:

| Mode (UI) | Role | Guide |
|---|---|---|
| **Data crawler** | Collect data: fetch records from NCBI, split by species, dedup | [DATA_CRAWLER_GUIDE.md](DATA_CRAWLER_GUIDE.md) |
| **Run pipeline** | Core processing: standardise names + assign annotation, export TSV/FASTA | [PIPELINE_RUNNER_GUIDE.md](PIPELINE_RUNNER_GUIDE.md) |
| **Alias manager** | Manage the gene-name dictionary (canonical/alias) per virus | [ALIAS_MANAGER_GUIDE.md](ALIAS_MANAGER_GUIDE.md) |

Deep technical docs about the code: [CODEBASE_GUIDE.md](CODEBASE_GUIDE.md). The Data crawler engine: [gbcrawler/README.md](gbcrawler/README.md). Internal validation & paper notes (Vietnamese): [docs/](docs/README.md).

## The problem

GenBank records for the same gene are often labelled very differently, across many qualifier fields (`/gene`, `/product`, `/note`...):

| Virus | Same gene, many naming styles |
|---|---|
| PRRSV | `GP5` · `gp5` · `glycoprotein 5` · `major envelope glycoprotein` · `ORF5` |
| FMDV | `3Cpro` · `3C` · `protease 3C` · `3C protease` |

Depending on the record, this creates two situations:

```text
Case 1 — Query already annotated, but names non-standard
  -> Standardise names via alias lookup (fast, no alignment)

Case 2 — Query not annotated
  -> Lift coordinates from the reference via tblastn, then name
```

Run pipeline handles **both** automatically in a single run.

## Pipeline overview

The three modes form one seamless flow:

```text
        ┌──────────────┐
        │ Data crawler │   fetch from NCBI, split by species
        └──────┬───────┘
               │  <species>.gb files (full feature-bearing GenBank)
               v
        ┌──────────────┐      ┌───────────────┐
        │ Run pipeline │◄─────│ Alias manager │  gene-name dictionary
        └──────┬───────┘      └───────────────┘  (canonical/alias)
               │  TSV + FASTA (standardised names + coordinates)
               v
        downstream primer-pipeline steps (e.g. extract ORF5)
```

- **Data crawler** prepares clean input (correct format, split by species).
- **Run pipeline** is the core processing, using a reference + alias config.
- **Alias manager** maintains the gene-name dictionary that Run pipeline looks up; all three share one virus registry (`app/config/virus_alias_registry.json`).

## Install & run

### Option 1 — Docker (recommended)

```bash
cd viralift
cp .env.example .env  # optional: fill OPENAI_API_KEY to enable LLM alias review
docker compose up        # first build ~2 min
```

Open **http://localhost:8501**. Rebuild after code changes: `docker compose up --build`.

### Option 2 — Run locally

```bash
cd viralift
python -m venv .venv && source .venv/bin/activate
pip install -r ui/requirements.txt
# BLAST+ (tblastn, makeblastdb) must be installed via the OS package manager, not pip
streamlit run ui/streamlit_app.py
```

> **BLAST+** is only needed for the tblastn path (assigning annotation to records that lack it). The name-standardisation (alias) path does not need BLAST+.

### Optional LLM alias review

ViraLift can automatically review low-confidence bootstrap alias suggestions with an LLM. The core pipeline still works without this.

```bash
cp .env.example .env
# edit .env:
#   VIRALIFT_LLM_ENABLED=1
#   OPENAI_API_KEY=...
docker compose up --build
```

Only uncertain alias rows are sent for review; sequence data and full GenBank records are not included.

## Output

| File | Meaning |
|---|---|
| `extracted_cds.tsv` | Result table: standardised gene names + coordinates + mapping status |
| `extracted_cds.fasta` | Extracted gene sequences, filtered by coverage/status |
| `manifest.csv` (Data crawler) | Table of fetched records: accession, organism, status |

## Directory structure

```text
viralift/
├── app/                # pipeline engine (parse, lift, alias, io)
│   └── config/         # registry + per-virus alias config
├── gbcrawler/          # Data crawler engine (NCBI fetch + split)
├── ui/                 # Streamlit app (3 modes)
├── docs/               # internal validation & paper notes (Vietnamese)
├── README.md           # this file
├── DATA_CRAWLER_GUIDE.md
├── PIPELINE_RUNNER_GUIDE.md
├── ALIAS_MANAGER_GUIDE.md
└── CODEBASE_GUIDE.md
```
