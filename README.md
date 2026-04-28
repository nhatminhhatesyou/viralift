# 🧬 ViraLift

> **Reference-guided viral gene name standardisation** — lift CDS coordinates from a well-annotated reference onto unannotated (or poorly annotated) genomes, then normalise gene names to a canonical vocabulary.

---

## 🧩 The Problem

GenBank viral genome submissions are inconsistent. The same gene across different submissions can be labelled:

| Virus | Same gene, different labels |
|-------|-----------------------------|
| PRRSV | `GP5` · `gp5` · `glycoprotein 5` · `major envelope glycoprotein` · `ORF5` |
| FMDV  | `3Cpro` · `3C` · `protease 3C` · `3C protease` |

This creates two distinct problems depending on the submission:

```
┌─────────────────────────────────────────────────────────────┐
│  Case 1 — Query is already annotated                        │
│  Gene names exist but are non-standard                      │
│  → Normalise names via alias lookup (fast, no alignment)    │
├─────────────────────────────────────────────────────────────┤
│  Case 2 — Query has no annotation                           │
│  No CDS features at all                                     │
│  → Lift coordinates from reference via tblastn, then name   │
└─────────────────────────────────────────────────────────────┘
```

ViraLift handles **both cases automatically** in a single run.

---

## ⚙️ How It Works

```
Reference GenBank (.gb)          Query GenBank (.gb)
   [well-annotated]                [1 – N records]
         │                               │
         ▼                               │
  Parse CDS features                     │
  Translate → protein                    │
         │                               │
         ▼                               ▼
  ┌──────────────────────────────────────────────────┐
  │          Smart Router (per query record)         │
  │                                                  │
  │  Already annotated + matching granularity?       │
  │  ├── YES → Direct alias normalisation            │
  │  └── NO  → tblastn protein-guided lifting        │
  └──────────────────────────────────────────────────┘
         │
         ▼
  Alias normalisation
  (raw name → canonical name via config)
         │
         ▼
  TSV / FASTA output  (coordinates + canonical names + status)
```

### 🔬 Why tblastn?

Proteins are ~3–4× more conserved than nucleotides across serotypes and lineages. `tblastn` searches a protein query against a nucleotide subject, making it robust even for divergent strains. All reference proteins are searched in **a single batched call** per genome.

---

## 🖥️ Web UI (recommended)

ViraLift ships with a **Streamlit web interface** served via Docker. No Python setup needed on your machine.

### Start the UI

```bash
# Clone the repo
git clone <repo-url>
cd viralift

# Start (builds image on first run, ~2 min)
docker compose up
```

Open **http://localhost:8501** in your browser.

To rebuild after code changes:

```bash
docker compose up --build
```

To run in the background:

```bash
docker compose up -d
docker compose down   # stop
```

---

### UI Walkthrough

#### Stage 1 — Upload 📂

1. Upload a **Reference GenBank** file (single well-annotated record) on the left.
2. Upload a **Query GenBank** file (one or more records) on the right.
3. Optionally expand **Advanced options** to tune lifting thresholds:
   - **Min coverage** — minimum fraction of the reference protein that must be covered by a tblastn hit (default `0.5`)
   - **Min identity** — minimum protein identity to accept a hit (default `0.3`)
   - **E-value** — tblastn significance threshold (default `1e-5`)
   - **Rescue window** — how many bp to scan around the lifted start position looking for ATG if start codon is missing (default `50`)
4. Toggle **"Use ref gene names as output"** if you want the output to use the reference's original gene names (e.g. `Lab`) instead of the canonical alias key (e.g. `Lpro`).
5. Click **▶ Run ViraLift**.

> The alias config for the detected virus is auto-selected from the registry. You don't need to specify anything manually.

---

#### Stage 2 — Resolve ⚠️ *(shown only when unknown names are detected)*

**Ref-side panel** (top) — shown if the reference has genes not yet in the alias database:

- Lists each unrecognised ref gene name.
- Check **➕ Add as new canonical** for any names you want to register.
- Click **💾 Save selected ref names** to write them permanently to the alias config.
- You can skip this and continue — lifting still works regardless (tblastn uses protein sequence, not name).

**Query-side resolver** (below) — shown if query records contain names not in the alias database:

| Column | What to do |
|--------|-----------|
| Gene name | The raw name found in the query file |
| Dropdown | Pick the canonical name to map it to, or leave as `-- ignore (keep raw name) --` |
| 💾 Save | Check to permanently add this mapping to the alias config for future runs |

Click **▶ Continue with these decisions** when done.

> Decisions are logged to `logs/viralift.log` even if you don't save them permanently — so you always have a trace of what was mapped in each run.

---

#### Stage 3 — Results 📊

**Summary badges** at the top show the status distribution across all records:

| Badge | Meaning |
|-------|---------|
| 🟢 OK | Lifted successfully, valid start + stop codon |
| 🟡 Rescued | Start codon was missing but recovered nearby |
| 🟠 Invalid boundary | Could not fix start or stop codon |
| 🟠 Low coverage | tblastn hit found but protein coverage below threshold |
| 🔴 No hit | tblastn found no alignment for this gene |
| 🔴 Translation fail | Reference feature could not be translated to protein |

**Per-record expanders** show a table for each query record with:
- Status, canonical name, raw original name, start/end coordinates, coverage, identity, method (direct/tblastn)

---

#### Export ⬇️

**TSV tab:**
- Download the full results table with canonical names or raw names.

**FASTA extraction tab:**
1. Select which genes to extract using the multiselect.
2. Choose output format: one FASTA per gene, or all genes in a single file.
3. Set quality filters (min coverage, min identity, include/exclude rescued).
4. Click **⬇ Generate & Download FASTA**.

> FASTA headers follow the format: `>{record_id}|{gene}|{start}|{end}|{strand}`

---

### Data Persistence

Three directories are volume-mounted so data survives Docker rebuilds:

| Directory | What persists |
|-----------|--------------|
| `app/config/` | Alias config files — new mappings saved via the UI |
| `logs/` | Run history and alias audit trail |
| `output/` | Any files you generate outside the UI |

---

## 💻 CLI

For scripting or batch runs without the UI.

### Installation

**Requirements:** Python ≥ 3.10, BLAST+ (`tblastn` in PATH)

```bash
git clone <repo-url>
cd viralift
python3 -m venv venv
source venv/bin/activate
pip install -r app/requirements.txt
```

Verify:

```bash
tblastn -version
```

### Examples

**Minimal run — auto-detect alias config**
```bash
python -m app.src.main \
  --reference app/data/PRRS_ref_test.gb \
  --query     app/data/PRRSV_test.gb \
  --output    output/prrs_run
```

**Explicit alias config**
```bash
python -m app.src.main \
  --reference    app/data/PRRS_ref_test.gb \
  --query        app/data/PRRSV_test.gb \
  --output       output/prrs_run \
  --alias-config app/config/prrsv_alias.json
```

**FMDV example**
```bash
python -m app.src.main \
  --reference    app/data/FMD_ref_test.gb \
  --query        app/data/FMD_test.gb \
  --output       output/fmd_run \
  --alias-config app/config/fmdv_alias.json
```

**Strict thresholds**
```bash
python -m app.src.main \
  --reference    app/data/PRRS_ref_test.gb \
  --query        app/data/PRRSV_test.gb \
  --output       output/prrs_strict \
  --min-coverage 0.7 \
  --min-identity 0.5 \
  --evalue       1e-10
```

### All Options

| Flag | Default | Description |
|------|---------|-------------|
| `--reference` | *(required)* | Reference GenBank file (single annotated record) |
| `--query` | *(required)* | Query GenBank file (one or more records) |
| `--output` | `output/run` | Output directory (created if missing) |
| `--alias-config` | auto-detect | Path to virus-specific alias JSON config |
| `--alias-registry` | `app/config/virus_alias_registry.json` | Registry mapping virus keywords → alias configs |
| `--min-coverage` | `0.5` | Minimum protein coverage for a lifted feature to be accepted |
| `--min-identity` | `0.3` | Minimum protein identity for a lifted feature to be accepted |
| `--evalue` | `1e-5` | tblastn E-value threshold |
| `--rescue-window` | `50` | Window (bp) to scan for nearby ATG if start codon is missing |
| `--quiet` | off | Suppress per-record detail; still shows summary |

### CLI Output

Results are written to `<output>/extracted_cds.tsv`:

```
record_id   name    source_name  start   end     strand  status       coverage  identity  method
AF331831.1  GP5     ORF5         13880   14482   +       ok           0.99      95.2      tblastn
AF331831.1  M       ORF6         14467   14991   +       ok           1.00      97.1      tblastn
AF331831.1  N       ORF7         14981   15352   +       ok_rescued   0.98      94.4      tblastn
AF331831.1  ORF1ab  polyprotein  190     12173   +       ok           1.00      91.3      direct
```

- **name** — canonical name from the alias config key
- **source_name** — original raw name in the query annotation (direct extracts only)

---

## 🗂️ Project Structure

```
viralift/
├── app/
│   ├── config/
│   │   ├── prrsv_alias.json           # PRRSV canonical name aliases
│   │   ├── fmdv_alias.json            # FMDV canonical name aliases
│   │   └── virus_alias_registry.json  # Keyword → alias config registry
│   ├── data/
│   │   ├── PRRS_ref_test.gb           # PRRSV reference genome
│   │   ├── FMD_ref_test.gb            # FMDV reference genome
│   │   └── cross_check/               # Validation datasets
│   └── src/
│       ├── main.py                    # 🚪 CLI entry point + core pipeline
│       ├── lifting/
│       │   ├── base.py                # 📦 LiftedFeature dataclass
│       │   └── tblastn_lifter.py      # 🔬 Protein-guided coordinate lifting
│       ├── annotation/
│       │   ├── gene_alias.py          # 📖 Alias lookup and normalisation
│       │   ├── alias_registry.py      # 🔍 Auto-detection of alias config
│       │   └── annotation_strategy.py # 🧭 Feature type detection (CDS / mat_peptide / None)
│       └── io/
│           ├── genbank_parser.py      # 📂 GenBank parsing utilities
│           └── run_logger.py          # 📝 Run and alias audit logging
├── ui/
│   ├── streamlit_app.py               # 🖥️ Web UI (4-stage Streamlit app)
│   └── requirements.txt               # UI-specific Python dependencies
├── logs/                              # 📋 Runtime logs (gitignored, volume-mounted)
├── output/                            # 📤 CLI output files (gitignored, volume-mounted)
├── Dockerfile
├── docker-compose.yml
└── README.md
```

---

## 🧪 Alias Config Format

Alias configs live in `app/config/`. Each file maps a canonical gene name to all known raw name variants:

```json
{
  "virus": "PRRSV",
  "ignored_names": ["polyprotein"],
  "canonical_names": {
    "GP5": [
      "GP5",
      "gp5",
      "ORF5",
      "glycoprotein 5",
      "major envelope glycoprotein"
    ],
    "N": [
      "N",
      "ORF7",
      "nucleocapsid protein",
      "N protein"
    ]
  }
}
```

- **`canonical_names`** — the key is the canonical output name; the list is every alias that maps to it
- **`ignored_names`** — names to exclude from alias scanning entirely (e.g. FMD's `polyprotein`, which is the whole-genome parent feature and carries no gene-level information)

### Adding a new virus

1. Create `app/config/<virus>_alias.json` following the format above
2. Add an entry to `app/config/virus_alias_registry.json`:

```json
{
  "virus_name": "MyVirus",
  "keywords": ["my virus", "myv", "myvirus strain"],
  "alias_config": "config/myvirus_alias.json"
}
```

Keywords are matched case-insensitively against the GenBank record's organism name, description, and accession.

---

## 📋 Logging

ViraLift writes an append-only audit log to `logs/viralift.log` (rotates at 5 MB, keeps 3 backups).

```
2026-04-28 10:23:41 | INFO  | RUN_START        | ref=AY150564.1 | queries=12 | virus=FMDV | alias=fmdv_alias.json | min_cov=0.50 ...
2026-04-28 10:25:03 | INFO  | RUN_COMPLETE     | ref=AY150564.1 | records=12 | features=154 | ok=140 rescued=8 invalid=3 ...
2026-04-28 10:26:11 | INFO  | SESSION_DECISION | "3C protease" → 3Cpro  [saved to config]
2026-04-28 10:26:11 | INFO  | ALIAS_ADDED      | file=fmdv_alias.json | "3C protease" → 3Cpro
2026-04-28 10:26:12 | INFO  | SESSION_DECISION | "VP-1" → VP1  [session only, not saved]
2026-04-28 10:26:13 | INFO  | CANONICAL_ADDED  | file=fmdv_alias.json | new_canonical="VP0"
```

| Event | When it fires |
|-------|--------------|
| `RUN_START` | Pipeline begins (ref ID, query count, all thresholds) |
| `RUN_COMPLETE` | Pipeline finishes (full status distribution) |
| `SESSION_DECISION` | User maps an unknown name in the resolver — logged whether or not 💾 Save was checked |
| `ALIAS_ADDED` | A raw→canonical mapping is permanently written to the alias config |
| `CANONICAL_ADDED` | A new canonical key is added to the alias config |
| `ERROR` | Any exception during record processing |

---

## 📊 Validation Results

Validated on 20 fully-annotated PRRSV and FMDV records from GenBank. Coordinate accuracy measured by IoU ≥ 0.90 against ground-truth annotations:

| Virus | Method | Accuracy |
|-------|--------|----------|
| PRRSV | 🥇 **tblastn** | **93.9%** |
| PRRSV | minimap2 | 43.3% |
| FMDV  | 🥇 **tblastn** | **100%** |
| FMDV  | minimap2 | 40.0% |

Alias name coverage (unique raw names resolved to canonical):

| Virus | Coverage |
|-------|----------|
| PRRSV | 95.2% (80 / 84 unique names) |
| FMDV  | 100% (28 / 28 unique names) |

> 💬 The 4 unmapped PRRSV names (`unknown protein`, `non-structural protein`, `proteinase`, `unglycosylated envelope protein`) are genuinely ambiguous and intentionally left as raw fallback.
