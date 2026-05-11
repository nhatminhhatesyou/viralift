# ViraLift — Codebase Guide

A walkthrough of every module: what it does, when it's called, and how data moves from input files to final output.

---

## Big Picture

Two entry points, one shared pipeline.

```
UI (Streamlit)          CLI (main.py)
      │                       │
      └──────────┬────────────┘
                 ▼
         shared pipeline
   load files → prepare ref → route per query → lift/extract → validate → output
```

The pipeline always does the same thing regardless of entry point. The UI just wraps it with file upload, an interactive resolver for unknown gene names, and a results viewer.

---

## Module Map

```
app/src/
├── io/
│   ├── genbank_parser.py       load .gb files, parse features into dicts
│   └── result_writer.py        write TSV output and print summary
├── features/
│   ├── annotation_strategy.py  decide which feature type a record uses
│   ├── ref_loader.py           orchestrate reference preparation
│   └── direct_extractor.py     extract + rename annotated features (no alignment)
├── alias/
│   ├── alias_registry.py       auto-detect the right alias config for a virus
│   └── gene_alias.py           build alias lookup, normalize raw names → canonical
└── lifting/
    ├── base.py                 LiftedFeature dataclass (output from any path)
    ├── tblastn_lifter.py       protein-guided coordinate lifting via tblastn
    └── validator.py            check / rescue start and stop codons
ui/
└── streamlit_app.py            4-stage web interface
```

---

## Step-by-Step Pipeline

### Step 1 — Load files

**`io/genbank_parser.py`**

Turns `.gb` files into Biopython `SeqRecord` objects and then into plain Python dicts.

| Function | What it does |
|---|---|
| `load_single_genbank(path)` | Loads exactly one record. Throws if the file has 0 or more than 1. Used for the **reference**. |
| `load_genbank_records(path)` | Loads all records from a file. Used for the **query** (may have dozens of genomes). |
| `parse_cds_features(record)` | Extracts `CDS` features into dicts. |
| `parse_mat_peptides(record)` | Same but for `mat_peptide` features — FMDV stores gene names here instead of in CDS. |
| `_feature_to_scored_dict(feature, ...)` | Internal. Converts one Biopython feature into a dict with all qualifier fields extracted: `name, gene, product, label, standard_name, locus_tag, note, start, end, strand, ...` |

All coordinates are **1-based inclusive** throughout the codebase, matching the GenBank convention.

**`_LOOKUP_QUALIFIER_KEYS`** is a module-level constant defining which qualifier fields get extracted and in what priority order:
```python
["gene", "product", "label", "standard_name", "locus_tag", "note"]
```
This same list is imported by `gene_alias.py` for consistent alias lookup across all field types.

---

### Step 2 — Detect feature type

**`features/annotation_strategy.py`**

```python
get_feature_type(record) → "CDS" | "mat_peptide" | None
get_strategy(query_record, ref_feature_type) → "direct" | "tblastn"
```

`get_feature_type` asks: does this record have `mat_peptide` features? If yes, use those. Otherwise, does it have `CDS`? If neither, return `None`.

`get_strategy` decides how to process a query record:
- **`"direct"`** — record has usable gene-level annotation, extract coordinates without alignment
- **`"tblastn"`** — record has no annotation, or only a single shell polyprotein CDS — must lift coordinates from the reference

Special case: a record with exactly one CDS whose name contains `"polyprotein"` or is blank is treated as unannotated and routed to tblastn, because it has no individual gene coordinates to extract from.

---

### Step 3 — Prepare reference

**`features/ref_loader.py` → `prepare_reference_features()`**

Called once per run. Orchestrates everything needed to get the reference ready:

```
ref_record
  ├─ get_feature_type()          → which feature type the ref uses
  ├─ parse_cds_features()        → raw feature dicts
  ├─ detect_alias_config_for_record()  → find the right alias JSON
  ├─ load_alias_lookup()         → build flat normalized lookup dict
  └─ apply_alias_to_features()   → normalize all ref feature names

returns: (ref_features, alias_config_path, virus_name, alias_lookup)
```

After this step, every `ref_feature["name"]` is a canonical key like `"ORF5"` or `"Lpro"`. Names not found in the alias config keep their raw value and get `name_source == "raw"`.

---

### Step 4 — Auto-detect alias config

**`alias/alias_registry.py`**

```python
detect_alias_config_for_record(record, registry_path) → Path | None
```

Reads `app/config/virus_alias_registry.json`. Each entry maps a set of keywords to an alias config file. The function builds a searchable string from the record's organism, description, and accession fields, then checks if any registered keyword appears in it.

```json
{
  "virus_name": "PRRSV",
  "keywords": ["porcine reproductive and respiratory syndrome", "prrsv"],
  "alias_config": "config/prrsv_alias.json"
}
```

This is what lets you upload any PRRSV GenBank file and automatically load `prrsv_alias.json` without specifying anything. Returns `None` for unregistered viruses — names just pass through unchanged.

---

### Step 5 — Build alias lookup

**`alias/gene_alias.py`**

```python
load_alias_lookup(config_path) → {normalized_alias: canonical_key}
```

Reads the alias JSON config (e.g. `prrsv_alias.json`) and builds a **flat dict** where every alias — including the canonical key itself — maps to the canonical key.

```json
"ORF5": ["GP5", "gp5", "glycoprotein 5", "major envelope glycoprotein"]
```
→
```python
{
  "orf5":                     "ORF5",   # canonical maps to itself
  "gp5":                      "ORF5",
  "glycoprotein5":            "ORF5",
  "majorenvelopeglycoprotein":"ORF5",
  ...
}
```

`normalize_text()` is applied to every key before insertion and every lookup — strips whitespace, lowercases, removes spaces/hyphens/underscores. So `"GP-5 Protein"` and `"gp5protein"` hit the same entry.

---

### Step 6 — Resolve names

**`alias/gene_alias.py` → `apply_alias_to_feature(feature, alias_lookup)`**

The name resolution logic for a single feature dict. Instead of checking only one field, it iterates over **all qualifier fields** in priority order and collects every field that hits the alias lookup:

```
for each field in [gene, product, label, standard_name, locus_tag, note]:
    if normalize(feature[field]) is in alias_lookup → collect hit

0 hits     → name_source = "raw",                   keep original name
1 hit      → name_source = "alias",                 use canonical
multiple hits, same canonical → name_source = "alias",  unanimous, use it
multiple hits, different canonicals → name_source = "alias_conflict_resolved",
                                      use the hit from the highest-priority field
```

This multi-field strategy handles common GenBank inconsistencies — for example, a record with `/product="envelope protein"` (generic, not in alias) and `/note="ORF5"` (specific, in alias) will correctly resolve to `ORF5` via the `note` field fallback.

Always writes back to the feature dict:
- `raw_name` — original display name before resolution
- `name` — canonical key (or raw if unresolved)
- `name_source` — `"alias"` | `"alias_conflict_resolved"` | `"raw"` | `"ignored"`

---

### Step 7A — Direct extraction (fast path)

**`features/direct_extractor.py` → `direct_extract_with_alias()`**

Used when the query record already has usable annotation. No alignment needed.

```
parse features from query record
  → apply_alias_to_features()       normalize names using the shared alias lookup
  → for each resolved feature:
      slice query sequence at [start:end]
      reverse-complement if strand == "-"
      → LiftedFeature(method="direct", status="ok", coverage=1.0)
```

`source_name` is set to the original raw name whenever `name_source` indicates a resolution happened (`"alias"`, `"alias_conflict_resolved"`). Features with `name_source == "ignored"` are skipped entirely.

---

### Step 7B — tblastn lifting (slow path)

**`lifting/tblastn_lifter.py`**

Used when the query record has no usable annotation. This is the core of the project.

#### Translate reference proteins

```python
translate_feature(feature, ref_record)
```

Slices the ref genome at each CDS's coordinates, reverse-complements if on `−` strand, translates to protein. Returns `None` if translation fails or the protein is shorter than 10 aa.

#### Batch tblastn

```python
run_tblastn_batch(proteins, query_genome, tmp_dir, evalue)
```

All reference proteins are written into a **single multi-FASTA query file** and searched against the query genome in **one tblastn subprocess call**. This is significantly faster than calling tblastn once per gene because the genome is indexed only once internally.

Returns `{query_id: [HSPs]}` for each protein.

#### Merge HSPs → genomic coordinates

```python
merge_hsps(hsps) → (start, end, strand, coverage, identity, bit_score)
```

tblastn often returns multiple High-Scoring Pairs for one gene (gaps, divergent regions). Merge strategy:
- **Strand**: majority vote by aligned length
- **Coordinates**: `min(all starts)` to `max(all ends)`
- **Identity**: weighted average by aligned length
- **Coverage**: unique query protein positions covered ÷ total query protein length

Returns 1-based genome coordinates.

#### Validate codons

**`lifting/validator.py`**

```python
validate_cds_boundaries(seq)       → {valid, has_start_codon, has_stop_codon}
rescue_stop_codon(record, ...)     → scan forward codon-by-codon (up to 30 codons)
rescue_start_codon(record, ...)    → scan ±N bp around lifted start for nearest ATG
```

tblastn aligns protein sequence, so the stop codon is not included in the HSP end coordinate — `rescue_stop_codon` always runs to extend the end to the actual stop. If the start codon is also missing, `rescue_start_codon` scans a window (default ±50 bp) upstream of the lifted start.

---

### Step 8 — Output

**`lifting/base.py` → `LiftedFeature`**

Both paths produce `LiftedFeature` dataclass instances. This is the standard data object passed between pipeline layers and written to output.

```python
name            # canonical key, e.g. "ORF5"
source_name     # raw name before alias resolution (direct path only, else None)
ref_start       # 1-based start on reference
ref_end         # 1-based end on reference
strand          # "+" or "-"
query_start     # 1-based start on query genome
query_end       # 1-based end on query genome
sequence        # extracted nucleotide string
coverage        # fraction of ref protein covered (0.0–1.0)
status          # see status codes below
method          # "direct" | "tblastn"
identity        # % identity from BLAST (tblastn only)
score           # bit score (tblastn only)
has_start_codon # bool (tblastn only)
has_stop_codon  # bool (tblastn only)
rescue_offset   # bp offset used to fix start codon (or None)
```

**`io/result_writer.py`**

`write_results_tsv()` flattens all `LiftedFeature` objects to rows in a TSV file. `summarize_counts()` tallies status codes for the console/UI summary.

---

## Status Codes

| Status | Meaning |
|---|---|
| `ok` | Valid ATG start + in-frame stop codon |
| `ok_rescued` | Start codon was missing but found nearby within rescue window |
| `invalid_boundaries` | Lifted but could not fix start or stop codon |
| `low_coverage` | tblastn hit found but protein coverage below threshold |
| `no_hit` | tblastn returned no alignment for this gene |
| `translation_fail` | Reference feature could not be translated to protein |

---

## Alias Config Format

```json
{
  "virus": "PRRSV",
  "ignored_names": ["polyprotein"],
  "canonical_names": {
    "ORF5": [
      "GP5", "gp5", "orf5",
      "glycoprotein 5", "major envelope glycoprotein"
    ],
    "ORF7": [
      "N", "n", "orf7",
      "nucleocapsid protein", "nucleocapsid protein n"
    ]
  }
}
```

- **`canonical_names`** — key is the output canonical name; list contains every known raw name variant that maps to it. The canonical key itself is also automatically included in the lookup.
- **`ignored_names`** — names excluded from alias scanning entirely (e.g. `"polyprotein"` for FMDV, which is a whole-genome wrapper CDS with no useful gene-level information).

PRRSV canonical naming convention: structural proteins use **ORF names** (`ORF2a`, `ORF2b`, `ORF3`–`ORF7`), replicase ORFs use `ORF1a` / `ORF1b` / `ORF1ab`. Individual nonstructural proteins processed from the polyprotein (NSP2–NSP12) retain **NSP names** since they have no individual ORF designation.

---

## Web UI — `ui/streamlit_app.py`

4-stage state machine. All state lives in `st.session_state` and persists across Streamlit reruns.

### Stage 1 — Upload

User uploads ref + query `.gb` files and sets optional thresholds (`min_coverage`, `min_identity`, `evalue`, `rescue_window`).

On submit:
1. Load and parse both files
2. Run `prepare_reference_features()` for the ref
3. Scan query records for gene names not in the alias lookup (`_scan_unknown_names`)
4. Scan ref features for names with `name_source == "raw"` (`_scan_unknown_ref_names`)
5. Any unknowns → go to **Stage 2**. All known → skip to **Stage 3**.

`_scan_unknown_names` checks all qualifier fields (`_LOOKUP_QUALIFIER_KEYS`), not just `gene`/`product`. A feature is only flagged as unknown when **every** field misses the alias lookup — mirroring `apply_alias_to_feature` logic exactly.

### Stage 2 — Resolve

Shown when there are unrecognised names.

**Ref panel** — lists ref names with `name_source == "raw"`. User can add them as new canonical entries (empty alias list) to the config file. Newly added names immediately appear in the query resolver dropdowns.

**Query resolver** — for each unknown feature group:
- Shows **all candidate qualifier values** (e.g. `` `envelope protein` `ORF4` ``) so the user has full context
- Dropdown to pick a canonical key (or ignore)
- 💾 Save checkbox — if checked, **all candidate values** in the group are written to the alias config, not just the representative shown. This ensures future records using any variant of that name are resolved automatically.

On Continue: all candidates are expanded into the session resolver dict, which is merged into the base alias lookup via `_build_effective_lookup()` for the current run.

### Stage 3 — Running (transient)

Immediately processes all query records and advances to Stage 4. For each record:
- `get_strategy()` → `"direct"` or `"tblastn"`
- Direct → `direct_extract_with_alias()`
- tblastn → `process_one_query_record()`

### Stage 4 — Results

Displays status badges, per-record expandable tables, and export options.

Export:
- **TSV** — with canonical names or raw names
- **FASTA** — gene multiselect, quality filters (min coverage/identity, include/exclude rescued), one file per gene or combined

---

## Experimental / Dead Code

| File | Status |
|---|---|
| `alignment/minimap_runner.py` | Experimental — was Phase 1 lifting engine using minimap2. Not called by the current pipeline. minimap2 showed ~40% accuracy vs tblastn's ~94% on PRRSV, so it was replaced. |
| `alignment/sam_lifter.py` | Experimental — parsed SAM output from minimap2. Not used. |

These files are kept for reference but can be safely ignored. The current lifting engine is exclusively tblastn.
