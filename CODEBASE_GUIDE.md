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
│   ├── result_writer.py        write TSV output, compute status summary
│   └── run_logger.py           rotating audit log for alias decisions and runs
├── features/
│   ├── annotation_strategy.py  select feature type; decide direct vs tblastn
│   ├── ref_loader.py           orchestrate reference preparation
│   └── direct_extractor.py     extract + rename annotated features (no alignment)
├── alias/
│   ├── alias_registry.py       auto-detect the right alias config for a virus
│   ├── gene_alias.py           build alias lookup, normalize raw names → canonical
│   └── alias_payload.py        build JSON payloads for LLM-assisted alias mapping (planned)
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
| `get_record_metadata(record)` | Returns `{id, name, description, organism, length}` for display/logging. |
| `_feature_to_scored_dict(feature, ...)` | Internal. Converts one Biopython feature into a dict with all qualifier fields extracted: `name, gene, product, label, standard_name, locus_tag, note, start, end, strand, length, rel_start, rel_end, order` |

All coordinates are **1-based inclusive** throughout the codebase, matching the GenBank convention.

**`_LOOKUP_QUALIFIER_KEYS`** is a module-level constant defining which qualifier fields get extracted and in what priority order:
```python
["gene", "product", "label", "standard_name", "locus_tag", "note"]
```
This same list is imported by `gene_alias.py` for consistent alias lookup across all field types.

---

### Step 2 — Detect feature type + decide routing

**`features/annotation_strategy.py`**

```python
select_feature_type(record, alias_lookup=None) → "CDS" | "mat_peptide" | None
get_strategy(query_record, alias_lookup=None)  → ("direct" | "tblastn", feature_type | None)
```

**`select_feature_type`** is the single entry point for feature type selection. It combines two behaviors:

- **With `alias_lookup`**: scores both `CDS` and `mat_peptide` by how many features can be alias-resolved. Scoring formula: `(resolved * 100) + raw - ignored_or_ambiguous`. Returns the level with the higher score, or `None` if neither level is informative.
- **Without `alias_lookup`**: simple existence check with a polyprotein-shell guard — returns `None` if the only CDS is a whole-genome `"polyprotein"` placeholder (no individual gene coordinates to extract).

**`get_strategy`** calls `select_feature_type` once and returns both the routing decision and the selected feature type as a tuple, so the caller does not need to call `select_feature_type` again:

```python
strategy, feature_type = get_strategy(query_record, alias_lookup)
# strategy   → "direct" | "tblastn"
# feature_type → "CDS" | "mat_peptide" | None
```

- **`"direct"`** — query has usable gene-level annotation; extract coordinates without alignment.
- **`"tblastn"`** — query has no useful annotation; lift coordinates from the reference. `feature_type` is `None`.

---

### Step 3 — Prepare reference

**`features/ref_loader.py` → `prepare_reference_features()`**

Called once per run. Returns a 5-tuple — all downstream code unpacks this and does not re-derive the feature type independently:

```python
ref_features, ref_feature_type, alias_config_path, virus_name, alias_lookup = (
    prepare_reference_features(ref_record, alias_config_arg, alias_registry_arg)
)
```

Internal processing order (single pass, no redundant parsing):

```
1. Resolve alias config path:
   │  user-provided --alias-config
   │  → auto-detect from virus_alias_registry.json
   │  → no config (raw names kept)
   ▼
2. Load alias lookup from config (or {} if no config)
   ▼
3. select_feature_type(ref_record, alias_lookup)
   → "CDS" or "mat_peptide"   (raises ValueError if None — ref must be usable)
   ▼
4. parse_cds_features() / parse_mat_peptides()
   ▼
5. apply_alias_to_features()   → normalize all ref feature names
```

After this step, every `ref_feature["name"]` is a canonical key like `"ORF5"` or `"Lpro"`. Names not in the alias config keep their raw value with `name_source == "raw"`.

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
  "alias_config": "app/config/prrsv_alias.json"
}
```

Returns `None` for unregistered viruses — names just pass through unchanged.

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
  "orf5":                     "ORF5",
  "gp5":                      "ORF5",
  "glycoprotein5":            "ORF5",
  "majorenvelopeglycoprotein":"ORF5",
  ...
}
```

`normalize_text()` is applied to every key before insertion and every lookup — strips whitespace, lowercases, removes spaces/hyphens/underscores. So `"GP-5 Protein"` and `"gp5protein"` hit the same entry.

Two special sentinel values in the lookup:
- **`"__ignored__"`** — feature should be skipped entirely (e.g. `"polyprotein"` wrapper CDS)
- **`"__ambiguous__"`** — name is shared across multiple genes; user must resolve manually

---

### Step 6 — Resolve names

**`alias/gene_alias.py` → `apply_alias_to_feature(feature, alias_lookup)`**

The name resolution logic for a single feature dict. Instead of checking only one field, it iterates over **all qualifier fields** in priority order and collects every field that hits the alias lookup:

```
for each field in [gene, product, label, standard_name, locus_tag, note]:
    if normalize(feature[field]) is in alias_lookup → collect hit

0 hits     → name_source = "raw",                  keep original name
1 hit      → name_source = "alias",                use canonical
multiple hits, same canonical → name_source = "alias",  unanimous, use it
multiple hits, different canonicals → name_source = "alias_conflict_resolved",
                                      use the hit from the highest-priority field
hit → "__ignored__"   → name_source = "ignored"
hit → "__ambiguous__" → name_source = "ambiguous"
```

This multi-field strategy handles common GenBank inconsistencies — for example, a record with `/product="envelope protein"` (generic, not in alias) and `/note="ORF5"` (specific, in alias) will correctly resolve to `ORF5` via the `note` field fallback.

Always writes back to the feature dict:
- `raw_name` — original display name before resolution
- `name` — canonical key (or raw if unresolved)
- `name_source` — `"alias"` | `"alias_conflict_resolved"` | `"raw"` | `"ignored"` | `"ambiguous"`

---

### Step 7A — Direct extraction (fast path)

**`features/direct_extractor.py` → `direct_extract_with_alias()`**

Used when the query record already has usable annotation. No alignment needed.

```
parse features from query record (CDS or mat_peptide)
  → apply_alias_to_features()       normalize names using the shared alias lookup
  → for each feature:
      "ignored"    → skip entirely
      "ambiguous"  → LiftedFeature(status="ambiguous_name")
      "raw"        → LiftedFeature(status="unresolved_name")
      resolved, not in ref → LiftedFeature(status="not_in_reference")
      resolved, in ref     → slice query sequence at [start:end]
                             reverse-complement if strand == "-"
                             → LiftedFeature(method="direct", status="ok", coverage=1.0)
```

`source_name` is set to the original raw name when alias resolution occurred (`name_source` is `"alias"` or `"alias_conflict_resolved"`).

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
- **Coverage**: unique query protein positions covered ÷ total reference protein length

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

### tblastn path

| Status | Meaning |
|---|---|
| `ok` | Valid ATG start + in-frame stop codon |
| `ok_rescued` | Start codon was missing but found nearby within rescue window |
| `invalid_boundaries` | Lifted but could not fix start or stop codon |
| `low_coverage` | tblastn hit found but protein coverage below threshold |
| `no_hit` | tblastn returned no alignment for this gene |
| `translation_fail` | Reference feature could not be translated to protein |

### direct path

| Status | Meaning |
|---|---|
| `ok` | Name resolved and coordinates extracted |
| `unresolved_name` | Name not found in alias lookup |
| `ambiguous_name` | Name is known-ambiguous; user must disambiguate manually |
| `not_in_reference` | Name resolved but gene is absent from the chosen reference |

---

## Alias Config Format

```json
{
  "virus": "PRRSV",
  "ignored_names": ["polyprotein", "nonstructural protein"],
  "ambiguous_names": ["envelope protein", "glycoprotein"],
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
- **`ignored_names`** — names excluded from alias scanning entirely (e.g. `"polyprotein"`, a whole-genome wrapper CDS with no useful gene-level information).
- **`ambiguous_names`** — names shared across multiple genes. Features matching these are flagged as `ambiguous_name` and shown to the user for manual disambiguation.

PRRSV canonical naming convention: structural proteins use **ORF names** (`ORF2a`, `ORF2b`, `ORF3`–`ORF7`), replicase ORFs use `ORF1a` / `ORF1b` / `ORF1ab`. Individual nonstructural proteins processed from the polyprotein (NSP2–NSP12) retain **NSP names** since they have no individual ORF designation.

---

## Web UI — `ui/streamlit_app.py`

4-stage state machine. All state lives in `st.session_state` and persists across Streamlit reruns.

### Stage 1 — Upload

User uploads ref + query `.gb` files and sets optional thresholds (`min_coverage`, `min_identity`, `evalue`, `rescue_window`).

On submit:
1. Load and parse both files
2. Run `prepare_reference_features()` for the ref → returns `(ref_features, ref_feature_type, alias_config_path, virus_name, alias_lookup)`
3. Scan query records for gene names not in the alias lookup (`_scan_unknown_names`)
4. Scan ref features for names with `name_source == "raw"` (`_scan_unknown_ref_names`)
5. Any unknowns → go to **Stage 2**. All known → skip to **Stage 3**.

`_scan_unknown_names` checks all qualifier fields (`_LOOKUP_QUALIFIER_KEYS`), not just `gene`/`product`. A feature is only flagged as unknown when **every** field misses the alias lookup — mirroring `apply_alias_to_feature` logic exactly.

### Stage 2 — Resolve

Shown when there are unrecognised names.

**Ref panel** — lists ref names with `name_source == "raw"`. User can add them as new canonical entries (empty alias list) to the config file. Newly added names immediately appear in the query resolver dropdowns.

**Query resolver** — for each unknown or ambiguous feature group:
- Shows **all candidate qualifier values** (e.g. `` `envelope protein` `ORF4` ``) so the user has full context
- Dropdown to pick a canonical key (or ignore)
- 💾 Save checkbox — if checked, **all candidate values** in the group are written to the alias config, not just the representative shown. This ensures future records using any variant of that name are resolved automatically.

On Continue: all candidates are expanded into the session resolver dict, which is merged into the base alias lookup via `_build_effective_lookup()` for the current run.

### Stage 3 — Running (transient)

Immediately processes all query records and advances to Stage 4. For each record:

```python
strategy, query_feature_type = get_strategy(qrec, effective_lookup)
if strategy == "direct":
    results = direct_extract_with_alias(qrec, query_feature_type, ref_features, effective_lookup)
else:
    results = process_one_query_record(ref_record, qrec, ref_features, ref_feature_type, ...)
```

### Stage 4 — Results

Displays status badges, per-record expandable tables, and export options.

Export:
- **TSV** — with canonical names or raw names
- **FASTA** — gene multiselect, quality filters (min coverage/identity, include/exclude rescued), one file per gene or combined

---

## Planned: LLM-Assisted Alias Building

**`alias/alias_payload.py`**

When a run encounters unresolved gene names, the pipeline can build a structured JSON payload containing the unresolved features (names, lengths, qualifier fields — no sequence data). Two task types:

- **`"map_aliases"`** — virus already in registry; some raw names not yet covered. LLM maps them to existing canonical keys.
- **`"build_alias_map"`** — virus not in registry at all. LLM builds a canonical name list from scratch.

The payload structure and builder functions are complete. The LLM API call and response parser are the remaining pieces to integrate. Currently, the manual resolver in Stage 2 of the UI serves this role interactively.

---

## Audit Log — `io/run_logger.py`

Appends structured events to `logs/viralift.log` (rotates at 5 MB, 3 backups).

| Event | When |
|---|---|
| `RUN_START` | Pipeline begins — ref ID, query count, all thresholds |
| `RUN_COMPLETE` | Pipeline finishes — full status distribution |
| `SESSION_DECISION` | User maps a name in the resolver — logged whether or not 💾 Save is checked |
| `ALIAS_ADDED` | A raw→canonical mapping is permanently written to the alias config |
| `CANONICAL_ADDED` | A new canonical key is added to the alias config |
| `ERROR` | Any exception during record processing |
