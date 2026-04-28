# ViraLift — Codebase Guide

> Complete walkthrough of every essential module: what it does, when it is called,
> and how data flows from file upload all the way to exported results.

---

## The Big Picture

Two entry points, one shared pipeline.

```
UI (Streamlit)          CLI (main.py)
      │                       │
      └──────────┬────────────┘
                 ▼
         shared pipeline
   prepare ref → route per query → lift → validate → output
```

---

## Layer 1 — I/O: `app/src/io/genbank_parser.py`

Everything starts here. Turns `.gb` files into Python objects.

| Function | Role |
|---|---|
| `load_single_genbank(path)` | Loads **exactly one** GenBank record. Strict — throws if 0 or 2+. Used for the **ref**. |
| `load_genbank_records(path)` | Loads **all** records from a file. Used for **query** (may contain dozens of genomes). |
| `parse_cds_features(record)` | Extracts `CDS` features → list of plain dicts: `{name, gene, product, start, end, strand, length, ...}` |
| `parse_mat_peptides(record)` | Same but for `mat_peptide` features — used for FMDV where all gene names live in mat_peptide sub-features instead of CDS. |

Raw coordinates are **1-based inclusive** throughout the codebase (same as GenBank).

---

## Layer 2 — Strategy Decision: `app/src/annotation/annotation_strategy.py`

```python
choose_strategy(record) → ("direct" | "minimap", "mat_peptide" | "CDS" | None)
```

Called **once for the ref** and **once per each query record**. It simply asks:
does this genome already have usable annotation?

| Has `mat_peptide`? | Has `CDS`? | Result |
|---|---|---|
| ✅ | — | `("direct", "mat_peptide")` |
| ❌ | ✅ | `("direct", "CDS")` |
| ❌ | ❌ | `("minimap", None)` |

`mat_peptide` takes priority over `CDS` because viruses like FMDV encode a single
whole-genome polyprotein `CDS` which is useless for naming — the real gene names
live in `mat_peptide` sub-features.

### ⚠️ Why does it say "minimap" when we use tblastn?

`"minimap"` is a **legacy label from Phase 1** of the project when minimap2 was
the lifting engine. The code that reads this return value only checks for `"direct"`:

```python
# main.py
granularity_matches = (
    query_strategy == "direct"          # ← only "direct" is special-cased
    and query_feature_type == ref_feature_type
)

if granularity_matches:
    results = direct_extract_with_alias(...)  # fast path
else:
    results = process_one_query_record(...)   # tblastn path — catches "minimap" too
```

So `"minimap"` effectively means **"needs lifting"**, and the actual engine used
is **always tblastn**. The minimap2 code in `alignment/minimap_runner.py` and
`alignment/sam_lifter.py` is dead code — leftover from Phase 1, never called
by the current pipeline.

**The correct label would be `"tblastn"` or `"lift"`**, but renaming it is a
cosmetic refactor with no functional impact.

---

## Layer 3 — Alias System: `app/src/annotation/gene_alias.py`

The naming brain. Standardises raw annotation names into canonical keys.

### `normalize_text(text) → str`

Applied before every lookup. Strips, lowercases, removes spaces / hyphens / underscores.

```
"GP-5 protein"  →  "gp5protein"
"ORF 5"         →  "orf5"
```

### `build_alias_lookup(config) → {normalized_alias: canonical_key}`

Built once from a JSON config file. Every alias (including the canonical key itself)
maps to the canonical key.

```json
"GP5": ["GP5", "ORF5", "glycoprotein 5", "glycoprotein5"]
```
→ `{"gp5": "GP5", "orf5": "GP5", "glycoprotein5": "GP5", ...}`

### `apply_alias_to_feature(feature, alias_lookup) → feature`

Core resolver for one feature dict. Try in order:

1. Normalize `feature["name"]` → lookup → canonical. `name_source = "alias"`
2. Fallback: normalize `feature["product"]` → lookup. `name_source = "product_alias"`
3. No match: keep raw name unchanged. `name_source = "raw"`

Always writes:
- `feature["raw_name"]` = original name before any change
- `feature["name"]` = canonical key (or raw if no match)
- `feature["name_source"]` = `"alias"` | `"product_alias"` | `"raw"`

`"raw"` is the flag that means *"this name is not in the alias DB"* — the UI reads
this field to detect unknown ref names.

---

## Layer 4 — Virus Auto-Detection: `app/src/annotation/alias_registry.py`

```python
detect_alias_config_for_record(record, registry_path) → Path | None
```

Reads `config/virus_alias_registry.json`. Each entry has a virus name, a list of
keywords, and a path to the alias config. The function builds a searchable string
from the record's organism / description / id fields, then checks if any keyword
appears in that string.

```json
{
  "virus_name": "PRRSV",
  "keywords": ["porcine reproductive and respiratory syndrome", "prrsv"],
  "alias_config": "config/prrsv_alias.json"
}
```

This is what lets you upload any PRRSV GenBank file and automatically get
`prrsv_alias.json` without specifying anything. Returns `None` for unregistered viruses
(names just pass through unchanged).

---

## Layer 5 — Reference Preparation: `main.prepare_reference_features()`

First real pipeline function called after loading files. Orchestrates layers 2–4
for the **ref only**:

```
ref_record
    │
    ├─ choose_strategy()           → pick CDS or mat_peptide
    ├─ parse_cds_features()        → raw feature dicts
    ├─ detect_alias_config()       → find the right alias JSON via registry
    ├─ load_alias_lookup()         → build normalized lookup table
    └─ apply_alias_to_features()   → normalize all ref names

returns: (ref_features, alias_config_path, virus_name, alias_lookup)
```

After this, every `ref_feature["name"]` is a canonical key like `"GP5"` or `"Lpro"`.
If a ref name was not in the alias DB, `name` keeps the raw name and `name_source == "raw"`.

---

## Layer 6 — Routing Per Query Record

For every query record, check whether to use the **fast path** or the **slow path**:

```python
query_strategy, query_feature_type = choose_strategy(query_record)

granularity_matches = (
    query_strategy == "direct"
    and query_feature_type == ref_feature_type  # both CDS, or both mat_peptide
)
```

| Query genome | Route |
|---|---|
| Has annotation, same feature type as ref | **Direct extract** (fast, no alignment) |
| Has annotation but different type (e.g. ref=mat_peptide, query=CDS) | **tblastn** (slow) |
| No annotation at all | **tblastn** (slow) |

The "different type" case matters: if ref uses mat_peptide and query only has a
single polyprotein CDS, you cannot name individual genes by direct lookup — you
need protein-guided lifting to find each gene's coordinates independently.

---

## Fast Path: `main.direct_extract_with_alias()`

Used when query already has its own annotation coordinates. **No alignment at all.**

```
query_record
    │
    ├─ parse_cds/mat_peptides()        raw query features
    ├─ apply_alias_to_features()       normalize names using same alias lookup
    └─ for each query feature:
        ├─ name match against ref      → get ref_start/ref_end for the record
        ├─ slice query sequence        → no lifting, direct coordinate slice
        └─ LiftedFeature(
               method="direct",
               status="ok",
               coverage=1.0,
               source_name=raw_name if alias resolved else None
           )
```

`source_name` is populated only when the alias actually resolved something
(`name_source == "alias"`). This is why direct extracts show a "raw name" column
in the UI — you can always trace back to the original annotation.

---

## Slow Path: `main.process_one_query_record()` → `lifting/tblastn_lifter.py`

Used when query has no annotation (or incompatible annotation type).
This is the core scientific contribution of the project.

### Step 1 — `translate_feature(feature, ref_record)`

Slices the ref genome at CDS coordinates, reverse-complements if on `−` strand,
translates to protein with Biopython. Returns a protein string or `None` if
translation fails or the protein is shorter than 10 aa.

### Step 2 — `run_tblastn(protein, query_genome, tmp_dir)`

- Writes protein to a temp FASTA file
- Writes query genome to another temp FASTA file
- Shells out to BLAST+:
  ```
  tblastn -query prot.fa -subject genome.fa -outfmt 5 -evalue 1e-5 -seg no
  ```
- Parses XML result, returns the list of **HSPs** (High-Scoring Pairs) from the
  best alignment, or `None` if no hit.

Why not nucleotide BLAST? Protein is ~3–4× more conserved than nucleotide sequence
across serotypes and lineages. Each gene is searched independently so overlapping
genes (common in dense viral genomes) don't interfere.

### Step 3 — `merge_hsps(hsps)`

tblastn often returns multiple HSPs for one gene (divergent regions, gaps).
Merge strategy:
- **Strand**: majority vote by aligned length
- **Coordinates**: `merged_start = min(all starts)`, `merged_end = max(all ends)`
- **Identity**: weighted average across HSPs by aligned length
- **Coverage**: unique query protein positions covered ÷ total query protein length

Returns `(start, end, strand, coverage, identity, bit_score)` in 1-based genome coords.

### Step 4 — Coordinate filtering

- `coverage < min_coverage` → `status = "low_coverage"`, no sequence extracted
- Otherwise: slice query genome at merged coords, reverse-complement if `−` strand

### Step 5 — `annotation/validator.py` — Codon validation + rescue

**`validate_cds_boundaries(seq)`** — checks:
- `seq[:3] == "ATG"` → has start codon
- `seq[-3:] in {TAA, TAG, TGA}` → has stop codon

**`rescue_start_codon()`** — if start codon missing:
Expands search ±1, ±2, ... ±N bp from the lifted start position, upstream first
(most common fix is a small frameshift). Returns the nearest ATG and the offset used.

**`rescue_stop_codon()`** — if stop codon missing:
Scans forward from `query_end` codon by codon (up to 30 codons = 90 bp) looking
for the next in-frame TAA/TAG/TGA.

Final status codes:

| Status | Meaning |
|---|---|
| `ok` | Valid ATG start + stop codon |
| `ok_rescued` | Start codon was missing but was found nearby |
| `invalid_boundaries` | Could not fix start or stop |
| `low_coverage` | Protein coverage below threshold |
| `no_hit` | tblastn returned no alignment |
| `translation_fail` | Ref protein could not be translated |

---

## The Output Object: `lifting/base.py` — `LiftedFeature`

Every path (direct or tblastn) produces `LiftedFeature` dataclass instances.
This is the standard currency between all pipeline layers.

```python
name          # canonical key, e.g. "GP5"           ← always set
source_name   # raw original name, e.g. "ORF5"      ← direct extracts only

ref_start     # 1-based start on reference genome
ref_end       # 1-based end on reference genome
strand        # "+" or "-"

query_start   # 1-based start on query genome
query_end     # 1-based end on query genome
sequence      # extracted nucleotide string

coverage      # fraction of ref protein covered (0.0–1.0)
status        # see table above
method        # "direct" | "tblastn"

identity      # % identity from BLAST (tblastn only, else None)
score         # bit score (tblastn only, else None)
has_start_codon   # bool (tblastn only)
has_stop_codon    # bool (tblastn only)
rescue_offset     # int, bp offset used to fix start (or None)
```

`.to_dict()` flattens it to a plain dict for TSV export.

---

## UI Layer: `ui/streamlit_app.py`

4-stage state machine driven by `st.session_state.stage`.
All state persists across Streamlit reruns inside `st.session_state`.

### Stage 1 — `stage_upload()`

User uploads ref + query `.gb` files and configures options.

On "Run ViraLift" click:
```
ref file  → load_single_genbank()
query file → load_genbank_records()
ref       → prepare_reference_features()   → ref_features, alias_lookup, virus_name
ref       → choose_strategy()              → ref_feature_type
ref_features → _scan_unknown_ref_names()   → list of ref names with name_source=="raw"
query records → _scan_unknown_names()      → {raw_name: [record_ids]} for unmapped names
```

Routing after scan:
- Any unknowns (ref or query) → **Stage 2 (resolve)**
- All names known → **Stage 3 (running)**

Options set here:
- `min_coverage`, `min_identity`, `evalue`, `rescue_window` (advanced thresholds)
- `use_ref_names` toggle (show ref's raw gene names in output instead of canonical keys)

### Stage 2 — `stage_resolve()`

Shown when there are unrecognised names in the ref or query.

**Ref-side panel** (top, shown only if ref has unknowns):
- Lists ref feature names that have `name_source == "raw"`
- Per-name checkbox to add as a new canonical entry to the alias JSON
- Save button calls `_add_new_canonicals_to_config()` → writes `"new_gene": []`
  to the alias config file (empty alias list, can be filled in later)
- Newly added names immediately appear in the query resolver dropdowns

**Query-side resolver** (below):
- Per-name selectbox: pick a canonical key or "ignore (keep raw name)"
- Per-name 💾 Save checkbox: if checked, the mapping is also written permanently to
  the alias JSON (appends raw_name to the canonical's alias list) via `_save_to_alias_config()`
- Decisions are stored in `st.session_state.resolver`

On "Continue": resolver decisions are merged into the alias lookup for this run
via `_build_effective_lookup()`.

### Stage 3 — `stage_running()` (transient)

Auto-advances immediately to Stage 4.

```
_build_effective_lookup(base_alias_lookup, resolver)
    → effective_lookup (base + user decisions for this run)

_run_pipeline(ref, queries, ref_features, effective_lookup, thresholds)
    → for each query record:
        choose_strategy()
        if direct and type matches → direct_extract_with_alias()
        else                       → process_one_query_record()  [tblastn]
    → [(query_id, [LiftedFeature]), ...]
```

Stores results to `st.session_state.all_results`.

### Stage 4 — `stage_results()`

```
ref_name_map = _canonical_to_ref_map(ref_features)   # if use_ref_names toggle is ON
    → {"Lpro": "Lab", "3Cpro": "3C", ...}

_results_to_df(all_results, ref_name_map)
    → DataFrame (name column uses ref names if map provided)

Display:
    → summary badges (ok / rescued / invalid / low_cov / no_hit)
    → per-record expandable tables (canonical name, raw name, coords, coverage)

Export tabs:
    → TSV: canonical names | raw names
    → FASTA: gene multiselect, quality filter, per-gene or all-in-one
```

The FASTA export uses `ref_name_map` for both the gene list (multiselect options)
and the FASTA headers, so the naming is consistent with whatever the toggle is set to.

---

## Config Files

```
app/config/virus_alias_registry.json
    which keyword → which alias config file

app/config/prrsv_alias.json
app/config/fmdv_alias.json
    {
      "virus": "PRRSV",
      "ignored_names": ["polyprotein"],      ← excluded from alias scanning
      "canonical_names": {
        "GP5": ["GP5", "ORF5", "glycoprotein 5", ...],
        "Lpro": ["Lab", "leader protease", ...]
      }
    }
```

`ignored_names` lets you exclude features like FMD's `"polyprotein"` (the parent
CDS that wraps the entire genome) from alias scanning and coverage statistics,
without affecting the lifting itself.

---

## Dead Code

| File | Status |
|---|---|
| `alignment/minimap_runner.py` | ☠️ Unused. Was the Phase 1 lifting engine (minimap2). |
| `alignment/sam_lifter.py` | ☠️ Unused. Parsed SAM output from minimap2 into a position map. |
| `annotation/extractor.py` | ⚠️ Was used with the minimap path. Still imported in tests but not called by the live pipeline. |
| `annotation/feature_matcher.py` | ⚠️ Structural scoring (length/position/order). Not used in the main pipeline but available as a utility. |
| `annotation/direct_extractor.py` | ⚠️ Simple sequence slicing from existing coords. The current `direct_extract_with_alias()` in `main.py` reimplements this inline with alias support. |

---

## One-Line Summary Per File

| File | Does what |
|---|---|
| `io/genbank_parser.py` | Loads `.gb` files, parses CDS/mat_peptide → raw dicts |
| `annotation/annotation_strategy.py` | Decides direct vs lift per record |
| `annotation/alias_registry.py` | Keyword-matches ref record → finds the right alias JSON |
| `annotation/gene_alias.py` | Normalizes gene names, builds lookup, applies canonical naming |
| `annotation/validator.py` | Checks ATG/stop codons, rescues missing start or stop |
| `annotation/extractor.py` | ☠️ Legacy — coordinate lifting via position map (minimap path) |
| `annotation/direct_extractor.py` | ⚠️ Legacy — simple sequence slice, superseded by `direct_extract_with_alias` |
| `annotation/feature_matcher.py` | ⚠️ Structural scoring utility (length/position/order), not used in pipeline |
| `lifting/base.py` | `LiftedFeature` dataclass — standard output from any lifting engine |
| `lifting/tblastn_lifter.py` | Translate → BLAST → merge HSPs → validate codons → `LiftedFeature` |
| `main.py` | Orchestrates everything: prepare ref, route per query, write TSV/FASTA output |
| `ui/streamlit_app.py` | 4-stage web UI: upload → resolve → run → results/export |
