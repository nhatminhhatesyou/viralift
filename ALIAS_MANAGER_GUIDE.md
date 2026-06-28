# ViraLift Alias Manager Guide

> Navigation: [README](README.md) · [Pipeline Runner](PIPELINE_RUNNER_GUIDE.md) · [Data Crawler](DATA_CRAWLER_GUIDE.md) · **Alias Manager**

This document explains the Alias Manager module: what it is for, when to use it, and how to read/edit the alias map for a new virus or one that already has a config.

## Table of contents

- [What is an alias map?](#what-is-an-alias-map)
- [Why do you need the Alias Manager?](#why-do-you-need-the-alias-manager)
- [Key files](#key-files)
- [Flow for a known virus](#flow-for-a-known-virus)
- [Flow for a new virus](#flow-for-a-new-virus)
- [Reading the Alias Manager UI](#reading-the-alias-manager-ui)
- [What the name types mean](#what-the-name-types-mean)
- [How the tool suggests aliases for a new virus](#how-the-tool-suggests-aliases-for-a-new-virus)
- [Granularity mismatch](#granularity-mismatch)
- [Worked example](#worked-example)
- [Common warnings](#common-warnings)
- [Best practices](#best-practices)

## What is an alias map?

In GenBank, the same gene can be written under many different names:

| Canonical name | Aliases you might see |
|---|---|
| `ORF5` | `GP5`, `ORF5 protein`, `major envelope glycoprotein` |
| `N` | `nucleoprotein`, `N protein` |
| `S` | `spike protein`, `S protein` |
| `ORF1ab` | `ORF1a/1b`, `ORF1a/b`, `Pol1`, `polyprotein 1ab` |

ViraLift needs to standardise these names to a single canonical name, called the `canonical name`.

Example:

```text
GP5 -> ORF5
ORF5 protein -> ORF5
major envelope glycoprotein -> may be ignored if too descriptive/generic
```

The alias map is per virus, not shared across all viruses. A name like `envelope protein` may be too vague for one virus but map reliably to `E` in another when reference/query evidence is clear.

## Why do you need the Alias Manager?

The Alias Manager lets you:

- See which viruses currently have an alias config.
- Fix canonical names and aliases when the tool maps incorrectly.
- Delete aliases added by mistake.
- Manage `ignored_names` and `ambiguous_names`.
- Edit the keywords used to auto-detect a virus.
- Create an alias config for a new virus from the reference + query annotation.

In short: this is where you manage each virus's "gene-name dictionary".

## Key files

| File | Role |
|---|---|
| `app/config/virus_alias_registry.json` | Registry listing viruses, detection keywords, and alias config paths |
| `app/config/*_alias.json` | Per-virus alias config |
| `app/src/alias/alias_manager.py` | Functions to read/edit/validate the alias config |
| `app/src/alias/alias_bootstrap.py` | Creates an alias config and suggests aliases for a new virus |
| `ui/streamlit_app.py` | The Alias Manager UI and the alias-seed flow |

Each time an alias config is saved via the UI, the tool writes a backup in:

```text
app/config/backups/
```

This backup folder is a runtime artifact and is git-ignored. It exists only to recover from accidental edits to the alias config.

## Flow for a known virus

When the user uploads a reference, the tool reads its metadata, e.g.:

```text
organism = Porcine reproductive and respiratory syndrome virus
```

It then compares this against the `keywords` in `virus_alias_registry.json`.

Example registry entry:

```json
{
  "virus_name": "PRRSV",
  "keywords": [
    "porcine reproductive and respiratory syndrome virus",
    "prrsv",
    "prrs virus"
  ],
  "alias_config": "app/config/prrsv_alias.json"
}
```

If the metadata contains a matching keyword, the tool auto-selects that alias config.

The pipeline then runs normally:

1. Parse features from the reference.
2. Standardise names with the alias config.
3. Parse the query.
4. If the query lacks useful annotation, use tblastn lifting.
5. Output the result with standardised names.

## Flow for a new virus

If the reference matches no virus in the registry, the tool enters the new-virus flow:

1. **Virus review**  
   The tool shows metadata taken from the reference such as `organism`, `description`, `record id`.

2. **The user picks one of two paths**

   - Choose an existing alias config if this is really a known virus whose keyword is simply missing.
   - Create a new alias config if this is a new virus.

3. **Alias seed**

   For a new virus, ViraLift takes the reference feature names as the initial canonical names.

   For example, a PED reference has:

   ```text
   ORF1a, ORF1b, S, ORF3, E, M, N
   ```

   So the new alias config has these canonicals:

   ```json
   {
     "canonical_names": {
       "ORF1a": [],
       "ORF1b": [],
       "S": [],
       "ORF3": [],
       "E": [],
       "M": [],
       "N": []
     }
   }
   ```

4. **Generate suggestions**

   The tool runs tblastn from the reference onto each query record, then compares the coordinates against the real annotation in the query to find which names should be added as aliases.

5. **User approval**

   Tick `save` for any reasonable suggestion. Leave `ignore` for any that are generic or wrong.

6. **Save config**

   The tool saves the new alias config into `app/config/` and adds the virus to the registry.

## Reading the Alias Manager UI

In the sidebar select:

```text
Alias manager
```

The main tabs:

### Registry

Used to edit:

- `virus_name`: the display name shown to the user.
- `keywords`: the keywords used to auto-detect the virus.

Note: auto-detect relies mainly on `keywords`, not just `virus_name`.

### Canonical aliases

Each canonical name has its own panel.

Example:

```text
E · 4 alias(es)
```

Inside is the list of aliases for `E`.

The user can:

- Tick one or more aliases.
- Click `Delete selected` to remove them.
- Tick `Delete canonical` to remove the whole canonical name.
- Add a new alias in the `Add canonical / alias` form.

### Ignored names

Holds names that should not be used to map a gene.

Example:

```text
protein
glycoprotein
unknown protein
replicase polyprotein
```

These names are usually too generic. If put into the alias map, the tool could mis-map many different genes.

Note: `envelope protein`, `membrane protein`, `nucleocapsid protein` do not always have to be ignored. For PED, these names have clear evidence and map respectively to `E`, `M`, `N`.

### Ambiguous names

Holds names that could map to several different genes.

Example:

```text
envelope protein
glycosylated membrane protein
```

If a name could resemble ORF2, ORF5, or ORF6 depending on the virus/record, keep it ambiguous or send it to manual review.

For example, PED has a raw gene `mp` in some records. `mp` could mean the ORF3 accessory membrane protein, but is also easily confused with membrane protein. So keeping `mp` in `ambiguous_names` is safer when it appears on its own.

### Raw JSON

Shows the raw alias config file. Useful for quick debugging.

## What the name types mean

### Canonical name

The final standardised name ViraLift wants to output.

Example:

```text
ORF5
```

### Alias

A different name that is specific enough to map to a canonical.

Example:

```text
GP5 -> ORF5
ORF5 protein -> ORF5
```

### Ignored name

A name to skip because it carries too little information.

Example:

```text
protein
glycoprotein
unknown protein
replicase polyprotein
```

### Ambiguous name

A name that carries information but is not certain enough to map to a single canonical.

Example:

```text
glycosylated membrane protein
```

This name could point to different genes depending on the virus or annotation convention.

## How the tool suggests aliases for a new virus

When you click `Generate suggestions`, the tool does the following:

1. Picks useful feature types in the query, e.g. `CDS` or `mat_peptide`.
2. Takes the fields that may contain a gene name:

   ```text
   gene, product, note, label, standard_name, locus_tag
   ```

3. Runs tblastn reference -> query.
4. Compares the query annotation coordinates against the tblastn-lifted coordinates using IoU.
5. If the IoU is high enough, the tool treats it as evidence that the query feature corresponds to the canonical from the reference.
6. Scores each raw name independently.

For example, the same ORF5 feature might have:

```json
{
  "gene": "GP5",
  "product": "major envelope glycoprotein",
  "note": "ORF5 protein"
}
```

If the query coordinates overlap the tblastn ORF5, the tool might suggest:

| Raw name | Field | Canonical | Action |
|---|---|---|---|
| `GP5` | `gene` | `ORF5` | `save_alias` |
| `ORF5 protein` | `note` | `ORF5` | `save_alias` |
| `major envelope glycoprotein` | `product` | `ORF5` | `ignore` or `manual_review` |

Reason: `GP5` and `ORF5 protein` are specific gene names. `major envelope glycoprotein` describes the protein but is not necessarily a safe alias for every record.

### The scoring formula

The decision `save_alias` / `manual_review` / `ignore` for each raw name is based on **a single accumulated score**. The code is in `app/src/alias/alias_classifier.py` (`classify_alias_candidate`). There are two layers:

**Layer 1 — coordinate evidence.** tblastn lifts the reference feature onto the query; the query's existing annotation is matched to the lifted feature via **IoU** (the overlap of coordinate ranges, 1-based inclusive). Only pairs with `IoU ≥ min_iou` (default `0.90`) are considered — this is the condition for a raw name to count as evidence for a canonical.

**Layer 2 — scoring each raw name.** Each qualifier string is scored independently with additions/subtractions:

| Factor | Points |
|---|---|
| IoU ≥ 0.95 | +5 |
| IoU ≥ 0.90 (and < 0.95) | +4 |
| Same strand as the lifted feature | +1 |
| **Strong** field: `gene`, `label`, `standard_name` | +3 |
| **Weak** field: `product`, `note`, `locus_tag` | −1 |
| Raw name is an **exact match** of the canonical name | +5 |
| Raw name **contains** the canonical name | +4 |
| Short gene symbol matches the numeric part of the canonical (e.g. `GP5` ↔ `ORF5`) | +3 |
| Generic name: `protein`, `glycoprotein`, `polyprotein`, `envelope protein`… | −8 |
| Biological descriptor word (`polymerase`, `capsid`, `nucleocapsid`…) **without** a specific gene name | −4 |
| Biological descriptor word **with** a specific gene name | +1 |
| Looks like a locus tag (e.g. `ABC_001234`) | −6 |

> Note: the three items "exact match / contains / short symbol matches number" are mutually exclusive — only the first satisfied one is added. A name containing a digit (e.g. `ORF3`, `ORF1a`) is **not** treated as generic, even if it contains a word like `orf`.

**Decision thresholds** from the total score:

| Total score | Action | Confidence | Saved by default |
|---|---|---|---|
| ≥ 8 | `save_alias` | high | ✅ yes |
| 3 – 7 | `manual_review` | medium | ❌ no (needs manual approval) |
| < 3 | `ignore` | low | ❌ no |

**Worked scoring example** (ORF5 feature, IoU = 1.0, same strand):

```text
GP5  (field=gene):
  +4  IoU ≥ 0.90
  +1  same strand
  +3  strong field (gene)
  +3  short symbol matches the number "5" of ORF5
  = 11  → ≥ 8 → save_alias (high)

major envelope glycoprotein  (field=product):
  +4  IoU ≥ 0.90
  +1  same strand
  −1  weak field (product)
  −8  generic name ("envelope protein"/"glycoprotein")
  = −4  → < 3 → ignore (low)
```

Design rationale: coordinate evidence (IoU) confirms *which feature* corresponds to which canonical, but **the name string itself** is still scored separately to avoid saving overly generic names (like `glycoprotein`) as global aliases — because a generic name can appear on many different genes in other records.

## Granularity mismatch

Some viruses have different annotation conventions between the reference and the query. The Alias Manager only standardises names; it does not split/merge genes itself.

PED example:

```text
Reference: ORF1a + ORF1b separate
Query:     ORF1ab is a single merged feature
```

In this case you should not map:

```text
ORF1ab -> ORF1a
ORF1a/1b -> ORF1a
Pol1 -> ORF1a
```

Because that would label the query's merged region with the single-gene name `ORF1a`.

The better approach is a separate canonical:

```json
{
  "canonical_names": {
    "ORF1a": ["ORF1A", "ORF1a protein"],
    "ORF1b": ["ORF1B", "ORF1b polyprotein"],
    "ORF1ab": [
      "ORF1",
      "ORF1a/1b",
      "ORF1a/b",
      "ORF 1a/1b",
      "ORF1ab polyprotein",
      "polyprotein 1ab",
      "Pol1",
      "POL1"
    ]
  }
}
```

When the reference has no `ORF1ab`, the output `ORF1ab` may be flagged `not_in_reference`. This is the correct signal: the query and reference differ in annotation granularity — it is not a wrong alias.

## Worked example

Suppose the query has annotation:

```json
{
  "raw_query_names": {
    "gene": "GP5",
    "product": "major envelope glycoprotein",
    "note": "ORF5 protein"
  },
  "query_coords": {
    "start": 13788,
    "end": 14390,
    "strand": "+"
  },
  "best_tblastn_match": {
    "canonical_name": "ORF5",
    "start": 13788,
    "end": 14390,
    "strand": "+",
    "iou": 1.0,
    "coverage": 1.0,
    "identity": 0.94
  }
}
```

Reading the result:

- `IoU = 1.0`: the query annotation and tblastn prediction coordinates overlap completely.
- `coverage = 1.0`: tblastn covers the full reference protein.
- `identity = 0.94`: the protein sequences are very similar.

Reasonable conclusion:

```text
GP5 -> ORF5: should save_alias
ORF5 protein -> ORF5: should save_alias
major envelope glycoprotein: consider ignore/manual_review
```

PED example after reviewing 100 records:

```text
envelope protein       -> E
membrane protein       -> M
nucleocapsid protein   -> N
accessory protein 3a   -> ORF3
ORF1a/1b, Pol1, ORF1ab -> ORF1ab
HNZK1                  -> ignore
mp                     -> ambiguous
```

Here `HNZK1` is a strain/isolate prefix that appears on many different genes, so it should not go into the alias.

## Common warnings

### `X maps to multiple canonicals`

Means the same alias is being mapped to multiple canonicals.

Example:

```text
HNZK1 maps to multiple canonicals: M, N, ORF3, S.
```

This is usually because `HNZK1` is not a gene name but a strain/isolate prefix. Remove it from the alias.

### `X is both ignored and an alias`

Means a name is both in `ignored_names` and in a canonical's aliases.

Example:

```text
ORF3 is both ignored and an alias for ORF3.
```

How to handle:

- If `ORF3` is a real gene name: remove it from ignored.
- If the name is too generic: remove it from the alias.

### No suggestions after clicking Generate suggestions

Possible causes:

- The query has no useful annotation.
- The query annotation has no gene/product/note name field.
- tblastn could not lift the feature.
- The IoU between the query annotation and the tblastn prediction is below threshold.
- The feature type in the query is unsuitable.

Check the diagnostics section in the UI to see which step skipped each record.

## Best practices

- Only put specific-enough names into aliases, e.g. `GP5`, `ORF5 protein`, `N protein`.
- Do not put overly generic names like `protein`, `glycoprotein`, `replicase polyprotein` into aliases if they could appear on multiple genes.
- If the query uses a merged gene like `ORF1ab` but the reference splits `ORF1a`/`ORF1b`, create a separate canonical `ORF1ab` instead of forcing it into `ORF1a`.
- For a new virus, review the suggestions before saving the config.
- If an alias makes the tool map incorrectly, go to the Alias Manager and delete that alias immediately.
- If virus auto-detect picks the wrong config, go to the `Registry` tab and fix the keyword.
- After editing an alias config, re-run a few representative queries to check the output.
- Don't be afraid of over-editing: every save via the UI creates a backup in `app/config/backups/`.
