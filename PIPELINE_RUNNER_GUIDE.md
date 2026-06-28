# ViraLift Pipeline Runner Guide

> Navigation: [README](README.md) · **Pipeline Runner** · [Data Crawler](DATA_CRAWLER_GUIDE.md) · [Alias Manager](ALIAS_MANAGER_GUIDE.md)

This document explains ViraLift's main processing stage: from the moment the user feeds a reference + query GenBank into the tool until it returns the TSV/FASTA output.

In the current UI this part is called **Run pipeline**. That label is fine for the button/tab. In technical docs it can be referred to more explicitly as the **Pipeline Runner** or the **ViraLift Pipeline**.

## Table of contents

- [What is the Pipeline Runner?](#what-is-the-pipeline-runner)
- [Input and output](#input-and-output)
- [Overall flow](#overall-flow)
- [Stages in the Web UI](#stages-in-the-web-ui)
- [How does the tool decide direct vs tblastn?](#how-does-the-tool-decide-direct-vs-tblastn)
- [How does the alias config affect the pipeline?](#how-does-the-alias-config-affect-the-pipeline)
- [What the thresholds mean](#what-the-thresholds-mean)
- [Boundary rescue](#boundary-rescue)
- [What the statuses mean](#what-the-statuses-mean)
- [Reading the results](#reading-the-results)
- [Running from the CLI](#running-from-the-cli)
- [Common cases](#common-cases)
- [Best practices](#best-practices)

## What is the Pipeline Runner?

The Pipeline Runner is ViraLift's core processing component. It takes:

```text
1 reference GenBank + many query GenBank records
```

and produces output consisting of:

```text
standardised gene names + coordinates on the query + sequence + mapping status
```

The pipeline handles two kinds of query automatically:

```text
Query already has useful annotation  -> direct extraction
Query lacks useful annotation        -> tblastn lifting from the reference
```

## Input and output

### Input

| Input | Meaning |
|---|---|
| Reference GenBank | A single, well-annotated record used as the source of canonical genes |
| Query GenBank | One or more records that need name standardisation / annotation lifting |
| Alias config | JSON file that standardises gene names for the matching virus |
| Thresholds | Conditions for accepting a tblastn hit: coverage, identity, e-value, rescue window |

### Output

| Output | Meaning |
|---|---|
| TSV canonical | Result table using canonical names (or reference names, depending on the setting) |
| TSV raw | Result table that prefers the original names when the query has annotation |
| FASTA | Extracted gene sequences, filterable by coverage/status |
| Run summary | Counts of features that are OK, rescued, no-hit, or need review |

## Overall flow

```text
Upload reference + query
        |
        v
Detect virus via registry keywords
        |
        v
Load the matching alias config
        |
        v
Parse reference features
        |
        v
Standardise reference names with the alias config
        |
        v
For each query record:
    |
    +-- Has useful annotation?
    |       |
    |       +-- Yes -> direct extraction
    |       |
    |       +-- No  -> tblastn lifting
    |
    v
Validate coordinate / codon / coverage
        |
        v
Show results + export TSV/FASTA
```

## Stages in the Web UI

### 1. Upload

The user uploads:

- Reference GenBank file: a single record.
- Query GenBank file: one or more records.

Advanced options are also available here:

- `min_coverage`
- `min_identity`
- `evalue`
- `rescue_window`
- `Use ref gene names as output`

If the tool recognises the virus, the pipeline continues. If it cannot recognise the virus, the UI moves to the review / new-alias-config stage.

### 2. Virus review

Shown only when the reference matches no keyword in the registry.

The user chooses:

- Use an existing alias config if this is a known virus with unusual metadata.
- Create a new alias config if this is a new virus.

### 3. Alias seed

Shown only when creating a new virus.

The tool takes the reference feature names as the initial canonical names, then can run tblastn to suggest aliases from the query annotation.

### 4. Resolve

Shown only when the query or reference contains names not yet in the alias config.

The user can:

- Map an unknown name to an existing canonical.
- Add a new name as a canonical.
- Save the mapping to the alias config so it is recognised automatically next time.
- Ignore a name if it is too generic or untrustworthy.

### 5. Run

The tool runs for real on each query record:

- record with good annotation -> direct
- record lacking annotation -> tblastn

The UI shows per-record progress.

### 6. Review

Displays:

- Total records processed.
- Total features found.
- Pass rate.
- Number of features needing review.
- A detailed per-record table.
- TSV/FASTA export.

## How does the tool decide direct vs tblastn?

ViraLift uses a strategy function to pick the processing path for each query record.

### Direct extraction

Used when the query already has useful gene-level annotation.

For example, a query with a CDS like:

```text
gene = ORF5
product = major envelope glycoprotein
```

If the alias config can resolve this name to a canonical, the tool needs no alignment. It only:

1. Parses the coordinates already present in the query annotation.
2. Standardises the name with the alias config.
3. Extracts the sequence directly from the query.

Advantages:

- Fast.
- Preserves the query's original annotation.
- Ideal when the query is already well annotated.

### tblastn lifting

Used when the query has no useful annotation.

The tool will:

1. Take the gene from the reference.
2. Translate the reference gene to protein.
3. Run `tblastn` with the reference protein against the query genome.
4. Merge HSPs to infer the gene coordinates on the query.
5. Extract the sequence.
6. Validate the start/stop codon where appropriate.

Advantages:

- Works on queries that are not yet annotated.
- Protein is more conserved than nucleotide, so it suits divergent lineages/serotypes.
- Handles short or variable genes better than minimap2.

## How does the alias config affect the pipeline?

The alias config decides which canonical a raw GenBank name is standardised to.

PED example:

```text
envelope protein       -> E
membrane protein       -> M
nucleocapsid protein   -> N
accessory protein 3a   -> ORF3
ORF1a/1b, Pol1, ORF1ab -> ORF1ab
```

If the alias config is missing a name:

- The query can still be extracted.
- But the name may come out as `unresolved_name`.
- The UI will ask the user which canonical to map it to.

If the alias config is wrong:

- Direct extraction may call the wrong gene.
- The validation / pipeline result gets noisy.
- Fix it in the Alias Manager.

## What the thresholds mean

| Threshold | Meaning | Default |
|---|---|---|
| `min_coverage` | Fraction of the reference protein that tblastn must cover | `0.5` |
| `min_identity` | Minimum protein identity of a hit | `0.3` |
| `evalue` | Statistical-significance threshold for a tblastn hit | `1e-5` |
| `rescue_window` | bp window around the start to search for an ATG when the boundary is off | `200` |

Tips:

- If the viruses are very close: the defaults are usually fine.
- If the query is distant from the reference: you may need to lower `min_identity`.
- If there are many noisy hits: raise `min_coverage` or check the reference.

## Boundary rescue

After `tblastn` finds the gene region on the query, ViraLift still has to check the boundary, because a tblastn HSP is a local alignment and does not always capture the exact CDS start/end.

For features of type `CDS`, the pipeline adds a codon-validation step:

```text
the sequence must:
1. start with ATG
2. end with a stop codon: TAA/TAG/TGA
3. have a length divisible by 3
```

If this check fails, the tool may rescue the boundary.

### Start codon rescue

Start rescue only triggers when the current sequence **does not start with ATG**.

It then searches for ATGs around `q_start` within `rescue_window`.

Key point: the tool does not simply pick the nearest ATG. It scores candidates in this order:

```text
1. prefer in-frame sequences
2. prefer the length closest to the reference CDS length
3. if still tied, prefer the ATG closer to q_start
4. if still tied, prefer upstream over downstream
```

The expected reference length is computed as:

```text
expected_length = protein_length * 3 + 3
```

Because the translated protein excludes the stop codon, we add `+3` for the stop codon.

Example:

```text
Reference ORF5 length: 603 bp
tblastn lifted span: 600 bp, does not start with ATG
rescue_window: 200 bp
```

The tool searches for ATGs around the start. If one ATG makes the sequence close to `603 bp` and still in-frame, that candidate is preferred over a nearer ATG that distorts the length a lot.

If the rescue succeeds and the rescued sequence is valid:

```text
status = ok_rescued
rescue_offset = how many bp the new start differs from the original start
```

### Stop codon rescue

Stop rescue runs before start validation. It is used when the HSP/tblastn span is missing the stop codon at the end.

The tool scans forward in-frame from `q_end`:

```text
q_end + 3
q_end + 6
q_end + 9
...
```

up to 30 codons. If it hits `TAA`, `TAG`, or `TGA`, it may update `q_end`.

Stop rescue no longer takes the first stop codon mechanically. If a reference length is available, it scores stop-codon candidates and prefers the stop that makes the CDS length closest to the reference:

```text
expected_length = protein_length * 3 + 3
```

If start rescue changes `q_start`, the tool retries stop rescue with the new start, because changing the start can change the frame and invalidate the old stop codon.

If there is no reference length, stop rescue falls back to safer behaviour: pick the nearest in-frame stop codon downstream.

### What `start`, `stop`, `frame` mean in the boundary check

In the UI, each tblastn feature can have a `boundary_check` column:

```text
start:yes, stop:no, frame:yes
```

Meaning:

- `start:yes`: the CDS starts with the start codon `ATG`.
- `stop:yes`: the CDS ends with a valid stop codon: `TAA`, `TAG`, or `TGA`.
- `frame:yes`: the total CDS length is divisible by 3, i.e. the sequence can be read as whole codons.

Example:

```text
ATG AAA GGG TAA
```

This is 12 bp, divisible by 3, starts with `ATG`, ends with `TAA`, so:

```text
start:yes, stop:yes, frame:yes
```

If the length is not divisible by 3 — for example missing 1-2 bp at either end — the UI shows:

```text
frame:no
```

If `start:no`, start rescue searches for an `ATG` around the start within `rescue_window`. The current default is `200 bp` because some viral records have boundaries off by more than 50 bp.

### Terminal extrapolation

Terminal extrapolation is a mechanism distinct from codon rescue.

It uses the HSP query coordinates to learn how many amino acids the protein alignment is missing at each end, then extends the boundary by:

```text
missing_aa * 3 bp
```

This mechanism is used for the non-codon-validated path, e.g. `mat_peptide`, because a mat_peptide does not necessarily have its own ATG/stop codon like a CDS.

### Why can ORF1b still have an off start?

For genes like PED/PRRSV `ORF1b`, the start boundary can involve the `ORF1a/ORF1b` region and a frameshift/annotation convention. So:

```text
the end boundary is usually more stable
the start boundary is more prone to drift
```

If the tblastn span already starts with a valid ATG, start rescue will **not trigger**, even when the GenBank truth uses a different start coordinate.

If start rescue does trigger, it still picks the ATG based on in-frame + length close to the reference. For ORF1b, the exact GenBank start sometimes reflects an annotation convention rather than a clear biological start codon. So you may see:

```text
coord_correct = True
exact_match = False
delta_end = 0
delta_start off by a few bp or many bp
```

This pattern should be read as a boundary ambiguity / granularity issue before concluding the tool lifted incorrectly.

## What the statuses mean

| Status | Meaning | Needs review? |
|---|---|---|
| `ok` | Feature found and boundary valid | No |
| `ok_rescued` | Feature found; boundary adjusted and passed start/stop/frame checks | Usually no, but worth noting if frequent |
| `direct` | Feature taken directly from the query annotation | No |
| `invalid_boundaries` | Hit found but boundary/codon invalid | Yes |
| `low_coverage` | Hit coverage below the threshold | Yes |
| `no_hit` | No tblastn hit found | Yes |
| `translation_fail` | Reference feature could not be translated to protein | Yes |
| `unresolved_name` | Query has a name that the alias config could not map | Yes |
| `ambiguous_name` | Name is in the ambiguous list; the user must decide | Yes |
| `not_in_reference` | Query resolved to a canonical name, but the selected reference has no such gene | No; still counts as mapped/pass, but worth noting |

Note: `not_in_reference` does not necessarily mean the tool is wrong. For example, the query has `ORF1ab` but the reference only has `ORF1a` and `ORF1b` — this is an annotation-granularity mismatch.

## Reading the results

The main result has these columns:

| Column | Meaning |
|---|---|
| `query_id` | Query record |
| `name` | Gene name after standardisation |
| `source_name` | Original raw name, if any |
| `ref_start`, `ref_end` | Feature coordinates on the reference |
| `start`, `end` | Feature coordinates on the query |
| `strand` | Gene strand |
| `method` | `direct` or `tblastn` |
| `status` | Mapping status |
| `coverage` | Fraction of the reference protein covered |
| `identity` | Protein identity of the tblastn hit |
| `has_start_codon`, `has_stop_codon` | Codon boundary checks |
| `rescue_offset` | Start-codon offset if start rescue moved the start |

### TSV canonical vs TSV raw

TSV canonical uses the standardised names after alias mapping.

TSV raw prefers the original names from the query when available — useful when you want to audit the original annotation.

### FASTA export

FASTA export lets you:

- Select which genes to export.
- Choose a single combined file or one file per gene.
- Filter by coverage/identity.
- Include or exclude `ok_rescued`.

FASTA header:

```text
>{record_id}|{gene}|{start}|{end}|{strand}
```

## Running from the CLI

### Auto-detect alias config

```bash
venv/bin/python -m app.src.main \
  --reference app/data/PED/PED_ref_1.gb \
  --query app/data/PED/PED_100seqs.gb \
  --output output/ped_run
```

### Specify the alias config manually

```bash
venv/bin/python -m app.src.main \
  --reference app/data/PED/PED_ref_1.gb \
  --query app/data/PED/PED_100seqs.gb \
  --output output/ped_run \
  --alias-config app/config/porcine_epidemic_diarrhea_virus_alias.json
```

### Raise the coverage threshold

```bash
venv/bin/python -m app.src.main \
  --reference app/data/PED/PED_ref_1.gb \
  --query app/data/PED/PED_100seqs.gb \
  --output output/ped_strict \
  --min-coverage 0.9
```

## Common cases

### Query has annotation but the name is non-standard

Example:

```text
gene = GP5
product = major envelope glycoprotein
```

If the alias config has `GP5 -> ORF5`, the pipeline goes direct and outputs `ORF5`.

### Query has no annotation

The pipeline uses tblastn to lift the whole gene from the reference onto the query.

The output name comes from the reference canonical.

### Query has a merged gene but the reference splits it

Example:

```text
query: ORF1ab
ref:   ORF1a + ORF1b
```

Do not force `ORF1ab` into `ORF1a`. Keep a separate canonical `ORF1ab`. If the reference has no `ORF1ab`, the status `not_in_reference` is reasonable.

### Virus not in the registry

The UI moves to Virus Review:

- Pick an existing config if auto-detect is missing a keyword.
- Or create a new alias config via Alias Seed.

### Many `invalid_boundaries`

This can be due to:

- The reference boundary differs from the query truth.
- A hard-to-determine start codon.
- A truncated gene.
- A tblastn local hit missing terminal amino acids.
- An insufficient rescue window or rescue logic picking a suboptimal boundary.

Look at it gene by gene, not just the total count.

## Best practices

- Use a reference with complete, trustworthy annotation.
- For a new virus, create the alias config before running in production.
- Do not map a merged gene onto a split gene just to make the score look nice.
- When the pass rate is low, read the status breakdown first before tweaking thresholds.
- `ok_rescued` is usually acceptable, but if it shows up en masse for the same gene, review the reference/query boundary.
- For FASTA export, filter to `ok` and `direct`, and include `ok_rescued` depending on your goal.
- If the UI auto-detects the wrong virus, fix the keyword in the Alias Manager instead of selecting manually each time.
