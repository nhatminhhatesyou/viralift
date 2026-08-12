# ViraLift — Validation Methodology (context)

How lifting accuracy is validated, and the exact conditions under which a prediction is
counted correct. This is the shared reference for every accuracy notebook and for the
tool-comparison notebooks.

---

## 0. Core principles

1. **Real tool functions only — no reimplementation.**
   Validation notebooks **import the shipped tool logic verbatim from `app.src`** (extraction,
   alias resolution, tblastn lifting, boundary rescue). The notebook never re-codes what the
   tool does; it only *calls* it. This guarantees the validation measures the real tool, not a
   look-alike. The **metric** functions (IoU, truth comparison, per-gene roll-up) live once in
   `validation/_shared/validation_utils.py` and are shared by all notebooks, so scoring is
   identical everywhere.

   Real tool functions used (imported from `app.src`):
   - `app.src.lifting.tblastn_lifter.lift_all_tblastn` — protein-guided tblastn lift + boundary rescue
   - `app.src.features.ref_loader.prepare_reference_features`, `annotation_strategy.get_strategy / select_feature_type`
   - `app.src.io.genbank_parser` — `load_genbank_records`, `parse_cds_features`, `parse_mat_peptides`
   - `app.src.alias.*` — `detect_alias_config_for_record`, `load_alias_lookup`, `apply_alias_to_features`
   - `app.src.pipeline.run_pipeline` / `PipelineConfig` (production defaults are reused, not re-set)

2. **Leakage-free.** No tool sees the truth. ViraLift lifts from the **reference proteins**;
   external tools (Liftoff/GATU) lift from the **reference annotation** onto the target
   **sequence with its own annotation stripped**. The record's own annotation is read **only**
   by the scorer.

3. **Generic.** No `if virus == X` and no hard-coded gene names anywhere in the lift or the
   metric. Everything (target gene set, codon-checkable set, etc.) is derived from the
   reference and the alias config at run time.

4. **Engine forced.** To measure *lifting* quality (not routing), the tblastn engine is forced
   on every record — including records that already carry annotation.

---

## 1. Truth

- Validate on records that **carry annotation**. That record's own GenBank annotation, after
  **name normalisation through the per-virus alias config**, is the **truth**.
- Coordinates are 1-based inclusive `[start, end]` (`parse_cds_features`: `start =
  location.start + 1`, `end = location.end`) — the same convention as GFF, so no conversion is
  needed when comparing to external tools.

---

## 2. The R ∩ truth restriction (evaluable set)

Accuracy is scored **only on genes whose name is in the reference gene set `R`** (i.e. the
reference actually carries that gene), intersected with the genes present in each record's
truth.

**Why:** a reference-based lifter — ViraLift, Liftoff, GATU, LiftOn — can only ever produce a
gene the reference carries. A truth gene whose name is **outside `R`** is unliftable *by
construction*, not a tool error, so it must not sit in the denominator. Concretely, these are
excluded:

- merge / granularity names: `ORF1ab` when the reference splits `ORF1a` / `ORF1b`; a single
  `3B` vs a split `3B1 / 3B2 / 3B3`;
- non-gene or catch-all labels from inconsistent lab submissions: `polyprotein`,
  `hypothetical protein`, `non-structural protein`, strain-specific ORFs such as `HNZK1`.

Implementation: `parse_truth_features(..., filter_nested=False, target_names=R,
keep_extra_names=[])`, with `R = {feature['name'] for feature in reference_features}`. This is
the same restriction the paper's accuracy harness applies via `parse_validation_truth_features`
→ `truth_target_names` → `should_use_target_truth_filter` (denominator = `R ∩ query-truth`).

**Important separation of concerns.** Name inconsistency (the "alias gap") is **not** folded
into the lifting metric. In lifting we *normalise names first* (via the alias config) and then
score **coordinates** on `R ∩ truth`. The **name-standardisation capability itself** is
validated in a **separate** harness (alias coverage / config reconstruction), never conflated
with lifting accuracy.

---

## 3. When is a gene "coordinate-correct"?

A gene is correct if its lifted coordinates satisfy:

> **IoU(pred, truth) ≥ 0.90**  **OR**  **both boundaries within ±6 bp** of truth

- **IoU** = overlap / union of `[start, end]` vs truth (1.0 = identical span).
- The **±6 bp** absolute tolerance rescues short features where IoU is over-sensitive (a
  1–2 codon annotation-convention difference leaves a long feature at IoU ≈ 0.99 but pushes a
  very short one below 0.90).
- **`exact`** (`start == truth_start and end == truth_end`) is a strict subset, tracked for
  reference.

---

## 4. Codon handling

The lifter computes `has_start_codon` / `has_stop_codon` / `in_frame` for CDS lifts. How these
are used depends on the harness:

- **ViraLift standalone accuracy:** hold a gene to the codon check **only if its reference CDS
  is itself codon-clean** (valid start + stop + frame), computed by
  `ref_codon_checkable_genes(bundle)`. Genes whose reference is partial / frameshift /
  overlapping — e.g. **PRRSV ORF1b** (`<7687..12072`, a −1 ribosomal frameshift with no
  independent ATG start) — are **exempt**: an "invalid codon" there is biology, not a tool
  error. Non-CDS reference types (mat_peptide) carry no codon check at all.
- **Cross-tool comparison:** the codon check is **disabled for both tools**
  (`codon_required_names=set()`), because external tools (Liftoff/GATU) do not report codon
  validity — holding only ViraLift to it would be unfair, and it would wrongly fail ORF1b even
  when coordinates are exact. The comparison is therefore **coordinate-only**.

---

## 5. Denominator

**Truth-anchored, per gene.** For each `(record, gene ∈ R ∩ truth)`, the gene is correct if the
tool produced a same-name prediction that is coordinate-correct. A gene the tool **fails to
lift at all counts as incorrect** — a tool is never rewarded for skipping hard genes.

`accuracy = correct / (reference genes present in truth)`.

---

## 6. Two harnesses (do not conflate their numbers)

| Harness | Purpose | Codon check | Denominator | Entry point |
|---|---|---|---|---|
| **ViraLift standalone accuracy** | headline lifting accuracy of ViraLift | per-gene exemption (`ref_codon_checkable_genes`) | `R ∩ truth`, per gene | `run_lifting_accuracy` → `run_tblastn_against_truth` → `summarize_per_gene` / `summarize_overall` |
| **Tool comparison** | fair head-to-head vs Liftoff / GATU | **off for both tools** | `R ∩ truth`, truth-anchored | `05_tool_comparison/liftoff_compare.ipynb`, `gatu_5case_validate.ipynb` |

Both harnesses use the **same lift** (`lift_all_tblastn`) and the same coordinate rule
(§3); they differ only in codon handling and reporting. The comparison number is therefore
*stricter* than, and must not be quoted as, the standalone headline.

---

## 7. What is NOT measured here

- **Name standardisation / alias reconciliation** → separate harness (alias coverage,
  config reconstruction). Lifting validation assumes names already normalised.
- **Routing** (annotated vs unannotated) → the engine is forced, so routing is out of scope
  here; recovery of unannotated records is covered by its own notebook.

---

## 8. Datasets

| Virus | Reference | Query set | Feature level | Reference genes `R` |
|---|---|---|---|---|
| PRRSV | `PRRS_ref_test.gb` (PQ623173.1) | `PRRS_100seq_anno.gb` (100) | CDS | ORF1a, ORF1b, ORF2a, ORF2b, ORF3, ORF4, ORF5, ORF6, ORF7 |
| FMDV | `FMD_ref_test.gb` (FJ175661.1) | `FMD_100seq_anno.gb` (100) | mat_peptide | Lpro, VP4, VP2, VP3, VP1, 2A, 2B, 2C, 3A, 3B, 3Cpro, 3Dpol |
| PEDV | `PED_ref_1.gb` (PZ105934.1) | `PED_100seqs.gb` (100) | CDS | ORF1a, ORF1b, S, ORF3, E, M, N |

---

## 9. Worked note — why one gene can show three different numbers

The **same** ViraLift lift of PRRSV can read as ~99.8% (standalone accuracy), ~95% (comparison,
`R ∩ truth`, codon-off), or ~73% (an earlier, flawed comparison run). The lift never changed —
only the scoring did:

- **codon check on every gene** (no exemption) → penalised ORF1b (no ATG) and any boundary-
  offset gene even with correct coordinates → collapses PRRSV;
- **denominator = predictions emitted** (not `R ∩ truth`) → counts a tool's not-in-truth extra
  predictions as failures, and lets a tool's *unlifted* genes vanish from its own denominator
  (this inflated Liftoff on FMDV, which silently drops ~1/4 of the mature peptides it cannot
  place);
- **including non-`R` truth names** (`polyprotein`, `ORF1ab`, `HNZK1`, `3B1/2/3`) → 0 % for
  both tools, dead weight that depresses absolutes equally.

The methodology above (codon exemption / codon-off, truth-anchored `R ∩ truth` denominator, real
imported lift) is what makes the numbers correct and comparable.
