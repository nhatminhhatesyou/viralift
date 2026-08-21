# Unit 05 — consistency of the three comparison notebooks

`liftoff_compare.ipynb`, `gatu_50case_validate.ipynb` and `gatu_5case_validate.ipynb` are now
aligned with each other and with the accuracy harness
(`../02_lifting_accuracy/` → `run_lifting_accuracy` → `run_tblastn_against_truth`).

## Changes applied

| # | Change | Applies to | Effect on numbers |
|---|---|---|---|
| 1 | **Production config passed explicitly.** `CFG = PipelineConfig()` and `min_coverage / min_identity / evalue / rescue_window` handed to `lift_all_tblastn`, exactly as `run_tblastn_against_truth` does, instead of relying on the function signature's defaults. | all three | none today (values coincide) — but the notebooks no longer silently diverge if `PipelineConfig` is edited |
| 2 | **exact / coord_only / failed buckets** reported, plus per-gene exact rate, matching the accuracy harness's reporting. `coverage_rows` now also records `exact`. | all three | additive |
| 3 | **`dedupe_truth_by_name`** applied to truth before scoring — see below. | all three | PRRSV only: denominator 797 → **796** |
| 4 | **Per-prediction detail + disagreements + failure-mode tables** for both tools (`liftoff_per_feature.tsv`, `liftoff_disagreements.tsv`, `liftoff_failure_detail.tsv`, `liftoff_failure_modes.tsv`), matching what the GATU notebook produces. | liftoff | additive; also makes the README's promised `liftoff_failure_modes.tsv` real |

### Change 3 in detail — a real (small) metric bug

`compare_predictions_to_truth` resolves truth **by name, longest wins** (its `truth_by_name`
loop): at most one truth feature per name can ever be matched. But the comparison notebooks built
their truth-anchored denominator from the raw truth list. When a record annotates the same gene
twice, the duplicate became a denominator slot **no tool could ever fill** — an automatic miss by
construction, not by performance.

One record in the whole corpus does this: PRRSV **`AF331831.1`**, which labels *both* ORF1a and
ORF1b "RNA polymerase", so ORF1b resolves twice. It cost both tools one gene.

`dedupe_truth_by_name()` was added to `_shared/validation_utils.py` (additive; nothing else
changed) and is now applied in all three notebooks. It affects both tools equally, so no
head-to-head conclusion moves — only the absolute PRRSV figures, which now agree with the harness.

## Result — ViraLift now matches the accuracy harness exactly

| Virus | | truth | exact | coord_only | failed | accuracy |
|---|---|---|---|---|---|---|
| PRRSV | accuracy harness | 796 | 706 | 88 | 2 | 99.75 % |
| PRRSV | **liftoff_compare** | **796** | **706** | **88** | **2** | **99.75 %** |
| PEDV | accuracy harness | 505 | 500 | 4 | 1 | 99.80 % |
| PEDV | **liftoff_compare** | **505** | **500** | **4** | **1** | **99.80 %** |
| FMDV | accuracy harness | 1171 | 1137 | 32 | 2 | 99.83 % |
| FMDV | **liftoff_compare** | **1172** | **1138** | 32 | 2 | **99.83 %** |

PRRSV and PEDV are now identical to the harness. FMDV differs by exactly **one truth feature** —
see below — and the accuracy is identical to two decimals.

## Remaining differences, all intentional

1. **Codon check.** Comparison: off for both tools (`codon_required_names=set()`). Harness:
   per-gene exemption via `ref_codon_checkable_genes`. Methodology §4/§6 — external tools do not
   report codon validity, so holding only ViraLift to it would be unfair. Documented, not a drift.

2. **FMDV truth, 1172 vs 1171.** `should_use_target_truth_filter` returns `feature_type == "CDS"`,
   so for FMDV's `mat_peptide` reference it is **False** and the harness takes
   `parse_truth_features(record, alias_lookup, feature_type)` with defaults — i.e. **no `R`
   restriction, but `filter_nested=True`**, which runs `filter_subfeatures` (a feature fully
   inside another and shorter than 80 % of it is dropped as a sub-feature). The comparison
   notebooks always pass `target_names=R`, and in that branch `parse_truth_features` returns
   early, so **no sub-feature filtering happens at all** (the `filter_nested=False` argument the
   notebooks pass is inert on this path — it only matters when `target_names is None`).

   The residual feature is `VP4` in `AY687334.1`. That record annotates `VP2 = 1639..2592` and
   `VP4 = 1684..1938`; VP4 sits entirely inside VP2 at 27 % of its length, so `filter_subfeatures`
   drops it. The nesting is really a **truth error** — biologically VP4 and VP2 are adjacent
   cleavage products, and that record's VP2 span wrongly swallows VP4. Both tools lift VP4 exactly
   (IoU 1.000); it is ViraLift's *VP2* (predicted 1939..2592, the biologically correct span) that
   the record's truth marks wrong.

   Impact: 1 gene in 1172, accuracy unchanged at 99.83 %. Left as-is deliberately — switching the
   comparison to the harness's path would drop the `R` restriction that methodology §2 requires
   for a fair cross-tool denominator, and would change how Liftoff and GATU are scored too, for
   one record.

   **Consequence to respect:** the FMDV comparison denominator is not the FMDV harness
   denominator. Do not quote the two FMDV numbers as if they came from the same count.

3. **Query passed to the lifter.** The GATU notebooks pass the stripped FASTA; liftoff_compare and
   the harness pass the annotated record. Verified equivalent — 450/450 identical coordinates on
   the 50-case set — which is an empirical proof that the lifter never reads query annotation.

4. **Alias normalisation of the external tool's output.** GATU output goes through
   `apply_alias_to_features`; Liftoff does not need it because it inherits canonical names from
   the reference GFF3 this notebook writes. Same end state: names are comparable to truth before
   scoring, so name inconsistency is never charged to either tool.

5. **Unmapped handling.** Both notebooks count a gene a tool never emitted as incorrect (§5).
   `liftoff_failure_modes.tsv` now splits that out — it is the difference between "lifted the
   wrong span" and "produced nothing", and for Liftoff on FMDV it is the whole story:
   **266 of 1172 truth genes (22.7 %) were never emitted**, versus only 6 lifted to wrong
   coordinates. ViraLift emitted every gene in all three viruses.

## Head-to-head result after alignment

| Virus | Liftoff | ViraLift | Δ |
|---|---|---|---|
| FMDV | 76.79 % | 99.83 % | +23.04 |
| PRRSV | 95.35 % | 99.75 % | +4.40 |
| PEDV | 99.21 % | 99.80 % | +0.59 |

Note the shape, not just the gap: on PEDV the two tools are near-equivalent; the FMDV gap is
almost entirely Liftoff failing to emit mature peptides at all, not misplacing them. Both points
belong in the write-up.

## Bonus — all 5 ViraLift misses across the 300 records

From `outputs/liftoff_failure_detail.tsv`:

| Virus | Record | Gene | ViraLift | Truth | IoU |
|---|---|---|---|---|---|
| FMD | MG372730.1 | 3A | 5313..5630 | 5313..5741 | 0.741 |
| FMD | AY687334.1 | VP2 | 1939..2592 | 1639..2592 | 0.686 |
| PED | KX550281.1 | ORF1b | 12811..20637 | 293..20637 | 0.385 |
| PRRS | AF331831.1 | ORF1b | 7684..12069 | 191..7699 | 0.001 |
| PRRS | AF331831.1 | ORF2b | 12076..12297 | 13786..14388 | 0.000 |

Four of the five are **truth-annotation problems, not lifting errors**:

- `AY687334.1` VP2 — the record's VP2 span swallows VP4 (see above); ViraLift's 1939..2592 is the
  biologically correct post-VP4 span.
- `KX550281.1` ORF1b — truth is `293..20637`, i.e. the **whole ORF1ab** filed under the name
  ORF1b. ViraLift returns the actual ORF1b portion.
- `AF331831.1` ORF1b and ORF2b — the record labels both ORF1a and ORF1b "RNA polymerase" and calls
  GP5 "envelope protein E"; adjudication assigns `blame = ref_truth`.

Only `MG372730.1` 3A is a genuine boundary miss (3' end short by 111 bp).

This is worth a line in the paper, but state it as an *observation with the adjudication rule
behind it*, not as "our failures don't count" — a reviewer will read the second version as
special pleading.
