# 05 — Tool comparison: LiftOn + Liftoff + GATU

Three external annotation-transfer tools are compared against ViraLift, all scored with the
identical shared metric (coordinate-only, IoU >= 0.90 or +-6 bp, codon check off for every tool,
truth-anchored `R ∩ truth` denominator). `CONSISTENCY_unit05.md` records how the notebooks were
aligned with each other and with `../02_lifting_accuracy/`; `LIFTON_research.md` covers the LiftOn
positioning and its first PRRSV results.

| Tool | Engine | Status |
|---|---|---|
| **LiftOn** (Chao & Salzberg 2025) | Liftoff + miniprot + protein-maximisation | strongest modern baseline — 3 viruses x 100 records, `lifton_compare.ipynb` |
| **Liftoff** (Shumate & Salzberg 2021) | minimap2, DNA only | 3 viruses x 100 records — `liftoff_compare.ipynb` |
| **GATU** (Tcherepanov 2006) | reference-protein BLAST, interactive Java GUI | 50 PRRSV records run by hand — `outputs/gatu_50case/` |

Headline: **ViraLift and LiftOn are tied above ~70 % protein identity** (FMDV 99.83 % vs 99.83 %,
PEDV 99.80 % vs 99.41 %, PRRSV 99.75 % vs 95.35 %); below it LiftOn drops to **0 %** and ViraLift
holds **100 %** (4 PRRSV-1 records at ~65 % identity, where LiftOn emits nothing). LiftOn *does*
handle FMDV `mat_peptide` — there is no capability gap there. Report this as **operating range**,
not average accuracy. Liftoff (DNA only) is the weaker engine baseline: FMDV 76.79 %, of which
22.7 % of truth genes are never emitted. GATU: 88.05 % on 50 PRRSV records.

---


Two external annotation-transfer tools are compared against ViraLift here.

## A. Liftoff — quantitative head-to-head (`liftoff_compare.ipynb`)

Question: **on the same records, same truth and same coordinate metric, how accurately does
Liftoff recover the target gene compared with ViraLift's protein-guided tblastn lifting?**

Liftoff (Shumate & Salzberg, 2021) lifts annotation with **minimap2** nucleotide alignment;
ViraLift lifts with **tblastn** (protein-guided). The notebook runs both on PRRSV, FMDV and PEDV
(100 records each) and scores them with the identical shared harness
(`compare_predictions_to_truth`, IoU ≥ 0.90 or ±6 bp), so any gap is engine/strategy, not metric.

### Run
1. In the project venv: `pip install liftoff`.
2. Ensure `minimap2` is on PATH (the engine-comparison notebook already uses it).
3. Run all cells.

### Main outputs
| File | Use |
|---|---|
| `outputs/liftoff_summary.tsv` | Coordinate-correct % by virus × method |
| `outputs/liftoff_summary_wide.tsv` | Table-friendly (ViraLift − Liftoff) |
| `outputs/liftoff_accuracy.png` | Main figure (bar) |
| `outputs/liftoff_per_feature.tsv` | Per-prediction detail |
| `outputs/liftoff_failure_modes.tsv` | Failure-mode + unmapped breakdown |

### Fairness notes
- The reference is written to GFF3 with **canonical** gene names, so Liftoff inherits the same
  names ViraLift assigns (name matching is not what is tested).
- Liftoff predictions are **not** held to the CDS codon check (benefit of the doubt); only
  coordinate correctness is compared.
- tblastn is scored with the same per-prediction metric as Liftoff, so it is the apples-to-apples
  reference point. The paper's headline per-gene accuracy comes from `run_lifting_accuracy`.

Related: `../04_engine_comparison/` isolates the engine (minimap2 vs tblastn) and predicts this
result — Liftoff is the named, published tool that instantiates the minimap2 strategy.

## B. GATU — external baseline (`gatu_compare.ipynb`, `outputs/gatu_inputs/`)

GATU (Tcherepanov et al., 2006) transfers annotation by BLASTing reference proteins against the
target — conceptually similar to tblastn — but ships as an interactive **Java 8 Web Start** GUI
that processes **one genome at a time** with human confirmation. It cannot be scripted over a
batch of records, which is the surveillance/primer-design setting ViraLift targets. This folder
holds representative FASTA inputs and a manifest; GATU results (where obtainable) go in
`outputs/gatu_output/`.

Report: which cases were selected, why representative, GATU's manual/one-at-a-time setup, and how
ViraLift's batch reference-guided output differs.
