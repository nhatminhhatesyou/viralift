# Running the VAPiD comparison (`vapid_compare.ipynb`)

VAPiD is added as a fourth external annotation-transfer tool, scored with the **identical shared
metric** as Liftoff / LiftOn / GATU (coordinate-only, IoU >= 0.90 or +-6 bp, codon check OFF for
both tools, truth-anchored `R ∩ truth` denominator). ViraLift is the real `lift_all_tblastn` engine.

VAPiD hard-requires **MAFFT**, which could not be installed in the automation sandbox (no root;
`mafft.cbrc.jp` and conda/bioconda are network-blocked there). So the VAPiD step is run **locally**,
exactly as GATU was run by hand. Everything else (input generation, ViraLift, scoring) is in the
notebook and runs anywhere.

## Prerequisites

- The repo's Python env (biopython, pandas, matplotlib) — same one the other notebooks use.
- **BLAST+ (`tblastn`)** on PATH — ViraLift's engine (`conda install -c bioconda blast` / `brew install blast`).
- **MAFFT** on PATH — VAPiD's aligner (`conda install -c bioconda mafft` / `brew install mafft`).
- **Internet** — VAPiD does a (discarded) NCBI Entrez lookup for the `--r` accession on every record.

## Steps

1. Clone VAPiD and point the notebook at it:

   ```bash
   git clone https://github.com/rcs333/VAPiD
   export VAPID_DIR=/absolute/path/to/VAPiD      # or edit VAPID_DIR in cell 1
   ```

2. Install MAFFT (see above) and confirm `mafft --version` works.

3. Open `vapid_compare.ipynb` and **Run All**. It will:
   - regenerate inputs in `outputs/vapid_inputs/` (`<virus>_ref.gbk`, `<virus>_targets.fasta`, `<virus>_meta.csv`);
   - run VAPiD once per virus into `outputs/vapid_runs/<virus>/<record>/<record>.tbl`
     (skipped automatically if the `.tbl` files already exist — set `SKIP_VAPID_IF_DONE = False` to force);
   - score both tools and write `outputs/vapid_summary.tsv`, `vapid_summary_wide.tsv`,
     `vapid_coverage_per_gene.tsv`, `vapid_per_feature.tsv`, `vapid_failure_detail.tsv`,
     `vapid_not_emitted.tsv`, and `vapid_accuracy.png`.

### Equivalent CLI (if you prefer to run VAPiD outside the notebook)

VAPiD writes strain folders into the **current directory**, so run it from a per-virus work dir:

```bash
cd outputs/vapid_runs/PRRS
python3 "$VAPID_DIR/vapid3.py" ../../vapid_inputs/PRRS_targets.fasta "$VAPID_DIR/example.sbt" \
    --r PQ623173.1 --f ../../vapid_inputs/PRRS_ref.gbk --metadata_loc ../../vapid_inputs/PRRS_meta.csv
```

Accessions: **PRRS = PQ623173.1, FMD = FJ175661.1, PED = PZ105934.1** (see `outputs/vapid_inputs/manifest.csv`).
Then re-run the notebook's parse/score cells (they will find the `.tbl` files).

## Why user-reference mode (methodology disclosure)

To compare fairly, VAPiD must lift from **ViraLift's reference**, not a database-selected one. VAPiD's
default auto-selects a best-BLAST-hit reference per genome; here we instead force ViraLift's reference:

- `--f <virus>_ref.gbk` — ViraLift's reference sequence + features, so VAPiD maps from the same
  reference and copies the same canonical gene names (mirrors how Liftoff/LiftOn were given the
  reference GFF).
- `--r <accession>` — needed only so VAPiD skips its bundled BLAST-database step; the record VAPiD
  fetches for that accession is immediately deleted and replaced by the `--f` file.

State this in the methods: *"VAPiD was run in user-reference mode (`--f`) with ViraLift's reference,
rather than its default database-selected reference, so the comparison isolates lifting accuracy."*

## Caveats to expect

- **Nucleotide aligner (MAFFT).** VAPiD aligns at the DNA level, so expect Liftoff-like behaviour:
  strong at high identity, degrading across genotypes (watch the PRRSV cross-genotype set).
- **FMD mat_peptides.** FMD truth is `mat_peptide`; VAPiD's main annotation path writes CDS. If VAPiD
  emits few/no mat_peptides, that is a real tool limitation, not a scoring artefact — report it as such.
- **Per-record Entrez lookup** (from `--r`) is slow due to NCBI throttling (VAPiD sleeps 1 s/record).
- Records with IDs > 23 chars trigger a VAPiD "gbf will be corrupted" warning — harmless here; we only
  read the `.tbl`, which is written before that step.

## Verified in-sandbox (before handoff)

- `.tbl` parser (`vapid_tbl_to_rows`) — correct on synthetic multi-interval + reverse-strand features.
- `build_reference_gbk` — GenBank round-trips with the right feature counts (PRRS 9 / FMD 12 / PED 7).
- ViraLift + `coverage_rows` scoring — byte-identical to the working `liftoff_compare.ipynb`
  (not re-run here only because BLAST+ is absent from the sandbox; it runs on your machine).
