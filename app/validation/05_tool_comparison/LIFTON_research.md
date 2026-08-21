# LiftOn — what it is, and a full head-to-head on FMDV, PRRSV and PEDV

## Short answer

Yes, LiftOn is newer, and it is the **strongest baseline in this space** — much stronger than
Liftoff or GATU. It is by the same group as Liftoff (Kuan-Hao Chao, Steven Salzberg, JHU),
published in **Genome Research 35(2):311** (accepted 19 Dec 2024; bioRxiv May 2024). Current
release **v1.0.11**; v1.0.9 was June 2026, so it is actively maintained.

## How it works, and why that matters here

LiftOn's whole thesis is the one this project independently arrived at: **DNA alignment alone is
not enough, use the protein.**

1. Runs **Liftoff** (minimap2, DNA) *and* **miniprot** (protein-to-genome) on the same locus.
2. A **chaining algorithm** picks the best protein alignment per gene locus from the two sources.
3. An **ORF search** step fixes truncated transcripts by testing alternative reading frames for the
   longest match to the reference protein.
4. Resolves overlapping loci and finds extra gene copies.

So LiftOn = Liftoff + miniprot + protein-maximisation. ViraLift = tblastn (protein) + codon
validation + start/stop rescue. **Different engines, same core insight** — which is why LiftOn is
the baseline a reviewer will ask about, and why beating Liftoff/GATU alone is not enough.

Benchmarked in the paper on human, mouse, honeybee, rice, *Arabidopsis* (within-species) and
mouse↔rat, *D. melanogaster*↔*D. erecta* (cross-species). **No viral benchmark** — which is the
opening for this paper.

## It runs on viral genomes — verified

Installed v1.0.11 and ran it on all three viruses (`lifton_compare.ipynb`) (same reference `PQ623173.1`,
same stripped target FASTA, same GFF3 the Liftoff notebook writes, same metric). It works: ~2 s
per 15 kb genome, all 9 genes lifted on a typical record.

Install notes (Python 3.11): `pip install --no-deps lifton`, then `intervaltree`, `duckdb`,
`pyarrow`, `gffutils`; `miniprot` must be built from source (`github.com/lh3/miniprot`, v0.18);
`minimap2` from apt. The pinned `numpy==1.21.0` in its metadata will not build, hence `--no-deps`.

## Result — all three viruses, 100 records each, identical metric

`lifton_compare.ipynb`. Same reference, same stripped-FASTA input, same metric as
`liftoff_compare.ipynb` and `gatu_50case_validate.ipynb`.

| Virus | Tool | exact | coord_only | failed | total | accuracy |
|---|---|---|---|---|---|---|
| FMDV | LiftOn | 1134 | 36 | 2 | 1172 | **99.83 %** |
| FMDV | ViraLift | 1138 | 32 | 2 | 1172 | **99.83 %** |
| PEDV | LiftOn | 497 | 5 | 3 | 505 | 99.41 % |
| PEDV | ViraLift | 500 | 4 | 1 | 505 | 99.80 % |
| PRRSV | LiftOn | 672 | 87 | 37 | 796 | 95.35 % |
| PRRSV | ViraLift | 706 | 88 | 2 | 796 | **99.75 %** |

### The hypothesised FMDV capability gap does NOT exist

The earlier draft of this note guessed that LiftOn, being protein-coding-gene centric, might not be
able to express `mat_peptide` features. **It can.** Given the reference mature peptides as GFF3
gene/mRNA/CDS quartets, LiftOn lifted them at **99.83 % — an exact tie with ViraLift**, with only
**1 of 1172** truth genes not emitted. There is no FMDV capability argument to make. Good that this
was measured rather than asserted.

### The real result: a divergence cliff

| mean protein identity to reference | LiftOn | ViraLift | truth genes |
|---|---|---|---|
| **< 70 %** | **0.00 %** | **100.00 %** | 35 |
| 70–80 % | 100.00 % | 100.00 % | 60 |
| 80–90 % | 100.00 % | 100.00 % | 95 |
| 90–95 % | 99.81 % | 99.91 % | 1077 |
| >= 95 % | 99.59 % | 99.67 % | 1206 |

The two tools are **statistically indistinguishable everywhere above ~70 % protein identity**, and
below it LiftOn goes to **zero** while ViraLift stays at **100 %**. Not degraded — zero: LiftOn
emits no output at all for those records.

The cliff is exactly 4 PRRSV records, all **PRRSV-1 (European genotype)** against a PRRSV-2
reference:

| Record | mean identity | LiftOn | ViraLift |
|---|---|---|---|
| `EU076704.1` | 64.9 % | 0/9 | 9/9 |
| `AY366525.1` | 64.9 % | 0/9 | 9/9 |
| `DQ864705.1` | 65.0 % | 0/8 | 8/8 |
| `DQ489311.1` | 65.0 % | 0/9 | 9/9 |

This also explains why the gap shows up on PRRSV only: it is the sole dataset that *contains*
cross-genotype records. Minimum mean identity is **64.9 % for PRRSV**, but **79.5 % for FMDV** and
**91.4 % for PEDV** — neither of the other two reaches the cliff.

Mechanistically this is the tblastn docstring's claim demonstrated: minimap2's nucleotide seeding
fails at that divergence, LiftOn's miniprot rescue does not recover the locus, and protein-guided
tblastn still finds every gene.

### LiftOn independently confirms the ORF1b convention

PRRSV ORF1b `delta_start` distribution, 100 records:

| delta_start | −75 | −12 | −9 | −3 | +6 |
|---|---|---|---|---|---|
| **LiftOn** | 1 | 40 | 0 | 32 | 10 |
| **ViraLift** | 1 | 40 | 2 | 32 | 12 |

Near-identical, and `delta_end = 0` for both on essentially every record. GATU, by contrast, sits
~600 bp downstream (IoU 0.86–0.88, just under the 0.90 cutoff, in 44 of 50 records).

Two independent protein-aware tools converging on the same frameshift start convention is strong
external evidence that GATU's 12 % ORF1b figure is a **convention artefact**, not a detection
failure — far better than asserting it ourselves.

## How to position this in the paper

The honest, and stronger, framing:

> ViraLift matches the state of the art (LiftOn) on within-genotype annotation transfer across three
> viruses, and remains at 100 % where LiftOn produces no output at all — records at ~65 % protein
> identity to the reference, i.e. a different genotype of the same virus.

Do **not** write "more accurate than LiftOn". On 2378 of the 2473 truth genes the two are within
0.1 pp of each other. The claim is **operating range**, not accuracy — and operating range is the
right claim for a surveillance/primer-design pipeline, where a PRRSV-1 sequence entering a
PRRSV-2-referenced workflow is routine.

The other two contributions remain untouched by this comparison, and are what LiftOn does not do at
all: **gene-name standardisation** across inconsistent GenBank submissions, and **routing** between
direct extraction and lifting.

> 📌 Suggested figure: accuracy vs. protein identity, one line per tool
> (`outputs/lifton_accuracy.png`, right panel). It makes the argument without a word of spin —
> flat and overlapping until the cliff, then a vertical separation.

## Why the cliff happens — root cause, and it is not tuning

The obvious objection ("you ran LiftOn at defaults") was tested and **the cliff is architectural,
not a parameter choice.**

Flags tried on all 4 cross-genotype records: `--miniprot-rescue`, `--miniprot-candidate`,
`--miniprot-cross-locus-rescue`, `--no-adaptive-rescue-floor`, and Liftoff's own thresholds down to
`-a 0.05 -s 0.05` (effectively no threshold). **Every combination returns 0 CDS.**

The log says why:

```
[ERROR] Liftoff alignment failed: Liftoff completed all alignment and recovery passes but
        lifted zero features.
[ERROR] LiftOn cannot proceed without a valid Liftoff baseline annotation.
```

**LiftOn hard-requires a Liftoff (minimap2, DNA) baseline and aborts without one.** At ~65 %
protein identity the nucleotide alignment finds nothing, so LiftOn stops — *before* any of its
protein logic can contribute.

And the protein evidence was there the whole time. LiftOn runs miniprot anyway and writes the
result to `lifton_output*miniprot/miniprot.gff3`; that file contains **all 9 genes on all 4
records**. Scored against truth with the same metric:

| | 4 cross-genotype records (35 truth genes) |
|---|---|
| **miniprot alone** (LiftOn's own intermediate file) | **33/35 = 94.3 %** |
| **LiftOn** (which produced that file) | **0/35 = 0 %** |
| **ViraLift** | **35/35 = 100 %** |

miniprot alone would have recovered nearly everything. LiftOn discards it because its DNA-first
pipeline has already given up.

Control: on a within-genotype record (`AB811785.1`, ~92 % identity) Liftoff succeeds, and LiftOn
returns all 9 genes normally. The failure is specific to the DNA-baseline step.

**How to phrase this in the paper.** Not "LiftOn is less accurate" and not "LiftOn's protein
alignment is weaker" — both are false. The accurate statement is about **composition order**:

> LiftOn combines DNA and protein evidence, but requires the DNA-based Liftoff pass to succeed
> first; where nucleotide alignment fails outright — as it does between PRRSV genotypes at ~65 %
> protein identity — LiftOn produces no annotation, even though its own miniprot pass recovers
> 33 of 35 genes. ViraLift's protein-first design has no such dependency.

This is a constructive, verifiable observation (and arguably a useful bug report for the LiftOn
authors), not a knock on the tool's alignment quality.

## Caveats

- The above was reproduced on LiftOn **v1.0.11**. Worth re-checking against a newer release before
  submission, since a fallback-to-miniprot path would remove the effect entirely.
- Only 4 records sit below 70 % identity, all PRRSV. The cliff is sharp and total, but it rests on
  a small n — deliberately sampling more divergent records (PRRSV-1 strains, or a cross-species
  reference) would make it publication-grade.
- Version pinned: LiftOn **v1.0.11**, miniprot **v0.18**, minimap2 **v2.26-r1175**, BLAST+ **2.12.0+**.

Sources: LiftOn paper (Genome Research), the LiftOn GitHub repository, and the author's project page.
