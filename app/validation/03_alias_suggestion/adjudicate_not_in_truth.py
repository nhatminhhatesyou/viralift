"""
Adjudicate the `not_in_truth` rows and fold the verdicts into the gold config.

`not_in_truth` means the curated config has no entry for a name the tool saved --
neither as an alias nor in `excluded_names`. Silence is ambiguous, and the two
readings are opposite:

    "helicase"                          gold holds "helicase; zinc-finger protein"
                                        -> ORF1b, so it plainly treats this family
                                        as ORF1b and simply never listed the bare
                                        form. Gold is incomplete; the tool is right.

    "viral attachment protein"           attachment is carried out by the GP5-M
                                        heterodimer, so the phrase fits ORF5 and
                                        ORF6 equally. Mapping it to ORF5 would
                                        make M records resolve to GP5. Tool wrong.

Because the metric cannot tell these apart, PRRSV precision is only bounded:
73.7% (every silent row counted against) to 95.1% (every silent row ignored).
Deciding each row once, here, replaces that range with one defensible number and
makes the gold config genuinely better for users at the same time.

PROVENANCE WARNING. These names were surfaced BY the tool. Folding them into gold
makes gold partly derived from tool output, so a reconstruction score measured
against the updated config is no longer independent for these 30 names. The
current config is therefore frozen first (`*.gold_frozen_<date>.json`) and that
frozen copy is what the paper's reconstruction number must be measured against.
Report the updated-config score, if at all, as a separate "after curation" figure.

Every verdict below is a judgement about the biology and about what a curator
would list -- not "the tool proposed it, so it is right".

Run from the viralift/ project root:
    python app/validation/03_alias_suggestion/adjudicate_not_in_truth.py --dry-run
    python app/validation/03_alias_suggestion/adjudicate_not_in_truth.py --apply
"""
import argparse
import csv
import json
import shutil
import sys
from datetime import date
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT))

from app.src.alias.alias_manager import save_validated_alias_config  # noqa: E402
from app.src.alias.gene_alias import (  # noqa: E402
    build_alias_lookup,
    lookup_field_value,
    normalize_text,
)

OUT_DIR = Path(__file__).resolve().parent / "outputs"
GOLD = REPO_ROOT / "app/config/prrsv_alias.json"

# --------------------------------------------------------------------------- #
# The 30 PRRSV `not_in_truth` rows, one verdict each.
#
#   ("raw name", canonical, verdict, reason)
#
# verdict:
#   alias    -> add under that canonical
#   exclude  -> add to excluded_names
#   omit     -> not enough evidence for either; NOT written to gold at all
#
# Verdicts marked AMBIGUOUS in the reason were checked against the literature
# (see the Sources note at the bottom of this file), not decided by intuition.
#
# The test is whether the string fits exactly ONE gene in this virus. Nothing
# else -- not how descriptive it reads, not whether it is a noun phrase.
# --------------------------------------------------------------------------- #
#
# REVISION, after checking each verdict against the existing config.
# A first pass used the rule "a predicate (contains X / encodes X / disulfide
# linked to X) is a remark, not a label, so exclude it". That rule is wrong, and
# gold already showed it: gold maps `contains RNA-dependent RNA polymerase` to
# ORF1b. It is right to. If a record's only ORF1b-spanning feature carries that
# note, resolving the note is the only way to extract ORF1b, so a curator working
# from purpose rather than grammar maps it. The operative test is the same one
# used everywhere else in this tool -- does the string fit exactly ONE gene in
# this virus -- and grammatical shape is irrelevant to it.
#
#   contains RNA-dependent RNA polymerase   only ORF1b has RdRp             -> alias
#   contains 1b replicase protein           only ORF1b                      -> alias
#   contains neutralizing epitope           GP5, M and GP4 all carry one    -> exclude
#   viral attachment protein                the GP5-M heterodimer attaches  -> exclude
#   zinc-finger protein                     ZF in nsp1a AND nsp10           -> exclude
#
# A related claim in the report was also wrong: "disulfide linked to the M
# protein" -> ORF5 and "disulfide linked to GP5" -> ORF6 were called a conflict.
# They are different strings, each fitting one gene, and the GP5-M disulfide pair
# is unique in PRRSV. There is no collision.
#
# A row with no adequate basis for either verdict is marked `omit` and is written
# to gold in neither direction -- putting a guess into the ground truth would make
# the score partly a measure of the guess. It stays `not_in_truth`, reported as
# undecided.
#
ADJUDICATION = [
    # -- ORF1a: protease activities of the nsp's processed from this polyprotein.
    # Gold already lists "papain-like cysteine protease" -> ORF1a, so a single
    # named activity at ORF1a coordinates is consistent. Contrast the excluded
    # "cysteine/serine protease": two unrelated activities joined names neither.
    ("cysteine protease", "ORF1a", "alias",
     "single named protease activity of an ORF1a nsp; gold already maps papain-like cysteine protease here"),
    ("serine protease", "ORF1a", "alias",
     "3C-like serine protease is nsp4, inside ORF1a"),
    ("cysteine proteas", "ORF1a", "alias",
     "misspelling of cysteine protease present verbatim in a record; mapping it resolves that record"),
    ("serine protease protease", "ORF1a", "alias",
     "duplicated-word typo of serine protease, present verbatim in a record"),
    ("polypeptide precursor of non-structural proteins NSP1-NSP8", "ORF1a", "alias",
     "NSP1-NSP8 is exactly the ORF1a product; this is a full description of the gene, not of a part"),
    ("contains 1a replicase protein", "ORF1a", "alias",
     "only ORF1a is the 1a replicase; 'contains' does not change which gene it fits"),

    # -- ORF1b: replicase domains. Gold maps RdRp and helicase+zinc-finger here.
    ("helicase", "ORF1b", "alias",
     "gold maps 'helicase; zinc-finger protein' to ORF1b; the bare form was simply never listed"),
    ("zinc-finger protein", "ORF1b", "exclude",
     "AMBIGUOUS: nsp1alpha (ORF1a) has a zinc-finger domain required for IFN-beta "
     "antagonism, and nsp10 (ORF1b) has an N-terminal RING-like + treble-clef "
     "zinc-binding domain. Two different genes, so the bare name picks out neither. "
     "Gold's compound 'helicase; zinc-finger protein' is fine because 'helicase' pins ORF1b"),
    ("nidovirus-like domain", "ORF1b", "alias",
     "nidovirus-specific domain of the ORF1b replicase; unambiguous within PRRSV"),
    ("nidovirus-like domain protein", "ORF1b", "alias",
     "same name with a 'protein' suffix"),
    ("1b replicase protein", "ORF1b", "alias",
     "carries the 1b designation; names the gene directly"),
    ("polypeptide precursor of non-structural proteins NSP9-NSP12, RNA dependent RNA polymerase/helicase",
     "ORF1b", "alias",
     "NSP9-NSP12 is the ORF1b product; a full description of the gene"),
    ("contains 1b replicase protein", "ORF1b", "alias",
     "only ORF1b is the 1b replicase"),
    ("contains RNA-dependent RNA polymerase", "ORF1b", "alias",
     "already an ORF1b alias in gold; only ORF1b carries the RdRp"),
    ("encodes polymerase, helicase and zinc finger", "ORF1b", "alias",
     "that combination of products is ORF1b and nothing else in PRRSV"),
    ("encodes polymerase, helicase, zinc finger and nidovirus motifs", "ORF1b", "alias",
     "same combination, longer wording"),
    ("encodes RNA dependent-RNA polymerase, helicase, zinc finger and conserved nidovirus motifs",
     "ORF1b", "alias",
     "same combination, longer wording"),
    ("ORF1B expressed by a ribosomal frameshifting mechanism", "ORF1b", "alias",
     "names ORF1B explicitly; the mechanism clause is extra text"),

    # -- ORF2b (E): a property fragment vs a real short name.
    ("protein 2b", "ORF2b", "alias",
     "conventional short form of the ORF2b product"),
    ("non-glycosylated", "ORF2b", "exclude",
     "AMBIGUOUS: E (ORF2b), M (ORF6) and N (ORF7) are all unglycosylated -- the PRRSV "
     "reference itself labels ORF6 'unglycosylated membrane protein M'"),

    # -- ORF5 / ORF6. The GP5-M disulfide bond is a well-known unique pair, so each
    # direction of the remark fits exactly one gene. The epitope and attachment
    # wordings are the genuinely ambiguous ones.
    ("disulfide linked to M matrix protein", "ORF5", "alias",
     "the GP5-M disulfide pair is unique in PRRSV; this wording only fits GP5"),
    ("disulfide linked to the M protein", "ORF5", "alias",
     "same bond, same direction"),
    ("disulfide linked to Membrane protein coded for by ORF6", "ORF5", "alias",
     "same bond, names ORF6 as the partner so the subject is GP5"),
    ("disulfide linked to GP5", "ORF6", "alias",
     "same bond from the M side; only M pairs with GP5 this way"),
    ("disulfide linked to the GP5 protein", "ORF6", "alias",
     "same bond, same direction"),

    ("contains neutralizing epitope", "ORF5", "exclude",
     "AMBIGUOUS: neutralizing epitopes are mapped in GP5, in M, and in GP4 for "
     "genotype 1 -- so the remark fits several genes"),
    ("has neutralizing epitope", "ORF5", "exclude",
     "AMBIGUOUS: same as above, different wording"),
    ("viral attachment protein", "ORF5", "exclude",
     "AMBIGUOUS: attachment is carried out by the GP5-M heterodimer, so the phrase "
     "fits ORF5 and ORF6 equally (GP2-GP3-GP4 mediate the later entry step via CD163)"),
    ("major viral attachment protein", "ORF5", "exclude",
     "AMBIGUOUS: attachment is the GP5-M heterodimer; 'major' does not single out GP5 over M"),
    ("presumed cell attachment protein", "ORF5", "exclude",
     "AMBIGUOUS: same heterodimer problem, and the wording is hedged"),

    # Left out entirely -- neither added as an alias nor excluded. The two clauses
    # point different ways: "putative viral attachment protein" is ambiguous between
    # ORF5 and ORF6, while "disulfide linked to membrane protein" does pin GP5. Not
    # enough to justify either verdict, so it stays out of the ground truth.
    ("putative viral attachment protein, disulfide linked to membrane protein", "ORF5", "omit",
     "CONFLICTING SIGNALS: ambiguous attachment claim combined with a GP5-specific "
     "disulfide clause; insufficient basis for either verdict"),
]


def freeze(path: Path) -> Path:
    frozen = path.with_suffix(f".gold_frozen_{date.today():%Y%m%d}.json")
    if not frozen.exists():
        shutil.copy2(path, frozen)
    return frozen


def apply_adjudication(config: dict) -> dict:
    canonical_names = config.setdefault("canonical_names", {})
    excluded = config.setdefault("excluded_names", [])
    added_alias = added_excl = skipped = 0

    for raw, canonical, verdict, _reason in ADJUDICATION:
        if verdict in ("review", "omit"):
            # Not enough evidence to justify either verdict. Writing a guess into
            # the ground truth the tool is scored against would make the score
            # partly a measure of the guess, so these stay out of gold entirely --
            # they remain `not_in_truth` and are reported as undecided.
            continue
        if verdict == "alias":
            entry = canonical_names.setdefault(canonical, [])
            if not isinstance(entry, list):
                entry = entry.setdefault("aliases", [])
            if raw not in entry:
                entry.append(raw)
                added_alias += 1
            else:
                skipped += 1
        else:
            if raw not in excluded:
                excluded.append(raw)
                added_excl += 1
            else:
                skipped += 1

    config["excluded_names"] = sorted(set(excluded), key=normalize_text)
    for name, entry in canonical_names.items():
        if isinstance(entry, list):
            canonical_names[name] = sorted(set(entry), key=normalize_text)
    note = (config.get("notes") or "").strip()
    stamp = (f"Adjudicated {len(ADJUDICATION)} names surfaced by the alias-suggestion "
             f"reconstruction on {date.today():%Y-%m-%d}; see "
             f"app/validation/03_alias_suggestion/adjudicate_not_in_truth.py")
    config["notes"] = f"{note} {stamp}".strip() if stamp not in note else note
    return {"added_alias": added_alias, "added_excluded": added_excl, "already_present": skipped}


def write_table() -> Path:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    path = OUT_DIR / "not_in_truth_adjudication.tsv"
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["raw_value", "canonical", "verdict", "reason"])
        w.writerows(ADJUDICATION)
    return path


def main():
    ap = argparse.ArgumentParser()
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--dry-run", action="store_true", help="report only, touch nothing")
    g.add_argument("--apply", action="store_true", help="freeze gold, then update it")
    args = ap.parse_args()

    from collections import Counter
    tally = Counter(r[2] for r in ADJUDICATION)
    print(f"{len(ADJUDICATION)} names: {tally['alias']} -> alias, "
          f"{tally['exclude']} -> excluded, "
          f"{tally['review'] + tally['omit']} -> LEFT OUT (insufficient evidence)")

    before = build_alias_lookup(json.loads(GOLD.read_text()))
    unseen = [r for r in ADJUDICATION if lookup_field_value(r[0], before) is None]
    print(f"{len(unseen)} of them are currently absent from gold (expected: all)")

    table = write_table()
    print(f"table -> {table}")

    if args.dry_run:
        print("\ndry run, gold untouched. Verdicts:")
        for raw, canonical, verdict, reason in ADJUDICATION:
            print(f"  {verdict:<8} {raw[:56]:<56} {canonical:<7} {reason[:60]}")
        return

    frozen = freeze(GOLD)
    print(f"frozen  -> {frozen.name}   (measure the paper's number against THIS)")

    config = json.loads(GOLD.read_text())
    stats = apply_adjudication(config)
    save_validated_alias_config(GOLD, config)   # tool's own validator + backup
    print(f"updated -> {GOLD.name}: +{stats['added_alias']} aliases, "
          f"+{stats['added_excluded']} excluded, {stats['already_present']} already present")

    after = build_alias_lookup(json.loads(GOLD.read_text()))
    still = [r[0] for r in ADJUDICATION if lookup_field_value(r[0], after) is None]
    print(f"unresolved after update: {len(still)} {still}")
    print(f"lookup size {len(before)} -> {len(after)}")


if __name__ == "__main__":
    main()

# --------------------------------------------------------------------------- #
# Sources consulted for the AMBIGUOUS verdicts
#
#   zinc finger in both ORF1a and ORF1b
#     Xue et al., "The Zinc-Finger Domain Was Essential for PRRSV nsp1-alpha to
#     Inhibit the Production of Interferon-beta" — PMC3665301
#     Deng et al., "Structural basis for the regulatory function of a complex
#     zinc-binding domain in a replicative arterivirus helicase" —
#     Nucleic Acids Research 42(5):3464
#     "A Complex Zinc Finger Controls the Enzymatic Activities of Nidovirus
#     Helicases" — J Virol 79(2):696
#
#   nsp11 NendoU is nidovirus-specific and lies in ORF1b
#     Nedialkova et al., "Biochemical Characterization of Arterivirus
#     Nonstructural Protein 11..." — J Virol 83(11):5671
#     "Research Progress on NSP11 of PRRSV" — PMC10384725
#
#   GP5-M disulfide-linked heterodimer, and attachment vs entry roles
#     ICTV Arteriviridae report — ictv.global/report_9th/RNApos/Nidovirales/Arteriviridae
#     "Palmitoylation of the envelope membrane proteins GP5 and M of PRRSV..." —
#     PMC8099100
#
#   neutralizing epitopes in GP5, M and GP4
#     "Linear epitopes of PRRSV-1 envelope proteins ectodomains are not correlated
#     with broad neutralization" — Porcine Health Management (2024)
# --------------------------------------------------------------------------- #
