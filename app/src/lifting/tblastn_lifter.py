"""
tblastn_lifter.py

Protein-guided coordinate lifting using tblastn.

Pipeline per feature:
    1. Translate ref CDS to protein (Biopython)
    2. Write protein as FASTA query, genome as FASTA subject
    3. Run tblastn → parse HSPs
    4. Merge overlapping HSPs to get full gene coordinates
    5. Extract nucleotide sequence from query genome
    6. Validate start/stop codons (with rescue for missing ATG)

Advantages over nucleotide coordinate transfer:
    - Protein is ~3-4x more conserved than nucleotide
    - Works across serotypes and lineages
    - Each gene is searched independently → no interference between overlapping genes

Limitations:
    - Frameshift genes (ORF1b in PRRSV) still have ambiguous start
    - Requires BLAST+ installed (tblastn in PATH)
    - Requires an external BLAST+ binary and is slower than direct extraction
"""

import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.SeqRecord import SeqRecord

from app.src.lifting.validator import (
    validate_cds_boundaries,
    rescue_start_codon,
    rescue_stop_codon,
    first_inframe_stop_end,
)
from app.src.lifting.base import LiftedFeature
from app.src.io.run_logger import log_error


# tblastn binary name; kept as a constant so the pre-flight check and the
# subprocess calls can never drift apart.
TBLASTN_BIN = "tblastn"


class BlastNotInstalledError(RuntimeError):
    """Raised when the tblastn executable cannot be found on PATH."""


def ensure_tblastn_available() -> str:
    """
    Verify that the tblastn executable is on PATH before any lifting starts.

    Returns the resolved path to tblastn. Raises BlastNotInstalledError with a
    clear, actionable message if it is missing — this turns an obscure
    FileNotFoundError deep inside subprocess into a single understandable error
    at the start of the run.
    """
    path = shutil.which(TBLASTN_BIN)
    if path is None:
        raise BlastNotInstalledError(
            "BLAST+ is required for tblastn lifting but 'tblastn' was not found "
            "on PATH. Install NCBI BLAST+ (e.g. `apt-get install ncbi-blast+`, "
            "`conda install -c bioconda blast`, or `brew install blast`) and "
            "ensure 'tblastn' is on PATH."
        )
    return path


# ---------------------------------------------------------------------------
# Translation
# ---------------------------------------------------------------------------

def ref_starts_with_atg(feature: Dict, ref_record: SeqRecord) -> bool:
    """Does the reference feature's own CDS begin with an ATG start codon?

    False for features the reference annotates without a canonical start (partial /
    -1 PRF frameshift genes such as PRRSV ORF1b). Such genes have no ATG to rescue, so
    the lifter keeps the homology-derived start instead of hunting one. Generic and
    data-driven: no gene or virus name is referenced.
    """
    start, end = feature["start"] - 1, feature["end"]
    nuc = ref_record.seq[start:end]
    if feature.get("strand", "+") == "-":
        nuc = nuc.reverse_complement()
    return str(nuc[:3]).upper() == "ATG"


def translate_feature(feature: Dict, ref_record: SeqRecord) -> Optional[str]:
    """
    Translate a ref CDS feature to protein sequence.

    Uses the feature's start/end/strand to extract from ref genome,
    then translates. Stops at first stop codon (to table=1).

    Returns protein string, or None if translation fails.
    """
    start = feature["start"] - 1  # convert to 0-based
    end = feature["end"]
    strand = feature.get("strand", "+")

    nuc = ref_record.seq[start:end]
    if strand == "-":
        nuc = nuc.reverse_complement()

    try:
        protein = str(nuc.translate(to_stop=True, table=1))
    except Exception as exc:
        # Biopython raises TranslationError (a ValueError subclass) plus a few
        # others for bad sequence; keep the broad catch but record the reason so
        # a translation_fail is never a silent black box.
        log_error(f"translate_feature:{feature.get('name')}", exc)
        return None

    if len(protein) < 10:
        return None

    return protein


# ---------------------------------------------------------------------------
# HSP merging → genome coordinates
# ---------------------------------------------------------------------------

def _hsp_strand(hsp) -> str:
    """Subject (genome) strand of a tblastn HSP.

    Biopython's NCBIXML reports sbjct_start <= sbjct_end even for minus-strand
    subject hits — the coordinates are always low->high, and the real strand lives
    in hsp.frame == (query_frame, subject_frame), where a negative subject frame
    means the hit is on the genome's minus strand. Inferring strand from the order
    of sbjct_start/sbjct_end therefore mislabels EVERY minus-strand hit as "+",
    which silently breaks reverse-complement-deposited genomes (all genes come out
    "+", extracted antisense → internal stops → invalid boundaries). Fall back to
    coordinate order only if the frame attribute is unavailable.
    """
    frame = getattr(hsp, "frame", None)
    if isinstance(frame, (tuple, list)) and len(frame) >= 2 and frame[1] is not None:
        return "-" if frame[1] < 0 else "+"
    return "+" if hsp.sbjct_start <= hsp.sbjct_end else "-"


def _hsp_to_genome_coords(hsp) -> Tuple[int, int, str]:
    """
    Convert a tblastn HSP to 1-based genome coordinates and strand.

    Coordinates are the low->high genomic span. Strand comes from the HSP frame
    (see _hsp_strand), NOT from the order of sbjct_start/sbjct_end — Biopython always
    reports those low->high, so coordinate order cannot distinguish strand.
    """
    start = min(hsp.sbjct_start, hsp.sbjct_end)
    end = max(hsp.sbjct_start, hsp.sbjct_end)
    return start, end, _hsp_strand(hsp)


# An HSP implies a gene origin: subtract its offset along the reference protein
# from its genomic position. Real HSPs of one gene agree on that origin; a
# chance hit elsewhere in the genome does not. Tolerance is proportional to the
# gene, because a long gene can carry proportionally larger indels. A -1
# ribosomal frameshift inside a single reference protein (ORF1ab) shifts the
# origin by exactly 1 bp, far inside this tolerance.
_ORIGIN_TOLERANCE_FRACTION = 0.10
_MIN_ORIGIN_TOLERANCE_BP = 150

# Backstop. Viruses have no introns, so a gene occupies roughly ref_len bases.
_MAX_SPAN_FACTOR = 2.0


def _count_n_term_missing_aa(hsps: List, strand: str) -> Optional[int]:
    """
    How many N-terminal reference residues no HSP managed to align.

    Diagnostic only. Returns None rather than raising when the HSP objects do
    not carry the expected coordinate attributes.
    """
    try:
        relevant = [
            h for h in hsps
            if _hsp_strand(h) == strand
        ] or list(hsps)
        return max(0, min(min(h.query_start, h.query_end) for h in relevant) - 1)
    except (AttributeError, TypeError, ValueError):
        return None


def _implied_origin(hsp, strand: str) -> int:
    """Genomic position this HSP implies for the start of the reference protein."""
    start, end, _ = _hsp_to_genome_coords(hsp)
    offset = (min(hsp.query_start, hsp.query_end) - 1) * 3
    return start - offset if strand == "+" else end + offset


def _select_collinear_hsps(hsps: List, protein_length: int, strand: str) -> List:
    """
    Drop HSPs that are not collinear with the main alignment.

    tblastn is a local aligner, so one gene normally yields several HSPs. It
    also reports short chance matches elsewhere in the genome. Merging by
    min(start)/max(end) swallows those, producing a span many times the gene
    length; the protein-axis coverage check cannot see this because the real
    HSPs already cover the whole protein.

    Grouping by implied origin separates the two: genuine HSPs agree, outliers
    do not. Uses only HSP coordinates, bit scores and the reference protein
    length -- never the query's own annotation.
    """
    if len(hsps) <= 1:
        return hsps

    tolerance = max(
        _MIN_ORIGIN_TOLERANCE_BP,
        int(_ORIGIN_TOLERANCE_FRACTION * max(1, protein_length * 3)),
    )
    ordered = sorted(((_implied_origin(h, strand), h) for h in hsps), key=lambda p: p[0])

    clusters = [[ordered[0]]]
    for origin, hsp in ordered[1:]:
        if origin - clusters[-1][-1][0] <= tolerance:
            clusters[-1].append((origin, hsp))
        else:
            clusters.append([(origin, hsp)])

    if len(clusters) == 1:
        return hsps

    def cluster_rank(cluster):
        # Rank by reference protein actually explained, then by total bit
        # score. Counting HSPs would be wrong: one real gene may produce two
        # HSPs while noise produces three tiny ones.
        covered = set()
        for _, hsp in cluster:
            lo, hi = sorted((hsp.query_start, hsp.query_end))
            covered.update(range(lo, hi + 1))
        return len(covered), sum(hsp.bits for _, hsp in cluster)

    return [hsp for _, hsp in max(clusters, key=cluster_rank)]


def _enforce_span_limit(hsps: List, protein_length: int) -> List:
    """
    Last-resort guard: shed edge HSPs while the span stays implausibly long.

    Collinearity filtering handles the usual case; this catches leftovers whose
    implied origins happen to fall within tolerance yet still stretch the span.
    """
    limit = max(1, int(protein_length * 3 * _MAX_SPAN_FACTOR))
    kept = list(hsps)
    while len(kept) > 1:
        coords = [_hsp_to_genome_coords(h) for h in kept]
        low, high = min(c[0] for c in coords), max(c[1] for c in coords)
        if high - low + 1 <= limit:
            break
        edge = [h for h, c in zip(kept, coords) if c[0] == low or c[1] == high]
        kept.remove(min(edge, key=lambda h: h.bits))
    return kept


def merge_hsps(hsps: List, protein_length: int) -> Tuple[int, int, str, float, float, float]:
    """
    Merge a list of HSPs into a single genomic span.

    Strategy:
        - Determine strand from majority vote
        - Discard HSPs not collinear with the main alignment, then any
          remaining outlier that keeps the span implausibly long
        - Take min(start) and max(end) across the surviving HSPs
        - Compute weighted average identity (by aligned length)

    Args:
        hsps: tblastn HSPs for one protein query.
        protein_length: Full length of the reference protein query in amino acids.

    Returns: (merged_start, merged_end, strand, coverage_fraction, identity, bit_score)
    where coverage_fraction is relative to the full query protein length.
    """
    if not hsps:
        raise ValueError("No HSPs to merge")

    # Strand vote (from HSP frame, not coordinate order — see _hsp_strand)
    plus_len = sum(abs(h.sbjct_end - h.sbjct_start) + 1 for h in hsps if _hsp_strand(h) == "+")
    minus_len = sum(abs(h.sbjct_end - h.sbjct_start) + 1 for h in hsps if _hsp_strand(h) == "-")
    strand = "+" if plus_len >= minus_len else "-"

    relevant = [h for h in hsps if _hsp_strand(h) == strand]
    if not relevant:
        relevant = hsps

    relevant = _select_collinear_hsps(relevant, protein_length, strand) or relevant
    relevant = _enforce_span_limit(relevant, protein_length) or relevant

    coords = [_hsp_to_genome_coords(h) for h in relevant]
    merged_start = min(c[0] for c in coords)
    merged_end = max(c[1] for c in coords)

    # Weighted identity
    total_aligned = sum(h.align_length for h in relevant)
    identity = (
        sum(h.identities / h.align_length * h.align_length for h in relevant) / total_aligned
        if total_aligned > 0 else 0.0
    )

    # Coverage: unique query aa covered / full reference protein length.
    # Do not use max(h.query_end) as the denominator; a hit covering only the
    # N-terminus would otherwise look like 100% coverage.
    query_positions = set()
    for h in relevant:
        qs = min(h.query_start, h.query_end)
        qe = max(h.query_start, h.query_end)
        query_positions.update(range(qs, qe + 1))
    coverage = (
        min(1.0, len(query_positions) / protein_length)
        if protein_length > 0 else 0.0
    )

    bit_score = max(h.bits for h in relevant)

    return merged_start, merged_end, strand, coverage, identity, bit_score


def extrapolate_terminal_boundaries(
    start: int,
    end: int,
    strand: str,
    hsps: List,
    protein_length: int,
    genome_length: int,
    coverage: float,
    min_coverage: float = 0.90,
    max_missing_aa: int = 10,
) -> Tuple[int, int, int, int]:
    """
    Conservatively extend tblastn boundaries for protein terminals not covered
    by the HSP span.

    tblastn is a local aligner: terminal amino acids can be omitted from an HSP
    even when the full feature is present in the query genome. HSP query
    coordinates tell us how many reference-protein amino acids are missing at
    the N/C termini. For non-CDS features such as mat_peptide, codon rescue is
    not applicable, so we can extend by missing_aa * 3 bp when coverage is high.

    Returns:
        (new_start, new_end, n_terminal_extended_bp, c_terminal_extended_bp)
    """
    if coverage < min_coverage or protein_length <= 0 or not hsps:
        return start, end, 0, 0

    relevant = [h for h in hsps if _hsp_strand(h) == strand]
    if not relevant:
        relevant = hsps

    min_query_start = min(min(h.query_start, h.query_end) for h in relevant)
    max_query_end = max(max(h.query_start, h.query_end) for h in relevant)

    missing_n_aa = max(0, min_query_start - 1)
    missing_c_aa = max(0, protein_length - max_query_end)

    if missing_n_aa > max_missing_aa or missing_c_aa > max_missing_aa:
        return start, end, 0, 0

    n_extension = missing_n_aa * 3
    c_extension = missing_c_aa * 3

    if strand == "+":
        new_start = max(1, start - n_extension)
        new_end = min(genome_length, end + c_extension)
    else:
        new_start = max(1, start - c_extension)
        new_end = min(genome_length, end + n_extension)

    return new_start, new_end, n_extension, c_extension


# ---------------------------------------------------------------------------
# Sequence extraction + validation
# ---------------------------------------------------------------------------

def extract_sequence(query_record: SeqRecord, start: int, end: int, strand: str) -> str:
    """Extract 1-based inclusive sequence from query genome."""
    seq = query_record.seq[start - 1: end]
    if strand == "-":
        seq = seq.reverse_complement()
    return str(seq)


def _threshold_failure_status(
    coverage: float,
    identity: float,
    min_coverage: float,
    min_identity: float,
) -> Optional[str]:
    """Return the precise quality-gate failure status, if thresholds fail."""
    low_coverage = coverage < min_coverage
    low_identity = identity < min_identity
    if low_coverage and low_identity:
        return "low_coverage_and_identity"
    if low_coverage:
        return "low_coverage"
    if low_identity:
        return "low_identity"
    return None


# ---------------------------------------------------------------------------
# Batched lifter — single tblastn call for all proteins per genome
# ---------------------------------------------------------------------------

def run_tblastn_batch(
    proteins: List[Tuple[str, str]],
    query_genome: SeqRecord,
    tmp_dir: Path,
    evalue: float = 1e-5,
) -> Dict[str, List]:
    """
    Run tblastn for a batch of proteins against one genome, in a single call.

    Args:
        proteins: List of (query_id, protein_seq) tuples
        query_genome: Genome to search
        tmp_dir: Temp directory for FASTA + XML files
        evalue: E-value threshold

    Returns:
        Dict mapping query_id -> list of HSPs from best alignment
        (empty list if no hits for that query)
    """
    if not proteins:
        return {}

    # Multi-FASTA query containing all proteins
    prot_path = tmp_dir / "all_proteins.fa"
    with open(prot_path, "w") as f:
        for qid, seq in proteins:
            f.write(f">{qid}\n{seq}\n")

    # Genome FASTA written once
    genome_path = tmp_dir / "genome.fa"
    SeqIO.write(query_genome, str(genome_path), "fasta")

    result_path = tmp_dir / "all_blast.xml"
    cmd = [
        TBLASTN_BIN,
        "-query", str(prot_path),
        "-subject", str(genome_path),
        "-evalue", str(evalue),
        "-outfmt", "5",
        "-out", str(result_path),
        "-seg", "no",
        "-soft_masking", "false",
    ]

    try:
        subprocess.run(cmd, check=True, capture_output=True)
    except FileNotFoundError as exc:
        # tblastn binary not on PATH. This is an environment error, not a
        # biological "no hit" — surface it loudly instead of silently turning
        # every gene into no_hit.
        raise BlastNotInstalledError(
            "tblastn executable not found on PATH while running BLAST. "
            "Install NCBI BLAST+ and ensure 'tblastn' is on PATH."
        ) from exc
    except subprocess.CalledProcessError as exc:
        log_error("run_tblastn_batch", exc)
        return {qid: [] for qid, _ in proteins}

    if not result_path.exists() or result_path.stat().st_size == 0:
        return {qid: [] for qid, _ in proteins}

    hsps_by_id: Dict[str, List] = {qid: [] for qid, _ in proteins}
    with open(result_path) as f:
        try:
            for record in NCBIXML.parse(f):
                # record.query is the FASTA header (without ">")
                qid = record.query.split()[0]
                if record.alignments:
                    hsps_by_id[qid] = record.alignments[0].hsps
        except (ValueError, SyntaxError) as exc:
            # Malformed/truncated BLAST XML. Log it so a parse failure is not
            # mistaken for a genuine absence of hits; return whatever was
            # parsed before the error.
            log_error("run_tblastn_batch:NCBIXML.parse", exc)

    return hsps_by_id


# Tolerated number of unaligned N-terminal reference residues before we suspect the
# lifted start sits on an internal ATG and try to recover the true upstream start.
_N_TERM_TRUST_AA = 3


def _build_lifted_from_hsps(
    feature: Dict,
    hsps: List,
    query_record: SeqRecord,
    protein_length: int,
    min_coverage: float,
    min_identity: float,
    rescue_window: int,
    validate_codons: bool,
    ref_has_atg_start: bool = True,
) -> LiftedFeature:
    """
    Build a LiftedFeature from already-computed HSPs.

    This is the single post-tblastn path: merge HSPs → rescue stop codon (or
    extrapolate terminals for non-CDS) → coverage/identity gate → extract
    sequence → validate boundaries → rescue start codon. Used by the batched
    lifter for every feature.
    """
    base = dict(
        name=feature["name"],
        source_name=None,  # tblastn: query has no annotation, source_name not applicable
        ref_start=feature["start"],
        ref_end=feature["end"],
        strand=feature.get("strand", "+"),
        method="tblastn",
    )

    if not hsps:
        return LiftedFeature(
            **base,
            query_start=None, query_end=None,
            sequence=None, coverage=0.0,
            status="no_hit",
        )

    q_start, q_end, strand, coverage, identity, score = merge_hsps(
        hsps,
        protein_length=protein_length,
    )
    # Report the strand of the QUERY match, not the reference feature's strand.
    # They differ when the query genome is a reverse-complement deposit: the gene
    # is then on the query's minus strand even though the reference annotates it "+".
    base["strand"] = strand
    # Diagnostics: raw merged-HSP coords (pre-rescue) + unaligned N-terminal ref
    # residues. Purely informational, so never let it break a real lift: HSP
    # objects come from an external parser and may not expose these attributes.
    base["raw_start"] = q_start
    base["raw_end"] = q_end
    base["n_term_missing_aa"] = _count_n_term_missing_aa(hsps, strand)
    n_term_extension = 0
    c_term_extension = 0
    stop_rescue_offset_bp = None

    if validate_codons and q_end is not None:
        rescued_stop = rescue_stop_codon(
            query_record,
            q_start,
            q_end,
            strand,
            max_codons=30,
            expected_length=protein_length * 3 + 3,
        )
        if rescued_stop:
            q_start, q_end, _, codons_extended = rescued_stop
            stop_rescue_offset_bp = codons_extended * 3
        else:
            if strand == "+":
                q_end = min(q_end + 3, len(query_record.seq))
            else:
                q_start = max(1, q_start - 3)
    elif q_start is not None and q_end is not None:
        q_start, q_end, n_term_extension, c_term_extension = extrapolate_terminal_boundaries(
            start=q_start,
            end=q_end,
            strand=strand,
            hsps=hsps,
            protein_length=protein_length,
            genome_length=len(query_record.seq),
            coverage=coverage,
        )

    threshold_failure = _threshold_failure_status(
        coverage,
        identity,
        min_coverage,
        min_identity,
    )
    if threshold_failure:
        return LiftedFeature(
            **base,
            query_start=q_start, query_end=q_end,
            sequence=None, coverage=round(coverage, 4),
            status=threshold_failure,
            identity=round(identity * 100, 1),
            score=round(score, 1),
        )

    seq_str = extract_sequence(query_record, q_start, q_end, strand)

    if not validate_codons:
        status = "ok_extrapolated" if (n_term_extension or c_term_extension) else "ok"
        return LiftedFeature(
            **base,
            query_start=q_start, query_end=q_end,
            sequence=seq_str, coverage=round(coverage, 4),
            status=status,
            identity=round(identity * 100, 1),
            score=round(score, 1),
        )

    # A real CDS ends at its FIRST in-frame stop. If the lifted end over-ran that stop
    # (e.g. ORF1a read past the ORF1a/ORF1b frameshift stop, following a longer reference),
    # the sequence contains internal stops. Trim the end back to the first in-frame stop.
    # No-op for genes that already stop correctly, so it never regresses valid CDS.
    # Only trim when the start is already a valid ATG, so the reading frame is anchored and
    # the "first in-frame stop" is meaningful (genes whose start still needs rescue are
    # handled by the start-rescue path below).
    trim_len = first_inframe_stop_end(seq_str) if seq_str[:3] == "ATG" else None
    if trim_len is not None and trim_len < len(seq_str):
        if strand == "+":
            q_end = q_start + trim_len - 1
        else:
            q_start = q_end - trim_len + 1
        seq_str = extract_sequence(query_record, q_start, q_end, strand)

    validation = validate_cds_boundaries(seq_str)

    if not ref_has_atg_start:
        # The reference feature itself has no ATG start (a -1 PRF / frameshift gene such as
        # PRRSV ORF1b, annotated partial). There is no canonical start codon to rescue, so
        # keep the homology-derived (HSP) start -- its coordinates follow the reference's
        # convention -- instead of hunting an ATG and grabbing a wrong internal one. The end
        # is still trimmed to the first in-frame stop (a real ORF end). Generic: triggered
        # only by the reference annotation, no gene/virus name.
        fs_trim = first_inframe_stop_end(seq_str)
        if fs_trim is not None and fs_trim < len(seq_str):
            if strand == "+":
                q_end = q_start + fs_trim - 1
            else:
                q_start = q_end - fs_trim + 1
            seq_str = extract_sequence(query_record, q_start, q_end, strand)
            validation = validate_cds_boundaries(seq_str)
        return LiftedFeature(
            **base,
            query_start=q_start, query_end=q_end,
            sequence=seq_str, coverage=round(coverage, 4),
            status="ok_no_start_codon",
            has_start_codon=validation["has_start_codon"],
            has_stop_codon=validation["has_stop_codon"],
            in_frame=validation["in_frame"],
            identity=round(identity * 100, 1),
            score=round(score, 1),
        )

    # N-terminal reference residues the HSP failed to align. Derived ONLY from the
    # reference-vs-query alignment (HSP query coordinates are positions on the reference
    # protein), never from the truth annotation -- so this stays leakage-free. A large
    # value means the lifted start may sit on an INTERNAL in-frame ATG: a truncated CDS
    # that still validates. In that case we try to recover the true upstream start.
    missing_n_aa = base.get("n_term_missing_aa") or 0

    n_term_recovery_bp = None
    if validation["valid"] and missing_n_aa > _N_TERM_TRUST_AA:
        recovered = rescue_start_codon(
            query_record, q_start, q_end, strand,
            max_window=max(rescue_window, missing_n_aa * 3 + 30),
            expected_length=protein_length * 3 + 3,
        )
        if recovered:
            r_start, r_end, r_seq, r_offset = recovered
            r_val = validate_cds_boundaries(r_seq)
            expected = protein_length * 3 + 3
            # adopt the recovered start only if it is a valid CDS AND lands closer to the
            # reference-implied length (i.e. it genuinely restores the missing N-terminus,
            # rather than swapping one wrong start for another)
            if r_val["valid"] and abs(len(r_seq) - expected) < abs(len(seq_str) - expected):
                q_start, q_end, seq_str, validation = r_start, r_end, r_seq, r_val
                n_term_recovery_bp = r_offset

    if validation["valid"]:
        rescued = (stop_rescue_offset_bp is not None) or (n_term_recovery_bp is not None)
        actions = []
        if n_term_recovery_bp is not None:
            actions.append(f"start {n_term_recovery_bp:+d} bp (N-term)")
        if stop_rescue_offset_bp is not None:
            actions.append(f"stop +{stop_rescue_offset_bp} bp")
        return LiftedFeature(
            **base,
            query_start=q_start, query_end=q_end,
            sequence=seq_str, coverage=round(coverage, 4),
            status="ok_rescued" if rescued else "ok",
            has_start_codon=True, has_stop_codon=True,
            in_frame=validation["in_frame"],
            rescue_offset=(n_term_recovery_bp if n_term_recovery_bp is not None else stop_rescue_offset_bp),
            rescue_target=(
                "start+stop" if (n_term_recovery_bp is not None and stop_rescue_offset_bp is not None)
                else "start" if n_term_recovery_bp is not None
                else "stop" if stop_rescue_offset_bp is not None else None
            ),
            rescue_action="; ".join(actions) if actions else None,
            identity=round(identity * 100, 1),
            score=round(score, 1),
        )

    if not validation["has_start_codon"]:
        # When the HSP left the N-terminus unaligned (a divergent front), the true start
        # is ~missing_n_aa codons UPSTREAM of the lifted start. Anchor the ATG search
        # there so we recover it, instead of grabbing a downstream internal ATG (which
        # is what happens if we search around the HSP start). Uses only the ref-vs-query
        # alignment (missing_n_aa), never the truth annotation -> leakage-free.
        anchor_start, anchor_end = q_start, q_end
        search_window = rescue_window
        if missing_n_aa > _N_TERM_TRUST_AA:
            # Anchor on the RAW HSP coordinates. A stop-rescue ran earlier from the
            # (truncated) HSP start and may have extended the end downstream; scoring the
            # start against that extended end biases toward the internal ATG. The raw HSP
            # end reflects the true ORF extent, so the upstream true start wins. After the
            # start is recovered, the stop is re-checked/rescued below.
            raw_s = base.get("raw_start") or q_start
            raw_e = base.get("raw_end") or q_end
            if strand == "+":
                anchor_start = max(1, raw_s - missing_n_aa * 3)
                anchor_end = raw_e
            else:
                anchor_end = min(len(query_record.seq), raw_e + missing_n_aa * 3)
                anchor_start = raw_s
            search_window = missing_n_aa * 3 + 30
        rescued = rescue_start_codon(
            query_record,
            anchor_start,
            anchor_end,
            strand,
            max_window=search_window,
            expected_length=protein_length * 3 + 3,
        )
        if rescued:
            new_start, new_end, new_seq, offset = rescued
            revalidation = validate_cds_boundaries(new_seq)
            rescue_parts = [f"start {offset:+d} bp"]
            if not revalidation["has_stop_codon"]:
                rescued_stop = rescue_stop_codon(
                    query_record,
                    new_start,
                    new_end,
                    strand,
                    max_codons=30,
                    expected_length=protein_length * 3 + 3,
                )
                if rescued_stop:
                    new_start, new_end, new_seq, codons_extended = rescued_stop
                    rescue_parts.append(f"stop +{codons_extended * 3} bp")
                    revalidation = validate_cds_boundaries(new_seq)
            status = "ok_rescued" if revalidation["valid"] else "invalid_boundaries"
            return LiftedFeature(
                **base,
                query_start=new_start, query_end=new_end,
                sequence=new_seq, coverage=round(coverage, 4),
                status=status,
                has_start_codon=True,
                has_stop_codon=revalidation["has_stop_codon"],
                in_frame=revalidation["in_frame"],
                rescue_offset=offset,
                rescue_target="start+stop" if len(rescue_parts) > 1 else "start",
                rescue_action="; ".join(rescue_parts),
                identity=round(identity * 100, 1),
                score=round(score, 1),
            )

    return LiftedFeature(
        **base,
        query_start=q_start, query_end=q_end,
        sequence=seq_str, coverage=round(coverage, 4),
        status="invalid_boundaries",
        has_start_codon=validation["has_start_codon"],
        has_stop_codon=validation["has_stop_codon"],
        in_frame=validation["in_frame"],
        identity=round(identity * 100, 1),
        score=round(score, 1),
    )


def _fill_contiguous_gaps(
    results: List[LiftedFeature],
    ref_features: List[Dict],
    query_record: SeqRecord,
) -> List[LiftedFeature]:
    """Recover a feature that failed to lift by placing it in the gap between its nearest
    successfully-lifted neighbours -- but only when the REFERENCE annotates the features
    contiguously (a polyprotein: each feature touches the next, ref_end == next ref_start - 1).

    Mature peptides of a polyprotein are contiguous with no gaps, so a peptide too short or
    divergent for tblastn to hit (e.g. FMDV 2A, 18 aa) is bounded exactly by its neighbours.
    Uses only lifted-neighbour coordinates + the reference ordering -- never the truth
    annotation (leakage-free). Generic: any contiguous polyprotein, no gene/virus name.
    """
    FAIL = {"no_hit", "low_coverage", "low_identity", "low_coverage_and_identity", "translation_fail"}
    by_ref_end = {rf.get("end"): i for i, rf in enumerate(ref_features)}
    by_ref_start = {rf.get("start"): i for i, rf in enumerate(ref_features)}
    for i, lf in enumerate(results):
        if lf.query_start is not None and lf.query_end is not None:
            continue
        if lf.status not in FAIL:
            continue
        rf = ref_features[i]
        if rf.get("strand", "+") != "+":
            continue  # + strand polyproteins; - strand handled by the normal path
        li = by_ref_end.get(rf.get("start", 0) - 1)          # ref neighbour ending just before
        ri = by_ref_start.get(rf.get("end", 0) + 1)          # ref neighbour starting just after
        if li is None or ri is None:
            continue
        left, right = results[li], results[ri]
        if left.query_end is None or right.query_start is None:
            continue
        new_start, new_end = left.query_end + 1, right.query_start - 1
        if new_start > new_end:
            continue
        ref_len = rf["end"] - rf["start"] + 1
        if not (0.5 * ref_len <= (new_end - new_start + 1) <= 1.5 * ref_len):
            continue  # gap implausibly sized vs the reference peptide -> leave as failed
        lf.query_start, lf.query_end = new_start, new_end
        lf.sequence = extract_sequence(query_record, new_start, new_end, "+")
        lf.coverage = 1.0
        lf.status = "ok_gap_filled"
    return results


def lift_all_tblastn(
    ref_features: List[Dict],
    ref_record: SeqRecord,
    query_record: SeqRecord,
    min_coverage: float = 0.5,
    min_identity: float = 0.3,
    evalue: float = 1e-5,
    rescue_window: int = 200,
    validate_codons: bool = True,
) -> List[LiftedFeature]:
    """
    Lift all ref features onto query using tblastn — batched implementation.

    Translates all ref features to protein, runs ONE tblastn call with a
    multi-FASTA query against the genome, then dispatches HSPs back to
    each feature for merge + validation.

    This is significantly faster than per-feature tblastn calls because:
        - Genome FASTA written once (not N times)
        - tblastn subprocess started once (not N times)
        - tblastn indexes the genome once internally
    """
    results: List[LiftedFeature] = []

    # 0. Fail fast with a clear message if BLAST+ is not installed, rather than
    #    letting a FileNotFoundError surface deep inside subprocess.
    ensure_tblastn_available()

    # 1. Translate every feature; track which features were translatable.
    proteins: List[Tuple[str, str]] = []
    qid_to_protein_length: Dict[str, int] = {}
    failed: List[Dict] = []  # features that failed to translate

    for idx, feature in enumerate(ref_features):
        protein = translate_feature(feature, ref_record)
        if not protein:
            failed.append(feature)
            continue
        qid = f"feat_{idx}"
        proteins.append((qid, protein))
        qid_to_protein_length[qid] = len(protein)

    # 2. Single batched tblastn call for all translated proteins.
    with tempfile.TemporaryDirectory() as tmp:
        tmp_dir = Path(tmp)
        hsps_by_id = run_tblastn_batch(
            proteins, query_record, tmp_dir, evalue=evalue
        )

    # 3. Build LiftedFeature for each feature, preserving input order.
    failed_set = {id(f) for f in failed}
    qid_iter = iter(proteins)

    for feature in ref_features:
        if id(feature) in failed_set:
            base = dict(
                name=feature["name"],
                source_name=None,  # tblastn: query has no annotation, source_name not applicable
                ref_start=feature["start"],
                ref_end=feature["end"],
                strand=feature.get("strand", "+"),
                method="tblastn",
            )
            results.append(LiftedFeature(
                **base,
                query_start=None, query_end=None,
                sequence=None, coverage=0.0,
                status="translation_fail",
            ))
            continue

        qid, _ = next(qid_iter)
        hsps = hsps_by_id.get(qid, [])
        results.append(_build_lifted_from_hsps(
            feature=feature,
            hsps=hsps,
            query_record=query_record,
            protein_length=qid_to_protein_length[qid],
            min_coverage=min_coverage,
            min_identity=min_identity,
            rescue_window=rescue_window,
            validate_codons=validate_codons,
            ref_has_atg_start=ref_starts_with_atg(feature, ref_record),
        ))

    # Post-process: recover features that failed to lift (e.g. peptides too short for
    # tblastn) from the gap between contiguous, successfully-lifted neighbours in a
    # polyprotein. No-op unless the reference annotates features contiguously.
    results = _fill_contiguous_gaps(results, ref_features, query_record)
    return results


def process_one_query_record(
    ref_record: SeqRecord,
    query_record: SeqRecord,
    ref_cds: List[Dict],
    ref_feature_type: str,
    min_coverage: float,
    min_identity: float = 0.3,
    evalue: float = 1e-5,
    rescue_window: int = 200,
    quiet: bool = False,
) -> List[LiftedFeature]:
    """
    Lift all reference features onto one query genome via tblastn.

    Args:
        ref_record:       Reference genome record.
        query_record:     Query genome record.
        ref_cds:          Parsed reference features (alias-normalized).
        ref_feature_type: "CDS" or "mat_peptide" — controls codon validation.
        min_coverage:     Minimum accepted protein coverage.
        min_identity:     Minimum accepted protein identity.
        evalue:           E-value threshold for tblastn.
        rescue_window:    Window size (bp) for start codon rescue.
        quiet:            Suppress per-record console output.

    Returns:
        List of LiftedFeature objects.
    """
    # mat_peptide boundaries don't require ATG/stop codon validation
    validate_codons = (ref_feature_type == "CDS")

    return lift_all_tblastn(
        ref_features=ref_cds,
        ref_record=ref_record,
        query_record=query_record,
        min_coverage=min_coverage,
        min_identity=min_identity,
        evalue=evalue,
        rescue_window=rescue_window,
        validate_codons=validate_codons,
    )


def lift_feature_tblastn(
    feature: Dict,
    ref_record: SeqRecord,
    query_record: SeqRecord,
    tmp_dir: Optional[Path] = None,
    min_coverage: float = 0.5,
    min_identity: float = 0.3,
    evalue: float = 1e-5,
    rescue_window: int = 200,
    validate_codons: bool = True,
) -> LiftedFeature:
    """
    Lift a single ref feature onto the query genome via tblastn.

    Backward-compatible wrapper kept for external callers (e.g. validation
    notebooks). Internally this now delegates to the batched path, which is the
    single source of truth for the lift logic. `tmp_dir` is accepted but ignored
    — the batched path manages its own temporary directory.
    """
    return lift_all_tblastn(
        ref_features=[feature],
        ref_record=ref_record,
        query_record=query_record,
        min_coverage=min_coverage,
        min_identity=min_identity,
        evalue=evalue,
        rescue_window=rescue_window,
        validate_codons=validate_codons,
    )[0]
