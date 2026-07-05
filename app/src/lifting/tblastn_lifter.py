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

from app.src.lifting.validator import validate_cds_boundaries, rescue_start_codon, rescue_stop_codon
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

def _hsp_to_genome_coords(hsp) -> Tuple[int, int, str]:
    """
    Convert a tblastn HSP to 1-based genome coordinates and strand.

    tblastn sbjct coords are 1-based. Strand determined by sbjct_start vs sbjct_end.
    """
    if hsp.sbjct_start <= hsp.sbjct_end:
        return hsp.sbjct_start, hsp.sbjct_end, "+"
    else:
        return hsp.sbjct_end, hsp.sbjct_start, "-"


def merge_hsps(hsps: List, protein_length: int) -> Tuple[int, int, str, float, float, float]:
    """
    Merge a list of HSPs into a single genomic span.

    Strategy:
        - Determine strand from majority vote
        - Take min(start) and max(end) across all HSPs on that strand
        - Compute weighted average identity (by aligned length)

    Args:
        hsps: tblastn HSPs for one protein query.
        protein_length: Full length of the reference protein query in amino acids.

    Returns: (merged_start, merged_end, strand, coverage_fraction, identity, bit_score)
    where coverage_fraction is relative to the full query protein length.
    """
    if not hsps:
        raise ValueError("No HSPs to merge")

    # Strand vote
    plus_len = sum(abs(h.sbjct_end - h.sbjct_start) + 1 for h in hsps if h.sbjct_start <= h.sbjct_end)
    minus_len = sum(abs(h.sbjct_end - h.sbjct_start) + 1 for h in hsps if h.sbjct_start > h.sbjct_end)
    strand = "+" if plus_len >= minus_len else "-"

    relevant = [h for h in hsps if (h.sbjct_start <= h.sbjct_end) == (strand == "+")]
    if not relevant:
        relevant = hsps

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

    relevant = [h for h in hsps if (h.sbjct_start <= h.sbjct_end) == (strand == "+")]
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


def _build_lifted_from_hsps(
    feature: Dict,
    hsps: List,
    query_record: SeqRecord,
    protein_length: int,
    min_coverage: float,
    min_identity: float,
    rescue_window: int,
    validate_codons: bool,
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

    validation = validate_cds_boundaries(seq_str)

    if validation["valid"]:
        rescued_by_stop = stop_rescue_offset_bp is not None
        return LiftedFeature(
            **base,
            query_start=q_start, query_end=q_end,
            sequence=seq_str, coverage=round(coverage, 4),
            status="ok_rescued" if rescued_by_stop else "ok",
            has_start_codon=True, has_stop_codon=True,
            in_frame=validation["in_frame"],
            rescue_offset=stop_rescue_offset_bp,
            rescue_target="stop" if rescued_by_stop else None,
            rescue_action=(
                f"stop +{stop_rescue_offset_bp} bp"
                if rescued_by_stop else None
            ),
            identity=round(identity * 100, 1),
            score=round(score, 1),
        )

    if not validation["has_start_codon"]:
        rescued = rescue_start_codon(
            query_record,
            q_start,
            q_end,
            strand,
            max_window=rescue_window,
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
        ))

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
