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

Advantages over minimap:
    - Protein is ~3-4x more conserved than nucleotide
    - Works across serotypes and lineages
    - Each gene is searched independently → no interference between overlapping genes

Limitations:
    - Frameshift genes (ORF1b in PRRSV) still have ambiguous start
    - Requires BLAST+ installed (tblastn in PATH)
    - Slower than minimap (~2-5s per genome vs ~0.1s)
"""

import subprocess
import tempfile
from io import StringIO
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from app.src.annotation.validator import validate_cds_boundaries, rescue_start_codon, rescue_stop_codon
from app.src.lifting.base import LiftedFeature


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
    except Exception:
        return None

    if len(protein) < 10:
        return None

    return protein


# ---------------------------------------------------------------------------
# BLAST runner
# ---------------------------------------------------------------------------

def run_tblastn(
    protein_seq: str,
    query_genome: SeqRecord,
    tmp_dir: Path,
    feature_name: str = "feature",
    evalue: float = 1e-5,
) -> Optional[List]:
    """
    Run tblastn of one protein against one genome.

    Returns list of Biopython HSP objects from the best alignment,
    or None if no hit found.
    """
    # Write protein FASTA
    prot_path = tmp_dir / f"{feature_name}_prot.fa"
    prot_path.write_text(f">{feature_name}\n{protein_seq}\n")

    # Write genome FASTA
    genome_path = tmp_dir / "genome.fa"
    SeqIO.write(query_genome, str(genome_path), "fasta")

    # Run tblastn
    result_path = tmp_dir / f"{feature_name}_blast.xml"
    cmd = [
        "tblastn",
        "-query", str(prot_path),
        "-subject", str(genome_path),
        "-evalue", str(evalue),
        "-outfmt", "5",       # XML output
        "-out", str(result_path),
        "-seg", "no",         # disable low-complexity filter for short viral proteins
        "-soft_masking", "false",
    ]

    try:
        subprocess.run(cmd, check=True, capture_output=True)
    except subprocess.CalledProcessError:
        return None

    if not result_path.exists() or result_path.stat().st_size == 0:
        return None

    with open(result_path) as f:
        try:
            records = list(NCBIXML.parse(f))
        except Exception:
            return None

    if not records or not records[0].alignments:
        return None

    best_alignment = records[0].alignments[0]
    return best_alignment.hsps


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


def merge_hsps(hsps: List) -> Tuple[int, int, str, float, float]:
    """
    Merge a list of HSPs into a single genomic span.

    Strategy:
        - Determine strand from majority vote
        - Take min(start) and max(end) across all HSPs on that strand
        - Compute weighted average identity (by aligned length)

    Returns: (merged_start, merged_end, strand, coverage_fraction, identity)
    where coverage_fraction is relative to the query protein length (in aa → *3 for nt).
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

    # Coverage: total query aa covered / query length
    query_positions = set()
    for h in relevant:
        qs = min(h.query_start, h.query_end)
        qe = max(h.query_start, h.query_end)
        query_positions.update(range(qs, qe + 1))
    query_length = max(h.query_end for h in relevant) if relevant else 1
    coverage = len(query_positions) / query_length

    bit_score = max(h.bits for h in relevant)

    return merged_start, merged_end, strand, coverage, identity, bit_score


# ---------------------------------------------------------------------------
# Sequence extraction + validation
# ---------------------------------------------------------------------------

def extract_sequence(query_record: SeqRecord, start: int, end: int, strand: str) -> str:
    """Extract 1-based inclusive sequence from query genome."""
    seq = query_record.seq[start - 1: end]
    if strand == "-":
        seq = seq.reverse_complement()
    return str(seq)


# ---------------------------------------------------------------------------
# Main lifter
# ---------------------------------------------------------------------------

def lift_feature_tblastn(
    feature: Dict,
    ref_record: SeqRecord,
    query_record: SeqRecord,
    tmp_dir: Path,
    min_coverage: float = 0.5,
    min_identity: float = 0.3,
    evalue: float = 1e-5,
    rescue_window: int = 50,
    validate_codons: bool = True,
) -> LiftedFeature:
    """
    Lift one ref CDS feature onto query genome using tblastn.

    Returns a LiftedFeature with status:
        ok               — lifted, valid start/stop
        ok_rescued       — lifted, ATG rescued within rescue_window
        invalid_boundaries — lifted but missing start/stop
        low_coverage     — hit found but protein coverage too low
        no_hit           — tblastn found no hit
        translation_fail — ref feature could not be translated
    """
    base = dict(
        name=feature["name"],
        canonical_name=feature.get("canonical_name"),
        ref_start=feature["start"],
        ref_end=feature["end"],
        strand=feature.get("strand", "+"),
        method="tblastn",
    )

    # 1. Translate ref feature
    protein = translate_feature(feature, ref_record)
    if not protein:
        return LiftedFeature(
            **base,
            query_start=None, query_end=None,
            sequence=None, coverage=0.0,
            status="translation_fail",
        )

    # 2. Run tblastn
    safe_name = feature["name"].replace(" ", "_").replace("/", "_")[:30]
    hsps = run_tblastn(protein, query_record, tmp_dir,
                       feature_name=safe_name, evalue=evalue)

    if not hsps:
        return LiftedFeature(
            **base,
            query_start=None, query_end=None,
            sequence=None, coverage=0.0,
            status="no_hit",
        )

    # 3. Merge HSPs → genomic coordinates
    q_start, q_end, strand, coverage, identity, score = merge_hsps(hsps)

    # tblastn aligns protein → stop codon not included in HSP end.
    # Use stop codon rescue: scan forward in-frame from q_end to find
    # the actual stop codon. Falls back gracefully if not found within window.
    if validate_codons and q_end is not None:
        rescued_stop = rescue_stop_codon(
            query_record, q_start, q_end, strand, max_codons=30
        )
        if rescued_stop:
            q_end, _, _ = rescued_stop
        else:
            # Fallback: at least add 3bp for stop codon
            q_end = min(q_end + 3, len(query_record.seq))

    if coverage < min_coverage or identity < min_identity:
        return LiftedFeature(
            **base,
            query_start=q_start, query_end=q_end,
            sequence=None, coverage=round(coverage, 4),
            status="low_coverage",
            identity=round(identity * 100, 1),
            score=round(score, 1),
        )

    # 4. Extract sequence
    seq_str = extract_sequence(query_record, q_start, q_end, strand)

    # 5. Skip codon validation if not requested (e.g. mat_peptide features)
    if not validate_codons:
        return LiftedFeature(
            **base,
            query_start=q_start, query_end=q_end,
            sequence=seq_str, coverage=round(coverage, 4),
            status="ok",
            identity=round(identity * 100, 1),
            score=round(score, 1),
        )

    # 6. Validate boundaries
    validation = validate_cds_boundaries(seq_str)

    if validation["valid"]:
        return LiftedFeature(
            **base,
            query_start=q_start, query_end=q_end,
            sequence=seq_str, coverage=round(coverage, 4),
            status="ok",
            has_start_codon=True, has_stop_codon=True,
            identity=round(identity * 100, 1),
            score=round(score, 1),
        )

    # 7. Try start codon rescue if ATG missing
    if not validation["has_start_codon"]:
        rescued = rescue_start_codon(
            query_record, q_start, q_end, strand, max_window=rescue_window
        )
        if rescued:
            new_start, new_seq, offset = rescued
            revalidation = validate_cds_boundaries(new_seq)
            status = "ok_rescued" if revalidation["has_stop_codon"] else "invalid_boundaries"
            return LiftedFeature(
                **base,
                query_start=new_start, query_end=q_end,
                sequence=new_seq, coverage=round(coverage, 4),
                status=status,
                has_start_codon=True,
                has_stop_codon=revalidation["has_stop_codon"],
                rescue_offset=offset,
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
    rescue_window: int = 50,
    validate_codons: bool = True,
) -> List[LiftedFeature]:
    """
    Lift all ref features onto query using tblastn.
    Uses a single shared temp directory per call.
    """
    results = []
    with tempfile.TemporaryDirectory() as tmp:
        tmp_dir = Path(tmp)
        for feature in ref_features:
            lifted = lift_feature_tblastn(
                feature, ref_record, query_record, tmp_dir,
                min_coverage=min_coverage,
                min_identity=min_identity,
                evalue=evalue,
                rescue_window=rescue_window,
                validate_codons=validate_codons,
            )
            results.append(lifted)
    return results
