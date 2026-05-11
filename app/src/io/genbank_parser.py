from pathlib import Path
from typing import Dict, List, Optional

from Bio import SeqIO
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord


"""
Module: genbank_parser.py

Purpose:
    Parse GenBank records into structured Python dictionaries for downstream use.

Main responsibilities:
    - Load one or multiple GenBank records
    - Extract genome-level metadata
    - Parse CDS and mat_peptide features
    - Convert features into normalized dictionaries for matching / scoring

Notes:
    - This module is responsible for extracting raw annotation information.
    - It does NOT perform alias resolution or canonical biological naming.
    - Alias mapping should be added as a separate naming layer upstream or downstream.
"""


# ---------------------------------------------------------------------
# File loading utilities
# ---------------------------------------------------------------------

def load_single_genbank(path: Path) -> SeqRecord:
    """
    Load exactly one GenBank record from a file.

    Args:
        path: Path to a GenBank file

    Returns:
        A single SeqRecord

    Raises:
        ValueError: If the file contains zero or multiple records
    """
    records = list(SeqIO.parse(str(path), "genbank"))

    if not records:
        raise ValueError(f"No records found in {path}")

    if len(records) > 1:
        raise ValueError(f"Expected 1 record, found {len(records)} in {path}")

    return records[0]


def load_genbank_records(path: Path) -> List[SeqRecord]:
    """
    Load all GenBank records from a file.

    Args:
        path: Path to a GenBank file

    Returns:
        List of SeqRecord objects
    """
    return list(SeqIO.parse(str(path), "genbank"))


def get_record_sequence(record: SeqRecord) -> str:
    """
    Return the full nucleotide sequence of a record as a string.

    Args:
        record: Input SeqRecord

    Returns:
        Genome sequence as string
    """
    return str(record.seq)


# Qualifier keys to extract from each feature for alias lookup.
# Order matches _get_feature_name priority so downstream lookup tries
# the most-specific field first when resolving canonical names.
_LOOKUP_QUALIFIER_KEYS = ["gene", "product", "label", "standard_name", "locus_tag", "note"]


# ---------------------------------------------------------------------
# Qualifier helpers
# ---------------------------------------------------------------------

def _get_first_qualifier(feature: SeqFeature, key: str) -> Optional[str]:
    """
    Return the first qualifier value for a given key, if present.

    Args:
        feature: Biopython SeqFeature
        key: Qualifier key (e.g. "gene", "product")

    Returns:
        First qualifier value, or None if missing
    """
    values = feature.qualifiers.get(key)
    if not values:
        return None
    return values[0]


def _get_feature_name(feature: SeqFeature, fallback_index: Optional[int] = None) -> str:
    """
    Select a stable display name for a feature.

    Priority:
        gene > product > label > standard_name > locus_tag > note > fallback

    Args:
        feature: Biopython SeqFeature
        fallback_index: Optional index used to build a fallback name

    Returns:
        A human-readable feature name
    """
    for key in ["gene", "product", "label", "standard_name", "locus_tag", "note"]:
        value = _get_first_qualifier(feature, key)
        if value:
            return value

    if fallback_index is not None:
        return f"{feature.type}_{fallback_index:03d}"

    return "unknown"


# ---------------------------------------------------------------------
# Feature conversion helpers
# ---------------------------------------------------------------------

def _feature_to_basic_dict(feature: SeqFeature, index: int) -> Dict:
    """
    Convert a feature into a minimal dictionary for extraction workflows.

    Fields:
        - name
        - all keys in _LOOKUP_QUALIFIER_KEYS (gene, product, label, standard_name, locus_tag, note)
        - start
        - end
        - strand

    Args:
        feature: Biopython SeqFeature
        index: Feature index used for fallback naming

    Returns:
        Dictionary containing basic feature information
    """
    start = int(feature.location.start) + 1
    end = int(feature.location.end)

    strand = "+"
    if feature.location.strand == -1:
        strand = "-"

    result = {
        "name": _get_feature_name(feature, index),
        "start": start,
        "end": end,
        "strand": strand,
    }

    for key in _LOOKUP_QUALIFIER_KEYS:
        result[key] = _get_first_qualifier(feature, key)

    return result


def _feature_to_scored_dict(feature: SeqFeature, genome_length: int, order_index: int) -> Dict:
    """
    Convert a feature into a richer dictionary for structural matching/scoring.

    Fields:
        - type
        - name
        - all keys in _LOOKUP_QUALIFIER_KEYS (gene, product, label, standard_name, locus_tag, note)
        - start / end
        - strand
        - length
        - rel_start / rel_end
        - order

    Args:
        feature: Biopython SeqFeature
        genome_length: Total genome length
        order_index: Positional order of the feature in the genome

    Returns:
        Dictionary used for heuristic feature matching
    """
    start = int(feature.location.start) + 1
    end = int(feature.location.end)
    strand = "+" if feature.location.strand == 1 else "-"
    length = end - start + 1

    result = {
        "type":      feature.type,
        "name":      _get_feature_name(feature, order_index),
        "start":     start,
        "end":       end,
        "strand":    strand,
        "length":    length,
        "rel_start": start / genome_length,
        "rel_end":   end / genome_length,
        "order":     order_index,
    }

    for key in _LOOKUP_QUALIFIER_KEYS:
        result[key] = _get_first_qualifier(feature, key)

    return result


# ---------------------------------------------------------------------
# Feature parsers
# ---------------------------------------------------------------------

def parse_cds_features(record: SeqRecord) -> List[Dict]:
    """
    Parse CDS features into dictionaries for structural comparison / matching.

    Args:
        record: Input SeqRecord

    Returns:
        List of CDS feature dictionaries
    """
    items: List[Dict] = []
    genome_length = len(record.seq)

    cds_features = [feature for feature in record.features if feature.type == "CDS"]

    for i, feature in enumerate(cds_features, start=1):
        items.append(_feature_to_scored_dict(feature, genome_length, i))

    return items


def parse_mat_peptides(record: SeqRecord) -> List[Dict]:
    """
    Parse mat_peptide features into dictionaries for structural comparison / matching.

    Args:
        record: Input SeqRecord

    Returns:
        List of mat_peptide feature dictionaries
    """
    items: List[Dict] = []
    genome_length = len(record.seq)

    mat_features = [feature for feature in record.features if feature.type == "mat_peptide"]

    for i, feature in enumerate(mat_features, start=1):
        items.append(_feature_to_scored_dict(feature, genome_length, i))

    return items


def parse_feature_levels(record: SeqRecord) -> Dict[str, List[Dict]]:
    """
    Parse multiple feature levels from a record.

    Args:
        record: Input SeqRecord

    Returns:
        Dictionary grouping parsed features by type
    """
    return {
        "CDS": parse_cds_features(record),
        "mat_peptide": parse_mat_peptides(record),
    }


def parse_cds_features_basic(record: SeqRecord) -> List[Dict]:
    """
    Parse CDS features into minimal dictionaries for extraction workflows.

    This is useful when only coordinate and naming metadata are needed,
    without relative-position scoring fields.

    Args:
        record: Input SeqRecord

    Returns:
        List of CDS feature dictionaries
    """
    items: List[Dict] = []

    cds_features = [feature for feature in record.features if feature.type == "CDS"]

    for i, feature in enumerate(cds_features, start=1):
        items.append(_feature_to_basic_dict(feature, i))

    return items


# ---------------------------------------------------------------------
# Record metadata
# ---------------------------------------------------------------------

def get_record_metadata(record: SeqRecord) -> Dict:
    """
    Extract basic metadata from a GenBank record.

    Args:
        record: Input SeqRecord

    Returns:
        Dictionary containing high-level record metadata
    """
    organism = None

    for feature in record.features:
        if feature.type == "source":
            organism = _get_first_qualifier(feature, "organism")
            break

    return {
        "id": record.id,
        "name": record.name,
        "description": record.description,
        "organism": organism,
        "length": len(record.seq),
    }