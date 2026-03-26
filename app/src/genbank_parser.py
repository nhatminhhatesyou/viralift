from pathlib import Path
from typing import Dict, List, Optional

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature


def load_single_genbank(gb_path: Path) -> SeqRecord:
    """Load a single GenBank record."""
    return SeqIO.read(str(gb_path), "genbank")


def get_record_sequence(record: SeqRecord) -> str:
    """Return the full nucleotide sequence as a string."""
    return str(record.seq)


def _get_first_qualifier(feature: SeqFeature, key: str) -> Optional[str]:
    """Return the first qualifier value if it exists."""
    values = feature.qualifiers.get(key)
    if not values:
        return None
    return values[0]


def _choose_feature_name(feature: SeqFeature, index: int) -> str:
    """
    Pick a stable name for a feature.
    Priority: gene > product > locus_tag > fallback name.
    """
    gene = _get_first_qualifier(feature, "gene")
    if gene:
        return gene

    product = _get_first_qualifier(feature, "product")
    if product:
        return product

    locus_tag = _get_first_qualifier(feature, "locus_tag")
    if locus_tag:
        return locus_tag

    return f"CDS_{index:03d}"


def parse_cds_features(record: SeqRecord) -> List[Dict]:
    """
    Extract CDS features from a GenBank record.

    Returns a list of dictionaries with:
    - name
    - gene
    - product
    - start
    - end
    - strand
    """
    cds_list: List[Dict] = []

    cds_index = 1
    for feature in record.features:
        if feature.type != "CDS":
            continue

        start = int(feature.location.start) + 1
        end = int(feature.location.end)
        strand = "+" if feature.location.strand == 1 else "-"

        gene = _get_first_qualifier(feature, "gene")
        product = _get_first_qualifier(feature, "product")
        name = _choose_feature_name(feature, cds_index)

        cds_list.append(
            {
                "name": name,
                "gene": gene,
                "product": product,
                "start": start,
                "end": end,
                "strand": strand,
            }
        )
        cds_index += 1

    return cds_list


def get_record_metadata(record: SeqRecord) -> Dict:
    """Extract simple metadata from a GenBank record."""
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
    
def load_genbank_records(gb_path: Path) -> List[SeqRecord]:
    """Load all GenBank records from a file."""
    return list(SeqIO.parse(str(gb_path), "genbank"))


def load_single_genbank(gb_path: Path) -> SeqRecord:
    """Load a single GenBank record."""
    return SeqIO.read(str(gb_path), "genbank")


def get_record_sequence(record: SeqRecord) -> str:
    """Return the full nucleotide sequence as a string."""
    return str(record.seq)


def _get_first_qualifier(feature: SeqFeature, key: str) -> Optional[str]:
    """Return the first qualifier value if it exists."""
    values = feature.qualifiers.get(key)
    if not values:
        return None
    return values[0]


def _choose_feature_name(feature: SeqFeature, index: int) -> str:
    """Pick a stable name for a feature."""
    gene = _get_first_qualifier(feature, "gene")
    if gene:
        return gene

    product = _get_first_qualifier(feature, "product")
    if product:
        return product

    locus_tag = _get_first_qualifier(feature, "locus_tag")
    if locus_tag:
        return locus_tag

    return f"CDS_{index:03d}"


def parse_cds_features(record: SeqRecord) -> List[Dict]:
    """Extract CDS features from a GenBank record."""
    cds_list: List[Dict] = []

    cds_index = 1
    for feature in record.features:
        if feature.type != "CDS":
            continue

        start = int(feature.location.start) + 1
        end = int(feature.location.end)
        strand = "+" if feature.location.strand == 1 else "-"

        gene = _get_first_qualifier(feature, "gene")
        product = _get_first_qualifier(feature, "product")
        name = _choose_feature_name(feature, cds_index)

        cds_list.append(
            {
                "name": name,
                "gene": gene,
                "product": product,
                "start": start,
                "end": end,
                "strand": strand,
            }
        )
        cds_index += 1

    return cds_list


def get_record_metadata(record: SeqRecord) -> Dict:
    """Extract simple metadata from a GenBank record."""
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