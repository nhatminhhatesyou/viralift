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

def load_single_genbank(path):
    records = list(SeqIO.parse(str(path), "genbank"))
    if not records:
        raise ValueError(f"No records found in {path}")
    if len(records) > 1:
        raise ValueError(f"Expected 1 record, found {len(records)}")
    return records[0]


def load_genbank_records(path):
    return list(SeqIO.parse(str(path), "genbank"))


def _get_feature_name(feature):
    qualifiers = feature.qualifiers

    for key in ["gene", "product", "label", "standard_name", "note"]:
        values = qualifiers.get(key)
        if values:
            return values[0]

    return "unknown"


def _feature_to_dict(feature, genome_length, order_index):
    start = int(feature.location.start) + 1
    end = int(feature.location.end)
    strand = feature.location.strand
    length = end - start + 1

    return {
        "type": feature.type,
        "name": _get_feature_name(feature),
        "start": start,
        "end": end,
        "strand": strand,
        "length": length,
        "rel_start": start / genome_length,
        "rel_end": end / genome_length,
        "order": order_index,
    }


def parse_cds_features(record):
    items = []
    genome_length = len(record.seq)

    cds_features = [f for f in record.features if f.type == "CDS"]

    for i, feature in enumerate(cds_features, 1):
        items.append(_feature_to_dict(feature, genome_length, i))

    return items


def parse_mat_peptides(record):
    items = []
    genome_length = len(record.seq)

    mat_features = [f for f in record.features if f.type == "mat_peptide"]

    for i, feature in enumerate(mat_features, 1):
        items.append(_feature_to_dict(feature, genome_length, i))

    return items


def parse_feature_levels(record):
    return {
        "CDS": parse_cds_features(record),
        "mat_peptide": parse_mat_peptides(record),
    }