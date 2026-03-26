from typing import Dict, List

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def extract_feature_sequence(record: SeqRecord, feature: Dict) -> Seq:
    """Extract one feature sequence from a query record."""
    seq = record.seq[feature["start"] - 1: feature["end"]]

    if feature["strand"] == "-":
        seq = seq.reverse_complement()

    return seq


def extract_all_direct(record: SeqRecord, features: List[Dict]) -> List[Dict]:
    """Extract all CDS features directly from query coordinates."""
    results = []

    for feature in features:
        seq = extract_feature_sequence(record, feature)

        results.append(
            {
                "name": feature["name"],
                "gene": feature.get("ref_gene"),
                "product": feature.get("ref_product"),
                "start": feature["start"],
                "end": feature["end"],
                "strand": feature["strand"],
                "sequence": str(seq),
                "method": "direct_annotation",
                "status": "ok",
            }
        )

    return results