from typing import Dict, List


def rename_query_cds_by_reference_order(ref_cds: List[Dict], query_cds: List[Dict]) -> List[Dict]:
    """Rename query CDS features using the order of reference CDS features."""
    renamed_features = []

    for ref_feature, query_feature in zip(ref_cds, query_cds):
        new_feature = query_feature.copy()
        new_feature["name"] = ref_feature["name"]
        new_feature["ref_name"] = ref_feature["name"]
        new_feature["ref_gene"] = ref_feature.get("gene")
        new_feature["ref_product"] = ref_feature.get("product")
        renamed_features.append(new_feature)

    return renamed_features