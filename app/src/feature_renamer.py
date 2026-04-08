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

def compute_length_score(q, r):
    return min(q["length"], r["length"]) / max(q["length"], r["length"])


def compute_position_score(q, r):
    pos_diff = abs(q["rel_start"] - r["rel_start"]) + abs(q["rel_end"] - r["rel_end"])
    score = 1.0 - pos_diff
    return max(0.0, score)


def compute_order_score(q, r):
    order_diff = abs(q["order"] - r["order"])
    return 1.0 / (1.0 + order_diff)


def compute_strand_score(q, r):
    return 1.0 if q["strand"] == r["strand"] else 0.0


def score_feature_pair(q, r, use_strand=False):
    length_score = compute_length_score(q, r)
    position_score = compute_position_score(q, r)
    order_score = compute_order_score(q, r)

    score = (
        0.4 * length_score +
        0.4 * position_score +
        0.2 * order_score
    )

    if use_strand:
        strand_score = compute_strand_score(q, r)
        score = 0.35 * length_score + 0.35 * position_score + 0.15 * order_score + 0.15 * strand_score

    return {
        "score": round(score, 4),
        "length_score": round(length_score, 4),
        "position_score": round(position_score, 4),
        "order_score": round(order_score, 4),
        "strand_score": round(compute_strand_score(q, r), 4),
    }
    
def match_one_feature(query_feature, ref_features, threshold=0.75):
    if not ref_features:
        return {
            "query": query_feature,
            "best_ref": None,
            "score_details": None,
            "status": "no_reference_features",
        }

    candidates = []

    for ref_feature in ref_features:
        details = score_feature_pair(query_feature, ref_feature)

        candidates.append(
            {
                "ref_feature": ref_feature,
                "score_details": details,
            }
        )

    candidates.sort(key=lambda x: x["score_details"]["score"], reverse=True)
    best = candidates[0]

    best_score = best["score_details"]["score"]

    if best_score >= threshold:
        status = "matched"
    elif best_score >= 0.5:
        status = "ambiguous"
    else:
        status = "unresolved"

    return {
        "query": query_feature,
        "best_ref": best["ref_feature"],
        "score_details": best["score_details"],
        "status": status,
    }

def match_features(query_features, ref_features, threshold=0.75):
    results = []

    for query_feature in query_features:
        result = match_one_feature(query_feature, ref_features, threshold=threshold)
        results.append(result)

    return results