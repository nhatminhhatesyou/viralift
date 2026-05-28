from __future__ import annotations

import sys
from collections import Counter
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


ROOT = Path(__file__).resolve().parents[4]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from app.validation._shared.validation_utils import (  # noqa: E402
    CROSS_CHECK,
    DATA,
    load_genbank_records,
    load_reference_bundle,
    parse_truth_features,
)


HERE = Path(__file__).resolve().parents[1]
TIERS_DIR = HERE / "tiers"

DATASETS = [
    ("PRRS", DATA / "PRRS_ref_test.gb", CROSS_CHECK / "PRRS_100seq_anno.gb"),
    ("FMD", DATA / "FMD_ref_test.gb", CROSS_CHECK / "FMD_100seq_anno.gb"),
]

TIER_NAMES = ("strict_clean", "usable_clean", "raw_realworld")
CANONICAL_SOURCES = {"alias", "alias_conflict_resolved"}


def truth_quality(record: SeqRecord, bundle: Dict) -> Dict:
    truth, subfeatures = parse_truth_features(
        record,
        bundle["alias_lookup"],
        bundle["feature_type"],
    )
    bad_features = [
        feature
        for feature in truth
        if feature.get("name_source") not in CANONICAL_SOURCES
    ]
    name_counts = Counter(feature.get("name") for feature in truth)
    duplicate_names = sorted(
        name
        for name, count in name_counts.items()
        if name and count > 1
    )
    return {
        "feature_type": bundle["feature_type"],
        "truth": truth,
        "subfeatures": subfeatures,
        "truth_count": len(truth),
        "subfeature_count": len(subfeatures),
        "canonical_count": len(truth) - len(bad_features),
        "bad_features": bad_features,
        "duplicate_names": duplicate_names,
    }


def keep_for_tier(tier: str, quality: Dict) -> Tuple[bool, str]:
    if quality["truth_count"] == 0:
        return False, "no_truth_features"

    if tier == "raw_realworld":
        return True, "kept_raw"

    if quality["bad_features"]:
        statuses = sorted({
            feature.get("name_source", "unknown")
            for feature in quality["bad_features"]
        })
        return False, "noncanonical_truth:" + ",".join(statuses)

    if tier == "strict_clean" and quality["duplicate_names"]:
        return False, "duplicate_canonical_names"

    return True, "kept_clean"


def manifest_row(virus: str, tier: str, record: SeqRecord, quality: Dict) -> Dict:
    keep, reason = keep_for_tier(tier, quality)
    bad_values = [
        f"{feature.get('raw_name') or feature.get('name')}[{feature.get('name_source')}]"
        for feature in quality["bad_features"]
    ]
    return {
        "virus": virus,
        "tier": tier,
        "record_id": record.id,
        "keep": keep,
        "reason": reason,
        "feature_type": quality["feature_type"],
        "truth_count": quality["truth_count"],
        "canonical_count": quality["canonical_count"],
        "bad_count": len(quality["bad_features"]),
        "duplicate_names": ";".join(quality["duplicate_names"]),
        "bad_truth_values": ";".join(bad_values),
    }


def write_records(path: Path, records: Sequence[SeqRecord]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(list(records), str(path), "genbank")


def build_tiers() -> None:
    all_manifest_rows: List[Dict] = []
    tier_records: Dict[Tuple[str, str], List[SeqRecord]] = {
        (tier, virus): []
        for tier in TIER_NAMES
        for virus, _, _ in DATASETS
    }

    for virus, ref_path, query_path in DATASETS:
        bundle = load_reference_bundle(ref_path)
        records = load_genbank_records(query_path)
        for record in records:
            quality = truth_quality(record, bundle)
            for tier in TIER_NAMES:
                row = manifest_row(virus, tier, record, quality)
                all_manifest_rows.append(row)
                if row["keep"]:
                    tier_records[(tier, virus)].append(record)

    manifest = pd.DataFrame(all_manifest_rows)

    for tier in TIER_NAMES:
        tier_dir = TIERS_DIR / tier
        inputs_dir = tier_dir / "inputs"
        outputs_dir = tier_dir / "outputs"
        inputs_dir.mkdir(parents=True, exist_ok=True)
        outputs_dir.mkdir(parents=True, exist_ok=True)
        (outputs_dir / ".gitkeep").touch()

        tier_manifest = manifest[manifest["tier"] == tier].copy()
        tier_manifest.to_csv(inputs_dir / "manifest.tsv", sep="\t", index=False)

        summary = (
            tier_manifest.groupby(["virus", "reason"], dropna=False)
            .agg(records=("record_id", "size"), kept=("keep", "sum"))
            .reset_index()
            .sort_values(["virus", "records"], ascending=[True, False])
        )
        summary.to_csv(inputs_dir / "summary.tsv", sep="\t", index=False)

        for virus, _, _ in DATASETS:
            write_records(inputs_dir / f"{virus}_records.gb", tier_records[(tier, virus)])

    manifest.to_csv(HERE / "tier_manifest_all.tsv", sep="\t", index=False)


if __name__ == "__main__":
    build_tiers()
