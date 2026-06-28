# Auto-split from the original monolithic streamlit_app.py.
import streamlit as st
from pathlib import Path
from app.src.lifting.base import PASS_STATUSES, STATUS_META


ROOT = Path(__file__).parent.parent


REGISTRY_PATH  = ROOT / "app/config/virus_alias_registry.json"


CONFIG_DIR = ROOT / "app/config"


GOOD_STATUSES = set(PASS_STATUSES)


REVIEW_STATUSES = {
    status for status, meta in STATUS_META.items()
    if meta.get("category") == "review"
}


STATUS_LABEL = {
    status: str(meta.get("label", status))
    for status, meta in STATUS_META.items()
}


STATUS_TONE = {
    status: str(meta.get("category", "review"))
    for status, meta in STATUS_META.items()
}


def _init_state():
    defaults = dict(
        stage="upload",
        app_mode="Run pipeline",
        tmp=None,               # tempfile.TemporaryDirectory object
        ref_record=None,
        query_records=None,
        ref_features=None,
        use_ref_names=False,    # if True, output uses ref's raw names instead of alias keys
        ref_feature_type=None,
        alias_lookup={},
        alias_config_path=None,
        virus_name=None,
        canonical_list=[],
        unknown_ref_names=[],   # [raw_name] ref features not in alias DB
        unknown_names={},       # {raw_name: [record_id, ...]}
        resolver={},            # {raw_name: canonical or "-- ignore --"}
        all_results=None,       # [(query_id, [LiftedFeature])]
        run_errors=[],          # [{record_id, error}]
        min_coverage=0.5,
        min_identity=0.3,
        evalue=1e-5,
        rescue_window=200,
        bootstrap_alias_config=None,
        bootstrap_alias_config_path=None,
        bootstrap_suggestions=[],
        bootstrap_diagnostics={},
        bootstrap_virus_name="",
        bootstrap_keywords="",
        alias_manager_config_path=None,
        virus_review_metadata=[],
        virus_review_guess="",
        approved_export_features=[],
    )
    for k, v in defaults.items():
        if k not in st.session_state:
            st.session_state[k] = v


def _reset():
    for k in list(st.session_state.keys()):
        del st.session_state[k]
    _init_state()
