"""
Module: run_logger.py

Purpose:
    Centralised logger for ViraLift.

    Two categories of events are recorded:

    ALIAS AUDIT — permanent record of every change to alias config files.
        Adding a new mapping (raw_name → canonical) or a new canonical key
        is a scientific decision that should be traceable over time.

    RUN SUMMARY — one entry per pipeline execution.
        Records which ref/query was used, thresholds, and the final
        status distribution (ok / rescued / failed / etc.).

Log file:
    logs/viralift.log  (rotating, max 5 MB × 3 backups)
    The directory is volume-mounted in Docker so logs survive rebuilds.

Format:
    2026-04-28 10:23:45 | INFO  | ALIAS_ADDED     | file=prrsv_alias.json | "GP-5 protein" → GP5
    2026-04-28 10:25:01 | INFO  | RUN_COMPLETE    | ref=AY150564.1 | records=12 | ok=84 rescued=3 ...
"""

import logging
from logging.handlers import RotatingFileHandler
from pathlib import Path
from typing import Dict, List, Optional


# ── configuration ────────────────────────────────────────────────────────────

LOG_DIR  = Path(__file__).parents[4] / "logs"   # repo_root/logs/
LOG_FILE = LOG_DIR / "viralift.log"

_MAX_BYTES  = 5 * 1024 * 1024   # 5 MB per file
_BACKUP_COUNT = 3                # keep 3 rotated files


# ── internal setup ───────────────────────────────────────────────────────────

def _setup_logger() -> logging.Logger:
    """
    Build and return the ViraLift logger.
    Safe to call multiple times — handlers are only added once.
    """
    logger = logging.getLogger("viralift")

    if logger.handlers:           # already initialised
        return logger

    logger.setLevel(logging.DEBUG)

    # ensure log dir exists (fails silently if can't create — e.g. read-only FS)
    try:
        LOG_DIR.mkdir(parents=True, exist_ok=True)
    except OSError:
        pass

    fmt = logging.Formatter(
        "%(asctime)s | %(levelname)-5s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    # rotating file handler
    try:
        fh = RotatingFileHandler(
            LOG_FILE,
            maxBytes=_MAX_BYTES,
            backupCount=_BACKUP_COUNT,
            encoding="utf-8",
        )
        fh.setFormatter(fmt)
        logger.addHandler(fh)
    except OSError:
        # if we can't write (e.g. container permission issue), fall back to stderr
        sh = logging.StreamHandler()
        sh.setFormatter(fmt)
        logger.addHandler(sh)

    return logger


_logger = _setup_logger()


# ── public API ───────────────────────────────────────────────────────────────

def log_alias_added(alias_config_name: str, raw_name: str, canonical: str) -> None:
    """
    Record that a raw_name → canonical mapping was added to an alias config.

    Called when the user saves a resolver decision permanently.
    """
    _logger.info(
        "ALIAS_ADDED | file=%s | %r → %s",
        alias_config_name, raw_name, canonical,
    )


def log_session_decisions(
    decisions: Dict[str, Optional[str]],
    saved_names: List[str],
) -> None:
    """
    Record every resolver decision made in the resolve stage.

    Logs both session-only decisions (not saved to alias config) and
    permanent ones (saved to alias config), so there is always a trace
    of what was mapped in a given run regardless of whether 💾 Save
    was checked.

    Args:
        decisions:   {raw_name: canonical_or_None} — full resolver output.
                     None means the user chose to ignore the name.
        saved_names: List of raw_names whose mappings were also written
                     permanently to the alias config file.
    """
    saved_set = set(saved_names)
    for raw_name, canonical in decisions.items():
        if canonical is None:
            _logger.info("SESSION_DECISION | %r → (ignored, session only)", raw_name)
        elif raw_name in saved_set:
            _logger.info(
                "SESSION_DECISION | %r → %s  [saved to config]", raw_name, canonical
            )
        else:
            _logger.info(
                "SESSION_DECISION | %r → %s  [session only, not saved]",
                raw_name, canonical,
            )


def log_canonical_added(alias_config_name: str, canonical: str) -> None:
    """
    Record that a brand-new canonical key was added to an alias config.

    Called when a ref gene name not in the DB is saved as a new canonical.
    """
    _logger.info(
        "CANONICAL_ADDED | file=%s | new_canonical=%r",
        alias_config_name, canonical,
    )


def log_run_start(
    ref_id: str,
    n_queries: int,
    min_coverage: float,
    min_identity: float,
    evalue: float,
    rescue_window: int,
    virus_name: Optional[str] = None,
    alias_config: Optional[str] = None,
) -> None:
    """
    Record the start of a pipeline run with its configuration.
    """
    virus_str  = virus_name  or "unknown"
    config_str = alias_config or "none"
    _logger.info(
        "RUN_START | ref=%s | queries=%d | virus=%s | alias=%s "
        "| min_cov=%.2f min_id=%.2f evalue=%s rescue=%d",
        ref_id, n_queries, virus_str, config_str,
        min_coverage, min_identity, evalue, rescue_window,
    )


def log_run_complete(
    ref_id: str,
    n_queries: int,
    summary: Dict[str, int],
) -> None:
    """
    Record the final status distribution after a pipeline run completes.

    Args:
        ref_id:    Reference record ID.
        n_queries: Number of query records processed.
        summary:   Dict from summarize_counts() — {status: count}.
    """
    total = sum(summary.values())
    _logger.info(
        "RUN_COMPLETE | ref=%s | records=%d | features=%d "
        "| ok=%d rescued=%d invalid=%d low_cov=%d no_hit=%d translation_fail=%d",
        ref_id, n_queries, total,
        summary.get("ok", 0),
        summary.get("ok_rescued", 0),
        summary.get("invalid_boundaries", 0),
        summary.get("low_coverage", 0),
        summary.get("no_hit", 0),
        summary.get("translation_fail", 0),
    )


def log_error(context: str, error: Exception) -> None:
    """
    Record an unexpected error with context description.
    """
    _logger.error("ERROR | %s | %s: %s", context, type(error).__name__, error)


def log_warning(message: str) -> None:
    """
    Record a non-fatal warning.
    """
    _logger.warning("WARNING | %s", message)
