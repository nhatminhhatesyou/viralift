import os
from dataclasses import dataclass
from typing import Mapping, Optional


def _env_bool(value: Optional[str], default: bool = False) -> bool:
    if value is None:
        return default
    return str(value).strip().lower() in {"1", "true", "yes", "on"}


@dataclass(frozen=True)
class LLMConfig:
    """Runtime config for optional LLM review."""

    enabled: bool = False
    api_key: Optional[str] = None
    model: str = "gpt-5.4-mini"
    fallback_model: Optional[str] = "gpt-5.4"
    max_rows: int = 20
    timeout_seconds: int = 45
    endpoint: str = "https://api.openai.com/v1/responses"

    @classmethod
    def from_env(cls, env: Optional[Mapping[str, str]] = None) -> "LLMConfig":
        values = env or os.environ
        return cls(
            enabled=_env_bool(values.get("VIRALIFT_LLM_ENABLED"), default=False),
            api_key=values.get("OPENAI_API_KEY") or values.get("VIRALIFT_OPENAI_API_KEY"),
            model=values.get("VIRALIFT_LLM_MODEL", "gpt-5.4-mini"),
            fallback_model=values.get("VIRALIFT_LLM_FALLBACK_MODEL", "gpt-5.4") or None,
            max_rows=max(1, int(values.get("VIRALIFT_LLM_MAX_ROWS", "20"))),
            timeout_seconds=max(5, int(values.get("VIRALIFT_LLM_TIMEOUT_SECONDS", "45"))),
            endpoint=values.get(
                "VIRALIFT_LLM_ENDPOINT",
                "https://api.openai.com/v1/responses",
            ),
        )

    @property
    def available(self) -> bool:
        return bool(self.enabled and self.api_key and not _is_placeholder_key(self.api_key))


def _is_placeholder_key(value: str) -> bool:
    normalized = str(value or "").strip().lower()
    return (
        not normalized
        or "replace_with" in normalized
        or "your_new_key" in normalized
        or normalized in {"sk-...", "..."}
    )
