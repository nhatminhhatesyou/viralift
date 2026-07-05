import json
import urllib.error
import urllib.request
from typing import Dict, Iterable, List, Optional

from app.src.llm.config import LLMConfig


class LLMProviderError(RuntimeError):
    """Raised when an LLM provider cannot return a usable review."""


class NoopLLMProvider:
    """Provider used when LLM review is disabled or unavailable."""

    def review_alias_suggestions(self, payload: Dict) -> Dict:
        return {"reviews": []}


class MockLLMProvider:
    """Tiny test provider that returns precomputed reviews."""

    def __init__(self, reviews: Iterable[Dict]):
        self.reviews = list(reviews)
        self.calls: List[Dict] = []

    def review_alias_suggestions(self, payload: Dict) -> Dict:
        self.calls.append(payload)
        return {"reviews": self.reviews}


class OpenAILLMProvider:
    """
    Minimal OpenAI Responses API client using stdlib HTTP.

    The feature is optional, so ViraLift avoids adding a hard SDK dependency.
    If the API shape changes, VIRALIFT_LLM_ENDPOINT can point to a compatible
    wrapper without changing the alias review code.
    """

    def __init__(self, config: LLMConfig):
        self.config = config

    def review_alias_suggestions(self, payload: Dict) -> Dict:
        models = [self.config.model]
        if self.config.fallback_model and self.config.fallback_model not in models:
            models.append(self.config.fallback_model)

        last_error: Optional[Exception] = None
        for model in models:
            try:
                return self._request(model, payload)
            except LLMProviderError as exc:
                last_error = exc
                continue

        raise LLMProviderError(str(last_error) if last_error else "LLM review failed")

    def _request(self, model: str, payload: Dict) -> Dict:
        if not self.config.api_key:
            raise LLMProviderError("OPENAI_API_KEY is not configured")

        request_body = {
            "model": model,
            "input": [
                {
                    "role": "system",
                    "content": (
                        "You review viral gene alias suggestions. Return only JSON "
                        "matching the schema. Do not invent canonicals outside the "
                        "provided available_canonicals list. You are advisory only; "
                        "the user will approve or skip each suggestion."
                    ),
                },
                {"role": "user", "content": json.dumps(payload, ensure_ascii=False)},
            ],
            "text": {
                "format": {
                    "type": "json_schema",
                    "name": "viralift_alias_reviews",
                    "schema": _alias_review_schema(),
                    "strict": True,
                }
            },
        }

        data = json.dumps(request_body).encode("utf-8")
        request = urllib.request.Request(
            self.config.endpoint,
            data=data,
            method="POST",
            headers={
                "Authorization": f"Bearer {self.config.api_key}",
                "Content-Type": "application/json",
            },
        )
        try:
            with urllib.request.urlopen(request, timeout=self.config.timeout_seconds) as response:
                response_data = json.loads(response.read().decode("utf-8"))
        except urllib.error.HTTPError as exc:
            body = exc.read().decode("utf-8", errors="replace")
            raise LLMProviderError(f"OpenAI HTTP {exc.code}: {body[:500]}") from exc
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as exc:
            raise LLMProviderError(f"OpenAI request failed: {exc}") from exc

        output_text = _extract_output_text(response_data)
        if not output_text:
            raise LLMProviderError("OpenAI response did not contain output text")

        try:
            parsed = json.loads(output_text)
        except json.JSONDecodeError as exc:
            raise LLMProviderError(f"OpenAI response was not valid JSON: {exc}") from exc

        if not isinstance(parsed, dict) or not isinstance(parsed.get("reviews"), list):
            raise LLMProviderError("OpenAI response did not contain a reviews list")
        return parsed


def _extract_output_text(response_data: Dict) -> str:
    if isinstance(response_data.get("output_text"), str):
        return response_data["output_text"]

    chunks = []
    for item in response_data.get("output", []) or []:
        for content in item.get("content", []) or []:
            if isinstance(content.get("text"), str):
                chunks.append(content["text"])
            elif content.get("type") == "output_text" and isinstance(content.get("text"), str):
                chunks.append(content["text"])
    return "".join(chunks).strip()


def _alias_review_schema() -> Dict:
    return {
        "type": "object",
        "additionalProperties": False,
        "properties": {
            "reviews": {
                "type": "array",
                "items": {
                    "type": "object",
                    "additionalProperties": False,
                    "properties": {
                        "review_id": {"type": "string"},
                        "recommendation": {
                            "type": "string",
                            "enum": ["save_alias", "ignore", "skip", "move_to_ambiguous"],
                        },
                        "canonical_name": {"type": ["string", "null"]},
                        "confidence": {
                            "type": "string",
                            "enum": ["low", "medium", "high"],
                        },
                        "reason": {"type": "string"},
                    },
                    "required": [
                        "review_id",
                        "recommendation",
                        "canonical_name",
                        "confidence",
                        "reason",
                    ],
                },
            }
        },
        "required": ["reviews"],
    }

