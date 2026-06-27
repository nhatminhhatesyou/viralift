# Translation/text accessor for the Streamlit UI.
# Display strings live in ui/ui_text.py (UI_TEXT); this module exposes the
# lookup helper _t so callers don't depend on the dict shape directly.
from ui.ui_text import UI_TEXT


def _t(key: str, **kwargs) -> str:
    template = UI_TEXT["en"].get(key, key)
    return template.format(**kwargs) if kwargs else template
