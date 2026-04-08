from src.annotation.gene_alias_data import GENE_ALIAS


def normalize_name(name: str) -> str:
    """Normalize gene name for alias lookup."""
    return (
        name.strip()
        .lower()
        .replace("-", "")
        .replace("_", "")
        .replace(" ", "")
    )


def resolve_gene_name(raw_name: str) -> str:
    """Map raw gene name to canonical name using alias table."""
    if not raw_name:
        return raw_name

    key = normalize_name(raw_name)

    return GENE_ALIAS.get(key, raw_name)