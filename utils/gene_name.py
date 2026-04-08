def normalize_text(s: str) -> str:
    if not s:
        return ""
    return (
        s.lower()
        .replace("-", "")
        .replace("_", "")
        .replace(" ", "")
        .strip()
    )


ALIASES = {
    "PRRSV": {
        "GP5": ["GP5", "ORF5", "glycoprotein 5", "glycoprotein GP5"],
        "M": ["M", "ORF6", "membrane protein"],
        "N": ["N", "ORF7", "nucleocapsid protein"],
    },
    "FMDV": {
        "VP1": ["VP1", "1D", "capsid protein VP1"],
        "VP2": ["VP2", "1B", "capsid protein VP2"],
        "VP3": ["VP3", "1C", "capsid protein VP3"],
    }
}


def resolve_alias(raw_name: str, virus: str):
    norm = normalize_text(raw_name)

    virus_aliases = ALIASES.get(virus, {})
    for canonical, alias_list in virus_aliases.items():
        for alias in alias_list:
            if normalize_text(alias) == norm:
                return canonical
    return None