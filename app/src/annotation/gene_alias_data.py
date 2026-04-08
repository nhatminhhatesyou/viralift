"""
Gene alias mapping (canonical naming system)

All keys must be normalized (lowercase, no spaces, no hyphens).
"""

GENE_ALIAS = {
    # PRRSV
    "orf5": "GP5",
    "gp5": "GP5",
    "glycoprotein5": "GP5",

    "orf4": "GP4",
    "gp4": "GP4",

    "orf3": "GP3",
    "gp3": "GP3",

    # special cases
    "orf5a": "ORF5a",
}