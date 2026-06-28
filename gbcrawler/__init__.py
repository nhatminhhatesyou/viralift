"""gbcrawler — NCBI GenBank crawler that feeds ViraLift.

Crawl sequences from NCBI nuccore (by search query and/or accession list),
download them in GenBank format *with features* (``gbwithparts``), and split
the result into one ``.gb`` file per virus species — ready to drop into
ViraLift (each species + its reference).
"""

__version__ = "0.1.0"
