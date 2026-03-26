# ViraLift

Reference-guided viral feature transfer and extraction pipeline using minimap2.

## Input
- reference GenBank file
- query GenBank file

## Current scope
- CDS-only transfer
- one reference, one query
- minimap2-based coordinate transfer

## Build Docker image

```bash
docker build -t viralift .