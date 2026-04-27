# Methods summary (minimal)

## Core quantities
- **B**: species-level transcript/genomic budget proxy (V7 uses `intron_over_cds_prior` when per-gene GTF is unavailable).
- **I**: topological load from a Cα contact-span graph (operational proxy; not Shannon information).
- **η**: η = B / I.

## Topology definition (V7-compatible)
- Keep Cα atoms with pLDDT ≥ 70 (AlphaFold B-factor field).
- Contacts: pairwise distance ≤ 8 Å.
- Non-locality: |i − j| > 4 in residue index.
- I = Σ|i − j| over contacts.

## Synthetic designed proteins
Main analyses use a minimal-budget proxy **B = 1** for synthetic proteins (see comments in `src/build_synthetic_dataset.py`).
