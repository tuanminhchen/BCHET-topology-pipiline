# Reproducibility notes

1. Raw AlphaFold structures are not redistributed.
2. Raw PDB structures can be downloaded using PDB IDs listed in `data/source_lists/synthetic_pdb_ids.txt`.
3. Processed tables in `data/processed/` are sufficient for the bundled lightweight plotting scripts.
4. Full reprocessing from raw structures requires downloading AlphaFold proteomes (large) and arranging them locally.
5. Synthetic proteins use **B = 1** as a minimal-budget proxy in the main synthetic analysis.
6. Null remapping randomizes B–I assignment while preserving empirical marginal distributions (see `src/null_remapping.py`).
