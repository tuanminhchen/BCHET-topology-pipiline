# BCHET Topology Pipeline

## 1. Overview
This repository contains the code and processed data used for the manuscript:
**"A Budget–Topology Structure Links Natural and Designed Proteins Across Species"**.

## 2. Main quantities
- **I**: species-level intron-derived budget proxy (V7-style: `intron_over_cds_prior` when per-gene GTF is unavailable).
- **T**: topological load from a Cα contact-span graph (**operational proxy**, not Shannon information).
- **η**: η = I / T.

Clarifications:
- **T** is an operational topological load proxy derived from AlphaFold coordinates under explicit contact rules.
- **I** is a **species-level prior**, not a per-gene causal measurement.

## 3. Data sources
- AlphaFold DB: https://alphafold.ebi.ac.uk/
- RCSB PDB: https://www.rcsb.org/
- Ensembl: https://www.ensembl.org/
- UniProt: https://www.uniprot.org/

## 4. Repository structure
- `data/source_lists/`: species lists, AlphaFold tar basenames (from local extraction metadata), synthetic PDB IDs.
- `data/processed/`: processed summary tables (+ packaged supplementary tables if available).
- `data/sample/`: small real subsamples for quick tests.
- `src/`: minimal computation + plotting entrypoints.
- `scripts/`: reproducibility shell wrappers.
- `docs/`: sources + methods + reproducibility notes.
- `figures/`: outputs from bundled plotting scripts (and optional copied upstream PNGs).

## 5. Installation
- Python >= 3.10
- `pip install -r requirements.txt`

## 6. Quick test
- `bash scripts/quick_test.sh`

## 7. Reproduce main figures
- `bash scripts/reproduce_main_figures.sh`

## 8. Reproduce supplementary figures
- `bash scripts/reproduce_supplementary_figures.sh`

## 9. Notes on raw data
Raw AlphaFold proteome tarballs and PDB/mmCIF archives are large and publicly available; therefore they are not redistributed here.

## 10. Citation
If you use this code, please cite the associated preprint (**TODO**: add DOI/biorxiv link when available).

## Known gaps (explicit)
- `species_list_46.csv:alphafold_url` is **NA** unless/until a verified canonical URL scheme is pinned for each proteome tarball.
- `species_budget_prior.csv:B` is **NA** because the packaged species summary snapshot does not provide a single scalar species-level B without extra assumptions.
