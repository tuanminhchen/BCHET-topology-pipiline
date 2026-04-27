# BCHET Topology Pipeline

## Overview
This repository contains code, processed tables, and figure scripts for the BCHET framework:

**Budget-Constraint Hypothesis of Emergent Topology (BCHET)**  
Exons encode protein material, while intron architecture provides a genome-level budget layer associated with non-local folding topology.

The current implementation is designed for falsifiable, cross-system tests rather than narrative-only interpretation.

## Core Quantities (Current Notation)
- **T_p**: protein topological load from non-local C-alpha contact spans (operational structural proxy).
- **I_{s,p}**: intron-derived budget proxy for protein `p` in species `s`.
- **eta_{s,p} = I_{s,p} / T_p**.

In this repository:
- `I_{s,p} = r_s * 3 * L_aa,p`
- `r_s` is the species-level intron-over-CDS prior.

### Important Clarifications
- `T_p` is **not** Shannon information and **not** a full thermodynamic model; it is an operational non-local closure load proxy.
- `I_{s,p}` is a **species-level prior-based proxy**, not a direct per-gene intron count in the current release.
- Results should be interpreted as evidence for a structured budget-topology association, with controls against simple algebraic or length-only explanations.

## Contact and Topology Rules
The BCHET topology pipeline uses:
- C-alpha residue coordinates (AlphaFold/PDB structures)
- confidence filter: `pLDDT >= 70` (for AlphaFold)
- contact: `d_ij <= 8 Angstrom`
- non-local: `|seq_idx_i - seq_idx_j| > 4`
- each unordered residue pair counted once

Then:
- `T_p = sum_{(i,j) in C_p} |seq_idx_i - seq_idx_j|`

## Data Sources
- AlphaFold DB: https://alphafold.ebi.ac.uk/
- RCSB PDB: https://www.rcsb.org/
- Ensembl: https://www.ensembl.org/
- UniProt: https://www.uniprot.org/

## Repository Layout
- `data/source_lists/`: species lists, AlphaFold source records, designed-protein PDB IDs
- `data/processed/`: processed summary tables and supplementary tables
- `data/sample/`: small real-data subsets for quick reproducibility checks
- `src/`: computation and plotting entrypoints
- `scripts/`: reproduction wrappers
- `docs/`: methods, provenance, reproducibility notes
- `figures/`: generated figures and packaged outputs

## Installation
- Python >= 3.10
- Install dependencies:
  - `pip install -r requirements.txt`

## Reproducibility Commands
- Quick test:
  - `bash scripts/quick_test.sh`
- Main figures:
  - `bash scripts/reproduce_main_figures.sh`
- Supplementary figures:
  - `bash scripts/reproduce_supplementary_figures.sh`

## What BCHET Tests
This release is organized around falsifiable checks:
- cross-lineage regime structure in `(I, T, eta)` space
- designed-vs-natural contrast
- scaling and length controls
- bootstrap stability
- null remapping controls (preserving `eta = I/T` form while breaking biological assignment)
- structure-level contact-map manifestation

## Scope and Limitations
- Raw AlphaFold proteome archives and full PDB/mmCIF collections are not redistributed here due to size and upstream availability.
- The present `I_{s,p}` uses species-level priors; gene-resolved intron-budget estimation is a future extension.
- Designed proteins are projected with a minimal intron-budget proxy for contrastive analysis, not biological intron attribution.

## Citation
If you use this code or processed tables, please cite the associated manuscript/preprint.
DOI/preprint link will be added once finalized.

## Known Gaps (Explicit)
- `species_list_46.csv:alphafold_url` may remain `NA` until a canonical per-proteome URL mapping is pinned.
- Some budget-prior fields are proxy-level summaries and should not be interpreted as gene-level ground truth.
