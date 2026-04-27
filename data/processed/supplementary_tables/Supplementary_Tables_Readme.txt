BCHET Supplementary Tables S1–S4 (auto-generated)

Inputs:
- /Volumes/扩展盘/蛋白质涌现研究/LinearHypo/PlotALL/data/bchet_v7_species_summary.csv
- /Volumes/扩展盘/蛋白质涌现研究/LinearHypo/PlotALL/data/bchet_v7_eta_quantiles_by_species.csv
- /Volumes/扩展盘/蛋白质涌现研究/LinearHypo/PlotALL/data/bchet_v7_species_meta.json
- /Volumes/扩展盘/蛋白质涌现研究/LinearHypo/PlotALL/data/synthetic_proteins_eta.csv
- /Volumes/扩展盘/蛋白质涌现研究/LinearHypo/PlotALL/output/FigS7_null_remapping_stats.csv
- /Volumes/扩展盘/蛋白质涌现研究/LinearHypo/PlotALL/output/FigS7_null_remapping_refined_note.txt
- /Volumes/扩展盘/蛋白质涌现研究/LinearHypo/PlotALL/data/per_species (used for Table S3 if no FigS5 stats CSV exists)

Outputs:
- PlotALL/output/tables/Table_S1_species_eta_summary.csv
- PlotALL/output/tables/Table_S2_synthetic_designed_proteins.csv
- PlotALL/output/tables/Table_S3_alternative_topology_definitions.csv
- PlotALL/output/tables/Table_S4_null_remapping_statistics.csv
- PlotALL/output/tables/BCHET_Supplementary_Tables.xlsx
- Mirrored copies under PlotALL/table/

Table S1 notes:
- `budget_prior` comes from `bchet_v7_species_meta.json` (`intron_over_cds_prior`).
- `B` is left NA at species-summary level because `bchet_v7_species_summary.csv` does not provide a species-level B column in this snapshot.
- `median_eta` is taken from quantiles `p50` after renaming to `q50`.
- `mean_eta` comes from `bchet_v7_species_summary.csv`.
- `n_proteins` maps from summary `n_total`; `n_structures` maps from `n_structures_attempted`.
- `fraction_positive_eta` = `n_positive_eta` / `n_total` when `n_total`>0.

Table S2 notes:
- Directly derived from `synthetic_proteins_eta.csv` with columns preserved/renamed minimally.
- Sorted by `log10_eta` ascending (NaNs last).

Table S3 notes:
No standalone FigS5 stats CSV was found under PlotALL/output; Table S3 was computed directly from `PlotALL/data/per_species/*.csv.gz` using the published Fig S5 alternative-I mapping (same as the plotting script).

Table S4 notes:
- Aggregates `FigS7_null_remapping_stats.csv` per `null_type`.
- Observed values for variance and between-lineage terms are parsed from `FigS7_null_remapping_refined_note.txt`.
- Observed Spearman rho reference is defined as 1.0 (perfect match to observed ordering), per Fig S7 panel B.

Warnings / missing fields:
- `synthetic_class` not present in synthetic CSV; left NA.
