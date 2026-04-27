#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

export MPLBACKEND=Agg

mkdir -p output

python3 src/compute_eta.py               --in-csv data/sample/sample_natural_proteins.csv               --out-csv output/sample_eta.csv

python3 src/null_remapping.py               --in-csv data/sample/sample_natural_proteins.csv               --meta-json data/bchet_v7_species_meta.json               --n-perm 100               --seed 42               --out-stats-csv output/sample_null_stats.csv
