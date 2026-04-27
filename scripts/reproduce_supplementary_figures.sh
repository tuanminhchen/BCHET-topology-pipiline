#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

export MPLBACKEND=Agg

mkdir -p figures/supplementary figures/_upstream_png
python3 src/plot_supplementary_figures.py --root "$ROOT"
