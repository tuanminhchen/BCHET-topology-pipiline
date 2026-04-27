#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

export MPLBACKEND=Agg

mkdir -p figures/main figures/_upstream_png
python3 src/plot_main_figures.py --root "$ROOT"
