#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Lightweight main figures from processed summary tables + optional copied upstream PNGs.
from __future__ import annotations

import argparse
import shutil
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from utils import repo_root


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=repo_root())
    args = ap.parse_args()

    root: Path = args.root
    proc = root / "data" / "processed"
    summ = pd.read_csv(proc / "bchet_species_summary.csv")
    qu = pd.read_csv(proc / "bchet_eta_quantiles_by_species.csv")
    meta = pd.read_csv(root / "data" / "source_lists" / "species_list_46.csv")

    m = summ.merge(qu, on="species_key", how="left", suffixes=("", "_qu"))
    m = m.merge(meta[["species_key", "domain", "lineage_group", "phylo_group"]], on="species_key", how="left")

    p50 = pd.to_numeric(m.get("p50"), errors="coerce")
    m["log10_median_eta"] = np.where(np.isfinite(p50) & (p50 > 0), np.log10(p50), np.nan)

    out_main = root / "figures" / "main"
    out_main.mkdir(parents=True, exist_ok=True)

    plt.figure(figsize=(7.2, 4.2), dpi=160)
    sns.scatterplot(data=m, x="phylo_group", y="log10_median_eta", hue="domain", s=35, alpha=0.85)
    plt.title("Species median log10(eta) vs phylo_group (processed tables)")
    plt.tight_layout()
    plt.savefig(out_main / "Fig2_like_species_median_eta_vs_phylo.png")
    plt.close()

    upstream = root / "figures" / "_upstream_png"
    if upstream.exists():
        for p in sorted(upstream.glob("*.png")):
            shutil.copy2(p, out_main / p.name)

    print(f"Wrote figures under {out_main}")


if __name__ == "__main__":
    main()
