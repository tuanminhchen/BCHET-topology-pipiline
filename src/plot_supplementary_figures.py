#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Lightweight supplementary figures + optional copied upstream FigS*.png exports.
from __future__ import annotations

import argparse
import shutil
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from utils import repo_root


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=repo_root())
    args = ap.parse_args()

    root: Path = args.root
    out_sup = root / "figures" / "supplementary"
    out_sup.mkdir(parents=True, exist_ok=True)

    s4 = root / "data" / "processed" / "supplementary_tables" / "Table_S4_null_remapping_statistics.csv"
    if s4.exists():
        df = pd.read_csv(s4)
        plt.figure(figsize=(8.2, 4.2), dpi=160)
        sns.barplot(data=df, x="statistic_name", y="empirical_p_value", hue="null_type")
        plt.xticks(rotation=20, ha="right")
        plt.title("Packaged Table S4: empirical p-values (as provided)")
        plt.tight_layout()
        plt.savefig(out_sup / "FigS7_like_tableS4_pvalues.png")
        plt.close()

    upstream = root / "figures" / "_upstream_png"
    upstream.mkdir(parents=True, exist_ok=True)
    if upstream.exists():
        for p in sorted(upstream.glob("FigS*.png")):
            shutil.copy2(p, out_sup / p.name)

    print(f"Wrote figures under {out_sup}")


if __name__ == "__main__":
    main()
