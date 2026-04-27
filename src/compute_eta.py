#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Compute eta=B/I and log10_eta for a protein-level table.
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--in-csv", type=Path, required=True)
    ap.add_argument("--out-csv", type=Path, required=True)
    args = ap.parse_args()

    df = pd.read_csv(args.in_csv)
    need = {"B", "I"}
    if not need <= set(df.columns):
        raise SystemExit(f"Missing columns: {sorted(need - set(df.columns))}")

    B = pd.to_numeric(df["B"], errors="coerce")
    I = pd.to_numeric(df["I"], errors="coerce")
    eta = B / I
    eta = eta.where(np.isfinite(eta) & (eta > 0))
    df["eta"] = eta
    df["log10_eta"] = np.where(np.isfinite(eta.to_numpy(dtype=float)), np.log10(eta.to_numpy(dtype=float)), np.nan)

    args.out_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out_csv, index=False)
    print(f"Wrote {args.out_csv}")


if __name__ == "__main__":
    main()
