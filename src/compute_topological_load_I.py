#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Compute BCHET topological load I (V7-compatible) from mmCIF.gz or a CA coordinate CSV.
from __future__ import annotations

import argparse
import gzip
from pathlib import Path

import numpy as np
import pandas as pd
from Bio.PDB import FastMMCIFParser
from scipy.spatial.distance import pdist, squareform

PLDDT_THRESHOLD = 70.0
DISTANCE_CUTOFF = 8.0
LOCAL_WINDOW = 4


def extract_high_confidence_ca(cif_path: Path, plddt_min: float) -> tuple[np.ndarray, np.ndarray]:
    parser = FastMMCIFParser(QUIET=True)
    with gzip.open(cif_path, "rt", encoding="utf-8", errors="ignore") as fh:
        structure = parser.get_structure("af", fh)
    model = next(structure.get_models())
    coords = []
    seq_idx = []
    for chain in model:
        for residue in chain:
            if residue.id[0] != " ":
                continue
            if "CA" not in residue:
                continue
            ca_atom = residue["CA"]
            bf = getattr(ca_atom, "bfactor", None)
            if bf is None or float(bf) < plddt_min:
                continue
            coords.append(ca_atom.coord)
            seq_idx.append(int(residue.id[1]))
    if len(coords) < 2:
        raise ValueError(f"Not enough high-confidence CA in {cif_path.name}")
    return np.asarray(coords, dtype=float), np.asarray(seq_idx, dtype=int)


def compute_I(coords: np.ndarray, seq_idx: np.ndarray) -> tuple[int, float]:
    n_res = len(coords)
    if n_res < 2:
        return n_res, float("nan")
    dist_matrix = squareform(pdist(coords, metric="euclidean"))
    iu, ju = np.triu_indices(n_res, k=1)
    spatial = dist_matrix[iu, ju] <= DISTANCE_CUTOFF
    non_local = np.abs(seq_idx[iu] - seq_idx[ju]) > LOCAL_WINDOW
    mask = spatial & non_local
    if not np.any(mask):
        return n_res, 0.0
    return n_res, float(np.sum(np.abs(seq_idx[iu[mask]] - seq_idx[ju[mask]])))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cif-gz", type=Path, default=None)
    ap.add_argument("--coords-csv", type=Path, default=None)
    ap.add_argument("--protein-id", type=str, default="unknown")
    ap.add_argument("--plddt-min", type=float, default=PLDDT_THRESHOLD)
    ap.add_argument("--out-csv", type=Path, required=True)
    args = ap.parse_args()

    if bool(args.cif_gz) == bool(args.coords_csv):
        raise SystemExit("Provide exactly one of --cif-gz or --coords-csv")

    if args.cif_gz is not None:
        coords, seq_idx = extract_high_confidence_ca(args.cif_gz, float(args.plddt_min))
        pid = args.protein_id
    else:
        df = pd.read_csv(args.coords_csv)
        need = {"residue_index", "x", "y", "z"}
        if not need <= set(df.columns):
            raise SystemExit(f"coords csv missing columns: {sorted(need - set(df.columns))}")
        seq_idx = pd.to_numeric(df["residue_index"], errors="coerce").astype(int).to_numpy()
        coords = df[["x", "y", "z"]].to_numpy(dtype=float)
        pid = args.protein_id

    L_aa, I_topo = compute_I(coords, seq_idx)
    out = pd.DataFrame([{"protein_id": pid, "I": float(I_topo), "length": int(L_aa)}])
    args.out_csv.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out_csv, index=False)
    print(f"Wrote {args.out_csv} ({pid} L={L_aa} I={I_topo})")


if __name__ == "__main__":
    main()
