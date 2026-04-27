#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Minimal null remapping controls (natural eukaryotes) mirroring Fig S7 intent.
#
# Inputs: long per-protein CSV with columns including:
#   species_key, domain, I_topo (or I), L_aa, B_intron_budget (or B), optional intron_over_cds_prior
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def between_lineage_var(med: pd.Series, lineage: pd.Series) -> float:
    x = med.groupby(lineage).mean().dropna()
    x = x[np.isfinite(x.to_numpy(dtype=float))]
    if len(x) < 2:
        return float("nan")
    return float(np.var(x.to_numpy(dtype=float), ddof=1))


def lineage_from_realm(realm: str) -> str:
    r = str(realm or "").lower()
    if "fung" in r:
        return "fungi"
    if "metaz" in r or "animal" in r:
        return "metazoa"
    if "plant" in r or "virid" in r:
        return "plants"
    if any(k in r for k in ("protist", "protozo", "amoeb", "alga", "chrom", "stramen", "apicom", "alveol")):
        return "protists"
    return "other"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--in-csv", type=Path, required=True)
    ap.add_argument("--meta-json", type=Path, default=None)
    ap.add_argument("--n-perm", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--out-stats-csv", type=Path, required=True)
    args = ap.parse_args()

    df = pd.read_csv(args.in_csv, low_memory=False)
    if "domain" not in df.columns:
        raise SystemExit("Missing domain column")

    df = df[df["domain"].astype(str).str.lower() == "eukaryota"].copy()
    if df.empty:
        raise SystemExit("No eukaryotic rows after filtering")

    Icol = "I_topo" if "I_topo" in df.columns else ("I" if "I" in df.columns else None)
    if Icol is None:
        raise SystemExit("Need I_topo or I")

    df[Icol] = pd.to_numeric(df[Icol], errors="coerce")
    if "L_aa" in df.columns:
        df["L_aa"] = pd.to_numeric(df["L_aa"], errors="coerce")
    elif "length" in df.columns:
        df["L_aa"] = pd.to_numeric(df["length"], errors="coerce")
    else:
        df["L_aa"] = np.nan
    df = df[np.isfinite(df[Icol]) & (df[Icol] > 0) & np.isfinite(df["L_aa"]) & (df["L_aa"] > 0)].copy()

    if "B_intron_budget" in df.columns and df["B_intron_budget"].notna().any():
        B = pd.to_numeric(df["B_intron_budget"], errors="coerce")
    elif "B" in df.columns and pd.to_numeric(df["B"], errors="coerce").notna().any():
        # Sample / simplified tables may ship a per-row species-level B proxy already.
        B = pd.to_numeric(df["B"], errors="coerce")
    else:
        prior = (
            pd.to_numeric(df["intron_over_cds_prior"], errors="coerce")
            if "intron_over_cds_prior" in df.columns
            else pd.Series(np.nan, index=df.index)
        )
        if bool(prior.notna().any()) is False:
            if args.meta_json is None:
                raise SystemExit("Need B (or B_intron_budget), intron_over_cds_prior column, or --meta-json")
            meta = json.loads(Path(args.meta_json).read_text(encoding="utf-8"))
            pmap = {k: float((v or {}).get("intron_over_cds_prior", np.nan)) for k, v in (meta.get("species") or {}).items()}
            prior = df["species_key"].astype(str).map(pmap)
        B = prior * (3.0 * df["L_aa"])

    eta = B / df[Icol]
    log10_eta = np.log10(eta.where(np.isfinite(eta) & (eta > 0)))
    obs_med = log10_eta.groupby(df["species_key"]).median().dropna()
    obs_var = float(np.var(obs_med.to_numpy(), ddof=1)) if len(obs_med) >= 2 else float("nan")

    realm = df.get("realm", pd.Series([""] * len(df))).astype(str)
    lineage = realm.groupby(df["species_key"]).agg(lambda s: s.dropna().iloc[0] if len(s.dropna()) else "").map(lineage_from_realm)
    obs_sep = between_lineage_var(obs_med, lineage.reindex(obs_med.index))

    rng = np.random.default_rng(int(args.seed))
    n = int(args.n_perm)

    # Null 1
    if "intron_over_cds_prior" in df.columns and pd.to_numeric(df["intron_over_cds_prior"], errors="coerce").notna().any():
        prior_row = pd.to_numeric(df["intron_over_cds_prior"], errors="coerce")
    else:
        if args.meta_json is None:
            raise SystemExit("Null 1 needs intron_over_cds_prior per row or --meta-json")
        meta = json.loads(Path(args.meta_json).read_text(encoding="utf-8"))
        pmap = {k: float((v or {}).get("intron_over_cds_prior", np.nan)) for k, v in (meta.get("species") or {}).items()}
        prior_row = df["species_key"].astype(str).map(pmap)

    base_term = (np.log10(3.0 * df["L_aa"]) - np.log10(df[Icol])).replace([np.inf, -np.inf], np.nan)
    base_term_med = base_term.groupby(df["species_key"]).median()

    species = sorted(set(obs_med.index.astype(str)) & set(base_term_med.index.astype(str)))
    pri_by_sp = prior_row.groupby(df["species_key"]).median().reindex(species).astype(float).to_numpy()
    base_by_sp = base_term_med.reindex(species).astype(float).to_numpy()
    obs_rank = obs_med.reindex(species).rank(ascending=True, method="average")

    var1 = np.empty(n)
    rho1 = np.empty(n)
    sep1 = np.empty(n)
    med_buf = np.empty(len(species), dtype=float)
    for t in range(n):
        pp = rng.permutation(pri_by_sp)
        med_buf[:] = np.log10(np.clip(pp, 1e-20, None)) + base_by_sp
        var1[t] = float(np.var(med_buf, ddof=1)) if len(med_buf) >= 2 else np.nan
        rho1[t] = float(
            obs_rank.corr(pd.Series(med_buf, index=species).rank(method="average"), method="pearson")
        )
        sep1[t] = between_lineage_var(pd.Series(med_buf, index=species), lineage.reindex(species))

    # Null 2
    Bv = pd.to_numeric(B, errors="coerce").to_numpy(dtype=float)
    Iv = df[Icol].to_numpy(dtype=float)
    sk = df["species_key"].astype(str).to_numpy()
    order = np.argsort(sk, kind="mergesort")
    sks = sk[order]
    Bs = Bv[order]
    Is = Iv[order]
    good = np.isfinite(Bs) & (Bs > 0) & np.isfinite(Is) & (Is > 0)
    sks, Bs, Is = sks[good], Bs[good], Is[good]
    _, start_idx, counts = np.unique(sks, return_index=True, return_counts=True)
    ends = start_idx + counts
    species2 = pd.Index(np.unique(sks))
    logB = np.log10(Bs)
    logI = np.log10(Is)
    obs2 = np.log10(Bs / Is)
    obs_med2 = pd.Series(obs2).groupby(sks).median().reindex(species2).dropna()
    obs_rank2 = obs_med2.rank(ascending=True, method="average")

    perm_idx = np.arange(len(Bs), dtype=int)
    diff = np.empty(len(Bs), dtype=float)
    med_buf2 = np.empty(len(species2), dtype=float)
    var2 = np.empty(n)
    rho2 = np.empty(n)
    sep2 = np.empty(n)
    lin2 = lineage.reindex(species2).fillna("other")
    for t in range(n):
        rng.shuffle(perm_idx)
        diff[:] = logB[perm_idx] - logI
        for i, (a, b) in enumerate(zip(start_idx, ends)):
            med_buf2[i] = float(np.median(diff[a:b]))
        med_s = pd.Series(med_buf2, index=species2)
        var2[t] = float(np.var(med_buf2, ddof=1)) if len(med_buf2) >= 2 else np.nan
        rho2[t] = float(obs_rank2.corr(med_s.rank(method="average"), method="pearson"))
        sep2[t] = between_lineage_var(med_s, lin2)

    rows_out = []
    for i in range(n):
        rows_out.append(
            {
                "null_type": "Null1_species_label_remap",
                "permutation_id": i,
                "var_species_medians": var1[i],
                "spearman_rho_to_observed": rho1[i],
                "between_lineage_variance": sep1[i],
            }
        )
        rows_out.append(
            {
                "null_type": "Null2_protein_B_shuffle",
                "permutation_id": i,
                "var_species_medians": var2[i],
                "spearman_rho_to_observed": rho2[i],
                "between_lineage_variance": sep2[i],
            }
        )

    out = pd.DataFrame(rows_out)
    args.out_stats_csv.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out_stats_csv, index=False)

    note_path = args.out_stats_csv.with_suffix(".note.txt")
    note_path.write_text(
        "\n".join(
            [
                "observed_var_species_medians={v}".format(v=obs_var),
                "observed_between_lineage_variance={s}".format(s=obs_sep),
                "observed_spearman_reference=1.0",
                "N_perm={n}".format(n=n),
                "seed={sd}".format(sd=int(args.seed)),
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    print(f"Wrote {args.out_stats_csv}")


if __name__ == "__main__":
    main()
