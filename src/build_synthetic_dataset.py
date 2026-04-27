#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
从 RCSB 拉取设计/合成相关 PDB，按 V7 相同规则算 I 与 η（B 为合成蛋白代理常数），
写出 PlotALL/data/synthetic_proteins_eta.csv。
"""

from __future__ import annotations

import argparse
import json
import logging
import math
import os
import sys
from pathlib import Path

_PKG = Path(__file__).resolve().parent
if str(_PKG) not in sys.path:
    sys.path.insert(0, str(_PKG))

import pandas as pd

from rcsb_client import (
    batch_download,
    build_design_query,
    configure_tls,
    fetch_entry_metadata,
    search_all_entry_ids,
)
from v7_topo import compute_I, extract_high_confidence_ca_pdb

_LINEAR_HYPO = _PKG.parent
_DEFAULT_OUT = _LINEAR_HYPO / "PlotALL" / "data" / "synthetic_proteins_eta.csv"
_DEFAULT_WORK = _PKG / "work"


def _safe_title(meta: dict) -> str:
    try:
        t = meta.get("struct", {}).get("title")
        return str(t).replace("\n", " ").strip() if t else ""
    except Exception:
        return ""


def _resolution(meta: dict) -> float | None:
    try:
        v = meta.get("rcsb_entry_info", {}).get("resolution_combined")
        if v is None:
            return None
        return float(v)
    except Exception:
        return None


def _engineered_hint(meta: dict) -> str:
    """尽力从 core entry 推断；无权威字段时留空。"""
    try:
        attrs = meta.get("rcsb_entry_container_identifiers", {})
        if attrs.get("entry_id"):
            pass
        # struct_keywords 列表
        kws = meta.get("struct_keywords", [])
        if isinstance(kws, list):
            flat = " ".join(str(x) for x in kws).lower()
            if "engineer" in flat or "design" in flat or "synthetic" in flat:
                return "yes"
    except Exception:
        pass
    return ""


def compute_rows_for_pdb(
    pdb_path: Path,
    pdb_id: str,
    *,
    B: float,
    B_epsilon: float,
    plddt_min: float,
    meta: dict | None,
) -> dict:
    coords, seq_idx = extract_high_confidence_ca_pdb(pdb_path, plddt_min=plddt_min)
    n_res, I_val = compute_I(coords, seq_idx)
    row: dict = {
        "protein_id": pdb_id.lower(),
        "length": int(n_res),
        "B": float(B),
        "I": float(I_val),
        "source": "synthetic",
    }
    if I_val > 0 and math.isfinite(I_val):
        eta = B / I_val
        row["eta"] = eta
        row["log10_eta"] = math.log10(eta) if eta > 0 else float("nan")
        eta_e = B_epsilon / I_val
        row["B_epsilon"] = float(B_epsilon)
        row["eta_epsilon"] = eta_e
        row["log10_eta_epsilon"] = math.log10(eta_e) if eta_e > 0 else float("nan")
    else:
        row["eta"] = float("nan")
        row["log10_eta"] = float("nan")
        row["B_epsilon"] = float(B_epsilon)
        row["eta_epsilon"] = float("nan")
        row["log10_eta_epsilon"] = float("nan")
    if meta is not None:
        row["title"] = _safe_title(meta)
        r = _resolution(meta)
        row["resolution"] = r if r is not None else float("nan")
        row["engineered_hint"] = _engineered_hint(meta)
    else:
        row["title"] = ""
        row["resolution"] = float("nan")
        row["engineered_hint"] = ""
    return row


def main(argv: list[str] | None = None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    ap = argparse.ArgumentParser(description="Build synthetic protein η=B/I dataset (V7 I).")
    ap.add_argument("--work-dir", type=Path, default=_DEFAULT_WORK)
    ap.add_argument("--out-csv", type=Path, default=_DEFAULT_OUT)
    ap.add_argument("--max-entries", type=int, default=500, help="最多处理的 PDB 条目数")
    ap.add_argument("--B", type=float, default=1.0, help="合成 B 代理（无内含子）主方案")
    ap.add_argument(
        "--B-epsilon",
        type=float,
        default=1e-12,
        help="敏感性分析用极小 B",
    )
    ap.add_argument(
        "--plddt-min",
        type=float,
        default=0.0,
        help="PDB 上默认 0=不过滤 B-factor；与 AF 的 pLDDT≥70 仅在此参数相同时可比",
    )
    ap.add_argument("--max-resolution", type=float, default=3.5)
    ap.add_argument("--skip-search", action="store_true", help="使用已有 id 列表 json")
    ap.add_argument("--id-list-json", type=Path, help="含 pdb id 列表的 JSON（键 entry_ids）")
    ap.add_argument("--skip-download", action="store_true")
    ap.add_argument("--skip-metadata", action="store_true")
    ap.add_argument(
        "--insecure-ssl",
        action="store_true",
        help="关闭 TLS 证书校验（仅当本机仍报 CERTIFICATE_VERIFY_FAILED 时使用）",
    )
    args = ap.parse_args(argv)

    env_insecure = os.environ.get("SYNTH_TLS_INSECURE", "").strip().lower() in (
        "1",
        "true",
        "yes",
    )
    configure_tls(insecure=bool(args.insecure_ssl or env_insecure))

    work = args.work_dir
    pdb_dir = work / "pdb"
    fasta_dir = work / "fasta"
    meta_dir = work / "metadata"
    for d in (pdb_dir, fasta_dir, meta_dir):
        d.mkdir(parents=True, exist_ok=True)

    if args.skip_search and args.id_list_json and args.id_list_json.is_file():
        data = json.loads(args.id_list_json.read_text(encoding="utf-8"))
        entry_ids = [str(x).lower() for x in data.get("entry_ids", [])]
        logging.info("Loaded %d ids from %s", len(entry_ids), args.id_list_json)
    else:
        q = build_design_query(max_resolution=args.max_resolution, protein_only=True)
        entry_ids = search_all_entry_ids(q)
        (work / "entry_ids.json").write_text(
            json.dumps({"entry_ids": entry_ids, "total": len(entry_ids)}, indent=2),
            encoding="utf-8",
        )
        logging.info("Search returned %d entry ids", len(entry_ids))

    entry_ids = entry_ids[: int(args.max_entries)]

    if not args.skip_download:
        ok, bad = batch_download(entry_ids, pdb_dir, fasta_dir, skip_existing=True)
        logging.info("Download ok=%d bad=%d", len(ok), len(bad))

    rows: list[dict] = []
    for pid in entry_ids:
        pl = pid.lower()
        p_pdb = pdb_dir / f"{pl}.pdb"
        p_cif = pdb_dir / f"{pl}.cif"
        if p_pdb.is_file():
            p = p_pdb
        elif p_cif.is_file():
            p = p_cif
        else:
            logging.warning("Missing coordinate file %s.pdb / %s.cif", pl, pl)
            continue
        meta_path = meta_dir / f"{pid.lower()}.json"
        meta: dict | None = None
        if not args.skip_metadata:
            if meta_path.is_file():
                meta = json.loads(meta_path.read_text(encoding="utf-8"))
            else:
                try:
                    meta = fetch_entry_metadata(pid)
                    meta_path.write_text(json.dumps(meta, indent=2), encoding="utf-8")
                except Exception as e:
                    logging.warning("[%s] metadata: %s", pid, e)
                    meta = None
        try:
            r = compute_rows_for_pdb(
                p,
                pid,
                B=float(args.B),
                B_epsilon=float(args.B_epsilon),
                plddt_min=float(args.plddt_min),
                meta=meta,
            )
            rows.append(r)
        except Exception as e:
            logging.warning("[%s] topology failed: %s", pid, e)

    args.out_csv.parent.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(rows)
    if df.empty:
        logging.error("No rows computed; check network / PDB parse errors.")
        return 1
    df.to_csv(args.out_csv, index=False)
    logging.info("Wrote %s (%d rows)", args.out_csv, len(df))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
