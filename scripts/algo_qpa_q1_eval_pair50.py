#!/usr/bin/env python3
"""Randomly evaluate Mix2Phase-QPA on pair1k_v3 rows using local standalone scripts."""
from __future__ import annotations

import argparse
import gzip
import json
import pickle
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

import algo_qpa_mixture_nnls as mixnnls
import algo_qpa_q1_cif_xy as q1
from mix2phase_common import (
    DEFAULT_MP20_LMDB,
    DEFAULT_PAIR_DIR,
    TASK0511_ROOT,
    composition_label,
    require_existing_dir,
    require_existing_file,
    safe_rel,
    write_mix_xy,
)

_MD_INTERPRETATION = """
## 解读：本评测在测什么？

### 主推：**查表（mp20）+ 混合 `.xy`**

- 与 **[`algo_qpa_mixture_nnls.py`](algo_qpa_mixture_nnls.py)** 一致：按 `src_idx_A/B` 从 **`mp20_data/test.lmdb`** 读 **`pxrd_x/y`** → 插值到统一 GRID → **max 归一** 得 \(y_A,y_B\)；混合谱由 `pxrd_y_mix` 写 `.xy` 再读入做 **max 归一** 得 \(y_{\mathrm{mix}}\)；**NNLS** 得 \(\hat w\)。
- 与 **`w1`（真值）**同源于 mp20 峰表管线，误差量级应接近数值/插值误差，用于验证 **「只解析混合谱、参考谱查表」** 的实现是否正确。

### 可选：`--run-cif`（备选 XRDC）

- 将同一结构写 **CIF** 后用 **pymatgen XRDCalculator** 重算参考谱；该路径与 mp20 存库峰表 **数值不一致** 时，\(\hat w\) 相对 `w1` 往往偏差更大，反映的是 **模拟器不一致**，不是查表主线失效。
""".strip().split("\n")


def load_pair_rows(split: str, pair_dir: Path) -> list:
    import lmdb

    path = pair_dir / f"{split}.lmdb"
    env = lmdb.open(str(path), subdir=False, readonly=True, lock=False)
    with env.begin() as txn:
        n = txn.stat()["entries"]
        rows = [pickle.loads(gzip.decompress(txn.get(str(i).encode()))) for i in range(n)]
    env.close()
    return rows


def entry_to_cif(entry: dict, path: Path) -> None:
    from pymatgen.core import Lattice, Structure

    lattice = np.asarray(entry["lattice_matrix"], dtype=np.float64)
    species = list(entry["atom_type"])
    frac = np.asarray(entry["atom_pos"], dtype=np.float64).reshape(-1, 3)
    structure = Structure(Lattice(lattice), species, frac, coords_are_cartesian=False)
    path.parent.mkdir(parents=True, exist_ok=True)
    structure.to(fmt="cif", filename=str(path))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=50, help="number of samples")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--split", type=str, default="test", choices=("train", "val", "test"))
    ap.add_argument("--pair-dir", type=str, default=str(DEFAULT_PAIR_DIR))
    ap.add_argument("--mp20-lmdb", type=str, default=str(DEFAULT_MP20_LMDB))
    ap.add_argument("--workdir", type=str, default=str(TASK0511_ROOT / "results/q1_eval50_run"))
    ap.add_argument("--md-out", type=str, default=str(TASK0511_ROOT / "results/_autogen/eval50_report.md"))
    ap.add_argument("--json-out", type=str, default=str(TASK0511_ROOT / "results/eval50_summary.json"))
    ap.add_argument("--run-cif", action="store_true", help="also compute CIF->XRDC comparison")
    cli = ap.parse_args()

    pair_dir = require_existing_dir(Path(cli.pair_dir), purpose="pair dataset directory")
    mp20_lmdb = require_existing_file(Path(cli.mp20_lmdb), purpose="mp20 lookup LMDB")
    require_existing_file(pair_dir / f"{cli.split}.lmdb", purpose=f"{cli.split} split LMDB")
    rows = load_pair_rows(cli.split, pair_dir)
    n_pop = len(rows)
    if cli.n > n_pop:
        raise SystemExit(f"split={cli.split} only has {n_pop} rows, cannot sample {cli.n}")

    rng = np.random.default_rng(cli.seed)
    chosen = sorted(rng.choice(n_pop, size=cli.n, replace=False).tolist())
    workdir = Path(cli.workdir)
    records = []

    for run_i, pair_idx in enumerate(chosen):
        row = rows[pair_idx]
        ia = int(row["src_idx_A"])
        ib = int(row["src_idx_B"])
        w_true = float(row["w1"])
        tag = f"{run_i:03d}_pair{pair_idx:04d}"
        d = workdir / tag
        mix_xy = d / "mix.xy"
        cif_a = d / "A.cif"
        cif_b = d / "B.cif"

        try:
            ea = mixnnls.load_mp20_entry(mp20_lmdb, ia)
            eb = mixnnls.load_mp20_entry(mp20_lmdb, ib)
            write_mix_xy(mix_xy, row["pxrd_y_mix"])
            y_a = mixnnls.spectrum_from_mp20_entry(ea, True)
            y_b = mixnnls.spectrum_from_mp20_entry(eb, True)
            y_m = q1.spectrum_from_mixture_path(mix_xy, True)
            w_lookup, coef_l, r_lookup = q1.nnls_fit_stack(y_m, y_a, y_b)
            rec = {
                "run": run_i + 1,
                "pair_idx": pair_idx,
                "status": "OK",
                "src_idx_A": ia,
                "src_idx_B": ib,
                "label_A": f"mp20[{ia}] {composition_label(ea['atom_type'])}",
                "label_B": f"mp20[{ib}] {composition_label(eb['atom_type'])}",
                "mix_xy_rel": safe_rel(mix_xy, TASK0511_ROOT),
                "w_hat_lookup_mp20": float(w_lookup),
                "w_true": w_true,
                "abs_error_lookup": float(abs(w_lookup - w_true)),
                "signed_error_lookup": float(w_lookup - w_true),
                "residual_l2_nnls_lookup": float(r_lookup),
                "nnls_coef_lookup": coef_l.tolist(),
            }
            if cli.run_cif:
                entry_to_cif(ea, cif_a)
                entry_to_cif(eb, cif_b)
                out_cif = q1.run_pipeline(cif_a, cif_b, mix_xy, "CuKa", True, None, emit=False)
                rec["w_hat_cif_xrdc"] = float(out_cif["w_A_hat"])
                rec["abs_error_cif"] = float(abs(rec["w_hat_cif_xrdc"] - w_true))
                rec["signed_error_cif"] = float(rec["w_hat_cif_xrdc"] - w_true)
                rec["residual_l2_nnls_cif"] = float(out_cif["residual_l2_nnls"])
                rec["cif_a_rel"] = safe_rel(cif_a, TASK0511_ROOT)
                rec["cif_b_rel"] = safe_rel(cif_b, TASK0511_ROOT)
            else:
                rec["w_hat_cif_xrdc"] = None
                rec["abs_error_cif"] = None
            records.append(rec)
        except Exception as e:
            records.append(
                {
                    "run": run_i + 1,
                    "pair_idx": pair_idx,
                    "status": "FAIL",
                    "error": str(e),
                    "src_idx_A": ia,
                    "src_idx_B": ib,
                    "w_true": w_true,
                    "run_cif": cli.run_cif,
                }
            )

    ok = [r for r in records if r["status"] == "OK"]
    w_true_arr = np.array([r["w_true"] for r in ok], dtype=np.float64)
    w_lu = np.array([r["w_hat_lookup_mp20"] for r in ok], dtype=np.float64)
    err_lu = np.array([r["abs_error_lookup"] for r in ok], dtype=np.float64)

    summary = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "script": "scripts/algo_qpa_q1_eval_pair50.py",
        "primary_pipeline": "mixture_xy + mp20_lookup",
        "run_cif": cli.run_cif,
        "split": cli.split,
        "seed": cli.seed,
        "requested_n": cli.n,
        "ok_n": len(ok),
        "fail_n": len(records) - len(ok),
        "lookup_mp20": {
            "mae_abs_w": float(err_lu.mean()) if len(ok) else None,
            "rmse_w": float(np.sqrt(np.mean((w_lu - w_true_arr) ** 2))) if len(ok) else None,
            "max_abs_error": float(err_lu.max()) if len(ok) else None,
        },
        "records": records,
    }

    if cli.run_cif:
        w_cif = np.array([r["w_hat_cif_xrdc"] for r in ok], dtype=np.float64)
        err_cif = np.array([r["abs_error_cif"] for r in ok], dtype=np.float64)
        summary["cif_xrdc_optional"] = {
            "mae_abs_w": float(err_cif.mean()) if len(ok) else None,
            "rmse_w": float(np.sqrt(np.mean((w_cif - w_true_arr) ** 2))) if len(ok) else None,
            "max_abs_error": float(err_cif.max()) if len(ok) else None,
        }

    json_out = Path(cli.json_out)
    json_out.parent.mkdir(parents=True, exist_ok=True)
    json_out.write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")

    lines = [
        "# Mix2Phase-QPA：随机 N 条 pair 评测",
        "",
        f"- **生成时间（UTC）**: {summary['generated_utc']}",
        f"- **数据**: `{safe_rel(pair_dir / (cli.split + '.lmdb'), TASK0511_ROOT)}`（总体 **{n_pop}** 条，无放回随机 **{cli.n}** 条，`seed={cli.seed}`）",
        f"- **主推流程**: `pxrd_y_mix` → **`mix.xy`**；`src_idx_A/B` → **mp20 `test.lmdb`** 查 **`pxrd_x/y`** → 统一插值与 **max 归一** → **NNLS**。",
        f"- **真值**: LMDB **`w1`**。",
        f"- **CIF/XRDC 对照**: {'已启用 `--run-cif`' if cli.run_cif else '**未启用**（默认）'}。",
        f"- **产物目录**: `{safe_rel(workdir, TASK0511_ROOT)}/<序号>_pairXXXX/mix.xy`" + ("；若 `--run-cif` 则同目录含 `A.cif`、`B.cif`。" if cli.run_cif else "。"),
        f"- **JSON**: `{safe_rel(json_out, TASK0511_ROOT)}`",
        "",
        "## 汇总指标（成功样本）",
        "",
        "| 成功 | 失败 | 指标 | MAE(|ŵ−w|) | RMSE | max |ŵ−w| |",
        "| --- | --- | --- | --- | --- | --- |",
        f"| {len(ok)} | {summary['fail_n']} | **查表（mp20）** | {summary['lookup_mp20']['mae_abs_w']:.6f} | {summary['lookup_mp20']['rmse_w']:.6f} | {summary['lookup_mp20']['max_abs_error']:.6f} |",
        "",
    ]
    if cli.run_cif:
        lines.append(
            f"| {len(ok)} | {summary['fail_n']} | **备选 CIF→XRDC** | {summary['cif_xrdc_optional']['mae_abs_w']:.6f} | {summary['cif_xrdc_optional']['rmse_w']:.6f} | {summary['cif_xrdc_optional']['max_abs_error']:.6f} |"
        )
        lines.append("")

    lines.extend(_MD_INTERPRETATION)
    lines.extend(["", "## 逐条明细", ""])
    if cli.run_cif:
        lines.extend(
            [
                "| # | pair | A | B | 混合 PXRD | ŵ 查表 | ŵ CIF | w_true | |e|查表 | |e|CIF |",
                "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |",
            ]
        )
    else:
        lines.extend(
            [
                "| # | pair | A | B | 混合 PXRD | ŵ 查表 | w_true | |e|查表 |",
                "| --- | --- | --- | --- | --- | --- | --- | --- |",
            ]
        )

    for r in records:
        if r["status"] != "OK":
            msg = r.get("error", "")
            esc = (msg[:80] + "…") if len(msg) > 80 else msg
            if cli.run_cif:
                lines.append(f"| {r['run']} | {r['pair_idx']} | — | — | — | — | — | {r['w_true']:.6f} | FAIL | `{esc}` |")
            else:
                lines.append(f"| {r['run']} | {r['pair_idx']} | — | — | — | — | {r['w_true']:.6f} | FAIL `{esc}` |")
            continue
        if cli.run_cif:
            lines.append(
                f"| {r['run']} | {r['pair_idx']} | {r['label_A']} | {r['label_B']} | `{r['mix_xy_rel']}` | "
                f"{r['w_hat_lookup_mp20']:.6f} | {r['w_hat_cif_xrdc']:.6f} | {r['w_true']:.6f} | {r['abs_error_lookup']:.6f} | {r['abs_error_cif']:.6f} |"
            )
        else:
            lines.append(
                f"| {r['run']} | {r['pair_idx']} | {r['label_A']} | {r['label_B']} | `{r['mix_xy_rel']}` | "
                f"{r['w_hat_lookup_mp20']:.6f} | {r['w_true']:.6f} | {r['abs_error_lookup']:.6f} |"
            )

    lines.extend(["", "## 说明", "", "- **混合 PXRD**：由 `pxrd_y_mix` 写成的 `.xy`，读入规则与推理脚本一致。", "- **查表**：mp20 索引即 pair 中 `src_idx_A` / `src_idx_B`，无需 CIF。", ""])
    md_path = Path(cli.md_out)
    md_path.parent.mkdir(parents=True, exist_ok=True)
    md_path.write_text("\n".join(lines), encoding="utf-8")

    print(f"[wrote] {md_path}", flush=True)
    print(f"[wrote] {json_out}", flush=True)
    extra = ""
    if cli.run_cif:
        extra = f" MAE_CIF={summary['cif_xrdc_optional']['mae_abs_w']}"
    print(f"[summary] ok={len(ok)} fail={summary['fail_n']} MAE_lookup={summary['lookup_mp20']['mae_abs_w']}{extra}", flush=True)


if __name__ == "__main__":
    main()
