#!/usr/bin/env python3
"""固定一对 A/B，扫描多种 w_A 比例，评测 Mix2Phase-QPA 的横向稳定性。"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

TASK0511_ROOT = Path(__file__).resolve().parents[1]

import algo_qpa_mixture_nnls as mixnnls  # noqa: E402
import algo_qpa_q1_cif_xy as q1  # noqa: E402
from mix2phase_common import require_existing_file


def composition_label(atom_types) -> str:
    from collections import Counter

    types = [str(x) for x in atom_types]
    uniq = sorted(set(types))
    c = Counter(types)
    parts = []
    for el in uniq:
        n = c[el]
        parts.append(f"{el}{n}" if n > 1 else el)
    return "".join(parts)


def write_mix_xy(path: Path, y_mix: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        for x, y in zip(q1.GRID, np.asarray(y_mix, dtype=np.float64).reshape(-1)):
            f.write(f"{x:.6f} {y:.8f}\n")


def parse_w_values(spec: str) -> list[float]:
    vals = [float(x.strip()) for x in spec.split(",") if x.strip()]
    if not vals:
        raise ValueError("empty --w-values")
    for w in vals:
        if not (0.0 <= w <= 1.0):
            raise ValueError(f"w must be in [0,1], got {w}")
    return vals


def main() -> None:
    ap = argparse.ArgumentParser(description="固定一对 A/B，在 10 个不同比例下生成混合谱并做 NNLS。")
    ap.add_argument("--ref-a-mp20", type=int, default=2155, help="A 相 mp20 cursor 下标")
    ap.add_argument("--ref-b-mp20", type=int, default=3497, help="B 相 mp20 cursor 下标")
    ap.add_argument("--pair-idx-note", type=int, default=73, help="仅用于文档说明的来源 pair_idx")
    ap.add_argument(
        "--w-values",
        type=str,
        default="0.05,0.15,0.25,0.35,0.45,0.55,0.65,0.75,0.85,0.95",
        help="逗号分隔的 w_A 列表",
    )
    ap.add_argument(
        "--workdir",
        type=str,
        default=str(TASK0511_ROOT / "results/fixed_pair_w10_run"),
        help="输出 mix.xy 的目录",
    )
    ap.add_argument(
        "--json-out",
        type=str,
        default=str(TASK0511_ROOT / "results/fixed_pair_w10_summary.json"),
        help="汇总 JSON 输出路径",
    )
    cli = ap.parse_args()

    w_values = parse_w_values(cli.w_values)
    workdir = Path(cli.workdir)
    json_out = Path(cli.json_out)

    require_existing_file(mixnnls.DEFAULT_MP20_LMDB, purpose="mp20 lookup LMDB")
    ea = mixnnls.load_mp20_entry(mixnnls.DEFAULT_MP20_LMDB, cli.ref_a_mp20)
    eb = mixnnls.load_mp20_entry(mixnnls.DEFAULT_MP20_LMDB, cli.ref_b_mp20)
    y_a = mixnnls.spectrum_from_mp20_entry(ea, True)
    y_b = mixnnls.spectrum_from_mp20_entry(eb, True)

    records = []
    for i, w_true in enumerate(w_values):
        y_mix = w_true * y_a + (1.0 - w_true) * y_b
        tag = f"{i:02d}_w{int(round(w_true * 1000)):03d}"
        mix_xy = workdir / tag / "mix.xy"
        write_mix_xy(mix_xy, y_mix)

        y_m = q1.spectrum_from_mixture_path(mix_xy, normalize_max=True)
        w_hat, coef, resid = q1.nnls_fit_stack(y_m, y_a, y_b)
        records.append(
            {
                "run": i + 1,
                "w_true": float(w_true),
                "w_hat": float(w_hat),
                "abs_error": float(abs(w_hat - w_true)),
                "signed_error": float(w_hat - w_true),
                "residual_l2_nnls": float(resid),
                "nnls_coef_alpha_beta": coef.tolist(),
                "mix_xy_rel": str(mix_xy.relative_to(TASK0511_ROOT)),
            }
        )

    w_true_arr = np.array([r["w_true"] for r in records], dtype=np.float64)
    w_hat_arr = np.array([r["w_hat"] for r in records], dtype=np.float64)
    err_arr = np.array([r["abs_error"] for r in records], dtype=np.float64)

    summary = {
        "script": "scripts/eval_fixed_pair_w10.py",
        "pipeline": "synthetic_mix_from_lookup_refs + same lookup refs + NNLS",
        "pair_idx_note_from_eval50": cli.pair_idx_note,
        "ref_a_mp20": cli.ref_a_mp20,
        "ref_b_mp20": cli.ref_b_mp20,
        "label_A": f"mp20[{cli.ref_a_mp20}] {composition_label(ea['atom_type'])}",
        "label_B": f"mp20[{cli.ref_b_mp20}] {composition_label(eb['atom_type'])}",
        "w_values": w_values,
        "metrics": {
            "mae_abs_w": float(err_arr.mean()),
            "rmse_w": float(np.sqrt(np.mean((w_hat_arr - w_true_arr) ** 2))),
            "max_abs_error": float(err_arr.max()),
        },
        "records": records,
    }

    json_out.parent.mkdir(parents=True, exist_ok=True)
    json_out.write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"[wrote] {json_out}")
    print(
        "[summary] "
        f"A={summary['label_A']} | B={summary['label_B']} | "
        f"MAE={summary['metrics']['mae_abs_w']:.6e} | "
        f"RMSE={summary['metrics']['rmse_w']:.6e} | "
        f"MAX={summary['metrics']['max_abs_error']:.6e}"
    )


if __name__ == "__main__":
    main()
