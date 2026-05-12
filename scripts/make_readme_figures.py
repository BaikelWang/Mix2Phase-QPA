#!/usr/bin/env python3
"""Generate the two figures embedded in README.md.

Outputs (PNG, 150 dpi):

  docs/figures/pxrd_demo_pair0008.png
      3-panel plot: A reference PXRD, B reference PXRD, mixed PXRD with
      true / predicted A weight and the NNLS-reconstructed mixture.

  docs/figures/accuracy_eval20.png
      2-panel plot summarising the 20-sample evaluation: parity scatter
      (w_true vs w_A_hat) plus per-sample absolute error bar chart.

Run from the repository root:

    python3 scripts/make_readme_figures.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
ROOT_DIR = SCRIPT_DIR.parent
sys.path.insert(0, str(SCRIPT_DIR))

from mix2phase_common import (  # noqa: E402  (sys.path mutation above)
    DEFAULT_MP20_LMDB,
    GRID,
    load_mp20_entry,
    require_existing_file,
)
from algo_qpa_mixture_nnls import spectrum_from_mp20_entry  # noqa: E402
import algo_qpa_q1_cif_xy as q1  # noqa: E402

FIG_DIR = ROOT_DIR / "docs" / "figures"
EVAL_JSON = ROOT_DIR / "results" / "eval20_summary.json"


def _strip_mp20_prefix(label: str) -> str:
    if "]" in label and label.startswith("mp20["):
        return label.split("]", 1)[1].strip()
    return label


def _short_record_label(record: dict) -> str:
    a = _strip_mp20_prefix(record.get("label_A", "A"))
    b = _strip_mp20_prefix(record.get("label_B", "B"))
    return f"{a} / {b}"


def make_pxrd_demo_figure(records: list[dict], out_path: Path) -> dict:
    target_pair = 8
    rec = next((r for r in records if r.get("pair_idx") == target_pair), records[0])
    mix_rel = rec["mix_xy_rel"]
    mix_path = require_existing_file(ROOT_DIR / mix_rel, purpose="example mixed PXRD")
    lmdb_path = require_existing_file(DEFAULT_MP20_LMDB, purpose="mp20 lookup LMDB")

    ea = load_mp20_entry(lmdb_path, rec["src_idx_A"])
    eb = load_mp20_entry(lmdb_path, rec["src_idx_B"])
    y_a = spectrum_from_mp20_entry(ea, normalize_max=True)
    y_b = spectrum_from_mp20_entry(eb, normalize_max=True)
    y_m = q1.spectrum_from_mixture_path(mix_path, normalize_max=True)

    w_true = float(rec["w_true"])
    w_hat = float(rec["w_hat_lookup_mp20"])
    alpha, beta = (float(v) for v in rec["nnls_coef_lookup"])
    y_recon = alpha * y_a + beta * y_b

    label_a = _strip_mp20_prefix(rec.get("label_A", "A"))
    label_b = _strip_mp20_prefix(rec.get("label_B", "B"))

    fig, axes = plt.subplots(3, 1, figsize=(9.5, 7.6), sharex=True)
    color_a = "#1f77b4"
    color_b = "#d62728"
    color_mix = "#2ca02c"
    color_recon = "#7f7f7f"

    axes[0].plot(GRID, y_a, color=color_a, lw=1.0)
    axes[0].set_title(f"Reference PXRD of phase A: {label_a}")
    axes[0].set_ylabel("Intensity (max=1)")
    axes[0].set_ylim(-0.05, 1.1)
    axes[0].grid(alpha=0.25)

    axes[1].plot(GRID, y_b, color=color_b, lw=1.0)
    axes[1].set_title(f"Reference PXRD of phase B: {label_b}")
    axes[1].set_ylabel("Intensity (max=1)")
    axes[1].set_ylim(-0.05, 1.1)
    axes[1].grid(alpha=0.25)

    axes[2].plot(GRID, y_m, color=color_mix, lw=1.1, label="mixed PXRD (input)")
    axes[2].plot(GRID, y_recon, color=color_recon, lw=0.9, ls="--",
                 label=f"NNLS reconstruction = {alpha:.3f}*A + {beta:.3f}*B")
    axes[2].set_title("Mixed PXRD and NNLS reconstruction")
    axes[2].set_xlabel(r"$2\theta$ (deg)")
    axes[2].set_ylabel("Intensity (max=1)")
    axes[2].set_ylim(-0.05, max(1.1, float(y_m.max()) * 1.05))
    axes[2].grid(alpha=0.25)
    axes[2].legend(loc="upper right", fontsize=9)

    annotation = (
        f"true   $w_A$ = {w_true:.4f},  $w_B$ = {1.0 - w_true:.4f}\n"
        f"pred  $\\hat w_A$ = {w_hat:.4f},  $\\hat w_B$ = {1.0 - w_hat:.4f}\n"
        f"abs error = {abs(w_hat - w_true):.2e}"
    )
    axes[2].text(
        0.012, 0.96, annotation,
        transform=axes[2].transAxes,
        ha="left", va="top",
        fontsize=9,
        bbox=dict(boxstyle="round,pad=0.35", facecolor="white", edgecolor="#bbbbbb"),
    )

    fig.suptitle(
        f"Mix2Phase-QPA on pair #{target_pair}: A = {label_a},  B = {label_b}",
        fontsize=12,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    return {
        "pair_idx": target_pair,
        "label_A": label_a,
        "label_B": label_b,
        "w_true": w_true,
        "w_hat": w_hat,
        "abs_error": abs(w_hat - w_true),
    }


def make_accuracy_figure(records: list[dict], out_path: Path) -> dict:
    n = len(records)
    w_true = np.array([float(r["w_true"]) for r in records])
    w_hat = np.array([float(r["w_hat_lookup_mp20"]) for r in records])
    abs_err = np.abs(w_hat - w_true)
    mae = float(abs_err.mean())
    rmse = float(np.sqrt(((w_hat - w_true) ** 2).mean()))
    max_err = float(abs_err.max())

    short_labels = [_short_record_label(r) for r in records]

    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.0), gridspec_kw={"width_ratios": [1, 1.3]})

    ax = axes[0]
    ax.plot([0, 1], [0, 1], color="#888888", lw=1.0, ls="--", label="ideal y = x")
    ax.scatter(w_true, w_hat, s=44, color="#1f77b4", edgecolor="white", zorder=3,
               label=f"20 test samples (MAE = {mae:.2e})")
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    ax.set_xlabel(r"True weight  $w_A$")
    ax.set_ylabel(r"Predicted weight  $\hat w_A$")
    ax.set_title("Parity plot: predicted vs true A weight")
    ax.grid(alpha=0.25)
    ax.legend(loc="upper left", fontsize=9)
    ax.set_aspect("equal", adjustable="box")

    ax2 = axes[1]
    xs = np.arange(n)
    ax2.bar(xs, abs_err * 1000.0, color="#2ca02c", alpha=0.85)
    ax2.axhline(mae * 1000.0, color="#d62728", lw=1.0, ls="--",
                label=f"MAE = {mae*1000:.3f} x 10$^{{-3}}$")
    ax2.axhline(max_err * 1000.0, color="#7f7f7f", lw=1.0, ls=":",
                label=f"max abs err = {max_err*1000:.3f} x 10$^{{-3}}$")
    ax2.set_xticks(xs)
    ax2.set_xticklabels(short_labels, rotation=75, fontsize=7, ha="right")
    ax2.set_ylabel(r"Absolute error  $|\hat w_A - w_A|$  ($\times 10^{-3}$)")
    ax2.set_title("Per-sample absolute error on 20 test pairs")
    ax2.grid(axis="y", alpha=0.25)
    ax2.legend(loc="upper right", fontsize=9)

    fig.suptitle(
        f"Mix2Phase-QPA accuracy on the 20-sample test set  "
        f"(MAE = {mae:.2e},  RMSE = {rmse:.2e},  max abs err = {max_err:.2e})",
        fontsize=12,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    return {"n": n, "mae": mae, "rmse": rmse, "max_abs_err": max_err}


def main() -> None:
    if not EVAL_JSON.is_file():
        raise SystemExit(f"missing eval summary: {EVAL_JSON}")
    summary = json.loads(EVAL_JSON.read_text(encoding="utf-8"))
    records = [r for r in summary.get("records", []) if r.get("status") == "OK"]
    if not records:
        raise SystemExit("no OK records in eval20 summary")

    fig1 = FIG_DIR / "pxrd_demo_pair0008.png"
    fig2 = FIG_DIR / "accuracy_eval20.png"

    info1 = make_pxrd_demo_figure(records, fig1)
    info2 = make_accuracy_figure(records, fig2)

    print(json.dumps({"pxrd_demo": {"path": str(fig1.relative_to(ROOT_DIR)), **info1},
                      "accuracy": {"path": str(fig2.relative_to(ROOT_DIR)), **info2}},
                     indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
