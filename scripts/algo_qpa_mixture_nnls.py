#!/usr/bin/env python3
"""Main Mix2Phase-QPA path: mixed PXRD + lookup/file references -> NNLS."""
from __future__ import annotations

import argparse
import json
import sys
import tempfile
from pathlib import Path

import numpy as np

import algo_qpa_q1_cif_xy as q1
from mix2phase_common import DEFAULT_MP20_LMDB, GRID, STEP, X_MAX, interpolate_pxrd, load_mp20_entry, require_existing_file


def spectrum_from_mp20_entry(entry: dict, normalize_max: bool) -> np.ndarray:
    _, py = interpolate_pxrd(
        X_MAX,
        np.asarray(entry["pxrd_x"], dtype=np.float64),
        np.asarray(entry["pxrd_y"], dtype=np.float64),
        num_fill=2,
        threshold=5,
        step=STEP,
    )
    y = np.asarray(py, dtype=np.float64)
    if normalize_max and y.max() > 0:
        y = y / y.max()
    return y


def spectrum_from_ref_path(path: Path, normalize_max: bool) -> np.ndarray:
    return q1.spectrum_from_mixture_path(path, normalize_max)


def run_mixture_nnls(
    mix_path: Path,
    y_a: np.ndarray,
    y_b: np.ndarray,
    *,
    normalize_max: bool,
    json_out: str | None,
    meta: dict,
    emit: bool = True,
) -> dict:
    y_m = q1.spectrum_from_mixture_path(mix_path, normalize_max)
    w_hat, coef, resid_nnls = q1.nnls_fit_stack(y_m, y_a, y_b)
    out = {
        "w_A_hat": w_hat,
        "w_B_hat": 1.0 - w_hat,
        "nnls_coef_alpha_beta": coef.tolist(),
        "residual_l2_nnls": resid_nnls,
        "residual_l2_convex_combo": q1.residual_l2_convex(y_m, y_a, y_b, w_hat),
        "grid": {"x_max": X_MAX, "step": STEP, "n": int(len(GRID))},
        "normalize_max": normalize_max,
        "pipeline": "mixture_nnls_lookup_refs",
        **meta,
    }
    if emit:
        print(json.dumps(out, indent=2, ensure_ascii=False))
    if json_out:
        Path(json_out).parent.mkdir(parents=True, exist_ok=True)
        Path(json_out).write_text(json.dumps(out, indent=2, ensure_ascii=False), encoding="utf-8")
        print(f"[wrote] {json_out}", file=sys.stderr)
    return out


def run_demo(normalize_max: bool, json_out: str | None) -> None:
    require_existing_file(DEFAULT_MP20_LMDB, purpose="mp20 lookup LMDB")
    e0 = load_mp20_entry(DEFAULT_MP20_LMDB, 0)
    e1 = load_mp20_entry(DEFAULT_MP20_LMDB, 1)
    y_a = spectrum_from_mp20_entry(e0, normalize_max)
    y_b = spectrum_from_mp20_entry(e1, normalize_max)
    w_true = 0.37
    y_m = w_true * y_a + (1.0 - w_true) * y_b + np.random.default_rng(0).normal(0, 0.005, y_a.shape)
    with tempfile.NamedTemporaryFile("w", suffix="_mix2phase_demo.xy", delete=False, encoding="utf-8") as tmp:
        mix_path = Path(tmp.name)
        for x, y in zip(GRID, y_m):
            tmp.write(f"{x:.6f} {y:.8f}\n")
    meta = {
        "refs": {"mode": "mp20_lmdb", "lmdb": str(DEFAULT_MP20_LMDB), "idx_A": 0, "idx_B": 1},
        "inputs": {"mixture": str(mix_path)},
        "demo_w_true": w_true,
    }
    run_mixture_nnls(
        mix_path,
        y_a,
        y_b,
        normalize_max=normalize_max,
        json_out=json_out,
        meta=meta,
        emit=True,
    )


def main() -> None:
    ap = argparse.ArgumentParser(description="Mix2Phase-QPA main path: mixed PXRD + lookup/file refs -> NNLS")
    ap.add_argument("--mixture", type=str, default=None, help="mixed pattern path (.xy/.csv/.npy/.npz)")
    ap.add_argument("--ref-a-mp20", type=int, default=None, help="A phase mp20 cursor index")
    ap.add_argument("--ref-b-mp20", type=int, default=None, help="B phase mp20 cursor index")
    ap.add_argument("--mp20-lmdb", type=str, default=str(DEFAULT_MP20_LMDB), help="mp20 lmdb with pxrd_x/y")
    ap.add_argument("--ref-a", type=str, default=None, help="A reference pattern file")
    ap.add_argument("--ref-b", type=str, default=None, help="B reference pattern file")
    ap.add_argument("--no-normalize", action="store_true", help="skip per-spectrum max=1")
    ap.add_argument("--json-out", type=str, default=None)
    ap.add_argument("--demo", action="store_true", help="synthesize a demo mixture from two mp20 entries")
    cli = ap.parse_args()
    normalize_max = not cli.no_normalize

    if cli.demo:
        run_demo(normalize_max, cli.json_out)
        return

    if not cli.mixture:
        ap.error("need --mixture (or use --demo)")

    has_mp20 = cli.ref_a_mp20 is not None or cli.ref_b_mp20 is not None
    has_files = cli.ref_a is not None or cli.ref_b is not None
    if has_mp20 and has_files:
        ap.error("do not mix mp20 lookup refs with file refs")
    if has_mp20:
        if cli.ref_a_mp20 is None or cli.ref_b_mp20 is None:
            ap.error("mp20 mode requires --ref-a-mp20 and --ref-b-mp20 together")
    elif has_files:
        if cli.ref_a is None or cli.ref_b is None:
            ap.error("file-ref mode requires --ref-a and --ref-b together")
    else:
        ap.error("please provide (--ref-a-mp20/--ref-b-mp20) or (--ref-a/--ref-b)")

    mix_path = Path(cli.mixture)
    if not mix_path.is_file():
        raise SystemExit(f"missing mixture file: {mix_path}")

    use_mp20 = has_mp20
    if use_mp20:
        lmdb_path = require_existing_file(Path(cli.mp20_lmdb), purpose="mp20 lookup LMDB")
        ea = load_mp20_entry(lmdb_path, cli.ref_a_mp20)
        eb = load_mp20_entry(lmdb_path, cli.ref_b_mp20)
        y_a = spectrum_from_mp20_entry(ea, normalize_max)
        y_b = spectrum_from_mp20_entry(eb, normalize_max)
        meta = {
            "refs": {
                "mode": "mp20_lmdb",
                "lmdb": str(lmdb_path),
                "idx_A": cli.ref_a_mp20,
                "idx_B": cli.ref_b_mp20,
            },
            "inputs": {"mixture": str(mix_path)},
        }
    else:
        ref_a = Path(cli.ref_a)
        ref_b = Path(cli.ref_b)
        if not ref_a.is_file() or not ref_b.is_file():
            raise SystemExit("missing --ref-a or --ref-b")
        y_a = spectrum_from_ref_path(ref_a, normalize_max)
        y_b = spectrum_from_ref_path(ref_b, normalize_max)
        meta = {
            "refs": {"mode": "file", "ref_a": str(ref_a), "ref_b": str(ref_b)},
            "inputs": {"mixture": str(mix_path)},
        }

    run_mixture_nnls(mix_path, y_a, y_b, normalize_max=normalize_max, json_out=cli.json_out, meta=meta)


if __name__ == "__main__":
    main()
