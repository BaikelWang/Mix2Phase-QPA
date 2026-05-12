#!/usr/bin/env python3
"""Fallback path: compute single-phase PXRD from CIF, then fit a mixed pattern with NNLS."""
from __future__ import annotations

import argparse
import json
import re
import sys
import tempfile
from pathlib import Path

import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import nnls

from mix2phase_common import (
    DEFAULT_MP20_LMDB,
    GRID,
    STEP,
    X_MAX,
    interpolate_pxrd,
    load_mp20_entry,
    require_existing_file,
)


def nnls_fit_stack(y_mix: np.ndarray, y_a: np.ndarray, y_b: np.ndarray) -> tuple[float, np.ndarray, float]:
    """NNLS on [y_a, y_b]; return w_A = alpha/(alpha+beta), coef, ||y-Xb||_2."""
    X = np.stack([y_a.reshape(-1), y_b.reshape(-1)], axis=1)
    coef, resid_norm = nnls(X, y_mix.reshape(-1))
    s = float(coef.sum())
    if s < 1e-15:
        w = 0.5
    else:
        w = float(np.clip(coef[0] / s, 0.0, 1.0))
    return w, coef, float(resid_norm)


def load_xy_two_column(path: Path) -> tuple[np.ndarray, np.ndarray]:
    xs, ys = [], []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = re.split(r"\s+", line.replace(",", " "))
        if len(parts) < 2:
            continue
        try:
            xs.append(float(parts[0]))
            ys.append(float(parts[1]))
        except ValueError:
            continue
    if len(xs) < 3:
        raise ValueError(f"Too few data points in {path}")
    xs = np.asarray(xs, dtype=np.float64)
    ys = np.asarray(ys, dtype=np.float64)
    order = np.argsort(xs)
    return xs[order], ys[order]


def _interp_theta_intensity(theta: np.ndarray, intensity: np.ndarray, normalize_max: bool) -> np.ndarray:
    f = interp1d(theta, intensity, kind="linear", bounds_error=False, fill_value=0.0)
    y = np.asarray(f(GRID), dtype=np.float64)
    if normalize_max and y.max() > 0:
        y = y / y.max()
    return y


def spectrum_from_xy(path: Path, normalize_max: bool) -> np.ndarray:
    theta, intensity = load_xy_two_column(path)
    return _interp_theta_intensity(theta, intensity, normalize_max)


def spectrum_from_mixture_path(path: Path, normalize_max: bool) -> np.ndarray:
    suf = path.suffix.lower()
    if suf == ".npy":
        arr = np.asarray(np.load(path, allow_pickle=False), dtype=np.float64)
        if arr.ndim != 1 or arr.size != len(GRID):
            raise ValueError(f"{path}: .npy mixture must be 1d with length {len(GRID)}")
        y = arr.copy()
        if normalize_max and y.max() > 0:
            y = y / y.max()
        return y

    if suf == ".npz":
        z = np.load(path)
        try:
            arr = None
            for key in ("y_mix", "pxrd_y_mix", "y", "intensity", "I"):
                if key in z.files:
                    arr = np.asarray(z[key], dtype=np.float64)
                    break
            if arr is None:
                for key in z.files:
                    arr = np.asarray(z[key], dtype=np.float64)
                    if arr.size >= 3:
                        break
                else:
                    raise ValueError("no usable array in .npz")
            if arr.ndim == 2 and arr.shape[1] >= 2:
                order = np.argsort(arr[:, 0])
                return _interp_theta_intensity(arr[order, 0].ravel(), arr[order, 1].ravel(), normalize_max)
            arr = arr.ravel()
            if arr.size != len(GRID):
                raise ValueError(f"{path}: 1d array in .npz must have length {len(GRID)}")
            y = arr.copy()
            if normalize_max and y.max() > 0:
                y = y / y.max()
            return y
        finally:
            z.close()

    return spectrum_from_xy(path, normalize_max)


def spectrum_from_cif(path: Path, wavelength: str, normalize_max: bool) -> np.ndarray:
    from pymatgen.analysis.diffraction.xrd import XRDCalculator
    from pymatgen.core import Structure

    structure = Structure.from_file(str(path))
    calc = XRDCalculator(wavelength=wavelength)
    pattern = calc.get_pattern(structure, two_theta_range=(0.0, X_MAX))
    _, py = interpolate_pxrd(
        X_MAX,
        np.asarray(pattern.x, dtype=np.float64),
        np.asarray(pattern.y, dtype=np.float64),
        num_fill=2,
        threshold=5,
        step=STEP,
    )
    y = np.asarray(py, dtype=np.float64)
    if normalize_max and y.max() > 0:
        y = y / y.max()
    return y


def residual_l2_convex(y_mix: np.ndarray, y_a: np.ndarray, y_b: np.ndarray, w: float) -> float:
    pred = w * y_a + (1.0 - w) * y_b
    return float(np.linalg.norm(y_mix - pred))


def run_pipeline(
    cif_a: Path,
    cif_b: Path,
    mixture: Path,
    wavelength: str,
    normalize_max: bool,
    json_out: str | None,
    *,
    emit: bool = True,
) -> dict:
    y_a = spectrum_from_cif(cif_a, wavelength, normalize_max)
    y_b = spectrum_from_cif(cif_b, wavelength, normalize_max)
    y_m = spectrum_from_mixture_path(mixture, normalize_max)
    w_hat, coef, resid_nnls = nnls_fit_stack(y_m, y_a, y_b)
    out = {
        "w_A_hat": w_hat,
        "w_B_hat": 1.0 - w_hat,
        "nnls_coef_alpha_beta": coef.tolist(),
        "residual_l2_nnls": resid_nnls,
        "residual_l2_convex_combo": residual_l2_convex(y_m, y_a, y_b, w_hat),
        "grid": {"x_max": X_MAX, "step": STEP, "n": int(len(GRID))},
        "normalize_max": normalize_max,
        "wavelength": wavelength,
        "inputs": {"cif_a": str(cif_a), "cif_b": str(cif_b), "mixture": str(mixture)},
        "pipeline": "cif_xrdc_nnls",
    }
    if emit:
        print(json.dumps(out, indent=2, ensure_ascii=False))
    if json_out:
        Path(json_out).parent.mkdir(parents=True, exist_ok=True)
        Path(json_out).write_text(json.dumps(out, indent=2, ensure_ascii=False), encoding="utf-8")
        if emit:
            print(f"[wrote] {json_out}", file=sys.stderr)
    return out


def run_demo(normalize_max: bool, json_out: str | None) -> None:
    from pymatgen.core import Lattice, Structure

    require_existing_file(DEFAULT_MP20_LMDB, purpose="mp20 lookup LMDB")
    e0 = load_mp20_entry(DEFAULT_MP20_LMDB, 0)
    e1 = load_mp20_entry(DEFAULT_MP20_LMDB, 1)

    def entry_to_cif(entry: dict, path: Path) -> None:
        lattice = np.asarray(entry["lattice_matrix"], dtype=np.float64)
        species = list(entry["atom_type"])
        frac = np.asarray(entry["atom_pos"], dtype=np.float64).reshape(-1, 3)
        structure = Structure(Lattice(lattice), species, frac, coords_are_cartesian=False)
        structure.to(fmt="cif", filename=str(path))

    with tempfile.TemporaryDirectory(prefix="mix2phase_demo_") as tmpdir:
        tmpdir = Path(tmpdir)
        p_a = tmpdir / "demo_A.cif"
        p_b = tmpdir / "demo_B.cif"
        entry_to_cif(e0, p_a)
        entry_to_cif(e1, p_b)

    _, ya = interpolate_pxrd(
        X_MAX,
        np.asarray(e0["pxrd_x"], dtype=np.float64),
        np.asarray(e0["pxrd_y"], dtype=np.float64),
        num_fill=2,
        threshold=5,
        step=STEP,
    )
    _, yb = interpolate_pxrd(
        X_MAX,
        np.asarray(e1["pxrd_x"], dtype=np.float64),
        np.asarray(e1["pxrd_y"], dtype=np.float64),
        num_fill=2,
        threshold=5,
        step=STEP,
    )
    ya = np.asarray(ya, dtype=np.float64)
    yb = np.asarray(yb, dtype=np.float64)
    if normalize_max:
        if ya.max() > 0:
            ya = ya / ya.max()
        if yb.max() > 0:
            yb = yb / yb.max()

    w_true = 0.37
    y_m = w_true * ya + (1.0 - w_true) * yb + np.random.default_rng(0).normal(0, 0.005, ya.shape)
        mix_path = tmpdir / "demo_mix.xy"
        with mix_path.open("w", encoding="utf-8") as f:
            for x, y in zip(GRID, y_m):
                f.write(f"{x:.6f} {y:.8f}\n")

        out = run_pipeline(
            p_a,
            p_b,
            mix_path,
            wavelength="CuKa",
            normalize_max=normalize_max,
            json_out=None,
            emit=False,
        )
        out["demo_w_true"] = w_true
        print(json.dumps(out, indent=2, ensure_ascii=False))
        if json_out:
            Path(json_out).parent.mkdir(parents=True, exist_ok=True)
            Path(json_out).write_text(json.dumps(out, indent=2, ensure_ascii=False), encoding="utf-8")
            print(f"[wrote] {json_out}", file=sys.stderr)


def main() -> None:
    ap = argparse.ArgumentParser(description="Mix2Phase-QPA fallback: 2x CIF + mixed PXRD -> NNLS")
    ap.add_argument("--cif-a", type=str, default=None)
    ap.add_argument("--cif-b", type=str, default=None)
    ap.add_argument("--mixture", type=str, default=None)
    ap.add_argument("--wavelength", type=str, default="CuKa", help="pymatgen XRDCalculator wavelength key")
    ap.add_argument("--no-normalize", action="store_true", help="skip per-spectrum max=1")
    ap.add_argument("--json-out", type=str, default=None)
    ap.add_argument("--demo", action="store_true", help="use two mp20 entries to synthesize a fake mixed PXRD")
    cli = ap.parse_args()
    normalize_max = not cli.no_normalize

    if cli.demo:
        run_demo(normalize_max, cli.json_out)
        return

    if not cli.cif_a or not cli.cif_b or not cli.mixture:
        ap.error("need --cif-a --cif-b --mixture (or use --demo)")

    p_a = Path(cli.cif_a)
    p_b = Path(cli.cif_b)
    p_m = Path(cli.mixture)
    for p in (p_a, p_b, p_m):
        if not p.is_file():
            raise SystemExit(f"missing file: {p}")

    run_pipeline(p_a, p_b, p_m, cli.wavelength, normalize_max, cli.json_out, emit=True)


if __name__ == "__main__":
    main()
