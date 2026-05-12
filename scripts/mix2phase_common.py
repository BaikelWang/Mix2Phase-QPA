#!/usr/bin/env python3
"""Shared helpers for the standalone Mix2Phase-QPA bundle."""
from __future__ import annotations

import gzip
import pickle
from pathlib import Path

import numpy as np
from scipy.interpolate import interp1d

TASK0511_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MP20_LMDB = TASK0511_ROOT / "data" / "mp20_data" / "test.lmdb"
DEFAULT_PAIR_DIR = TASK0511_ROOT / "data" / "pair1k_v3"

X_MAX = 120.0
STEP = 0.1
GRID = np.arange(0.0, X_MAX, STEP, dtype=np.float64)

_MP20_CURSOR_LIST_CACHE: dict[str, list] = {}


def random_fill(
    known_peaks: np.ndarray,
    fill_length: int,
    x_min: float,
    x_max: float,
    dist_margin: float = 0.2,
    candidate_multiplier: int = 50,
) -> np.ndarray:
    known_peaks = np.asarray(known_peaks, dtype=np.float64)
    num_to_fill = fill_length - len(known_peaks)
    if num_to_fill <= 0:
        return np.zeros(0, dtype=np.float64)
    candidate_count = num_to_fill * candidate_multiplier
    candidates = np.random.uniform(x_min, x_max, candidate_count)
    diffs = np.abs(candidates[:, None] - known_peaks[None, :])
    mask = np.all(diffs > dist_margin, axis=1)
    return np.asarray(candidates[mask][:num_to_fill], dtype=np.float64)


def interpolate_pxrd(
    pxrd_x_range_max: float,
    pxrd_x: np.ndarray,
    pxrd_y: np.ndarray,
    num_fill: int = 128,
    threshold: float = 5.0,
    step: float = STEP,
) -> tuple[np.ndarray, np.ndarray]:
    mask = np.asarray(pxrd_y, dtype=np.float64) > threshold
    pxrd_x = np.asarray(pxrd_x, dtype=np.float64)[mask]
    pxrd_y = np.asarray(pxrd_y, dtype=np.float64)[mask]
    pad_pxrd_x = random_fill(pxrd_x, num_fill, 0.0, pxrd_x_range_max)
    pxrd_x = np.concatenate([pxrd_x, pad_pxrd_x])
    pxrd_y = np.concatenate([pxrd_y, np.zeros(pad_pxrd_x.shape[0], dtype=np.float64)])
    sorted_index = np.argsort(pxrd_x)
    pxrd_x = pxrd_x[sorted_index]
    pxrd_y = pxrd_y[sorted_index]
    if pxrd_x[0] > 0.0:
        pxrd_x = np.insert(pxrd_x, 0, 0.0)
        pxrd_y = np.insert(pxrd_y, 0, 0.0)
    if pxrd_x[-1] < pxrd_x_range_max:
        pxrd_x = np.append(pxrd_x, pxrd_x_range_max)
        pxrd_y = np.append(pxrd_y, 0.0)
    new_x = np.arange(0.0, pxrd_x_range_max, step)
    interpolator = interp1d(pxrd_x, pxrd_y, kind="linear", fill_value="extrapolate")
    new_y = interpolator(new_x)
    return new_x.astype(np.float64), np.asarray(new_y, dtype=np.float64)


def load_mp20_cursor_list(lmdb_path: Path) -> list:
    import lmdb

    key = str(lmdb_path.resolve())
    if key in _MP20_CURSOR_LIST_CACHE:
        return _MP20_CURSOR_LIST_CACHE[key]
    env = lmdb.open(str(lmdb_path), subdir=False, readonly=True, lock=False)
    out = []
    with env.begin() as txn:
        for _, value in txn.cursor():
            try:
                out.append(pickle.loads(gzip.decompress(value)))
            except Exception:
                out.append(pickle.loads(value))
    env.close()
    _MP20_CURSOR_LIST_CACHE[key] = out
    return out


def load_mp20_entry(lmdb_path: Path, cursor_index: int) -> dict:
    entries = load_mp20_cursor_list(lmdb_path)
    if cursor_index < 0 or cursor_index >= len(entries):
        raise IndexError(f"{lmdb_path}: cursor_index {cursor_index} out of range [0,{len(entries)})")
    return entries[cursor_index]


def composition_label(atom_types) -> str:
    from collections import Counter

    types = [str(x) for x in atom_types]
    uniq = sorted(set(types))
    counter = Counter(types)
    parts = []
    for el in uniq:
        n = counter[el]
        parts.append(f"{el}{n}" if n > 1 else el)
    return "".join(parts)


def safe_rel(path: Path, base: Path) -> str:
    rp = path.resolve()
    rb = base.resolve()
    try:
        return str(rp.relative_to(rb))
    except ValueError:
        return str(rp)


def write_mix_xy(path: Path, y_mix: np.ndarray) -> None:
    y_mix = np.asarray(y_mix, dtype=np.float64).reshape(-1)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        for x, y in zip(GRID, y_mix):
            f.write(f"{x:.6f} {y:.8f}\n")


def require_existing_file(path: Path, *, purpose: str) -> Path:
    if not path.is_file():
        raise FileNotFoundError(
            f"Missing {purpose}: {path}\n"
            "If you are using the public repository, place the required runtime data under the "
            "expected `data/` layout or override the relevant CLI path argument."
        )
    return path


def require_existing_dir(path: Path, *, purpose: str) -> Path:
    if not path.is_dir():
        raise FileNotFoundError(
            f"Missing {purpose}: {path}\n"
            "If you are using the public repository, place the required runtime data under the "
            "expected `data/` layout or override the relevant CLI path argument."
        )
    return path
