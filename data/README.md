# Data Layout

This directory documents the **expected runtime data layout** for Mix2Phase-QPA.

## Expected Structure

| Path | Description |
| --- | --- |
| `pair1k_v3/{train,val,test}.lmdb` | Binary mixture evaluation dataset |
| `mp20_data/test.lmdb` | Single-phase lookup library used by the mp20 reference path |
| `mp20_data/test.csv` | Small metadata/index helper file |

At the moment, only lightweight files that are safe to ship with the repository should live here by default.  
Large LMDB assets should be added manually, distributed separately, or managed through an external release channel if redistribution is allowed.

## How Scripts Use This Directory

The scripts default to the following repository-relative paths:

- `data/pair1k_v3`
- `data/mp20_data/test.lmdb`

If your data lives elsewhere, override the CLI arguments such as:

- `--pair-dir`
- `--mp20-lmdb`

## Provenance and Redistribution

Some runtime data used during internal experiments originated from other local project directories.  
Before publishing this repository publicly, please verify:

1. whether those datasets may be redistributed;
2. whether they contain any fields that should be removed before release;
3. whether large binary artifacts should stay out of normal Git history.

If redistribution is not permitted, keep this repository code-only and document a separate data access path instead.
