# Mix2Phase-QPA

**Mix2Phase-QPA** stands for **Mixture 2-Phase Quantitative Phase Analysis**.  
它解决的是一个非常具体的 PXRD 定量问题：

> 当样品的主要物相已经缩小到 **A / B 两相**，并且这两相的单相参考谱可查或可生成时，如何用一条混合 PXRD 快速估计两相比例。

本项目采用的主线方法是：

`混合谱 + 两条参考谱 -> NNLS -> A/B 占比`

## Repository Contents

| Path | Description |
| --- | --- |
| [`scripts/`](scripts/) | Core inference, CIF fallback, and evaluation scripts |
| [`Mix2Phase-QPA_全文.md`](Mix2Phase-QPA_全文.md) | 中文完整说明 |
| [`Mix2Phase-QPA_full_EN.md`](Mix2Phase-QPA_full_EN.md) | English documentation |
| [`报告_Mix2Phase-QPA_项目总结.md`](报告_Mix2Phase-QPA_项目总结.md) | 中文汇报版总结 |
| [`results/`](results/) | Example generated artifacts and evaluation summaries |
| [`data/README.md`](data/README.md) | Data layout, provenance notes, and expected local paths |

## Quick Start

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

Daily inference from the repository root:

```bash
python3 scripts/algo_qpa_mixture_nnls.py \
  --mixture /path/to/mix.xy \
  --ref-a-mp20 <idx_a> \
  --ref-b-mp20 <idx_b>
```

Batch evaluation wrappers:

```bash
bash scripts/run_eval20.sh
bash scripts/run_fixed_pair_w10.sh
```

More script examples are in [`scripts/README.md`](scripts/README.md).

## Data Availability

This repository includes:

- scripts
- documentation
- small metadata files such as `data/mp20_data/test.csv`
- example result artifacts under `results/`

This repository does **not** currently redistribute the full runtime LMDB datasets by default.  
If you want to reproduce the bundled evaluations, you must place the required files at the expected local paths:

- `data/pair1k_v3/{train,val,test}.lmdb`
- `data/mp20_data/test.lmdb`

See [`data/README.md`](data/README.md) for the expected layout and provenance notes.  
Before publishing this repository publicly, make sure you have the right to redistribute any external data.

## Included Results

| Path | Description |
| --- | --- |
| [`results/eval20_summary.json`](results/eval20_summary.json) | 20-sample evaluation summary |
| [`results/q1_eval20_run/`](results/q1_eval20_run/) | Example mixed-pattern `.xy` files |
| [`results/fixed_pair_w10_summary.json`](results/fixed_pair_w10_summary.json) | Fixed-A/B sweep across 10 ratios |
| [`results/fixed_pair_w10_run/`](results/fixed_pair_w10_run/) | Synthetic mixed-pattern `.xy` files for the fixed-pair sanity check |

The fixed-pair 10-ratio experiment is a **sanity check** rather than an independent benchmark: the mixtures are synthesized from the same lookup references used during fitting.

## Notes for Public Release

- Choose and add a repository license before publishing.
- Confirm redistribution rights for any dataset copied from external sources.
- If full LMDB data is too large or cannot be redistributed, publish only scripts/docs/results and provide a separate data acquisition path.
