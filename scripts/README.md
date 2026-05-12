# Mix2Phase-QPA — 脚本入口

本目录包含 Mix2Phase-QPA 的核心脚本。

默认数据路径也已切换到 **`../data/`**：

- `../data/pair1k_v3/*.lmdb`
- `../data/mp20_data/test.lmdb`

| 脚本 | 作用 |
| --- | --- |
| [`algo_qpa_mixture_nnls.py`](algo_qpa_mixture_nnls.py) | 主推路径：混合谱 + mp20 查表 / 文件参考谱 → NNLS |
| [`algo_qpa_q1_cif_xy.py`](algo_qpa_q1_cif_xy.py) | 备选路径：A/B 的 CIF → XRDC → NNLS |
| [`algo_qpa_q1_eval_pair50.py`](algo_qpa_q1_eval_pair50.py) | 随机抽样 pair 数据集做批量评测 |
| [`mix2phase_common.py`](mix2phase_common.py) | 本项目共用常量、插值、mp20 读取与写谱辅助函数 |
| [`run_eval20.sh`](run_eval20.sh) | 重新生成 `results/eval20_summary.json` 与 `results/q1_eval20_run/*/mix.xy` |
| [`run_fixed_pair_w10.sh`](run_fixed_pair_w10.sh) | 固定一对 A/B，生成 10 个不同比例的横向对比结果 |
| [`make_readme_figures.py`](make_readme_figures.py) | 由 `results/eval20_summary.json` 生成 `docs/figures/pxrd_demo_pair0008.png` 与 `docs/figures/accuracy_eval20.png` |

**日常推理（混合谱 + mp20 查表）**：

```bash
cd /path/to/Mix2Phase-QPA
python3 scripts/algo_qpa_mixture_nnls.py \
  --mixture /path/to/mix.xy \
  --ref-a-mp20 <idx_a> \
  --ref-b-mp20 <idx_b>
```

索引语义（mp20 cursor 序）见仓库根文档 **[`../Mix2Phase-QPA_全文.md`](../Mix2Phase-QPA_全文.md)** 第三章。

**CIF 备选路径**：

```bash
cd /path/to/Mix2Phase-QPA
python3 scripts/algo_qpa_q1_cif_xy.py \
  --cif-a /path/to/A.cif --cif-b /path/to/B.cif --mixture /path/to/mix.xy
```

**横向对比（固定 A/B，扫描 10 个比例）**：

```bash
bash scripts/run_fixed_pair_w10.sh
```
