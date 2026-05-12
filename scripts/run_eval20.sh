#!/usr/bin/env bash
# Reproduce the 20-sample evaluation from the repository root.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

python3 "${ROOT_DIR}/scripts/algo_qpa_q1_eval_pair50.py" \
  --n 20 --seed 42 --split test \
  --workdir "${ROOT_DIR}/results/q1_eval20_run" \
  --json-out "${ROOT_DIR}/results/eval20_summary.json" \
  --md-out "${ROOT_DIR}/results/_autogen/eval20_report.md"
echo "[done] JSON -> ${ROOT_DIR}/results/eval20_summary.json"
echo "      See ${ROOT_DIR}/Mix2Phase-QPA_全文.md for the main narrative summary"
