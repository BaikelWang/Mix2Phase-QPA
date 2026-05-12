#!/usr/bin/env bash
# Generate the fixed-pair 10-ratio sanity-check results.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

python3 "${ROOT_DIR}/scripts/eval_fixed_pair_w10.py" "$@"
