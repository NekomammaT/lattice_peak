#!/usr/bin/env bash
# Scan biascoeff: seeds 0..NS-1 for each coefficient, N parallel single-threaded jobs.
# usage: ./run_bias_scan.sh [NS=200] [N=26] [bc list...]
set -u
cd "$(dirname "$0")"
NS=${1:-200}; N=${2:-26}; shift 2 2>/dev/null
BCS=${@:-4 5 6 7 8 10 12.5}
mkdir -p scan_bias logs_scan
for bc in $BCS; do for ((s=0;s<NS;s++)); do echo "$s $bc"; done; done | \
  xargs -P "$N" -L1 bash -c './Gaussian_LN_GB_scan $0 $1 > logs_scan/s$0_bc$1.log 2>&1'
