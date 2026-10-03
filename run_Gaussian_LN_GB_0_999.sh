#!/usr/bin/env bash
# Run seeds 0-999 of Gaussian_LN_GB with N parallel single-threaded processes (default 24).
set -u
cd "$(dirname "$0")"
N=${1:-24}
mkdir -p logs
seq 0 999 | xargs -P "$N" -I{} bash -c "./Gaussian_LN_GB {} > logs/Gaussian_LN_GB_seed{}.log 2>&1"
