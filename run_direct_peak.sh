#!/usr/bin/env bash
# Direct sampling with Gaussian_LN_peak (NL=256, As=5e-3, s2=0.1): seeds S0..S1-1, N parallel jobs
# usage: ./run_direct_peak.sh [S0=0] [S1=10] [N=10]
set -u
cd "$(dirname "$0")"
S0=${1:-0}; S1=${2:-10}; N=${3:-10}
mkdir -p direct_256 logs_direct
seq $S0 $((S1-1)) | xargs -P "$N" -I{} bash -c 'OMP_NUM_THREADS=1 ./Gaussian_LN_peak {} > logs_direct/s{}.log 2>&1'
