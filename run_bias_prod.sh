#!/usr/bin/env bash
# Production run at a fixed biascoeff: seeds S0..S1-1, N parallel single-threaded jobs.
# Appends to scan_bias/muk_bc<BC>.csv.   usage: ./run_bias_prod.sh BC S0 S1 [N=30]
set -u
cd "$(dirname "$0")"
BC=$1; S0=$2; S1=$3; N=${4:-30}
mkdir -p scan_bias logs_scan
seq "$S0" $((S1-1)) | xargs -P "$N" -I{} bash -c "./Gaussian_LN_GB_scan {} $BC > logs_scan/s{}_bc$BC.log 2>&1"
