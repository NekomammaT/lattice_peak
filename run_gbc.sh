#!/usr/bin/env bash
# compaction-peak IS: seeds S0..S1-1, N parallel single-threaded jobs, bias coefficient BC.
# Appends to gbc/muk_bc<BC>.csv.   usage: ./run_gbc.sh BC S0 S1 [N=24]
set -u
cd "$(dirname "$0")"
BC=$1; S0=$2; S1=$3; N=${4:-24}
mkdir -p gbc logs_gbc
seq "$S0" $((S1-1)) | xargs -P "$N" -I{} bash -c "arch -arm64 ./Gaussian_LN_GBC {} $BC > logs_gbc/s{}_bc$BC.log 2>&1"
