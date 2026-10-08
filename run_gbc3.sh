#!/usr/bin/env bash
# local compaction-peak IS at origin window (As=1e-2): usage ./run_gbc3.sh BC S0 S1 [N=24] -> gbc/ev_As1e-2_bc<BC>.csv
set -u
cd "$(dirname "$0")"
BC=$1; S0=$2; S1=$3; N=${4:-24}
mkdir -p gbc logs_gbc
seq "$S0" $((S1-1)) | xargs -P "$N" -I{} bash -c "arch -arm64 ./Gaussian_LN_GBC3 {} $BC > logs_gbc/ev_s{}_bc$BC.log 2>&1"
