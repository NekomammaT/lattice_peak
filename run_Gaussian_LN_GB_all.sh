#!/usr/bin/env bash
# Run the ten seed-range jobs (1--999) concurrently and wait for all of them.
set -u

cd "$(dirname "$0")"

if [[ ! -x ./Gaussian_LN_GB ]]; then
  echo "error: ./Gaussian_LN_GB is missing or not executable; run 'make' first." >&2
  exit 1
fi

mkdir -p logs

scripts=(Gaussian_LN_GB.sh)
for job in {1..9}; do
  scripts+=("Gaussian_LN_GB${job}.sh")
done

pids=()

stop_jobs() {
  echo
  echo "interrupt received: stopping running jobs..." >&2
  for pid in "${pids[@]}"; do
    # Each job is a Bash process whose current child is Gaussian_LN_GB.
    # Stop that child first so terminating the wrapper cannot leave it behind.
    pkill -TERM -P "$pid" 2>/dev/null || true
    kill -TERM "$pid" 2>/dev/null || true
  done
  for pid in "${pids[@]}"; do
    wait "$pid" 2>/dev/null || true
  done
  exit 130
}

trap stop_jobs INT TERM

for script in "${scripts[@]}"; do
  if [[ ! -f "$script" ]]; then
    echo "error: $script is missing." >&2
    exit 1
  fi
  log="logs/${script%.sh}.log"
  echo "start: $script (log: $log)"
  bash "$script" >"$log" 2>&1 &
  pids+=("$!")
done

status=0
for n in "${!pids[@]}"; do
  if wait "${pids[n]}"; then
    echo "done:  ${scripts[n]}"
  else
    echo "failed: ${scripts[n]} (see logs/${scripts[n]%.sh}.log)" >&2
    status=1
  fi
done

exit "$status"
