#!/bin/bash
# Thread sweep of the 200 bench fixtures: for each sweep, for each thread count, run `bench`
# (first call + 3 repetitions, best per fixture, memos emptied before every call) and record the
# load average before and after. Aggregate with tools/par_sweep.py. Shared machine: take minima.
#   tools/par_sweep.sh OUTDIR SWEEPS "1 2 4 8 16 32 64"
out=${1:-bench/parallel/runs}; sweeps=${2:-3}; ts=${3:-"1 2 4 8 16 32 64"}
mkdir -p "$out"
for s in $(seq 1 "$sweeps"); do
  for t in $ts; do
    if [ "$t" = 1 ]; then env="CCA2_PAR=0"; else env="CCA2_THREADS=$t"; fi
    l0=$(cut -d' ' -f1 /proc/loadavg)
    env $env ./target/release/bench ../python/bench/fixtures.json "$out/s${s}_t${t}.jsonl" 3 20 > /dev/null 2>&1
    echo "sweep $s T=$t load_before=$l0 load_after=$(cut -d' ' -f1 /proc/loadavg)" >> "$out/loads.txt"
  done
done
