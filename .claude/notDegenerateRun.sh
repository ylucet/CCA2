#!/usr/bin/env bash
# Runs .claude/notDegenerateProbe.m over the 8 KINK fixtures x {conjCPLQ, biconjCPLQ}, one MATLAB
# per job under `timeout`, J at a time. Every maxQuaPar:notDegenerate firing lands in
# $OUT/<job>.guard.log (the guard writes it when CCA2_GUARD_LOG is set), one log PER JOB so a
# firing names its fixture. Any nonempty log => REACHABLE, witnessed by that fixture.
#   bash .claude/notDegenerateRun.sh [J=3] [TIMEOUT=4h] [ONLY=regex over "shape c op"]
J=${1:-3}; TO=${2:-4h}; ONLY=${3:-.}
cd "$(dirname "$0")/.." || exit 1
OUT=.claude/notDegenerate-$(date +%Y%m%d-%H%M); mkdir -p "$OUT"
jobs_list=()
for op in conjCPLQ biconjCPLQ; do for shape in square2 poly2; do for c in -2 -0.5 0.5 2; do
  echo "$shape $c $op" | grep -Eq -- "$ONLY" && jobs_list+=("$shape $c $op"); done; done; done
run() { set -- $1; tag="$1_c$2_$3"
  export CCA2_GUARD_LOG="$(cygpath -w "$PWD/$OUT/$tag.guard.log" 2>/dev/null || echo "$PWD/$OUT/$tag.guard.log")"
  timeout "$TO" matlab -batch "addpath('.claude'); notDegenerateProbe('$1', $2, '$3')" > "$OUT/$tag.txt" 2>&1
  rc=$?; [ $rc -eq 124 ] && echo "OUTCOME TIMEOUT $TO" >> "$OUT/$tag.txt"
  echo "$(date +%T) done $tag rc=$rc: $(grep -m1 OUTCOME "$OUT/$tag.txt")" >> "$OUT/progress.txt"; }
export -f run; export OUT TO
printf '%s\n' "${jobs_list[@]}" | xargs -P "$J" -I{} bash -c 'run "{}"'
echo "FINISHED. Firings per job:" | tee -a "$OUT/progress.txt"
for f in "$OUT"/*.guard.log; do [ -s "$f" ] && echo "  $(basename "$f" .guard.log): $(grep -c . "$f")"; done | tee -a "$OUT/progress.txt"
