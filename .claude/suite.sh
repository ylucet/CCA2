#!/usr/bin/env bash
# Run the whole suite with a PER-SUITE TIMEOUT, one MATLAB process per suite.
#
#   .claude/suite.sh                  # default timeout per suite
#   CCA2_TEST_TIMEOUT=1800 .claude/suite.sh
#   CCA2_TEST_TIMEOUT=0 .claude/suite.sh          # no timeout at all
#   .claude/suite.sh --fast                       # fast bucket, one process, budget 5 min
#   .claude/suite.sh --normal                     # normal bucket, one process, budget 10 min
#   .claude/suite.sh --slow                       # slow (symbolic) bucket, one process per suite
#   .claude/suite.sh regionTest conjCPLQTest      # only these suites
#
# WHY A SHELL DRIVER RATHER THAN .claude/suite.m. MATLAB's unittest framework has no per-test
# timeout, and nothing inside a MATLAB process can interrupt a symbolic call that is not going to
# return -- `solve`/`simplify` are not interruptible from M-code. The only reliable timeout is a
# process the OS can kill, so each suite gets its own `matlab -batch` under `timeout`. The cost is
# one MATLAB startup per suite (~8 s); the benefit is that one wedged suite no longer costs the
# whole run.
#
# This matters here because the suite is genuinely slow and getting slower: conjCPLQTest is around
# an hour, of which indefiniteTriangleThreeConvexEdgesUsesStep3 alone is ~33 minutes. A timeout
# that is too tight will report TIMEOUT for a suite that is merely slow, which is why the default
# is generous and the override is the first thing documented above.
#
# Exit status: 0 if every suite passed, 1 if any failed, timed out, or errored.

set -u
CCA2DIR="${CCA2DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
TIMEOUT="${CCA2_TEST_TIMEOUT:-3600}"
cd "$CCA2DIR" || exit 1

# ---- BUCKETS ---------------------------------------------------------------------------------
# fast   : the whole bucket in ONE MATLAB, budget 5 minutes.
# normal : the whole bucket in ONE MATLAB, budget 10 minutes.
# slow   : one MATLAB per suite under `timeout`, no budget -- these are the symbolic ones.
#
# WHY THE FAST AND NORMAL BUCKETS SHARE A PROCESS. Measured per suite (18 suites, 6 concurrent):
# every fast suite came in at 41-61 s, and that is almost entirely MATLAB startup plus toolbox
# load, not test time -- 26 startups is minutes of pure overhead. Paying it once is what makes a
# 5-minute budget reachable at all. The slow bucket keeps one process per suite because that is
# the only thing the OS can kill: `solve`/`simplify` are not interruptible from M-code, so a
# wedged symbolic call can only be bounded by killing the process (see the note above).
#
# SLOW MEANS SYMBOLIC. Every member of the slow bucket runs the Symbolic Math Toolbox pipeline;
# every fast member is closed-form numerics. That is the whole distribution -- so moving a case
# from the symbolic path to closed form is also what moves its test out of the slow bucket.
FAST_SUITES=(PLQVCTest QuaParTest RatParTest RatPolTest addQuaParTest addQuaPolTest
             clipArcByHalfPlaneTest conjPieceCPLQTest convEnvCPLQTest infConvTest
             lasryLionsTest maxQuaParTest moreauTest proxAverageTest regionTest
             testSymbolicFunction)
NORMAL_SUITES=(cplqAdapterTest functionNDomainTest testfunctionNDomain)
SLOW_SUITES=(biconjCPLQTest biconjugateTest conjCPLQTest testMaxMultiRegion testRegion
             testcPLQ unboundedFaceTest)

ONEPROC=0
case "${1:-}" in
    --fast)   SUITES=("${FAST_SUITES[@]}");   ONEPROC=1; shift ;;
    --normal) SUITES=("${NORMAL_SUITES[@]}"); ONEPROC=1; shift ;;
    --slow)   SUITES=("${SLOW_SUITES[@]}");            shift ;;
    *)
        if [ "$#" -gt 0 ]; then
            SUITES=("$@")
        else
            mapfile -t SUITES < <(ls *Test.m test*.m 2>/dev/null | sed 's/\.m$//' | sort -u)
        fi
        ;;
esac

printf '=== suite: %s (timeout %ss per suite) ===\n' "$CCA2DIR" "${TIMEOUT:-none}"
tp=0; tf=0; ti=0; tt=0; bad=0
run_started=$(date +%s)

if [ "$ONEPROC" -eq 1 ]; then
    # One process for the whole bucket: build a cell array of suite names for runtests.
    list=$(printf "'%s'," "${SUITES[@]}"); list="{${list%,}}"
    script="r = runtests($list);
            fprintf('COUNTS %d %d %d\\n', sum([r.Passed]), sum([r.Failed]), sum([r.Incomplete]));
            for k = 1:numel(r), if r(k).Failed, fprintf('FAILED %s\\n', r(k).Name); end, end"
    if [ "$TIMEOUT" = "0" ]; then
        out=$(matlab -batch "$script" 2>&1); rc=$?
    else
        out=$(timeout "$TIMEOUT" matlab -batch "$script" 2>&1); rc=$?
    fi
    elapsed=$(( $(date +%s) - run_started ))
    if [ "$rc" -eq 124 ]; then
        printf 'BUCKET TIMEOUT after %ss\n' "$TIMEOUT"; exit 1
    fi
    counts=$(printf '%s\n' "$out" | grep -m1 '^COUNTS ')
    if [ -z "$counts" ]; then
        printf 'BUCKET ERRORED (no result; exit %s)\n' "$rc"; exit 1
    fi
    read -r _ np nf ni <<<"$counts"
    printf '%s\n' "$out" | grep '^FAILED ' | sed 's/^FAILED /  FAIL /'
    printf 'TOTALS pass=%d fail=%d incomplete=%d over %d suites in %ds\n' \
           "$np" "$nf" "$ni" "${#SUITES[@]}" "$elapsed"
    [ "$nf" -eq 0 ] && [ "$ni" -eq 0 ] && exit 0
    exit 1
fi

for s in "${SUITES[@]}"; do
    s_started=$(date +%s)
    script="r = runtests('$s');
            fprintf('COUNTS %d %d %d\\n', sum([r.Passed]), sum([r.Failed]), sum([r.Incomplete]));
            for k = 1:numel(r), if r(k).Failed, fprintf('FAILED %s\\n', r(k).Name); end, end"
    if [ "$TIMEOUT" = "0" ]; then
        out=$(matlab -batch "$script" 2>&1); rc=$?
    else
        out=$(timeout "$TIMEOUT" matlab -batch "$script" 2>&1); rc=$?
    fi

    elapsed=$(( $(date +%s) - s_started ))
    if [ "$rc" -eq 124 ]; then
        printf 'SUITE %-28s TIMEOUT after %ss\n' "$s" "$TIMEOUT"
        tt=$((tt+1)); bad=1
        continue
    fi
    counts=$(printf '%s\n' "$out" | grep -m1 '^COUNTS ')
    if [ -z "$counts" ]; then
        printf 'SUITE %-28s ERRORED (no result; exit %s)\n' "$s" "$rc"
        bad=1
        continue
    fi
    read -r _ np nf ni <<<"$counts"
    printf 'SUITE %-28s pass=%3d fail=%3d incomplete=%3d  %5ds\n' "$s" "$np" "$nf" "$ni" "$elapsed"
    printf '%s\n' "$out" | grep '^FAILED ' | sed 's/^FAILED /  FAIL /'
    tp=$((tp+np)); tf=$((tf+nf)); ti=$((ti+ni))
    [ "$nf" -eq 0 ] && [ "$ni" -eq 0 ] || bad=1
done

printf 'TOTALS pass=%d fail=%d incomplete=%d timeout=%d over %d suites\n' \
       "$tp" "$tf" "$ti" "$tt" "${#SUITES[@]}"
exit "$bad"
