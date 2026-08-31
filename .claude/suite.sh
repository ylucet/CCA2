#!/usr/bin/env bash
# Run the whole suite with a PER-SUITE TIMEOUT, one MATLAB process per suite.
#
#   .claude/suite.sh                  # default timeout per suite
#   CCA2_TEST_TIMEOUT=1800 .claude/suite.sh
#   CCA2_TEST_TIMEOUT=0 .claude/suite.sh          # no timeout at all
#   .claude/suite.sh --fast                       # fast bucket, one process, budget 5 min
#   .claude/suite.sh --fast --coverage            # ...plus a Cobertura report under
#                                                  #    .claude/coverage/cobertura.xml
#   .claude/suite.sh --normal                     # normal bucket, one process, budget 10 min
#   .claude/suite.sh --slow                       # slow (symbolic) bucket, one process per suite
#   .claude/suite.sh --slow -j 4                  # ...four at a time
#   .claude/suite.sh --verylong -j 4              # the multi-minute pipeline suites
#   CCA2_TEST_JOBS=4 .claude/suite.sh --slow      # same, via the environment
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
# fast     : the whole bucket in ONE MATLAB, budget 5 minutes.
# normal   : the whole bucket in ONE MATLAB, budget 10 minutes.
# slow     : one MATLAB per suite under `timeout`, no budget -- the symbolic ones.
# verylong : the same, for suites where a SINGLE test runs for many minutes.
#
# WHY VERYLONG EXISTS (2026-08-19). Measured: the slow bucket was 7207 s, and 6333 s of that was
# two suites -- testcPLQ and testMaxMultiRegion -- whose individual tests run 15 to 53 minutes
# each. They are whole-pipeline runs (triangulate -> envelope -> conjugate -> Step 3 ->
# biconjugate) on fixtures big enough that Step 3's symbolic max dominates. Keeping them in
# `--slow` meant the bucket you run after touching region.m cost two hours, so it was not run --
# and that is exactly how a stale expectation in conjCPLQTest went unnoticed for a day.
#
# Split out, `--slow` is the four suites carrying real assertions and finishes in minutes;
# `--verylong` is the pipeline endurance run, for before a tag or overnight.
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
#
# ---- PARALLELISM (-j, slow bucket only) ------------------------------------------------------
# The per-suite loop was always SERIAL: one process per suite bought KILLABILITY, never
# concurrency. Measured 2026-08-19, serial: 7207 s over seven suites, of which testcPLQ 3562 s and
# testMaxMultiRegion 2935 s -- 90% of the run in two suites. So fanning out whole SUITES floors
# the wall clock at 59 minutes, and the only way past that floor is to split those two.
#
# HOW THE JOB LIST IS BUILT. Every suite in SHARD_SUITES is expanded to one job PER TEST
# (`runtests('suite/test')`); every other suite stays one job. Jobs are then run at most
# CCA2_TEST_JOBS at a time, so wall clock is about max(longest single job, total / jobs).
#
# WHY 4 AND NOT 16. Each job is a MATLAB plus its own MuPAD engine, and each holds a LICENCE from
# the campus server; this machine also runs other projects' work. Four gets nearly all of the
# available speed-up because the distribution is so top-heavy -- past the first few workers there
# is nothing left to hand out -- and leaves the machine usable. Raise it with -j if the licences
# and the cores are genuinely idle.
#
# WHY FAST AND NORMAL ARE NOT FANNED OUT. They are 87 s and 215 s of mostly MATLAB startup; more
# processes would cost more than they save. -j is accepted and ignored for those buckets.
#
# ORDERING. Parallel jobs finish out of order, so results are buffered per job and printed in JOB
# ORDER at the end -- a run's output does not depend on how the scheduler happened to interleave
# it. Per-job elapsed times are printed, which is what tells you where to shard next.
FAST_SUITES=(PLQVCTest QuaParTest RatParTest RatPolTest addQuaParTest addQuaPolTest
             clipArcByHalfPlaneTest conicMeetTest conjConvexPolygonTest conjEdgeLowerBoundTest
             conjAffinePLQTest conjPieceCPLQTest conjSymFreeTest
             convEnvCPLQTest
             exactQTest frameAndFanTest infConvTest
             lasryLionsTest maxQuaParTest mergeSameQuadFacesTest meshPredicateTest
             moreauTest proxAverageTest
             ratQTest regionTest regionUnitTest testSymbolicFunction)
NORMAL_SUITES=(cplqAdapterTest functionNDomainTest regionMinusTest testfunctionNDomain)
# Split PER TEST when running in parallel: with tests this long, one job per SUITE floors the
# wall clock at the suite's whole runtime.
SHARD_SUITES=(testcPLQ testMaxMultiRegion)
SLOW_SUITES=(biconjCPLQTest biconjugateTest conjCPLQTest testRegion unboundedFaceTest)
# Individual tests here run for many minutes; see the VERYLONG note above.
VERYLONG_SUITES=(testcPLQ testMaxMultiRegion)

run_job () {
    # $1 = job name ("suite" or "suite/test"), $2 = index. Writes ONE result file and prints
    # nothing, so that concurrent jobs cannot interleave their output.
    local job="$1" idx="$2" started out rc elapsed counts
    started=$(date +%s)
    local script="r = runtests('$job');
            fprintf('COUNTS %d %d %d\n', sum([r.Passed]), sum([r.Failed]), sum([r.Incomplete]));
            for k = 1:numel(r), if r(k).Failed, fprintf('FAILED %s\n', r(k).Name); end, end"
    if [ "$TIMEOUT" = "0" ]; then
        out=$(matlab -batch "$script" 2>&1); rc=$?
    else
        out=$(timeout "$TIMEOUT" matlab -batch "$script" 2>&1); rc=$?
    fi
    elapsed=$(( $(date +%s) - started ))
    {
        if [ "$rc" -eq 124 ]; then
            printf 'STATUS timeout %s
' "$elapsed"
        else
            counts=$(printf '%s
' "$out" | grep -m1 '^COUNTS ')
            if [ -z "$counts" ]; then
                printf 'STATUS errored %s %s
' "$elapsed" "$rc"
            else
                printf 'STATUS ok %s %s
' "$elapsed" "$counts"
                printf '%s
' "$out" | grep '^FAILED ' || true
            fi
        fi
    } > "$RESDIR/$(printf '%04d' "$idx")"
}

# ---- SINGLE-JOB MODE (internal) ----------------------------------------------------------
# `suite.sh --run-one "<idx><TAB><job>"` runs exactly one job and writes its result file into
# CCA2_RESDIR. Used only by the xargs workers below; nothing else should call it.
if [ "${1:-}" = "--run-one" ]; then
    RESDIR="${CCA2_RESDIR:?--run-one needs CCA2_RESDIR}"
    printf '%s' "$2" | { IFS=$'	' read -r _idx _job; run_job "$_job" "$_idx"; }
    exit 0
fi

ONEPROC=0
JOBS="${CCA2_TEST_JOBS:-1}"
case "${1:-}" in
    --fast)   SUITES=("${FAST_SUITES[@]}");   ONEPROC=1; shift ;;
    --normal) SUITES=("${NORMAL_SUITES[@]}"); ONEPROC=1; shift ;;
    --slow)   SUITES=("${SLOW_SUITES[@]}");            shift ;;
    --verylong) SUITES=("${VERYLONG_SUITES[@]}");     shift ;;
    --all)    SUITES=("${FAST_SUITES[@]}" "${NORMAL_SUITES[@]}" "${SLOW_SUITES[@]}"                       "${VERYLONG_SUITES[@]}");       shift ;;
    *)
        if [ "$#" -gt 0 ] && [ "${1:-}" != "-j" ]; then
            SUITES=("$@")
            set --
        else
            mapfile -t SUITES < <(ls *Test.m test*.m 2>/dev/null | sed 's/\.m$//' | sort -u)
        fi
        ;;
esac
# -j N anywhere after the bucket flag. Any remaining NAMES narrow the bucket to those suites --
# silently dropping them (which an earlier version of this parsing did) turns
# `--slow -j 4 regionTest` into a two-hour run when the user asked for one suite.
COVERAGE=0
EXPLICIT=()
while [ "$#" -gt 0 ]; do
    case "$1" in
        -j) JOBS="${2:-1}"; shift 2 ;;
        -j*) JOBS="${1#-j}"; shift ;;
        --coverage) COVERAGE=1; shift ;;
        -*) printf 'suite.sh: unknown option "%s"
' "$1" >&2; exit 2 ;;
        *) EXPLICIT+=("$1"); shift ;;
    esac
done
if [ "${#EXPLICIT[@]}" -gt 0 ]; then
    SUITES=("${EXPLICIT[@]}")
fi
case "$JOBS" in
    ''|*[!0-9]*) printf 'suite.sh: -j needs a positive integer, got "%s"
' "$JOBS" >&2; exit 2 ;;
esac
[ "$JOBS" -lt 1 ] && JOBS=1

printf '=== suite: %s (timeout %ss per suite) ===\n' "$CCA2DIR" "${TIMEOUT:-none}"
tp=0; tf=0; ti=0; tt=0; bad=0
run_started=$(date +%s)

if [ "$COVERAGE" -eq 1 ] && [ "$ONEPROC" -ne 1 ]; then
    printf 'suite.sh: --coverage only works with --fast or --normal (one MATLAB process)\n' >&2
    exit 2
fi

if [ "$ONEPROC" -eq 1 ]; then
    # One process for the whole bucket: build a cell array of suite names for runtests.
    list=$(printf "'%s'," "${SUITES[@]}"); list="{${list%,}}"
    if [ "$COVERAGE" -eq 1 ]; then
        # Same suites, run through matlab.unittest.TestRunner instead of the plain `runtests`
        # convenience function, so a CodeCoveragePlugin can be attached. Coverage is measured over
        # every .m file directly under the repo root (one level, not test files or .claude/) --
        # that is where the production code lives; test classes measuring their own coverage
        # would just show 100% and dilute the number that matters.
        mkdir -p "$CCA2DIR/.claude/coverage"
        # `pwd` INSIDE matlab, not bash's $CCA2DIR: bash's $(pwd) is a Git-Bash-style path
        # (/c/Users/...) that native Windows MATLAB does not understand as a string literal --
        # measured directly, CoberturaFormat('/tmp/...') silently wrote nothing. matlab -batch
        # inherits the OS-level cwd bash already `cd`'d into above, which IS the correct native
        # path, so let MATLAB read its own pwd rather than re-embedding bash's version of it.
        script="import matlab.unittest.TestSuite
                import matlab.unittest.TestRunner
                import matlab.unittest.plugins.CodeCoveragePlugin
                import matlab.unittest.plugins.codecoverage.CoberturaFormat
                suite = [$(printf "TestSuite.fromClass(?%s), " "${SUITES[@]}" | sed 's/, $//')];
                runner = TestRunner.withTextOutput;
                runner.addPlugin(CodeCoveragePlugin.forFolder(pwd, ...
                    'IncludingSubfolders', false, ...
                    'Producing', CoberturaFormat(fullfile(pwd,'.claude','coverage','cobertura.xml'))));
                r = runner.run(suite);
                fprintf('COUNTS %d %d %d\\n', sum([r.Passed]), sum([r.Failed]), sum([r.Incomplete]));
                for k = 1:numel(r), if r(k).Failed, fprintf('FAILED %s\\n', r(k).Name); end, end"
    else
        script="r = runtests($list);
                fprintf('COUNTS %d %d %d\\n', sum([r.Passed]), sum([r.Failed]), sum([r.Incomplete]));
                for k = 1:numel(r), if r(k).Failed, fprintf('FAILED %s\\n', r(k).Name); end, end"
    fi
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

# ---- BUILD THE JOB LIST ------------------------------------------------------------------
# A job is either a whole suite or ONE TEST of a suite. Sharding only happens when running in
# parallel: serially it would just add MATLAB startups to the same total.
JOBS_LIST=()
for s in "${SUITES[@]}"; do
    shard=0
    if [ "$JOBS" -gt 1 ]; then
        for ss in "${SHARD_SUITES[@]}"; do [ "$s" = "$ss" ] && shard=1; done
    fi
    if [ "$shard" -eq 1 ]; then
        # Ask MATLAB for the test names once; if that fails for any reason, fall back to running
        # the suite whole rather than losing it.
        mapfile -t tests < <(matlab -batch "t = matlab.unittest.TestSuite.fromClass(meta.class.fromName('$s')); for k = 1:numel(t), fprintf('TEST %s\n', t(k).ProcedureName); end" 2>/dev/null | sed -n 's/^TEST //p')
        if [ "${#tests[@]}" -gt 0 ]; then
            for t in "${tests[@]}"; do JOBS_LIST+=("$s/$t"); done
        else
            printf 'suite.sh: could not enumerate %s, running it whole\n' "$s" >&2
            JOBS_LIST+=("$s")
        fi
    else
        JOBS_LIST+=("$s")
    fi
done

RESDIR=$(mktemp -d)
trap 'rm -rf "$RESDIR"' EXIT

# ---- RUN, AT MOST $JOBS AT A TIME --------------------------------------------------------
# THE THROTTLE IS `xargs -P`, NOT `wait -n`, and that is a measured decision rather than a style
# one. The first version of this used the textbook `job & ; if running >= JOBS; then wait -n; fi`
# loop. Under this bash it did not hold the cap: asked for 4, it had EIGHT MATLABs live and
# climbing (measured 2026-08-19 by process start time and CPU -- 8 sessions at ~1.7 GB each). The
# failure mode is silent and it is the expensive kind: it oversubscribes a shared machine while
# reporting the -j the user asked for. `xargs -P` does the same thing in one well-tested
# primitive, so the cap is not this script's invention to get wrong.
#
# Each worker re-invokes THIS script with --run-one, which is why the mode exists; it is not a
# public entry point.
if [ "$JOBS" -gt 1 ]; then
    printf 'running %d jobs, %d at a time
' "${#JOBS_LIST[@]}" "$JOBS"
    for idx in "${!JOBS_LIST[@]}"; do
        printf '%s	%s
' "$idx" "${JOBS_LIST[$idx]}"
    done | CCA2DIR="$CCA2DIR" CCA2_TEST_TIMEOUT="$TIMEOUT" CCA2_RESDIR="$RESDIR"            xargs -P "$JOBS" -d '
' -I{} bash "$0" --run-one {}
else
    for idx in "${!JOBS_LIST[@]}"; do
        run_job "${JOBS_LIST[$idx]}" "$idx"
    done
fi

# ---- REPORT, IN JOB ORDER ----------------------------------------------------------------
for idx in "${!JOBS_LIST[@]}"; do
    job="${JOBS_LIST[$idx]}"
    f="$RESDIR/$(printf '%04d' "$idx")"
    if [ ! -f "$f" ]; then
        printf 'SUITE %-40s NO RESULT\n' "$job"; bad=1; continue
    fi
    read -r _ st elapsed rest < "$f"
    case "$st" in
        timeout)
            printf 'SUITE %-40s TIMEOUT after %ss\n' "$job" "$TIMEOUT"
            tt=$((tt+1)); bad=1 ;;
        errored)
            printf 'SUITE %-40s ERRORED (no result; exit %s)\n' "$job" "${rest#* }"
            bad=1 ;;
        ok)
            read -r _ np nf ni <<<"$rest"
            printf 'SUITE %-40s pass=%3d fail=%3d incomplete=%3d  %5ds\n' \
                   "$job" "$np" "$nf" "$ni" "$elapsed"
            grep '^FAILED ' "$f" | sed 's/^FAILED /  FAIL /' || true
            tp=$((tp+np)); tf=$((tf+nf)); ti=$((ti+ni))
            [ "$nf" -eq 0 ] && [ "$ni" -eq 0 ] || bad=1 ;;
        *)
            printf 'SUITE %-40s UNREADABLE RESULT\n' "$job"; bad=1 ;;
    esac
done

wall=$(( $(date +%s) - run_started ))
printf 'TOTALS pass=%d fail=%d incomplete=%d timeout=%d over %d jobs in %ds (-j %d)\n' \
       "$tp" "$tf" "$ti" "$tt" "${#JOBS_LIST[@]}" "$wall" "$JOBS"
exit "$bad"
