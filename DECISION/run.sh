#!/usr/bin/env bash
# ==================================================================================================
# CCA2 language-decision benchmark -- runner
#   1) bench.py     -> results/cca2_py.json   (PYTHON side; runs where Python+sympy exist)
#   2) bench_jl.jl  -> results/cca2_jl.json   (JULIA side; runs ONLY where a Julia+Nemo toolchain exists,
#                                            otherwise marked PENDING -- we never fabricate a result)
#   3) decide.py    -> results/decision.json  (unions both, applies the G->A->B->F rule)
# ==================================================================================================
set -uo pipefail
cd "$(dirname "$0")"
PY=${python3:-python3}
RES=results
mkdir -p "$RES"
# best-effort: put juliaup's bin on PATH (default install location) if it exists
if [ -d "$HOME/.juliaup/bin" ] && ! command -v julia >/dev/null 2>&1; then
  export PATH="$HOME/.juliaup/bin:$PATH"
fi

echo "==============================================================================="
echo " CCA2 Python-vs-Julia language decision   (real runs in, honest 'pending' out)"
echo "==============================================================================="

# ---- 1) PYTHON side (real on this host) -------------------------------------
echo
echo ">>> [1/3] PYTHON side  (bench.py: sympy + mpmath)"
if command -v "$PY" >/dev/null 2>&1 && "$PY" -c 'import sympy, mpmath' >/dev/null 2>&1; then
  "$PY" bench.py
else
  echo "    ! python3/sympy/mpmath not available on this host -> Python side PENDING"
  printf '{"host":"python","status":"PENDING (no python/sympy on this host)"}\n' > "$RES/cca2_py.json"
fi

# ---- 2) JULIA side (only if a real toolchain exists) ------------------------
echo
echo ">>> [2/3] JULIA side  (bench_jl.jl: Nemo)"
if command -v julia >/dev/null 2>&1; then
  echo "    julia found: $(julia --version 2>/dev/null | head -1)"
  if julia bench_jl.jl; then
    :
  else
    echo "    ! julia benchmark raised -> recording as PENDING/absent (no fabricated result)"
    printf '{"host":"julia","status":"PENDING (julia present but bench_jl.jl did not run clean)"}' > "$RES/cca2_jl.json"
  fi
else
  echo "    ! no 'julia' binary on this host -> Julia side PENDING (toolchain absent, not a measured loss)"
  printf '{"host":"julia","status":"PENDING (no Julia/Nemo toolchain on this host)","backends":{}}' > "$RES/cca2_jl.json"
fi

# ---- 3) decision rule (unions both) -----------------------------------------
echo
echo ">>> [3/3] DECISION  (decide.py: G gate -> A coverage/cost -> B wash -> F tiebreak)"
"$PY" decide.py

echo
echo "artifacts:"
echo "    results/cca2_py.json   results/cca2_jl.json   results/decision.json"
echo "==============================================================================="
