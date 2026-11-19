#!/usr/bin/env bash
# Build the CCA2 manuscript and keep the project .bib in step with the master.
#
#   cd manuscript && ./build.sh
#
# Steps:
#   1. Keep manuscript/biblio.bib in step with the master
#      (../biblio/biblio.bib): append every cited entry that is in the master
#      but not yet in the project file, report cited keys found in neither.
#      (Runs BEFORE the build so the first build already has its .bib; it
#      reads the cited keys from main.aux when present, else from the \cite
#      commands in main.tex.)
#   2. latexmk (pdflatex + bibtex, as needed) -- the journal class (mpc.cls)
#      is not available in this environment, so this uses `article` (see the
#      note at the top of main.tex).
#   3. Report undefined references / citations, if any.
set -u
cd "$(dirname "$0")"

# The shared bibliography lives at <AI root>/biblio, next to the project.
# Under the AI root that is ../biblio; from a worktree elsewhere (e.g. a
# Kandev task folder) point AI_ROOT at the AI root explicitly.
AIROOT="${AI_ROOT:-..}"
BIBSUBSET="$AIROOT/biblio/tools/bibsubset.py"
if [ ! -f "$BIBSUBSET" ]; then
    echo "cannot find the master bibliography at $BIBSUBSET"
    echo "set AI_ROOT to the AI root (the folder holding biblio/) and re-run"
    exit 1
fi

echo "== 1/3 biblio in step with the master =="
python3 "$BIBSUBSET" main.tex || {
    echo "bibsubset reported a problem (cited keys missing from the master)"
    exit 1
}

echo "== 2/3 latexmk =="
latexmk -pdf -interaction=nonstopmode main.tex || {
    echo "latexmk failed; see main.log"
    exit 1
}

echo "== 3/3 undefined references / citations =="
if grep -E 'LaTeX Warning: (Reference|Citation) .* undefined' main.log; then
    echo "UNDEFINED REFERENCES OR CITATIONS -- see above"
    exit 1
fi
echo "clean: no undefined references or citations"
