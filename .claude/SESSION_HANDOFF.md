# Session Handoff

_Last updated: 2026-11-05T18:53-07:00_

## What happened this session
The overnight branch was merged to main. Every biconj result is now an AlgCon (EnvelopeMin joined
as a no-subdivision leaf). The bridge producer's 2^53 overflow refusals are gone: the
witnessOverflow, branchOverflow and bigZContent refusals now finish in BigInteger (0 hits in 2220
pairs), and polynomialOverflow was removed after Psi got a BigInteger path. Large data (entries
1e3..1e6) now builds faces that match the definition to <= 2.3e-13. DECISIONS 2026-11-05 (third) to (fifth).

## Where things stand
- Branch: `main` @ `ccf7282` -- "handoff: polynomialOverflow at large data fixed"
  (work commit 5b99d37, before that 973efcf)
- Pushed: no (main is ahead of origin/main; never pushed this session)
- Tests: fast 831/0, normal 35/0, slow 128/0 (+2 known-incomplete), verylong Ex2Subdivision 3/0, guards 33/0
- Known reds: none. Blocked: nothing.

## Next steps (proposed to the user, in this order; not yet chosen)
1. Test oracle: `tieMax`/`tiePoints` in bridgeFaceAlgTest.m sweep s1 in [-25,25] and keep |S| < 1e4,
   so they miss ties centred at |c0| ~ 1e3 (7 false disagreements at large data). Scale the window to
   the tie centre and size. A data-scaled version (`tieMaxW`) is in the scratch probe
   (scratchpad po/or/pocheck.m, may be gone).
2. Release/paper: state the MATLAB + Symbolic Math Toolbox dependency; cite `convexdb` by its Zenodo DOI.
3. Example 2 is an EnvelopeMin by default (overlay budget 3000); a full subdivision takes ~hours. Decide
   whether a faster subdivision is worth the effort.
4. Low value: two zero-area faces in the c = -2 AlgCrv; "two bridge faces owning one cell" refusal
   (never observed); the legacy 'cplq' engine's issues (parked).
- Still open but by design: bridgeFaceAlg:noWitness (the tie has no usable rational point).
- Decided: leave Example 2 pair 2's non-winning faces.

## Relevant files
- `bridgeFaceAlg.m` -- implicitPoly -> implicitPolyDouble / implicitPolyBig, discBigRaw, bigTail
  (big z-content division), imageConic overflow -> row omitted, isSingular
- `bridgePhi.m` -- 'big' mode: Phi big form; U and det(M_mu) as bigZ when past 2^53
- `bridgeWitness.m` -- memberExact / memberExactArcBig (BigInteger Sturm), detNonzero
- `ratQ.m` -- countRootsInHalfOpenBig, isolateRootsBig, signAtBig, bigGcdRow, bigDivLinear; bigRow reads bigZ
- `bridgeFaceAlgTest.m` -- theFormerOVERFLOWRefusalsNowBuildTheFace, aWITNESSQuartic...,
  LARGEDataBuildsTheFaceWherePsiItselfPasses2to53; oracle helpers at the bottom (tieMax etc.)
- `guardExerciseTest.m` -- guard inventory (remainder list)
- `DECISIONS.md` 2026-11-05 (third)-(fifth); `SUPPORT_MATRIX.md` 3.2; `TODO.md` 'Open now'
