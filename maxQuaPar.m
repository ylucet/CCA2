function g = maxQuaPar(g1, g2)
% maxQuaPar  Step 3 of the 'cplq' pipeline: pointwise maximum of two full-domain QuaPar objects.
%
% objective: h(s) = max(g1(s), g2(s)) for all s in R^2, returned as one QuaPar. This is the
%   missing primitive Step 3 needs: when a nonconvex triangle piece's original quadratic has more
%   than one "sub-piece" after Step 1 (convEnvCPLQ splits it), Step 2 (conjPieceCPLQ) conjugates
%   each sub-piece into its OWN full-domain QuaPar, and the true conjugate of the original piece
%   is the pointwise max of those. See DESIGN.md II.5.1 and conjCPLQ.m's conjSingleTriangle.
%
% [input]  g1, g2 : QuaPar, each finite everywhere (domain = R^2) -- i.e. every conjPieceCPLQ
%                    output. NOT a general "max of any two piecewise quadratics" utility: g1, g2
%                    must be conjugates of adjacent sub-pieces of the SAME originally-nonconvex
%                    domain (see the scoping caveat below).
% [output] g      : QuaPar, h = max(g1,g2), finite everywhere.
%
% SCOPING CAVEAT (load-bearing, not a generic fact): this only works because, for THIS pipeline's
%   adjacent sub-pieces, g1.f(k,:)-g2.f(l,:) is always a DEGENERATE conic (its full 3x3
%   discriminant Delta vanishes), so it factors into a parabola, a single line, or -- the case
%   this file's first version got wrong, see below -- a pair of two distinct straight lines, of
%   which only one is ever active inside any one comparison cell. That is a theorem about pieces
%   produced by convEnvCPLQ's own triangle-splitting, NOT a generic fact: two arbitrary PD
%   quadratics CAN have a genuinely irreducible (ellipse/hyperbola) difference (e.g. Hessians
%   diag(1,4) and diag(4,1) give a difference-Hessian diag(-3,3) with no affine correction that
%   makes Delta vanish), which QuaPar cannot represent as an edge. This is checked at runtime
%   (splitCell asserts Delta~=0) rather than assumed silently, so a violating input errors loudly
%   instead of producing a wrong/unrepresentable result.
%
%   HISTORY: an earlier version of this check tested only the quadratic PART's discriminant
%   delta=b^2-4ac (equivalently rank(Qd) for the 2x2 submatrix [a b/2;b/2 c]), which decides the
%   conic's TYPE (hyperbolic/parabolic/elliptic) but not its IRREDUCIBILITY. On the f(x,y)=xy,
%   three-convex-edge-triangle example (see /home/ylucet/CCA2/3-edge.tex), the two comparison
%   boundaries have delta>0 (so the old check rejected them as "not degenerate") but the full
%   discriminant Delta=0 exactly: both are pairs of straight lines, not hyperbolas. Fixed below to
%   test Delta directly, and to detect -- per split, from the two actual boundary-crossing points,
%   not from a symbolic factorization -- whether the connecting curve between them is straight
%   (the generic case for this pipeline, per the above) or a genuine parabola.
%
%   HISTORY (later session): with the Delta fix above in place, full end-to-end assembly of
%   maxQuaPar(g1,g2) on the same f(x,y)=xy example still failed with "no matching neighbour" on a
%   plain, cleanly-decided cell (g1 face 1 vs g2 face 4) nowhere near the hyperbola cells -- a
%   separate face-clipping topology gap, not a conic degeneracy issue. Four distinct bugs were
%   found and fixed, each uncovered only after fixing the previous one let the computation get
%   further:
%     1. clipPolyHalfPlane's bounded-polygon "keep the far/wrapped arc" branch computed
%        keepIdx = mod((p2):(p1-1), nv) + 1, which is ALWAYS empty in MATLAB (p1<p2 always here,
%        so p1-1<p2, and the colon operator does not wrap) -- silently turning "keep the far arc"
%        into "keep nothing" and collapsing the cell to empty. Fixed by adding +nv before modding:
%        mod((p2):(p1-1+nv), nv) + 1. The identical bug, identically fixed, was in splitCell's
%        analogous restIdx.
%     2. Once #1 actually returned vertices, they were wired up wrong: that same far-arc branch
%        built [X1; kept-vertices; X2], but X1 and X2 are always adjacent via the new cut edge
%        regardless of which arc survives, so for the far-arc case the correct order is
%        [X2; kept-vertices; X1] -- the old order produced a self-intersecting "bowtie" cell.
%     3. assemblePieces encodes both rays of an unbounded piece apex-first in Ep (matching
%        QuaPar's own E-matrix convention: "column 1 is always the finite apex" -- see facePoly's
%        header comment), which, unlike segments, is NOT walk-order. So two adjacent pieces
%        sharing one physical ray both encode it as the SAME (a,b) pair, not swapped -- but the
%        half-edge matching loop was uniformly searching for a swapped pair (correct for segments,
%        wrong for rays), so no ray ever found its neighbour. Fixed by branching the search: rays
%        match on identical (a,b), segments on swapped (a,b).
%     4. g2 face 1's own three real vertices happen to be exactly collinear, so its two
%        consecutive boundary edges (to face 2, then to face 3) clip g1 face 2 by the SAME
%        half-plane twice -- the second clip is a geometric no-op, so the shared vertex between
%        those two collinear edges (where the true neighbouring face changes from face 2 to
%        face 3) never becomes an explicit vertex of the clipped cell, leaving one straight cell
%        edge that silently spans two different neighbours. Fixed by adding
%        insertPassthroughVertices, called at the end of clipByFace, which re-inserts any
%        polyK/polyL vertex lying in the open interior of one of the clipped cell's edges.
%   After all four fixes, maxQuaPar(g1,g2) on this example fully assembles and matches ground
%   truth at 6 of 7 sample points; the 7th (s=(-3,2)) was wrong due to two further, separate bugs,
%   both now fixed (later session):
%     5. QuaPar.m's orderEdges (NOT in this file) recomputed its boundary-walk pivot vertex from
%        scratch every iteration via a left/right rule based only on the current edge's own
%        orientation; that rule doesn't know which vertex the walk just arrived at, and for a
%        segment edge entered from a ray, could pick the vertex the walk came FROM instead of the
%        one it arrived AT, duplicating one boundary edge and dropping its true neighbour. Fixed
%        by tracking the pivot vertex as the previous iteration's iNext instead of recomputing it.
%     6. This file's assemblePieces assigned ray edges' F(:,1)/F(:,2) (left/right face) by
%        processing order ([HE(h).piece, HE(opp).piece], whichever half-edge was enumerated
%        first) -- fine for segments, where the CCW walk direction is naturally reversed between
%        the two adjacent pieces and so encodes which is on the left, but WRONG for rays, which
%        both adjacent pieces encode identically (apex-first), carrying no left/right information
%        via processing order. This let one piece's cell silently claim territory that geometrically
%        belonged to its neighbour across a ray boundary. Fixed by deriving left/right from which
%        end of its own CCW boundary each piece uses for that ray: the piece for which the ray is
%        OUTGOING (walked apex->direction, matching the stored a->b order) is on the left, same as
%        segments; the piece for which it is INCOMING (walked direction->apex, i.e. b->a) is on the
%        right.
%   After fixes 5 and 6, maxQuaPar(g1,g2) on this example matches ground truth at all 7 sample
%   points, including s=(-3,2) -- see maxQuaParTest.m's
%   maxQuaParResolvesBothHyperbolaCellsWithoutMisclassifying.
%
%   HISTORY (later session): splitCell's degeneracy guard rejected ANY pair whose full 3x3
%   discriminant Delta was nonzero, mislabelling the rejection "a genuine ellipse/hyperbola" --
%   but Delta~=0 with delta=b^2-4ac==0 (parabolic TYPE) is a genuine, representable PARABOLA, not
%   a hyperbola/ellipse; the guard only needs to reject delta~=0 (hyperbolic/elliptic type) cases
%   whose Delta is also nonzero. Found via a randomized stress test across many f(x,y)=xy
%   triangles (not just the one hard-coded 3-edge.tex example), e.g.
%   T=(0,0),(7.02,0.67),(8.43,7.63): a legitimate g1-face/g2-face boundary there has delta~0 but
%   Delta~0.063, and maxQuaPar(g1,g2) crashed on it even though the rest of splitCell (the
%   isStraight/edgeEc construction) already builds curved parabolic edges correctly once let
%   through. Fixed by conditioning the guard on delta as well as Delta. See
%   maxQuaParTest.splitCellAcceptsGenuineNonDegenerateParabola.
%
%   HISTORY (later session): assemblePieces used to merge ALL pieces' vertices into one global
%   list via a SINGLE coordinate-distance tolerance, then matched half-edges by global vertex
%   INDEX equality -- unsound for near-degenerate inputs, since a genuine cross-arithmetic-noise
%   gap (~1e-5, two different formulas computing "the same" cell corner) and a genuine distinct
%   nearby-vertex gap (as small as ~1e-5 too, for a tight fan of many pieces meeting near one
%   point) can overlap, so no single tolerance can separate them: loosening it enough to fix one
%   "no matching neighbour" crash could make some piece merge two of its OWN vertices instead,
%   causing a DIFFERENT crash (QuaPar.orderEdges rejecting a self-touching face). A prior session
%   tried exactly that blanket widen, found it traded one crash for MORE crashes overall on a
%   broader stress test, and reverted it (see the session handoff this was picked up from). Fixed
%   properly by never relying on coordinate-distance vertex identity at all: half-edges are now
%   matched directly by GEOMETRY (matchHalfEdges), comparing pairs of DIFFERENT pieces' edges (a
%   segment's two endpoints, or a ray's apex+direction), so two vertices belonging to the SAME
%   piece are never compared against each other. Global vertex identity (buildGlobalVertices) is
%   then derived via union-find driven purely by those confirmed half-edge matches, so it can
%   never unify two of one piece's own vertices either. See assemblePieces' own HISTORY comment
%   for the full derivation. Verified on the exact reproduction case from the handoff
%   (T=(6.0365,4.9504),(9.8960,6.3015),(1.4908,3.3753), see
%   maxQuaParTest.assemblePiecesResolvesNearDuplicateApexCluster) and via a ~5000-triangle
%   randomized stress test: the maxQuaPar:internal crash rate dropped from ~4/800 valid samples to
%   ~1/800 (the residual case is a genuinely ambiguous 3-way vertex cluster ~2e-4 apart, which
%   needs vertex PROVENANCE -- tracking which original g1/g2 face boundary each vertex came from --
%   to resolve correctly; still open), with ZERO new wrong-answer regressions (every "wrong
%   answer" case found by the stress test was independently confirmed to reproduce identically on
%   the unmodified pre-fix code, i.e. pre-existing and unrelated to this change).
%
% STATUS (incremental implementation -- see DESIGN.md II.5.1 and the session plan):
%   * IMPLEMENTED: g1, g2 purely polyhedral (all-zero Ec, i.e. every edge a line/ray/segment) --
%     exactly what conjPSDRank1QuadTriangle/conjPSDRank1QuadTriangleTie/conjLinearTriangle produce.
%     The face FUNCTIONS may still be genuinely quadratic; only the DOMAIN boundaries must be
%     polyhedral for this version.
%   * IMPLEMENTED (this session): exactly ONE of g1/g2 carrying curved (parabolic) input edges --
%     conjBilinearXYoneCE / conjIndefiniteQuadTriangle-with-1-convex-edge output. Since max is
%     symmetric, the two operands are SWAPPED below if needed so the curved one is always g1: the
%     (k,l) loop only ever clips polyK (a g1 face) against polyL's (a g2 face's) half-planes, so
%     with that normalization the only new geometric primitive needed is arc-vs-half-plane
%     clipping (clipArcByHalfPlane.m), never polygon-vs-conic-side clipping. A curved cell edge is
%     threaded through clipByFace by clipPolyHalfPlaneCurved (the polyhedral path,
%     clipPolyHalfPlaneStraight, is left untouched -- see its own dense HISTORY below for why that
%     matters).
%   * IMPLEMENTED (later session): splitting a cell that ALREADY carries an arc -- see splitCell.
%     It needs NO general conic-conic solver, and the older wording here (and on the guards) that
%     said it did was wrong. Every curved edge in play is a PARABOLA (QuaPar.assertParabolic
%     rejects b^2-4ac~=0 outright; the irreducible ellipse/hyperbola case is unreachable, see the
%     scoping caveat above), so restricting the splitting conic to the arc via
%     parabolaArcFrame.conicCoeffs gives one univariate QUARTIC in the frame's global monotone
%     parameter u, and roots() finishes it.
%     ONE ARC PER FACE remains the invariant -- QuaPar stores one conic per edge -- and it is
%     maintained by SUBDIVIDING with a straight chord (splitTwoArcPiece), never by giving a piece a
%     second curve slot. Measured on the named fixture plus a 395-quadrilateral sweep: of 22
%     curved-cell splits, the splitting curve crossed two STRAIGHT edges every time and left the
%     existing arc either untouched (19) or tangent (3), never cut -- the C1 tangency structure
%     this file already documents. Assembled results went 58 -> 76 of 395, every one of the 77
%     (sweep + fixture) exact against the closed-form sup (worst 2.8e-14) and violation-free under
%     maxQuaParTest.arrangementViolations.
%   * TODO: BOTH inputs curved -- not implemented, errors clearly. Also not implemented, and
%     likewise erroring clearly rather than guessing: a clip line that cuts one arc TWICE (the arc
%     bulging across it between its own two endpoints -- the result is either disconnected or has
%     two separate arcs, neither representable as one QuaPar face; see clipPolyHalfPlaneCurved);
%     the splitting curve genuinely CROSSING a cell's arc (splitCell checks for this rather than
%     assuming the tangency structure); and splitting an UNBOUNDED cell that carries an arc. None
%     of the three was produced by the sweep.
%     Residual, and NOT caused by the split-a-curved-cell work: 3 of 395 now reach assemblePieces
%     and fail its half-edge matching ('no matching neighbour'). That is the pre-existing ambiguous
%     vertex-cluster limitation documented at assemblePieces (~1/800 on the polyhedral stress
%     test); those 3 errored before too, just earlier and under a different identifier. One of them
%     fails with no subdivision performed at all, which is what rules the new code out as a cause.
%
%   VALIDATION of the curved case (randomized sweep, same generator as the counts above; ground
%   truth = the closed-form sup of the bilinear objective over the quadrilateral):
%     - Of 115 sampled splits: 85 assembled, 30 hit the split-a-curved-cell guard, and ZERO hit
%       either the arc-bulge guard or maxQuaPar:internal. (The arc-bulge case is therefore a
%       defensive guard, not an observed configuration -- the sweep never produced one, consistent
%       with the tangency structure described below.)
%     - Those 85 results were checked at every result VERTEX, every straight-edge MIDPOINT, and 60
%       generic interior points each. Exact (to ~1e-14) at every edge midpoint (340/340) and every
%       interior point (5100/5100).
%     - CAVEAT ON EVERY NUMBER IN THIS BLOCK: the sweep that produced them was a throwaway
%       script, not committed, and its seed was not recorded, so none of these figures can
%       be re-derived. Re-running any of this means writing a NEW seeded, committed sweep --
%       see SUPPORT_MATRIX.md section 0.1, which now requires that of any quoted measurement.
%     - THAT SWEEP NOW EXISTS: sweepMaxQuaParCurvedSplit.m, seeded and committed. Prefer its
%       numbers to the ones above, which are kept only as a record of what was once claimed.
%       sweepMaxQuaParCurvedSplit(20260802, 200) gives: 131 sampled splits, 30 assembled,
%       112/112 straight-edge midpoints and 1800/1800 interior points exact, and 0 of 369
%       result VERTICES disagreeing (the ~0.8% below no longer reproduces -- the QuaPar.eval
%       tolerance appears to have closed the curved half as well as the polyhedral one).
%       It also breaks the non-assembling cases out by error identifier, which the old count
%       did not: 80 of the 131 never reach maxQuaPar at all, because Step 2 refuses the
%       triangle with conjPieceCPLQ:notImplemented, and only 21 hit maxQuaPar's own guard. So
%       "85 of 115 assembled" was measuring a pre-filtered population, and the dominant
%       obstacle to this configuration is UPSTREAM of this file.
%     - Exactly AT a result vertex, ~0.8% (9/1105) disagree -- QuaPar.eval's exact, no-tolerance
%       point location can leave a corner unclaimed (+Inf) or claimed by an adjacent face. This is
%       PRE-EXISTING, not introduced here: the same sweep restricted to purely POLYHEDRAL pairs
%       (untouched code path) disagrees at ~1.4% (5/356) of result vertices, i.e. slightly more
%       often, while likewise being exact at every edge midpoint and interior point. Probing rings
%       of radius as small as 1e-8 around every disagreeing vertex gives the right answer to
%       ~1e-15, confirming the geometry is right and only the exactly-at-a-corner tie-break is at
%       issue. Tracked separately; see maxQuaParTest.curvedSamplePoints.
%     - ARRANGEMENT VALIDITY was checked directly on the assembled results, independently of this
%       file's own piece bookkeeping (curveAfter, half-edge matching, F) -- using only the final
%       V/E/Ec geometry: every pair of distinct edges must meet in the empty set or in a proper
%       face, i.e. only at a shared endpoint, with arcs sampled along the parabola itself rather
%       than along their chords. Zero violations across all 109 results the sweep assembled. The
%       INPUT-side counterpart of the same invariant -- a face vertex landing in the open interior
%       of an arc, which would force that arc to be split into two sub-arcs -- likewise never
%       occurred (0 occurrences), and is now DETECTED rather than silently ignored; see
%       insertPassthroughVertices. Regression test: maxQuaParTest.maxQuaParResultsAreValidArrangements.
%     - A clip that cuts an arc STRICTLY BETWEEN its endpoints never arose in the supported regime:
%       there, the other operand's face boundaries are tangent lines to the parabola (a conjugate
%       is C1 where its pieces join), so a clip either keeps the arc whole or reduces it to a
%       point. Such cuts DO occur for unrelated (non-adjacent) triangle pairs, but those violate
%       this file's scoping caveat and are already wrong on the pre-existing polyhedral path too
%       (14 of 15 such polyhedral pairs give wrong answers), so they are not usable as evidence
%       either way. That branch of clipPolyHalfPlaneCurved is therefore currently covered only by
%       clipArcByHalfPlaneTest's own direct 'cut' tests, not end to end.

    curved1 = hasCurvedEdge(g1); curved2 = hasCurvedEdge(g2);
    if false %#ok<UNRCH> -- kept for the reasoning it records; see the note inside
        % NOT "because conic-conic intersection is hard" -- splitCell now does exactly that
        % intersection, in closed form, via parabolaArcFrame.conicCoeffs. What is missing here is
        % the clipping side: clipByFace only ever clips a g1 face against g2's HALF-PLANES, so two
        % curved operands would need arc-vs-arc clipping, which clipPolyHalfPlaneCurved has no
        % case for.
        %
        % THAT REASONING WAS TOO COARSE, and the refusal is gone. It is a statement about a FACE
        % PAIR, not about the operands. Face intersection is SYMMETRIC, so a pair in which only
        % one of the two faces carries an arc can always be clipped in the direction that puts
        % the curved face first -- clipByFace performs that swap itself. Only the pair in which
        % BOTH faces are curved genuinely needs conic-vs-conic clipping, and that is now detected
        % per pair, which is a far smaller claim than refusing the whole operand. A QuaPar
        % conjugate carries its arcs on a minority of its faces, so most pairs never reach it.
    end
    if curved2 && ~curved1
        tmp = g1; g1 = g2; g2 = tmp;   % max is symmetric; prefer the curved operand as g1
    end
    assertCurvedEdgesAreArcs(g1);
    assertCurvedEdgesAreArcs(g2);
    assertFullDomain(g1, 'g1');
    assertFullDomain(g2, 'g2');

    pieces = struct('V', {}, 'dirIn', {}, 'dirOut', {}, 'dirInSign', {}, 'dirOutSign', {}, ...
        'curveAfter', {}, 'curveEc', {}, 'f', {});
    for k = 1:g1.nf
        polyK = facePoly(g1, k);
        for l = 1:g2.nf
            polyL = facePoly(g2, l);
            cells = clipByFace(polyK, polyL);   % a LIST: an arc-vs-arc clip can yield two cells
            f1row = g1.f(k,:); f2row = g2.f(l,:);
            for ci = 1:numel(cells)
                cell = cells{ci};
                if isempty(cell), continue; end
                [decided, winRow] = decideWinner(cell, f1row, f2row);
                if decided
                    cell.f = winRow;   % cell carries its own curveAfter/curveEc (0/[] if none)
                    pieces(end+1) = cell; %#ok<AGROW>
                    continue
                end
                newCells = splitCell(cell, f1row, f2row);
                for z = 1:numel(newCells)
                    pieces(end+1) = newCells(z); %#ok<AGROW>
                end
            end
        end
    end
    pieces = dedupPieces(pieces);
    pieces = dropSubsumedPieces(pieces);
    g = assemblePieces(pieces);
end

% ============================================================================================
% ----- collapsing duplicate cells produced by two different (k,l) pairs ----------------------
function pieces = dedupPieces(pieces)
% Two DIFFERENT (g1-face k, g2-face l) pairs can produce GEOMETRICALLY IDENTICAL cells: this is not
% a rare fluke but a structurally expected occurrence when g1 and g2 are conjugates of two
% sub-pieces glued along a shared PRIMAL seam (see splitTwoConvexEdges/convEnvCPLQ.m) -- a point
% common to both primal pieces (e.g. a seam endpoint) puts a genuinely shared feature into BOTH
% g1's and g2's own face structure there, and several of g2's OWN faces can meet at the dual image
% of that point (a "fan"), any one of which a touching g1 face can appear to align with exactly.
% This differs from convEnvCPLQ's splitThreeConvex case (always a C1-smooth seam), where no such
% dual-side fan alignment was ever observed -- it is a genuinely new configuration this session's
% fix to convEnvCPLQ (COAP Appendix A.4) exposes, common enough (about half of random cases, per
% this session's stress test) to need handling here rather than at the primal level.
%
% Left unhandled, assemblePieces sees each duplicate's boundary edges compete for the SAME single
% neighbour, orphaning one copy's edges (maxQuaPar:internal, "no matching neighbour").
%
% Fix: collapse each group of identical-geometry pieces (same vertex count, same real vertices as
% a SET, same dirIn/dirOut, no curved edge) into ONE. If every member of a group agrees on `f`,
% simply drop the duplicates. If they DISAGREE (each drawn from a different, only PARTIALLY
% applicable candidate g2/g1 face -- see header), reconcile by evaluating every candidate row at a
% point verifiably interior to the shared cell and keeping the largest: both g1's and g2's own rows
% are independently exact everywhere (conjPieceCPLQ's per-piece output, verified against ground
% truth to machine precision), so whichever row is larger AT an interior point of this cell is
% provably the correct one for the WHOLE cell (decideWinner's own premise -- a single row wins
% across an entire convex cell -- already relies on exactly this).
    n = numel(pieces);
    if n == 0, return; end
    groupOf = zeros(1,n);
    nextGroup = 0;
    for i = 1:n
        if groupOf(i) ~= 0, continue; end
        nextGroup = nextGroup + 1;
        groupOf(i) = nextGroup;
        for j = i+1:n
            if groupOf(j) == 0 && samePieceGeometry(pieces(i), pieces(j))
                groupOf(j) = nextGroup;
            end
        end
    end
    keep = false(1,n);
    for gi = 1:nextGroup
        idx = find(groupOf == gi);
        keep(idx(1)) = true;
        if numel(idx) == 1, continue; end
        pt = interiorSample(pieces(idx(1)));
        best = -inf; bestF = pieces(idx(1)).f;
        for jj = idx
            val = QuaPar.evalPoly(pieces(jj).f, pt);
            if val > best, best = val; bestF = pieces(jj).f; end
        end
        pieces(idx(1)).f = bestF;
    end
    pieces = pieces(keep);
end

function tf = samePieceGeometry(a, b)
    tol = 1e-6;
    tf = false;
    if a.curveAfter ~= 0 || b.curveAfter ~= 0, return; end
    if isempty(a.dirIn) ~= isempty(b.dirIn), return; end
    if size(a.V,1) ~= size(b.V,1), return; end
    if ~isempty(a.dirIn)
        if norm(a.dirIn/norm(a.dirIn) - b.dirIn/norm(b.dirIn)) > tol, return; end
        if norm(a.dirOut/norm(a.dirOut) - b.dirOut/norm(b.dirOut)) > tol, return; end
    end
    Va = sortrows(a.V); Vb = sortrows(b.V);
    tf = all(abs(Va(:) - Vb(:)) < tol);
end

function pt = interiorSample(piece)
% A point verifiably interior to `piece` (not on its own boundary), used to compare candidate
% winner rows -- see dedupPieces. For an unbounded piece with >=2 real vertices, stepping along the
% (shared, since only this configuration is ever deduplicated as unbounded -- see samePieceGeometry
% requiring matching dirIn/dirOut) ray direction from the midpoint of the real vertices stays
% equidistant from both bounding rays, hence interior regardless of which side is "left"; for
% exactly 1 real vertex (a cone), nudge by the sum of the two ray directions (the angle bisector's
% direction, up to scale); for a bounded piece (>=3 vertices), the centroid is always interior.
    nv = size(piece.V,1);
    if isempty(piece.dirIn)
        pt = mean(piece.V,1);
    elseif nv >= 2
        pt = mean(piece.V,1) + piece.dirIn;
    else
        pt = piece.V(1,:) + piece.dirIn + piece.dirOut;
    end
end

% ============================================================================================
% ----- dropping pieces geometrically SUBSUMED by another piece with the SAME winner -----------
function pieces = dropSubsumedPieces(pieces)
% Two DIFFERENT (k,l) pairs can produce two pieces that agree on the winning row `f` (same
% quadratic/affine formula) yet are NOT geometrically identical (so dedupPieces, which only
% collapses exact-geometry groups, leaves both) yet still OVERLAP -- specifically, one can be a
% strict subset of the other. This happens when a g1 (or g2) face k is cut into several (k,l)
% sub-pieces by different opposing faces l, l': every such sub-piece independently inherits face
% k's own boundary edges/rays as part of ITS OWN boundary, so two sub-pieces that both end up
% deciding "f1row (face k) wins" can each be a geometrically valid but differently-sized portion
% of face k's own territory -- and if one is wholly contained in the other, the smaller one is
% pure redundant duplication of already-covered territory (both, individually, are provably
% correct: every point in either does have `f` as its true winning row), not a distinct region of
% the final partition. Left in, the SMALLER one's boundary competes for neighbours that the
% LARGER one's containing edges already legitimately claim, starving genuine neighbours elsewhere
% of a match (maxQuaPar:internal in assemblePieces) -- see matchHalfEdges' oppositeSides for the
% sibling bug this complements (that one rejects a spurious PAIRING between two such pieces; this
% one removes the redundant piece outright when one wholly contains the other, which oppositeSides
% alone cannot fix since rejecting the pairing merely leaves both orphaned).
%
% Fix: for every pair of pieces sharing the same winning row (within tolerance), drop whichever is
% a subset of the other (by testing every real vertex, plus -- for an unbounded piece -- both
% recession directions, against the other's own half-plane/ray constraints via polyConstraints).
% Curved (parabolic) pieces are left untouched (curveAfter~=0 skipped): the source of this specific
% redundancy is confined to plain polyhedral sub-pieces of one polyhedral parent face.
    n = numel(pieces);
    if n <= 1, return; end
    tolF = 1e-6;
    drop = false(1,n);
    for i = 1:n
        if drop(i) || pieces(i).curveAfter ~= 0, continue; end
        for j = 1:n
            if i == j || drop(j) || pieces(j).curveAfter ~= 0, continue; end
            if norm(pieces(i).f - pieces(j).f) > tolF*(1+norm(pieces(j).f)), continue; end
            if isSubsumed(pieces(i), pieces(j))
                drop(i) = true;
                break
            end
        end
    end
    pieces = pieces(~drop);
end

function tf = isSubsumed(a, b)
% True iff a's whole region is contained in b's (same convention/shape as facePoly's poly struct:
% .V real CCW vertices, .dirIn/.dirOut ray directions or [] if bounded).
    tf = false;
    if ~isempty(a.dirIn) && isempty(b.dirIn), return; end   % unbounded a can't fit in bounded b
    cons = polyConstraints(b);
    scale = 1 + max(abs(b.V(:)));
    tol = 1e-6 * scale;
    val = a.V * cons(:,1:2)' - cons(:,3)';
    if any(val(:) > tol), return; end
    if ~isempty(a.dirIn)
        % Homogeneous (direction-only) check: no positional scale, just a small absolute tolerance.
        if any(cons(:,1:2)*a.dirIn' > 1e-6) || any(cons(:,1:2)*a.dirOut' > 1e-6)
            return
        end
    end
    tf = true;
end

function tf = hasCurvedEdge(g)
    tf = ~isempty(g.Ec) && any(g.Ec(:) ~= 0);
end

function tf = pieceIsCurved(p)
% Whether a poly/cell/piece actually carries a PARABOLIC edge. Deliberately NOT just
% `p.curveAfter ~= 0`: boundedPiece tags every piece it builds with curveAfter = "the closing
% edge", including the isStraight case where curveEc is all zeros, and two long-standing
% behaviours depend on that tag staying set (samePieceGeometry and dropSubsumedPieces both use
% curveAfter~=0 to opt splitCell's pieces out of dedup/subsumption entirely). So the tag alone
% means "this piece came from a split", and only a nonzero conic means "this edge is curved".
    tf = p.curveAfter ~= 0 && ~isempty(p.curveEc) && any(p.curveEc ~= 0);
end

function assertCurvedEdgesAreArcs(g)
% A curved edge must be a bounded ARC (E(j,3)==1). QuaPar has no representation for an unbounded
% curved ray (see QuaPar.m's header and splitCell's own note), and every downstream step here --
% facePoly's CCW curve index, clipPolyHalfPlaneCurved's arc clipping, localEdgeLists' Ec row --
% assumes the curve sits on a segment between two finite vertices.
    if isempty(g.Ec), return; end
    for j = 1:size(g.Ec,1)
        if any(g.Ec(j,:) ~= 0) && g.E(j,3) == 0
            error('maxQuaPar:notImplemented', ...
                'edge %d is a curved RAY; QuaPar curved edges must be bounded arcs.', j);
        end
    end
end

function assertFullDomain(g, name)
    if g.nv == 0 && g.nf == 1, return; end   % bare full-domain quadratic
    if any(g.F(:) == 0)
        error('maxQuaPar:notFullDomain', '%s is not finite everywhere (F has a 0 entry).', name);
    end
end

% ============================================================================================
% ----- extracting a face's boundary as {V (CCW finite vertices), dirIn, dirOut} --------------
function poly = facePoly(obj, k)
% Boundary of face k as an ORDERED list of finite vertices (poly.V, CCW, walked via obj.P{k},
% which orderEdges documents as clockwise -- so we REVERSE it here to get CCW), plus, for an
% unbounded face, the two ray directions poly.dirIn (ray ending at V(1,:), arriving from infinity)
% and poly.dirOut (ray leaving from V(end,:) to infinity); both empty for a bounded face.
% A ray edge in QuaPar always stores its finite apex in column 1 of E and a (non-unit) direction
% point in column 2, REGARDLESS of the sign it carries in P -- the sign only says whether, walking
% the boundary, this ray is traversed apex->direction (outgoing, last in the chain) or
% direction->apex (incoming, first in the chain); by the "0 or 2 rays" invariant (orderEdges),
% when rays are present they are always exactly the first and last chain elements.
%
% Vertex extraction walks the chain and, for each element, appends the vertex the walk ARRIVES AT
% (its "trailing" vertex): a forward/reversed segment arrives at V(b)/V(a) respectively; the
% incoming ray (t==1) arrives at its apex; the outgoing ray (t==n) arrives nowhere new -- its apex
% was ALWAYS already appended, either as the trailing vertex of the preceding segment (n>2), or by
% t==1's own append (n==2, a pure two-ray wedge sharing one apex, no segments) -- so it is never
% re-appended, regardless of n. See this session's inspect_P.m/test_facePoly.m investigation
% (verified against real conjPSDRank1QuadTriangle output by hand; an earlier draft of this
% function wrongly re-appended the apex when n>2, producing a duplicate vertex).
%
% CURVED EDGE: at most one boundary edge of face k may carry a nonzero conic Ec(j,:) (a parabola
% arc). It is reported as poly.curveAfter = the index i of the CCW boundary edge poly.V(i) ->
% poly.V(i+1) that IS the arc (wrapping to poly.V(nv)->poly.V(1) when the face is bounded), 0 when
% the face is purely polyhedral -- the SAME "vertex index the curve starts at" convention
% splitCell/boundedPiece already use for the pieces they build. Mapping the walk position t of the
% curved element to that CCW index: the CW walk's element t is the edge V_cw(t-1)->V_cw(t) and
% flipud reverses the list, giving i = n+1-t for a bounded face (n elements, nv=n vertices; t==1 is
% the closing edge, i=n) and i = n-t for an unbounded one (n elements, nv=n-1 vertices, since the
% two rays contribute one and zero vertices respectively -- see above). Its conic is returned as
% poly.curveEc, SIGN-NORMALIZED so that evalConic(poly.curveEc,.) > 0 on face k's own interior:
% QuaPar's orientation invariant already guarantees Ec(j,:) > 0 on the LEFT of V(E(j,1))->V(E(j,2)),
% i.e. on face F(j,1), so this is an exact index test (F(j,1)==k ? keep : flip), with no need for
% an interior sample point. clipPolyHalfPlaneCurved relies on that convention to know which side of
% the arc is inside.
    if obj.nv == 0   % bare full-domain quadratic: the "face" IS all of R^2, no boundary at all
        poly.V = zeros(0,2); poly.dirIn = []; poly.dirOut = [];
        poly.dirInSign = []; poly.dirOutSign = [];
        poly.curveAfter = 0; poly.curveEc = [];
        return
    end
    Pk = obj.P{k};
    n = numel(Pk);
    V = zeros(0,2); dirIn = []; dirOut = []; dirInSign = []; dirOutSign = [];
    curveT = 0; curveEc = [];
    for t = 1:n
        j = abs(Pk(t)); s = sign(Pk(t));
        a = obj.E(j,1); b = obj.E(j,2);
        if ~isempty(obj.Ec) && any(obj.Ec(j,:) ~= 0)
            if curveT ~= 0
                error('maxQuaPar:notImplemented', ...
                    ['face %d has more than one curved edge; a piece here carries at most one ' ...
                     '(QuaPar''s own per-piece Ec slot).'], k);
            end
            curveT = t;
            if obj.F(j,1) == k, curveEc = obj.Ec(j,:); else, curveEc = -obj.Ec(j,:); end
        end
        if obj.E(j,3) == 0   % ray: column 1 is always the finite apex
            apex = obj.V(a,:); d = obj.V(b,:) - apex; d = d/norm(d);
            if t == 1
                V(end+1,:) = apex; dirIn = d; dirInSign = s; %#ok<AGROW>
            else   % t == n, the only other ray: apex already appended, never re-add it
                dirOut = d; dirOutSign = s;
            end
        else
            if s > 0, V(end+1,:) = obj.V(b,:); else, V(end+1,:) = obj.V(a,:); end %#ok<AGROW>
        end
    end
    % obj.P{k} is CLOCKWISE; reverse to CCW, which swaps the roles of "first ray" (incoming when
    % reversed becomes outgoing) and "last ray" -- the ray direction VECTORS themselves (apex to
    % infinity) don't change, only which end of the (now-reversed) vertex list they attach to. The
    % ray's own sign(Pk(t)) (dirInSign/dirOutSign) is an intrinsic property of THAT specific P{k}
    % entry (not of its role), so it is carried along with the direction it was computed from,
    % swapped identically -- see polyConstraints' HISTORY for why this is needed.
    poly.V = flipud(V);
    poly.dirIn = dirOut;    % swapped, see comment above
    poly.dirOut = dirIn;
    poly.dirInSign = dirOutSign;
    poly.dirOutSign = dirInSign;
    if curveT == 0
        poly.curveAfter = 0; poly.curveEc = [];
    elseif isempty(poly.dirIn)
        poly.curveAfter = n + 1 - curveT;   % bounded: nv = n (see the CURVED EDGE note above)
        poly.curveEc = curveEc;
    else
        poly.curveAfter = n - curveT;       % unbounded: nv = n-1
        poly.curveEc = curveEc;
    end
end

% ============================================================================================
% ----- clipping one convex poly by every boundary constraint of another (face intersection) --
function cells = clipByFace(polyK, polyL)
% polyK intersected with the convex region bounded by polyL. Returns a CELL ARRAY of clipped
% cells -- usually one, but two when the survivor carries two arcs and must be subdivided, which
% is the invariant QuaPar's one-Ec-per-edge representation forces. Empty when nothing survives.
%
% polyConstraints turns polyL's boundary into HALF-PLANES, which a parabolic edge is not, so the
% arc is handled separately at the end:
%   * if only ONE of the two faces carries an arc, SWAP so the curved one is polyK -- face
%     intersection is symmetric, so polyK n polyL is the same set either way, and only the two
%     operands' f rows distinguish them, which the caller passes separately and this never
%     touches;
%   * if BOTH carry one, polyL's arc becomes a CURVED cut applied after the straight ones by
%     clipPolyByConic. Deferring it is deliberate: a cell the straight clips already empty needs
%     no curved surgery, and the curved step gets the smallest cell available.
    if polyL.curveAfter ~= 0 && polyK.curveAfter == 0
        tmp = polyK; polyK = polyL; polyL = tmp;
    end
    cutConic = []; cutX0 = []; cutX1 = [];
    if polyL.curveAfter ~= 0
        cutConic = polyL.curveEc;
        % polyL's boundary is the bounded ARC of that conic, not the whole conic. Its two
        % endpoints are carried through so the cut can be restricted to the arc's own span --
        % without that restriction the extended conic cuts the cell where polyL's boundary does
        % not, which shows up as a single unpaired crossing.
        nvL = size(polyL.V,1);
        cutX0 = polyL.V(polyL.curveAfter,:);
        cutX1 = polyL.V(mod(polyL.curveAfter, nvL)+1,:);
    end

    cell = polyK;
    cons = polyConstraints(polyL);
    for i = 1:size(cons,1)
        cell = clipPolyHalfPlane(cell, cons(i,1:2), cons(i,3));
        if isempty(cell), cells = {}; return; end
    end
    % Two consecutive real vertices of polyL (or polyK) can be exactly collinear with a third --
    % e.g. this pipeline's own faces sometimes have 3 real vertices on one straight line -- in
    % which case polyConstraints emits the SAME half-plane twice (once per collinear edge) and the
    % second clip is a geometric no-op. The vertex BETWEEN those two collinear edges is where the
    % true neighbouring face changes (e.g. from g2's face2 to face3), but it never becomes a vertex
    % of `cell` since no half-plane clip actually cut there -- leaving a straight cell edge that
    % silently spans TWO different neighbours, so assemblePieces can never find a match for either
    % sub-portion (see maxQuaPar.m header HISTORY). Explicitly re-insert any polyK/polyL vertex
    % that lies in the open interior of one of cell's own edges to restore that missing corner.
    cell = insertPassthroughVertices(cell, [polyK.V; polyL.V]);
    % A BOUNDED cell normally needs 3 vertices to have positive area -- but a bounded cell with a
    % curved edge can legitimately be a two-vertex "lens" (one parabola arc plus one straight
    % chord/cut edge between the same two endpoints), which does.
    minV = 3;
    if pieceIsCurved(cell), minV = 2; end
    if size(cell.V,1) < 1 || (isempty(cell.dirIn) && size(cell.V,1) < minV)
        cells = {}; return
    end
    if isempty(cutConic)
        cells = {cell}; return
    end
    cells = clipPolyByConic(cell, cutConic, cutX0, cutX1);
end

function cell = insertPassthroughVertices(cell, pts)
% Subdivide cell's straight boundary edges (segments, plus the two rays if unbounded) at any point
% of `pts` that lies in the OPEN interior of an existing edge and isn't already a vertex. See
% clipByFace's call site for why this is needed.
%
% BUGFIX (found via the still-open silent-wrong-answer issue from the prior session's handoff):
% `pts` is g1's/g2's ORIGINAL face vertex coordinates, while `cell.V` was built by a DIFFERENT
% arithmetic path (clipPolyHalfPlane's crossingPoint formula, chained across possibly several
% clips) -- the same cross-arithmetic-noise situation documented in assemblePieces'/
% matchHalfEdges' HISTORY, just one step earlier in the pipeline. The "already a vertex" check
% below used to share `tol`=1e-7 with onOpenSegment/onOpenRay's own matching tolerance, which was
% tighter than that noise floor (observed ~3.1e-5 on a real repro, T=(7.8665,4.6784),
% (2.6908,1.9477),(0.3892,0.7130)): a `pts` entry that geometrically IS an already-present cell
% vertex narrowly failed this check and fell through to onOpenSegment, which (correctly, by ITS
% OWN tight tolerance) still recognized p as within a hair of the edge's endpoint... except the
% version of this code before this fix inserted p as a brand-new vertex regardless, creating a
% near-zero-length sliver edge whose line equation is dominated by floating-point noise in the tiny
% direction vector, wrongly excluding a real region of the plane from its face in QuaPar.eval's
% exact (no-tolerance) membership test -- see
% maxQuaParTest.insertPassthroughVerticesDropsNearDuplicateCrossingPoint.
%
% FIX: widen ONLY this "already a vertex" pre-check (tolSnap), leaving onOpenSegment/onOpenRay's own
% matching tolerance (tol, still 1e-7) completely untouched. This was deliberately chosen over
% widening `tol` itself (tried first): `tol` also controls onOpenSegment/onOpenRay's "is p actually
% on this edge, and not too close to either of ITS OWN endpoints" tests, which are answering a
% DIFFERENT question (whether a genuinely-distinct point elsewhere on this cell's boundary should
% split it) than tolSnap (whether p merely coincides with a vertex the cell ALREADY has). Widening
% `tol` broadly perturbs the former and was observed to wrongly absorb a genuinely DISTINCT nearby
% vertex (~7.9e-4 away, a real corner, not noise) into the wrong edge on
% maxQuaParTest.dedupHitsMergesCrossingsAtACellCorner's T=(0,0),(2.11,1.43),(8.84,4.50), while
% narrower values still broke maxQuaParTest.checkOrphanHalfEdgesDropsProvablyDegenerateOrphanEdges's
% near-degenerate triangles (whose genuine small-scale features sit right around 1e-4, per that
% test's own header) -- no single value of the SHARED tolerance separated all three cases.
% Decoupling them resolves this: tolSnap only ever causes p to be treated as coincident with a
% vertex the cell construction already produced (never changes whether a genuinely new point gets
% inserted), so it cannot manufacture the kind of wrong topology `tol` did.
%
% CURVED EDGE: the cell's arc (cell.curveAfter, if any) is SKIPPED by the onOpenSegment scan. That
% test asks whether p lies on the straight segment between two consecutive vertices, which for an
% arc is its CHORD, not the boundary -- a `pts` entry on that chord is strictly inside the cell (or
% strictly outside it), never a missing corner of it, so splitting there would invent a vertex that
% is not on the boundary at all. Every insertion BEFORE the arc's own start vertex shifts it by
% one, so cell.curveAfter is updated alongside cell.V.
%
% ...but skipping the CHORD must not mean ignoring the ARC. This whole function exists to maintain
% the arrangement invariant that any two edges meet in the empty set or in a proper face (a shared
% vertex / a shared whole edge) -- wherever that fails, the missing vertex has to be added. That
% applies to an arc exactly as it does to a segment: if a `pts` point lies in the OPEN INTERIOR of
% the arc, the neighbouring face differs on either side of it, so the arc must be split there into
% two sub-arcs (of the same parabola). A piece carries only ONE curve slot (curveAfter/curveEc), so
% it cannot represent the two sub-arcs that split would produce; rather than let the arc silently
% span two different neighbours -- precisely the silent-wrong-adjacency failure mode the collinear
% straight-edge case above was fixed for -- this is DETECTED and raised.
%
% Measured, not assumed: over a randomized sweep of 109 assembled curved results, NO vertex of
% either input ever landed in the open interior of an arc (the conjugate's arc meets its
% neighbours tangentially, and the other operand's face vertices sit off it), so this guard never
% fires in the supported regime -- see maxQuaPar.m's header VALIDATION note. Lifting it means
% generalizing a piece to carry several arcs (curveAfter becoming a set), which also breaks
% clipPolyHalfPlaneCurved's "at most 2 crossings" invariant -- deliberately left out of this step.
    if isempty(pts), return; end
    tol = 1e-7;
    tolSnap = 1e-4;
    ciSkip = 0;
    if pieceIsCurved(cell), ciSkip = cell.curveAfter; end
    for pi = 1:size(pts,1)
        p = pts(pi,:);
        again = true;
        while again
            again = false;
            nv = size(cell.V,1);
            if nv == 0 || any(all(abs(cell.V - p) < tolSnap, 2)), break; end
            if ciSkip ~= 0 && onOpenArc(cell, p, tolSnap)
                error('maxQuaPar:notImplemented', ...
                    ['insertPassthroughVertices: a face vertex lies in the open interior of this ' ...
                     'cell''s parabolic edge, so that arc borders two different neighbours and must ' ...
                     'be split there into two sub-arcs. A piece carries only one curve slot, so ' ...
                     'this is not representable yet (see this function''s header).']);
            end
            if isempty(cell.dirIn)
                for i = 1:nv
                    if i == ciSkip, continue; end            % the arc: not a straight edge
                    j = mod(i,nv)+1;
                    if onOpenSegment(cell.V(i,:), cell.V(j,:), p, tol)
                        cell = insertVertexAt(cell, i, p);
                        again = true; break
                    end
                end
            else
                for i = 1:nv-1
                    if i == ciSkip, continue; end            % the arc: not a straight edge
                    if onOpenSegment(cell.V(i,:), cell.V(i+1,:), p, tol)
                        cell = insertVertexAt(cell, i, p);
                        again = true; break
                    end
                end
                if ~again && onOpenRay(cell.V(1,:), cell.dirIn, p, tol)
                    cell = insertVertexAt(cell, 0, p); again = true;
                elseif ~again && onOpenRay(cell.V(end,:), cell.dirOut, p, tol)
                    cell.V = [cell.V; p]; again = true;   % appended last: shifts no earlier index
                end
            end
        end
    end
end

function cell = insertVertexAt(cell, i, p)
% Insert p into cell.V right after row i (i==0 prepends), keeping cell.curveAfter -- an index INTO
% cell.V -- pointing at the same arc.
    cell.V = [cell.V(1:i,:); p; cell.V(i+1:end,:)];
    if cell.curveAfter > i, cell.curveAfter = cell.curveAfter + 1; end
end

function tf = onOpenArc(cell, p, tolSnap)
% Does p lie ON the cell's parabola arc, strictly between the arc's own two endpoints? Two
% conditions, both needed: p satisfies the conic (to a scale-relative tolerance), AND its frame
% parameter u lies strictly inside the arc's own u-span. Since u is a global monotonic parameter
% along the whole conic (parabolaArcFrame), the second test is exactly "between the endpoints
% along the arc" -- no chord/arc confusion. The endpoint margin reuses tolSnap, matching the
% caller's own "p is already one of this cell's vertices" pre-check, so a p that coincides with an
% arc ENDPOINT is never reported as interior to it.
    [X0, X1] = curveEndpoints(cell);
    ec = cell.curveEc;
    if abs(QuaPar.evalConic(ec, p)) > 1e-7*(1 + max(abs(ec))), tf = false; return; end
    fr = parabolaArcFrame(ec, 'maxQuaPar');
    u0 = fr.uOf(X0); u1 = fr.uOf(X1); up = fr.uOf(p);
    lo = min(u0,u1); hi = max(u0,u1);
    if hi - lo < tolSnap, tf = false; return; end
    tf = up > lo + tolSnap && up < hi - tolSnap;
end

function tf = onOpenSegment(a, b, p, tol)
    d = b - a; L = norm(d);
    if L < tol, tf = false; return; end
    t = dot(p-a, d) / L^2;
    if t <= tol/L || t >= 1 - tol/L, tf = false; return; end
    tf = norm(a + t*d - p) < tol;
end

function tf = onOpenRay(apex, dir, p, tol)
    dn = dir/norm(dir);
    t = dot(p-apex, dn);
    if t <= tol, tf = false; return; end
    tf = norm(apex + t*dn - p) < tol;
end

function cons = polyConstraints(poly)
% Outward half-plane constraints [n1 n2 c] (inside iff n*x'<=c) for every boundary edge of poly,
% in CCW order: interior is on the LEFT of each edge's direction of travel, so the outward normal
% is rot90cw(direction) = (dy,-dx) (same convention as conjPieceCPLQ.m's unitNormal).
%
% BUGFIX (found while implementing addQuaPar.m): a BOUNDED poly (dirIn empty) is a closed cycle
% of nv edges (1,2),...,(nv-1,nv),(nv,1) -- the old "for i=1:nv-1" loop always dropped the
% closing edge (nv,1), so clipByFace never enforced that one constraint of polyL (silent
% under-constraining of the clipped cell -- not caught by the existing test suite, which never
% happened to need that specific constraint). An UNBOUNDED poly's nv-1 real-vertex edges
% (1,2),...,(nv-1,nv) do NOT wrap (the two ends connect to rays, not each other), so that case is
% unaffected: `last` below is nv-1 for it, exactly as before.
    nv = size(poly.V,1);
    cons = zeros(0,3);
    last = nv - 1;
    if isempty(poly.dirIn), last = nv; end   % bounded: also include the closing edge (nv,1)
    for i = 1:last
        if poly.curveAfter == i
            % This edge is a parabolic ARC, not a segment. Emitting the half-plane of its CHORD
            % would clip against the wrong boundary -- on the bulging side it cuts away part of
            % the face, on the other it keeps territory that is not the face's. The arc is
            % applied separately, as a curved cut, by clipByFace/clipPolyByConic.
            continue
        end
        jn = mod(i,nv) + 1;
        d = poly.V(jn,:) - poly.V(i,:);
        n = [d(2), -d(1)];
        cons(end+1,:) = [n, n*poly.V(i,:)']; %#ok<AGROW>
    end
    if ~isempty(poly.dirIn)
        % BUGFIX: a ray edge's outward normal is NOT always rot90cw(-dirIn)/rot90cw(dirOut) --
        % that fixed-role formula silently assumes dirIn's own P{k} sign is always +1 and dirOut's
        % is always -1, which is exactly what obj.P{k}'s "-j=left/+j=right" convention (see
        % orderEdges) produces for a TYPICAL unbounded face (real vertices in between the two
        % rays constrain the walk enough that this pattern always holds), but is NOT guaranteed
        % for a face whose ENTIRE boundary is just the two rays sharing one apex (poly.V has a
        % single row, dirIn and dirOut anchored at the identical point) -- there, orderEdges' own
        % ray-selection logic (picking the start ray via a DIFFERENT rule, see its own comments)
        % can legitimately give BOTH rays the SAME sign. Found via a coverage/sampling stress test
        % on the paper's own V=[2 1;0 0;1 0] example (maxQuaPar:internal on g1 face 3 vs g2's
        % faces): g2 face 3's P{3}=[3 2] (both entries positive) is exactly this single-apex case,
        % and the fixed-role formula silently computed a ~12-degree sliver as "inside" when the
        % true (QuaPar.eval-confirmed ground truth) face is the ~170-degree complementary wedge.
        %
        % FIX: use the TRUE sign each ray carries in its own original P{k} entry (captured by
        % facePoly as dirInSign/dirOutSign, swapped alongside dirIn/dirOut through the CW->CCW
        % reversal) instead of assuming a fixed role-based sign. This is a strict generalization:
        % outward normal = sign * rot90ccw(direction), which reduces to EXACTLY the old formula
        % when dirInSign==+1 and dirOutSign==-1 (rot90cw(-dirIn)=+1*rot90ccw(dirIn) and
        % rot90cw(dirOut)=-1*rot90ccw(dirOut)) -- i.e. every previously-passing case (where that
        % sign pattern held) is completely unaffected; only the previously-mishandled
        % same-sign/single-apex case now differs, correctly. See
        % maxQuaParTest.matchHalfEdgesRejectsSameSideRayPairingAndDropsSubsumedPieces.
        n = poly.dirInSign * [-poly.dirIn(2), poly.dirIn(1)];
        cons(end+1,:) = [n, n*poly.V(1,:)']; %#ok<AGROW>
        n = poly.dirOutSign * [-poly.dirOut(2), poly.dirOut(1)];
        cons(end+1,:) = [n, n*poly.V(end,:)']; %#ok<AGROW>
    end
end

function s = sign2(v, tol)
    s = zeros(size(v));
    s(v > tol) = 1; s(v < -tol) = -1;
end

function poly2 = clipPolyHalfPlane(poly, nrm, c)
% Clip poly by the half-plane {nrm*x'<=c}, dispatching on whether poly carries a parabolic edge.
%
% The two paths are kept SEPARATE deliberately. clipPolyHalfPlaneStraight below is the original,
% purely polyhedral implementation, left byte-identical: its own HISTORY comment documents several
% sessions' worth of subtle, non-crashing WRONG-ANSWER bugs (a wraparound-index bug, an
% endpoint-ordering "bowtie" bug, a degenerate parallel-ray case), so generalizing it in place --
% rather than adding a parallel curve-aware path that the polyhedral case never enters -- would
% have put every one of those hard-won fixes back at risk for no benefit. Once a clip removes a
% cell's arc entirely, the result is plain polyhedral again and subsequent clips go back through
% the straight path.
    if pieceIsCurved(poly)
        poly2 = clipPolyHalfPlaneCurved(poly, nrm, c);
        return
    end
    poly2 = clipPolyHalfPlaneStraight(poly, nrm, c);
    if ~isempty(poly2), poly2.curveAfter = 0; poly2.curveEc = []; end
end

function polys = clipPolyByConic(poly, Ecut, cutX0, cutX1)
% Clip poly by the region { evalConic(Ecut,.) >= 0 } -- the CURVED-cut analogue of
% clipPolyHalfPlane, and the missing half of arc-vs-arc face clipping. Returns a LIST of polys
% (0, 1 or 2): two when the survivor would carry two arcs and has to be subdivided, which is the
% invariant QuaPar's single Ec slot per edge forces.
%
% SIGN CONVENTION. The half-plane paths keep {nrm*x' - c <= 0}; here the kept side is
% {evalConic >= 0}, which is facePoly's normalization (a face's curveEc is oriented > 0 on that
% face's own interior). To reuse the same branch shapes the value is NEGATED into `st`, so
% st <= 0 means "keep" exactly as it does there.
%
% CROSSINGS. Each boundary element restricts the cutting conic to a univariate polynomial:
%   * a straight edge or a ray  -> conicAlongRay's quadratic in t (roots in the element's range);
%   * poly's OWN parabolic edge -> clipArcByConic, which solves the quartic in the arc's frame.
% A conic can cross ONE straight edge twice, unlike a line, so "2 crossings on 2 distinct edges"
% is checked rather than assumed -- see the error below.
    nv = size(poly.V,1);
    unbounded = ~isempty(poly.dirIn);
    if nv == 0, polys = {poly}; return, end

    val = zeros(nv,1);
    for i = 1:nv, val(i) = -QuaPar.evalConic(Ecut, poly.V(i,:)); end
    sc  = max(1, max(abs(Ecut)) * max(1, max(abs(poly.V(:))))^2);
    tol = 1e-9 * sc;

    ci = poly.curveAfter;
    hasArc = (ci ~= 0);

    % ---- boundary elements, in walk order, with their crossings --------------------------
    % Node/pair numbering identical to clipPolyHalfPlaneCurved's, so the branch shapes below can
    % be read against it: for an unbounded cell node 1 is the dirIn marker, nodes 2..nv+1 are the
    % vertices and node m is the dirOut marker, so pair p connects nodes p and p+1 -- pair 1 is
    % the incoming ray, pairs 2..nv the straight edges V(p-1)->V(p), pair nv+1 the outgoing ray.
    if unbounded
        stIn  = signAtInfinity(Ecut, poly.V(1,:),    poly.dirIn);
        stOut = signAtInfinity(Ecut, poly.V(end,:),  poly.dirOut);
        st = [-stIn; sign2(val,tol); -stOut];
        edges = [(2:nv)'-1, (2:nv)'];            % straight edges, as (fromVertex, toVertex)
        edgePair = (2:nv)';                       % their pair indices
        cePair = ci + 1;                          % the arc's own pair index
    else
        st = sign2(val,tol);
        edges = [(1:nv)', mod((1:nv),nv)'+1];
        edgePair = (1:nv)';
        cePair = ci;
    end

    hits = zeros(0,3);                            % [pairIndex, x, y]
    for pe = 1:size(edges,1)
        i = edges(pe,1); j = edges(pe,2);
        p = edgePair(pe);
        if hasArc && p == cePair
            [status, Xn] = clipArcByConic(poly.curveEc, poly.V(i,:), poly.V(j,:), Ecut);
            switch status
                case 'cut'
                    Xc = Xn(1,:);
                    if norm(Xc - poly.V(i,:)) < 1e-9, Xc = Xn(2,:); end
                    hits(end+1,:) = [p, Xc]; %#ok<AGROW>
                case 'twice'
                    error('maxQuaPar:notImplemented', ...
                        ['clipPolyByConic: the cutting conic crosses this cell''s own arc TWICE, ' ...
                         'so the survivor is disconnected or bounded by two sub-arcs of the same ' ...
                         'parabola; neither is one QuaPar face.']);
            end
            continue
        end
        d = poly.V(j,:) - poly.V(i,:);
        ts = conicRootsAlong(Ecut, poly.V(i,:), d, 0, 1, tol);
        for k = 1:numel(ts)
            hits(end+1,:) = [p, poly.V(i,:) + ts(k)*d]; %#ok<AGROW>
        end
    end
    if unbounded
        tsIn = conicRootsAlong(Ecut, poly.V(1,:), poly.dirIn, 0, inf, tol);
        for k = 1:numel(tsIn)
            hits(end+1,:) = [1, poly.V(1,:) + tsIn(k)*poly.dirIn]; %#ok<AGROW>
        end
        tsOut = conicRootsAlong(Ecut, poly.V(end,:), poly.dirOut, 0, inf, tol);
        for k = 1:numel(tsOut)
            hits(end+1,:) = [nv+1, poly.V(end,:) + tsOut(k)*poly.dirOut]; %#ok<AGROW>
        end
        hits = sortrows(hits, 1);
    end

    % ---- restrict the cut to the CUTTING ARC's own span -----------------------------------
    % Ecut is a whole conic, but polyL's boundary is only the bounded arc of it between cutX0 and
    % cutX1. A crossing of the extended conic that lies outside that span is not a crossing of
    % polyL's boundary at all -- polyL is bounded there by its straight edges, which have already
    % been applied -- so keeping it would cut the cell along a curve that is not a boundary. It
    % shows up as a single unpaired crossing, which no keep/cut branch can act on.
    if ~isempty(cutX0)
        frC = parabolaArcFrame(Ecut, 'maxQuaPar');
        ua = frC.uOf(cutX0); ub = frC.uOf(cutX1);
        ulo = min(ua,ub); uhi = max(ua,ub); uspan = uhi - ulo;
        keep = false(size(hits,1),1);
        for k = 1:size(hits,1)
            uk = frC.uOf(hits(k,2:3));
            keep(k) = (uk >= ulo - 1e-7*max(1,uspan)) && (uk <= uhi + 1e-7*max(1,uspan));
        end
        hits = hits(keep,:);
    end

    % ---- the cutting arc may END inside the cell ------------------------------------------
    % polyL's arc is BOUNDED. If it enters polyK through one edge and terminates at one of its own
    % endpoints, that endpoint is where polyL's boundary hands over from the arc to a straight
    % edge -- and that straight edge's half-plane has already been applied, so the endpoint lies
    % ON the clipped cell's boundary. It is then the second crossing, and without it the clip
    % reports a single unpaired hit that no keep/cut branch can act on.
    if ~isempty(cutX0) && size(hits,1) == 1
        for X = {cutX0, cutX1}
            pIdx = pairOnBoundary(poly, X{1}, edges, edgePair, unbounded, nv);
            if pIdx > 0 && pIdx ~= hits(1,1)
                hits(end+1,:) = [pIdx, X{1}]; %#ok<AGROW>
                break
            end
        end
        hits = sortrows(hits, 1);
    end

    % ---- no crossing: the whole cell is on one side --------------------------------------
    if isempty(hits)
        % No crossing WITHIN the arc's span. The cell therefore lies wholly on one side of
        % polyL's curved boundary; decide it at an interior sample rather than from the vertex
        % signs, which are signs of the EXTENDED conic and can disagree away from the arc.
        smp = interiorSample(poly);
        if isempty(smp), polys = {poly}; return, end
        if QuaPar.evalConic(Ecut, smp) >= -tol, polys = {poly}; else, polys = {}; end
        return
    end

    if unbounded
        if size(hits,1) ~= 2 || hits(1,1) == hits(2,1)
            error('maxQuaPar:notImplemented', ...
                ['clipPolyByConic: expected 2 crossings on 2 distinct boundary elements of an ' ...
                 'unbounded cell, found %d on %s.'], size(hits,1), mat2str(hits(:,1)'));
        end
        p1 = hits(1,1); p2 = hits(2,1);
        Xa = hits(1,2:3); Xb = hits(2,2:3);
        if st(1) <= tol
            % Both ray ends are on the KEPT side: the cut removes a middle bulge, and both rays
            % survive. Same surgery as the straight path's corresponding branch, except that the
            % new Xa->Xb edge is an ARC of the cutting conic rather than a segment.
            Vnew = [poly.V(1:p1-1,:); Xa; Xb; poly.V(p2:nv,:)];
            cutEdge = p1;                          % the new edge follows the kept head
            out = mkPoly(Vnew, cutEdge, orientConicInto(Ecut, [Xa; Xb; poly.V(:,:)]), ...
                         poly.dirIn, poly.dirInSign, poly.dirOut, poly.dirOutSign);
            qCurve = 0;
            if hasArc
                if cePair <= p1,      qCurve = cePair - 1;
                elseif cePair >= p2,  qCurve = cePair - p2 + p1 + 1;
                end
                if qCurve < 1 || qCurve > size(Vnew,1)-1, qCurve = 0; end
            end
        else
            % Both ray ends are on the REMOVED side: only the middle survives, and the result is
            % BOUNDED -- closed by the arc from Xb back to Xa.
            Vnew = [Xa; poly.V(p1:p2-1,:); Xb];
            cutEdge = size(Vnew,1);
            out = mkPoly(Vnew, cutEdge, orientConicInto(Ecut, Vnew), [], [], [], []);
            qCurve = 0;
            if hasArc && cePair >= p1 && cePair <= p2
                qCurve = cePair - p1 + 1;
                if qCurve > size(Vnew,1)-1, qCurve = 0; end
            end
        end
        if size(out.V,1) < 2, polys = {}; return, end
        if qCurve == 0, polys = {out}; return, end
        polys = num2cell(splitTwoArcPiece(out, qCurve, poly.curveEc));
        return
    end
    if size(hits,1) ~= 2 || hits(1,1) == hits(2,1)
        error('maxQuaPar:notImplemented', ...
            ['clipPolyByConic: expected 2 crossings on 2 distinct edges, found %d on edges %s. ' ...
             'A conic can cut one straight edge twice, unlike a line, so this configuration is ' ...
             'real and needs its own branch.'], size(hits,1), mat2str(hits(:,1)'));
    end

    % ---- bounded, two crossings: keep the chain that is inside ----------------------------
    p1 = hits(1,1); p2 = hits(2,1);
    Xa = hits(1,2:3); Xb = hits(2,2:3);
    midIdx = mod((p1):(p2-1), nv) + 1;
    if isempty(midIdx) || all(st(midIdx) <= tol)
        Vnew = [Xa; poly.V(midIdx,:); Xb];
        pStart = p1;
    else
        keepIdx = mod((p2):(p1-1+nv), nv) + 1;
        Vnew = [Xb; poly.V(keepIdx,:); Xa];
        pStart = p2;
    end
    qCurve = 0;
    if hasArc
        qCurve = mod(ci - pStart, nv) + 1;
        if qCurve > numel(Vnew)-1, qCurve = 0; end   % the inherited arc did not survive
    end
    cutEdge = numel(Vnew);                            % the closing edge IS the new cut arc

    out.V = Vnew;
    out.dirIn = []; out.dirOut = []; out.dirInSign = []; out.dirOutSign = [];
    out.curveAfter = cutEdge;
    out.curveEc = orientConicInto(Ecut, Vnew);
    out.f = [];
    if size(out.V,1) < 2, polys = {}; return, end

    if qCurve == 0
        polys = {out};                                % one arc only: the new cut
        return
    end
    % Two arcs -- the inherited one and the cut -- which one QuaPar face cannot hold. This is
    % exactly what splitTwoArcPiece exists for, and is the same subdivide-never-widen rule
    % splitCell already follows.
    pieces = splitTwoArcPiece(out, qCurve, poly.curveEc);
    polys = num2cell(pieces);
end

function p = pairOnBoundary(poly, X, edges, edgePair, unbounded, nv)
% The pair index of the boundary element X lies on, or 0. Straight elements are tested by
% collinearity plus range; rays by direction plus a nonnegative parameter. Used only to recognise
% a cutting arc's own endpoint as a crossing, so an exact-ish tolerance is what is wanted.
    p = 0;
    tolG = 1e-7 * max(1, max(abs(poly.V(:))));
    for pe = 1:size(edges,1)
        A = poly.V(edges(pe,1),:); B = poly.V(edges(pe,2),:);
        d = B - A; L = norm(d);
        if L < 1e-12, continue, end
        t = dot(X - A, d) / L^2;
        if t < -1e-9 || t > 1 + 1e-9, continue, end
        if norm(X - (A + t*d)) <= tolG, p = edgePair(pe); return, end
    end
    if unbounded
        for which = [1 2]
            if which == 1, base = poly.V(1,:);   dir = poly.dirIn;  pi0 = 1;
            else,          base = poly.V(nv,:);  dir = poly.dirOut; pi0 = nv+1;
            end
            t = dot(X - base, dir);
            if t < -1e-9, continue, end
            if norm(X - (base + t*dir)) <= tolG, p = pi0; return, end
        end
    end
end

function p = mkPoly(V, curveAfter, curveEc, dirIn, dirInSign, dirOut, dirOutSign)
% Assemble one clipPolyByConic result, running the same consecutive-duplicate removal with index
% tracking that finishCurved uses so curveAfter keeps pointing at the same edge.
    [V, curveAfter] = dedupConsecutiveTracked(V, curveAfter);
    p.V = V;
    p.dirIn = dirIn; p.dirOut = dirOut;
    p.dirInSign = dirInSign; p.dirOutSign = dirOutSign;
    p.curveAfter = curveAfter;
    if curveAfter == 0, p.curveEc = []; else, p.curveEc = curveEc; end
    p.f = [];        % filled by the caller once the winner is decided; splitTwoArcPiece and
                     % subPiece both copy this field through, so it must exist by then
end

function s = signAtInfinity(Ecut, apex, dir)
% Sign of the cutting conic far along apex + t*dir, from the leading nonzero coefficient of its
% restriction -- the ray analogue of evaluating at a vertex.
    [A,B,C] = conicAlongRayLocal(Ecut, apex, dir);
    if abs(A) > 1e-12, s = sign(A); elseif abs(B) > 1e-12, s = sign(B); else, s = sign(C); end
end

function ts = conicRootsAlong(Ecut, base, dir, tlo, thi, tol)
% Roots of t -> evalConic(Ecut, base + t*dir) strictly inside (tlo, thi), sign-change only.
    [A,B,C] = conicAlongRayLocal(Ecut, base, dir);
    ts = [];
    if abs(A) <= 1e-14
        if abs(B) > 1e-14, r = -C/B; else, return, end
    else
        disc = B^2 - 4*A*C;
        if disc < 0, return, end
        sq = sqrt(disc);
        r = [(-B - sq)/(2*A), (-B + sq)/(2*A)];
    end
    for k = 1:numel(r)
        t = r(k);
        if t <= tlo + 1e-9 || t >= thi - 1e-9, continue, end
        h = 1e-7 * max(1, abs(t));
        v0 = polyval([A B C], t - h); v1 = polyval([A B C], t + h);
        if sign2(v0,tol)*sign2(v1,tol) < 0, ts(end+1) = t; end %#ok<AGROW>
    end
    ts = sort(ts);
end

function [A,B,C] = conicAlongRayLocal(ecRow, apex, dir)
% Coefficients of t -> evalConic(ecRow, apex + t*dir) = A t^2 + B t + C. Same three lines as
% addQuaPar's own conicAlongRay, duplicated because that one is file-local there.
    a=ecRow(1); b=ecRow(2); c=ecRow(3); d=ecRow(4); e=ecRow(5);
    A = a*dir(1)^2 + b*dir(1)*dir(2) + c*dir(2)^2;
    B = 2*a*apex(1)*dir(1) + b*(apex(1)*dir(2)+apex(2)*dir(1)) + 2*c*apex(2)*dir(2) ...
        + d*dir(1) + e*dir(2);
    C = QuaPar.evalConic(ecRow, apex);
end

function Ec = orientConicInto(Ec, V)
% facePoly's convention: a face's stored curveEc must be > 0 on that face's OWN interior. The cut
% keeps {evalConic(Ecut,.) >= 0}, so Ecut already has that sign -- but the check is cheap and the
% convention is the one clipPolyHalfPlaneCurved and assignSide both rely on, so it is enforced
% here rather than assumed.
    if size(V,1) < 3, return, end
    ctr = mean(V,1);
    if QuaPar.evalConic(Ec, ctr) < 0, Ec = -Ec; end
end

function poly2 = clipPolyHalfPlaneCurved(poly, nrm, c)
% Clip a poly carrying exactly ONE curved (parabolic) boundary edge by the half-plane {nrm*x'<=c}.
%
% Same node/pair indexing and the same four keep/cut branches as clipPolyHalfPlaneStraight (see its
% header for the derivation of each); what is added here is (a) crossings that land on the arc are
% computed ON THE ARC (clipArcByHalfPlane) instead of on its chord, and (b) the arc's position is
% threaded through each branch's vertex-list surgery.
%
% Threading the arc: poly.curveAfter is a VERTEX index -- the arc runs poly.V(i) -> poly.V(i+1)
% (wrapping to poly.V(nv)->poly.V(1) when bounded) -- the same convention splitCell/boundedPiece
% already use. Internally it is first converted to this function's own boundary-EDGE/node-pair
% numbering (identical to the straight path's `pairs`): for a bounded poly, edge p is
% V(p)->V(p mod nv+1), so cE = i; for an unbounded one, node 1 is the dirIn ray marker and node j+1
% is V(j), so edge p connects V(p-1)->V(p) and cE = i+1. Each branch below then maps original edge
% indices to output edge indices in closed form (see the per-branch comments), and the result is
% converted back to a vertex index. Rays are never curved (assertCurvedEdgesAreArcs), so cE is
% never 1 or m.
%
% NOT IMPLEMENTED, and detected rather than silently mishandled: the clip line cutting the arc
% TWICE while both of the arc's own endpoints stay on the same side (the arc bulging across it in
% between). Unlike a convex polygon, a face bounded on one side by a parabola is not convex, so
% "at most 2 boundary crossings, on 2 distinct edges" -- the invariant the straight path is built
% on -- can genuinely fail. The result of such a clip is either DISCONNECTED (two components) or a
% single region whose boundary carries TWO separate arcs; neither is representable as one QuaPar
% face, so this errors loudly (arcBulgesAcross below).
    nv = size(poly.V,1);
    tol = 1e-9*(1+abs(c)+norm(nrm));
    val = poly.V*nrm' - c;
    unbounded = ~isempty(poly.dirIn);
    ci = poly.curveAfter;
    iA = ci; iB = mod(ci, nv) + 1;               % the arc's own two endpoints, in walk order
    X0 = poly.V(iA,:); X1 = poly.V(iB,:);
    if unbounded, cE = ci + 1; else, cE = ci; end

    if arcBulgesAcross(poly.curveEc, X0, X1, nrm, c)
        error('maxQuaPar:notImplemented', ...
            ['clipPolyHalfPlaneCurved: the clip line cuts this cell''s parabolic edge TWICE (both ' ...
             'arc endpoints on the same side, the arc bulging across in between). The clipped cell ' ...
             'is then either disconnected or bounded by two separate arcs, neither of which is ' ...
             'representable as a single QuaPar face.']);
    end

    if unbounded
        inSign = sign2(poly.dirIn*nrm', tol);
        if inSign == 0, inSign = sign2(val(1), tol); end
        outSign = sign2(poly.dirOut*nrm', tol);
        if outSign == 0, outSign = sign2(val(end), tol); end
        st = [inSign; sign2(val,tol); outSign];
    else
        st = sign2(val,tol);
    end
    m = numel(st);

    % Sound now that arcBulgesAcross has ruled out the arc leaving and re-entering: with every node
    % inside (resp. outside) and no bulge, the WHOLE boundary -- arc included -- is inside (outside).
    if all(st <= 0), poly2 = poly; return; end
    if all(st >= 0) && any(st > 0), poly2 = []; return; end

    pairs = [(1:m-1)', (2:m)'];
    if ~unbounded, pairs = [pairs; m 1]; end
    cross = find((st(pairs(:,1)) > 0) ~= (st(pairs(:,2)) > 0));
    xpt = @(p) crossingPointCurved(poly, val, nrm, c, p, m, unbounded, cE, X0, X1, tol);

    if ~unbounded
        if numel(cross) ~= 2
            error('maxQuaPar:internal', ...
                'clipPolyHalfPlaneCurved: expected 2 crossings on a bounded cell, found %d.', numel(cross));
        end
        p1 = cross(1); p2 = cross(2); Xa = xpt(p1); Xb = xpt(p2);
        midIdx = mod((p1):(p2-1), nv) + 1;
        if isempty(midIdx) || all(st(midIdx) <= 0)
            Vnew = [Xa; poly.V(midIdx,:); Xb];   % kept arc runs forward from Xa; cut edge closes it
            pStart = p1;
        else
            keepIdx = mod((p2):(p1-1+nv), nv) + 1;
            Vnew = [Xb; poly.V(keepIdx,:); Xa];  % kept arc is the wrapped complement (see straight path)
            pStart = p2;
        end
        % Output edge q (q=1..numel(Vnew)-1) is the surviving part of original edge
        % mod(pStart+q-2,nv)+1; output edge numel(Vnew) is the brand-new straight cut edge.
        qCurve = mod(cE - pStart, nv) + 1;
        if qCurve > numel(Vnew)-1, qCurve = 0; end
        poly2 = finishCurved(Vnew, qCurve, poly, [], [], [], [], true);
        return
    end

    if numel(cross) == 2
        p1 = cross(1); p2 = cross(2); Xa = xpt(p1); Xb = xpt(p2);
        if st(1) <= 0   % both ray ends inside: discard the middle bulge, keep both rays
            keepBefore = poly.V(1:p1-1,:);
            keepAfter  = poly.V(p2:nv,:);
            Vnew = [keepBefore; Xa; Xb; keepAfter];
            % Output edges: 1..p1-1 are original edges 1..p1-1; p1 is the surviving tail of edge
            % p1; p1+1 is the new cut edge; p1+2 is the surviving head of edge p2; and p1+2+r is
            % original edge p2+r.
            if cE <= p1
                qE = cE;
            elseif cE >= p2
                qE = cE - p2 + p1 + 2;
            else
                qE = 0;
            end
            poly2 = finishCurved(Vnew, max(qE-1,0), poly, poly.dirIn, poly.dirInSign, ...
                poly.dirOut, poly.dirOutSign, false);
        else            % both ray ends outside: keep only the middle; the result becomes bounded
            Vnew = [Xa; poly.V((p1):(p2-1),:); Xb];
            % Output edge q is original edge p1+q-1 for q=1..p2-p1+1; the last is the cut edge.
            if cE >= p1 && cE <= p2, qE = cE - p1 + 1; else, qE = 0; end
            poly2 = finishCurved(Vnew, qE, poly, [], [], [], [], true);
        end
        return
    end
    if numel(cross) ~= 1
        error('maxQuaPar:internal', ...
            'clipPolyHalfPlaneCurved: expected 1 or 2 crossings on an unbounded cell, found %d.', numel(cross));
    end
    % Exactly one end outside: keep the inside end up to the single crossing and replace the
    % outside end's ray with one running along the clip line (same construction, and the same
    % default replaced-ray signs, as the straight path -- see its own comment for the derivation).
    p = cross(1); X = xpt(p);
    if st(1) <= 0   % dirIn side inside, dirOut side outside: keep 1..p, replace dirOut
        Vnew = [poly.V(1:min(p-1,nv),:); X];
        % Output edge k is original edge k for k=1..p (edge p surviving as its tail V(p-1)->X).
        if cE <= p, qE = cE; else, qE = 0; end
        poly2 = finishCurved(Vnew, max(qE-1,0), poly, poly.dirIn, poly.dirInSign, ...
            [-nrm(2), nrm(1)], -1, false);
    else            % dirOut side inside, dirIn side outside: keep p+1..m, replace dirIn
        Vnew = [X; poly.V(max(p,1):nv,:)];
        % Output edge 1 is the new ray; output edge 2 is the surviving head X->V(p) of original
        % edge p; output edge 2+r is original edge p+r.
        if cE >= p, qE = cE - p + 2; else, qE = 0; end
        poly2 = finishCurved(Vnew, max(qE-1,0), poly, [nrm(2), -nrm(1)], 1, ...
            poly.dirOut, poly.dirOutSign, false);
    end
end

function poly2 = finishCurved(Vnew, curveAfter, poly, dirIn, dirInSign, dirOut, dirOutSign, bounded)
% Assemble one clipPolyHalfPlaneCurved result: run the same consecutive-duplicate removal the
% straight path uses, keeping curveAfter (an index INTO Vnew) pointing at the same arc, then apply
% the minimum-vertex test. A BOUNDED result normally needs 3 vertices, but only 2 when it keeps the
% arc (a "lens": one parabola arc closed by one straight cut edge -- see clipByFace's own note).
    [Vnew, curveAfter] = dedupConsecutiveTracked(Vnew, curveAfter);
    if isempty(Vnew), poly2 = []; return; end
    minV = 3;
    if ~bounded, minV = 1; elseif curveAfter ~= 0, minV = 2; end
    if size(Vnew,1) < minV, poly2 = []; return; end
    poly2.V = Vnew;
    poly2.dirIn = dirIn; poly2.dirOut = dirOut;
    poly2.dirInSign = dirInSign; poly2.dirOutSign = dirOutSign;
    poly2.curveAfter = curveAfter;
    if curveAfter == 0, poly2.curveEc = []; else, poly2.curveEc = poly.curveEc; end
end

function [V, q] = dedupConsecutiveTracked(V, q)
% dedupConsecutive with index tracking: dropping row d (because it duplicates row d-1) shifts every
% index above d down by one, and maps d itself onto d-1 -- both covered by subtracting the number
% of dropped rows at or below q.
%
% DEGENERATE ARC: the arc's own two endpoint rows can themselves collapse together, so this cannot
% just shift q blindly. That happens whenever the clip line is TANGENT to the parabola at one of
% the arc's own endpoints -- not a contrived case but a structural one here, since a conjugate's
% curved edge meets its neighbouring straight edges exactly tangentially (the conjugate is C1 where
% its pieces join), and the OTHER operand's face boundaries frequently run along that same tangent
% line. The clip then keeps only a single point of the arc, so the arc is GONE from the result and q
% must become 0. Shifting it instead silently moved the "this edge is a parabola" label onto the
% neighbouring STRAIGHT edge, which then failed to pair with that edge's genuine straight neighbour
% (maxQuaPar:internal, "no matching neighbour") -- found on the very first curved end-to-end
% fixture, maxQuaParTest.maxQuaParCombinesOneCurvedInputWithAPolyhedralOne.
    tol = sqrt(eps);
    n = size(V,1);
    keep = true(n,1);
    for i = 2:n
        if norm(V(i,:)-V(i-1,:)) < tol, keep(i) = false; end
    end
    if n > 1 && norm(V(1,:)-V(end,:)) < tol, keep(n) = false; end
    if q > 0
        qb = mod(q, n) + 1;                       % the arc runs row q -> row qb
        if norm(V(q,:) - V(qb,:)) < tol
            q = 0;                                % zero-length arc: it did not survive the clip
        else
            q = q - sum(~keep(1:min(q,n)));
            if q < 1, q = 0; end
        end
    end
    V = V(keep,:);
end

function tf = arcBulgesAcross(Ec, X0, X1, nrm, c)
% True iff the arc from X0 to X1 leaves the half-plane {nrm*x'<=c} and comes back (or vice versa)
% strictly BETWEEN its own two endpoints -- i.e. both endpoints agree on which side they are on,
% but the arc's own extremum of nrm*x'-c along the arc lies in between and disagrees. Restricted
% to the conic, nrm*x'-c is the explicit quadratic A2*u^2+A1*u+A0 in the frame's parameter u
% (parabolaArcFrame), so the extremum is at u=-A1/(2*A2) in closed form; a genuine parabola arc can
% cross a line at most twice, so checking that single stationary point is exhaustive.
    tolS = 1e-8;
    v0 = nrm*X0' - c; v1 = nrm*X1' - c;
    if (v0 <= tolS) ~= (v1 <= tolS), tf = false; return; end   % opposite sides: one crossing, fine
    fr = parabolaArcFrame(Ec, 'maxQuaPar');
    A = fr.lineCoeffs(nrm, c);
    if abs(A(1)) <= 1e-12*(1+abs(A(2))+abs(A(3))), tf = false; return; end   % affine in u: no bulge
    ustar = -A(2)/(2*A(1));
    u0 = fr.uOf(X0); u1 = fr.uOf(X1);
    lo = min(u0,u1); hi = max(u0,u1);
    if ustar <= lo || ustar >= hi, tf = false; return; end
    vstar = A(1)*ustar^2 + A(2)*ustar + A(3);
    inEnds = (v0 <= tolS);
    marg = 1e-9*(1 + abs(v0) + abs(v1));
    tf = (inEnds && vstar > marg) || (~inEnds && vstar < -marg);
end

function pt = crossingPointCurved(poly, val, nrm, c, pairIdx, m, unbounded, cE, X0, X1, tol)
% crossingPoint, except that a crossing landing on the cell's parabolic edge is solved ON THE ARC
% (clipArcByHalfPlane) rather than on the arc's chord.
    if pairIdx ~= cE
        pt = crossingPoint(poly, val, nrm, c, pairIdx, m, unbounded);
        return
    end
    x0Inside = (nrm*X0' - c) <= tol;
    [status, Xnew] = clipArcByHalfPlane(poly.curveEc, X0, X1, nrm, c);
    if strcmp(status, 'cut')
        if x0Inside, pt = Xnew(2,:); else, pt = Xnew(1,:); end
        return
    end
    % clipArcByHalfPlane classified the arc as wholly inside/outside using its own endpoint
    % tolerance, while this function's caller (sign2, on a scaled tolerance) saw the endpoints
    % straddle the line. The two only disagree when the sign flip sits within tolerance of one of
    % the arc's own endpoints -- in which case that endpoint IS the crossing, to available precision.
    if x0Inside, pt = X1; else, pt = X0; end
end

function poly2 = clipPolyHalfPlaneStraight(poly, nrm, c)
% Clip poly by the half-plane {nrm*x'<=c}. A line crosses a convex region's boundary in at most 2
% points, so only 0, 1 (unbounded only, mixed end-status), or 2 crossings are possible.
%
% Node indexing: bounded (closed cycle, m=nv) -- node i = poly.V(i,:), pairs (1,2)..(nv,1) (wraps).
% Unbounded (open chain, m=nv+2) -- node 1 = "the dirIn ray" (no fixed point: status = sign of the
% ASYMPTOTIC value along that ray, i.e. sign(dirIn*nrm')); node i for i=2..nv+1 = poly.V(i-1,:);
% node m = "the dirOut ray" (status = sign(dirOut*nrm')). Pairs (1,2)..(m-1,m), NOT wrapped -- an
% open chain has two distinct ends, not a cycle. A crossing on pair (1,2) or (m-1,m) lies ON the
% corresponding ray itself (solved by the ray's own parametrization); every other pair is an
% ordinary segment-vs-line intersection.
%
% DEGENERATE CASE (ray parallel to the clip boundary, dir*nrm'~0): the "asymptotic trend" sign is
% meaningless then -- the constraint value along that ray is CONSTANT (equal to its value at the
% ray's own finite anchor vertex), not trending to +-inf, so sign2 correctly returning 0 for the
% trend must NOT be fed into the crossing-detection logic as a genuine tie (that wrongly flags a
% "crossing" between the ray and its neighbouring vertex whenever their true, constant signs
% differ, and crossingPoint's ray formula then divides by this near-zero denom and returns a
% huge/NaN "intersection" that is not actually on the boundary -- see maxQuaPar.m header
% HISTORY). Fall back to the anchor vertex's own (non-asymptotic) sign in that case.
    nv = size(poly.V,1);
    tol = 1e-9*(1+abs(c)+norm(nrm));
    val = poly.V*nrm' - c;
    unbounded = ~isempty(poly.dirIn);

    if unbounded
        inSign = sign2(poly.dirIn*nrm', tol);
        if inSign == 0, inSign = sign2(val(1), tol); end
        outSign = sign2(poly.dirOut*nrm', tol);
        if outSign == 0, outSign = sign2(val(end), tol); end
        st = [inSign; sign2(val,tol); outSign];
    else
        st = sign2(val,tol);
    end
    m = numel(st);

    if all(st <= 0)
        poly2 = poly; return               % fully inside (incl. all-tied-to-boundary case)
    end
    if all(st >= 0) && any(st > 0)
        poly2 = []; return                  % fully outside
    end

    pairs = [(1:m-1)', (2:m)'];
    if ~unbounded, pairs = [pairs; m 1]; end
    cross = find((st(pairs(:,1)) > 0) ~= (st(pairs(:,2)) > 0));

    xpt = @(p) crossingPoint(poly, val, nrm, c, p, m, unbounded);

    if ~unbounded
        if numel(cross) ~= 2
            error('maxQuaPar:internal', ...
                'clipPolyHalfPlane: expected 2 crossings on a bounded poly, found %d.', numel(cross));
        end
        p1 = cross(1); p2 = cross(2); X1 = xpt(p1); X2 = xpt(p2);
        % Closed polygon: keep whichever arc (p1+1..p2, or its wrapped complement) is inside; close
        % with a chord between the two crossing points either way. p1<p2 always (cross is found in
        % ascending pair-index order), so the wrapped complement p2..p1-1 must add nv before
        % modding -- MATLAB's colon operator does not wrap on its own, and p2:(p1-1) (p2>p1-1) is
        % simply empty rather than the intended wrapped range (see maxQuaPar.m header HISTORY: this
        % silently turned "keep the far arc" into "keep nothing," producing a spurious empty cell).
        midIdx = mod((p1):(p2-1), nv) + 1;
        if isempty(midIdx) || all(st(midIdx) <= 0)
            % Kept arc runs forward from X1 to X2 (through the mid vertices): X1, mid..., X2, then
            % the new cut edge closes X2 back to X1.
            Vnew = dedupConsecutive([X1; poly.V(midIdx,:); X2]);
        else
            % Kept arc is the WRAPPED complement, running forward from X2 (through p2+1..p1) back
            % to X1: X2, rest..., X1, then the new cut edge closes X1 back to X2. Swapped relative
            % to the mid case above -- X1 and X2 are always adjacent via the cut edge regardless of
            % which arc is kept, so whichever crossing begins the kept arc must come first (see
            % maxQuaPar.m header HISTORY: putting X1 first here produced a self-intersecting
            % "bowtie" cell, since X1 would then be wrongly wired to the FAR end of the kept arc).
            keepIdx = mod((p2):(p1-1+nv), nv) + 1;
            Vnew = dedupConsecutive([X2; poly.V(keepIdx,:); X1]);
        end
        if size(Vnew,1) < 3, poly2 = []; return; end
        poly2.V = Vnew; poly2.dirIn = []; poly2.dirOut = [];
        poly2.dirInSign = []; poly2.dirOutSign = [];
        return
    end

    if numel(cross) == 2
        p1 = cross(1); p2 = cross(2); X1 = xpt(p1); X2 = xpt(p2);
        % st(1)==st(m) is forced whenever there are exactly 2 crossings (each crossing flips the
        % running status once; 2 flips return to the starting status) -- see header derivation.
        vIdxMid = (p1):(p2-1);   % real V-indices strictly between the two crossings
        if st(1) <= 0   % both ray ends inside: discard the middle bulge, keep both rays (possibly
                         % with a new, closer apex if a crossing landed ON a ray -- automatic, since
                         % an empty keepBefore/keepAfter just leaves X1/X2 as the new end vertex).
            keepBefore = poly.V(1:p1-1,:);
            keepAfter  = poly.V(p2:nv,:);
            poly2.V = dedupConsecutive([keepBefore; X1; X2; keepAfter]);
            poly2.dirIn = poly.dirIn; poly2.dirOut = poly.dirOut;
            poly2.dirInSign = poly.dirInSign; poly2.dirOutSign = poly.dirOutSign;
        else            % both ray ends outside: keep only the middle; result becomes bounded
            Vnew = dedupConsecutive([X1; poly.V(vIdxMid,:); X2]);
            if size(Vnew,1) < 3, poly2 = []; return; end
            poly2.V = Vnew; poly2.dirIn = []; poly2.dirOut = [];
            poly2.dirInSign = []; poly2.dirOutSign = [];
        end
        return
    end
    if numel(cross) ~= 1
        error('maxQuaPar:internal', ...
            'clipPolyHalfPlane: expected 1 or 2 crossings on an unbounded poly, found %d.', numel(cross));
    end
    % Exactly one end is outside (mixed end-status): keep the inside end's real vertices up to the
    % single crossing, and REPLACE the outside end's ray with a new one running along the clip line
    % itself (direction chosen so the retained CCW region has this new edge's interior on its left,
    % i.e. its outward normal is exactly nrm -- see header derivation: rot90ccw(nrm)=(-nrm2,nrm1)
    % for a trailing/outgoing replacement, rot90cw(nrm)=(nrm2,-nrm1) for a leading/incoming one).
    % The REPLACED ray's sign is fixed at its "default" value (dirInSign=+1, dirOutSign=-1) -- the
    % one value for which polyConstraints' sign-aware normal (see its own HISTORY) reconstructs
    % exactly nrm itself for this brand-new edge, by construction; the KEPT ray's sign passes
    % through unchanged from `poly` (see polyConstraints' HISTORY -- this is what lets a piece
    % that never has either of its original rays clipped away keep its true, possibly
    % non-default, sign all the way through to dropSubsumedPieces' isSubsumed).
    p = cross(1); X = xpt(p);
    if st(1) <= 0   % dirIn side inside, dirOut side outside: keep 1..p, replace dirOut
        keepV = poly.V(1:min(p-1,nv),:);
        poly2.V = dedupConsecutive([keepV; X]);
        poly2.dirIn = poly.dirIn; poly2.dirInSign = poly.dirInSign;
        poly2.dirOut = [-nrm(2), nrm(1)]; poly2.dirOutSign = -1;
    else            % dirOut side inside, dirIn side outside: keep p+1..m, replace dirIn
        keepV = poly.V(max(p,1):nv,:);
        poly2.V = dedupConsecutive([X; keepV]);
        poly2.dirIn = [nrm(2), -nrm(1)]; poly2.dirInSign = 1;
        poly2.dirOut = poly.dirOut; poly2.dirOutSign = poly.dirOutSign;
    end
    if isempty(poly2.V), poly2 = []; end
end

function pt = crossingPoint(poly, val, nrm, c, pairIdx, m, unbounded)
% Crossing point of the half-plane boundary {nrm*x'=c} with the boundary pair (pairIdx,pairIdx+1).
    if unbounded && pairIdx == 1
        v0 = poly.V(1,:)*nrm' - c; denom = poly.dirIn*nrm';
        pt = poly.V(1,:) - (v0/denom)*poly.dirIn;
        return
    end
    if unbounded && pairIdx == m-1
        vn = poly.V(end,:)*nrm' - c; denom = poly.dirOut*nrm';
        pt = poly.V(end,:) - (vn/denom)*poly.dirOut;
        return
    end
    if unbounded
        i = pairIdx - 1; j = pairIdx;    % node(pairIdx)=V(pairIdx-1), node(pairIdx+1)=V(pairIdx)
        vA = val(i); vB = val(j);
        pt = poly.V(i,:) + (vA/(vA-vB))*(poly.V(j,:)-poly.V(i,:));
    else
        i = pairIdx; j = mod(pairIdx, m) + 1;
        vA = val(i); vB = val(j);
        pt = poly.V(i,:) + (vA/(vA-vB))*(poly.V(j,:)-poly.V(i,:));
    end
end

function V = dedupConsecutive(V)
    tol = sqrt(eps);
    keep = true(size(V,1),1);
    for i = 2:size(V,1)
        if norm(V(i,:)-V(i-1,:)) < tol, keep(i) = false; end
    end
    if size(V,1) > 1 && norm(V(1,:)-V(end,:)) < tol, keep(end) = false; end
    V = V(keep,:);
end

% ============================================================================================
% ----- deciding whether one candidate face-function dominates a whole cell -------------------
function [decided, winRow] = decideWinner(cell, f1row, f2row)
% True iff f1row or f2row is >= the other over the WHOLE cell: check every finite vertex, plus (for
% an unbounded cell) the closed-form asymptotic behaviour along each ray (leading quadratic term,
% then linear, then value at the apex -- whichever is first nonzero decides the sign as t->inf).
    diffRow = f1row - f2row;
    vals = [];
    if size(cell.V,1) > 0
        vals = QuaPar.evalPoly(diffRow, cell.V);
    end
    if pieceIsCurved(cell)
        % A curved cell's boundary bulges away from the chord between the arc's two endpoints, so
        % the vertices alone can miss a sign change that happens along the arc. Add the arc's own
        % midpoint (a genuine point of the closed cell, and not one of its vertices) as an extra
        % sample; this can only make the domination test STRICTER, never wrongly "decide" a cell.
        vals(end+1,1) = QuaPar.evalPoly(diffRow, arcMidpoint(cell)); %#ok<AGROW>
    end
    if ~isempty(cell.dirIn)
        vals(end+1,1) = asymptoticSign(diffRow, cell.V(1,:), cell.dirIn); %#ok<AGROW>
        vals(end+1,1) = asymptoticSign(diffRow, cell.V(end,:), cell.dirOut); %#ok<AGROW>
    end
    tol = 1e-9*(1+max(abs(vals)));
    if all(vals >= -tol)
        decided = true; winRow = f1row; return
    end
    if all(vals <= tol)
        decided = true; winRow = f2row; return
    end
    decided = false; winRow = [];
end

function [X0, X1] = curveEndpoints(piece)
% The parabolic edge's own two endpoints, in the piece's CCW walk order: piece.V(i) -> piece.V(i+1)
% for i = piece.curveAfter (wrapping to piece.V(nv)->piece.V(1) for a bounded piece; an unbounded
% piece's arc always sits strictly between two real vertices, so the wrap never triggers there).
    nv = size(piece.V,1);
    i = piece.curveAfter;
    X0 = piece.V(i,:);
    X1 = piece.V(mod(i,nv)+1,:);
end

function pt = arcMidpoint(piece)
% The point halfway (in the parabola frame's own parameter u) along the piece's arc -- a genuine
% boundary point of the piece that is not one of its vertices.
    [X0, X1] = curveEndpoints(piece);
    fr = parabolaArcFrame(piece.curveEc, 'maxQuaPar');
    pt = fr.point(0.5*(fr.uOf(X0) + fr.uOf(X1)));
end

function pt = insideArcSample(piece)
% A point strictly INSIDE the piece, just off the midpoint of its arc. By the CCW convention every
% piece here obeys (polyConstraints' "interior is on the LEFT of each edge's direction of travel"),
% the interior is on the left of the arc walked X0 -> X1, so stepping from the arc's midpoint along
% the left normal of its tangent lands inside -- for ANY piece shape, including a two-vertex lens
% whose centroid sits on its own boundary. Used to fix a curved edge's conic SIGN, which is the one
% thing a vertex (always ON the arc, where the conic vanishes) cannot decide.
    [X0, X1] = curveEndpoints(piece);
    fr = parabolaArcFrame(piece.curveEc, 'maxQuaPar');
    u0 = fr.uOf(X0); u1 = fr.uOf(X1);
    um = 0.5*(u0+u1);
    t = fr.tangent(um);
    if u1 < u0, t = -t; end          % orient the tangent along the walk, from X0 towards X1
    t = t/norm(t);
    h = 1e-4*(1 + norm(X1-X0));
    pt = fr.point(um) + h*[-t(2), t(1)];
end

function s = asymptoticSign(diffRow, apex, dir)
% Sign of diffRow(apex+t*dir) as t->+inf: sign of the leading nonzero coefficient among the
% quadratic, linear, then constant term (exact, via the ray's own quadratic-in-t coefficients).
    [A,B,C] = quadAlongRay(diffRow, apex, dir);
    tol = 1e-9*(1+abs(A)+abs(B)+abs(C));
    if abs(A) > tol, s = sign(A)*1e6; return; end   % large sentinel: dominates any finite vals()
    if abs(B) > tol, s = sign(B)*1e6; return; end
    s = C;   % genuinely constant along this ray to available precision: use the apex value
end

function [A,B,C] = quadAlongRay(diffRow, apex, dir)
% Exact coefficients of t -> diffRow(apex+t*dir) = A t^2 + B t + C (apex, dir: 1x2 row vectors).
    Q = [diffRow(5), diffRow(6); diffRow(6), diffRow(7)];
    L = [diffRow(8); diffRow(9)]; K = diffRow(10);
    A = 0.5*(dir*Q*dir');
    B = apex*Q*dir' + L'*dir';
    C = 0.5*(apex*Q*apex') + L'*apex' + K;
end

% ============================================================================================
% ----- splitting a cell by the (degenerate) curve where f1row and f2row are equal ------------
function newPieces = splitCell(cell, f1row, f2row)
% Split cell along {f1row=f2row} and return the resulting pieces (1 or more), each tagged with the
% row that wins on it and each carrying AT MOST ONE curved edge. See the file header scoping
% caveat: this REQUIRES f1row-f2row to be a degenerate conic (full 3x3 discriminant Delta==0),
% asserted here, and requires exactly 2 boundary crossings (asserted too) -- both are theorems
% about this pipeline's own adjacent sub-pieces, not generic facts, so a violation errors loudly.
%
% A cell that ALREADY carries an arc is handled too (it used to error). One arc per face stays the
% invariant; it is maintained by SUBDIVIDING, never by giving a piece a second curve slot -- see
% splitTwoArcPiece. What makes that enough is measured, not assumed: over the named fixture plus a
% 395-quadrilateral randomized sweep, every one of the 22 curved-cell splits had the splitting
% curve cross exactly two STRAIGHT edges, and the existing arc was either untouched (19) or met
% TANGENTIALLY (3) -- never crossed. That is the tangency structure this file already documents
% (a conjugate is C1 where its pieces join, so the other operand's face boundaries are tangent to
% the parabola), so the arc always survives whole, inside exactly one of the two halves, and that
% half is the only one that ends up with two curves.
    arcPos0 = 0; arcEc0 = [];
    if pieceIsCurved(cell)
        if ~isempty(cell.dirIn)
            error('maxQuaPar:notImplemented', ...
                ['maxQuaPar:splitCell: splitting an UNBOUNDED cell that already carries an arc ' ...
                 'is not implemented (never observed: every curved cell in the sweep was ' ...
                 'bounded). The bounded case is handled.']);
        end
        arcPos0 = cell.curveAfter;
        arcEc0  = cell.curveEc;
    end
    diffRow = f1row - f2row;
    a = diffRow(5)/2; b = diffRow(6); c = diffRow(7)/2;
    d = diffRow(8); e = diffRow(9); f = diffRow(10);
    % delta=b^2-4ac decides the quadratic part's TYPE (>0 hyperbolic, ==0 parabolic, <0 elliptic).
    % Delta (the full 3x3 discriminant) decides IRREDUCIBILITY (~=0 <=> a genuine curve of that
    % type; ==0 <=> degenerates to a point, a line, or a pair of lines).
    %
    % Representability by QuaPar depends on BOTH: delta==0 (parabolic type, including its `a line`
    % degeneracy -- see QuaPar.isParabola) is ALWAYS representable regardless of Delta, since a
    % genuine non-degenerate parabola is exactly the curved edge QuaPar's Ec was built for (see
    % the isStraight/edgeEc construction below, which already handles it once this guard lets it
    % through). Only delta~=0 (hyperbolic/elliptic type) NEEDS Delta==0 (degenerating to one or two
    % real straight lines, handled via the isStraight branch below) to be representable; if
    % Delta~=0 there it is a genuine irreducible ellipse/hyperbola, which QuaPar cannot represent.
    %
    % HISTORY: an earlier version of this guard rejected whenever Delta~=0, with no delta check at
    % all -- correct for the file's original bug (see file header HISTORY, which only needed to
    % rule out hyperbolic-type false positives) but wrong in the OTHER direction: it also rejected
    % the delta==0/Delta~=0 case, i.e. a genuine, representable parabola, mislabelling it "a
    % genuine ellipse/hyperbola" and erroring on legitimate input. Found via a randomized stress
    % test across many f(x,y)=xy triangles (not just the one hard-coded 3-edge.tex example): e.g.
    % T=(0,0),(7.02,0.67),(8.43,7.63) produces a g1-face/g2-face boundary with delta~6.8e-49 (~0,
    % parabolic type) but Delta~0.063 (genuinely non-degenerate) -- a real parabola arc that the
    % rest of this function already has the machinery to build. See
    % maxQuaParTest.splitCellAcceptsGenuineNonDegenerateParabola.
    delta = b^2 - 4*a*c;
    tolDelta = 1e-6 * (abs(a) + abs(b) + abs(c))^2;
    isParabolicType = abs(delta) <= max(tolDelta, 1e-12);
    M = [a, b/2, d/2; b/2, c, e/2; d/2, e/2, f];
    sc = max(1e-9, norm(M, 'fro'));
    % det(M) is a degree-3-homogeneous function of M's entries, so a relative tolerance on those
    % entries bounds det(M) by tolRel*sc^3 (linear in tolRel), NOT (tolRel*sc)^3 (which is what an
    % earlier version of this line wrote -- cubing the tolerance itself made the threshold many
    % orders of magnitude tighter than floating-point noise, rejecting even exactly-degenerate
    % inputs; e.g. on the 3-edge.tex example, det(M)~-1.7e-17 against a true zero, but the old
    % threshold (1e-6*sc)^3~3e-18 was smaller still).
    if ~isParabolicType && abs(det(M)) > 1e-6*sc^3
        error('maxQuaPar:notDegenerate', ...
            ['maxQuaPar:splitCell: the difference of the two candidate quadratics is a genuine ' ...
             'irreducible ellipse/hyperbola (delta=%g ~= 0 and the full 3x3 discriminant is nonzero), ' ...
             'which QuaPar cannot represent. This should never happen for two conjugates of adjacent ' ...
             'sub-pieces of the same originally-nonconvex domain (see maxQuaPar.m header); the input ' ...
             'pair violates that.'], delta);
    end

    edges = cellEdgeList(cell);
    hits = struct('edge', {}, 't', {}, 'pt', {});
    for i = 1:numel(edges)
        if i == arcPos0
            % cellEdgeList reports this edge as the CHORD between the arc's endpoints, but the
            % boundary here is the ARC. Crossings of it are found on the conic itself below;
            % solving along the chord would both miss real ones and invent ones the boundary
            % does not have.
            continue
        end
        [A,B,C] = quadAlongRay(diffRow, edges(i).apex, edges(i).dir);
        for r = solveQuad(A,B,C)
            if r >= -1e-6 && r <= edges(i).tMax + 1e-6
                % Snap a root near EITHER endpoint exactly onto it (not just clamp values that
                % overshoot past the endpoint): at the singular point of a degenerate two-line
                % conic (where both lines cross), every direction through it is a double root of
                % the local quadratic, but the two solveQuad roots straddle the true value by a
                % little floating-point noise -- one lands just past tMax (already clamped below)
                % and the other just short of it (previously left un-snapped), so their computed
                % points differed by up to ~1e-7*edgeLength, just over dedupHits' tolerance, and
                % were wrongly kept as two separate hits (see maxQuaPar.m header HISTORY; observed
                % noise on the 3-edge.tex example ranged up to ~1.5e-7, hence the 1e-6 margin here
                % rather than exactly matching that noise floor).
                if abs(r) < 1e-6, r = 0; end
                if abs(r - edges(i).tMax) < 1e-6, r = edges(i).tMax; end
                t = min(max(r,0), edges(i).tMax);
                hits(end+1) = struct('edge', i, 't', t, 'pt', edges(i).apex + t*edges(i).dir); %#ok<AGROW>
            end
        end
    end
    hits = dedupHits(hits);
    if arcPos0 ~= 0
        % The splitting curve restricted to the arc's own parabola is a quartic in the frame's
        % global monotone parameter u (parabolaArcFrame.conicCoeffs) -- this is the whole of the
        % "conic-conic intersection" this case was once thought to need, and no general conic
        % solver is involved: every curved edge here is a parabola (QuaPar.assertParabolic), so
        % one univariate quartic settles it.
        %
        % Only a genuine CROSSING would cut the arc in two and take the split beyond two pieces.
        % A tangency does not: the arc stays whole, and the touch point needs no vertex because
        % the winner does not change across it. Erroring on a real crossing keeps that assumption
        % honest rather than assumed -- the sweep produced 3 tangencies and no crossing.
        [crosses, XcArc] = arcHasStrictCrossing(cell, diffRow);
        if crosses
            % The splitting curve CROSSES the arc, so the arc is cut into two sub-arcs that end up
            % in different halves. That used to be refused as "never observed" -- true only while a
            % cell could not carry the OTHER operand's arc, which arc-vs-arc clipping now makes
            % routine. It needs no new machinery: the crossing is simply a THIRD kind of boundary
            % hit, on the arc's own edge (index arcPos0, which the loop above deliberately skips),
            % and the existing two-hit split then divides the arc along with everything else.
            hits(end+1) = struct('edge', arcPos0, 't', 0, 'pt', XcArc); %#ok<AGROW>
            [~, ord] = sort([hits.edge]);      % the split below assumes e1 < e2
            hits = hits(ord);
        end
    end
    if numel(hits) == 1
        % The cell only TOUCHES {diffRow=0} at a single point -- a tangency at the degenerate
        % conic's singular point, which happens to coincide with a cell vertex shared with other
        % face-pairs elsewhere in the arrangement (see maxQuaPar.m header HISTORY). decideWinner
        % flagged this cell "undecided" only because that one vertex evaluates to a tiny nonzero
        % residual (floating-point noise around the true value 0); the cell does not actually
        % split into two regions. Resolve the winner from the centroid (always strictly interior
        % for a convex cell, hence never the touch point itself) and return the WHOLE cell intact.
        cellA = cell;
        if arcPos0 == 0
            cellA.curveAfter = 0; cellA.curveEc = [];
        end   % a curved cell keeps its own arc: the cell is returned INTACT, arc included
        if QuaPar.evalPoly(diffRow, mean(cell.V,1)) >= 0
            cellA.f = f1row;
        else
            cellA.f = f2row;
        end
        newPieces = cellA;
        return
    end
    if numel(hits) ~= 2
        error('maxQuaPar:internal', ...
            ['maxQuaPar:splitCell: expected exactly 2 boundary crossings of the splitting curve, found ' ...
             '%d (see the file header scoping caveat -- this pipeline should never need more).'], numel(hits));
    end
    e1 = hits(1).edge; e2 = hits(2).edge; X1 = hits(1).pt; X2 = hits(2).pt;

    % Is the curve connecting X1 to X2 straight or genuinely curved? {diffRow=0} is a degenerate
    % conic (checked above), so it is a parabola, a single line, or a pair of two distinct lines;
    % in every case it can only touch a cell's interior along ONE connected branch (see
    % maxQuaPar.m header HISTORY note -- verified visually for the 3-edge.tex example, where only
    % one of the two lines cuts through each cell). Rather than reconstruct that branch by
    % symbolic factorization, test it directly: diffRow restricted to the line through X1,X2 is
    % A*t^2+B*t+C in the chord parameter t in [0,1]; X1,X2 are already roots (C=diffRow(X1)~=0,
    % A+B+C=diffRow(X2)~0), so A~0 iff the WHOLE chord lies on {diffRow=0}, i.e. the connecting
    % branch is exactly this straight segment. A~=0 means a genuinely curved (parabolic) branch,
    % handled by the pre-existing general-conic construction below.
    [Achord, ~, ~] = quadAlongRay(diffRow, X1, X2 - X1);
    scChord = max(1e-9, norm(diffRow(5:10), Inf)*max(1, norm(X2-X1))^2);
    isStraight = abs(Achord) <= 1e-8*scChord;

    if isStraight
        edgeEc = zeros(1,6);
    else
        edgeEc = [0.5*diffRow(5), diffRow(6), 0.5*diffRow(7), diffRow(8), diffRow(9), diffRow(10)];
        edgeEc = edgeEc / max(abs(edgeEc));   % a conic is scale-invariant; normalize (see pushforwardQuaParDual)
    end

    nv = size(cell.V,1);
    unbounded = ~isempty(cell.dirIn);
    m = numel(edges);   % number of boundary edges = number of node-pairs

    if ~unbounded
        midIdx = mod((e1):(e2-1), nv) + 1;    % real V-indices strictly between the two crossings
        % e1<e2 always (hits are collected walking edges in increasing order), so, exactly as in
        % clipPolyHalfPlane's analogous bounded-arc split (see that function's HISTORY comment),
        % the wrapped complement e2..e1-1 needs +nv before modding or MATLAB's colon operator
        % yields an empty range instead of wrapping.
        restIdx = mod((e2):(e1-1+nv), nv) + 1;
        cellMidB = boundedPiece(X1, cell.V(midIdx,:), X2, edgeEc);
        cellRestB = boundedPiece(X2, cell.V(restIdx,:), X1, edgeEc);
        cellA = assignSide(cellMidB, diffRow, f1row, f2row);
        cellB = assignSide(cellRestB, diffRow, f1row, f2row);
        if arcPos0 == 0
            newPieces = [cellA, cellB];
            return
        end
        % The original arc survives whole (checked above) inside exactly ONE of the two halves,
        % which therefore now carries two curves: its inherited arc and the splitting curve
        % boundedPiece has just tagged as its closing edge. Restore the inherited arc -- which
        % boundedPiece knows nothing about and would otherwise silently flatten to a chord -- and
        % subdivide that half so the one-arc-per-face invariant still holds.
        [aX0, aX1] = arcEndpointsOf(cell, arcPos0);
        % When the arc was CROSSED, neither half holds the whole arc: each holds a SUB-arc running
        % from the crossing point to one of the original endpoints. Both candidates are tried, so
        % the same code covers the crossed and uncrossed cases.
        cand = {[aX0; aX1]};
        if exist('XcArc','var') && ~isempty(XcArc)
            cand = {[XcArc; aX0], [XcArc; aX1]};
        end
        newPieces = [];
        for half = [cellA, cellB]
            p = 0;
            for cIdx = 1:numel(cand)
                p = findArcPosition(half, cand{cIdx}(1,:), cand{cIdx}(2,:));
                if p ~= 0, break, end
            end
            if p == 0 || p == half.curveAfter
                newPieces = [newPieces, half]; %#ok<AGROW>
            else
                newPieces = [newPieces, splitTwoArcPiece(half, p, arcEc0)]; %#ok<AGROW>
            end
        end
        return
    end

    % Unbounded: node(i) = dirIn-marker (i==1), V(i-1) (2<=i<=nv+1), dirOut-marker (i==m+1==nv+2).
    % Edge i connects node(i),node(i+1); e1<e2 in 1..m. The two candidate pieces are the "middle"
    % (edges e1+1..e2-1, i.e. real V-indices e1..e2-1) and the "rest" (edges 1..e1-1 and e2+1..m,
    % wrapping through node 1's dirIn ray and node end's dirOut ray -- NOT reconnectable except via
    % the SAME curve, so if the "rest" piece keeps BOTH original rays it stays unbounded with the
    % curve as a middle bridge; there is no way for a QuaPar edge to be an unbounded curved ray, so
    % if the retained piece would need the curve itself to go to infinity, that's a hard limitation
    % (not expected here -- see header caveat -- and will surface via the hit-count assertion above
    % or the isDomBounded-style checks in assemblePieces).
    vMidIdx = (e1):(e2-1);
    cellMid = boundedPiece(X1, cell.V(vMidIdx,:), X2, edgeEc);

    keepBefore = cell.V(1:e1-1,:);
    keepAfter  = cell.V(e2:nv,:);
    cellRest.V = dedupConsecutive([keepBefore; X1; X2; keepAfter]);
    cellRest.dirIn = cell.dirIn; cellRest.dirOut = cell.dirOut;
    cellRest.dirInSign = cell.dirInSign; cellRest.dirOutSign = cell.dirOutSign;
    cellRest.curveAfter = size(keepBefore,1) + 1; cellRest.curveEc = edgeEc;

    % Determine which of cellMid/cellRest is "f1 wins" by evaluating f1row-f2row at one interior
    % sample of each (a point strictly on one side; the cell's own vertices adjacent to the curve,
    % other than the two crossing points, are exactly such samples -- fall back to the crossing
    % midpoint if a piece has no other vertex, e.g. a bounded middle with e2==e1+1).
    cellA = assignSide(cellMid, diffRow, f1row, f2row);
    cellB = assignSide(cellRest, diffRow, f1row, f2row);
    newPieces = [cellA, cellB];
end

function [X0, X1] = arcEndpointsOf(piece, i)
% Endpoints of the boundary edge piece.V(i) -> piece.V(i+1), wrapping for a bounded piece.
% curveEndpoints does this for a piece's OWN curveAfter; splitCell needs it for the arc index it
% saved before the piece was rebuilt.
    nv = size(piece.V,1);
    X0 = piece.V(i,:);
    X1 = piece.V(mod(i,nv)+1,:);
end

function pos = findArcPosition(piece, X0, X1)
% Index i with piece.V(i) -> piece.V(i+1) being the edge whose endpoints are X0,X1 (in either
% order), or 0 if this piece does not have that edge. Used to locate the parent cell's arc in a
% half that splitCell has just rebuilt, since boundedPiece reindexes the vertices.
    pos = 0;
    nv = size(piece.V,1);
    if nv < 2, return, end
    tol = 1e-6 * (1 + max(abs(piece.V(:))));
    if isempty(piece.dirIn), last = nv; else, last = nv-1; end
    for i = 1:last
        A = piece.V(i,:); B = piece.V(mod(i,nv)+1,:);
        if (norm(A-X0) < tol && norm(B-X1) < tol) || (norm(A-X1) < tol && norm(B-X0) < tol)
            pos = i; return
        end
    end
end

function out = splitTwoArcPiece(piece, arcPos, arcEc)
% Cut a bounded piece that carries TWO curved edges -- the inherited arc at edge arcPos (conic
% arcEc) and the splitting curve at edge piece.curveAfter -- into two pieces with ONE each, using a
% single straight chord. One arc per face is the invariant this file (and QuaPar's single Ec slot
% per edge) is built on; it is maintained by subdividing, never by widening the representation.
%
% The chord runs between the two arcs' facing endpoints, so each half keeps one whole arc and the
% chord closes it. Both diagonals that separate the arcs are tried, and the chord is ACCEPTED only
% if its midpoint is genuinely inside the piece: a face bounded by a parabola is not convex on that
% side, so a diagonal is not automatically interior. If neither works the piece is returned
% unsplit, which is not silently wrong -- it still carries its splitting curve, and the assembly's
% own arrangement check (maxQuaParTest.maxQuaParResultsAreValidArrangements) is what would catch a
% dropped arc.
    nv = size(piece.V,1);
    c  = piece.curveAfter;
    cands = [mod(arcPos, nv)+1, mod(c, nv)+1; ...      % arc end -> curve end
             arcPos,            c          ];          % arc start -> curve start
    for k = 1:size(cands,1)
        a = cands(k,1); b = cands(k,2);
        if a == b, continue, end
        chainA = cycIdx(b, a, nv);        % walk b -> ... -> a, closed by the chord a -> b
        chainB = cycIdx(a, b, nv);        % walk a -> ... -> b, closed by the chord b -> a
        if numel(chainA) < 3 || numel(chainB) < 3, continue, end
        mid = 0.5*(piece.V(a,:) + piece.V(b,:));
        if ~insideStraightHull(piece, arcPos, arcEc, mid), continue, end
        pA = subPiece(piece, chainA, arcPos, arcEc);
        pB = subPiece(piece, chainB, arcPos, arcEc);
        if isempty(pA) || isempty(pB), continue, end
        out = [pA, pB];
        return
    end
    out = piece;
end

function idx = cycIdx(from, to, nv)
% The cyclic index walk from..to inclusive, wrapping through nv.
    idx = from;
    while idx(end) ~= to
        idx(end+1) = mod(idx(end), nv) + 1; %#ok<AGROW>
        if numel(idx) > nv, break, end
    end
end

function tf = insideStraightHull(piece, arcPos, arcEc, pt)
% Is pt inside the piece? Tested as: on the inner side of every STRAIGHT boundary edge, and on the
% interior side of both conics. Both curved edges carry the sign convention assignSide/facePoly
% establish (evalConic > 0 on the piece's own interior), so the conic tests are exact; the straight
% edges use the CCW "interior is on the left" convention polyConstraints documents.
    nv = size(piece.V,1);
    tf = false;
    for i = 1:nv
        if i == arcPos || i == piece.curveAfter, continue, end
        A = piece.V(i,:); B = piece.V(mod(i,nv)+1,:);
        d = B - A;
        if (d(1)*(pt(2)-A(2)) - d(2)*(pt(1)-A(1))) < -1e-9*(1+norm(d)), return, end
    end
    if QuaPar.evalConic(arcEc, pt) <= 0, return, end
    if any(piece.curveEc ~= 0) && QuaPar.evalConic(piece.curveEc, pt) <= 0, return, end
    tf = true;
end

function p = subPiece(piece, idx, arcPos, arcEc)
% One side of the chord: the vertices piece.V(idx) in walk order, closed by the chord from the last
% back to the first. Exactly one of the parent's two curved edges survives in it, and which one
% decides the new curveAfter/curveEc.
    p = [];
    V = piece.V(idx,:);
    nvOld = size(piece.V,1);
    arcHere = [];  ecHere = [];
    for t = 1:numel(idx)-1
        e = idx(t);                                   % parent edge idx(t) -> idx(t+1)
        if mod(e, nvOld)+1 ~= idx(t+1), continue, end % not a parent edge (should not happen)
        if e == arcPos,            arcHere = t; ecHere = arcEc; end
        if e == piece.curveAfter,  arcHere = t; ecHere = piece.curveEc; end
    end
    if isempty(arcHere), return, end
    p.V = V;
    p.dirIn = []; p.dirOut = []; p.dirInSign = []; p.dirOutSign = [];
    p.curveAfter = arcHere;
    p.curveEc = ecHere;
    p.f = piece.f;
    if pieceIsCurved(p) && QuaPar.evalConic(p.curveEc, insideArcSample(p)) < 0
        p.curveEc = -p.curveEc;
    end
end

function [tf, Xc] = arcHasStrictCrossing(cell, diffRow)
% Does {diffRow=0} genuinely CROSS this cell's arc strictly between its endpoints (as opposed to
% missing it, or touching it tangentially)? The splitting curve restricted to the arc's parabola is
% a quartic in the frame's global monotone parameter u; for each root strictly inside the arc's
% u-range, the sign of diffRow just before and just after decides crossing vs tangency.
    fr = parabolaArcFrame(cell.curveEc, 'maxQuaPar');
    [X0, X1] = curveEndpoints(cell);
    u0 = fr.uOf(X0); u1 = fr.uOf(X1);
    ulo = min(u0,u1); uhi = max(u0,u1); span = uhi - ulo;
    % Same 1/2 weights on the x^2/y^2 slots that QuaPar.evalPoly applies and that splitCell's own
    % edgeEc uses -- a QuaPar row is NOT a conic row.
    q = fr.conicCoeffs([0.5*diffRow(5), diffRow(6), 0.5*diffRow(7), ...
                        diffRow(8), diffRow(9), diffRow(10)]);
    sc = max(abs(q));
    tf = false; Xc = [];
    if sc == 0, return, end
    r = roots(q/sc);
    r = real(r(abs(imag(r)) < 1e-7*(1+abs(real(r)))));
    tolU = 1e-6*(1+span);
    h = 1e-5*(1+span);
    for z = 1:numel(r)
        if r(z) <= ulo + tolU || r(z) >= uhi - tolU, continue, end
        if sign(QuaPar.evalPoly(diffRow, fr.point(r(z)-h))) ~= ...
           sign(QuaPar.evalPoly(diffRow, fr.point(r(z)+h)))
            tf = true; Xc = fr.point(r(z)); return
        end
    end
end

function piece = boundedPiece(Xstart, Vmid, Xend, ecRow)
    piece.V = dedupConsecutive([Xstart; Vmid; Xend]);
    piece.dirIn = []; piece.dirOut = [];
    piece.dirInSign = []; piece.dirOutSign = [];
    piece.curveAfter = size(piece.V,1);   % curve is the CLOSING edge (V(end) back to V(1))
    piece.curveEc = ecRow;
end

function piece = assignSide(piece, diffRow, f1row, f2row) %#ok<INUSD>
% Tag which row (f1row or f2row) wins on this piece, using any vertex strictly interior to it
% (not one of the two curve endpoints) if available, else the curve-chord midpoint. Also normalizes
% the piece's own curve orientation to the convention facePoly establishes and
% clipPolyHalfPlaneCurved relies on: evalConic(piece.curveEc,.) > 0 on the piece's OWN interior.
% (splitCell builds edgeEc once from f1row-f2row, so it is positive where f1 wins -- i.e. correct
% for one of the two pieces it splits into and inverted for the other.)
    if size(piece.V,1) > 2
        samplePt = piece.V(2,:);   % V(1)=one curve endpoint, V(end)=the other; V(2) is a real vertex
    else
        samplePt = mean(piece.V,1);
    end
    d = QuaPar.evalPoly(diffRow, samplePt);
    if d >= 0, piece.f = f1row; else, piece.f = f2row; end
    if pieceIsCurved(piece) && QuaPar.evalConic(piece.curveEc, insideArcSample(piece)) < 0
        piece.curveEc = -piece.curveEc;
    end
end

function edges = cellEdgeList(cell)
% Boundary edges of cell in walk order, as {apex, dir, tMax} triples (tMax=1 for a bounded segment
% V(i)->V(i+1), tMax=Inf for a ray). Matches the node/edge indexing used throughout this file: edge
% i connects node(i) to node(i+1).
    nv = size(cell.V,1);
    edges = struct('apex', {}, 'dir', {}, 'tMax', {});
    if isempty(cell.dirIn)
        for i = 1:nv
            j = mod(i,nv)+1;
            edges(end+1) = struct('apex', cell.V(i,:), 'dir', cell.V(j,:)-cell.V(i,:), 'tMax', 1); %#ok<AGROW>
        end
    else
        edges(end+1) = struct('apex', cell.V(1,:), 'dir', cell.dirIn, 'tMax', Inf);
        for i = 1:nv-1
            edges(end+1) = struct('apex', cell.V(i,:), 'dir', cell.V(i+1,:)-cell.V(i,:), 'tMax', 1); %#ok<AGROW>
        end
        edges(end+1) = struct('apex', cell.V(end,:), 'dir', cell.dirOut, 'tMax', Inf);
    end
end

function r = solveQuad(A,B,C)
    tol = 1e-9*(1+abs(A)+abs(B)+abs(C));
    if abs(A) <= tol
        if abs(B) <= tol, r = []; return; end
        r = -C/B; return
    end
    disc = B^2 - 4*A*C;
    if disc < -tol, r = []; return; end
    disc = max(disc,0);
    r = [(-B+sqrt(disc))/(2*A), (-B-sqrt(disc))/(2*A)];
end

function hits = dedupHits(hits)
% Merge boundary crossings that are the SAME physical point computed via two different cell
% edges' independent quadratic-root arithmetic. sqrt(eps)~1.5e-8 is too tight for this: two
% genuinely-coincident hits (e.g. the curve crossing exactly at a cell corner shared by two
% adjacent boundary edges) can disagree by ~1e-7 between the two arithmetic paths -- the same
% cross-arithmetic noise floor already documented (and handled with a 1e-6 absolute tolerance) in
% assemblePieces' global vertex merge, see its HISTORY comment. Without this, e.g. a hit pair
% ~7e-8 apart at a triangle corner was wrongly kept as 2 separate hits, inflating a genuine 2-hit
% split into 3 and tripping the "expected exactly 2 boundary crossings" assertion below on
% legitimate input (found via a randomized triangle stress test, T=(0,0),(2.11,1.43),(8.84,4.50)).
    if numel(hits) < 2, return; end
    tol = 1e-6;
    keep = true(1,numel(hits));
    for i = 1:numel(hits)
        if ~keep(i), continue; end
        for j = i+1:numel(hits)
            if keep(j) && norm(hits(i).pt-hits(j).pt) < tol, keep(j) = false; end
        end
    end
    hits = hits(keep);
end

% ============================================================================================
% ----- final reassembly: merge all (V,dirIn,dirOut,curveAfter,curveEc,f) pieces into a QuaPar -
function g = assemblePieces(pieces)
% Generalizes convEnvCPLQ.m's assembleTriangles/assembleTwoTriangles (coordinate vertex dedup +
% half-edge pairing) to: arbitrary edge counts, unbounded faces (rays, represented the same way
% QuaPar itself does -- apex plus a direction-only vertex), and at most one curved (Ec) edge per
% piece. Since every input to maxQuaPar is full-domain, every edge produced here MUST pair with
% exactly one neighbour; an unpaired edge is treated as an internal-consistency error, not a valid
% "boundary of the domain" (there is no domain boundary -- the result is finite everywhere).
%
% HISTORY: an earlier version of this function matched half-edges by first collapsing ALL pieces'
% vertices into one global list via a single coordinate-distance tolerance (a cell corner can
% arrive here computed via two different arithmetic paths -- e.g. clipByFace's intersection
% formula vs. splitCell's ray-quadratic root -- agreeing only to ~1e-5, so the tolerance had to be
% at least that loose), then matched half-edges by GLOBAL VERTEX INDEX equality. That approach is
% fundamentally unsound for near-degenerate inputs: a single piece can have two of its OWN,
% genuinely distinct corners only ~1e-5 apart (a thin sliver edge, not arithmetic noise -- a real
% small-scale feature of a fan of many pieces meeting near one point), which a tolerance loose
% enough to fix the cross-arithmetic-noise case will ALSO merge, silently making that one piece
% self-touching (one vertex touched >2 times within a single face's own boundary, later rejected by
% QuaPar.orderEdges: "expected 2 but got 3/4"). No single distance tolerance can separate these two
% cases, since the noise floor and the genuine small-feature scale overlap (both ~1e-5): merging
% two of a piece's own vertices is not the same failure as failing to notice two different pieces
% share a boundary, but a coordinate-clustering approach conflates them because both are ultimately
% "is this point equal to that point".
%
% FIX: sidestep vertex identity entirely for the purpose of finding each edge's neighbour. Instead,
% match whole HALF-EDGES directly by geometry (matchHalfEdges below): a segment matches another
% piece's segment iff BOTH of its endpoints coincide (swapped order) with the candidate's, and a ray
% matches another piece's ray iff their apexes coincide AND their (unit) directions agree -- always
% comparing across two DIFFERENT pieces, never within one piece's own edge list. Only once every
% edge has a confirmed neighbour is a global vertex numbering derived, via union-find, restricted
% to EXACTLY the vertex identifications implied by those confirmed matches (buildGlobalVertices
% below) -- so two vertices of the SAME piece are provably never unified: every union relates one
% piece's vertex to a genuinely-matched DIFFERENT piece's vertex, never two vertices of one piece to
% each other, however close together they happen to be. This is the vertex-PROVENANCE approach the
% session handoff called for, in place of reconciling near-duplicates by raw coordinate distance.
%
% HISTORY (later session): the fix above dropped the crash rate on a randomized stress test from
% ~4/800 to ~1/800; the residual ~1/800 is a genuinely AMBIGUOUS 3-way (or more) vertex cluster --
% e.g. three pieces whose edges are ALL mutually within matchHalfEdges' tolPos, at a scale
% (~1e-4) too large to be cross-arithmetic noise (tolPos's intended target) but too small to be a
% separate genuine feature -- for which matchHalfEdges' best-first greedy matching can only pair
% up 2 of the 3, unavoidably orphaning the third half-edge regardless of matching strategy (this
% is a topological tie, not a tolerance-tuning problem: full vertex PROVENANCE, tagging each edge
% with which original g1/g2 face-pair boundary produced it, would resolve it in principle, but
% turned out to be unnecessary -- see below). Diagnosed via checkOrphanHalfEdges's namesake
% investigation: for every such orphaned half-edge found (6 randomly-sampled repro triangles), its
% own two endpoints ALWAYS resolved to the very same global vertex once the OTHER, successfully
% -matched half-edges on its own piece's boundary were accounted for -- i.e. the orphaned edge is
% provably zero-length in the resolved geometry, so dropping it (emitting no edge for it, which
% buildFinalEdgesAndFaces already does for any unmatched half-edge) is exactly correct, not a
% guess. Fixed by moving the "no matching neighbour" error out of matchHalfEdges and into the new
% checkOrphanHalfEdges (called after buildGlobalVertices, so global vertex identity is available):
% it only raises the error for an orphan whose endpoints do NOT already coincide globally (a
% genuine, still-unresolved topology gap) or for an orphaned RAY (no evidence yet that rays can be
% legitimately degenerate this way, so they keep the original strict behaviour). See
% checkOrphanHalfEdges's own header for the full argument and
% maxQuaParTest.residualVertexClusterCrash* for the regression tests.
    n = numel(pieces);
    [allNV, allE, allEc] = localEdgeLists(pieces);
    HE = buildHalfEdgeList(n, allNV, allE, allEc);
    opp = matchHalfEdges(pieces, HE);
    [V, rootOf] = buildGlobalVertices(pieces, allNV, HE, opp);
    checkOrphanHalfEdges(HE, opp, rootOf);
    [V, E, Ec, F] = buildFinalEdgesAndFaces(pieces, HE, opp, V, rootOf);

    f = zeros(n,10);
    for p = 1:n, f(p,:) = pieces(p).f; end
    g = QuaPar(V, E, Ec, f, F);
end

function [allNV, allE, allEc] = localEdgeLists(pieces)
% Per-piece edge list using LOCAL real-vertex indices (piece.V rows) only. A ray edge's "b" column
% is 0 (no local vertex): unlike the old synthetic apex+dir "vertex", ray direction is looked up
% on demand from the piece's own dirIn/dirOut wherever needed (see rayDirAt), never manufactured
% as a point to be matched by coordinate distance -- see assemblePieces' HISTORY.
    n = numel(pieces);
    allNV = zeros(1,n); allE = cell(1,n); allEc = cell(1,n);
    for p = 1:n
        piece = pieces(p);
        nv = size(piece.V,1);
        allNV(p) = nv;
        if isempty(piece.dirIn)
            Ep = zeros(nv,3);
            for i = 1:nv, Ep(i,:) = [i, mod(i,nv)+1, 1]; end
            Ecp = zeros(nv,6);
        else
            ne = nv+1;
            Ep = zeros(ne,3); Ecp = zeros(ne,6);
            Ep(1,:) = [1, 0, 0];      % incoming ray: apex = local vertex 1
            for i = 1:nv-1, Ep(i+1,:) = [i, i+1, 1]; end
            Ep(ne,:) = [nv, 0, 0];    % outgoing ray: apex = local vertex nv
        end
        if piece.curveAfter > 0
            % piece.curveAfter is a VERTEX index: the arc runs piece.V(i) -> piece.V(i+1) (wrapping
            % to V(nv)->V(1) when bounded). For a BOUNDED piece that is Ep row i directly, since
            % row i is exactly the edge (i, i mod nv + 1). For an UNBOUNDED one it is row i+1: Ep's
            % first row is the incoming RAY, so the edge V(i)->V(i+1) sits one row further down.
            % (Getting this wrong put the conic on the edge PRECEDING the arc -- a latent bug, since
            % no test in this repo produced a curved piece before curved INPUTS were supported: the
            % only other producer, splitCell's isStraight==false branch, was never reached by the
            % suite. Verified by instrumenting splitCell across maxQuaParTest and cplqAdapterTest.)
            if isempty(piece.dirIn)
                Ecp(piece.curveAfter,:) = piece.curveEc;
            else
                Ecp(piece.curveAfter+1,:) = piece.curveEc;
            end
        end
        allE{p} = Ep; allEc{p} = Ecp;
    end
end

function HE = buildHalfEdgeList(n, allNV, allE, allEc) %#ok<INUSD>
    HE = struct('piece', {}, 'aLoc', {}, 'bLoc', {}, 'isSeg', {}, 'ec', {}, 'rayOut', {});
    for p = 1:n
        Ep = allE{p}; Ecp = allEc{p};
        for e = 1:size(Ep,1)
            rayOut = ~Ep(e,3) && e == size(Ep,1);   % Ep's last row is always the OUTGOING ray
            HE(end+1) = struct('piece', p, 'aLoc', Ep(e,1), 'bLoc', Ep(e,2), ...
                'isSeg', Ep(e,3), 'ec', Ecp(e,:), 'rayOut', rayOut); %#ok<AGROW>
        end
    end
end

function pt = vertexAt(pieces, p, loc)
    pt = pieces(p).V(loc,:);
end

function d = rayDirAt(pieces, he)
    if he.rayOut, d = pieces(he.piece).dirOut; else, d = pieces(he.piece).dirIn; end
    d = d/norm(d);
end

function v = raySideVector(pieces, he)
% A representative vector, from the ray's apex, pointing toward the REST of he's own piece --
% used by oppositeSides to tell a genuine twin ray (shared by two ADJACENT pieces, one on each
% side) from two DIFFERENT pieces that merely inherit the SAME physical ray from a shared
% ancestor face (see oppositeSides' header) without actually bordering each other across it.
    piece = pieces(he.piece);
    nv = size(piece.V,1);
    if nv >= 2
        if he.aLoc == 1, other = piece.V(2,:); else, other = piece.V(nv-1,:); end
        v = other - piece.V(he.aLoc,:);
    else   % a pure 2-ray cone (nv==1): use the OTHER ray's direction as the representative side
        if he.rayOut, v = piece.dirIn; else, v = piece.dirOut; end
    end
end

function tf = oppositeSides(pieces, he, he2, d)
% True iff he and he2 -- two DIFFERENT pieces' rays with the SAME apex and direction (already
% confirmed by the caller) -- are genuinely on OPPOSITE sides of that shared ray, i.e. are the
% two true neighbours meeting along it.
%
% BUGFIX: matching purely on (apex,direction) is NOT sufficient. When a g1 (or g2) face k is cut
% into several (k,l) sub-pieces by different g2 (or g1) faces l, EVERY one of those sub-pieces
% independently inherits face k's own boundary rays as part of ITS OWN boundary -- not because
% each is adjacent to some other sub-piece THERE, but simply because they all descend from the
% same parent face. The old code paired up the first two such inheritors it found as if they were
% mutual neighbours, when in fact both sit on the SAME side of the ray (their true neighbour, if
% any, is a DIFFERENT piece descending from whichever face lies on the ray's other side) --
% producing a piece pair that geometrically OVERLAPS rather than tiles, which then starves the
% genuine neighbour of a match and crashes checkOrphanHalfEdges elsewhere (see maxQuaPar.m header
% HISTORY -- diagnosed via a coverage stress test on the paper's own V=[2 1;0 0;1 0] example,
% confirming pieces matched this way overlap over a positive area, not just share a boundary).
%
% Fix: each piece's interior lies, by the CCW convention, on a definite side of its own directed
% ray; take a point representative of "the rest of the piece" from the apex (raySideVector) and
% test its side via the 2D cross product with the shared ray direction d. Two pieces are genuine
% opposite neighbours iff their representative points fall on OPPOSITE sides (cross products of
% opposite sign) -- same-sign means both inherited the ray from the same side of the original
% parent face and must NOT be paired.
    v1 = raySideVector(pieces, he);
    v2 = raySideVector(pieces, he2);
    c1 = d(1)*v1(2) - d(2)*v1(1);
    c2 = d(1)*v2(2) - d(2)*v2(1);
    tf = (c1 * c2) < 0;
end

function opp = matchHalfEdges(pieces, HE)
% Pair every half-edge with its neighbour by direct geometry (see assemblePieces' HISTORY for why
% this replaces coordinate-clustering-then-index-equality). tolPos matches two DIFFERENT pieces'
% shared vertex/apex (loose: the ~1e-5 cross-arithmetic noise floor documented above); tolDir
% matches two DIFFERENT pieces' shared ray direction (tight: directions are unit vectors, and
% genuinely different rays sharing a nearby apex must stay distinguishable).
%
% A near-degenerate cluster of many pieces meeting close together can have several tiny sliver
% edges whose geometry falls within tolPos of MORE THAN ONE candidate (the true match, at
% ~1e-13-1e-7 cross-arithmetic noise, plus one or more spurious near-misses at up to ~1e-3). A
% purely LOCAL choice -- either "first candidate found in array order" or "closest candidate for
% THIS half-edge, considered alone" -- can still misassign, because accepting a mediocre match for
% one half-edge can consume the vertex/direction that a DIFFERENT, more tightly-matching half-edge
% pair actually needed. Instead, collect every candidate PAIR globally, sort by match quality
% (ascending distance), and accept pairs greedily best-first: the truly-corresponding pairs have
% distances at the cross-arithmetic noise floor (order 1e-13-1e-7), far below any spurious
% same-tolerance-band candidate (order 1e-3), so they are accepted first and never contested by a
% worse pairing that touched one of the same half-edges.
    tolPos = 1e-3;
    tolDir = 1e-6;
    m = numel(HE);
    cand = zeros(0,3);   % [h, h2, score]
    for h = 1:m
        if HE(h).isSeg
            Ah = vertexAt(pieces, HE(h).piece, HE(h).aLoc);
            Bh = vertexAt(pieces, HE(h).piece, HE(h).bLoc);
        else
            Ah = vertexAt(pieces, HE(h).piece, HE(h).aLoc);
            dh = rayDirAt(pieces, HE(h));
        end
        for h2 = h+1:m
            if HE(h2).piece == HE(h).piece || HE(h2).isSeg ~= HE(h).isSeg, continue; end
            % A parabolic edge and a straight one can share both endpoints (the arc and the chord
            % of a two-vertex "lens" cell, or an arc and the cut edge that closed the cell it was
            % split from) while being completely different boundaries, so endpoint coincidence
            % alone must not pair them: a curved half-edge only ever matches another curved one.
            if any(HE(h2).ec ~= 0) ~= any(HE(h).ec ~= 0), continue; end
            if HE(h).isSeg
                % Segment: each piece walks it in its own CCW order, necessarily reversed between
                % two pieces sharing it -- so endpoints match SWAPPED.
                A2 = vertexAt(pieces, HE(h2).piece, HE(h2).aLoc);
                B2 = vertexAt(pieces, HE(h2).piece, HE(h2).bLoc);
                dA = norm(Ah-B2); dB = norm(Bh-A2);
                if dA < tolPos && dB < tolPos
                    cand(end+1,:) = [h, h2, max(dA,dB)]; %#ok<AGROW>
                end
            else
                % Ray: both pieces sharing one physical ray have the SAME apex and direction.
                A2 = vertexAt(pieces, HE(h2).piece, HE(h2).aLoc);
                d2 = rayDirAt(pieces, HE(h2));
                dA = norm(Ah-A2); dD = norm(dh-d2);
                if dA < tolPos && dD < tolDir && oppositeSides(pieces, HE(h), HE(h2), dh)
                    cand(end+1,:) = [h, h2, dA]; %#ok<AGROW>
                end
            end
        end
    end
    [~, ord] = sort(cand(:,3));
    cand = cand(ord,:);

    opp = zeros(1,m);
    used = false(1,m);
    for r = 1:size(cand,1)
        h = cand(r,1); h2 = cand(r,2);
        if used(h) || used(h2), continue; end
        used(h) = true; used(h2) = true;
        opp(h) = h2; opp(h2) = h;
    end
    % Any still-zero entry of opp is left for assemblePieces' checkOrphanHalfEdges to resolve,
    % once global vertex identity is available -- see its header for why some of these are
    % provably safe to drop rather than treat as errors (a genuinely ambiguous 3-way vertex
    % cluster, not fixable by this function's own local view of the candidate list alone).
end

function [root, parent] = findRoot(parent, x)
    root = x;
    while parent(root) ~= root, root = parent(root); end
    while parent(x) ~= root, nx = parent(x); parent(x) = root; x = nx; end
end

function parent = unionKeys(parent, x, y)
    [rx, parent] = findRoot(parent, x);
    [ry, parent] = findRoot(parent, y);
    if rx ~= ry, parent(rx) = ry; end
end

function [V, rootOf] = buildGlobalVertices(pieces, allNV, HE, opp)
% Union-find over REAL (piece,localVertex) slots ONLY, driven purely by the confirmed half-edge
% matches in `opp` -- never by raw coordinate clustering -- so two vertices of the SAME piece can
% never end up unified (every union relates one piece's vertex to a DIFFERENT, already-matched
% piece's vertex). See assemblePieces' HISTORY.
    n = numel(pieces);
    offset = zeros(1,n+1);
    for p = 1:n, offset(p+1) = offset(p) + allNV(p); end
    total = offset(end);
    parent = 1:total;
    key = @(p,loc) offset(p) + loc;

    for h = 1:numel(HE)
        h2 = opp(h);
        if h2 < h, continue; end   % process each matched pair once
        if HE(h).isSeg
            parent = unionKeys(parent, key(HE(h).piece,HE(h).aLoc), key(HE(h2).piece,HE(h2).bLoc));
            parent = unionKeys(parent, key(HE(h).piece,HE(h).bLoc), key(HE(h2).piece,HE(h2).aLoc));
        else
            parent = unionKeys(parent, key(HE(h).piece,HE(h).aLoc), key(HE(h2).piece,HE(h2).aLoc));
        end
    end

    rootOf = struct('offset', offset, 'parent', parent, 'globalOf', zeros(1,total));
    V = zeros(0,2);
    for p = 1:n
        for loc = 1:allNV(p)
            [r, rootOf.parent] = findRoot(rootOf.parent, key(p,loc));
            if rootOf.globalOf(r) == 0
                V(end+1,:) = pieces(p).V(loc,:); %#ok<AGROW>
                rootOf.globalOf(r) = size(V,1);
            end
        end
    end
end

function checkOrphanHalfEdges(HE, opp, rootOf)
% Every half-edge matchHalfEdges left unpaired (opp==0) is normally a genuine topology bug (the
% error below) -- EXCEPT one provably safe case: a genuinely AMBIGUOUS 3-way (or more) vertex
% cluster, where several pieces meeting near one point each independently compute a slightly
% different position (order ~1e-4, a real feature of a near-degenerate/thin input triangle -- NOT
% the ~1e-5-1e-7 cross-arithmetic noise floor tolPos is tuned for) for what is mathematically ONE
% single vertex. matchHalfEdges' greedy best-first matching can then pair up only 2 of the 3
% (or more) mutually-close half-edges -- there is no valid 1-1 pairing that covers all of them,
% however the candidates are chosen -- leaving exactly one tiny sliver edge (connecting two of
% the near-duplicate points) without a partner.
%
% That leftover edge is safe to simply DROP (emit nothing for it -- buildFinalEdgesAndFaces
% already does this for any h with opp(h)==0) precisely WHEN its own two endpoints have already
% been identified as the SAME global vertex via the OTHER, independently-confirmed half-edge
% matches elsewhere on this piece's own boundary: that means the "gap" this edge would have
% bridged has zero length once the surrounding topology is resolved, so omitting it changes
% nothing geometrically -- the piece's boundary still closes properly through the shared vertex.
% If the two endpoints do NOT resolve to the same global vertex, this is still a genuine
% unresolved topology gap and raises the original error.
%
% Verified on 6 randomly-found reproduction triangles (all near-degenerate/thin, e.g.
% T=(8.5697,2.6142),(5.0151,1.8051),(1.3296,0.9185)) that used to throw maxQuaPar:internal here --
% in every case the orphaned edge's two endpoints resolved to one shared global vertex, and
% dropping it produced a QuaPar matching ground truth at all sample points (see
% maxQuaParTest.residualVertexClusterCrash*). See also matchHalfEdges' HISTORY.
    for h = 1:numel(HE)
        if opp(h) ~= 0, continue; end
        if ~HE(h).isSeg
            error('maxQuaPar:internal', ...
                ['assemblePieces: a boundary ray of piece %d has no matching neighbour -- inputs ' ...
                 'should be full-domain (finite everywhere), so every edge must pair with exactly ' ...
                 'one other.'], HE(h).piece);
        end
        gA = globalVertexIndex(rootOf, HE(h).piece, HE(h).aLoc);
        gB = globalVertexIndex(rootOf, HE(h).piece, HE(h).bLoc);
        if gA ~= gB
            error('maxQuaPar:internal', ...
                ['assemblePieces: a boundary edge of piece %d has no matching neighbour -- inputs ' ...
                 'should be full-domain (finite everywhere), so every edge must pair with exactly ' ...
                 'one other.'], HE(h).piece);
        end
        % else: gA==gB -- a zero-length orphan edge, safe to drop (see header).
    end
end

function gIdx = globalVertexIndex(rootOf, p, loc)
    key = rootOf.offset(p) + loc;
    [r, rootOf.parent] = findRoot(rootOf.parent, key); %#ok<NASGU>
    gIdx = rootOf.globalOf(r);
end

function [V, E, Ec, F] = buildFinalEdgesAndFaces(pieces, HE, opp, V, rootOf)
    E = zeros(0,3); Ec = zeros(0,6); F = zeros(0,2);
    for h = 1:numel(HE)
        h2 = opp(h);
        if h2 < h, continue; end   % emit each matched pair once
        aG = globalVertexIndex(rootOf, HE(h).piece, HE(h).aLoc);
        if HE(h).isSeg
            bG = globalVertexIndex(rootOf, HE(h).piece, HE(h).bLoc);
        else
            % Ray direction marker: never shared/merged across pieces (only used locally to encode
            % direction via V(b,:)-V(a,:), see assemblePieces' HISTORY), so give it its own row.
            V(end+1,:) = V(aG,:) + rayDirAt(pieces, HE(h)); %#ok<AGROW>
            bG = size(V,1);
        end
        E(end+1,:) = [aG, bG, HE(h).isSeg]; %#ok<AGROW>
        ecRow = HE(h).ec; if all(ecRow==0), ecRow = HE(h2).ec; end
        if any(ecRow ~= 0)
            % QuaPar's orientation invariant requires evalConic(Ec(j,:),.) > 0 on the LEFT of the
            % stored directed edge V(E(j,1))->V(E(j,2)). Each piece normalizes its own curveEc to
            % be positive on its OWN interior (facePoly / assignSide), and HE(h).piece is on the
            % left of (aG,bG) by construction, so this is normally already correct -- but ecRow may
            % have been taken from HE(h2) (the RIGHT-hand piece) when HE(h)'s own row is zero, so
            % re-establish it here rather than assume. The sample must be strictly INSIDE the piece
            % and off the arc: a vertex will not do, since a curved piece's vertices can all lie ON
            % the arc, where the conic vanishes and its sign carries no information.
            if any(pieces(HE(h).piece).curveEc ~= 0)
                samplePt = insideArcSample(pieces(HE(h).piece));
            else
                Vp = pieces(HE(h).piece).V;
                if size(Vp,1) > 2, samplePt = Vp(2,:); else, samplePt = mean(Vp,1); end
            end
            if QuaPar.evalConic(ecRow, samplePt) < 0
                ecRow = -ecRow;
            end
        end
        Ec(end+1,:) = ecRow; %#ok<AGROW>
        if HE(h).isSeg
            F(end+1,:) = [HE(h).piece, HE(h2).piece]; %#ok<AGROW>
        else
            % Ray: both pieces sharing one physical ray encode it apex-first, carrying no
            % left/right information via processing order -- derive it instead from which end each
            % piece uses: the piece for which the ray is OUTGOING (walked apex->direction) is on the
            % left, matching segments; the piece for which it is INCOMING is on the right.
            if HE(h).rayOut
                F(end+1,:) = [HE(h).piece, HE(h2).piece]; %#ok<AGROW>
            else
                F(end+1,:) = [HE(h2).piece, HE(h).piece]; %#ok<AGROW>
            end
        end
    end
end
