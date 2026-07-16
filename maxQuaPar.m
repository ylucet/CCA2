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
%   * TODO: g1 or g2 with a curved (parabolic) input edge (conjBilinearXYoneCE /
%     conjIndefiniteQuadTriangle-with-1-convex-edge output) -- needs face-vs-conic clipping
%     (clipByConicSide) and possibly conic-conic intersection; not implemented, errors clearly.

    if (~isempty(g1.Ec) && any(g1.Ec(:)~=0)) || (~isempty(g2.Ec) && any(g2.Ec(:)~=0))
        error('maxQuaPar:notImplemented', ...
            ['maxQuaPar currently requires both inputs to have purely polyhedral domains (no ' ...
             'curved/parabolic edges); at least one input has a nonzero Ec row.']);
    end
    assertFullDomain(g1, 'g1');
    assertFullDomain(g2, 'g2');

    pieces = struct('V', {}, 'dirIn', {}, 'dirOut', {}, 'curveAfter', {}, 'curveEc', {}, 'f', {});
    for k = 1:g1.nf
        polyK = facePoly(g1, k);
        for l = 1:g2.nf
            polyL = facePoly(g2, l);
            cell = clipByFace(polyK, polyL);
            if isempty(cell), continue; end
            f1row = g1.f(k,:); f2row = g2.f(l,:);
            [decided, winRow] = decideWinner(cell, f1row, f2row);
            if decided
                cell.curveAfter = 0; cell.curveEc = []; cell.f = winRow;
                pieces(end+1) = cell; %#ok<AGROW>
                continue
            end
            [cellA, cellB] = splitCell(cell, f1row, f2row);
            if ~isempty(cellA), pieces(end+1) = cellA; end %#ok<AGROW>
            if ~isempty(cellB), pieces(end+1) = cellB; end %#ok<AGROW>
        end
    end
    pieces = dedupPieces(pieces);
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
    if obj.nv == 0   % bare full-domain quadratic: the "face" IS all of R^2, no boundary at all
        poly.V = zeros(0,2); poly.dirIn = []; poly.dirOut = [];
        return
    end
    Pk = obj.P{k};
    n = numel(Pk);
    V = zeros(0,2); dirIn = []; dirOut = [];
    for t = 1:n
        j = abs(Pk(t)); s = sign(Pk(t));
        a = obj.E(j,1); b = obj.E(j,2);
        if obj.E(j,3) == 0   % ray: column 1 is always the finite apex
            apex = obj.V(a,:); d = obj.V(b,:) - apex; d = d/norm(d);
            if t == 1
                V(end+1,:) = apex; dirIn = d; %#ok<AGROW>
            else   % t == n, the only other ray: apex already appended, never re-add it
                dirOut = d;
            end
        else
            if s > 0, V(end+1,:) = obj.V(b,:); else, V(end+1,:) = obj.V(a,:); end %#ok<AGROW>
        end
    end
    % obj.P{k} is CLOCKWISE; reverse to CCW, which swaps the roles of "first ray" (incoming when
    % reversed becomes outgoing) and "last ray" -- the ray direction VECTORS themselves (apex to
    % infinity) don't change, only which end of the (now-reversed) vertex list they attach to.
    poly.V = flipud(V);
    poly.dirIn = dirOut;    % swapped, see comment above
    poly.dirOut = dirIn;
end

% ============================================================================================
% ----- clipping one convex poly by every boundary constraint of another (face intersection) --
function cell = clipByFace(polyK, polyL)
% polyK intersected with the convex region bounded by polyL, by clipping polyK against every
% boundary half-plane of polyL in turn (Sutherland-Hodgman style). Returns [] if the result is
% empty or degenerates below a proper 2D cell.
    cell = polyK;
    cons = polyConstraints(polyL);
    for i = 1:size(cons,1)
        cell = clipPolyHalfPlane(cell, cons(i,1:2), cons(i,3));
        if isempty(cell), return; end
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
    if size(cell.V,1) < 1 || (isempty(cell.dirIn) && size(cell.V,1) < 3)
        cell = []; return
    end
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
    if isempty(pts), return; end
    tol = 1e-7;
    tolSnap = 1e-4;
    for pi = 1:size(pts,1)
        p = pts(pi,:);
        again = true;
        while again
            again = false;
            nv = size(cell.V,1);
            if nv == 0 || any(all(abs(cell.V - p) < tolSnap, 2)), break; end
            if isempty(cell.dirIn)
                for i = 1:nv
                    j = mod(i,nv)+1;
                    if onOpenSegment(cell.V(i,:), cell.V(j,:), p, tol)
                        cell.V = [cell.V(1:i,:); p; cell.V(i+1:end,:)];
                        again = true; break
                    end
                end
            else
                for i = 1:nv-1
                    if onOpenSegment(cell.V(i,:), cell.V(i+1,:), p, tol)
                        cell.V = [cell.V(1:i,:); p; cell.V(i+1:end,:)];
                        again = true; break
                    end
                end
                if ~again && onOpenRay(cell.V(1,:), cell.dirIn, p, tol)
                    cell.V = [p; cell.V]; again = true;
                elseif ~again && onOpenRay(cell.V(end,:), cell.dirOut, p, tol)
                    cell.V = [cell.V; p]; again = true;
                end
            end
        end
    end
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
        jn = mod(i,nv) + 1;
        d = poly.V(jn,:) - poly.V(i,:);
        n = [d(2), -d(1)];
        cons(end+1,:) = [n, n*poly.V(i,:)']; %#ok<AGROW>
    end
    if ~isempty(poly.dirIn)
        % travel direction into V(1,:) is -dirIn (arriving from infinity); outward normal =
        % rot90cw(-dirIn) = (-dirIn(2), dirIn(1)).
        n = [-poly.dirIn(2), poly.dirIn(1)];
        cons(end+1,:) = [n, n*poly.V(1,:)']; %#ok<AGROW>
        % travel direction out of V(end,:) is +dirOut; outward normal = rot90cw(dirOut).
        n = [poly.dirOut(2), -poly.dirOut(1)];
        cons(end+1,:) = [n, n*poly.V(end,:)']; %#ok<AGROW>
    end
end

function s = sign2(v, tol)
    s = zeros(size(v));
    s(v > tol) = 1; s(v < -tol) = -1;
end

function poly2 = clipPolyHalfPlane(poly, nrm, c)
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
        else            % both ray ends outside: keep only the middle; result becomes bounded
            Vnew = dedupConsecutive([X1; poly.V(vIdxMid,:); X2]);
            if size(Vnew,1) < 3, poly2 = []; return; end
            poly2.V = Vnew; poly2.dirIn = []; poly2.dirOut = [];
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
    p = cross(1); X = xpt(p);
    if st(1) <= 0   % dirIn side inside, dirOut side outside: keep 1..p, replace dirOut
        keepV = poly.V(1:min(p-1,nv),:);
        poly2.V = dedupConsecutive([keepV; X]);
        poly2.dirIn = poly.dirIn;
        poly2.dirOut = [-nrm(2), nrm(1)];
    else            % dirOut side inside, dirIn side outside: keep p+1..m, replace dirIn
        keepV = poly.V(max(p,1):nv,:);
        poly2.V = dedupConsecutive([X; keepV]);
        poly2.dirIn = [nrm(2), -nrm(1)];
        poly2.dirOut = poly.dirOut;
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
function [cellA, cellB] = splitCell(cell, f1row, f2row)
% Split cell into cellA (f1row wins) and cellB (f2row wins) along {f1row=f2row}. See the file
% header scoping caveat: this REQUIRES f1row-f2row to be a degenerate conic (full 3x3 discriminant
% Delta==0), asserted here, and requires exactly 2 boundary crossings (asserted too) -- both are
% theorems about this pipeline's own adjacent sub-pieces, not generic facts, so a violation errors
% loudly.
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
    if numel(hits) == 1
        % The cell only TOUCHES {diffRow=0} at a single point -- a tangency at the degenerate
        % conic's singular point, which happens to coincide with a cell vertex shared with other
        % face-pairs elsewhere in the arrangement (see maxQuaPar.m header HISTORY). decideWinner
        % flagged this cell "undecided" only because that one vertex evaluates to a tiny nonzero
        % residual (floating-point noise around the true value 0); the cell does not actually
        % split into two regions. Resolve the winner from the centroid (always strictly interior
        % for a convex cell, hence never the touch point itself) and return the WHOLE cell intact.
        cellA = cell; cellA.curveAfter = 0; cellA.curveEc = [];
        if QuaPar.evalPoly(diffRow, mean(cell.V,1)) >= 0
            cellA.f = f1row;
        else
            cellA.f = f2row;
        end
        cellB = [];
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
    cellRest.curveAfter = size(keepBefore,1) + 1; cellRest.curveEc = edgeEc;

    % Determine which of cellMid/cellRest is "f1 wins" by evaluating f1row-f2row at one interior
    % sample of each (a point strictly on one side; the cell's own vertices adjacent to the curve,
    % other than the two crossing points, are exactly such samples -- fall back to the crossing
    % midpoint if a piece has no other vertex, e.g. a bounded middle with e2==e1+1).
    cellA = assignSide(cellMid, diffRow, f1row, f2row);
    cellB = assignSide(cellRest, diffRow, f1row, f2row);
end

function piece = boundedPiece(Xstart, Vmid, Xend, ecRow)
    piece.V = dedupConsecutive([Xstart; Vmid; Xend]);
    piece.dirIn = []; piece.dirOut = [];
    piece.curveAfter = size(piece.V,1);   % curve is the CLOSING edge (V(end) back to V(1))
    piece.curveEc = ecRow;
end

function piece = assignSide(piece, diffRow, f1row, f2row) %#ok<INUSD>
% Tag which row (f1row or f2row) wins on this piece, using any vertex strictly interior to it
% (not one of the two curve endpoints) if available, else the curve-chord midpoint.
    if size(piece.V,1) > 2
        samplePt = piece.V(2,:);   % V(1)=one curve endpoint, V(end)=the other; V(2) is a real vertex
    else
        samplePt = mean(piece.V,1);
    end
    d = QuaPar.evalPoly(diffRow, samplePt);
    if d >= 0, piece.f = f1row; else, piece.f = f2row; end
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
            Ecp(piece.curveAfter,:) = piece.curveEc;
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
                if dA < tolPos && dD < tolDir
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
            % stored directed edge V(E(j,1))->V(E(j,2)). A curved edge's Ec row is built once in
            % splitCell from f1row-f2row and reused unchanged for both neighbouring pieces, so its
            % sign is tied to "f1 wins" globally, not to whichever piece ends up on the left of
            % (aG,bG) here -- flip if it doesn't evaluate positive on HE(h).piece's own interior
            % (which is always on the left of (aG,bG) by construction).
            Vp = pieces(HE(h).piece).V;
            if size(Vp,1) > 2, samplePt = Vp(2,:); else, samplePt = mean(Vp,1); end
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
