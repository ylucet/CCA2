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
    g = assemblePieces(pieces);
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
    if size(cell.V,1) < 1 || (isempty(cell.dirIn) && size(cell.V,1) < 3)
        cell = []; return
    end
end

function cons = polyConstraints(poly)
% Outward half-plane constraints [n1 n2 c] (inside iff n*x'<=c) for every boundary edge of poly,
% in CCW order: interior is on the LEFT of each edge's direction of travel, so the outward normal
% is rot90cw(direction) = (dy,-dx) (same convention as conjPieceCPLQ.m's unitNormal).
    nv = size(poly.V,1);
    cons = zeros(0,3);
    for i = 1:nv-1
        d = poly.V(i+1,:) - poly.V(i,:);
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
        % with a chord between the two crossing points either way.
        midIdx = mod((p1):(p2-1), nv) + 1;
        if isempty(midIdx) || all(st(midIdx) <= 0)
            keepIdx = midIdx;
        else
            keepIdx = mod((p2):(p1-1), nv) + 1;
        end
        Vnew = dedupConsecutive([X1; poly.V(keepIdx,:); X2]);
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
    % Delta is the standard conic invariant that decides irreducibility (Delta~=0 <=> a genuine
    % curve of type delta=b^2-4ac; Delta==0 <=> degenerates to a point, a line, or a pair of
    % lines). Checking only delta (as an earlier version of this file did, via rank of the 2x2
    % quadratic-part submatrix [a b/2;b/2 c]) cannot tell a genuine hyperbola/ellipse apart from
    % this last, still-representable case -- see the file header HISTORY note.
    M = [a, b/2, d/2; b/2, c, e/2; d/2, e/2, f];
    sc = max(1e-9, norm(M, 'fro'));
    % det(M) is a degree-3-homogeneous function of M's entries, so a relative tolerance on those
    % entries bounds det(M) by tolRel*sc^3 (linear in tolRel), NOT (tolRel*sc)^3 (which is what an
    % earlier version of this line wrote -- cubing the tolerance itself made the threshold many
    % orders of magnitude tighter than floating-point noise, rejecting even exactly-degenerate
    % inputs; e.g. on the 3-edge.tex example, det(M)~-1.7e-17 against a true zero, but the old
    % threshold (1e-6*sc)^3~3e-18 was smaller still).
    if abs(det(M)) > 1e-6*sc^3
        error('maxQuaPar:notDegenerate', ...
            ['maxQuaPar:splitCell: the difference of the two candidate quadratics is not (numerically) ' ...
             'a degenerate conic (its full 3x3 discriminant is nonzero) -- a genuine ellipse/hyperbola ' ...
             'boundary would be needed, which QuaPar cannot represent. This should never happen for two ' ...
             'conjugates of adjacent sub-pieces of the same originally-nonconvex domain (see maxQuaPar.m ' ...
             'header); the input pair violates that.']);
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
        restIdx = mod((e2):(e1-1), nv) + 1;
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
    if numel(hits) < 2, return; end
    tol = sqrt(eps);
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
% Global vertex merge tolerance: an absolute 1e-6, not sqrt(eps)~1.5e-8 -- the SAME physical
% vertex can arrive here computed via two different arithmetic paths (e.g. a cell corner from
% clipByFace's intersection formula vs. the same point recomputed by splitCell's ray-quadratic
% root, snapped to an edge endpoint) that agree only to ~1e-7, which sqrt(eps) is too tight to
% treat as equal, silently leaving two "different" vertices whose edges then fail to pair up in
% the half-edge matching below (see maxQuaPar.m header HISTORY).
    tol = 1e-6;
    n = numel(pieces);
    allV = cell(1,n); allE = cell(1,n); allEc = cell(1,n);
    for p = 1:n
        piece = pieces(p);
        nv = size(piece.V,1);
        if isempty(piece.dirIn)
            Vp = piece.V;
            Ep = zeros(nv,3); Ecp = zeros(nv,6);
            for i = 1:nv
                Ep(i,:) = [i, mod(i,nv)+1, 1];
            end
        else
            Vp = [piece.V; piece.V(1,:)+piece.dirIn; piece.V(end,:)+piece.dirOut];
            ne = nv+1;
            Ep = zeros(ne,3); Ecp = zeros(ne,6);
            Ep(1,:) = [1, nv+1, 0];
            for i = 1:nv-1
                Ep(i+1,:) = [i, i+1, 1];
            end
            Ep(ne,:) = [nv, nv+2, 0];
        end
        if piece.curveAfter > 0
            Ecp(piece.curveAfter,:) = piece.curveEc;
        end
        allV{p} = Vp; allE{p} = Ep; allEc{p} = Ecp;
    end

    V = zeros(0,2);
    globalIdx = cell(1,n);
    for p = 1:n
        Vp = allV{p}; gi = zeros(size(Vp,1),1);
        for r = 1:size(Vp,1)
            hit = 0;
            for q = 1:size(V,1)
                if norm(V(q,:)-Vp(r,:)) < tol, gi(r) = q; hit = 1; break; end
            end
            if ~hit, V(end+1,:) = Vp(r,:); gi(r) = size(V,1); end %#ok<AGROW>
        end
        globalIdx{p} = gi;
    end

    HE = struct('piece', {}, 'a', {}, 'b', {}, 'isSeg', {}, 'ec', {});
    for p = 1:n
        Ep = allE{p}; Ecp = allEc{p}; gi = globalIdx{p};
        for e = 1:size(Ep,1)
            HE(end+1) = struct('piece', p, 'a', gi(Ep(e,1)), 'b', gi(Ep(e,2)), ...
                'isSeg', Ep(e,3), 'ec', Ecp(e,:)); %#ok<AGROW>
        end
    end
    haVec = [HE.a]; hbVec = [HE.b]; hsVec = [HE.isSeg];
    used = false(1,numel(HE));
    E = zeros(0,3); Ec = zeros(0,6); F = zeros(0,2);
    for h = 1:numel(HE)
        if used(h), continue; end
        used(h) = true;
        opp = find(~used & hbVec==HE(h).a & haVec==HE(h).b & hsVec==HE(h).isSeg, 1);
        if isempty(opp)
            error('maxQuaPar:internal', ...
                ['assemblePieces: boundary edge (%d,%d) has no matching neighbour -- inputs should ' ...
                 'be full-domain (finite everywhere), so every edge must pair with exactly one other.'], ...
                HE(h).a, HE(h).b);
        end
        used(opp) = true;
        E(end+1,:) = [HE(h).a, HE(h).b, HE(h).isSeg]; %#ok<AGROW>
        ecRow = HE(h).ec; if all(ecRow==0), ecRow = HE(opp).ec; end
        Ec(end+1,:) = ecRow; %#ok<AGROW>
        F(end+1,:) = [HE(h).piece, HE(opp).piece]; %#ok<AGROW>
    end

    f = zeros(n,10);
    for p = 1:n, f(p,:) = pieces(p).f; end
    g = QuaPar(V, E, Ec, f, F);
end
