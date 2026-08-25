function g = conjConvexPolygon(W, dFirst, dLast, A, L, c)
% conjConvexPolygon  The conjugate of a CONVEX quadratic over ANY convex polygon, bounded or not,
%   in closed form and with no symbolic call and no triangulation.
%
% objective: remove the two reasons `conj` still reaches the symbolic Case C. Both are structural
%   rather than mathematical, and both disappear here:
%
%   1. AN UNBOUNDED DOMAIN went to Case C unconditionally, because every numeric route below
%      conjCPLQ needs a bounded TRIANGLE. A recession direction changes one branch of this
%      construction and nothing else, so the gate is not needed.
%   2. A POLYGON was fan-triangulated and the triangles' conjugates folded by Step 3. An
%      indefinite triangle conjugates to a PARABOLIC QuaPar, and folding two of those can produce
%      a cell whose arc is cut twice by a straight clip -- `maxQuaPar:notImplemented`,
%      clipPolyHalfPlaneCurved. A CONVEX piece never needed the triangulation at all: its
%      conjugate is polyhedral, so nothing curved ever enters Step 3 for it. ALGORITHM.md says
%      exactly this ("A convex piece never needed a triangle; splitting it only forces Step 3 to
%      glue back together what was never broken") and this is that sentence implemented.
%
% [input]  W      : n x 2 vertices of the polygon in CCW order (interior on the LEFT of the walk)
%          dFirst : 1 x 2 direction of the ray leaving W(1,:) to infinity, or [] if bounded
%          dLast  : 1 x 2 direction of the ray leaving W(end,:) to infinity, or [] if bounded
%          A, L, c: q(x) = 1/2 x'A x + L'x + c with A POSITIVE DEFINITE
% [output] g      : QuaPol -- the subdivision is POLYHEDRAL and every face function is a quadratic
%
% WHY THE RESULT IS A QuaPol AND NOT A QuaPar, which is the whole point. For a convex q the sup is
% attained at the unconstrained optimum when that is feasible and on the boundary otherwise, so
% the cells are indexed by the ACTIVE SET and every boundary between them is the affine condition
% "this multiplier changes sign". CONJ_FIELD_PROOF.md 7.3 states it as a theorem: the subdivision
% is polyhedral and rational, with no parabola and no surd anywhere. Parabolas begin with
% INDEFINITE pieces and with the cross-piece max, neither of which is here.
%
% THE CONSTRUCTION, one cell per face of the polygon's own face lattice:
%
%   F0   the interior branch. x* = A^{-1}(s - L) is feasible, and f* = 1/2 (s-L)'A^{-1}(s-L) - c.
%        Its cell is exactly A*P + L, the affine image of the polygon: x* in P iff s in A*P + L.
%   S_k  one per EDGE k. The maximiser is in the relative interior of edge k, so it maximises
%        along the LINE through the edge, and f*(s) = <s,v> - q(v) + (<s,e> - g0)^2 / (2 e'A e)
%        with v the edge's base vertex, e its direction and g0 = <grad q(v), e>.
%   C_i  one per VERTEX. The maximiser is v_i, so f*(s) = <s,v_i> - q(v_i) -- affine, on the
%        vertex's normal cone translated to the dual point s_i = grad q(v_i).
%
% The dual mesh is those cells glued exactly as the polygon's own boundary is: S_k meets F0 along
% the segment s_k s_{k+1}, and meets C_k and C_{k+1} along the two rays from those dual points in
% the direction of edge k's OUTWARD NORMAL. A ray edge of the polygon is the same construction
% with its far dual point at infinity, so S_k becomes a CONE at the one finite dual point it has,
% spanned by that normal and by A*d -- the image of the polygon's own recession direction.
%
% For n = 3 and a bounded polygon this is exactly conjPieceCPLQ's own conjConvexQuadTriangle, and
% conjConvexPolygonTest checks the two agree rather than assuming it.

    if nargin < 6, error('conjConvexPolygon:nargin', 'expects (W, dFirst, dLast, A, L, c).'); end
    W = double(W);  A = double(A);  L = double(L(:));  c = double(c);
    n = size(W,1);
    if size(W,2) ~= 2 || n < 1
        error('conjConvexPolygon:shape', 'W must be n x 2 with n >= 1; got %dx%d.', n, size(W,2));
    end
    A = (A + A.')/2;
    ev = eig(A);
    if min(ev) <= sqrt(eps) * max(1, max(ev))
        error('conjConvexPolygon:notPositiveDefinite', ...
            ['A must be POSITIVE DEFINITE (eigenvalues [%g %g]). A semidefinite or indefinite ' ...
             'piece has a different cell structure -- the sup can be infinite in a recession ' ...
             'direction -- and is handled elsewhere.'], ev(1), ev(2));
    end
    unbounded = ~isempty(dFirst) || ~isempty(dLast);
    if unbounded && (isempty(dFirst) || isempty(dLast))
        error('conjConvexPolygon:oneRay', ...
            ['an unbounded convex polygon has TWO rays; got one. A single ray bounds a ' ...
             'half-plane-like set whose boundary walk is not closed and which this ' ...
             'construction cannot orient.']);
    end

    % ---- the polygon's edges in walk order, with OUTWARD normals -----------------------------
    % Bounded: edges 1..n, edge k from W(k) to W(k+1 mod n).
    % Unbounded: edges 0..n. Edge 0 is the ray AT W(1) going to infinity along dFirst; as a
    % directed boundary element in walk order it is traversed INWARD, so its direction for the
    % normal is -dFirst. Edge n is the ray at W(n) leaving along dLast.
    [base, dir, isRay, sdir] = edgeList(W, dFirst, dLast, unbounded);
    m = numel(base);                              % number of polygon edges
    nrm = zeros(m, 2);
    for k = 1:m
        nrm(k,:) = outwardNormal(dir(k,:));
    end
    assertOrientation(W, dFirst, dLast, base, dir, nrm);

    % ---- dual points and the geometry that hangs off them -------------------------------------
    S = (A*W.' + L).';                            % s_i = grad q(v_i), one per polygon vertex
    scale = dualScale(S, A, dFirst, dLast);       % a length for placing ray marker points

    % ---- assemble the mesh ---------------------------------------------------------------------
    % Edges are declared with the PAIR of faces they separate, unoriented, and `finish` decides
    % which is left and which is right from each face's own representative interior point. Doing
    % it by reasoning instead was tried and got two of the six cases backwards -- the two rays of
    % one strip have OPPOSITE orientations relative to that strip, and the polygon's two recession
    % rays differ from each other as well, because one is walked inward and the other outward.
    % A representative point cannot be wrong about which side it is on.
    B = meshBuilder(scale);
    idxS = zeros(n,1);
    for i = 1:n, idxS(i) = B.vertex(S(i,:)); end

    F0 = B.face();                                % the interior branch
    Sf = zeros(m,1);  for k = 1:m, Sf(k) = B.face(); end

    % A "vertex" whose two edges are PARALLEL is not a corner of the polygon, and its normal cone
    % is a single ray of zero width. Building a face for it produces an empty cell that eval can
    % never reach, so it is dropped and the two strips meeting there become neighbours instead.
    % This is not an edge case to tolerate but the ONLY way to describe some real domains: a
    % HALF-PLANE has no corner at all, and is written here as one marker point with two opposite
    % rays. The two strips then carry the SAME function -- stripFun depends on the edge direction
    % only through e*e' and through g0*e, both even in e -- which is the arithmetic saying the
    % same thing.
    isCorner = true(n,1);
    for i = 1:n
        [kp, kn] = incidentEdges(i, n, m, unbounded);
        if norm(nrm(kp,:) - nrm(kn,:)) <= 1e-12 * (1 + norm(nrm(kp,:)))
            isCorner(i) = false;
        end
    end
    Cf = zeros(n,1);
    for i = 1:n
        if isCorner(i), Cf(i) = B.face(); end
    end

    pend = zeros(0,3);        % [dualVertexIdx, edgeIdx, vertexIdx] for edges at a DROPPED corner
    for k = 1:m
        if isRay(k)
            % A RAY edge of the polygon. Its strip cell is a CONE at the one finite dual point it
            % has: bounded by the ray along A*d -- which is F0's own recession direction, since
            % F0 = A*P + L -- and by the ray along the edge's outward normal. `sdir`, not `dir`:
            % F0 recedes along the image of the polygon's own recession direction, which points
            % AWAY from the base vertex on both rays.
            i = base(k);
            B.ray(idxS(i), (A * sdir(k,:).').', F0, Sf(k));
            if isCorner(i), B.ray(idxS(i), nrm(k,:), Sf(k), Cf(i)); else, pend(end+1,:) = [idxS(i), k, i]; end %#ok<AGROW>
        else
            i = base(k); j = nextVertex(i, n);
            B.segment(idxS(i), idxS(j), F0, Sf(k));
            if isCorner(i), B.ray(idxS(i), nrm(k,:), Sf(k), Cf(i)); else, pend(end+1,:) = [idxS(i), k, i]; end %#ok<AGROW>
            if isCorner(j), B.ray(idxS(j), nrm(k,:), Sf(k), Cf(j)); else, pend(end+1,:) = [idxS(j), k, j]; end %#ok<AGROW>
        end
    end

    % Each dropped corner leaves exactly TWO pending edges -- one per incident strip -- and they
    % are the same ray, so they become ONE edge separating the two strips.
    for i = find(~isCorner).'
        rows = find(pend(:,3) == i);
        if numel(rows) ~= 2
            error('conjConvexPolygon:internal', ...
                'dropped corner %d left %d pending edges; it must leave exactly 2.', i, numel(rows));
        end
        k1 = pend(rows(1),2); k2 = pend(rows(2),2);
        B.ray(idxS(i), nrm(k1,:), Sf(k1), Sf(k2));
    end

    % ---- one representative interior point per face, for the orientation ----------------------
    rep = zeros(B.count(), 2);
    rep(F0,:) = mean(S, 1);
    if unbounded
        rep(F0,:) = rep(F0,:) + scale * unitOr0(((A*(dFirst(:) + dLast(:))).'));
    end
    for k = 1:m
        if isRay(k)
            i = base(k);
            rep(Sf(k),:) = S(i,:) + scale * unitOr0(unitOr0((A*sdir(k,:).').') + nrm(k,:));
        else
            i = base(k); j = nextVertex(i, n);
            rep(Sf(k),:) = 0.5*(S(i,:) + S(j,:)) + 0.5*scale*nrm(k,:);
        end
    end
    for i = 1:n
        if ~isCorner(i), continue, end
        [kp, kn] = incidentEdges(i, n, m, unbounded);
        rep(Cf(i),:) = S(i,:) + 0.5*scale * unitOr0(nrm(kp,:) + nrm(kn,:));
    end

    % ---- the face functions -----------------------------------------------------------------
    f = zeros(B.count(), 6);
    f(F0,:) = interiorFun(A, L, c);
    for k = 1:m
        % sdir, not dir: the boundary WALK of the first ray goes inward (from infinity to W(1)),
        % but the strip's own parametrisation must run OUT along the ray from its base vertex.
        % Using the walk direction there flips the sign of g0 and silently mirrors the cell.
        f(Sf(k),:) = stripFun(W(base(k),:).', sdir(k,:).', A, L, c);
    end
    for i = 1:n
        if ~isCorner(i), continue, end
        f(Cf(i),:) = vertexFun(W(i,:).', A, L, c);
    end

    [Vd, E, F] = B.finish(rep);
    g = QuaPol(Vd, E, f, F);

    % An orientation or adjacency slip shows up as a face whose function is never reached, which
    % is silent in the values at most probe points and catastrophic at a few. The mesh is small,
    % so checking it is cheap and it is checked always rather than under a flag.
    assertEveryFaceReachable(g, size(f,1));
end

% ================================================================================================

function [base, dir, isRay, sdir] = edgeList(W, dFirst, dLast, unbounded)
% objective: the polygon's boundary elements in walk order.
% [output] base : index of the vertex each edge starts from
%          dir  : the edge's direction as WALKED, so the outward normal is a fixed rotation of it
%          isRay: true for the two unbounded elements
%          sdir : the direction the STRIP is parametrised along, from base towards the rest of the
%                 edge. Identical to `dir` except for the first ray, which is walked INWARD --
%                 keeping them as one array is how the sign of g0 gets flipped without anyone
%                 noticing, so they are two arrays.
    n = size(W,1);
    if ~unbounded
        base = (1:n).';
        dir  = W(mod((1:n), n) + 1, :) - W((1:n).', :);
        sdir = dir;
        isRay = false(n,1);
        return
    end
    fin  = W(2:n,:) - W(1:n-1,:);
    base  = [1; (1:n-1).'; n];
    dir   = [-dFirst(:).'; fin; dLast(:).'];
    sdir  = [ dFirst(:).'; fin; dLast(:).'];
    isRay = [true; false(n-1,1); true];
end

function j = nextVertex(i, n)
    j = mod(i, n) + 1;
end

function nOut = outwardNormal(e)
% objective: the OUTWARD unit normal of an edge walked in direction e with the interior on its
%            LEFT. Rotating by -90 degrees sends "left" to "in", so it sends e to "out".
    nOut = [e(2), -e(1)];
    nn = norm(nOut);
    if nn <= 0
        error('conjConvexPolygon:degenerateEdge', 'an edge has zero direction.');
    end
    nOut = nOut / nn;
end

function assertOrientation(W, dFirst, dLast, base, dir, nrm)
% objective: refuse a polygon that is not CCW, or whose ray directions do not point out of it,
%            rather than returning a mesh built on the wrong side.
%
% The test is the definition: an outward normal must not point towards any other vertex of the
% polygon. Checked against every vertex rather than against a centroid, because a centroid of an
% unbounded polygon's finite vertices can lie outside it.
    n = size(W,1);
    tol = 1e-9 * (1 + max(abs(W(:))));
    for k = 1:numel(base)
        p = W(base(k),:);
        for i = 1:n
            if (W(i,:) - p) * nrm(k,:).' > tol
                error('conjConvexPolygon:orientation', ...
                    ['edge %d''s outward normal points towards vertex %d, so the vertex list is ' ...
                     'not CCW with the interior on the left (or the polygon is not convex).'], k, i);
            end
        end
        for d = {dFirst, dLast}
            if ~isempty(d{1}) && d{1}(:).' * nrm(k,:).' > tol
                error('conjConvexPolygon:orientation', ...
                    ['a recession direction points out through edge %d''s outward normal, so the ' ...
                     'set described is not convex.'], k);
            end
        end
    end
    if n >= 3 && isempty(dFirst)
        a = 0;
        for i = 1:n
            j = nextVertex(i, n);
            a = a + W(i,1)*W(j,2) - W(j,1)*W(i,2);
        end
        if a <= 0
            error('conjConvexPolygon:orientation', ...
                'the bounded polygon has signed area %g; it must be positive (CCW).', a/2);
        end
    end
end

function d = dualScale(S, A, dFirst, dLast)
% objective: a length at which to place the ray MARKER points, so that the mesh's stored second
%            vertex of a ray is neither degenerate nor absurdly far from the rest of it.
%   The markers are only direction carriers (RatCon.m's `E`), but eval and plotting read them as
%   coordinates, so their magnitude should match the mesh's.
    d = 1;
    if size(S,1) >= 2
        d = 0.5 * max(pdistMax(S));
    end
    for dd = {dFirst, dLast}
        if ~isempty(dd{1}), d = max(d, 0.5 * norm((A*dd{1}(:)).')); end
    end
    if ~isfinite(d) || d <= 0, d = 1; end
end

function mx = pdistMax(S)
    n = size(S,1); mx = 0;
    for i = 1:n-1
        for j = i+1:n
            mx = max(mx, norm(S(i,:) - S(j,:)));
        end
    end
end

function B = meshBuilder(scale)
% objective: accumulate vertices, edges and per-edge face adjacency without the caller having to
%            keep three parallel index spaces straight. `finish` orients each edge's [left right]
%            from the faces the caller declared, which is the one thing that cannot be got right
%            by inspection.
    st.V = zeros(0,2);
    st.E = zeros(0,3);
    st.adj = zeros(0,2);       % [faceOnLeftOfWalk faceOnRightOfWalk] as DECLARED by the caller
    st.nf = 0;
    st.scale = scale;
    B = struct();
    B.state = st;
    B.vertex   = @(p)                addVertex(p);
    B.face     = @()                 addFace();
    B.segment  = @(i, j, fl, fr)     addEdge(i, j, 1, fl, fr);
    B.ray      = @(i, d, fl, fr)     addRay(i, d, fl, fr);
    B.finish   = @(rep)              finish(rep);
    B.count    = @()                 st.nf;

    function k = addVertex(p)
        st.V(end+1,:) = p; k = size(st.V,1);
    end
    function k = addFace()
        st.nf = st.nf + 1; k = st.nf;
    end
    function addEdge(i, j, isSeg, fl, fr)
        st.E(end+1,:) = [i, j, isSeg];
        st.adj(end+1,:) = [fl, fr];
    end
    function addRay(i, d, fl, fr)
        d = d(:).' / max(norm(d), eps);
        k = addVertex(st.V(i,:) + st.scale * d);
        addEdge(i, k, 0, fl, fr);
    end
    function [V, E, F] = finish(rep)
        % F(j,:) = [faceOnLEFT faceOnRIGHT] of the DIRECTED edge j, which is what RatCon's `F`
        % means. The caller declared the PAIR only; the side is decided here from each face's own
        % representative interior point, which is the one method that cannot be reasoned wrong.
        V = st.V; E = st.E;
        F = zeros(size(st.adj));
        for jj = 1:size(E,1)
            p = V(E(jj,1),:); q = V(E(jj,2),:);
            t = q - p;
            lft = [-t(2), t(1)];                     % +90 degrees: the left of the directed edge
            fa = st.adj(jj,1); fb = st.adj(jj,2);
            sa = (rep(fa,:) - p) * lft.';
            sb = (rep(fb,:) - p) * lft.';
            if abs(sa - sb) <= 1e-12 * (1 + abs(sa) + abs(sb))
                error('conjConvexPolygon:orientation', ...
                    ['edge %d cannot be oriented: the representative points of faces %d and %d ' ...
                     'lie on the same side of it, which means the declared adjacency is wrong.'], ...
                     jj, fa, fb);
            end
            if sa > sb, F(jj,:) = [fa fb]; else, F(jj,:) = [fb fa]; end
        end
    end
end

function [kp, kn] = incidentEdges(i, n, m, unbounded)
% objective: the two boundary elements meeting vertex i, as indices into the edge list. The
%   unbounded list is [ray0, e1..e(n-1), rayN] and is 1-based, so vertex i meets elements i and
%   i+1; the bounded list wraps.
    if unbounded
        kp = i; kn = i + 1;
    else
        kp = mod(i - 2, m) + 1; kn = i;
    end
end

function u = unitOr0(v)
    nv = norm(v);
    if nv <= 0, u = zeros(1, numel(v)); else, u = v / nv; end
end

function s6 = interiorFun(A, L, c)
% f*(s) = 1/2 (s-L)' A^{-1} (s-L) - c, in the stored weighted 6-basis [x^2 xy y^2 x y 1].
    M = inv(A); %#ok<MINV>
    g = -M*L;
    s6 = [M(1,1), M(1,2), M(2,2), g(1), g(2), 0.5*(L.'*M*L) - c];
end

function s6 = stripFun(v, e, A, L, c)
% f*(s) = <s,v> - q(v) + (<s,e> - g0)^2 / D  on the strip of the edge through v with direction e,
% where D = 2 e'Ae and g0 = <grad q(v), e>. This is the max along the LINE, which is what the
% cell's own definition makes correct there.
    D  = 2*(e.'*A*e);
    g0 = (A*v + L).' * e;
    qv = 0.5*(v.'*A*v) + L.'*v + c;
    e1 = e(1); e2 = e(2);
    px2 = e1^2/D;  pxy = 2*e1*e2/D;  py2 = e2^2/D;
    pxc = v(1) - 2*g0*e1/D;  pyc = v(2) - 2*g0*e2/D;  pc = -qv + g0^2/D;
    s6 = [2*px2, pxy, 2*py2, pxc, pyc, pc];
end

function s6 = vertexFun(v, A, L, c)
% f*(s) = <s,v> - q(v) on the vertex's normal cone: affine.
    qv = 0.5*(v.'*A*v) + L.'*v + c;
    s6 = [0, 0, 0, v(1), v(2), -qv];
end

function assertEveryFaceReachable(g, nf)
% objective: every face of the built mesh must be the one `eval` selects somewhere inside it.
%   A face that no point reaches is an adjacency or orientation slip, and it is silent -- the
%   values stay right at most probe points and are wrong on a region.
    miss = [];
    for i = 1:nf
        p = faceProbe(g, i);
        if isempty(p), miss(end+1) = i; continue, end %#ok<AGROW>
        [~, r] = g.eval(p);
        if ~(r == i || r == 0), miss(end+1) = i; end %#ok<AGROW>
    end
    if ~isempty(miss)
        error('conjConvexPolygon:unreachableFace', ...
            ['faces [%s] of the constructed conjugate are not reached by eval at their own ' ...
             'interior probe, which means the edge adjacency or the orientation is wrong.'], ...
            num2str(miss));
    end
end

function p = faceProbe(g, i)
% objective: a point in the interior of face i, built from the face's own boundary: the average of
%   its vertices, pushed along the average of its ray directions when it is unbounded.
    ej = find(any(g.F == i, 2));
    if isempty(ej), p = []; return, end
    vi = unique(g.E(ej,1:2));
    p = mean(g.V(vi,:), 1);
    rays = ej(g.E(ej,3) == 0);
    if ~isempty(rays)
        d = zeros(1,2);
        for t = 1:numel(rays)
            dd = g.V(g.E(rays(t),2),:) - g.V(g.E(rays(t),1),:);
            d = d + dd / max(norm(dd), eps);
        end
        if norm(d) > 1e-12
            p = p + 0.25 * max(1, max(abs(g.V(:)))) * d / norm(d);
        end
    end
end
