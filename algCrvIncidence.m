function [E, F, info] = algCrvIncidence(g, Vname, X)
% algCrvIncidence  The EDGE and FACE incidence of an AlgCrv: which vertices bound each edge, and
%   which faces each edge separates.
%
% objective: `algCrvVertices` names the vertices, which is what an object needs to have topology at
%   all; this is what a figure or an edge list needs on top of it. Without it the mesh is a set of
%   faces that happen to share curves, with nothing saying they are adjacent.
%
% [input]  g     : AlgCrv
%          Vname : as `algCrvVertices` returns it -- nv x 5 [edgeA tableA edgeB tableB rootIdx]
%          X     : optional, nv x 2, the vertices' coordinates -- `algCrvVertices`'s info.X. When
%                  omitted it is recomputed by calling `algCrvVertices` again (minutes).
% [output] E     : one row per ARC of the mesh, [edgeIndex, table, vA, vB, flag]. vA and vB index
%                  ROWS OF Vname and are 0 where the arc has no such endpoint; flag is 1 for a
%                  segment, 0 for a ray or line. Several rows may share edgeIndex: an edge is split
%                  into the arcs between consecutive vertices along it, an arc that bounds no face is
%                  left out, and one arc spelled by several rows of the table is listed once.
%          F     : one row per arc, the faces on each side: [k1 k2], with k2 = 0 where only one
%                  face has it as boundary.
%          info  : .edges (distinct edges), .arcs (rows of E), .segments, .rays, .oversubdivided
%                  (edges with more than two distinct vertex POINTS on them), .split (of those, split
%                  into arcs), .untraced (edges with a vertex whose tracing failed, left as one row
%                  read off the tables), .unowned (arcs that bound no face, dropped), .merged (arcs
%                  listed again under another spelling of the same curve, dropped), .shared,
%                  .overshared (arcs with more than two faces; the first two are kept)
%
% ==================================================================================================
% WHAT IS EXACT HERE AND WHAT IS COMBINATORIAL
%
% The VERTICES are exact names (`algCrvVertices`). Which vertices an arc runs between, and which
% faces it bounds, is decided in DOUBLES -- the same trade `algCrvVertices` makes for membership and
% `assembleQuaConCells` for a curved corner, for the same reason: an arc wrongly split changes the
% mesh's TOPOLOGY and never a face's VALUE, since `eval` reads the face functions and their sign
% conditions and neither mentions an arc. The exact alternative is comparing algebraic numbers of
% degree up to 16 along a quartic. What is decided numerically, and nothing else:
%
%   * THE BREAKPOINTS of an edge are every vertex of the mesh lying on its curve -- not only the
%     names that mention it, because a vertex named for one spelling of a curve lies on the others;
%   * THE ARCS: the curve is traced from each breakpoint along every branch leaving it (sign changes
%     of the curve on a small circle: two at a smooth point, more at a node) to the next breakpoint,
%     or off to infinity -- a ray. A doubled line l^2 = 0 is traced as l, since it never changes sign;
%   * THE FACES of an arc, by the DEFINITION of a boundary: just off the arc's midpoint on one side is
%     in the face and on the other side is not. So a row like x^2 >= 0, true on both sides of x = 0,
%     is not a boundary there. A thin face (a row of sign 0) bounds its own edge where its rows hold.
%
% An edge with no vertex on it at all has nothing to trace from; it stays one row read off the
% tables -- its endpoints 0, its faces those whose rows reference it -- which is what this routine
% returned for every edge before arcs were traced.
%
% A POINT NAMED TWICE -- three edges through one point give three names -- is ONE breakpoint. Names
% are clustered by coordinates, and an arc's endpoint is a name that mentions its edge when there is
% one.

    if nargin < 2 || nargin > 3
        error('algCrvIncidence:nargin', 'expects (g, Vname) or (g, Vname, X).');
    end
    if ~isa(g, 'AlgCrv')
        error('algCrvIncidence:type', 'expects an AlgCrv; got a %s.', class(g));
    end
    if ~isempty(Vname) && size(Vname,2) ~= 5
        error('algCrvIncidence:shape', ...
              'Vname must be nv x 5 [edgeA tableA edgeB tableB rootIdx]; got %d columns.', ...
              size(Vname,2));
    end
    if nargin < 3, X = []; end
    if ~isempty(X) && ~isequal(size(X), [size(Vname,1), 2])
        error('algCrvIncidence:shape', 'X must be nv x 2, row for row with Vname; got %s.', ...
              mat2str(size(X)));
    end

    % ---- every edge either table actually uses, in a stable order ---------------------------------
    used = zeros(0,2);
    for k = 1:g.nf
        rows = g.FC{k};
        for t = 1:size(rows,1)
            key = [rows(t,1), rowTable(rows, t)];
            if isempty(find(all(used == key, 2), 1)), used(end+1,:) = key; end %#ok<AGROW>
        end
    end

    ne = size(used,1);
    E = zeros(0, 5);
    F = zeros(0, 2);
    M = zeros(0, 2);                             % each row's traced midpoint, NaN where untraced
    Fall = {};                                   % each row's full face list, before truncation
    info = struct('edges', ne, 'arcs', 0, 'segments', 0, 'rays', 0, 'oversubdivided', 0, ...
                  'split', 0, 'untraced', 0, 'unowned', 0, 'merged', 0, 'shared', 0, ...
                  'overshared', 0);

    for e = 1:ne
        idx = used(e,1);  tbl = used(e,2);

        % ---- the vertices whose NAME mentions this edge -------------------------------------------
        onIt = [];
        for v = 1:size(Vname,1)
            if (Vname(v,1) == idx && Vname(v,2) == tbl) || (Vname(v,3) == idx && Vname(v,4) == tbl)
                onIt(end+1) = v; %#ok<AGROW>
            end
        end
        fs = facesUsing(g, idx, tbl);

        % ---- the breakpoints on this curve, clustered by POINT, then trace ----------------------------
        if ~isempty(Vname)
            if isempty(X)
                [V2, i2] = algCrvVertices(g);
                if ~isequal(V2, Vname)
                    error('algCrvIncidence:coords', ...
                          ['X was omitted and recomputing the vertices gave a different Vname, ' ...
                           'so the coordinates cannot be matched to the names. Pass X.']);
                end
                X = i2.X;
            end
            % THE BREAKPOINTS ARE EVERY MESH VERTEX ON THE CURVE, not only those named for it: the
            % diagonal x + y = 1 of the c = -2 object passes through (1/2, 1/2), a vertex named for
            % the pair of lines (x-1)^2 - y^2 that contains the diagonal, and not splitting there left
            % two faces' boundaries meeting at a point that one of them did not know was a corner.
            % Names that mention the edge come first, so a cluster is represented by one of them.
            C = curvePoly(g, idx, tbl);
            others = setdiff(1:size(Vname,1), onIt);
            others = others(arrayfun(@(v) onCurve(C, X(v,:)), others));
            cand = [onIt, others];
            [reps, P] = clusterPoints(X(cand,:));
            onIt = cand;
            if size(P,1) > 2, info.oversubdivided = info.oversubdivided + 1; end
            if ~isempty(cand)
                arcs = [];
                try
                    arcs = traceArcs(C, P);
                catch ME
                    if ~startsWith(ME.identifier, 'algCrvIncidence:'), rethrow(ME), end
                end
                if isempty(arcs)
                    info.untraced = info.untraced + 1;
                else
                    if size(P,1) > 2, info.split = info.split + 1; end
                    for a = 1:size(arcs,1)
                        own = [];
                        for k = fs
                            if arcBoundsFace(g, k, idx, tbl, arcs(a,3:4), arcs(a,5:6))
                                own(end+1) = k; %#ok<AGROW>
                            end
                        end
                        if isempty(own), info.unowned = info.unowned + 1; continue, end
                        vA = onIt(reps(arcs(a,1)));
                        vB = 0;
                        if arcs(a,2) > 0, vB = onIt(reps(arcs(a,2))); end
                        E(end+1,:) = [idx, tbl, vA, vB, double(vB > 0)]; %#ok<AGROW>
                        M(end+1,:) = arcs(a,3:4); %#ok<AGROW>
                        Fall{end+1} = own; %#ok<AGROW>
                    end
                    continue
                end
            end
            onIt = onIt(reps(1:min(2, numel(reps))));
        end

        % ---- nothing to trace from (or tracing failed): one row, read off the tables ----------------
        flag = 0;  vA = 0;  vB = 0;
        if numel(onIt) >= 1, vA = onIt(1); end
        if numel(onIt) >= 2, vB = onIt(2);  flag = 1; end
        E(end+1,:) = [idx, tbl, vA, vB, flag]; %#ok<AGROW>
        M(end+1,:) = [NaN NaN]; %#ok<AGROW>
        Fall{end+1} = fs; %#ok<AGROW>
    end

    % ---- ONE ARC OF THE MESH, SEVERAL ROWS OF THE TABLE -----------------------------------------------
    % The conic table holds the same zero set under different integer rows -- on the c = -2 object
    % y = 1 is there as y - 1, 1 - y, y^2 - 1 and (y - 1)^2 -- so one boundary arc can be traced once
    % per spelling. Two traced arcs with the same endpoint POINTS and the same midpoint are the same
    % arc: the first row is kept and the faces are united. A chord and a conic arc between the same two
    % points have different midpoints, so they are not merged.
    keep = true(size(E,1), 1);
    for r = 1:size(E,1)
        if ~keep(r) || any(isnan(M(r,:))), continue, end
        for s = r+1:size(E,1)
            if ~keep(s) || any(isnan(M(s,:))), continue, end
            if sameEnds(E(r,3:4), E(s,3:4), X) && norm(M(r,:) - M(s,:)) < 1e-6 * (1 + norm(M(r,:)))
                Fall{r} = unique([Fall{r}, Fall{s}]);
                keep(s) = false;
                info.merged = info.merged + 1;
            end
        end
    end
    E = E(keep,:);  Fall = Fall(keep);
    for r = 1:size(E,1)
        [F, info] = addFaces(F, info, E(r,5), Fall{r});
    end
end

function tf = sameEnds(a, b, X)
% Do the two endpoint pairs name the same two POINTS, in either order? 0 is the point at infinity.
    pa = endPts(a, X);  pb = endPts(b, X);
    tf = (samePt(pa(1,:), pb(1,:)) && samePt(pa(2,:), pb(2,:))) || ...
         (samePt(pa(1,:), pb(2,:)) && samePt(pa(2,:), pb(1,:)));
end

function P = endPts(v, X)
    P = inf(2,2);
    for c = 1:2, if v(c) > 0, P(c,:) = X(v(c),:); end, end
end

function tf = samePt(p, q)
    if any(isinf(p)) || any(isinf(q)), tf = all(isinf(p)) && all(isinf(q)); return, end
    tf = norm(p - q) <= 1e-7 * max(1, max(abs([p q])));
end

% ==================================================================================================
function [F, info] = addFaces(F, info, flag, fs)
    info.arcs = info.arcs + 1;
    if flag == 1, info.segments = info.segments + 1; else, info.rays = info.rays + 1; end
    if numel(fs) >= 2, info.shared = info.shared + 1; end
    if numel(fs) > 2, info.overshared = info.overshared + 1; end
    F(end+1,:) = [fs(1), 0];
    if numel(fs) >= 2, F(end,2) = fs(2); end
end

function tbl = rowTable(rows, t)
    tbl = 0;
    if size(rows,2) >= 4, tbl = rows(t,4); end
end

function fs = facesUsing(g, idx, tbl)
    fs = [];
    for k = 1:g.nf
        rows = g.FC{k};
        for t = 1:size(rows,1)
            if rows(t,1) == idx && rowTable(rows, t) == tbl, fs(end+1) = k; break, end %#ok<AGROW>
        end
    end
    fs = unique(fs);
end

function [reps, P] = clusterPoints(Xe)
% One entry per distinct POINT: reps(c) is the first row of Xe in cluster c, P(c,:) its coordinates.
    D = max(1, max(abs(Xe(:))));
    reps = [];  P = zeros(0,2);
    for i = 1:size(Xe,1)
        if isempty(P) || min(vecnorm(P - Xe(i,:), 2, 2)) > 1e-7 * D
            reps(end+1) = i; %#ok<AGROW>
            P(end+1,:) = Xe(i,:); %#ok<AGROW>
        end
    end
end

% --------------------------------------------------------------------------------------------------
function C = curvePoly(g, idx, tbl)
% The edge as a polyN in (x, y), C(i,j) the coefficient of x^(i-1) y^(j-1) -- `algCrvVertices`'s
% curveTable for one entry.
    if tbl == 1
        a = g.EcA(idx);
        C = bridgeQuadCut(a.FP, a.num, a.den);
    else
        c = g.EcQ(idx,:);
        % A DOUBLED LINE, l^2 = 0, has the zero set of l but never changes sign, so there is nothing
        % to trace; the c = -2 object carries x^2 and y^2 as rows. Its matrix has rank 1 -- every 2x2
        % minor vanishes, exactly, in integers -- and any nonzero row of it is a multiple of l.
        M2 = [2*c(1), c(2), c(4); c(2), 2*c(3), c(5); c(4), c(5), 2*c(6)];
        if rankOne(M2)
            [~, i] = max(abs(diag(M2)));
            l = M2(i,:);
            C = zeros(2,2);  C(2,1) = l(1);  C(1,2) = l(2);  C(1,1) = l(3);
            return
        end
        C = zeros(3,3);
        C(3,1) = c(1);  C(2,2) = c(2);  C(1,3) = c(3);
        C(2,1) = c(4);  C(1,2) = c(5);  C(1,1) = c(6);
    end
end

function tf = rankOne(M)
    tf = any(M(:) ~= 0);
    for i = 1:2
        for j = i+1:3
            for k = 1:2
                for l = k+1:3
                    if M(i,k)*M(j,l) - M(i,l)*M(j,k) ~= 0, tf = false; return, end
                end
            end
        end
    end
end

function tf = onCurve(C, x)
% Does x lie on C = 0, to the accuracy a vertex is located with? Relative to the size of the terms.
    n = size(C,1);  m = size(C,2);
    T = abs(C) .* (abs(x(1)).^(0:n-1)).' .* abs(x(2)).^(0:m-1);
    tf = abs(polyAndGrad(C, x)) <= 1e-7 * max(1, sum(T(:)));
end

function [v, gr] = polyAndGrad(C, x)
    n = size(C,1);  m = size(C,2);
    px = x(1).^(0:n-1);  py = x(2).^(0:m-1);
    dpx = [0, (1:n-1) .* x(1).^(0:n-2)];
    dpy = [0, (1:m-1) .* x(2).^(0:m-2)];
    v  = px * C * py.';
    gr = [dpx * C * py.', px * C * dpy.'];
end

function x = project(C, x)
% Newton onto C = 0 along the gradient, a few steps.
    for it = 1:8
        [v, gr] = polyAndGrad(C, x);
        n2 = gr * gr.';
        if n2 == 0, return, end
        dx = -v / n2 * gr;
        x = x + dx;
        if norm(dx) < 1e-14 * max(1, norm(x)), return, end
    end
end

% --------------------------------------------------------------------------------------------------
function arcs = traceArcs(C, P)
% Every arc of the curve C = 0 that leaves one of the points P, as rows [i, j, mx, my]: from P(i) to
% the NEXT point P(j) along the curve (j = 0 when the curve runs off to infinity first), with (mx, my)
% a point on the arc half-way along it. Each arc between two points is found from both ends and kept
% once. The branches leaving a point are found by sign changes of C on a small circle around it,
% which gives two at a smooth point and more at a node.
    np = size(P,1);
    c0 = mean(P, 1);
    D = max(1, max(vecnorm(P - c0, 2, 2)) * 2);
    hmin = 1e-5 * D;
    far = 50 * D;
    arcs = zeros(0,6);
    for i = 1:np
        others = setdiff(1:np, i);
        dnear = min([inf; vecnorm(P(others,:) - P(i,:), 2, 2)]);
        r0 = min(1e-3 * D, 0.3 * dnear);
        for x1 = branchStarts(C, P(i,:), r0)'
            [j, mid, tan] = walk(C, P, i, x1', hmin, far, c0);
            if isnan(j), continue, end
            dup = false;
            for a = 1:size(arcs,1)
                if ((arcs(a,1) == i && arcs(a,2) == j) || (arcs(a,1) == j && arcs(a,2) == i)) ...
                        && norm(arcs(a,3:4) - mid) < 1e-2 * D
                    dup = true; break
                end
            end
            if ~dup, arcs(end+1,:) = [i, j, mid, tan]; end %#ok<AGROW>
        end
    end
end

function S = branchStarts(C, p, r)
% Points on the circle |x - p| = r where C changes sign, each projected onto the curve.
    th = linspace(0, 2*pi, 721);  th(end) = [];
    Q = p + r * [cos(th(:)), sin(th(:))];
    v = zeros(numel(th),1);
    for s = 1:numel(th), v(s) = polyAndGrad(C, Q(s,:)); end
    S = zeros(0,2);
    for s = 1:numel(th)
        t = mod(s, numel(th)) + 1;
        if sign(v(s)) ~= sign(v(t)) || v(s) == 0
            lo = th(s);  hi = th(t);  if t == 1, hi = hi + 2*pi; end
            vlo = v(s);
            for it = 1:50                                   % bisect on the circle
                mdl = (lo + hi) / 2;
                vm = polyAndGrad(C, p + r * [cos(mdl), sin(mdl)]);
                if sign(vm) == sign(vlo), lo = mdl; vlo = vm; else, hi = mdl; end
            end
            q = p + r * [cos(lo), sin(lo)];
            if isempty(S) || min(vecnorm(S - q, 2, 2)) > 1e-3 * r, S(end+1,:) = q; end %#ok<AGROW>
        end
    end
end

function [j, mid, tan] = walk(C, P, i, x, hmin, far, c0)
% Follow the curve from x (just off P(i), heading away from it) until it reaches a point of P.
% Step length adapts to the distance to the nearest point and to the local turning, so the walk is
% cheap far from every vertex and careful near one.
    np = size(P,1);
    [~, gr] = polyAndGrad(C, x);
    t = [-gr(2), gr(1)];
    tan = [NaN NaN];
    if norm(t) == 0, j = NaN; mid = [NaN NaN]; return, end
    t = t / norm(t);
    if dot(t, x - P(i,:)) < 0, t = -t; end
    path = [P(i,:); x];
    len = norm(x - P(i,:));
    j = NaN;
    for step = 1:200000
        d = vecnorm(P - x, 2, 2);
        if len < 10 * hmin, d(i) = inf; end               % not back at the start yet
        h = min(max(0.3 * min(d), hmin), 0.05 * max(far/50, norm(x - c0)));
        while true
            y = project(C, x + h * t);
            [~, gy] = polyAndGrad(C, y);
            ty = [-gy(2), gy(1)];
            if norm(ty) == 0, ty = t; else, ty = ty / norm(ty); end
            if dot(ty, t) < 0, ty = -ty; end
            if (dot(ty, t) > cos(0.2) && abs(norm(y - x) - h) < 0.5 * h) || h <= hmin
                break
            end
            h = h / 2;
        end
        % did the step pass a point?
        seg = y - x;  L2 = seg * seg.';
        for q = 1:np
            if q == i && len < 10 * hmin, continue, end
            s = 0;
            if L2 > 0, s = min(1, max(0, ((P(q,:) - x) * seg.') / L2)); end
            if norm(x + s * seg - P(q,:)) < max(hmin, 1e-9)
                j = q;
                path(end+1,:) = P(q,:); %#ok<AGROW>
                break
            end
        end
        if ~isnan(j), break, end
        path(end+1,:) = y; %#ok<AGROW>
        len = len + norm(y - x);
        x = y;  t = ty;
        if norm(x - c0) > far, j = 0; break, end
    end
    if isnan(j), mid = [NaN NaN]; return, end
    % the point half-way along, by arc length -- on a ray, half-way to where the walk stopped, which
    % is past every vertex and so still inside the one arc the ray is
    if j == 0
        k = find(vecnorm(path - P(i,:), 2, 2) > 2 * max(vecnorm(P - P(i,:), 2, 2)) + hmin, 1);
        if isempty(k), k = size(path,1); end
        cut = path(1:k,:);
    else
        cut = path;
    end
    sl = [0; cumsum(vecnorm(diff(cut), 2, 2))];
    % interpolated inside the step, not snapped to a step end: steps grow to 5% of the scale far from
    % a vertex, and the SAME arc walked from its other end must land on the same point to be merged
    k = max(2, find(sl >= sl(end) / 2, 1));
    s = (sl(end)/2 - sl(k-1)) / max(eps, sl(k) - sl(k-1));
    mid = project(C, cut(k-1,:) + s * (cut(k,:) - cut(k-1,:)));
    % the tangent from the PATH, not the gradient: a doubled line such as y^2 = 0 has a zero
    % gradient all along itself, and it is exactly such rows the side test must see through
    tan = cut(k,:) - cut(k-1,:);
    tan = tan / max(eps, norm(tan));
end

% --------------------------------------------------------------------------------------------------
function tf = arcBoundsFace(g, k, idx, tbl, x, tan)
% Does face k have edge (idx, tbl) as boundary AT x, a point of the arc with tangent tan? By the
% DEFINITION of a boundary: just off the arc on one side is in face k and just off it on the other
% is not. A row that holds on BOTH sides -- x^2 >= 0 along x = 0, the doubled lines the c = -2 object
% is full of -- is therefore not a boundary there, which reading the row as an equality on the curve
% could not see. A THIN face (a row of sign 0 on this edge) is outside on both sides and is its own
% case: it bounds its own edge wherever its other rows hold.
    rows = g.FC{k};
    thin = any(rows(:,1) == idx & arrayfun(@(t) rowTable(rows, t), (1:size(rows,1)).') == tbl ...
               & rows(:,2) == 0);
    if thin
        tf = inFace(g, k, x, 1e-7, [idx, tbl]);
        return
    end
    n = [-tan(2), tan(1)];
    ep = 1e-6 * max(1, norm(x));
    tf = xor(inFace(g, k, x + ep * n, 0, []), inFace(g, k, x - ep * n, 0, []));
end

function tf = inFace(g, k, x, tol, skip)
% Face k's rows at x, each side * value >= -tol; an equality row as |value| <= tol. The row of edge
% `skip` ([idx tbl]) is not tested.
    tf = true;
    rows = g.FC{k};
    for t = 1:size(rows,1)
        ti = rows(t,1);  tt = rowTable(rows, t);
        if ~isempty(skip) && ti == skip(1) && tt == skip(2), continue, end
        if tt == 1
            a = g.EcA(ti);
            v = algBranch(a.FP, x, a.FB);
            if ~isfinite(v), tf = false; return, end
            r = ratQ.evalFace(a.num, a.den, x) - v;
        else
            c = g.EcQ(ti, :);
            br = 0;
            if size(rows,2) >= 3, br = rows(t,3); end
            if br ~= 0
                r = ratQ.pairLineEval(c, br, x);
            else
                r = c(1)*x(1)^2 + c(2)*x(1)*x(2) + c(3)*x(2)^2 + c(4)*x(1) + c(5)*x(2) + c(6);
            end
        end
        if rows(t,2) == 0
            if abs(r) > max(tol, 1e-9), tf = false; return, end
        elseif rows(t,2) * r < -tol
            tf = false; return
        end
    end
end
