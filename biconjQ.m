function h = biconjQ(obj)
% biconjQ  The biconjugate f** = cl co f, computed EXACTLY over the rationals.
%
% objective: the closed convex envelope of a QuaPol, with every coefficient an exact rational and
%   every decision made by exact integer arithmetic.
%
% [input]  obj : QuaPol, operable (degree <= 2) and EXACT (it must carry fN/fD -- see
%                QuaPol.assertExact for why computing exactly on inexact data is worse than not
%                computing at all).
% [output] h   : QuaCon.
%
% ------------------------------------------------------------------------------------------------
% DIRECT, NOT conj(conj(f)). ALGORITHM.md states the principle and the measurement behind it:
% CCA2 has a direct algorithm for the envelope, and computing it as a double conjugation is a
% detour that cost 436 s on an input whose answer was "return f unchanged". The double conjugation
% would also need the conjugate of a QuaCon, whose return type is an open question (the plan's
% Phase 6), so the direct route is both faster and better defined.
%
% WHAT IS IMPLEMENTED
%
%   * f CONVEX on its own domain -> f itself. Nothing to compute: co f = f. Decided exactly by
%     testing every piece's Hessian for positive semidefiniteness, which is Sylvester on integers.
%
%   * a CONCAVE or AFFINE piece on a bounded polygon -> the affine interpolant over the LOWER HULL
%     of the lifted vertices. This is the case the convex envelope is cleanest for, and the reason
%     is that a concave function on a polytope has its whole envelope determined by the VERTEX
%     VALUES: the graph of co f is the lower convex hull of {(v, q(v))}, so the computation is a
%     hull of m points in R^3 with rational coordinates -- no calculus, no calls to an optimiser,
%     and every predicate an integer sign test.
%
% WHY THE ANSWER IS A QuaCon AND NOT YET THE PLANNED AlgAlg. Both implemented cases have RATIONAL
% face coefficients: the convex case returns f, and the concave case returns affine functions
% interpolating rational vertex values. AlgAlg exists for the faces that CANNOT be written
% rationally -- an affine cell of a general f** is <p,x> - f*(p) with p a dual vertex of degree up
% to 4, so it is rational only over Q(p) and the cell has to NAME the dual vertex instead. Nothing
% here produces one, so introducing the type now would be speculative; the trigger is the
% INDEFINITE piece, which is refused by name below and is where the first such face appears.

    if ~isa(obj, 'QuaPol')
        error('PLQ:biconjQ:input', ...
            'biconjQ takes a QuaPol (quadratic on polyhedral); got a %s.', class(obj));
    end
    obj.assertOperable();
    obj.assertExact();

    nf = size(obj.fN, 1);

    % ---- a domain of dimension < 2, BEFORE the convex short-circuit -------------------------------------------------------------
    % co f is convex, so its domain is conv(dom f). For a needle or a single segment that is the
    % domain itself -- a point or a segment -- and both are THIN, which the H-form now expresses
    % with an EQUALITY side.
    %
    % BEFORE the convex short-circuit, not after, and the ordering is load-bearing: an affine or
    % PSD q on such a domain IS convex, so that branch would claim it and then read a face through
    % pieceGeometryB -- and these meshes have nf = 0 and an empty F, so there is no face to read.
    % Measured: the needle and the convex segment both raised noFace while the CONCAVE segment,
    % which falls past the short-circuit, worked.
    if obj.dom.dim < 2
        if size(obj.fN,1) > 1
            error('PLQ:biconjQ:ambiguousChain', ...
                ['a domain of dimension < 2 carrying %d different functions has no F to say ' ...
                 'which edge takes which.'], size(obj.fN,1));
        end
        h = thinEnvelope(obj);
        return
    end

    % ---- f CONVEX: there is nothing to compute -------------------------------------------------
    % The necessary-and-here-sufficient test for a SINGLE piece is that its Hessian is PSD. For a
    % multi-piece f that is necessary but not sufficient (the gradient jump across every shared
    % edge has to be consistent too -- ALGORITHM.md), so a multi-piece input is only accepted on
    % this branch when the caller has asserted convexity through fIsConvex, which is exactly what
    % that flag is for.
    allPSD = true;
    for k = 1:nf
        Hn = [obj.fN(k,5) obj.fN(k,6); obj.fN(k,6) obj.fN(k,7)];
        if ~ratQ.isPSD2(Hn), allPSD = false; break, end
    end
    if allPSD && (nf == 1 || isequal(obj.fIsConvex, true))
        h = quaPolAsQuaCon(obj);
        return
    end
    if isequal(obj.fIsConvex, true) && ~allPSD
        % The free necessary condition contradicts the caller's assertion. Refused LOUDLY rather
        % than trusted -- ALGORITHM.md's rule, because the failure it prevents is silent: biconj
        % would return a non-convex f as its own convex envelope.
        error('PLQ:biconjQ:notConvex', ...
            ['fIsConvex is asserted true but some piece has a Hessian that is not positive ' ...
             'semidefinite, which is a free necessary condition for convexity.']);
    end

    if nf ~= 1
        error('PLQ:biconjQ:notImplemented', ...
            ['the convex envelope of a %d-piece non-convex f COUPLES its pieces -- the convex ' ...
             'hull of a union is not determined piecewise -- so it is not a fold like the ' ...
             'conjugate is. Only a single piece is implemented.'], nf);
    end

    Hn = [obj.fN(1,5) obj.fN(1,6); obj.fN(1,6) obj.fN(1,7)];
    D  = ratQ.detExact(Hn);
    [~, ~, bounded] = pieceGeometryB(obj, 1);
    if ~bounded
        error('PLQ:biconjQ:unbounded', ...
            ['the envelope of a non-convex piece on an UNBOUNDED domain can be -infinity (the ' ...
             'function may be unbounded below on its own recession cone), which is a correct ' ...
             'answer with nowhere to be stored -- the same gap conjCPLQ records as ' ...
             'convEnvUnbounded:minusInfinity.']);
    end
    % ---- EDGE-CONCAVITY, which is weaker than concavity and is the real condition -------------
    % The envelope is the lower hull of the lifted vertices whenever q is concave along every EDGE
    % DIRECTION of the polygon -- H itself need not be negative semidefinite. Why the weaker
    % condition suffices, on a triangle and hence on each cell of the hull: q - L vanishes at the
    % three corners and is concave along each edge, so q >= L on the BOUNDARY; and with H not
    % positive semidefinite q - L has no interior minimum, so its minimum over the closed triangle
    % is on that boundary and q >= L throughout. (This is the classical vertex-polyhedral /
    % edge-concave criterion; it is verified here rather than taken on trust -- see
    % biconjQTest and .claude/envelope-sweep.m, which check the envelope's three defining
    % properties on indefinite edge-concave inputs.)
    %
    % IT MATTERS BECAUSE IT IS WHERE McCORMICK LIVES. q = xy on the unit square has H = [0 1; 1 0],
    % genuinely indefinite, and yet d'Hd = 0 along both edge directions -- so its envelope is the
    % lower hull of the four corner values, which is exactly max(0, x+y-1), the McCormick /
    % Al-Khayyal-Falk envelope. ALGORITHM.md lists that as a case worth having in closed form.
    edgeCurv = [];
    for j = find(any(obj.F == 1, 2)).'
        d = obj.VN(obj.E(j,2),:)/obj.VD(obj.E(j,2)) - obj.VN(obj.E(j,1),:)/obj.VD(obj.E(j,1));
        [dn, ~] = ratQ.fromDouble(d);
        edgeCurv(end+1) = ratQ.chk(dn * Hn * dn.', 'edge curvature'); %#ok<AGROW>
    end
    if any(edgeCurv > 0)
        error('PLQ:biconjQ:notImplemented', ...
            ['the envelope is implemented for a CONVEX f (nothing to compute) and for a piece ' ...
             'that is concave along every EDGE DIRECTION (the lower hull of the lifted ' ...
             'vertices). This piece curves UPWARD along %d of its %d edges, so its envelope is ' ...
             'not determined by the vertex values -- it needs [COAP] A.2-A.5, and is where the ' ...
             'first face that cannot be written rationally appears (see this file''s header on ' ...
             'AlgAlg).'], sum(edgeCurv > 0), numel(edgeCurv));
    end

    h = concaveEnvelope(obj);
end

% ==================================================================================================

function h = concaveEnvelope(obj)
% objective: co(q + I_P) for q CONCAVE or AFFINE on a bounded polygon P, exactly.
%
% THE MATHEMATICS IS ONE SENTENCE. A concave function on a polytope lies above every chord, so the
% lower convex hull of its graph is the lower convex hull of the finitely many LIFTED VERTICES
% {(v_i, q(v_i))} -- and that hull is the graph of co f. A convex piecewise-affine function is the
% MAX of the affine functions of its facets, so
%       co f (x) = max over lower-hull facets of  l_F(x),      x in P,
% and +infinity off P. Both halves are exactly the structure conjQ already builds: a max of affine
% functions, intersected with the polygon's own half-planes.
%
% EVERYTHING IS INTEGER. The lifted points are put over one denominator, a candidate facet is the
% plane through three of them (an integer cross product), and "is this a LOWER hull facet" is the
% sign of an integer for every other point. No hull library, no orientation predicate that could be
% read backwards, and no tolerance.
    [fN, fD] = obj.faceQ(1);
    Hn = [fN(5) fN(6); fN(6) fN(7)];
    Ln = [fN(8); fN(9)];
    kn = fN(10);

    [Vi, vd] = ratQ.combineDen(obj.VN, obj.VD);
    [vs, ~, ~] = pieceGeometryB(obj, 1);
    m = numel(vs);
    if m < 3
        error('PLQ:biconjQ:degenerate', 'a 2-D piece needs at least three vertices; got %d.', m);
    end

    % ---- the lifted vertices, over ONE denominator ---------------------------------------------
    % x and y are Vi/vd and the value is qv/(2*fD*vd^2), so the common denominator is 2*fD*vd^2 and
    % the x,y numerators scale up by 2*fD*vd.
    Dl = ratQ.chk(2 * fD * vd^2, 'lift denominator');
    Q3 = zeros(m, 3);
    for i = 1:m
        v = Vi(vs(i), :).';
        qv = ratQ.chk(v.'*Hn*v + 2*vd*(Ln.'*v) + 2*vd^2*kn, 'vertex value');
        Q3(i,:) = [ratQ.chk(2*fD*vd*v(1), 'lift x'), ratQ.chk(2*fD*vd*v(2), 'lift y'), qv];
    end

    % ---- the LOWER hull facets, by brute force over triples --------------------------------------
    % m is the number of vertices of one polygon -- three to a handful -- so O(m^4) is a few hundred
    % integer operations and needs no incremental hull. Every test below is exact.
    facets = zeros(0, 4);                       % rows [a b c d] with a X + b Y + c Z = d, c > 0
    for i = 1:m
        for j = i+1:m
            for k = j+1:m
                u = Q3(j,:) - Q3(i,:);
                w = Q3(k,:) - Q3(i,:);
                n = ratQ.chk(cross(u, w), 'facet normal');
                if n(3) == 0, continue, end     % vertical plane: not the graph of a function
                if n(3) < 0, n = -n; end        % orient upward, so "above" is n.p >= d
                d = ratQ.chk(n * Q3(i,:).', 'facet offset');
                below = false;
                for p = 1:m
                    if ratQ.chk(n * Q3(p,:).' - d, 'facet side') < 0, below = true; break, end
                end
                if below, continue, end          % some vertex lies UNDER it: not a lower facet
                % Reduce to a primitive integer row so that one plane found through several
                % different triples is ONE row and `unique` collapses it -- the same canonical-form
                % argument ratQ.conic makes for curves.
                g = gcd(gcd(abs(n(1)), abs(n(2))), gcd(abs(n(3)), abs(d)));
                if g == 0, g = 1; end
                facets(end+1, :) = [n, d] / g; %#ok<AGROW>
            end
        end
    end
    facets = unique(facets, 'rows');
    if isempty(facets)
        error('PLQ:biconjQ:noFacet', ...
            'no lower-hull facet was found, which cannot happen for three non-collinear vertices.');
    end

    % ---- each facet's affine function, and the cell where it is the largest ---------------------
    % From a X + b Y + c Z = d in the lifted coordinates, the function is
    %       l(x,y) = (d - a*Dl*x - b*Dl*y) / (Dl*c)  ... written over the denominator Dl*c.
    nfac = size(facets,1);
    % THE ENVELOPE'S DOMAIN IS conv(P), NOT P, and for a non-convex face those differ. co f must
    % be CONVEX, so its domain is convex: the closed convex hull of dom f. For a convex face the
    % two coincide and nothing changes; for a non-convex one, using the face's own edges produced a
    % region that was not even the face -- measured on an L, where biconjQ returned +Inf at (2,0),
    % (2,1), (1,2) and (1.5,0.5), all of which are IN the L and where the envelope is finite.
    %
    % The FACE FUNCTIONS needed no change: the lower hull of the lifted vertices is already the
    % right answer there, reflex vertices included. Checked against an LP over a dense sample of
    % the true L -- the hull and the LP agree at every probe, inside the L and in conv(L) beyond
    % it, so only the domain was wrong.
    dom = hullHalfPlanes(Vi, vd, vs);
    cells = struct('num', {}, 'den', {}, 'con', {});
    for t = 1:nfac
        a = facets(t,1);  b = facets(t,2);  c = facets(t,3);  d = facets(t,4);
        num = [0 0 0 0, 0, 0, 0, ratQ.chk(-a*Dl,'c8'), ratQ.chk(-b*Dl,'c9'), d];
        den = ratQ.chk(Dl * c, 'facet denominator');

        rows = dom;
        for s = 1:nfac
            if s == t, continue, end
            [dn, ~] = ratQ.sub(num, den, ...
                [0 0 0 0, 0, 0, 0, ratQ.chk(-facets(s,1)*Dl,'c8'), ...
                 ratQ.chk(-facets(s,2)*Dl,'c9'), facets(s,4)], ...
                ratQ.chk(Dl * facets(s,3), 'facet denominator'));
            if all(dn == 0), continue, end
            rows(end+1,:) = [ratQ.chk(2*dn(8),'d1'), ratQ.chk(2*dn(9),'d2'), ...
                             ratQ.chk(2*dn(10),'d0')]; %#ok<AGROW>
        end
        if ~ratQ.feasible2(rows, true), continue, end   % this facet never wins inside P
        con = zeros(size(rows,1), 7);
        for r = 1:size(rows,1)
            con(r,:) = [ratQ.conic([0 0 0, rows(r,:)]), sgnOfB(rows(r,:), +1)];
        end
        cells(end+1) = struct('num', num, 'den', den, 'con', con); %#ok<AGROW>
    end

    h = assembleQuaConCells(cells);
end

function rows = hullHalfPlanes(Vi, vd, vs)
% objective: the half-planes of the CONVEX HULL of the given vertices, exactly.
% [output] rows : h x 3, [a b c] meaning a*x + b*y + c >= 0 on the hull, in ACTUAL coordinates
%
% A pair of vertices spans a hull edge exactly when every other vertex lies weakly on one side of
% the line through them -- an integer sign test per vertex, so the hull is DECIDED. O(m^3) for a
% face's handful of corners, which needs no incremental hull algorithm and has no orientation
% convention to read backwards.
%
% Numerator coordinates decide the side; the emitted row is in actual coordinates (linear part
% scaled by vd). See conjQ's pieceShape for why mixing the two yields a parallel line.
    rows = zeros(0,3);
    m = numel(vs);
    for a = 1:m
        for b = 1:m
            if a == b, continue, end
            p = Vi(vs(a), :);  q = Vi(vs(b), :);
            n  = [q(2)-p(2), -(q(1)-p(1))];
            if all(n == 0), continue, end                 % coincident vertices
            c0 = ratQ.chk(-(n * p.'), 'offset');
            sgn = 0;  ok = true;
            for i = 1:m
                t = ratQ.chk(n(1)*Vi(vs(i),1) + n(2)*Vi(vs(i),2) + c0, 'hull side');
                if t == 0, continue, end
                if sgn == 0, sgn = sign(t);
                elseif sign(t) ~= sgn, ok = false; break
                end
            end
            if ~ok || sgn == 0, continue, end             % vertices on both sides: not a hull edge
            row = sgn * [ratQ.chk(vd*n(1), 'hull normal'), ...
                         ratQ.chk(vd*n(2), 'hull normal'), c0];
            row = row / gcdRow(row);
            if ~any(all(rows == row, 2)), rows(end+1,:) = row; end %#ok<AGROW>
        end
    end
end

function g = gcdRow(row)
    g = 0;
    for i = 1:numel(row)
        if row(i) ~= 0, g = gcd(g, abs(row(i))); end
    end
    if g == 0, g = 1; end
end

function rows = polygonHalfPlanes(obj, Vi, vd, vs, k)
% objective: face k's own constraints, as exact half-planes [a b c] meaning a*x + b*y + c >= 0 at
%            the ACTUAL coordinates (x,y).
%
% TWO COORDINATE SYSTEMS MEET HERE AND MUST NOT BE MIXED. `Vi` holds vertex NUMERATORS over the
% shared denominator vd, so a point's actual coordinates are X/vd. With n a normal to the edge and
% c0 = -n.p taken at numerators, one and the same line reads
%
%       n.X + c0 = 0        in numerator coordinates
%       vd*n.x + c0 = 0     in actual coordinates,   since X = vd*x
%
% so the SIDE test (done at numerators, where the vertices live) uses n, while the emitted ROW
% (read by eval at actual points) uses vd*n. Using one where the other belongs rescales the offset
% and yields a PARALLEL line: measured on a triangle with half-integer vertices, whose domain came
% out as 2x + 5y >= 3 where the edge is 2x + 5y = 3/2, so the function excluded two of its own
% three vertices. Invisible on integer vertices, which is what every hand-built fixture had.
%
% The side is fixed by requiring the piece's OWN vertices to satisfy it, so no orientation
% convention is consulted and none can be read backwards.
    rows = zeros(0,3);
    for j = find(any(obj.F == k, 2)).'
        p = Vi(obj.E(j,1), :);  q = Vi(obj.E(j,2), :);
        n  = [q(2)-p(2), -(q(1)-p(1))];
        c0 = ratQ.chk(-(n * p.'), 'offset');
        sgn = 0;
        for i = 1:numel(vs)
            v = Vi(vs(i), :);
            t = ratQ.chk(n(1)*v(1) + n(2)*v(2) + c0, 'side');
            if t ~= 0, sgn = sign(t); break, end
        end
        if sgn == 0, continue, end                  % every vertex on the line: degenerate
        rows(end+1,:) = sgn * [ratQ.chk(vd*n(1), 'edge normal'), ...
                               ratQ.chk(vd*n(2), 'edge normal'), c0]; %#ok<AGROW>
    end
end

function s = sgnOfB(row, want)
% The H-form side, after ratQ.conic may have negated the row -- the same guard conjQ's sgnOf is,
% duplicated because MATLAB gives a function file no way to share a local with another file.
    r = [0 0 0, row(1), row(2), row(3)];
    nz = find(r ~= 0, 1);
    if isempty(nz)
        error('ratQ:zeroConic', 'a constraint with all-zero coefficients names no half-plane.');
    end
    if r(nz) < 0, s = -want; else, s = want; end
end

function h = quaPolAsQuaCon(obj)
% objective: a convex QuaPol re-expressed as a QuaCon, its faces unchanged.
% co f = f when f is convex, so this is the whole computation for that branch.
    [Vi, vd] = ratQ.combineDen(obj.VN, obj.VD);
    cells = struct('num', {}, 'den', {}, 'con', {});
    for k = 1:size(obj.fN,1)
        [vs, ~, bounded] = pieceGeometryB(obj, k);
        if ~bounded
            error('PLQ:biconjQ:unbounded', ...
                'piece %d is unbounded; re-expressing it in H-form needs its recession cone.', k);
        end
        rows = polygonHalfPlanes(obj, Vi, vd, vs, k);
        % A convex f is returned UNCHANGED, so its pieces must be re-expressed exactly -- and a
        % face described by half-planes IS convex. A non-convex face cannot be, so re-expressing it
        % would silently replace the region by something else. Detected by comparing against the
        % face's own convex hull: for a convex face the two descriptions have the same rows.
        hull = hullHalfPlanes(Vi, vd, vs);
        if size(hull,1) ~= size(rows,1)
            error('PLQ:biconjQ:nonConvexPiece', ...
                ['piece %d is NOT CONVEX, and co f = f here, so the answer must carry that same ' ...
                 'non-convex region -- which an intersection of half-planes cannot express. ' ...
                 'Unlike the conjugate this cannot be fixed by splitting, because the convex ' ...
                 'hull of a union is not the union of hulls.'], k);
        end
        con = zeros(size(rows,1), 7);
        for r = 1:size(rows,1)
            con(r,:) = [ratQ.conic([0 0 0, rows(r,:)]), sgnOfB(rows(r,:), +1)];
        end
        cells(end+1) = struct('num', obj.fN(k,:), 'den', obj.fD(k), 'con', con); %#ok<AGROW>
    end
    h = assembleQuaConCells(cells);
end

function [vs, rays, bounded] = pieceGeometryB(obj, k)
% face k's true vertices and recession directions -- see conjQ's pieceGeometry, which this mirrors.
    own = find(any(obj.F == k, 2));
    if isempty(own)
        error('PLQ:biconjQ:noFace', 'face %d has no edges.', k);
    end
    E = obj.E(own, :);
    [Vi, ~] = ratQ.combineDen(obj.VN, obj.VD);
    vs = [];  rays = zeros(0,2);
    for j = 1:size(E,1)
        if E(j,3) ~= 0
            vs = [vs; E(j,1); E(j,2)]; %#ok<AGROW>
        else
            vs = [vs; E(j,1)]; %#ok<AGROW>
            rays(end+1,:) = ratQ.chk(Vi(E(j,2),:) - Vi(E(j,1),:), 'ray direction'); %#ok<AGROW>
        end
    end
    vs = unique(vs);
    bounded = isempty(rays);
end

function h = thinEnvelope(obj)
% objective: co f when the DOMAIN itself has dimension < 2 -- a needle or a single segment.
%
% co f is convex, so its domain is conv(dom f), and for these inputs that is the domain itself: a
% point or a segment. Both are THIN, which the H-form expresses with an EQUALITY side (a side of 0
% on a curve means ON it) -- the same machinery conj uses for a point- or line-supported conjugate.
%
%   NEEDLE. dom f is one point p, and a single point is convex, so co f = f: the value q(p) at p
%   and +infinity elsewhere. Two equalities and a constant face.
%
%   SEGMENT. The problem is ONE-dimensional in the segment parameter. Where q curves UP along the
%   segment (d'Hd >= 0) it is already convex there and co f = f; where it curves DOWN, co f is the
%   CHORD -- the affine interpolant of the two endpoint values, which is the 1-D lower hull. One
%   equality for the line, two inequalities for the ends.
%
% A CHAIN of edges is not handled here and is not thin: the convex hull of two non-collinear
% segments is a two-dimensional polygon, so that case belongs with the ordinary envelope, reading
% its vertices from the whole mesh rather than from a face.
    [fN, fD] = obj.faceQ(1);
    Hn = [fN(5) fN(6); fN(6) fN(7)];
    Ln = [fN(8); fN(9)];
    kn = fN(10);
    [Vi, vd] = ratQ.combineDen(obj.VN, obj.VD);

    if isempty(obj.E)
        % ---- a NEEDLE ------------------------------------------------------------------------
        p  = Vi(1,:).';
        qp = ratQ.chk(p.'*Hn*p + 2*vd*(Ln.'*p) + 2*vd^2*kn, 'needle value');   % over 2*fD*vd^2
        r1 = [ratQ.chk(vd,'e1'), 0, ratQ.chk(-p(1),'e0')];
        r2 = [0, ratQ.chk(vd,'e1'), ratQ.chk(-p(2),'e0')];
        con = [ratQ.conic([0 0 0, r1]), 0; ratQ.conic([0 0 0, r2]), 0];
        h = assembleQuaConCells(struct('num', [0 0 0 0 0 0 0 0 0 qp], ...
                                       'den', ratQ.chk(2*fD*vd^2,'needle denominator'), ...
                                       'con', con));
        return
    end

    if size(obj.E,1) ~= 1
        error('PLQ:biconjQ:chainNotImplemented', ...
            ['a CHAIN of %d edges is not thin: the convex hull of two non-collinear segments is ' ...
             'a two-dimensional polygon, so its envelope is the ordinary lower-hull case read ' ...
             'from the whole mesh rather than from a face.'], size(obj.E,1));
    end
    if obj.E(1,3) == 0
        error('PLQ:biconjQ:unbounded', ...
            ['the domain is a RAY. Where q curves down along it the envelope runs to -infinity, ' ...
             'which is a correct answer with nowhere to be stored -- the gap conjCPLQ records as ' ...
             'convEnvUnbounded:minusInfinity.']);
    end

    % ---- a SEGMENT --------------------------------------------------------------------------
    a = Vi(obj.E(1,1),:).';  b = Vi(obj.E(1,2),:).';
    d = ratQ.chk(b - a, 'segment direction');
    alpha = ratQ.chk(d.' * Hn * d, 'segment curvature');

    n  = [d(2); -d(1)];
    c0 = ratQ.chk(-(n.' * a), 'line offset');
    rows = [ratQ.chk(vd*n(1),'l1'), ratQ.chk(vd*n(2),'l2'), c0];          % ON the line
    lo   = [ratQ.chk(vd*d(1),'s1'), ratQ.chk(vd*d(2),'s2'), ratQ.chk(-(d.'*a),'s0')];
    hi   = [ratQ.chk(-vd*d(1),'t1'), ratQ.chk(-vd*d(2),'t2'), ratQ.chk(d.'*b,'t0')];
    con  = [ratQ.conic([0 0 0, rows]), 0; ...
            ratQ.conic([0 0 0, lo]), sgnOfB(lo, +1); ...
            ratQ.conic([0 0 0, hi]), sgnOfB(hi, +1)];

    if alpha >= 0
        % q is convex along the segment, so it is already its own envelope there.
        num = fN;  den = fD;
    else
        % q curves DOWN along the segment: the envelope is the CHORD, the affine interpolant of the
        % two endpoint values. Written over 2*fD*vd^2*(d.d), and checked to reproduce q(a) and q(b)
        % exactly at the ends.
        qa = ratQ.chk(a.'*Hn*a + 2*vd*(Ln.'*a) + 2*vd^2*kn, 'endpoint value');
        qb = ratQ.chk(b.'*Hn*b + 2*vd*(Ln.'*b) + 2*vd^2*kn, 'endpoint value');
        dd = ratQ.chk(d.'*d, 'segment length squared');
        num = [0 0 0 0, 0, 0, 0, ...
               ratQ.chk((qb-qa)*vd*d(1), 'c8'), ratQ.chk((qb-qa)*vd*d(2), 'c9'), ...
               ratQ.chk(qa*dd - (qb-qa)*(a.'*d), 'c10')];
        den = ratQ.chk(2*fD*vd^2*dd, 'chord denominator');
    end
    h = assembleQuaConCells(struct('num', num, 'den', den, 'con', con));
end
