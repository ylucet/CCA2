function out = sweepMaxQuaParCurvedSplit(seed, nCases, verbose)
% sweepMaxQuaParCurvedSplit  SEEDED, COMMITTED replacement for maxQuaPar's curved-split
%   validation sweep (SUPPORT_MATRIX.md section 0.1).
%
% objective: re-derive, reproducibly, the figures `maxQuaPar.m`'s header and `DESIGN.md` quote for
%   the "two adjacent sub-pieces of one domain" configuration -- how many random splits assemble,
%   how many hit each documented guard, and how exact the assembled results are at result
%   vertices, at straight-edge midpoints and at interior points.
%
% [input]  seed    : rng seed (REQUIRED)
%          nCases  : number of random quadrilaterals (default 120)
%          verbose : per-case line (default false)
% [output] out     : struct of counts, rates and the worst errors
%
%   >> out = sweepMaxQuaParCurvedSplit(20260802, 120)
%
% THE CONFIGURATION, which is not incidental. maxQuaPar's scoping caveat is that its operands must
% be the conjugates of two ADJACENT sub-pieces of one domain -- there, the other operand's face
% boundaries are tangent lines to the parabola (a conjugate is C1 where its pieces join), which is
% what makes clipping an arc either keep it whole or reduce it to a point. So each case here is a
% convex quadrilateral split by a diagonal into two triangles, f = x*y on both, and the operands
% are conjPieceCPLQ of each triangle -- exactly what Step 2 hands Step 3.
%
% GROUND TRUTH is closed-form, and owes nothing to the conjugate pipeline: for f(x) = <s,x> - x1*x2
% the Hessian is [0 -1; -1 0], indefinite, so f has no interior local maximum and its sup over a
% compact convex polygon is attained on the BOUNDARY. Restricted to an edge it is a univariate
% quadratic, maximized exactly at an endpoint or at its own stationary point. supOverQuad below is
% that, and nothing more.
%
% WHY THIS FILE EXISTS: the sweep behind "115 sampled splits, 85 assembled, 340/340 edge midpoints,
% 5100/5100 interior points" was a throwaway script whose seed was not recorded, so none of those
% numbers can be re-derived and section 0.1 forbids quoting them onward. This one is committed and
% seeded, so its numbers can be. They are NOT expected to equal the retired ones: the random
% quadrilateral population is regenerated here, not recovered.

    if nargin < 1 || isempty(seed)
        error('sweepMaxQuaParCurvedSplit:seed', ...
            'A seed is REQUIRED: an unseeded sweep is not a reproducible measurement.');
    end
    if nargin < 2 || isempty(nCases), nCases = 120; end
    if nargin < 3, verbose = false; end

    rng(seed);

    nSampled = 0; nAssembled = 0;
    guards = containers.Map('KeyType','char','ValueType','double');
    nVtx = 0; badVtx = 0; worstVtx = 0;
    nMid = 0; badMid = 0; worstMid = 0;
    nInt = 0; badInt = 0; worstInt = 0;

    for c = 1:nCases
        Q = randomConvexQuad();
        if isempty(Q), continue; end
        nSampled = nSampled + 1;

        T1 = Q([1 2 3],:);
        T2 = Q([1 3 4],:);
        E = [1 2 1; 2 3 1; 3 1 1]; F = [1 0; 1 0; 1 0];
        try
            g1 = conjPieceCPLQ(QuaPol(toCCW(T1), E, [0 1 0 0 0 0], F));
            g2 = conjPieceCPLQ(QuaPol(toCCW(T2), E, [0 1 0 0 0 0], F));
        catch ME
            guards = bump(guards, ['conj:' ME.identifier]);
            continue
        end
        try
            g = maxQuaPar(g1, g2);
        catch ME
            guards = bump(guards, ME.identifier);
            continue
        end
        nAssembled = nAssembled + 1;

        % (1) exactly AT every result vertex
        for iv = 1:g.nv
            v = g.V(iv,:);
            if any(~isfinite(v)), continue; end
            nVtx = nVtx + 1;
            e = relErr(g.eval(v), supOverQuad(Q, v));
            if e > 1e-9, badVtx = badVtx + 1; worstVtx = max(worstVtx, e); end
        end

        % (2) the midpoint of every STRAIGHT edge (an arc's midpoint is not on its chord)
        for j = 1:size(g.E,1)
            if g.E(j,3) ~= 1, continue; end
            if ~isempty(g.Ec) && any(g.Ec(j,:) ~= 0), continue; end
            A = g.V(g.E(j,1),:); B = g.V(g.E(j,2),:);
            if any(~isfinite([A B])), continue; end
            mp = (A + B)/2;
            nMid = nMid + 1;
            e = relErr(g.eval(mp), supOverQuad(Q, mp));
            if e > 1e-9, badMid = badMid + 1; worstMid = max(worstMid, e); end
        end

        % (3) generic interior points, over a box comfortably containing the arrangement
        fin = g.V(all(isfinite(g.V),2), :);
        if isempty(fin), fin = [0 0]; end
        rad = 1 + max(abs(fin(:)));
        Pts = (2*rand(60,2) - 1) * rad;
        for k = 1:size(Pts,1)
            nInt = nInt + 1;
            e = relErr(g.eval(Pts(k,:)), supOverQuad(Q, Pts(k,:)));
            if e > 1e-9, badInt = badInt + 1; worstInt = max(worstInt, e); end
        end

        if verbose
            fprintf('case %3d: assembled, nv=%d ne=%d\n', c, g.nv, size(g.E,1));
        end
    end

    out = struct();
    out.seed = seed; out.nCases = nCases;
    out.nSampled = nSampled; out.nAssembled = nAssembled;
    out.guards = guards;
    out.nVertexProbes = nVtx;   out.badVertex   = badVtx;   out.worstVertex   = worstVtx;
    out.nMidpointProbes = nMid; out.badMidpoint = badMid;   out.worstMidpoint = worstMid;
    out.nInteriorProbes = nInt; out.badInterior = badInt;   out.worstInterior = worstInt;

    fprintf(['sweepMaxQuaParCurvedSplit(seed=%d, nCases=%d)\n' ...
             '  sampled splits            : %d\n' ...
             '  assembled                 : %d\n'], seed, nCases, nSampled, nAssembled);
    ks = keys(guards);
    for i = 1:numel(ks)
        fprintf('  guard %-38s : %d\n', ks{i}, guards(ks{i}));
    end
    fprintf(['  AT a result vertex        : %d/%d disagree (%.2f%%), worst %.3g\n' ...
             '  straight-edge midpoints   : %d/%d disagree, worst %.3g\n' ...
             '  interior points           : %d/%d disagree, worst %.3g\n'], ...
        badVtx, nVtx, 100*badVtx/max(1,nVtx), worstVtx, ...
        badMid, nMid, worstMid, badInt, nInt, worstInt);
end

% ================================================================================================
function e = relErr(got, want)
    if ~isfinite(got)
        e = inf;
        return
    end
    e = abs(got - want) / max(1, abs(want));
end

% ================================================================================================
function v = supOverQuad(Q, s)
% sup over the convex polygon Q (CCW) of <s,x> - x1*x2, in closed form. The objective's Hessian is
% [0 -1; -1 0] -- indefinite -- so it has no interior local maximum and the sup over a compact
% convex set is attained on the boundary; on an edge it is a univariate quadratic.
    v = -inf;
    m = size(Q,1);
    for i = 1:m
        A = Q(i,:); B = Q(mod(i,m)+1,:); d = B - A;
        a2 = -d(1)*d(2);
        a1 = s(1)*d(1) + s(2)*d(2) - A(1)*d(2) - A(2)*d(1);
        a0 = s(1)*A(1) + s(2)*A(2) - A(1)*A(2);
        ts = [0 1];
        if a2 < 0
            tst = -a1/(2*a2);
            if tst > 0 && tst < 1, ts(end+1) = tst; end %#ok<AGROW>
        end
        for t = ts
            v = max(v, a2*t^2 + a1*t + a0);
        end
    end
end

% ================================================================================================
function Q = randomConvexQuad()
% Four points in convex position, CCW. Returns [] when the random draw's hull is not a
% quadrilateral (three points, or one inside the triangle of the others).
    P = randn(4,2);
    k = convhull(P(:,1), P(:,2));
    k(end) = [];
    if numel(k) ~= 4
        Q = [];
        return
    end
    Q = P(k,:);
    if signedArea(Q) < 0, Q = flipud(Q); end
end

% ================================================================================================
function T = toCCW(T)
    a = (T(2,1)-T(1,1))*(T(3,2)-T(1,2)) - (T(2,2)-T(1,2))*(T(3,1)-T(1,1));
    if a < 0, T = T([1 3 2],:); end
end

% ================================================================================================
function a = signedArea(V)
    n = size(V,1); a = 0;
    for i = 1:n
        j = mod(i,n)+1;
        a = a + V(i,1)*V(j,2) - V(j,1)*V(i,2);
    end
    a = a/2;
end

% ================================================================================================
function m = bump(m, key)
    if isKey(m, key), m(key) = m(key) + 1; else, m(key) = 1; end
end
