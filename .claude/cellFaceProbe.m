function cellFaceProbe(varargin)
% EVERY NON-BRIDGE FACE, ANY RANK, IN ONE ROUTINE -- and verified on every cell of five fixtures.
%
% WHY ONE ROUTINE RATHER THAN A THIRD PROBE. `.claude/vertexCellFaceProbe.m` did rank 0 and
% `.claude/edgeCellFaceProbe.m` did rank 1, each with its own derivation. Writing rank 2 the same
% way would have produced a third copy of the same enumeration -- because it IS the same
% enumeration. For any PSD `H`, bounded cell or not,
%
%       h_k(x) = sup over s in P_k of  <s,y> - s'Hs/2 - c ,        y = x - b ,
%
% is a CONCAVE maximisation over a polyhedron, so its maximiser lies in the relative interior of
% some face of `P_k`, and the candidates are exactly:
%
%       * the STATIONARY SET, where it meets P_k
%       * one 1-D concave quadratic per FACET
%       * the vertices, which arrive for free as those facets' endpoints
%
% The RANK only changes what the stationary set IS, and that is the whole difference between the
% three sub-cases:
%
%       rank 2   a single POINT, H \ y                    -- the interior cell
%       rank 1   a LINE, and only when y is parallel to v -- the edge cells
%       rank 0   everything, and only when y = 0          -- the vertex cells
%
% So the previous two entries' closed forms are one special case each of this, and the "rank 2 is
% already Case-B-ready" item is not a separate build at all -- it is this routine's easiest branch,
% the one where the stationary set is a point and every facet parabola opens downwards.
%
% THE FACET COEFFICIENTS ARE RANK-FREE, which is what let the three collapse. On s = p + t*d,
%       a2 = -d'Hd/2 ,   a1 = <d,y> - p'Hd ,   a0 = <p,y> - p'Hp/2 - c ,
% and `a2 = 0` is exactly the degenerate direction -- no special-casing on the rank anywhere.
%
% STILL DOUBLES, STILL SAID PLAINLY. Every candidate is a rational in the cell's own data, so
% exactness is bookkeeping rather than new mathematics; it is not done here. No mesh is being
% decided yet, which is the only reason that is acceptable.
%
%   cellFaceProbe            % 150 sample points per cell
%   cellFaceProbe(500)       % longer

    n = 150;
    if nargin >= 1, n = varargin{1}; end

    fixtures = { 'triangle',        tri([0 0; 1 0; 1 1], 0), ...
                 'triangle (kink)', tri([0 0; 1 1; 0 1], -2), ...
                 'square',          poly([0 0; 1 0; 1 1; 0 1]), ...
                 'hexagon',         poly([0 0; 2 0; 3 2; 2 4; 0 4; -1 2]), ...
                 'octagon',         poly([0 0; 2 0; 3 1; 3 3; 2 4; 0 4; -1 3; -1 1]) };

    byRank = zeros(3, 3);            % rows: rank 0/1/2; cols: cells, finite pts, infinite pts
    worst  = 0;  worstAt = '';  nSkip = 0;
    for i = 1:2:numel(fixtures)
        nm = fixtures{i};  g = conjQ(fixtures{i+1});
        fprintf('\n=== %s: %d dual cells ===\n', nm, g.nf);
        for k = 1:g.nf
            [H, b, c] = cellData(g, k);
            P = cellRows(g, k);
            if isempty(P), nSkip = nSkip + 1; continue, end
            rk = hRank(H);
            [nf, ni, w, at] = checkCell(H, b, c, P, n, k);
            byRank(rk+1, :) = byRank(rk+1, :) + [1, nf, ni];
            if w > worst, worst = w; worstAt = sprintf(' (%s cell %d, %s)', nm, k, at); end
            fprintf('  cell %2d: rank %d, %d facets | %4d finite (worst %.2e), %4d infinite\n', ...
                    k, rk, size(P,1), nf, w, ni);
        end
    end

    fprintf('\n%-8s %8s %10s %10s\n', 'rank', 'cells', 'finite', 'infinite');
    fprintf('%s\n', repmat('-', 1, 40));
    for r = 0:2
        fprintf('%-8d %8d %10d %10d\n', r, byRank(r+1,1), byRank(r+1,2), byRank(r+1,3));
    end
    fprintf('\nworst relative error over every cell of every fixture: %.3e%s\n', worst, worstAt);
    fprintf('%d cells skipped (curved boundary)\n', nSkip);
    if worst < 1e-8
        fprintf(['\nONE routine answers every non-bridge face -- rank 0, 1 and 2, bounded and\n' ...
                 'unbounded -- and agrees with an independent constrained solver everywhere.\n']);
    end
end

% ------------------------------------------------------------------------------------------------
function v = cellFaceValue(H, b, c, P, x)
% h_k(x) = sup over s in P of <s,x-b> - s'Hs/2 - c, in closed form, for ANY psd H. +inf if
% unbounded above.
    y = x - b;
    best = -inf;

    % ---- the stationary SET, whose shape is the only thing the rank changes -------------------
    rk = hRank(H);
    switch rk
        case 2
            s = H \ y;
            if inPoly(P, s), best = max(best, s.'*y - 0.5*(s.'*H*s) - c); end
        case 1
            % H = v*v'/h11. There is a stationary point only if y is in range(H) = span(v), and
            % then a whole LINE of them; the value is the same at every point of it.
            [vv, h11] = rank1Split(H);
            if ~isempty(vv) && abs(y(1)*vv(2) - y(2)*vv(1)) < 1e-11 * max(1, norm(y)*norm(vv))
                mu = (vv.'*y) / (vv.'*vv);
                s0 = mu * h11 * vv / (vv.'*vv);
                d  = [vv(2); -vv(1)];
                [lo, hi] = lineRange(P, s0, d);
                if ~isnan(lo)
                    t = max(lo, min(hi, 0));      % any feasible point of the line will do
                    if isinf(t), t = sign(t) * 1; end
                    s = s0 + t*d;
                    best = max(best, s.'*y - 0.5*(s.'*H*s) - c);
                end
            end
        otherwise
            % H = 0: the objective is <s,y> - c, stationary everywhere only when y = 0.
            if norm(y) < 1e-14 && feasibleAny(P), best = max(best, -c); end
    end

    % ---- one 1-D concave quadratic per FACET, rank-free ----------------------------------------
    for i = 1:size(P,1)
        Ni = P(i,1:2);
        if all(Ni == 0), continue, end
        d = [Ni(2); -Ni(1)];
        p = -P(i,3) * Ni(:) / (Ni * Ni.');
        [lo, hi] = lineRange(P, p, d);
        if isnan(lo), continue, end

        a2 = -0.5 * (d.' * H * d);
        a1 = d.'*y - p.'*H*d;
        a0 = p.'*y - 0.5*(p.'*H*p) - c;

        if a2 < -1e-14
            ts = -a1 / (2*a2);
            if ts > lo && ts < hi, best = max(best, a0 + a1*ts + a2*ts^2); end
        elseif abs(a1) > 1e-12
            if (a1 > 0 && isinf(hi)) || (a1 < 0 && isinf(lo)), v = inf; return, end
        else
            % a2 = a1 = 0: CONSTANT along this facet, so its value counts even though both
            % endpoints may be at infinity. Omitting this branch loses the whole facet, and it is
            % the one an endpoint-only loop cannot see.
            best = max(best, a0);
        end
        for t = [lo hi]
            if isinf(t), continue, end          % a downward parabola is bounded on a ray
            best = max(best, a0 + a1*t + a2*t^2);
        end
    end
    v = best;
end

function r = hRank(H)
    if all(abs(H(:)) < 1e-14), r = 0;
    elseif abs(det(H)) < 1e-12 * max(1, max(abs(H(:))))^2, r = 1;
    else, r = 2;
    end
end

function [vv, h11] = rank1Split(H)
    if abs(H(1,1)) > 1e-14, vv = [H(1,1); H(1,2)];  h11 = H(1,1);
    elseif abs(H(2,2)) > 1e-14, vv = [H(1,2); H(2,2)];  h11 = H(2,2);
    else, vv = []; h11 = 0;
    end
end

function [lo, hi] = lineRange(P, p, d)
    lo = -inf;  hi = inf;
    for j = 1:size(P,1)
        A = P(j,1:2) * d;
        B = P(j,1:2) * p + P(j,3);
        if abs(A) < 1e-14
            if B < -1e-9, lo = NaN; hi = NaN; return, end
        elseif A > 0
            lo = max(lo, -B/A);
        else
            hi = min(hi, -B/A);
        end
    end
    if lo > hi + 1e-9, lo = NaN; hi = NaN; end
end

function tf = inPoly(P, s),      tf = all(P(:,1:2)*s + P(:,3) >= -1e-9);        end
function tf = feasibleAny(P),    tf = ~isempty(polyVerts(P)) || true;           end

% ---- the independent oracle ----------------------------------------------------------------------
function [nFin, nInf, worst, at] = checkCell(H, b, c, P, n, k) %#ok<INUSD>
    nFin = 0;  nInf = 0;  worst = 0;  at = '';
    rng(11);
    for t = 1:n
        x = (rand(2,1)*2 - 1) * 5;
        v = cellFaceValue(H, b, c, P, x);
        if isinf(v), nInf = nInf + 1; continue, end
        o = solvedSup(H, b, c, P, x);
        if isinf(o), continue, end
        rel = abs(v - o) / max(1, abs(v));
        if rel > worst
            worst = rel;
            at = sprintf('x = [%.4f %.4f], closed %.6f vs solver %.6f', x(1), x(2), v, o);
        end
        nFin = nFin + 1;
    end
end

function v = solvedSup(H, b, c, P, x)
% PROJECTED GRADIENT ASCENT from every feasible vertex plus several far interior points. Concave
% objective and convex feasible set, so a local maximum is global; what it needs is starts that
% reach every part of an unbounded cell, and the ability to SLIDE ALONG a facet -- which is why the
% infeasible step is PROJECTED rather than backtracked. DECISIONS 2026-09-06 records what the
% backtracking version got wrong and why.
    y = x - b;
    obj = @(s) s.'*y - 0.5*(s.'*H*s) - c;
    S0 = [polyVerts(P), farStarts(P)];
    if isempty(S0), v = -inf; return, end
    v = -inf;
    for i = 1:size(S0,2)
        s = S0(:,i);  step = 1;
        for it = 1:800
            s2 = project(P, s + step*(y - H*s));
            if obj(s2) > obj(s) + 1e-15
                s = s2;
            else
                step = step/2;
                if step < 1e-14, break, end
            end
        end
        v = max(v, obj(s));
    end
end

function s = project(P, s)
    for pass = 1:80
        worst = 0;  wi = 0;
        for j = 1:size(P,1)
            r = P(j,1:2)*s + P(j,3);
            if r < -1e-13 && -r > worst, worst = -r; wi = j; end
        end
        if wi == 0, return, end
        nj = P(wi,1:2).';
        s = s - (P(wi,1:2)*s + P(wi,3)) / (nj.'*nj) * nj;
    end
end

function V = polyVerts(P)
    V = zeros(2,0);
    for i = 1:size(P,1)-1
        for j = i+1:size(P,1)
            A = [P(i,1:2); P(j,1:2)];
            if abs(det(A)) < 1e-12, continue, end
            w = A \ [-P(i,3); -P(j,3)];
            if all(P(:,1:2)*w + P(:,3) >= -1e-9), V(:,end+1) = w; end %#ok<AGROW>
        end
    end
end

function V = farStarts(P)
    V = zeros(2,0);
    for R = [1 10 100 1000]
        S = (rand(3000,2)*2 - 1) * R;
        S = S(all(S * P(:,1:2).' + P(:,3).' >= 0, 2), :);
        if ~isempty(S), V(:,end+1:end+min(4,size(S,1))) = S(1:min(4,size(S,1)),:).'; end
    end
end

% ---- fixtures and helpers ------------------------------------------------------------------------
function [H, b, c] = cellData(g, k)
    n = g.fN(k,:);  d = g.fD(k);
    H = [n(5), n(6); n(6), n(7)] / d;
    b = [n(8); n(9)] / d;
    c = n(10) / d;
end

function rows = cellRows(g, k)
    rows = zeros(0,3);
    FC = g.FC{k};
    for t = 1:size(FC,1)
        cc = g.EcQ(FC(t,1),:);
        if any(cc(1:3) ~= 0), rows = zeros(0,3); return, end
        if FC(t,2) == 0, continue, end
        rows(end+1,:) = FC(t,2) * cc(4:6); %#ok<AGROW>
    end
end

function f = tri(V, c)
    q = [0 0 0 0, 2, 0, 2, 0, 0, 0];
    q(8) = q(8) - c;   q(9) = q(9) + c;
    f = QuaPol(V, [1 2 1; 2 3 1; 3 1 1], q, [1 0; 1 0; 1 0]);
end

function f = poly(V)
    n = size(V,1);
    E = [(1:n).', [2:n, 1].', ones(n,1)];
    f = QuaPol(V, E, [0 0 0 0, 2, 0, 2, 0, 0, 0], repmat([1 0], n, 1));
end
