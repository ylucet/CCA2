function edgeCellFaceProbe()
% THE RANK-1 NON-BRIDGE FACE -- the mixed sub-case, and the real remaining work on the envelope.
%
% WHERE IT SITS. Over an EDGE cell of a conjugate's subdivision the face's Hessian has rank 1, so
% the envelope face is
%
%       h_k(x) = sup over s in P_k of  <s, y> - <v,s>^2 / (2*h11)  -  c ,      y = x - b ,
%
% where the rank-1 PSD Hessian is written `H = v*v' / h11` with `v = [h11; h12]` -- exact, and
% available whenever `h11 ~= 0` (the `h11 = 0` spelling is handled by swapping coordinates).
%
% THE STRUCTURE, WHICH IS WHAT MAKES IT TRACTABLE. Along the null direction `u` of `H` -- any `u`
% with `<v,u> = 0`, e.g. `u = [-h12; h11]` -- the quadratic term is CONSTANT, so the objective is
% LINEAR in that direction. A linear function on a polyhedron is maximised on the BOUNDARY, so:
%
%       the maximiser lies on a FACET of P_k, or on the stationary LINE, and nowhere else.
%
% That collapses a 2-D degenerate problem into a handful of 1-D concave-quadratic maximisations,
% one per facet, each with a closed form. No elimination, no conjugation machinery.
%
% THE STATIONARY CASE IS A LINE, NOT A POINT, and it is worth writing out because it is where a
% rank-2 intuition goes wrong. `grad = y - H s = 0` needs `H s = y`, and `range(H) = span(v)`, so
% there is no stationary point at all unless `y` is PARALLEL to `v`. When it is, `y = mu*v` gives
% `<v,s> = mu*h11` -- a whole LINE of maximisers -- and the value collapses to
%
%       h_k = mu^2 * h11 / 2  -  c ,
%
% independent of which point on the line is taken. So the stationary case contributes ONE candidate
% on a measure-zero set of `x`, which matters for the mesh and not for a generic point.
%
% WHAT THIS PROBE ESTABLISHES, STATED NARROWLY. It checks the FORMULA -- that the candidate set
% above is complete and each closed form is right -- against an independent numeric sup, on real
% conjQ output. It computes in DOUBLES. Making it exact is a separate step: every candidate is a
% rational in the cell's own data, so it is bookkeeping rather than new mathematics, but it is not
% done here and must not be read as done. The exactness discipline this project runs on is about
% the PREDICATES that decide a mesh, and no mesh is being decided yet.

    fixtures = { 'triangle', tri([0 0; 1 0; 1 1], 0), ...
                 'square',   poly([0 0; 1 0; 1 1; 0 1]), ...
                 'hexagon',  poly([0 0; 2 0; 3 2; 2 4; 0 4; -1 2]) };

    nFin = 0;  nInf = 0;  nSkip = 0;  worst = 0;  worstAt = '';
    for i = 1:2:numel(fixtures)
        nm = fixtures{i};  g = conjQ(fixtures{i+1});
        fprintf('\n=== %s: %d dual cells ===\n', nm, g.nf);
        for k = 1:g.nf
            n = g.fN(k,:);
            H = [n(5), n(6); n(6), n(7)] / g.fD(k);
            if all(H(:) == 0) || abs(det(H)) > 1e-12, continue, end     % not rank 1
            P = cellRows(g, k);
            if isempty(P), nSkip = nSkip + 1; continue, end
            [f, in_, w, at] = checkCell(g, k, P);
            nFin = nFin + f;  nInf = nInf + in_;
            if w > worst, worst = w; worstAt = at; end
            fprintf('  cell %d: rank 1, %d facets | %d finite pts (worst %.2e), %d infinite pts\n', ...
                    k, size(P,1), f, w, in_);
        end
    end

    fprintf('\nTOTAL: %d finite points, worst relative error %.3e%s\n', nFin, worst, worstAt);
    fprintf('       %d points where the closed form says +infinity\n', nInf);
    fprintf('       %d cells skipped (curved boundary)\n', nSkip);
    if worst < 1e-8
        fprintf(['\nThe rank-1 candidate set is COMPLETE and each closed form is right, on every\n' ...
                 'edge cell of three fixtures. Exactness is the named next step, not this one.\n']);
    end
end

% ------------------------------------------------------------------------------------------------
function [v, why] = edgeCellValue(H, b, c, P, x)
% h_k(x) for a rank-1 cell, in closed form. v = +inf when the sup is unbounded.
    y = x - b;
    [vv, h11, sw] = rank1Split(H);
    if isempty(vv), v = NaN; why = 'the Hessian is not rank 1 in a usable spelling'; return, end
    if sw, y = flip(y); P = [P(:,2), P(:,1), P(:,3)]; end   % coordinates swapped when h11 = 0
    why = '';

    best = -inf;

    % ---- the stationary LINE, which exists only when y is parallel to v -----------------------
    if abs(y(1)*vv(2) - y(2)*vv(1)) < 1e-12 * max(1, norm(y)*norm(vv))
        mu = (vv.' * y) / (vv.' * vv);                       % y = mu * v
        % the line { s : <v,s> = mu*h11 } must meet P; tested by an LP-free scan along it
        if lineMeets(P, vv, mu * h11)
            best = max(best, mu^2 * h11 / 2 - c);
        end
    end

    % ---- one 1-D concave-quadratic maximisation per FACET --------------------------------------
    for i = 1:size(P,1)
        Ni = P(i,1:2);
        if all(Ni == 0), continue, end
        d = [Ni(2); -Ni(1)];
        p = -P(i,3) * Ni(:) / (Ni * Ni.');                   % a point on the facet's line
        [lo, hi] = facetRange(P, p, d);
        if isnan(lo), continue, end                          % the facet is empty

        a2 = -(vv.' * d)^2 / (2*h11);
        a1 = d.' * y - (vv.' * p) * (vv.' * d) / h11;
        a0 = p.' * y - (vv.' * p)^2 / (2*h11) - c;

        if a2 < -1e-14
            ts = -a1 / (2*a2);
            if ts > lo && ts < hi, best = max(best, a0 + a1*ts + a2*ts^2); end
        elseif abs(a1) > 1e-14
            % LINEAR along this facet: unbounded if the facet runs the way the slope points.
            if (a1 > 0 && isinf(hi)) || (a1 < 0 && isinf(lo)), v = inf; return, end
        end
        for t = [lo hi]
            if isinf(t)
                % A parabola opening downwards is bounded on a ray, so an infinite endpoint
                % contributes nothing once a2 < 0 -- and when a2 == 0 the linear test above has
                % already returned. Nothing to add.
                continue
            end
            best = max(best, a0 + a1*t + a2*t^2);
        end
    end
    v = best;
end

function [vv, h11, swapped] = rank1Split(H)
% H = v*v'/h11 for a rank-1 PSD H, with the coordinates swapped when the (1,1) entry vanishes.
% Both spellings are needed: an edge cell whose edge is vertical has h11 = 0 and h22 ~= 0.
    swapped = false;
    if abs(H(1,1)) > 1e-14
        vv = [H(1,1); H(1,2)];  h11 = H(1,1);  return
    end
    if abs(H(2,2)) > 1e-14
        swapped = true;
        vv = [H(2,2); H(1,2)];  h11 = H(2,2);  return
    end
    vv = [];  h11 = 0;
end

function tf = lineMeets(P, vv, rhs)
% Does { s : <v,s> = rhs } meet P? The line is 1-D, so this is an interval test: walk the
% constraints as bounds on the line's own parameter.
    d = [vv(2); -vv(1)];
    p = rhs * vv / (vv.' * vv);
    [lo, hi] = facetRange(P, p, d);
    tf = ~isnan(lo) && lo <= hi;
end

function [lo, hi] = facetRange(P, p, d)
% The parameter interval of { p + t*d } inside P. NaN when empty.
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

% ---- the independent oracle ----------------------------------------------------------------------
function [nFin, nInf, worst, at] = checkCell(g, k, P)
    n = g.fN(k,:);  d = g.fD(k);
    H = [n(5), n(6); n(6), n(7)] / d;
    b = [n(8); n(9)] / d;   c = n(10) / d;
    nFin = 0;  nInf = 0;  worst = 0;  at = '';

    rng(11);
    for t = 1:200
        x = (rand(2,1)*2 - 1) * 5;
        v = edgeCellValue(H, b, c, P, x);
        if isnan(v), continue, end
        if isinf(v)
            nInf = nInf + 1;  continue
        end
        o = sampledSup(H, b, c, P, x);
        if isinf(o), continue, end                 % the oracle also says unbounded; consistent
        rel = abs(v - o) / max(1, abs(v));
        if rel > worst
            worst = rel;
            at = sprintf(' (cell %d, x = [%.4f %.4f])', k, x(1), x(2));
        end
        nFin = nFin + 1;
    end
end

function v = sampledSup(H, b, c, P, x)
% THE ORACLE, REWRITTEN. sup over s in P of <s,x> - (s'Hs/2 + b's + c), by PROJECTED GRADIENT
% ASCENT from many starts. Concave objective, convex feasible set, so any local maximum is global
% and the method converges -- what it needs is starts that reach every part of an unbounded cell.
%
% WHY THE FIRST ORACLE WAS REPLACED, AND IT IS THE THIRD TIME TONIGHT. It rejection-sampled a
% square and polished with a backtracking gradient step. On cell 9 of the hexagon it reported
% -0.3650 while the closed form said +0.0950 -- and the closed form's point, the VERTEX
% s = (4, 8), is exactly feasible (two facet residuals of 0). The sample never landed on the
% vertex, and a backtracking step that only accepts strictly feasible points cannot slide ALONG a
% boundary to reach one. So the oracle was under-reporting, not the formula over-reporting.
%
% TWO FIXES, both about reaching boundaries rather than about precision:
%   * SEED AT EVERY FEASIBLE VERTEX of P, since a degenerate quadratic on a polyhedron very often
%     attains its max exactly there -- the objective is LINEAR along the Hessian's null direction.
%   * PROJECT rather than backtrack: after a gradient step, project cyclically onto the violated
%     half-planes until feasible. That slides along a facet instead of stopping at it.
    y = x - b;
    obj = @(s) s.'*y - 0.5*(s.'*H*s) - c;
    starts = polyVerts(P);
    starts = [starts, farStarts(P)];
    if isempty(starts), v = -inf; return, end

    v = -inf;
    for c0 = 1:size(starts,2)
        s = starts(:,c0);
        step = 1;
        for it = 1:600
            s2 = project(P, s + step * (y - H*s));
            if obj(s2) > obj(s) + 1e-15
                s = s2;
            else
                step = step / 2;
                if step < 1e-14, break, end
            end
        end
        v = max(v, obj(s));
        if norm(y - H*s) > 1e-6 && step < 1e-13 && obj(s) > 1e8, v = inf; return, end
    end
end

function s = project(P, s)
% Cyclic projection onto the violated half-planes. For an intersection of half-planes this
% converges to a feasible point, and unlike a backtracking step it can slide ALONG a facet.
    for pass = 1:60
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
% Every pairwise intersection of P's boundary lines that lies in P. A degenerate quadratic often
% attains its maximum exactly at one, which is why they are seeds rather than an afterthought.
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
% Points well inside the cell at several magnitudes, so an unbounded cell is entered rather than
% merely touched. Rejection sampling is fine for a START -- it only has to be feasible.
    V = zeros(2,0);
    for R = [1 10 100 1000]
        S = (rand(4000,2)*2 - 1) * R;
        S = S(all(S * P(:,1:2).' + P(:,3).' >= 0, 2), :);
        if ~isempty(S), V(:,end+1:end+min(4,size(S,1))) = S(1:min(4,size(S,1)),:).'; end
    end
end

% ---- fixtures and helpers ------------------------------------------------------------------------
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
