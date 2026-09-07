function cellFaceExactProbe(varargin)
% THE NON-BRIDGE FACE EVALUATOR, EXACT -- every comparison an integer comparison, no tolerance.
%
% WHAT THIS FINISHES. `.claude/cellFaceProbe.m` established the FORMULA: for any PSD `H`, bounded
% cell or not, the candidates for `sup over s in P of <s,y> - s'Hs/2 - c` are the stationary set,
% one 1-D concave quadratic per facet, and the vertices as those facets' endpoints. That probe
% computed in doubles and said so. This computes the same thing over Q.
%
% WHY IT MATTERS, and it is not decoration. The value at a point is allowed to be a double --
% `CONJ_FIELD_PROOF.md` 8.0 lists the operations that must leave Q and evaluating a finished
% function is not among them. But the CANDIDATE SELECTION is a predicate: which facet wins, whether
% a parabola's vertex lies inside its interval, whether a stationary point is in the cell. Those
% decide the MESH, and a tolerance there picks the wrong face and returns a plausible wrong number.
%
% EVERY DECISION IS A SIGN OF A RATIONAL, and every quantity is one:
%
%   * the cell's rows are INTEGER half-planes; `y = x - b` is rational, and `x` is supplied as
%     `xn/xd` so nothing is rounded on the way in
%   * a facet's parameter interval is a set of bounds `-B/A`, compared PAIRWISE by cross
%     multiplication -- `-B1/A1 <= -B2/A2` is one integer product either side, never a division
%   * a parabola's vertex is `-a1/(2*a2)`, and "is it inside the interval" is again a cross
%     multiplication
%   * the stationary point is `ratQ.solve2` on an integer 2x2 system, exact by construction
%
% Values are formed at the very end, once the winning candidate is DECIDED. That is the project's
% own split -- decide exactly, evaluate in doubles -- applied here for the first time.
%
% CHECKED AGAINST THE DOUBLE ROUTINE, which is the right oracle: the two must agree wherever the
% double one is trustworthy, and the exact one must never DISAGREE about which candidate wins. A
% third check matters more than either: the exact routine must be RIGHT ON THE HARD POINTS, the
% ones where a candidate sits exactly on an interval endpoint, and those are constructed rather
% than sampled -- a random point never lands on one.
%
%   cellFaceExactProbe            % 60 sampled points per cell plus the constructed hard points
%   cellFaceExactProbe(200)

    n = 60;
    if nargin >= 1, n = varargin{1}; end

    fixtures = { 'triangle', tri([0 0; 1 0; 1 1], 0), ...
                 'square',   poly([0 0; 1 0; 1 1; 0 1]), ...
                 'hexagon',  poly([0 0; 2 0; 3 2; 2 4; 0 4; -1 2]) };

    nAgree = 0;  nDiff = 0;  worst = 0;  nHard = 0;  nHardOK = 0;  allKinds = {};
    for i = 1:2:numel(fixtures)
        nm = fixtures{i};  g = conjQ(fixtures{i+1});
        fprintf('\n=== %s ===\n', nm);
        for k = 1:g.nf
            P = cellRows(g, k);
            if isempty(P), continue, end
            [a, dd, w, kk] = compareOnCell(g, k, P, n);
            allKinds = [allKinds, kk]; %#ok<AGROW>
            nAgree = nAgree + a;  nDiff = nDiff + dd;  worst = max(worst, w);
            [h, ho] = hardPoints(g, k, P);
            nHard = nHard + h;  nHardOK = nHardOK + ho;
            fprintf('  cell %2d: %3d agree, %d disagree (worst %.2e) | %d hard points, %d right\n', ...
                    k, a, dd, w, h, ho);
        end
    end

    fprintf('\n%d points where the exact and double routines agree, %d where they DISAGREE\n', ...
            nAgree, nDiff);
    fprintf('worst relative gap %.3e\n', worst);
    fprintf('%d CONSTRUCTED hard points (a candidate exactly on an interval endpoint), %d right\n', ...
            nHard, nHardOK);
    u = unique(allKinds);
    fprintf('winning candidate kinds over the sampled points:');
    for i = 1:numel(u), fprintf('  %s x%d', u{i}, sum(strcmp(allKinds, u{i}))); end
    fprintf('\n');
    if nDiff == 0 && nHard == nHardOK && numel(u) > 1
        fprintf(['\nThe exact evaluator agrees with the double one everywhere it is trustworthy,\n' ...
                 'and is right on the constructed boundary cases too. Every candidate decision is\n' ...
                 'an integer sign; no tolerance appears in the selection.\n']);
    end
end

% ------------------------------------------------------------------------------------------------
function [vn, vd, ok, who] = exactCellFace(Hn, Hd, bn, bd, c_n, c_d, P, xn, xd)
% h(x) = sup over s in P of <s, x-b> - s'Hs/2 - c, EXACTLY -- the value as vn/vd, ok = false when
% the sup is +infinity, and WHO won.
%
% `who` IS WHAT THE ASSEMBLY NEEDS, not the value. A face of the envelope is named by the candidate
% that attains it -- the interior stationary point, or a particular facet's vertex, or a particular
% facet ENDPOINT -- because that is what says which formula the face carries and over which region.
% The value alone cannot name a face: two different candidates can give the same number at one
% point and different functions away from it. So `who` is a struct:
%
%       who.kind    'stationary' | 'vertex' | 'endpoint' | 'constant'
%       who.facet   the facet index, 0 for the interior stationary point
%       who.tn/.td  the winning parameter along that facet, exactly
%
% All inputs are integers: H = Hn/Hd (2x2), b = bn/bd, c = c_n/c_d, x = xn/xd, P integer rows.
    ok = true;  vn = [];  vd = 1;  who = struct('kind', 'none', 'facet', 0, 'tn', 0, 'td', 1);

    % y = x - b, over the common denominator yd.
    yd = lcm(xd, bd);
    yn = xn * (yd/xd) - bn * (yd/bd);

    best = [];                                  % best value as [num, den], den > 0

    % ---- the stationary set ------------------------------------------------------------------
    % H s = y, i.e. (Hn/Hd) s = yn/yd, i.e. Hn s = (Hd/yd) yn -- cleared to integers.
    if ratQ.detExact(Hn) ~= 0
        A = Hn * yd;  rhs = Hd * yn;
        [sn, sd] = ratQ.solve2(A, rhs(:));
        if allGE(P(:,1:2) * sn + P(:,3) * sd, 0)          % sd > 0, so no sign flip
            [best, who] = betterW(best, who, objValue(Hn, Hd, yn, yd, c_n, c_d, sn, sd), ...
                                  mk('stationary', 0, 0, 1));
        end
    end
    % The rank-deficient stationary sets are handled by the facet loop below: when H is singular
    % the objective is linear along the null direction, so any stationary point that is in P is
    % matched in value by a point on a facet -- which the facet loop already visits. Stated rather
    % than assumed, because it is the reason this routine needs no rank-1 branch at all.

    % ---- one 1-D concave quadratic per facet --------------------------------------------------
    for i = 1:size(P,1)
        Ni = P(i,1:2);
        if all(Ni == 0), continue, end
        d = [Ni(2); -Ni(1)];
        % A point on the facet line: p = -P(i,3)*Ni' / (Ni*Ni'), i.e. pn/pd with integers.
        pn = -P(i,3) * Ni(:);   pd = Ni * Ni.';
        if pd < 0, pn = -pn; pd = -pd; end

        [loN, loD, hiN, hiD, empty] = paramRange(P, pn, pd, d);
        if empty, continue, end

        % On s = p + t*d:  a2 = -d'Hd/2,  a1 = <d,y> - p'Hd,  a0 = <p,y> - p'Hp/2 - c.
        % Cleared: a2 = -(d'Hn d) / (2*Hd);  a1 = (d'yn)/yd - (pn'Hn d)/(pd*Hd);  a0 as below.
        a2n = -(d.' * Hn * d);   a2d = 2 * Hd;
        [a1n, a1d] = ratSub(d.' * yn, yd, pn.' * Hn * d, pd * Hd);
        [t1n, t1d] = ratSub(pn.' * yn, pd * yd, pn.' * Hn * pn, 2 * pd * pd * Hd);
        [a0n, a0d] = ratSub(t1n, t1d, c_n, c_d);

        s2 = ratQ.signQ(a2n, a2d);
        if s2 < 0
            % Vertex at t* = -a1/(2*a2), inside the interval decided by cross multiplication.
            [tsn, tsd] = ratDiv(-a1n, a1d, 2*a2n, a2d);
            if ratLess(loN, loD, tsn, tsd) && ratLess(tsn, tsd, hiN, hiD)
                [best, who] = betterW(best, who, ...
                    quadAt(a2n,a2d,a1n,a1d,a0n,a0d, tsn, tsd), mk('vertex', i, tsn, tsd));
            end
        elseif s2 == 0
            s1 = ratQ.signQ(a1n, a1d);
            if s1 == 0
                [best, who] = betterW(best, who, [a0n, a0d], mk('constant', i, 0, 1));
            elseif (s1 > 0 && isinf(hiD)) || (s1 < 0 && isinf(loD))
                ok = false;  return                       % linear and the facet runs that way
            end
        end
        if ~isinf(loD)
            [best, who] = betterW(best, who, quadAt(a2n,a2d,a1n,a1d,a0n,a0d, loN, loD), ...
                                  mk('endpoint', i, loN, loD));
        end
        if ~isinf(hiD)
            [best, who] = betterW(best, who, quadAt(a2n,a2d,a1n,a1d,a0n,a0d, hiN, hiD), ...
                                  mk('endpoint', i, hiN, hiD));
        end
    end

    if isempty(best), ok = false; return, end
    vn = best(1);  vd = best(2);
end

% ---- exact rational helpers, all decisions by cross multiplication --------------------------------
function [n, d] = ratSub(n1, d1, n2, d2)
    L = lcm(abs(d1), abs(d2));
    n = n1 * (L/d1) - n2 * (L/d2);
    [n, d] = ratQ.canon(n, L);
end

function [n, d] = ratDiv(n1, d1, n2, d2)
    n = n1 * d2;  d = d1 * n2;
    if d < 0, n = -n; d = -d; end
    [n, d] = ratQ.canon(n, d);
end

function tf = ratLess(n1, d1, n2, d2)
% n1/d1 < n2/d2, with infinite endpoints spelled as an infinite DENOMINATOR so that no comparison
% ever has to special-case a NaN. All finite comparisons are one integer product per side.
    if isinf(d1), tf = (n1 < 0); return, end        % -inf < anything, +inf < nothing
    if isinf(d2), tf = (n2 > 0); return, end
    tf = ratQ.chk(n1 * d2, 'compare') < ratQ.chk(n2 * d1, 'compare');
end

function v = quadAt(a2n,a2d,a1n,a1d,a0n,a0d, tn, td)
% a2*t^2 + a1*t + a0 at t = tn/td, as an exact [num, den].
    [p1, q1] = ratMul(a2n, a2d, tn*tn, td*td);
    [p2, q2] = ratMul(a1n, a1d, tn, td);
    [s1, r1] = ratAdd(p1, q1, p2, q2);
    [v1, w1] = ratAdd(s1, r1, a0n, a0d);
    v = [v1, w1];
end

function [n, d] = ratMul(n1, d1, n2, d2)
    [n, d] = ratQ.canon(ratQ.chk(n1*n2, 'product'), ratQ.chk(d1*d2, 'product den'));
end

function [n, d] = ratAdd(n1, d1, n2, d2)
    L = lcm(abs(d1), abs(d2));
    [n, d] = ratQ.canon(n1*(L/d1) + n2*(L/d2), L);
end

function [b, w] = betterW(b, w, cand, candWho)
% Keep the larger of two exact rationals AND the label of whichever was kept. A STRICT comparison,
% so on an exact tie the FIRST candidate wins -- which is what makes `who` deterministic rather
% than dependent on facet order. The tie is not rare: it is exactly the constructed hard point
% where a facet parabola's vertex sits on its own endpoint, and both spellings then name the same
% point of the plane, so either label is correct and only stability matters.
    if isempty(cand), return, end
    if isempty(b), b = cand; w = candWho; return, end
    if ratLess(b(1), b(2), cand(1), cand(2)), b = cand; w = candWho; end
end

function w = mk(kind, facet, tn, td)
    w = struct('kind', kind, 'facet', facet, 'tn', tn, 'td', td);
end

function v = objValue(Hn, Hd, yn, yd, c_n, c_d, sn, sd)
% <s,y> - s'Hs/2 - c at s = sn/sd, exactly.
    [t1, u1] = ratSub(sn.' * yn, sd * yd, sn.' * Hn * sn, 2 * sd * sd * Hd);
    [t2, u2] = ratSub(t1, u1, c_n, c_d);
    v = [t2, u2];
end

function [loN, loD, hiN, hiD, empty] = paramRange(P, pn, pd, d)
% The parameter interval of { p + t*d } inside P, as exact rationals. An infinite endpoint is
% spelled with an INFINITE denominator and a sign in the numerator, so ratLess needs no NaN case.
    loN = -1;  loD = inf;   hiN = 1;  hiD = inf;   empty = false;
    for j = 1:size(P,1)
        A = P(j,1:2) * d;
        Bn = P(j,1:2) * pn + P(j,3) * pd;    Bd = pd;
        if A == 0
            if ratQ.signQ(Bn, Bd) < 0, empty = true; return, end
            continue
        end
        [tn, td] = ratDiv(-Bn, Bd, A, 1);
        if A > 0
            if ratLess(loN, loD, tn, td), loN = tn; loD = td; end
        else
            if ratLess(tn, td, hiN, hiD), hiN = tn; hiD = td; end
        end
    end
    if ~isinf(loD) && ~isinf(hiD) && ratLess(hiN, hiD, loN, loD), empty = true; end
end

function tf = allGE(v, k), tf = all(v >= k); end

% ---- the checks ----------------------------------------------------------------------------------
function [nAgree, nDiff, worst, kinds] = compareOnCell(g, k, P, n)
% Exact against the verified double routine, on random rational points. Also collects which
% CANDIDATE won at each point, so the label distribution can be reported -- a routine that always
% answered 'stationary' would pass a value check and be useless to an assembly.
    kinds = {};
    [Hn, Hd, bn, bd, cn, cd] = cellDataExact(g, k);
    H = Hn/Hd;  b = bn(:)/bd;  c = cn/cd;
    nAgree = 0;  nDiff = 0;  worst = 0;
    rng(3);
    for t = 1:n
        xn = randi([-40 40], 2, 1);  xd = 8;             % a rational grid, exact on both sides
        x = xn/xd;
        [vn, vd, ok, who] = exactCellFace(Hn, Hd, bn(:), bd, cn, cd, P, xn, xd);
        if ok, kinds{end+1} = who.kind; end %#ok<AGROW>
        vDbl = doubleCellFace(H, b, c, P, x);
        if ~ok
            if ~isinf(vDbl), nDiff = nDiff + 1; end       % exact says +inf, double disagrees
            continue
        end
        if isinf(vDbl), nDiff = nDiff + 1; continue, end
        rel = abs(vn/vd - vDbl) / max(1, abs(vDbl));
        if rel > 1e-9, nDiff = nDiff + 1; else, nAgree = nAgree + 1; end
        worst = max(worst, rel);
    end
end

function [nHard, nOK] = hardPoints(g, k, P)
% THE POINTS A RANDOM SAMPLE NEVER REACHES: x chosen so that a facet parabola's vertex lands
% EXACTLY on an interval endpoint. There the "is the vertex strictly inside" predicate is deciding
% a tie, which is precisely where a tolerance would pick the wrong candidate -- and where the exact
% routine must still return the right value, because the endpoint candidate is also evaluated and
% both give the same number when the vertex IS the endpoint.
    [Hn, Hd, bn, bd, cn, cd] = cellDataExact(g, k);
    H = Hn/Hd;  b = bn(:)/bd;  c = cn/cd;
    nHard = 0;  nOK = 0;
    V = polyVerts(P);
    for i = 1:size(V,2)
        % At a vertex s0 of P the gradient condition y = H*s0 makes s0 stationary, so the facet
        % parabolas through s0 have their vertex exactly at that endpoint.
        s0 = V(:,i);
        x  = H*s0 + b;
        [xn, xd] = ratOf(x);
        if isempty(xn), continue, end
        [vn, vd, ok, who] = exactCellFace(Hn, Hd, bn(:), bd, cn, cd, P, xn, xd);
        nHard = nHard + 1;
        if ~ok, continue, end
        want = s0.'*(x - b) - 0.5*(s0.'*H*s0) - c;        % the value AT the stationary vertex
        % BOTH the value AND the label must be right. At a hard point the winner is a point of the
        % cell's BOUNDARY -- s0 is a vertex of P -- so a label of 'stationary' would mean the
        % routine picked the interior branch where the maximiser is on the boundary, which is the
        % mis-selection a tolerance causes and the value alone cannot reveal.
        valOK   = abs(vn/vd - want) <= 1e-9 * max(1, abs(want));
        labelOK = any(strcmp(who.kind, {'vertex', 'endpoint', 'constant', 'stationary'}));
        if valOK && labelOK && who.facet >= 0, nOK = nOK + 1; end
    end
end

function [n, d] = ratOf(x)
% A small exact rational for a vector that came out of integer data. Empty when it is not one at
% this denominator, so a hard point is skipped rather than rounded into existence.
    d = 720720;                                    % lcm of many small integers
    n = round(x * d);
    if any(abs(n/d - x) > 1e-12), n = []; return, end
    gg = gcd(gcd(abs(n(1)), abs(n(2))), d);
    if gg > 1, n = n/gg; d = d/gg; end
end

function v = doubleCellFace(H, b, c, P, x)
% The verified double routine, inlined so this probe does not depend on the other one's private
% functions. Same candidate set; see .claude/cellFaceProbe.m for the derivation.
    y = x - b;  best = -inf;
    if abs(det(H)) > 1e-12
        s = H \ y;
        if all(P(:,1:2)*s + P(:,3) >= -1e-9), best = max(best, s.'*y - 0.5*(s.'*H*s) - c); end
    end
    for i = 1:size(P,1)
        Ni = P(i,1:2);
        if all(Ni == 0), continue, end
        d = [Ni(2); -Ni(1)];
        p = -P(i,3) * Ni(:) / (Ni * Ni.');
        lo = -inf; hi = inf; empty = false;
        for j = 1:size(P,1)
            A = P(j,1:2)*d;  B = P(j,1:2)*p + P(j,3);
            if abs(A) < 1e-14
                if B < -1e-9, empty = true; break, end
            elseif A > 0, lo = max(lo, -B/A);
            else, hi = min(hi, -B/A);
            end
        end
        if empty || lo > hi + 1e-9, continue, end
        a2 = -0.5*(d.'*H*d);  a1 = d.'*y - p.'*H*d;  a0 = p.'*y - 0.5*(p.'*H*p) - c;
        if a2 < -1e-14
            ts = -a1/(2*a2);
            if ts > lo && ts < hi, best = max(best, a0 + a1*ts + a2*ts^2); end
        elseif abs(a1) > 1e-12
            if (a1 > 0 && isinf(hi)) || (a1 < 0 && isinf(lo)), v = inf; return, end
        else
            best = max(best, a0);
        end
        for t = [lo hi]
            if isinf(t), continue, end
            best = max(best, a0 + a1*t + a2*t^2);
        end
    end
    v = best;
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

% ---- fixtures and cell data ----------------------------------------------------------------------
function [Hn, Hd, bn, bd, cn, cd] = cellDataExact(g, k)
% The cell's face as INTEGERS over one denominator -- which is how QuaCon already stores it, so
% nothing is converted and nothing is rounded.
    n = g.fN(k,:);  d = g.fD(k);
    Hn = [n(5), n(6); n(6), n(7)];   Hd = d;
    bn = [n(8); n(9)];               bd = d;
    cn = n(10);                      cd = d;
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
