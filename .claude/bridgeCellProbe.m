function bridgeCellProbe()
% CLIPPING. The elimination so far was written per-instance, with the pieces' conjugates hard-coded
% as unclipped paraboloids. The real `biconjQ` case is not that: each piece lives on a POLYTOPE, so
% its conjugate is piecewise -- one quadratic per cell of the normal-cone subdivision, AFFINE on the
% cells where the argmax is pinned at a vertex. A bridge is a pair of CELLS, and there is one
% elimination per pair.
%
% THE GENERALISATION THAT MAKES THAT POSSIBLE is small, because the tangency point IS the gradient
% of the conjugate. Everything the elimination needs about a cell is
%
%       g_k(s) = ( s'A_k s / 2 + c_k's + d_k ) / D_k ,      y_k(s) = grad g_k = (A_k s + c_k)/D_k,
%
% so ONE routine parametrised by (A_k, c_k, d_k, D_k) covers every cell pair and clipping stops
% being a special case:
%
%   * interior cell -> A_k = inv(H_k), the unclipped paraboloid;
%   * EDGE cell     -> A_k of rank 1;
%   * VERTEX cell   -> A_k = 0 and c_k the vertex, so g_k is affine and y_k is CONSTANT.
%
% AND CLIPPING TO A VERTEX MAKES THE PROBLEM EASIER, NOT HARDER. D = det[x - y1, y2 - y1] is
% generically quadratic in s; its quadratic part is det[(-A1)s, (D1*A2 - D2*A1)s], which vanishes
% exactly when A1 and D1*A2 - D2*A1 act in parallel. With A2 = 0 that second matrix is -D2*A1, which
% is ALWAYS parallel to A1 -- so a vertex cell paired with anything is linear in s2 and can be
% SUBSTITUTED rather than resultanted. That is the cheap route, and it is what makes clipping
% tractable at all.
%
% Cleared to integers the three equations are
%   C:  D2*(s'A1 s + 2c1's + 2d1) - D1*(s'A2 s + 2c2's + 2d2) = 0
%   D:  det[ D1*x - A1 s - c1 ,  D1*(A2 s + c2) - D2*(A1 s + c1) ] = 0
%   Z:  s'A1 s + 2c1's + 2d1 - 2*D1*<s,x> + 2*D1*z = 0

    a = [13; -7]/10;  c = -1/2;
    fprintf('%-28s %-9s %-16s %-12s %s\n', 'cell pair', 'time', 'degrees(x,y,z)', 'on face', 'off face');
    fprintf('%s\n', repmat('-', 1, 86));

    % 1. isotropic, both interior. g1 = |s|^2/2; g2 = <s,a> + |s|^2/4 - c, cleared by 20.
    run1('isotropic (int/int)',   eye(2), [0;0], 0, 1,  10*eye(2), 20*a, -20*c, 20, a, c, 1);

    % 2. anisotropic, both interior. inv(H2) = diag(1/2,1/4), cleared by 40.
    run1('anisotropic (int/int)', eye(2), [0;0], 0, 1,  diag([20 10]), 40*a, -40*c, 40, a, c, 0);

    % 3. VERTEX-CLIPPED: piece 2's argmax pinned at v, so g2(s) = <s,v> - q2(v) is affine.
    % q2(v) must be small enough for the two conjugates to MEET at all: the tie is the circle
    % |s - v|^2 = |v|^2 - 2*q2(v), so q2(v) < |v|^2/2 = 2.5 here. The first attempt used 3 and
    % the tie was empty -- which the oracle correctly reported as 'no tie points', and which
    % would have been easy to misread as the elimination failing.
    run1('vertex-clipped (int/vtx)', eye(2), [0;0], 0, 1,  zeros(2), [1;2], 0, 1, a, c, 0);

    % 4. EDGE-CLIPPED: A2 of rank 1 -- the argmax slides along one edge.
    run1('edge-clipped (int/edge)', eye(2), [0;0], 0, 1,  [2 0; 0 0], [1;2], -3, 2, a, c, 0);
end

% ------------------------------------------------------------------------------------------------
function run1(name, A1, c1, d1, D1, A2, c2, d2, D2, ~, ~, ~)
    t0 = tic;
    try
        p = eliminateCellPair(A1, c1, d1, D1, A2, c2, d2, D2);
    catch e
        fprintf('%-28s %-9s %s\n', name, sprintf('%.1f s', toc(t0)), e.identifier);
        return
    end
    el  = toc(t0);
    deg = sprintf('(%d,%d,%d)', size(p,3)-1, size(p,4)-1, size(p,5)-1);

    % ORACLE: points ON the bridge by construction -- s on the tie curve, x the midpoint of the two
    % tangency points, z the face value. It uses none of the elimination.
    rng(20260906);
    worstOn = 0;  worstOff = inf;  n = 0;
    for t = 1:400
        s0 = 5*(rand(2,1) - 0.5);
        [ok, s] = onTie(s0, A1, c1, d1, D1, A2, c2, d2, D2);
        if ~ok, continue, end
        y1 = (A1*s + c1)/D1;   y2 = (A2*s + c2)/D2;
        x  = 0.5*y1 + 0.5*y2;
        z  = s.'*x - (s.'*A1*s/2 + c1.'*s + d1)/D1;
        worstOn  = max(worstOn,  relAt(p, x(1), x(2), z));
        worstOff = min(worstOff, relAt(p, x(1), x(2), z + 0.05));
        n = n + 1;
        if n >= 60, break, end
    end
    if n == 0
        fprintf('%-28s %-9.1f %-16s (no tie points found)\n', name, el, deg);
        return
    end
    fprintf('%-28s %-9s %-16s %-12.2e %.2e  (%d pts)\n', ...
            name, sprintf('%.1f s', el), deg, worstOn, worstOff, n);
end

function [ok, s] = onTie(s0, A1, c1, d1, D1, A2, c2, d2, D2)
% A point on the tie curve g1 = g2, found by SOLVING rather than bisecting.
%
% The first version walked a ray from s0 and bisected on a sign change. That failed to bracket a
% root most of the time -- the tie is a conic and a fixed direction often misses it -- and reported
% "no tie points found" for cell pairs whose elimination was in fact perfectly correct. An oracle
% that silently finds nothing looks exactly like an elimination that silently fails.
%
% Fixing s1 leaves a QUADRATIC in s2, which is solved in closed form.
    ok = false;  s = s0;
    s1 = s0(1);
    al = (A1(2,2)/D1 - A2(2,2)/D2) / 2;
    be = (A1(1,2)/D1 - A2(1,2)/D2)*s1 + (c1(2)/D1 - c2(2)/D2);
    ga = (A1(1,1)*s1^2/2 + c1(1)*s1 + d1)/D1 - (A2(1,1)*s1^2/2 + c2(1)*s1 + d2)/D2;
    if abs(al) < 1e-14
        if abs(be) < 1e-14, return, end
        s2 = -ga/be;
    else
        disc = be^2 - 4*al*ga;
        if disc < 0, return, end
        s2 = (-be + sqrt(disc)) / (2*al);
    end
    if ~isfinite(s2), return, end
    s  = [s1; s2];
    ok = true;
end

function r = relAt(p, x, y, z)
    v = 0;  s = 0;
    for j = 1:size(p,3)
        for k = 1:size(p,4)
            for l = 1:size(p,5)
                cf = p(1,1,j,k,l);
                if cf == 0, continue, end
                t = cf * x^(j-1) * y^(k-1) * z^(l-1);
                v = v + t;  s = s + abs(t);
            end
        end
    end
    r = abs(v)/max(s,1);
end

% ------------------------------------------------------------------------------------------------
function p = eliminateCellPair(A1, c1, d1, D1, A2, c2, d2, D2)
% Variables (s1, s2, x, y, z). D is SOLVED for s2 when it is linear there -- see the header for when
% that holds and why a vertex cell always does -- and only then is one resultant taken in s1.
    nv = 5;
    m  = @(e) ratQ.monoN(e, 1, nv);
    S1 = m([1 0 0 0 0]);  S2 = m([0 1 0 0 0]);
    X  = m([0 0 1 0 0]);  Y  = m([0 0 0 1 0]);  Zz = m([0 0 0 0 1]);
    one = m([0 0 0 0 0]);
    ad = @(varargin) foldAdd(varargin);
    mu = @(A,B) ratQ.polyMulN(A,B);
    qd = @(A) ad(mu(S1,S1)*A(1,1), mu(S1,S2)*(2*A(1,2)), mu(S2,S2)*A(2,2));

    C = ratQ.trimN(ratQ.contentFree(ad(qd(A1)*D2, S1*(2*c1(1)*D2), S2*(2*c1(2)*D2), one*(2*d1*D2), ...
                            qd(A2)*(-D1), S1*(-2*c2(1)*D1), S2*(-2*c2(2)*D1), one*(-2*d2*D1))));

    u1 = ad(X*D1, S1*(-A1(1,1)), S2*(-A1(1,2)), one*(-c1(1)));
    u2 = ad(Y*D1, S1*(-A1(2,1)), S2*(-A1(2,2)), one*(-c1(2)));
    w1 = ad(S1*(D1*A2(1,1)-D2*A1(1,1)), S2*(D1*A2(1,2)-D2*A1(1,2)), one*(D1*c2(1)-D2*c1(1)));
    w2 = ad(S1*(D1*A2(2,1)-D2*A1(2,1)), S2*(D1*A2(2,2)-D2*A1(2,2)), one*(D1*c2(2)-D2*c1(2)));
    % trimN BEFORE asking the degree: contentFree does not trim, so trailing zero slices would
    % make a linear D look quadratic -- which is exactly what it did, on all four cell pairs.
    Dp = ratQ.trimN(ratQ.contentFree(ad(mu(u1,w2), mu(u2,w1)*(-1))));

    Zp = ratQ.trimN(ratQ.contentFree(ad(qd(A1), S1*(2*c1(1)), S2*(2*c1(2)), one*(2*d1), ...
                             mu(S1,X)*(-2*D1), mu(S2,Y)*(-2*D1), Zz*(2*D1))));

    if size(Dp,2) > 2
        error('ratQ:notLinearInS2', ...
              ['D is genuinely quadratic in s2, so it cannot be substituted. That needs two ' ...
               'resultants and a 16x16 Sylvester, which cofactor expansion cannot do.']);
    end

    % D = As*s2 + Bs. The two share a scale, so any common factor comes off TOGETHER -- reducing
    % them separately would move the root -Bs/As, which is a mistake this project has already paid
    % for once (DECISIONS.md 2026-09-06).
    As = sliceS2(Dp, 2, nv);   Bs = sliceS2(Dp, 1, nv);
    g  = gcd(ratQ.contentOf(As), ratQ.contentOf(Bs));
    if g > 1, As = As/g;  Bs = Bs/g; end
    As2 = mu(As,As);  Bs2 = mu(Bs,Bs);  AB = mu(As,Bs);

    Ct = subst(C,  As, Bs, As2, Bs2, AB, nv);
    Zt = subst(Zp, As, Bs, As2, Bs2, AB, nv);
    Ct = ratQ.trimN(ratQ.contentFree(Ct));   Zt = ratQ.trimN(ratQ.contentFree(Zt));
    if getenv('CELLPROBE_VERBOSE')
        fprintf('      C %s  D %s  Z %s | As %s Bs %s | Ct %s Zt %s\n', ...
                sz(C), sz(Dp), sz(Zp), sz(As), sz(Bs), sz(Ct), sz(Zt));
    end

    p = ratQ.trimN(ratQ.fromBig(ratQ.bigContentFree(ratQ.resultantVarBig(Ct, Zt, 1, nv))));
end

function T = subst(P, As, Bs, As2, Bs2, AB, nv)
% substitute s2 = -Bs/As into P and clear As^2. P has degree <= 2 in s2.
    P0 = sliceS2(P, 1, nv);
    P1 = sliceS2(P, 2, nv);
    P2 = sliceS2(P, 3, nv);
    T = ratQ.polyAddN(ratQ.polyMulN(P0, As2), ...
        ratQ.polyAddN(ratQ.polyMulN(P1, AB)*(-1), ratQ.polyMulN(P2, Bs2)));
end

function S = sliceS2(P, k, nv)
% the coefficient of s2^(k-1), as an array with a singleton in the s2 slot
    sz = [size(P), ones(1, nv - ndims(P))];
    if k > sz(2), S = 0; return, end
    idx = repmat({':'}, 1, nv);
    idx{2} = k;
    S = P(idx{:});
end

function A = foldAdd(items)
    A = 0;
    for i = 1:numel(items), A = ratQ.polyAddN(A, items{i}); end
end

function s = sz(A)
    d = size(A);
    s = sprintf('%dx%dx%dx%dx%d', [d, ones(1,5-numel(d))]);
end
