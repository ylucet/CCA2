function bridgeRegionSides()
% ITEM 2, THE LAST PIECE: which SIDE of each boundary conic the bridge region occupies.
%
% WHAT WAS ALREADY SETTLED. The region's two boundary conics are constructed exactly
% (`.claude/bridgeRegionBuild.m`), and there is no envelope arc, so the region really is the
% intersection of TWO SIGN CONDITIONS (`.claude/bridgeRegionCross.m`). An H-form face needs those
% two signs, and nothing so far supplies them.
%
% WHY NOT JUST EVALUATE AT A SAMPLED INTERIOR POINT. Because the sign would then rest on a
% tolerance, and it would be deciding WHICH SIDE OF A FACE THE ENVELOPE LIVES ON -- the one kind of
% question this project has never let a tolerance answer. A point near the boundary would flip it
% silently.
%
% THE ROUTE, and it is the project's own idiom: SCREEN WITH A TOLERANCE, DECIDE EXACTLY.
%
%   1. Propose a candidate interior point in floating point (any classifier will do).
%   2. VERIFY exactly that it is in the region.
%   3. Read the two signs there, exactly, by evaluating integer conics at a rational point.
%
% Step 2 is the part that needed a new idea, and it turns out to be a Sturm count. A point x is in
% the region exactly when the sweeping segment through it exists, i.e. when
%
%       x = M_t s + v_t   for some t in [0,1] and some s on the tie,
%       M_t = (1-t) A1 + t A2 ,  v_t = (1-t) b1 + t b2 ,
%
% and eliminating s gives, after clearing det(M_t),
%
%       Phi(x,t) = q2( adj(M_t)(x - v_t) ) + D_t * q1( adj(M_t)(x - v_t) ) + D_t^2 * q0 ,
%       D_t = det(M_t) ,
%
% with q2/q1/q0 the quadratic, linear and constant parts of the tie conic. `adj(M_t)` is AFFINE in
% t and `D_t` is quadratic, so Phi is a QUARTIC in t with rational coefficients once x is rational.
%
%       x is in the region  <=>  Phi(x, .) has a root in [0,1]
%
% and `ratQ.countRootsIn` decides exactly how many roots a rational polynomial has in an interval,
% by Sturm, with no root finder anywhere. So the floating point only ever PROPOSES; every decision
% is integer.
%
% This also hands over something item 2 needed anyway: an exact membership test for the region.

    fprintf('\n=== isotropic ===\n');
    a = [1; 0.5];  c = -0.25;
    g1 = struct('A', eye(2),     'b', [0;0], 'k', 0);
    g2 = struct('A', 0.5*eye(2), 'b', a,     'k', -c);
    sides(g1, g2, 2*a);

    fprintf('\n=== anisotropic ===\n');
    sides(g1, struct('A', [0.5 0.1; 0.1 0.9], 'b', a, 'k', -c), 2*a);
end

% ------------------------------------------------------------------------------------------------
function sides(g1, g2, ctr)
    q  = diffConic(g1, g2);
    C1 = imageConic(q, g1);   C2 = imageConic(q, g2);
    if isempty(C1) || isempty(C2)
        fprintf('  a boundary degenerates; this routine covers the two-conic case only.\n');  return
    end

    % ---- 1. PROPOSE. Any interior point will do; this walks out from the centre until the
    %         classifier says "bridge", then steps to the middle of that interval.
    [x0, gd] = proposeInterior(g1, g2, ctr);
    if isempty(x0), fprintf('  no interior point proposed.\n'); return, end
    [xn, xd] = rationalise(x0, gd);
    fprintf('  proposed interior point %s (exactly %d/%d, %d/%d)\n', ...
            mat2str(x0, 6), xn(1), xd, xn(2), xd);

    % ---- 2. VERIFY, exactly.
    [inRegion, nRoots] = memberExact(q, g1, g2, xn, xd);
    fprintf('  Sturm on the quartic Phi(x,.): %d root(s) in [0,1]  =>  in the region: %d\n', ...
            nRoots, inRegion);
    if ~inRegion
        fprintf('  the proposal was NOT verified, so no side is claimed. That is the guard working.\n');
        return
    end

    % ---- 3. READ the signs, exactly.
    s1 = signAtRational(C1, xn, xd);
    s2 = signAtRational(C2, xn, xd);
    fprintf('  H-form: region = { x : %d*C1(x) >= 0  and  %d*C2(x) >= 0 }\n', s1, s2);
    fprintf('    C1 = [%s]\n    C2 = [%s]\n', num2str(C1), num2str(C2));

    % ---- cross-check MEMBERSHIP itself, which is what the earlier scaling bug slipped past.
    % The sides check below evaluates C1 and C2 and never touches Phi, so it cannot see a wrong
    % Phi. This one does: the exact Sturm membership test against the numeric classifier, on
    % points chosen on a coarse rational grid so the Sturm chain stays small.
    okM = 0;  badM = 0;
    for a1 = -6:2:10
        for a2 = -6:2:10
            xg = [a1; a2] / 2;
            [xg_n, xg_d] = rationalise(xg, 2);
            exact = memberExact(q, g1, g2, xg_n, xg_d);
            if exact == onTie(g1, g2, xg), okM = okM + 1; else, badM = badM + 1; end
        end
    end
    fprintf('  exact Sturm membership vs the classifier on %d grid points: %d agree, %d disagree\n', ...
            okM + badM, okM, badM);

    % ---- cross-check the H-form against the independent classifier, on fresh points.
    agree = 0;  tested = 0;  wrong = 0;
    rng(11);
    for i = 1:4000
        x = ctr + (rand(2,1)*2 - 1) * 8;
        hform = (s1 * evalConic(C1, x.') >= 0) && (s2 * evalConic(C2, x.') >= 0);
        cls   = onTie(g1, g2, x);
        tested = tested + 1;
        if hform == cls, agree = agree + 1; else, wrong = wrong + 1; end
    end
    fprintf('  H-form vs the independent classifier on %d random points: %d agree, %d disagree\n', ...
            tested, agree, wrong);
end

function [inR, n] = memberExact(q, g1, g2, xn, xd)
% Does the sweeping segment through xn/xd exist? Build Phi(t) over the rationals and Sturm it.
    [A1n, A2n, b1n, b2n, den] = intifyPair(g1, g2);
    % Everything is scaled to a common integer denominator `den`, and x is xn/xd. Phi is formed
    % with polynomial arithmetic in t; a common positive factor never changes a root count, so the
    % scalings are carried rather than tracked.
    dM = [A2n(1,1)-A1n(1,1), A1n(1,1); A2n(1,2)-A1n(1,2), A1n(1,2); ...
          A2n(2,1)-A1n(2,1), A1n(2,1); A2n(2,2)-A1n(2,2), A1n(2,2)];   % each entry as [slope, const]
    m11 = dM(1,:); m12 = dM(2,:); m21 = dM(3,:); m22 = dM(4,:);        % M_t entries, degree 1 in t

    adj11 =  m22;  adj12 = -m12;  adj21 = -m21;  adj22 =  m11;         % adj(M_t), degree 1
    D     = padd(pmul(m11, m22), -pmul(m12, m21));                     % det(M_t), degree 2

    % x - v_t, where v_t = (1-t)b1 + t b2, all over the same den; x is xn/xd, so everything is put
    % over den*xd and the whole of Phi is scaled by a positive power of it.
    w1 = padd([0, den*xn(1)], -pmul([b2n(1)-b1n(1), b1n(1)], [0, xd]));
    w2 = padd([0, den*xn(2)], -pmul([b2n(2)-b1n(2), b1n(2)], [0, xd]));

    u1 = padd(pmul(adj11, w1), pmul(adj12, w2));                       % degree 2
    u2 = padd(pmul(adj21, w1), pmul(adj22, w2));

    % q's parts, with the scalings that make every term the same degree in the carried factors.
    % THE SCALING, WRITTEN OUT, because getting it wrong is invisible. With
    %       U = u_computed = den^2 * xd * adj(M_t) * (x - v_t)   and   Dc = den^2 * det(M_t),
    % the true argument is
    %       s = adj(M_t)(x - v_t) / det(M_t) = U / (xd * Dc) ,
    % so with lambda = xd * Dc, clearing q(s) = 0 by lambda^2 gives
    %       Phi = q2(U) + lambda * q1(U) + lambda^2 * q0 .
    % An earlier version used `den*xd` for lambda instead of `xd`, i.e. it was off by `den` in the
    % linear term and `den^2` in the constant. The ISOTROPIC instance still reported a root and
    % still produced the right sides, because the sides come from evaluating C1 and C2 rather than
    % from Phi -- so the 4000/4000 agreement did NOT validate Phi at all. Only the anisotropic
    % instance, once the BigInteger Sturm path let it run, showed the error by rejecting a witness
    % the classifier accepted.
    lam  = pmul(D, [0, xd]);                          % lambda = xd * Dc, degree 2 in t
    quad = padd(padd(pmul([q(1)], pmul(u1,u1)), pmul([q(2)], pmul(u1,u2))), pmul([q(3)], pmul(u2,u2)));
    linr = pmul(lam, padd(pmul([q(4)], u1), pmul([q(5)], u2)));
    cons = pmul([q(6)], pmul(lam, lam));
    Phi  = padd(padd(quad, linr), cons);

    Phi = Phi / gcdAll(Phi);
    n   = ratQ.countRootsIn(ratQ.chk(Phi, 'Phi'), 0, 1);
    % countRootsIn is half-open (0,1]; t = 0 is the boundary curve C1 and is handled by evaluating.
    if n == 0 && polyval(Phi, 0) == 0, n = 1; end
    inR = n > 0;
end

% ---- small dense-polynomial helpers, highest degree first ---------------------------------------
function c = pmul(a, b), c = conv(a, b); end
function c = padd(a, b)
    n = max(numel(a), numel(b));
    c = [zeros(1, n-numel(a)), a] + [zeros(1, n-numel(b)), b];
end
function g = gcdAll(p)
    g = 0;
    for v = p, g = gcd(g, abs(round(v))); end
    if g == 0, g = 1; end
end

function [A1n, A2n, b1n, b2n, den] = intifyPair(g1, g2)
    [X, den] = intify([g1.A(:); g2.A(:); g1.b; g2.b].');
    A1n = reshape(X(1:4), 2, 2);  A2n = reshape(X(5:8), 2, 2);
    b1n = X(9:10).';              b2n = X(11:12).';
end

function s = signAtRational(C, xn, xd)
% The sign of the integer conic C at the rational point xn/xd, decided in integers.
    v = C(1)*xn(1)^2 + C(2)*xn(1)*xn(2) + C(3)*xn(2)^2 + ...
        C(4)*xn(1)*xd  + C(5)*xn(2)*xd  + C(6)*xd^2;
    s = sign(ratQ.chk(v, 'conic value'));
end

function [n, d] = rationalise(x, d)
    n = round(x * d);
end

function [x0, den] = proposeInterior(g1, g2, ctr)
% A candidate interior point ON A COARSE RATIONAL GRID, and the coarseness is the point.
%
% THE DENOMINATOR IS NOT COSMETIC. Phi's coefficients carry (den*xd)^2, and the Sturm chain's
% pseudo-remainder squares its leading coefficient at every step -- ratQ.sturmChain's own header
% records a quartic with six-digit entries reaching 1e23 in a few steps. A point rationalised at
% 1e6 overflowed at 4.0e30 on the first instance tried. The verification does not need a precise
% point, only a verifiably interior one, so the search runs over a grid of small rationals and
% takes the first the classifier accepts with room on both sides.
    x0 = [];  den = 0;
    for d = [1 2 4 8 16 32]
        best = [];  bestMargin = 0;
        lo = ctr - 8;  hi = ctr + 8;
        for a = ceil(lo(1)*d):floor(hi(1)*d)
            for b = ceil(lo(2)*d):floor(hi(2)*d)
                x = [a; b] / d;
                if ~onTie(g1, g2, x), continue, end
                % Keep the candidate furthest from a classifier flip, so a boundary-adjacent grid
                % point never becomes the witness.
                % Prefer a point well clear of the boundary AND of SMALL MAGNITUDE. Phi's
                % coefficients grow with |x|^2, and the Sturm chain squares its leading coefficient
                % each step, so the witness's size is one of the few handles on the overflow.
                m = margin(g1, g2, x, 1/d);
                if m > bestMargin || (m == bestMargin && ~isempty(best) && norm(x) < norm(best))
                    bestMargin = m;  best = x;
                end
            end
        end
        if ~isempty(best), x0 = best; den = d; return, end
    end
end

function m = margin(g1, g2, x, h)
% How many of the four neighbours at distance h the classifier also calls interior.
    m = 0;
    for dxy = [h 0; -h 0; 0 h; 0 -h].'
        if onTie(g1, g2, x + dxy), m = m + 1; end
    end
end

function tf = onTie(g1, g2, x)
    s1 = g1.A \ (x - g1.b);
    if gval(g1,s1) >= gval(g2,s1) - 1e-12*max(1, s1.'*s1), tf = false; return, end
    s2 = g2.A \ (x - g2.b);
    if gval(g2,s2) >= gval(g1,s2) - 1e-12*max(1, s2.'*s2), tf = false; return, end
    tf = true;
end
function v = gval(g, s), v = 0.5*(s.'*g.A*s) + g.b.'*s + g.k; end

function cimg = imageConic(q, gk)
    A = gk.A;  b = gk.b;
    if abs(det(A)) < 1e-12, cimg = []; return, end
    [An, ad] = intify(A);  [bn, bd] = intify(b);
    L = lcm(ad, bd);  An = An * (L/ad);  bn = bn * (L/bd);
    adjA = [An(2,2), -An(1,2); -An(2,1), An(1,1)];
    dAn  = ratQ.detExact(An);
    [n2, ~] = ratQ.substAffine([0 0 0 0, 2*q(1), q(2), 2*q(3), q(4), q(5), q(6)], 1, ...
                               ratQ.chk(L*adjA), ratQ.chk(-adjA*bn), dAn);
    cimg = ratQ.conic([n2(5), 2*n2(6), n2(7), 2*n2(8), 2*n2(9), 2*n2(10)]);
end

function q = diffConic(g1, g2)
    A = g1.A - g2.A;  b = g1.b - g2.b;  k = g1.k - g2.k;
    [q, ~] = intify([A(1,1)/2, A(1,2), A(2,2)/2, b(1), b(2), k]);
end

function [N, d] = intify(X)
    d = 1;
    for v = X(:).'
        [~, dv] = rat(v, 1e-12);
        d = lcm(d, abs(dv));
    end
    N = round(X * d);
end

function v = evalConic(c, X)
    v = c(1)*X(:,1).^2 + c(2)*X(:,1).*X(:,2) + c(3)*X(:,2).^2 + ...
        c(4)*X(:,1)      + c(5)*X(:,2)        + c(6);
end
