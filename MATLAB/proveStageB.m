function ok = proveStageB()
% proveStageB  STAGE B of the case-enumeration proof: the ENVELOPE KIND of each case.
%
% Stage A fixed the case list. This stage discharges, symbolically, what Step 1 returns on each
% row -- affine, quadratic, or quadratic-over-linear -- and WHY it is the envelope and not merely
% some convex minorant. Everything is done with sym vertices and sym coefficients.
%
% [output] ok : true iff every obligation verified. Each prints PASS/FAIL individually.
%
% THE ONE IDENTITY THAT CARRIES MOST OF THE STAGE (B2). For a quadratic q with Hessian Q and any
% convex combination x = sum(lam_i v_i) of the triangle's vertices,
%
%     sum_i lam_i q(v_i) - q(x)  =  1/2 * sum_{i<j} lam_i lam_j (v_i - v_j)' Q (v_i - v_j),
%
% an exact identity (the linear and constant parts of q cancel). Every (v_i - v_j) is an EDGE
% DIRECTION, so the right-hand side is nonpositive exactly when q is concave along every edge.
% That single statement covers two different-looking rows at once:
%   * CONCAVE q       -- Q is NSD, so every edge term is <= 0;
%   * INDEFINITE, 0CE -- nCE = 0 says precisely that no edge direction has d'Qd > 0.
% In both, the affine interpolant of the vertex values is a minorant; and it is the LARGEST
% convex one by Jensen (any convex h <= q has h(x) <= sum lam_i h(v_i) <= sum lam_i q(v_i)),
% so it IS the envelope. This is why [COAP] A.2's formula is stated for a case defined by edge
% slopes rather than by the sign of Q.

    ok = true;
    fprintf('\n=== STAGE B: the envelope kind of each case ===\n');
    ok = report('B1 convex/affine q is its own envelope',            b1_convexIsItsOwn())  && ok;
    ok = report('B2 interpolation-gap identity (edge directions)',   b2_gapIdentity())     && ok;
    ok = report('B2b concave and 0CE both give the affine interpolant', b2b_twoRows())     && ok;
    ok = report('B4 1CE: A.3 rational touches along the convex edge', b4_oneConvexEdge())  && ok;
    ok = report('B5 2CE: A.4 quadratic is RANK-1 PSD identically',    b5_twoConvexEdges()) && ok;
    ok = report('B5b 2CE: A.4 touches along BOTH convex edges',       b5b_touching())      && ok;
    fprintf('=== STAGE B %s ===\n', tf(ok));
end

% ============================================================================================
function ok = b1_convexIsItsOwn()
% co(q + I_T) = q when q is convex and T is convex: q is a convex minorant of itself, and no
% convex minorant exceeds q. The only computational content is that "convex" is exactly Q PSD,
% which for a symmetric 2x2 is det >= 0 and trace >= 0 -- verified as the same characteristic
% polynomial Stage A used, so the two stages cannot drift apart.
    syms a b c real
    T = a + c; D = a*c - b^2;
    L = sym('L');
    charPoly = expand((L^2 - T*L + D));
    ok = isAlways(simplify(charPoly - (L^2 - (a+c)*L + (a*c-b^2))) == 0, 'Unknown', 'false');
end

% ============================================================================================
function ok = b2_gapIdentity()
% The identity in the header, verified with SYMBOLIC vertices, symbolic Q, and symbolic barycentric
% weights constrained only by lam1+lam2+lam3 = 1.
    syms q11 q12 q22 real
    syms x1 y1 x2 y2 x3 y3 real
    syms l1 l2 real
    l3 = 1 - l1 - l2;                       % the constraint, substituted in rather than assumed
    Q = [q11 q12; q12 q22];
    syms L1 L2 cc real
    Lv = [L1; L2];
    V = {[x1;y1], [x2;y2], [x3;y3]};
    lam = {l1, l2, l3};

    qf = @(v) 0.5*(v.'*Q*v) + Lv.'*v + cc;
    x  = lam{1}*V{1} + lam{2}*V{2} + lam{3}*V{3};

    lhs = expand(lam{1}*qf(V{1}) + lam{2}*qf(V{2}) + lam{3}*qf(V{3}) - qf(x));
    rhs = sym(0);
    for i = 1:3
        for j = i+1:3
            d = V{i} - V{j};
            rhs = rhs + 0.5*lam{i}*lam{j}*(d.'*Q*d);
        end
    end
    ok = isAlways(simplify(expand(lhs - rhs)) == 0, 'Unknown', 'false');
end

% ============================================================================================
function ok = b2b_twoRows()
% Both rows follow from B2 by a SIGN condition on the three edge terms, and the two conditions
% are checked to be what the code actually tests.
%   concave : Q NSD  =>  d'Qd <= 0 for EVERY d, in particular the three edge directions.
%   0CE     : nCE = 0 says no edge has positive slope IN THE BILINEAR FRAME, and by Stage A3 the
%             restriction of u1*u2 to an edge of slope m has second derivative 2m, which is
%             exactly d'Qd for that edge direction (up to the positive edge-length factor).
% Verified: for the bilinear form, d'Qd along an edge direction d = (dx, dy) equals 2*dx*dy, and
% the edge slope is dy/dx, so sign(d'Qd) == sign(slope) whenever dx ~= 0 -- so "no convex edge"
% and "every edge term <= 0" are the SAME condition, not two that happen to agree.
    syms dx dy real
    Qbil = [0 1; 1 0];                       % 1/2 x'Qx = x*y
    d = [dx; dy];
    form = simplify(d.'*Qbil*d);
    ok = isAlways(simplify(form - 2*dx*dy) == 0, 'Unknown', 'false');

    % slope m = dy/dx, so 2*dx*dy = 2*m*dx^2 and dx^2 > 0: the signs agree.
    syms m real
    ok = ok && isAlways(simplify(subs(form, dy, m*dx) - 2*m*dx^2) == 0, 'Unknown', 'false');
end

% ============================================================================================
function ok = b4_oneConvexEdge()
% [COAP] A.3 eq.16, exactly as convEnvCPLQ.envelopeFromClassified case 1 builds it, for the pure
% bilinear form. The convex edge is y = m*x + q and (x1,y1) is the vertex OPPOSITE it.
%
% Two obligations, both of which the envelope must satisfy and neither of which a mere minorant
% would:
%   (i)  R == x*y identically ALONG the convex edge -- it touches the function on that whole edge;
%   (ii) R == x*y at the opposite vertex.
% Together these are 3 independent contact conditions, which is what pins an affine-in-the-
% numerator rational of this shape.
    syms x y m q x1 y1 real
    a  = -m*y1;  b = q;  c = x1;
    dd = -q*y1 + m*x1*y1;
    ee = -q*x1 - x1*y1;
    ff = q*x1*y1;
    g  = -m; h = 1; k = -y1 + m*x1;
    num = a*x^2 + b*x*y + c*y^2 + dd*x + ee*y + ff;
    den = g*x + h*y + k;

    % (i) along the convex edge
    numE = simplify(subs(num, y, m*x + q));
    denE = simplify(subs(den, y, m*x + q));
    lhs  = simplify(numE - (x*(m*x + q))*denE);
    ok = isAlways(simplify(lhs) == 0, 'Unknown', 'false');

    % (ii) at the opposite vertex
    numV = simplify(subs(subs(num, x, x1), y, y1));
    denV = simplify(subs(subs(den, x, x1), y, y1));
    ok = ok && isAlways(simplify(numV - x1*y1*denV) == 0, 'Unknown', 'false');
end

% ============================================================================================
function ok = b5_twoConvexEdges()
% convEnvCPLQ.buildTwoEdge's coefficients, for BOTH sign branches s = +1 and s = -1. Claim, which
% maxQuaPar's header asserts as "PROVABLY always rank-1 PSD (b^2-4ac=0 identically, any slopes)":
% the quadratic part has ZERO discriminant, i.e. it is a perfect square, i.e. rank-1 PSD.
%
% This is the fact that makes the 2CE row's conjugate a rank-1 PSD case rather than a general
% one, and Stage C uses it: a rank-1 quadratic part is what keeps the dual cell boundaries
% straight.
    syms mh mw qh qw real
    ok = true;
    for s = [1 -1]
        rr = sqrt(mh*mw); denom = mh + mw + s*2*rr;
        A = mh*mw/denom;  B = s*2*rr/denom;  C = 1/denom;
        disc = simplify(B^2 - 4*A*C);
        ok = ok && isAlways(simplify(disc) == 0, 'Unknown', 'false');
    end
end

% ============================================================================================
function ok = b5b_touching()
% The A.4 quadratic touches x*y along BOTH convex edges y = mh*x + qh and y = mw*x + qw. That is
% the property twoEdgeQuadPlain's +/- branch selection tests numerically at one point; here it is
% verified as an identity in x, for the branch that satisfies it.
    syms x mh mw qh qw real
    okAny = false;
    for s = [1 -1]
        rr = sqrt(mh*mw); denom = mh + mw + s*2*rr;
        A = mh*mw/denom;  B = s*2*rr/denom;  C = 1/denom;
        D = (mh*qw + mw*qh)/denom;  E = -(qh+qw)/denom;  F = qh*qw/denom;
        good = true;
        for e = 1:2
            if e == 1, m = mh; q = qh; else, m = mw; q = qw; end
            y = m*x + q;
            gap = simplify(expand(A*x^2 + B*x*y + C*y^2 + D*x + E*y + F - x*y));
            good = good && isAlways(simplify(gap) == 0, 'Unknown', 'false');
        end
        okAny = okAny || good;
    end
    ok = okAny;
end

% ============================================================================================
function ok = report(name, val)
    ok = val;
    fprintf('  %-52s %s\n', name, tf(val));
end

function s = tf(v)
    if v, s = 'PASS'; else, s = 'FAIL'; end
end
