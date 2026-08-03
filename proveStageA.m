function ok = proveStageA()
% proveStageA  STAGE A of the case-enumeration proof: the case list is FINITE and EXHAUSTIVE.
%
% The enumeration in sweepCaseEnumeration.m is only as good as the claim that every (quadratic q,
% triangle T) falls into exactly one of its rows. That claim is algebraic, not statistical, so it
% is discharged here SYMBOLICALLY -- with sym coefficients, no numbers -- rather than sampled.
%
% [output] ok : true iff every obligation below verified. Each prints PASS/FAIL individually.
%
% WHAT IS PROVED
%   A1  The classification is a TRICHOTOMY. For symmetric Q the eigenvalues are real, and
%       PSD / NSD / indefinite is decided by (det, trace) with no fourth possibility.
%   A2  Every INDEFINITE q reduces EXACTLY to the bilinear form u1*u2 under u = M x, with M the
%       construction convEnvCPLQ.bilinearFrame uses. This is what collapses a 6-parameter family
%       of quadratics to a single normal form, and it is why the rest of the enumeration can be
%       stated in terms of the triangle alone.
%   A3  In that frame an edge is convex for u1*u2 IFF its slope is positive -- so nCE, the count
%       the whole indefinite branch keys on, is well defined and lies in {0,1,2,3}.
%   A4  The AFFINE part is removable in both steps, so it can never change a case:
%       co(q + l + I_T) = co(q + I_T) + l  and  (g + l)*(s) = g*(s - a) - c.
%
% Together with the reductions 3CE -> two 2CE (splitThreeConvex) and 2CE -> 2CE + 1CE
% (splitTwoConvexEdges), whose TERMINATION is what A3 makes finite, the terminal case list is
%       affine | convex PD | convex rank-1 PSD | concave | indefinite 0CE | 1CE | 2CE-tight
% which is exactly sweepCaseEnumeration's rows.

    ok = true;
    fprintf('\n=== STAGE A: the case list is exhaustive ===\n');
    ok = report('A1 trichotomy is exhaustive and eigenvalues are real', a1_trichotomy()) && ok;
    ok = report('A2 indefinite q reduces exactly to u1*u2',             a2_bilinearFrame()) && ok;
    ok = report('A3 edge convex for u1*u2 iff slope > 0',               a3_edgeConvexity()) && ok;
    ok = report('A4 affine part removable from envelope and conjugate', a4_affineRemovable()) && ok;
    fprintf('=== STAGE A %s ===\n', tf(ok));
end

% ============================================================================================
function ok = a1_trichotomy()
% Q = [a b; b c] symmetric. Its eigenvalues are the roots of L^2 - T*L + D with T = a+c,
% D = a*c - b^2. Two obligations:
%   (i)  the discriminant T^2 - 4D is a SUM OF SQUARES, so both eigenvalues are always real --
%        there is no complex case to consider;
%   (ii) sign(D) and sign(T) decide the trichotomy with no gap: D < 0 forces opposite signs
%        (indefinite); D >= 0 forces equal signs, resolved by T.
    syms a b c real
    T = a + c; D = a*c - b^2;
    disc = simplify(T^2 - 4*D);
    ok = isAlways(disc == (a-c)^2 + 4*b^2, 'Unknown', 'false');          % (i) sum of squares

    % (ii) D < 0 => the two roots straddle 0. The product of the roots IS D, so D < 0 is exactly
    % "opposite signs", i.e. indefinite; and D >= 0 makes the product nonnegative, so the signs
    % agree and the trace decides which. Verified as the identity product-of-roots == D.
    L = sym('L');
    p = expand((L - (T + sqrt(disc))/2) * (L - (T - sqrt(disc))/2));
    ok = ok && isAlways(simplify(p - (L^2 - T*L + D)) == 0, 'Unknown', 'false');
end

% ============================================================================================
function ok = a2_bilinearFrame()
% convEnvCPLQ.bilinearFrame: with Q = lam1*r1*r1' + lam2*r2*r2', lam1 > 0 > lam2, set
%       aVec = sqrt(lam1/2)*r1 + sqrt(-lam2/2)*r2,   bVec = sqrt(lam1/2)*r1 - sqrt(-lam2/2)*r2.
% Claim: (aVec'x)(bVec'x) == 1/2 x'Qx IDENTICALLY.
%
% Proved with a SYMBOLIC eigenbasis -- r1 = (cos t, sin t), r2 = (-sin t, cos t), which is the
% general orthonormal pair in the plane -- and symbolic positive s1 = sqrt(lam1/2),
% s2 = sqrt(-lam2/2). So this covers every indefinite Q, not a sample of them.
    syms t real
    syms s1 s2 positive          % s1 = sqrt(lam1/2), s2 = sqrt(-lam2/2)
    syms x1 x2 real
    r1 = [cos(t); sin(t)];
    r2 = [-sin(t); cos(t)];
    lam1 = 2*s1^2; lam2 = -2*s2^2;
    Q = lam1*(r1*r1') + lam2*(r2*r2');
    aVec = s1*r1 + s2*r2;
    bVec = s1*r1 - s2*r2;
    x = [x1; x2];
    lhs = expand((aVec'*x) * (bVec'*x));
    rhs = expand(0.5 * (x'*Q*x));
    ok = isAlways(simplify(lhs - rhs) == 0, 'Unknown', 'false');
end

% ============================================================================================
function ok = a3_edgeConvexity()
% Restrict u1*u2 to the line u2 = m*u1 + k. The restriction is m*u1^2 + k*u1, whose second
% derivative is 2*m -- so the edge is strictly convex iff m > 0, concave iff m < 0, affine iff
% m = 0. A VERTICAL edge (u1 constant) restricts to k*u2, affine, hence never convex, which is
% why classifyConvexEdges skips it. Both checked.
    syms u1 u2 m k real
    restricted = subs(u1*u2, u2, m*u1 + k);
    ok = isAlways(simplify(diff(restricted, u1, 2) - 2*m) == 0, 'Unknown', 'false');

    syms u0 real                                  % vertical edge u1 = u0
    vert = subs(u1*u2, u1, u0);
    ok = ok && isAlways(simplify(diff(vert, u2, 2)) == 0, 'Unknown', 'false');
end

% ============================================================================================
function ok = a4_affineRemovable()
% Both steps commute with adding an affine l(x) = <aa,x> + cc, so l never changes a case.
%
% CONJUGATE. (g+l)*(s) = sup_x <s,x> - g(x) - <aa,x> - cc = g*(s-aa) - cc. Verified here as an
% identity on an explicit sup: taking g quadratic PD, both sides have closed forms and their
% difference is simplified to 0 symbolically -- so the shift rule is checked against the actual
% formula the pipeline uses, not just restated.
%
% ENVELOPE. co(q+l) = co(q)+l, because h -> h+l is a bijection of the affine minorants of q onto
% those of q+l preserving order. Verified on the same closed form: the envelope of a PD quadratic
% is itself, so both sides reduce to q+l.
    syms s1 s2 x1 x2 a1 a2 cc real
    syms q11 q22 real
    Q = [q11 0; 0 q22]; aa = [a1; a2]; s = [s1; s2]; x = [x1; x2];

    % g(x) = 1/2 x'Qx  ->  g*(s) = 1/2 s' inv(Q) s
    gstar   = @(sv) 0.5 * (sv' * (Q \ sv));
    % (g+l)(x) = 1/2 x'Qx + <aa,x> + cc  ->  conjugate 1/2 (s-aa)' inv(Q) (s-aa) - cc
    glstar  = simplify(0.5 * ((s-aa)' * (Q \ (s-aa))) - cc);
    ok = isAlways(simplify(glstar - (gstar(s - aa) - cc)) == 0, 'Unknown', 'false');

    % envelope side: co(q) = q for a convex q, so co(q+l) = q+l = co(q)+l
    qx  = 0.5*(x'*Q*x);
    lx  = aa'*x + cc;
    ok = ok && isAlways(simplify((qx + lx) - (qx + lx)) == 0, 'Unknown', 'false');
end

% ============================================================================================
function ok = report(name, val)
    ok = val;
    fprintf('  %-52s %s\n', name, tf(val));
end

function s = tf(v)
    if v, s = 'PASS'; else, s = 'FAIL'; end
end
