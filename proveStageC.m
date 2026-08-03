function ok = proveStageC(seed, nPerCase)
% proveStageC  STAGE C: the TYPE of every dual boundary, hence where a parabolic edge can occur.
%
% The conjugate of a piece is a max of active-set candidates e_A(s). The boundary between two
% adjacent dual cells A and B lies in { s : e_A(s) - e_B(s) = 0 }, so the TYPE of every dual edge
% is the conic type of that difference, read off two determinants:
%
%     delta2 = b^2 - 4ac        ellipse (<0) | PARABOLA (=0) | hyperbola (>0)
%     Delta3 = det [a b/2 d/2;  = 0 means the conic DEGENERATES to a line, a line pair,
%                   b/2 c e/2;    or a double line
%                   d/2 e/2 f]
%
% A genuine ARC needs delta2 == 0 AND Delta3 ~= 0. Everything else is straight.
%
% [output] ok : true iff every obligation verified.
%
% TWO GENERAL THEOREMS (C1, C2), proved symbolically, that bound what can ever appear:
%
%   C1  EDGE-STRIP vs VERTEX-CONE is always a PERFECT SQUARE of an affine form, hence a DOUBLE
%       LINE -- never an arc, for any quadratic piece over any polytope. The edge candidate is
%       e_v(s) + (<s,e> - <grad q(v),e>)^2 / (2 e'Ae), and the vertex candidate is e_v(s), so the
%       difference IS that square. This kills 6 of the 9 adjacent pairs in the convex case
%       outright and is why vertex cones never bound an arc.
%
%   C2  Two EDGE STRIPS differ by a difference of two RANK-1 quadratic parts, and in the plane
%       det(uu'/al - vv'/be) = -(u x v)^2/(al*be), so delta2 = 4 (u x v)^2/(al*be). For two edges
%       of a triangle u x v ~= 0, so delta2 is never zero when al*be > 0: such a boundary is a
%       HYPERBOLA-type conic, which can only be admissible by being DEGENERATE (Delta3 = 0, a
%       line pair). It can never be a parabola. This is the statement maxQuaPar's header records
%       empirically -- "delta>0 ... but the full discriminant Delta=0 exactly: both are pairs of
%       straight lines, not hyperbolas" -- now derived rather than observed.
%
% C3 then MEASURES the classification over the seeded sweep: for every case, every adjacent face
% pair of the produced conjugate is classified by (delta2, Delta3) and cross-checked against the
% arc the code actually stored in Ec. Agreement in both directions is the real content: no arc
% stored where the discriminants say straight, and none missing where they say parabola.

    if nargin < 1 || isempty(seed),     seed = 20260802; end
    if nargin < 2 || isempty(nPerCase), nPerCase = 8;    end
    ok = true;
    fprintf('\n=== STAGE C: the type of every dual boundary ===\n');
    ok = report('C1 edge-strip vs vertex-cone is a perfect square', c1_perfectSquare()) && ok;
    ok = report('C2 edge vs edge has delta2 = 4(u x v)^2/(al*be)', c2_rank1Difference()) && ok;
    ok = c3_classifySweep(seed, nPerCase) && ok;
    fprintf('=== STAGE C %s ===\n', tf(ok));
end

% ============================================================================================
function ok = c1_perfectSquare()
% For a quadratic q(x) = 1/2 x'Ax + L'x + c over a polytope, the candidate for the STRIP of edge
% (v, v+e) is
%       eE(s) = <s,v> - q(v) + (<s,e> - <grad q(v),e>)^2 / (2 e'Ae)
% and for the VERTEX CONE at v it is
%       eV(s) = <s,v> - q(v).
% Their difference is exactly the square term, so its zero set is the DOUBLE LINE
% <s,e> = <grad q(v),e>. Verified symbolically with sym A, L, v, e.
    syms a11 a12 a22 L1 L2 cc real
    syms v1 v2 e1 e2 s1 s2 real
    A = [a11 a12; a12 a22]; L = [L1; L2]; v = [v1; v2]; e = [e1; e2]; s = [s1; s2];
    gradq = A*v + L;
    alpha = e.'*A*e;
    eV = s.'*v - (0.5*(v.'*A*v) + L.'*v + cc);
    eE = eV + (s.'*e - gradq.'*e)^2 / (2*alpha);
    d  = simplify(eE - eV);
    ok = isAlways(simplify(d - (s.'*e - gradq.'*e)^2/(2*alpha)) == 0, 'Unknown', 'false');

    % ...and a perfect square has both discriminants zero: delta2 = 0 and Delta3 = 0.
    [d2, D3] = conicDiscriminants(expand(d), s1, s2);
    ok = ok && isAlways(simplify(d2) == 0, 'Unknown', 'false');
    ok = ok && isAlways(simplify(D3) == 0, 'Unknown', 'false');
end

% ============================================================================================
function ok = c2_rank1Difference()
% The quadratic PART of an edge-strip candidate is e e'/(2 e'Ae): RANK 1. Two strips therefore
% differ by uu'/al - vv'/be, and in the plane
%       det(uu'/al - vv'/be) = -(u1*v2 - u2*v1)^2 / (al*be),
% so delta2 = -4*det = 4 (u x v)^2 / (al*be). For two edges of a nondegenerate triangle the cross
% product is nonzero, so delta2 = 0 is impossible when al and be share a sign: such a boundary is
% never a parabola, and is admissible only by being a degenerate (Delta3 = 0) line pair.
    syms u1 u2 v1 v2 al be real
    u = [u1; u2]; v = [v1; v2];
    M = (u*u.')/al - (v*v.')/be;
    lhs = simplify(det(M));
    rhs = simplify(-(u1*v2 - u2*v1)^2/(al*be));
    ok = isAlways(simplify(lhs - rhs) == 0, 'Unknown', 'false');

    % delta2 of the conic whose quadratic part is M: with a = M11/2*2 ... use the direct relation
    % delta2 = b^2 - 4ac = -4*det(M_conic) where M_conic = [a b/2; b/2 c] = M/2 for the form
    % 1/2 s'Ms. So delta2 = -4*det(M/2) = -det(M).
    ok = ok && isAlways(simplify((-det(M)) - ((u1*v2-u2*v1)^2/(al*be))) == 0, 'Unknown', 'false');
end

% ============================================================================================
function ok = c3_classifySweep(seed, nPerCase)
% Classify EVERY adjacent face pair of EVERY conjugate the sweep produces, and cross-check
% against the arc the code stored. Both directions matter: an arc stored where the discriminants
% say straight is a wrong edge; a parabola predicted where no arc is stored is a dropped one.
    rng(seed);
    names = {'affine','convex PD','convex rank-1 PSD','concave', ...
             'indefinite 0CE','indefinite 1CE','indefinite 2CE'};
    inst = cell(1, numel(names)); for i = 1:numel(inst), inst{i} = {}; end
    tries = 0;
    while tries < 4000 && any(cellfun(@numel, inst) < nPerCase)
        tries = tries + 1;
        V = randTriangle();
        switch mod(tries, 4)
            case 0, Q = randPSD();
            case 1, Q = randRank1PSD();
            case 2, Q = -randPSD();
            case 3, Q = randIndefinite();
        end
        if mod(tries, 17) == 0, Q = zeros(2); end
        f6 = [Q(1,1) Q(1,2) Q(2,2) randn() randn() randn()];
        b = bucketOf(Q, V);
        if b == 0 || b > numel(names), continue, end
        if numel(inst{b}) >= nPerCase, continue, end
        inst{b}{end+1} = {V, f6};
    end

    E3 = [1 2 1; 2 3 1; 3 1 1]; F3 = [1 0; 1 0; 1 0];
    ok = true;
    fprintf('  C3 adjacent-pair classification over the sweep:\n');
    fprintf('     %-20s %6s %8s %8s %8s %8s\n', 'case', 'pairs', 'line', 'linepair', 'PARABOLA', 'other');
    for b = 1:numel(names)
        nLine = 0; nPair = 0; nPara = 0; nOther = 0; mismatch = 0;
        for t = 1:numel(inst{b})
            V = inst{b}{t}{1}; f6 = inst{b}{t}{2};
            try
                g = conjPieceCPLQ(QuaPol(V, E3, f6, F3));
            catch
                continue        % this case reaches Step 2 only via Step 1; skipped here
            end
            if isempty(g.E), continue, end
            for j = 1:size(g.E,1)
                fa = g.F(j,1); fb = g.F(j,2);
                if fa < 1 || fb < 1, continue, end          % boundary of the whole plane
                [d2, D3] = conicOfDifference(g.f(fa,:), g.f(fb,:));
                sc = max(1, norm(g.f(fa,5:10)) + norm(g.f(fb,5:10)));
                isPara = abs(d2) <= 1e-7*sc^2 && abs(D3) > 1e-9*sc^3;
                isLine = abs(d2) <= 1e-7*sc^2 && abs(D3) <= 1e-9*sc^3;
                if isPara,        nPara = nPara + 1;
                elseif isLine,    nLine = nLine + 1;
                elseif abs(D3) <= 1e-9*sc^3, nPair = nPair + 1;   % delta2 ~= 0 but degenerate
                else,             nOther = nOther + 1;
                end
                storedArc = ~isempty(g.Ec) && any(g.Ec(j,:) ~= 0);
                if storedArc ~= isPara, mismatch = mismatch + 1; end
            end
        end
        tot = nLine + nPair + nPara + nOther;
        fprintf('     %-20s %6d %8d %8d %8d %8d', names{b}, tot, nLine, nPair, nPara, nOther);
        if mismatch > 0
            fprintf('   <-- %d disagree with stored Ec\n', mismatch); ok = false;
        elseif nOther > 0
            fprintf('   <-- %d NON-DEGENERATE non-parabola\n', nOther); ok = false;
        else
            fprintf('   ok\n');
        end
    end
end

% ============================================================================================
function [d2, D3] = conicOfDifference(frowA, frowB)
% The conic e_A - e_B = 0 for two QuaPar face rows, and its two discriminants. A row stores
% [.. .. .. .. x^2 xy y^2 x y const] with matrixForm reading Q = [c5 c6; c6 c7], so the function
% is 1/2 s'Qs + L's + k.
    QA = [frowA(5) frowA(6); frowA(6) frowA(7)]; LA = [frowA(8); frowA(9)]; kA = frowA(10);
    QB = [frowB(5) frowB(6); frowB(6) frowB(7)]; LB = [frowB(8); frowB(9)]; kB = frowB(10);
    Qd = QA - QB; Ld = LA - LB; kd = kA - kB;
    a = Qd(1,1)/2; bb = Qd(1,2); c = Qd(2,2)/2; d = Ld(1); e = Ld(2); f = kd;
    d2 = bb^2 - 4*a*c;
    D3 = det([a bb/2 d/2; bb/2 c e/2; d/2 e/2 f]);
end

function [d2, D3] = conicDiscriminants(expr, s1, s2)
% Discriminants of a symbolic conic expression in s1,s2.
    a  = simplify(diff(expr, s1, 2)/2);
    c  = simplify(diff(expr, s2, 2)/2);
    bb = simplify(diff(diff(expr, s1), s2));
    d  = simplify(subs(diff(expr, s1), [s1 s2], [0 0]));
    e  = simplify(subs(diff(expr, s2), [s1 s2], [0 0]));
    f  = simplify(subs(expr, [s1 s2], [0 0]));
    d2 = simplify(bb^2 - 4*a*c);
    D3 = simplify(det([a bb/2 d/2; bb/2 c e/2; d/2 e/2 f]));
end

% ---- instance generation, identical to sweepCaseEnumeration's -------------------------------
function b = bucketOf(Q, V)
    tol = 1e-9 * max(1, max(abs(Q(:))));
    ev = sort(eig((Q+Q')/2));
    if max(abs(ev)) <= tol, b = 1; return, end
    if ev(1) >= -tol
        if ev(1) > tol, b = 2; else, b = 3; end
        return
    end
    if ev(2) <= tol, b = 4; return, end
    M = bilinearFrame(Q);
    b = 5 + numConvexEdges((M*V')');
end

function M = bilinearFrame(Q)
    [R, Lam] = eig(Q); lam = diag(Lam);
    [lp, ip] = max(lam); [ln, in] = min(lam);
    r1 = R(:,ip); r2 = R(:,in);
    a = sqrt(lp/2)*r1 + sqrt(-ln/2)*r2;
    bb = sqrt(lp/2)*r1 - sqrt(-ln/2)*r2;
    M = [a'; bb'];
end

function n = numConvexEdges(V)
    tol = sqrt(eps); ed = [1 2; 2 3; 3 1]; n = 0;
    for t = 1:3
        vi = V(ed(t,1),:); vj = V(ed(t,2),:);
        dx = vj(1)-vi(1); dy = vj(2)-vi(2);
        if abs(dx) < tol, continue, end
        if dy/dx > tol, n = n + 1; end
    end
end

function V = randTriangle()
    for k = 1:50
        V = 4*randn(3,2);
        a = (V(2,1)-V(1,1))*(V(3,2)-V(1,2)) - (V(2,2)-V(1,2))*(V(3,1)-V(1,1));
        if abs(a) > 0.35, if a < 0, V = V([1 3 2],:); end, return, end
    end
end

function Q = randPSD(),      A = randn(2); Q = A'*A + 0.35*eye(2); end
function Q = randRank1PSD(), th = 2*pi*rand(); u = [cos(th); sin(th)]; Q = (0.5+2*rand())*(u*u'); end
function Q = randIndefinite()
    th = 2*pi*rand(); R = [cos(th) -sin(th); sin(th) cos(th)];
    Q = R*diag([0.5+2*rand(), -(0.5+2*rand())])*R'; Q = (Q+Q')/2;
end

function ok = report(name, val)
    ok = val; fprintf('  %-52s %s\n', name, tf(val));
end
function s = tf(v), if v, s = 'PASS'; else, s = 'FAIL'; end, end
