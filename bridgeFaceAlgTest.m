classdef bridgeFaceAlgTest < matlab.unittest.TestCase
% Tests for bridgeFaceAlg -- the producer that fills AlgCon.FP/.FB.
%
% BUCKET: fast. Everything here is closed-form numerics on 2x2 data; no symbolic call, no solver.
%
% THE ORACLE IS THE DEFINITION, AND IT SHARES NO ARITHMETIC WITH THE SUBJECT. The face is
%
%       z(x) = max { <s,x> - g1(s) : g1(s) = g2(s) } ,
%
% so the oracle maximises exactly that, over a parameterisation of the tie curve built from the
% quadratics themselves -- a LINE solved in closed form when the two Hessians agree, an ELLIPSE
% swept in theta and refined with fminbnd otherwise. bridgeFaceAlg reaches the same number through
% a discriminant in the Lagrange multiplier and a Sturm count; the two meet only at the answer.
%
% WHAT IS PINNED HERE, beyond "it runs":
%   - p VANISHES on constructed face points and does NOT vanish just off them. A polynomial that
%     vanishes everywhere would pass the first half alone, which is why the second half exists.
%   - the REGION equals the set where the maximiser is on the tie, against an independent
%     classifier that maximises each piece separately.
%   - the value is the TOP branch, and -- the finding that forced the negative index -- the
%     ASCENDING index of that same value is NOT constant on the face. That last test is the one
%     that goes red if anyone "tidies" FB back to a positive number.

    methods (Test)

        % ---- the polynomial ----------------------------------------------------------------

        function theImplicitPolynomialVanishesOnTheFaceAndNotBesideIt(testCase)
        % THE TWO HALVES ARE ASKED IN DIFFERENT CURRENCIES, and that is deliberate.
        %   ON the face, p is EVALUATED: the residual relative to the sum of its own term
        %   magnitudes is the sharpest statement available, and it is about p alone.
        %   OFF the face, a residual proves nothing -- the leading coefficient in z can be small
        %   and then |p| is small far from the surface too (measured: 1.1e-9 a full 0.1 above it).
        %   So that half asks about ROOTS: is z + 0.1 clear of every root of p(x,y,.)? It is a
        %   question about the surface, so it is asked of the surface.
        % The root DISTANCE is not used for the first half because it would then be measuring
        % `roots`, whose accuracy on a quartic with 1e6 coefficients is about 1e-8 -- the
        % root-finder's error, not the polynomial's.
            for f = fixtures()
                fx = f{1}; %#ok<*PFBNS>
                face = bridgeFaceAlg(fx.n1, fx.d1, fx.n2, fx.d2);
                [onFace, offFace] = faceResiduals(fx, face.FP, 0.1);
                testCase.verifyLessThan(onFace, 1e-12, ...
                    sprintf('%s: p must vanish at every constructed face point', fx.name));
                testCase.verifyGreaterThan(offFace, 0.05, ...
                    sprintf('%s: 0.1 above the face must be clear of every root', fx.name));
            end
        end

        function theDegenerateTieTakesTheQUADRATICPathInTheMultiplier(testCase)
        % Two pieces with the SAME Hessian tie along a LINE, not a conic, and the mu^3 coefficient
        % is half the determinant of the tie's own 3x3 matrix -- zero exactly then. This is the
        % common case, not an exotic one: any two translates of one quadratic are such a pair, so
        % the quadratic branch of the discriminant has to work.
            fx = fixtureByName('sameHessian');
            [~, info] = bridgeFaceAlg(fx.n1, fx.d1, fx.n2, fx.d2);
            testCase.verifyEqual(info.muDegree, 2, ...
                'a degenerate tie must drop Psi to degree 2 in mu');
        end

        % ---- the region --------------------------------------------------------------------

        function theRegionIsWhereTheMaximiserSitsOnTheTie(testCase)
            for f = fixtures()
                fx = f{1};
                face = bridgeFaceAlg(fx.n1, fx.d1, fx.n2, fx.d2);
                [agree, disagree] = regionAgreement(fx, face, 1500);
                testCase.verifyGreaterThan(agree, 100, ...
                    sprintf('%s: too few points to conclude anything', fx.name));
                testCase.verifyEqual(disagree, 0, ...
                    sprintf('%s: the H-form and the maximiser disagree somewhere', fx.name));
            end
        end

        % ---- the value, and the branch ------------------------------------------------------

        function theValueIsTheTopBranchOfTheImplicitPolynomial(testCase)
            for f = fixtures()
                fx = f{1};
                face = bridgeFaceAlg(fx.n1, fx.d1, fx.n2, fx.d2);
                testCase.verifyEqual(face.FB, -1, ...
                    sprintf('%s: a supremum is the top branch', fx.name));
                [worst, n] = valueAgreement(fx, face);
                testCase.verifyGreaterThan(n, 50, ...
                    sprintf('%s: too few region points to conclude', fx.name));
                testCase.verifyLessThan(worst, 1e-9, ...
                    sprintf('%s: algBranch(FP,-1) must reproduce the definition', fx.name));
            end
        end

        function theASCENDINGIndexCannotNameThisFace(testCase)
        % THE MEASUREMENT THAT FORCED THE NEGATIVE INDEX, pinned as a test rather than left in a
        % document. Two real roots can collide and LEAVE the real line inside the face, which drops
        % the count by two and shifts every index above the departed pair. So the value's ascending
        % index takes more than one value on ONE face, while its distance from the TOP never
        % changes. Anything that re-spells FB as a positive index makes this red.
            fx = fixtureByName('anisotropic');
            face = bridgeFaceAlg(fx.n1, fx.d1, fx.n2, fx.d2);
            idx = ascendingIndexOfValue(fx, face);
            testCase.verifyGreaterThan(numel(unique(idx)), 1, ...
                'the ascending index must be seen to VARY -- that is the whole reason FB is -1');
            testCase.verifyEqual(unique(idx(:)) > 0, true(numel(unique(idx)),1), ...
                'and every one of them is a real index, not a miss');
        end

        % ---- the deliverable: it fills an AlgCon -----------------------------------------------

        function theFaceFillsAnAlgConAndEvaluatesFromItAlone(testCase)
        % What the item was FOR. AlgCon has carried FP/FB since 2026-09-13 with nothing to put in
        % them; this is the first object built from a produced face, and its value is read back
        % from the payload alone -- no formula, no provenance.
            fx = fixtureByName('isotropic');
            face = bridgeFaceAlg(fx.n1, fx.d1, fx.n2, fx.d2);
            a = AlgCon();
            a.FP = {face.FP};
            a.FB = face.FB;
            testCase.verifyEqual(a.FB, -1);
            x = face.wit(1:2) / face.wit(3);
            got  = algBranch(a.FP{1}, x(:).', a.FB);
            want = maxOverTie(fx, x(:));
            testCase.verifyEqual(got, want, 'RelTol', 1e-9, ...
                'the value must come back from FP/FB alone');
        end

        % ---- what it refuses, by name ------------------------------------------------------

        function aVERTEXCellIsRefusedByName(testCase)
        % An AFFINE face means grad g is constant: its end of the sweep is a POINT, the region is a
        % fan rather than a band, and this routine's two-conic shape does not hold it.
            testCase.verifyError(@() bridgeFaceAlg([0 0 0 0 1 0 1 0 0 0], 1, ...
                                                   [0 0 0 0 0 0 0 1 2 3], 1), ...
                'bridgeFaceAlg:vertexCell');
        end

        function anEDGECellIsRefusedByName(testCase)
        % A rank-1 Hessian collapses the plane onto a line, so that boundary is not a conic either.
            testCase.verifyError(@() bridgeFaceAlg([0 0 0 0 1 0 1 0 0 0], 1, ...
                                                   [0 0 0 0 1 0 0 1 2 3], 1), ...
                'bridgeFaceAlg:edgeCell');
        end

        function TWOFacesThatCannotTIEAreRefusedByName(testCase)
        % Same Hessian AND same linear part: g1 - g2 is a CONSTANT, so { g1 = g2 } is empty or the
        % whole plane. There is no tie curve and so no bridge face, and Psi degenerates in mu.
            testCase.verifyError(@() bridgeFaceAlg([0 0 0 0 1 0 1 2 3 0], 1, ...
                                                   [0 0 0 0 1 0 1 2 3 5], 1), ...
                'bridgeFaceAlg:degreeInMu');
        end

        function aTIEWithNoREALPointsIsRefusedByName(testCase)
        % Two concentric quadratics that never take the same value: the tie conic is imaginary, so
        % there is no curve to sweep and no bridge. The pair filter would not have sent this here,
        % which is why the refusal names the caller's mistake rather than guessing a region.
            testCase.verifyError(@() bridgeFaceAlg([0 0 0 0 1 0 1 0 0 0], 1, ...
                                                   [0 0 0 0 2 0 2 0 0 5], 1), ...
                'bridgeFaceAlg:noWitness');
        end

        function AN2POWER53CeilingIsANamedREFUSALAndNotARawOverflow(testCase)
        % The one input class this producer cannot answer, and it must say so in its own words: the
        % branch certification substitutes the witness into FP and clears the denominator to the
        % total degree, so the Sturm chain can pass 2^53 while everything else succeeded. The
        % refusal names the ceiling; it does not let ratQ:overflow escape as if the face were
        % malformed.
            testCase.verifyError(@() bridgeFaceAlg([0 0 0 0 2 0 3 1 -1 0], 2, ...
                                                   [0 0 0 0 5 -1 4 -3 2 1], 3), ...
                'bridgeFaceAlg:branchOverflow');
        end

        function aCUBICFaceIsRefusedByName(testCase)
        % A conjugate face is quadratic. A cubic here means the caller passed something that is not
        % a conjugate cell, and that is worth saying at the boundary.
            testCase.verifyError(@() bridgeFaceAlg([1 0 0 0 1 0 1 0 0 0], 1, ...
                                                   [0 0 0 0 2 0 2 1 1 1], 1), ...
                'bridgeFaceAlg:cubicFace');
        end

        function theARGUMENTSAreChecked(testCase)
            testCase.verifyError(@() bridgeFaceAlg([0 0 0 0 1 0 1 0 0 0], 1), ...
                'bridgeFaceAlg:nargin');
            testCase.verifyError(@() bridgeFaceAlg([1 2 3], 1, [0 0 0 0 1 0 1 0 0 0], 1), ...
                'bridgeFaceAlg:shape');
            testCase.verifyError(@() bridgeFaceAlg([0 0 0 0 1 0 1 0 0 0], 0, ...
                                                   [0 0 0 0 2 0 2 1 1 1], 1), ...
                'bridgeFaceAlg:denominator');
        end
    end
end

% ==================================================================================================
% fixtures
% ==================================================================================================
function F = fixtures()
    F = {fixtureByName('isotropic'), fixtureByName('anisotropic'), ...
         fixtureByName('bothShifted'), fixtureByName('sameHessian')};
end

function fx = fixtureByName(name)
    switch name
        case 'isotropic'                        % A2 = 2*A1: the concentric-circle instance
            fx = mk(name, [0 0 0 0 1 0 1 0 0 0], 1, [0 0 0 0 4 0 4 -4 2 3], 2);
        case 'anisotropic'                      % a genuinely different Hessian: degree 4 in z
            fx = mk(name, [0 0 0 0 1 0 1 0 0 0], 1, [0 0 0 0 3 1 2 -4 2 3], 1);
        case 'bothShifted'                      % neither piece centred, both denominators > 1
            fx = mk(name, [0 0 0 0 4 1 3 2 5 -2], 3, [0 0 0 0 2 0 2 -1 1 1], 1);
        case 'sameHessian'                      % the tie is a LINE: the degenerate-conic path
            fx = mk(name, [0 0 0 0 2 0 2 0 0 0], 1, [0 0 0 0 2 0 2 -3 1 2], 1);
        otherwise
            error('bridgeFaceAlgTest:fixture', 'no fixture named %s', name);
    end
end

function fx = mk(name, n1, d1, n2, d2)
    fx = struct('name', name, 'n1', n1, 'd1', d1, 'n2', n2, 'd2', d2, ...
                'A1', [n1(5) n1(6); n1(6) n1(7)]/d1, 'b1', [n1(8); n1(9)]/d1, 'k1', n1(10)/d1, ...
                'A2', [n2(5) n2(6); n2(6) n2(7)]/d2, 'b2', [n2(8); n2(9)]/d2, 'k2', n2(10)/d2);
end

% ==================================================================================================
% the oracle: the tie curve, parameterised from the two quadratics
% ==================================================================================================
function [s, tangent] = tiePoint(fx, t)
% A point of { g1 = g2 } at parameter t, and the parameterisation's own derivative is not needed --
% the maximisations below are over t directly.
    [Dl, dl, kp] = tieData(fx);
    if all(Dl(:) == 0)
        % a LINE: dl's + kp = 0
        s0  = -kp * dl / (dl.'*dl);
        dir = [-dl(2); dl(1)];
        s = s0 + t*dir;  tangent = dir;
        return
    end
    [V, L, R] = ellipseFrame(fx);
    s = tieCentre(fx) + V*[sqrt(R/L(1))*cos(t); sqrt(R/L(2))*sin(t)];
    tangent = V*[-sqrt(R/L(1))*sin(t); sqrt(R/L(2))*cos(t)];
end

function [Dl, dl, kp] = tieData(fx)
    Dl = fx.A1 - fx.A2;  dl = fx.b1 - fx.b2;  kp = fx.k1 - fx.k2;
end

function c0 = tieCentre(fx)
    [Dl, dl, ~] = tieData(fx);
    c0 = -Dl \ dl;
end

function [V, L, R] = ellipseFrame(fx)
% The tie in centred form: with u = s - c0, the tie is u'Dl u = R. When Dl is definite and R has
% its sign, that is an ellipse and (V, L) is the eigenframe of the matrix that makes it one.
    [Dl, dl, kp] = tieData(fx);
    c0 = tieCentre(fx);
    R  = -2*((c0.'*Dl*c0)/2 + dl.'*c0 + kp);
    M  = Dl;
    if R < 0, M = -Dl;  R = -R; end
    [V, D] = eig(M);
    L = diag(D);
    if any(L <= 0)
        error('bridgeFaceAlgTest:notAnEllipse', ...
            ['the oracle parameterises a LINE or an ELLIPSE, and this tie is neither. A ' ...
             'hyperbolic tie needs its own parameterisation rather than a silent sample.']);
    end
end

function z = maxOverTie(fx, x)
% max over the tie of <s,x> - g1(s): the definition, maximised over the parameterisation above.
    g1 = @(s) (s.'*fx.A1*s)/2 + fx.b1.'*s + fx.k1;
    obj = @(t) -objAt(fx, g1, t, x);                       % fminbnd MINIMISES
    [Dl, ~, ~] = tieData(fx);
    if all(Dl(:) == 0)
        % a line: the objective is a concave quadratic in t and its maximiser is closed form
        [s0, dir] = deal(tiePoint(fx, 0), lineDir(fx));
        a = -(dir.'*fx.A1*dir)/2;
        b = dir.'*x - (s0.'*fx.A1*dir + fx.b1.'*dir);
        if a >= 0, z = inf; return, end                    % unbounded along the tie
        t = -b/(2*a);
        z = objAt(fx, g1, t, x);
        return
    end
    ts = linspace(0, 2*pi, 721);
    vals = arrayfun(@(t) objAt(fx, g1, t, x), ts);
    [~, i] = max(vals);
    lo = ts(max(i-1,1));  hi = ts(min(i+1, numel(ts)));
    tb = fminbnd(obj, lo, hi, optimset('TolX', 1e-14));
    z = max(objAt(fx, g1, tb, x), vals(i));
end

function d = lineDir(fx)
    [~, dl, ~] = tieData(fx);
    d = [-dl(2); dl(1)];
end

function v = objAt(fx, g1, t, x)
    s = tiePoint(fx, t);
    v = s.'*x - g1(s);
end

% ==================================================================================================
% the checks the tests read
% ==================================================================================================
function [onFace, offFace] = faceResiduals(fx, FP, h)
% CONSTRUCTED face points: at a tie point s the whole segment between the two gradients lies on the
% face, with value <s,x> - g1(s). Nothing here solves anything, so it cannot agree with the subject
% by accident.
% [output] onFace  : the worst RELATIVE residual |p| at a constructed face point
%          offFace : the smallest distance from (that value + h) to ANY root of p(x,y,.)
    g1 = @(s) (s.'*fx.A1*s)/2 + fx.b1.'*s + fx.k1;
    [Dl, ~, ~] = tieData(fx);
    if all(Dl(:) == 0), ts = linspace(-3, 3, 9); else, ts = linspace(0, 2*pi, 9); end
    onFace = 0;  offFace = inf;
    for t = ts
        s = tiePoint(fx, t);
        for lam = [0 0.25 0.5 0.75 1]
            x = (1-lam)*(fx.A1*s + fx.b1) + lam*(fx.A2*s + fx.b2);
            z = s.'*x - g1(s);
            [~, info] = algBranch(FP, x(:).', -1);
            rr = info.allReal{1};
            if isempty(rr)
                onFace = inf;  offFace = 0;  return         % no root at all: a definite failure
            end
            onFace  = max(onFace,  relEval(FP, x, z));
            offFace = min(offFace, min(abs(rr - (z + h))));
        end
    end
end

function [agree, disagree] = regionAgreement(fx, face, n)
    ctr = face.wit(1:2)/face.wit(3);
    rng(11);
    agree = 0;  disagree = 0;
    for i = 1:n
        x = ctr(:) + (rand(2,1)*2-1)*5;
        if inHForm(face, x) == maximiserIsOnTheTie(fx, x)
            if inHForm(face, x), agree = agree + 1; end
        else
            disagree = disagree + 1;
        end
    end
end

function [worst, n] = valueAgreement(fx, face)
    ctr = face.wit(1:2)/face.wit(3);
    rng(13);
    worst = 0;  n = 0;
    for i = 1:600
        x = ctr(:) + (rand(2,1)*2-1)*5;
        if ~inHForm(face, x), continue, end
        want = maxOverTie(fx, x);
        if ~isfinite(want), continue, end
        got = algBranch(face.FP, x(:).', face.FB);
        n = n + 1;
        worst = max(worst, abs(got - want)/max(1, abs(want)));
    end
end

function idx = ascendingIndexOfValue(fx, face)
    ctr = face.wit(1:2)/face.wit(3);
    rng(17);
    idx = [];
    for i = 1:800
        x = ctr(:) + (rand(2,1)*2-1)*5;
        if ~inHForm(face, x), continue, end
        want = maxOverTie(fx, x);
        if ~isfinite(want), continue, end
        [~, info] = algBranch(face.FP, x(:).', -1);
        rr = info.allReal{1};
        if isempty(rr), continue, end
        [gap, j] = min(abs(rr - want));
        if gap > 1e-8*max(1, abs(want)), continue, end
        idx(end+1) = j; %#ok<AGROW>
    end
end

function tf = inHForm(face, x)
    tf = true;
    for r = 1:size(face.con, 1)
        c = face.con(r,1:6);
        v = c(1)*x(1)^2 + c(2)*x(1)*x(2) + c(3)*x(2)^2 + c(4)*x(1) + c(5)*x(2) + c(6);
        if face.con(r,7)*v < 0, tf = false; return, end
    end
end

function tf = maximiserIsOnTheTie(fx, x)
% The independent classifier: maximise <s,x> - g_k(s) for each piece SEPARATELY and ask whether
% either unconstrained maximiser is feasible for its own piece. If neither is, the maximiser of the
% max is on the tie. No conic, no polynomial, no multiplier.
    s1 = fx.A1 \ (x - fx.b1);
    if gval(fx.A1, fx.b1, fx.k1, s1) >= gval(fx.A2, fx.b2, fx.k2, s1) - 1e-12
        tf = false;  return
    end
    s2 = fx.A2 \ (x - fx.b2);
    if gval(fx.A2, fx.b2, fx.k2, s2) >= gval(fx.A1, fx.b1, fx.k1, s2) - 1e-12
        tf = false;  return
    end
    tf = true;
end

function v = gval(A, b, k, s), v = (s.'*A*s)/2 + b.'*s + k; end

function r = relEval(P, x, z)
% |p| at (x,z), RELATIVE to the sum of the term magnitudes: a residual means nothing beside
% coefficients of size 1e6 unless it is read against the terms that produced it.
    sz = size(P);  sz(end+1:3) = 1;
    v = 0;  s = 0;
    for i = 1:sz(1)
        for j = 1:sz(2)
            for k = 1:sz(3)
                c = P(i,j,k);
                if c == 0, continue, end
                t = c * x(1)^(i-1) * x(2)^(j-1) * z^(k-1);
                v = v + t;  s = s + abs(t);
            end
        end
    end
    if s == 0, s = 1; end
    r = abs(v)/s;
end
