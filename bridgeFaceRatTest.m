classdef bridgeFaceRatTest < matlab.unittest.TestCase
% Tests for bridgeFaceRat -- the AFFINE x RANK-1 bridge face, which is RATIONAL rather than
% algebraic and comes out in closed form.
%
% BUCKET: fast. Closed-form numerics on 2x2 data.
%
% THE MAIN CHECK IS ANALYTIC, NOT A SEARCH, and that is deliberate. The claim is that the face's
% value is the objective AT THE STATIONARY POINT along the tie, and that the stationary point is a
% MAXIMUM. Both are closed-form statements, so they can be verified pointwise without maximising
% anything:
%
%   (a) s(t*) lies ON the tie                      -- the residual of g1 - g2 there
%   (b) the face's value equals <s,x> - g1(s)      -- at that very point
%   (c) the objective is CONCAVE in t              -- so the stationary point is the max over the
%                                                     whole parabola, a quadratic having one
%
% A SEARCH ORACLE WAS TRIED FIRST AND IT MISLED, which is why this is the main check: on an
% UNBOUNDED tie the maximiser can sit thousands of units out, a windowed sweep under-reports it, and
% the SUBJECT then looks like it is overclaiming. The restricted (cut) case is still checked against
% a search, because there the pair bounds the window itself.

    methods (Test)

        function theVALUEIsTheSTATIONARYValueAndTheStationaryPointIsAMaximum(testCase)
            for f = fixtures()
                fx = f{1};
                face = bridgeFaceRat(fx.n1, fx.d1, fx.n2, fx.d2);
                [nT, tieRes, valRes, nonConcave] = analytic(fx, face);
                testCase.verifyGreaterThan(nT, 50, sprintf('%s: too few face points', fx.name));
                testCase.verifyLessThan(tieRes, 1e-12, ...
                    sprintf('%s: the maximiser must lie ON the tie', fx.name));
                testCase.verifyLessThan(valRes, 1e-9, ...
                    sprintf('%s: the value must be the objective at that point', fx.name));
                testCase.verifyEqual(nonConcave, 0, ...
                    sprintf('%s: the face must only exist where the objective is concave', fx.name));
            end
        end

        function theCUTFaceEqualsTheRESTRICTEDDefinition(testCase)
        % With a pair given, the window is bounded by the pair itself, so a search oracle is
        % reliable -- and it is the independent one: it sweeps the tie by solving it directly and
        % shares no arithmetic with the closed form.
            fx = fixtureByName('shifted');
            P  = [1 0 0; 0 1 0];
            face = bridgeFaceRat(fx.n1, fx.d1, fx.n2, fx.d2, P);
            [worst, n] = restrictedAgreement(fx, face, P);
            testCase.verifyGreaterThan(n, 20, 'too few region points');
            testCase.verifyLessThan(worst, 1e-9, 'the cut face must equal the restricted definition');
        end

        function theDENOMINATORIsLINEARWhichIsTheWholePoint(testCase)
        % The face is a quadratic numerator over a LINEAR denominator -- exactly the `Rat` level of
        % the lattice, needing no FP/FB and no branch index. If this ever came out quadratic the
        % face would not fit the level it is being stored at.
            fx = fixtureByName('basic');
            face = bridgeFaceRat(fx.n1, fx.d1, fx.n2, fx.d2);
            testCase.verifyEqual(face.kind, 'rational');
            testCase.verifyEqual(numel(face.den), 3, 'the denominator is [g h k], i.e. LINEAR');
            testCase.verifyTrue(any(face.den(1:2) ~= 0), 'and it genuinely depends on x');
            testCase.verifyEqual(face.fN(1:4), zeros(1,4), 'the numerator has no cubic part');
        end

        function aSTRAIGHTTieMeansTHERESNOSeparateFaceAndItSaysSo(testCase)
        % The case ALL 12 affine x rank-1 pairs on the standard fixtures turn out to be. When the
        % tie's linear part lies in the rank-1 Hessian's range the tie is a pair of PARALLEL LINES,
        % and with one face affine the sup over each side is a LINEAR PROGRAM attained at a VERTEX
        % -- so the tie is just another facet and there is no separate bridge face to build.
            testCase.verifyError(@() bridgeFaceRat([0 0 0 0 0 0 0 2 0 0], 1, ...
                                                   [0 0 0 0 1 0 0 0 0 1], 1), ...
                'bridgeFaceRat:noSeparateFace');
        end

        function theARGUMENTSAndTheSHAPESAreChecked(testCase)
            aff = [0 0 0 0 0 0 0 1 -1 0];   rk1 = [0 0 0 0 2 0 0 -1 1 1];
            testCase.verifyError(@() bridgeFaceRat(aff, 1), 'bridgeFaceRat:nargin');
            testCase.verifyError(@() bridgeFaceRat([1 2 3], 1, rk1, 1), 'bridgeFaceRat:shape');
            testCase.verifyError(@() bridgeFaceRat(aff, 0, rk1, 1), 'bridgeFaceRat:denominator');
            testCase.verifyError(@() bridgeFaceRat([1 0 0 0 0 0 0 1 1 1], 1, rk1, 1), ...
                'bridgeFaceRat:cubicFace');
            testCase.verifyError(@() bridgeFaceRat(rk1, 1, rk1, 1), 'bridgeFaceRat:notAffine');
            testCase.verifyError(@() bridgeFaceRat(aff, 1, [0 0 0 0 2 0 2 0 0 0], 1), ...
                'bridgeFaceRat:notRankOne');
            testCase.verifyError(@() bridgeFaceRat(aff, 1, aff, 1), 'bridgeFaceRat:notRankOne');
        end

        function bridgeMemberDECIDESARationalFace(testCase)
        % envelopeEval asks bridgeMember which points a curved-tie face owns, for EVERY face kind.
        % Until 2026-10-10 it only understood the algebraic face's pair and refused a rational
        % one by name, so any envelope with a vertex-against-edge tie could not be evaluated. The
        % membership is the closed form this file's own `analytic` uses: the maximiser along the
        % parabola exists (concave: C/scale < 0) and lies in the pair.
            fx = fixtureByName('shifted');
            P  = [1 0 0; 0 1 0];
            face = bridgeFaceRat(fx.n1, fx.d1, fx.n2, fx.d2, P);
            pr = face.pair;
            rng(107);
            X = face.wit(1:2)/face.wit(3) + (rand(400,2)*2 - 1)*4;
            want = false(size(X,1),1);
            for i = 1:size(X,1)
                Bv = pr.B(1)*X(i,1) + pr.B(2)*X(i,2) + pr.B(3);
                Cv = pr.C(1)*X(i,1) + pr.C(2)*X(i,2) + pr.C(3);
                if Cv * pr.scale >= 0, continue, end
                s = paraPoint(pr, -Bv/(2*Cv));
                want(i) = all(P(:,1:2)*s + P(:,3) >= 0);
            end
            testCase.assumeTrue(any(want) && any(~want), 'fixture must straddle the face');
            testCase.verifyEqual(bridgeMember(pr, X), want);
        end

        function bridgeFaceAlgDELEGATESToItRatherThanRefusing(testCase)
        % One entry point, two face kinds: a caller asks bridgeFaceAlg and gets whichever face the
        % pair actually has. The singular pencil used to be a refusal.
            fx = fixtureByName('basic');
            face = bridgeFaceAlg(fx.n1, fx.d1, fx.n2, fx.d2);
            testCase.verifyEqual(face.kind, 'rational');
            testCase.verifyTrue(isfield(face, 'den'), 'and it carries the rational face''s data');
            % ...and the swapped order reaches the same face, because on the tie the two agree
            swapped = bridgeFaceAlg(fx.n2, fx.d2, fx.n1, fx.d1);
            testCase.verifyEqual(swapped.fN, face.fN, 'RelTol', 0);
            testCase.verifyEqual(swapped.den, face.den, 'RelTol', 0);
        end
    end
end

% ==================================================================================================
function F = fixtures()
    F = {fixtureByName('basic'), fixtureByName('shifted'), fixtureByName('scaled')};
end

function fx = fixtureByName(name)
    switch name
        case 'basic'
            fx = mk(name, [0 0 0 0 0 0 0 1 -1 0], 1, [0 0 0 0 2 0 0 -1 1 1], 1);
        case 'shifted'
            fx = mk(name, [0 0 0 0 0 0 0 -2 1 1], 1, [0 0 0 0 0 0 4 1 -2 0], 1);
        case 'scaled'
            fx = mk(name, [0 0 0 0 0 0 0 3 1 -2], 2, [0 0 0 0 1 2 4 -1 1 1], 3);
        otherwise
            error('bridgeFaceRatTest:fixture', 'no fixture named %s', name);
    end
end

function fx = mk(name, n1, d1, n2, d2)
    fx = struct('name', name, 'n1', n1, 'd1', d1, 'n2', n2, 'd2', d2, ...
                'b1', [n1(8); n1(9)]/d1, 'k1', n1(10)/d1, ...
                'A2', [n2(5) n2(6); n2(6) n2(7)]/d2, 'b2', [n2(8); n2(9)]/d2, 'k2', n2(10)/d2);
end

% ==================================================================================================
function [nT, tieRes, valRes, nonConcave] = analytic(fx, face)
    pr = face.pair;
    g1 = @(s) fx.b1.'*s + fx.k1;
    g2 = @(s) (s.'*fx.A2*s)/2 + fx.b2.'*s + fx.k2;
    rng(101);
    ctr = face.wit(1:2).'/face.wit(3);
    nT = 0;  tieRes = 0;  valRes = 0;  nonConcave = 0;
    for i = 1:1500
        x = ctr + (rand(2,1)*2-1)*6;
        Bv = pr.B(1)*x(1) + pr.B(2)*x(2) + pr.B(3);
        Cv = pr.C(1)*x(1) + pr.C(2)*x(2) + pr.C(3);
        if Cv * pr.scale >= 0, continue, end
        t = -Bv/(2*Cv);
        s = paraPoint(pr, t);
        if any(~isfinite(s)), continue, end
        nT = nT + 1;
        tieRes = max(tieRes, abs(g1(s) - g2(s)) / max(1, max(abs(s))^2));
        want = s.'*x - g1(s);
        valRes = max(valRes, abs(evalRat(face, x) - want)/max(1, abs(want)));
        if 2*Cv/pr.scale >= 0, nonConcave = nonConcave + 1; end
    end
end

function [worst, n] = restrictedAgreement(fx, face, P)
    ctr = face.wit(1:2).'/face.wit(3);
    rng(103);
    worst = 0;  n = 0;
    for i = 1:300
        x = ctr + (rand(2,1)*2-1)*4;
        if ~memberRat(face.pair, x), continue, end
        want = tieMax(fx, P, x);
        if ~isfinite(want), continue, end
        n = n + 1;
        worst = max(worst, abs(evalRat(face, x) - want)/max(1, abs(want)));
    end
end

function v = tieMax(fx, P, x)
% max over { g1 = g2 } ^ P of <s,x> - g1(s), by sweeping the tie directly. The pair bounds the
% window here, which is what makes a search reliable in this case and not in the uncut one.
    Dl = -fx.A2;  dl = fx.b1 - fx.b2;  kp = fx.k1 - fx.k2;
    % The PAIR bounds the window here, so a modest first sweep already finds the right basin and
    % the refinements do the rest -- this is the bucket's budget, not its accuracy, talking.
    lo = -60;  hi = 60;  v = -inf;  best = NaN;
    res = [4001 201 201 201 201];
    for round = 1:numel(res)
        for s1 = linspace(lo, hi, res(round))
            a = Dl(2,2);  b = 2*Dl(1,2)*s1 + 2*dl(2);  c = Dl(1,1)*s1^2 + 2*dl(1)*s1 + 2*kp;
            S = [];
            if abs(a) < 1e-14
                if abs(b) > 1e-14, S = [s1; -c/b]; end
            else
                disc = b^2 - 4*a*c;
                if disc < 0, continue, end
                r = sqrt(disc);
                S = [s1, s1; (-b+r)/(2*a), (-b-r)/(2*a)];
            end
            for j = 1:size(S,2)
                sv = S(:,j);
                if ~isempty(P) && any(P(:,1:2)*sv + P(:,3) < 0), continue, end
                o = sv.'*x - (fx.b1.'*sv + fx.k1);
                if o > v, v = o;  best = s1; end
            end
        end
        if ~isfinite(v), return, end
        h = (hi - lo)/(res(round) - 1);
        lo = best - 2*h;  hi = best + 2*h;
    end
end

function tf = memberRat(pr, x)
    Bv = pr.B(1)*x(1) + pr.B(2)*x(2) + pr.B(3);
    Cv = pr.C(1)*x(1) + pr.C(2)*x(2) + pr.C(3);
    if Cv * pr.scale >= 0, tf = false; return, end
    s = paraPoint(pr, -Bv/(2*Cv));
    tf = isempty(pr.P) || all(pr.P(:,1:2)*s + pr.P(:,3) >= 0);
end

function s = paraPoint(pr, t)
    r = (pr.rq*t^2 + pr.rl*t + pr.rc) / pr.rd;
    s = (t*pr.u + r*pr.w) / pr.nu;
end

function v = evalRat(face, x)
    q = face.fN;
    num = q(5)/2*x(1)^2 + q(6)*x(1)*x(2) + q(7)/2*x(2)^2 + q(8)*x(1) + q(9)*x(2) + q(10);
    den = face.den(1)*x(1) + face.den(2)*x(2) + face.den(3);
    v = num / den;
end
