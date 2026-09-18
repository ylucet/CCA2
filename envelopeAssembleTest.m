classdef envelopeAssembleTest < matlab.unittest.TestCase
% Tests for envelopeAssemble / envelopeEval -- the convex envelope, assembled from the pieces'
% conjugates and evaluated as the MAX over its faces.
%
% BUCKET: normal. The assembly itself is fast (0.3 s convex, ~30 s with bridge faces, dominated by
% the bridge witness search); what costs is the ORACLE for the non-convex fixtures, which maximises
% <s,x> - max(gA,gB)(s) numerically on a fine grid.
%
% TWO ORACLES, AND NEITHER SHARES ANYTHING WITH THE SUBJECT.
%
%   c = 0 is CONVEX, so f** = f. That is a PERFECT oracle -- the input itself -- and it is the one
%   that can be asserted at zero tolerance.
%
%   c < 0 is not, so the oracle is a numeric sup over s. Note which direction that is rigorous in:
%   ANY s gives a LOWER bound on f**, so the oracle can only UNDER-report. The assembly must
%   therefore never fall below it; and it cannot rise above the truth by construction, because every
%   face's value is attained at an s that is FEASIBLE for its own cell pair. Both directions are
%   still checked, the second at a tolerance that reflects the oracle's own accuracy.

    methods (Test)

        function theCONVEXFixtureReproducesFExactly(testCase)
        % The strongest statement available anywhere in this file: where the input is convex the
        % envelope IS the input, so the comparison is against f itself and the tolerance is ZERO.
        % It exercises the whole chain -- the pair walk, the tangential rule, the candidates, their
        % regions and the max -- with no bridge face involved.
            [fA, fB] = fixture(0);
            [env, info] = envelopeAssemble(conjQ(fA), conjQ(fB));
            testCase.verifyEqual(info.nCurved, 0, 'a convex fixture needs no bridge face');
            testCase.verifyGreaterThan(numel(env.quad), 10, 'and it does produce faces');
            rng(13);
            worst = 0;  n = 0;
            for i = 1:100
                x = [rand, rand];
                want = min(fA.eval(x), fB.eval(x));
                if ~isfinite(want), continue, end
                got = envelopeEval(env, x);
                testCase.verifyTrue(isfinite(got), 'the envelope is finite on the domain');
                n = n + 1;
                worst = max(worst, abs(got - want));
            end
            testCase.verifyGreaterThan(n, 50, 'the check was actually exercised');
            testCase.verifyEqual(worst, 0, 'AbsTol', 0, ...
                'a convex input is its own envelope, exactly');
        end

        function theNONCONVEXFixturesMatchTheDefinition(testCase)
        % With bridge faces in play. The oracle is a numeric sup, so it can only under-report --
        % "no undershoot" is the rigorous direction and is asserted tightly; the other direction is
        % asserted at the oracle's own accuracy.
            for c = [-1 -2]
                [fA, fB] = fixture(c);
                gA = conjQ(fA);  gB = conjQ(fB);
                [env, info] = envelopeAssemble(gA, gB);
                testCase.verifyGreaterThan(info.nCurved, 0, ...
                    sprintf('c = %g must need a bridge face', c));
                testCase.verifyEqual(numel(env.bridge), info.nCurved, ...
                    'every curved pair contributes exactly one bridge face');
                rng(17);
                worst = 0;  n = 0;
                for i = 1:8
                    x = [rand, rand];
                    want = supOracle(gA, gB, x);
                    if ~isfinite(want), continue, end
                    got = envelopeEval(env, x);
                    testCase.verifyTrue(isfinite(got), 'the envelope is finite on the domain');
                    testCase.verifyGreaterThanOrEqual(got, want - 1e-9*max(1,abs(want)), ...
                        sprintf('c = %g: the envelope may never fall BELOW a feasible s', c));
                    n = n + 1;
                    worst = max(worst, abs(got - want)/max(1, abs(want)));
                end
                testCase.verifyGreaterThan(n, 4, 'the check was actually exercised');
                testCase.verifyLessThan(worst, 1e-8, ...
                    sprintf('c = %g: the envelope must equal the definition', c));
            end
        end

        function aBRIDGEFaceIsUSEDAndNotJustBuilt(testCase)
        % The faces must actually ANSWER somewhere, or the fixture would be proving nothing about
        % them. `envelopeEval` reports which face supplied each value, so this asks it directly.
            [fA, fB] = fixture(-2);
            env = envelopeAssemble(conjQ(fA), conjQ(fB));
            rng(19);
            X = rand(400,2);
            [~, who] = envelopeEval(env, X);
            kinds = {who.kind};
            testCase.verifyGreaterThan(sum(strcmp(kinds, 'bridge')), 0, ...
                'some point must be answered by a BRIDGE face');
            % AND NOT "some by a quadratic one", which was asserted first and is FALSE here: at
            % c = -2 the kink is deep enough that the bridge faces answer the whole square, and
            % their values are right to 3e-12 against the definition. A quadratic face answering
            % somewhere is a property of the CONVEX fixture, where it is checked, not of this one.
            quadFixture = fixture(0);
            envQ = envelopeAssemble(conjQ(quadFixture), conjQ(second(0)));
            [~, whoQ] = envelopeEval(envQ, rand(50,2));
            testCase.verifyGreaterThan(sum(strcmp({whoQ.kind}, 'quad')), 0, ...
                'and a quadratic face must answer where there is no bridge at all');
        end

        function theARGUMENTSAreChecked(testCase)
            [fA, fB] = fixture(0);
            testCase.verifyError(@() envelopeAssemble(conjQ(fA)), 'envelopeAssemble:nargin');
            testCase.verifyError(@() envelopeAssemble(fA, conjQ(fB)), 'envelopeAssemble:type');
            env = envelopeAssemble(conjQ(fA), conjQ(fB));
            testCase.verifyError(@() envelopeEval(env), 'envelopeEval:nargin');
            testCase.verifyError(@() envelopeEval(struct('a',1), [0 0]), 'envelopeEval:shape');
            testCase.verifyError(@() envelopeEval(env, [0 0 0]), 'envelopeEval:points');
        end
    end
end

% ==================================================================================================
function [fA, fB] = fixture(c)
% The standard two-triangle fixture: one quadratic on each half of the unit square, perturbed by c
% times the linear form vanishing on the shared diagonal -- so the two pieces AGREE there for every
% c and the input is always a function. c = 0 is convex; c < 0 is a kink of increasing depth.
    fA = tri([0 0; 1 0; 1 1], 0);
    fB = tri([0 0; 1 1; 0 1], c);
end

function fB = second(c)
% the SECOND piece alone, for a caller that already has the first
    [~, fB] = fixture(c);
end

function f = tri(V, c)
    q = [0 0 0 0, 2, 0, 2, 0, 0, 0];
    q(8) = q(8) - c;
    q(9) = q(9) + c;
    f = QuaPol(V, [1 2 1; 2 3 1; 3 1 1], q, [1 0; 1 0; 1 0]);
end

% --------------------------------------------------------------------------------------------------
function v = supOracle(gA, gB, x)
% sup_s [ <s,x> - max(gA,gB)(s) ], numerically. A grid to find the basin, then Nelder-Mead from two
% starts -- the objective is a max of quadratics, so it is KINKED, and a single start can stall on
% the kink.
    obj = @(s) -( s(:).'*x(:) - max(ev(gA,s), ev(gB,s)) );
    best = inf;  bs = [0;0];
    for s1 = -8:0.1:8
        for s2 = -8:0.1:8
            o = obj([s1;s2]);
            if o < best, best = o;  bs = [s1;s2]; end
        end
    end
    op = optimset('TolX',1e-14,'TolFun',1e-14,'MaxFunEvals',20000,'Display','off');
    s1 = fminsearch(obj, bs, op);
    s2 = fminsearch(obj, bs + [1e-3; -1e-3], op);
    v = -min([obj(s1), obj(s2), best]);
end

function v = ev(g, s)
    v = g.eval(s(:).');
    if isempty(v) || ~isfinite(v), v = -inf; end
end
