classdef envelopeSubdivisionSlowTest < matlab.unittest.TestCase
% The conic subdivision of the CONVEX two-triangle fixture, at zero tolerance -- moved 2026-10-20 from
% envelopeAssembleTest (normal), where its ~94 s kept the bucket over budget (TODO item 5).
%
% BUCKET: slow.

    methods (Test)
        function theARRANGEMENTBuildsASubdivisionThatReproducesFExactly(testCase)
        % THE ROUTE THAT WORKS. Both fold orders are ruled out by measurement; building the
        % ARRANGEMENT once never forms a partial max, so nothing can blow up on the way. On the
        % convex fixture the subdivision must reproduce f itself, at zero tolerance.
            [fA, fB] = fixture(0);
            env = envelopePrune(envelopeAssemble(conjQ(fA), conjQ(fB)));
            [g, info] = envelopeSubdivision(env);
            testCase.verifyGreaterThan(g.nf, 0, 'the subdivision must have faces');
            testCase.verifyLessThan(info.cells, 500, ...
                'the arrangement must stay small -- that is the whole point of not folding');
            rng(43);
            worst = 0;  n = 0;
            for i = 1:150
                x = [rand, rand];
                want = min(fA.eval(x), fB.eval(x));
                if ~isfinite(want), continue, end
                got = g.eval(x);
                testCase.verifyTrue(isfinite(got), 'the subdivision is finite on the domain');
                n = n + 1;
                worst = max(worst, abs(got - want));
            end
            testCase.verifyGreaterThan(n, 80, 'the check was actually exercised');
            testCase.verifyEqual(worst, 0, 'AbsTol', 0, ...
                'a convex input is its own envelope, and the SUBDIVISION must say so exactly');
        end
    end
end

% --------------------------------------------------------------------------------------------------
function [fA, fB] = fixture(c)
% envelopeAssembleTest's standard two-triangle fixture: one quadratic on each half of the unit
% square, perturbed by c times the linear form vanishing on the shared diagonal.
    fA = tri([0 0; 1 0; 1 1], 0);
    fB = tri([0 0; 1 1; 0 1], c);
end

function f = tri(V, c)
    q = [0 0 0 0, 2, 0, 2, 0, 0, 0];
    q(8) = q(8) - c;
    q(9) = q(9) + c;
    f = QuaPol(V, [1 2 1; 2 3 1; 3 1 1], q, [1 0; 1 0; 1 0]);
end
