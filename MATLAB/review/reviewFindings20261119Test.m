classdef reviewFindings20261119Test < matlab.unittest.TestCase
% Findings of the 2026-11-19 code review (reviews/codereview-2026-11-19.md). Each test asserts the
% CORRECT behaviour, so it FAILS against the reviewed code; it is expected to pass once the finding
% is fixed. This file lives in MATLAB/review/, which suite.sh does not run (it lists *Test.m in
% MATLAB/ only). Run it explicitly:
%     cd MATLAB; addpath(pwd); runtests(fullfile('review','reviewFindings20261119Test.m'))

    methods (TestClassSetup)
        function onPath(tc) %#ok<MANU>
            addpath(fileparts(fileparts(mfilename('fullpath'))));
        end
    end

    methods (Test)
        function coneWithInteriorNullDirection(tc)
            % Finding 1 (CRITICAL). f(x,y) = y^2/2 on the cone {x >= |y|}. Along x = (t,0) the
            % objective <s,x> - f(x) = s1*t, so f*(s) = +inf whenever s1 > 0. conjQ only tests the
            % cone's two edge rays for zero curvature; the null direction (1,0) is strictly inside.
            f = QuaPol([0 0; 1 -1; 1 1], [1 2 0; 1 3 0], [0 0 0 0 0 0 1 0 0 0], [1 0; 0 1]);
            g = conjQ(f);
            v = g.eval([1 0; 0.5 0.3; 0.1 0; -1 0]);
            tc.verifyEqual(v(1:3).', [Inf Inf Inf], 'f*(s) must be +inf for s1 > 0');
            tc.verifyEqual(v(4), 0, 'AbsTol', 1e-12);
        end

        function cellInteriorMissesAFarVerticalStrip(tc)
            % Finding 2 (MAJOR). The cell {10 <= x <= 10.5} minus the disc |x-(100,0)| < 1 has
            % interior (e.g. (10.25, 0)). No two of its curves meet inside it, the random probe box
            % is around the origin, and the "complete" CAD fallback skips vertical lines, so the
            % cell is declared EMPTY. Its horizontal twin is (correctly) found.
            EcQ = [0 0 0 1 0 -10; 0 0 0 2 0 -21; 1 0 1 -200 0 9999];
            rows = [1 1; 2 -1; 3 1];
            tc.verifyTrue(cellHasInterior(EcQ, rows));
        end

        function cellInteriorMissesSharedLineFactorCell(tc)
            % Finding 2, second mechanism: four degenerate conics sharing x + y + 3 = 0; conicMeet
            % calls each pair degenerate and returns no points, so the CAD misses abscissae.
            % Reached in a real conj call; (-58.4253, -109.3822) satisfies every row strictly.
            EcQ = [101 76 -25 3378 3000 9225; 91 38 -53 -174 -606 -1341; 64 37 -27 2386 2113 6582;
                   1221 154 -1067 -42386 -49250 -138147; 0 0 0 3 -5 65; 0 0 0 3 -5 -31;
                   0 0 0 -18 8 -119; 0 0 0 9 -4 313; 0 0 0 5 -2 169];
            rows = [1 1; 2 -1; 3 -1; 4 1; 5 1; 6 1; 7 1; 8 1; 9 1];
            tc.verifyTrue(cellHasInterior(EcQ, rows));
        end

        function conjLShapeY2(tc)
            % Finding 8 (CRITICAL): f = y^2/2 on the L-shape; sup at s = (2.53125, -1.453125) is
            % 2*s1 = 5.0625 (attained at (2,0)). conjQ returns (2 s1 + s2)^2/2 = 6.5137939453125.
            V = [0 0; 2 0; 2 1; 1 1; 1 2; 0 2];
            E = [(1:6).' [2:6 1].' ones(6,1)];
            f = QuaPol(V, E, [0 0 0 0 0 0 1 0 0 0], repmat([1 0], 6, 1));
            g = conjQ(f);
            tc.verifyEqual(g.eval([2.53125 -1.453125]), 5.0625, 'AbsTol', 1e-12);
        end

        function biconjInfiniteOutsideTheSquare(tc)
            % Finding 9a (CRITICAL): f = (y^2 - x^2)/2 on the unit square; (-0.49, 0.01) is outside
            % the square, so f** = +inf there.
            f = QuaPol([0 0; 1 0; 1 1; 0 1], [1 2 1; 2 3 1; 3 4 1; 4 1 1], [0 0 0 0 -1 0 1 0 0 0], repmat([1 0], 4, 1));
            h = biconjQ(f);
            tc.verifyEqual(h.eval([-0.4921875 0.01171875]), Inf);
        end

        function biconjFenchelYoungLShape(tc)
            % Finding 9b (CRITICAL): f = -xy + y^2/2 on the L-shape. f**(x) >= <s,x> - f*(s) for
            % every s; with conjQ's own value at s the bound is 1.04606.
            V = [0 0; 2 0; 2 1; 1 1; 1 2; 0 2];
            E = [(1:6).' [2:6 1].' ones(6,1)];
            f = QuaPol(V, E, [0 0 0 0 0 -1 1 0 0 0], repmat([1 0], 6, 1));
            x = [0.2578125 1.76171875];  s = [-1.9484626037100141 1.6789473678974558];
            g = conjQ(f);  h = biconjQ(f);
            tc.verifyGreaterThanOrEqual(h.eval(x), s*x.' - g.eval(s) - 1e-9);
        end

        function biconjFiniteInsideTheDomainR3072(tc)
            % Finding 9c/9d (CRITICAL): two continuous pieces on a quadrilateral (R3-072). x is in
            % piece 1 (f(x) = -12.29), so f** is finite there; certified bracket
            % [-14.68169, -14.67993]. At (-2.49, 1.85) a convex combination of three domain points
            % gives -18.455, so f** <= that. The tie-arc screen parametrised the parabola tie
            % through a spurious 1e34 point and dropped the needed atom pair.
            f = QuaPol([-3 2; -1 0; 2 0; -2 3], [1 2 1; 2 3 1; 3 1 1; 3 4 1; 4 1 1], ...
                       [0 0 0 0 0 2 -2 2 1 1; 0 0 0 0 4 9 8 -4 -8 5], [1 0; 1 0; 1 2; 2 0; 2 0]);
            h = biconjQ(f);
            v = h.eval([-2.4921875 1.51171875; -2.49 1.85]);
            tc.verifyGreaterThanOrEqual(v(1), -14.681688881664497 - 1e-9);
            tc.verifyLessThanOrEqual(v(1), -14.679931640625);
            tc.verifyLessThanOrEqual(v(2), -18.455);
        end

        function biconjFiniteInsideTheDomainR6075(tc)
            % Finding 9c (CRITICAL), small-cell mechanism: two continuous pieces on a two-cone fan
            % (R6-075). x is in piece 2, so f** is finite; certified bracket [1.80731, 1.86682]. The
            % subdivision dropped the live unit-wide cell holding x (a grid over a 1e3 box).
            f = QuaPol([0 0; 1 1; 0 1; -1 2], [1 2 0; 1 3 0; 1 4 0], ...
                       [0 0 0 0 2 1 2 1 1 0; 0 0 0 0 0 -1 2 3 1 0], [1 0; 2 1; 0 2]);
            v = biconjQ(f).eval([-0.4088541666666667 1.1783854166666667]);
            tc.verifyGreaterThanOrEqual(v, 1.8073059721130218 - 1e-9);
            tc.verifyLessThanOrEqual(v, 1.8668169186760986);
        end

        function conjCPLQTinyPositiveEigenvalue(tc)
            % Finding 7 (legacy 'cplq' engine, MAJOR there). f = (x^2 + y^2)/2^31 is strictly convex,
            % f*(s) = 2^29 |s|^2 everywhere. conjCPLQ's absolute eig(Q) > sqrt(eps) test sends it
            % to the AFFINE branch and returns a needle: +inf at every s ~= 0.
            f = QuaPol([2^-30 0 2^-30 0 0 0]);
            g = conjCPLQ(f);
            v = g.eval([1 0; 0 1]);
            tc.verifyEqual(v(:).', [2^29 2^29], 'RelTol', 1e-9);
        end

        function conjCPLQSmallEigenvalueDoesNotCrash(tc)
            % Finding 7, second face. Q = diag(1, 2^-30) (strictly convex) and Q = diag(1, -2^-30)
            % (indefinite, f* = +inf everywhere) both fall into the rank-1 branch, where
            % null(Q) is EMPTY for a numerically full-rank Q: MATLAB:badsubscript, unnamed.
            g = conjCPLQ(QuaPol([1 0 2^-30 0 0 0]));
            v = g.eval([0 1]);
            tc.verifyEqual(v, 2^29, 'RelTol', 1e-9);
            g = conjCPLQ(QuaPol([1 0 -2^-30 0 0 0]));
            v = g.eval([0 1; 1 0]);
            tc.verifyEqual(v(:).', [Inf Inf]);
        end
    end
end
