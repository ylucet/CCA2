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
