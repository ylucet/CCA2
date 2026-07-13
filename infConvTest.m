classdef infConvTest < matlab.unittest.TestCase
    % Tests for infConv.m: (f box g) = conj(add(conj f,conj g)), valid for f,g both convex
    % (DESIGN.md II.6). Exercised here on full-domain strictly convex quadratics -- the one
    % end-to-end path currently closed under conjCPLQ (Step 3, max-of-conjugates over a
    % multi-face domain, is not implemented yet, so a bounded-triangle f,g pair cannot round-trip
    % all the way back through the final conj call).

    methods (Test)
        function isotropicEnergiesHalveTheCurvature(testCase)
            % f=g=0.5|x|^2 (default engine, nargin<3 branch): f box g = 0.25|x|^2, a classical
            % identity (ab/(a+b) with a=b=1 gives 1/2, halved again by the leading 0.5).
            f = QuaPoly([1 0 1 0 0 0]);
            g = QuaPoly([1 0 1 0 0 0]);
            h = infConv(f, g);
            testCase.verifyClass(h, 'QuaPoly');
            testCase.verifyEqual(h.nv, 0);
            testCase.verifyEqual(h.eval([2 -3]), 0.25*(2^2+3^2), 'AbsTol', 1e-10);
        end

        function anisotropicDiagonalMatchesClosedForm(testCase)
            % f(x)=x1^2+1.5x2^2 (A=diag(2,3)), g(x)=2x1^2+0.5x2^2 (B=diag(4,1)): for pure
            % quadratics f box g = 0.5 x'Cx with C = A(A+B)^-1 B = diag(2*4/6, 3*1/4).
            f = QuaPoly([2 0 3 0 0 0]);
            g = QuaPoly([4 0 1 0 0 0]);
            h = infConv(f, g, 'cplq');
            testCase.verifyClass(h, 'QuaPoly');
            c1 = (2*4)/(2+4); c2 = (3*1)/(3+1);
            testCase.verifyEqual(h.eval([1 1]), 0.5*(c1+c2), 'AbsTol', 1e-8);
            testCase.verifyEqual(h.eval([2 -1]), 0.5*(c1*4+c2*1), 'AbsTol', 1e-8);
        end

        function shiftedQuadraticsMatchIndependentMinimization(testCase)
            % General PD quadratics with linear/constant terms, cross-checked against the
            % independent closed-form minimizer z* = argmin_z f(z)+g(x0-z) (solved directly from
            % the stationarity condition, NOT via conj's own formula) rather than reusing the
            % code path under test.
            Af = diag([2 3]); pf = [1; -2]; af = 0.5;
            Ag = diag([1 4]); pg = [-1; 1]; ag = -0.3;
            f = QuaPoly([Af(1,1), Af(1,2), Af(2,2), pf(1), pf(2), af]);
            g = QuaPoly([Ag(1,1), Ag(1,2), Ag(2,2), pg(1), pg(2), ag]);
            h = infConv(f, g);

            x0 = [0.3; 0.7];
            zstar = (Af+Ag) \ (Ag*x0 + pg - pf);
            hx0 = 0.5*zstar'*Af*zstar + pf'*zstar + af ...
                + 0.5*(x0-zstar)'*Ag*(x0-zstar) + pg'*(x0-zstar) + ag;
            testCase.verifyEqual(h.eval(x0'), hx0, 'AbsTol', 1e-8);
        end
    end
end
