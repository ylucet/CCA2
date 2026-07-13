classdef proxAverageTest < matlab.unittest.TestCase
    % Tests for proxAverage.m: P = P_mu(f,g;lambda,1-lambda). Cross-checked against P's OWN
    % defining characterization -- e_mu*P = lambda*e_mu*f + (1-lambda)*e_mu*g -- via the
    % already-tested moreau.m, rather than a hand-derived closed form for P itself (avoids
    % re-deriving proxAverage's algebra by hand; exercises the actual mathematical property that
    % motivates the formula in DESIGN.md II.6).

    methods (Test)
        function isotropicQuadraticsSatisfyMoreauEnvelopeCharacterization(testCase)
            f = QuaPoly([2 0 2 0 0 0]);   % 0.5*2*|x|^2
            g = QuaPoly([3 0 3 0 0 0]);   % 0.5*3*|x|^2
            lambda = 0.4; mu = 0.5;
            h = proxAverage(f, g, lambda, mu);
            testCase.verifyClass(h, 'QuaPoly');

            pts = [1 2; -3 0.5];
            lhs = moreau(h, mu).eval(pts);
            rhs = lambda*moreau(f, mu).eval(pts) + (1-lambda)*moreau(g, mu).eval(pts);
            testCase.verifyEqual(lhs, rhs, 'AbsTol', 1e-8);
        end

        function anisotropicDiagonalQuadraticsSatisfyMoreauCharacterization(testCase)
            Af = diag([1 4]); Ag = diag([3 0.5]);
            f = QuaPoly([Af(1,1) 0 Af(2,2) 0 0 0]);
            g = QuaPoly([Ag(1,1) 0 Ag(2,2) 0 0 0]);
            lambda = 0.7; mu = 0.3;
            h = proxAverage(f, g, lambda, mu, 'cplq');
            testCase.verifyClass(h, 'QuaPoly');

            pts = [1 -1; 2 3];
            lhs = moreau(h, mu).eval(pts);
            rhs = lambda*moreau(f, mu).eval(pts) + (1-lambda)*moreau(g, mu).eval(pts);
            testCase.verifyEqual(lhs, rhs, 'AbsTol', 1e-8);
        end
    end
end
