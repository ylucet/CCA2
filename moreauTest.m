classdef moreauTest < matlab.unittest.TestCase
    % Tests for moreau.m (Moreau envelope via a single conjugate, [HIRIART-URRUTY-07]) and its
    % addQuadratic/addScaledEnergy prerequisites (QuaPoly/QuaPar per-face coefficient update).

    methods (Test)
        function isotropicEnergyMatchesClosedForm(testCase)
            % f(x)=0.5|x|^2 (a=1): e_mu f(x) = 0.5*a/(1+a*mu)*|x|^2 (default engine, nargin<3).
            f = QuaPoly([1 0 1 0 0 0]);
            mu = 2;
            h = moreau(f, mu);
            testCase.verifyClass(h, 'QuaPoly');
            coeff = 0.5/(1+mu);   % = 0.5*a/(1+a*mu) with a=1
            testCase.verifyEqual(h.eval([3 4]), coeff*(3^2+4^2), 'AbsTol', 1e-8);
        end

        function anisotropicDiagonalMatchesClosedForm(testCase)
            % f(x)=x1^2+1.5x2^2 (A=diag(2,3)): componentwise e_mu f_i = 0.5*a_i/(1+a_i*mu)*xi^2.
            A = diag([2 3]); mu = 0.5;
            f = QuaPoly([A(1,1) A(1,2) A(2,2) 0 0 0]);
            h = moreau(f, mu, 'cplq');
            c1 = 0.5*A(1,1)/(1+A(1,1)*mu); c2 = 0.5*A(2,2)/(1+A(2,2)*mu);
            testCase.verifyEqual(h.eval([2 -1]), c1*2^2 + c2*(-1)^2, 'AbsTol', 1e-8);
        end

        function nonconvexInputNeedsNoConvexity(testCase)
            % f(x)=0.5*(x1^2-x2^2) is genuinely NONCONVEX (A=diag(1,-1)): moreau's whole point
            % ([HIRIART-URRUTY-07]) is that this is still well-defined (single conj call, no
            % biconjugation) as long as mu*A+I stays positive definite (mu<1 here). Same
            % componentwise closed form as above, now with a negative eigenvalue.
            A = diag([1 -1]); mu = 0.5;
            f = QuaPoly([A(1,1) A(1,2) A(2,2) 0 0 0]);
            h = moreau(f, mu);
            c1 = 0.5*A(1,1)/(1+A(1,1)*mu); c2 = 0.5*A(2,2)/(1+A(2,2)*mu);
            testCase.verifyEqual(h.eval([2 3]), c1*2^2 + c2*3^2, 'AbsTol', 1e-8);
        end

        function shiftedConvexQuadraticMatchesIndependentMinimization(testCase)
            % General PD quadratic with linear/constant terms, cross-checked against the
            % independent closed-form minimizer z* = argmin_z f(z)+(1/(2mu))|x0-z|^2 (solved
            % directly from the stationarity condition, NOT via moreau's own conj call).
            A = diag([2 3]); p = [0.5; -1]; r = 0.2; mu = 0.3;
            f = QuaPoly([A(1,1), A(1,2), A(2,2), p(1), p(2), r]);
            h = moreau(f, mu);

            x0 = [1; -0.5];
            zstar = (A + eye(2)/mu) \ (x0/mu - p);
            hx0 = 0.5*zstar'*A*zstar + p'*zstar + r + (1/(2*mu))*sum((x0-zstar).^2);
            testCase.verifyEqual(h.eval(x0'), hx0, 'AbsTol', 1e-8);
        end
    end
end
