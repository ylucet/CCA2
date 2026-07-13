classdef lasryLionsTest < matlab.unittest.TestCase
    % Tests for lasryLions.m: h = -e_mu(-e_lambda(f)), a pure composition of moreau (no new
    % geometry). Cross-checked against the componentwise closed form for a diagonal-quadratic
    % input, computed independently in each test (not by calling moreau/lasryLions again).

    methods (Test)
        function isotropicQuadraticMatchesNestedClosedForm(testCase)
            % f(x)=0.5*a*|x|^2 (default engine, nargin<4): chain moreau's own closed-form
            % coefficient recursion c -> c/(1+c*t) through lambda, negate, mu, negate again.
            a = 2; lambda = 0.3; mu = 0.4;
            f = QuaPoly([a 0 a 0 0 0]);
            h = lasryLions(f, lambda, mu);
            testCase.verifyClass(h, 'QuaPoly');

            step = @(c,t) c/(1+c*t);
            c1 = step(a, lambda);        % = e_lambda f's own coefficient
            c2 = -step(-c1, mu);         % negate, moreau(mu), negate again
            testCase.verifyEqual(h.eval([3 -2]), 0.5*c2*(3^2+2^2), 'AbsTol', 1e-8);
        end

        function anisotropicMixedConvexityMatchesNestedClosedForm(testCase)
            % A=diag(2,-1): component 1 convex, component 2 CONCAVE -- lasryLions needs no
            % convexity on f (inherited from moreau, see moreau.m), only that each nested moreau
            % call's own positive-definiteness condition holds for the chosen lambda,mu.
            a1 = 2; a2 = -1; lambda = 0.2; mu = 0.3;
            f = QuaPoly([a1 0 a2 0 0 0]);
            h = lasryLions(f, lambda, mu, 'cplq');
            testCase.verifyClass(h, 'QuaPoly');

            step = @(c,t) c/(1+c*t);
            nested = @(a) -step(-step(a,lambda), mu);
            c1 = nested(a1); c2 = nested(a2);
            testCase.verifyEqual(h.eval([1 2]), 0.5*c1*1^2 + 0.5*c2*2^2, 'AbsTol', 1e-8);
        end
    end
end
