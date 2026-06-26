classdef conjCPLQTest < matlab.unittest.TestCase
    % Tests for the 'cplq' conjugate engine (conjCPLQ) and the QuaPoly operator interface.
    % Covers the currently implemented case (full-domain quadratics) plus the
    % validation/rejection behaviour. See DESIGN.md II.5.1.

    methods (Test)
        function energyIsSelfConjugate(testCase)
            % f = 1/2(x^2+y^2) is self-conjugate.
            p = QuaPoly.energy();
            q = p.conj('cplq');
            testCase.verifyEqual(q.f, p.f, 'AbsTol', 1e-12);
        end

        function defaultEngineIsCPLQ(testCase)
            % conj() with no engine argument defaults to 'cplq'.
            p = QuaPoly.energy();
            testCase.verifyEqual(p.conj().f, p.conj('cplq').f, 'AbsTol', 1e-12);
        end

        function biconjOfConvexQuadraticIsItself(testCase)
            % f** = f for a closed convex function.
            p = QuaPoly.energy();
            b = p.biconj('cplq');
            testCase.verifyEqual(b.f, p.f, 'AbsTol', 1e-12);
        end

        function generalPositiveDefiniteQuadratic(testCase)
            % f(x) = 1/2 x'Q x + L'x + k,  Q=[2 0;0 4], L=[1;-1], k=3.
            % f*(s) = 1/2 (s-L)' inv(Q) (s-L) - k.
            f6 = [2 0 4 1 -1 3];          % [x^2 xy y^2 x y const]
            p  = QuaPoly(f6);
            q  = p.conj('cplq');
            Q = [2 0; 0 4]; L = [1; -1]; k = 3;
            M = inv(Q); grad = -M*L; d = 0.5*(L'*M*L) - k;
            expf6 = [M(1,1) M(1,2) M(2,2) grad(1) grad(2) d];
            expq  = QuaPoly(expf6);
            % Compare both the coefficients and the values on sample dual points.
            testCase.verifyEqual(q.f, expq.f, 'AbsTol', 1e-12);
            S = [0 0; 1 1; -2 3; 0.5 -1; 4 -4];
            testCase.verifyEqual(q.eval(S), expq.eval(S), 'AbsTol', 1e-10);
        end

        function fenchelYoungAtMinimizer(testCase)
            % Sanity: for f above, f*(s) at s=L (gradient at the minimizer x=0... here x*=Q\(s-L))
            % equals -k when s=L, since the minimizer of f is at x=-Q\L and f*(L)= -inf? -- instead
            % check the duality value at a chosen point against the direct sup over a fine grid.
            f6 = [2 0 4 1 -1 3];
            p  = QuaPoly(f6);
            q  = p.conj('cplq');
            s  = [3; 2];
            % direct: sup_x <s,x> - f(x), maximizer x* = Q\(s-L)
            Q = [2 0; 0 4]; L = [1; -1]; k = 3;
            xstar = Q \ (s - L);
            fval  = 0.5*xstar'*Q*xstar + L'*xstar + k;
            direct = s'*xstar - fval;
            testCase.verifyEqual(q.eval(s'), direct, 'AbsTol', 1e-10);
        end

        function indefiniteQuadraticNotImplemented(testCase)
            % f = xy is indefinite; its conjugate is not a full-domain quadratic (QuaPar, TODO).
            p = QuaPoly([0 1 0 0 0 0]);
            testCase.verifyError(@() p.conj('cplq'), 'PLQ:conjCPLQ:notImplemented');
        end

        function cubicRejectedByOperators(testCase)
            % Cubic numerator is storable but rejected by operators (allowed for isConvex only).
            p = QuaPoly([1 0 0 0 0 0 0 0 0 0]);   % x^3/6 term present -> degree 3
            testCase.verifyEqual(p.degree, 3);
            testCase.verifyError(@() p.conj('cplq'), 'PLQ:op:unsupportedType');
        end

        function unimplementedEnginesError(testCase)
            p = QuaPoly.energy();
            testCase.verifyError(@() p.conj('pqp'),   'PLQ:conj:engine');
            testCase.verifyError(@() p.conj('graph'), 'PLQ:conj:engine');
        end

        function partialConjEngineRestriction(testCase)
            p = QuaPoly.energy();
            testCase.verifyError(@() p.partialConj(1,'graph'), 'PLQ:partialConj:engine');
        end

        function plqvcAliasStillWorks(testCase)
            % Backward compatibility: PLQVC is an alias of QuaPoly.
            p = PLQVC.energy();                 % inherited static factory
            testCase.verifyTrue(isa(p, 'QuaPoly'));
            V = [0 0;-1 0; 0 1;1 0;0 -1];
            E = [1 2 0;1 3 0;1 4 0;1 5 0];
            f = [1 0 1 0 0 0;1 0 2 0 0 0;2 0 2 0 0 0;2 0 1 0 0 0];
            F = [1 2;2 3;3 4;4 1];
            p2 = PLQVC(V,E,f,F);                % construct via the alias
            testCase.verifyTrue(isa(p2, 'PLQVC'));
            testCase.verifyEqual(p2.nf, 4);
        end
    end
end
