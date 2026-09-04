classdef conjQTest < matlab.unittest.TestCase
% Unit tests for conjQ, the EXACT conjugate: rational coefficients throughout, no tolerance in any
% decision, returning a QuaCon.
%
% BUCKET: fast (closed-form integer arithmetic; no symbolic call anywhere).
%
% The oracle is twofold, and deliberately so:
%   * the DEFINITION -- f*(s) = sup_x <s,x> - f(x) -- recomputed independently, which is what
%     CLAUDE.md section 7 asks for. For a strictly convex quadratic the sup is attained at the
%     stationary point, so the definition is available in closed form without a search.
%   * the LEGACY double pipeline (conjCPLQ's Case A), which is the thing being replaced judging the
%     replacement, CLAUDE.md section 6.

    methods (Test)

        function theExactConjugateAgreesWithTheDefinitionItself(testCase)
        % f*(s) = <s,x*> - f(x*) at the stationary point x* = H^-1 (s - L). Computed here from H
        % and L directly, in doubles, with no reference to conjQ's own algebra -- so agreement is
        % evidence about the mathematics and not about a shared expression.
            rng(20260903);
            for k = 1:40
                H = randi([-4 4], 2, 2);  H = H + H';       % symmetric integer
                if ~(H(1,1) > 0 && det(H) > 0.5), continue, end   % integer det: exact test
                L = randi([-5 5], 2, 1);
                kappa = randi([-5 5]);
                f = QuaPol([H(1,1), H(1,2), H(2,2), L(1), L(2), kappa]);

                g = conjQ(f);
                testCase.verifyEqual(g.kind(), 'QuaCon');

                S = [0 0; 1 0; 0 1; -2 3; 5 -4; randn(10,2)];
                got = g.eval(S);
                want = zeros(size(S,1),1);
                for i = 1:size(S,1)
                    s  = S(i,:)';
                    xs = H \ (s - L);                        % the maximiser
                    want(i) = s' * xs - (0.5 * xs' * H * xs + L' * xs + kappa);
                end
                testCase.verifyEqual(got, want, 'RelTol', 1e-9, 'AbsTol', 1e-9, ...
                    sprintf('case %d: conjQ disagrees with the definition of the sup', k));
            end
        end

        function theExactConjugateAgreesWithTheLegacyDoublePipeline(testCase)
        % Differential against conjCPLQ Case A. The legacy route returns a QuaPol of doubles; this
        % one returns a QuaCon of exact rationals. They must be the same function.
            cases = { [2 0 2 0 0 0], [1 0 1 3 -2 5], [4 1 3 0 0 -7], [1 0 9 -1 2 0] };
            for k = 1:numel(cases)
                f  = QuaPol(cases{k});
                ge = conjQ(f);
                gl = f.conj();
                S  = [0 0; 1 1; -3 2; 7 -5; 0.25 0.75];
                testCase.verifyEqual(ge.eval(S), gl.eval(S), 'RelTol', 1e-12, ...
                    sprintf('case %d: exact and legacy conjugates differ', k));
            end
        end

        function theConjugateCoefficientsAreExactRationalsNotRoundedDoubles(testCase)
        % f = (3x^2 + 2xy + 3y^2)/2 has H = [3 1; 1 3], det 8, so f* carries thirds and eighths --
        % values a double cannot hold exactly. The point of the class is that the STORED
        % coefficients are the exact rationals, whatever the double evaluation later rounds to.
            f = QuaPol([3 1 3 0 0 0]);
            g = conjQ(f);
            % H^-1 = [3 -1; -1 3]/8, so f*(s) = (3 s1^2 - 2 s1 s2 + 3 s2^2)/16 and the weighted
            % slots are c5 = 3/8, c6 = -1/8, c7 = 3/8.
            testCase.verifyEqual(g.fD, 8);
            testCase.verifyEqual(g.fN, [0 0 0 0, 3, -1, 3, 0, 0, 0]);
        end

        function theExactTestAcceptsAStrictlyConvexQuadraticEigWouldRefuse(testCase)
        % THE reason the exact test is worth having. Q = [1 16384; 16384 268435457] has determinant
        % exactly 1 and leading minor 1, so it is positive definite -- but its smaller eigenvalue is
        % about 3.7e-9, below conjCPLQ's fixed sqrt(eps) threshold, so the legacy route does not
        % take its strictly-convex branch. A fixed absolute tolerance on a scale-dependent quantity
        % is not a convexity test.
            k = 16384;  N = k^2 + 1;
            f = QuaPol([1, k, N, 0, 0, 0]);

            H = [1 k; k N];
            testCase.verifyEqual(det(H), 1, 'AbsTol', 0.5, 'the determinant really is 1');
            testCase.verifyLessThan(min(eig(H)), sqrt(eps), ...
                'and eig really does put the smaller eigenvalue below the legacy threshold');

            g = conjQ(f);                                  % exact: accepted
            testCase.verifyEqual(g.fN, [0 0 0 0, N, -k, 1, 0, 0, 0], ...
                'H^-1 = adj(H)/det = adj(H), exactly');
            testCase.verifyEqual(g.fD, 1);

            % and it is RIGHT, against the definition rather than against itself
            S = [1 0; 0 1; 1 1; 2 -3];
            want = zeros(size(S,1),1);
            for i = 1:size(S,1)
                s = S(i,:)';  xs = H \ s;
                want(i) = s' * xs - 0.5 * xs' * H * xs;
            end
            testCase.verifyEqual(g.eval(S), want, 'RelTol', 1e-9);

            % NOT asserted here, deliberately: what the LEGACY route does with this input. It does
            % not refuse -- it returns a QuaPar that evaluates to +Inf at all four points above,
            % i.e. it reports dom f* as empty when f* is finite everywhere. That is a silent wrong
            % answer in the path this one replaces, recorded with its reproducer in DECISIONS.md
            % and TODO.md 2026-09-03. Pinning it here would be a golden-value test of a defect,
            % which would then have to be edited when the defect goes.
        end

        function aNonStrictlyConvexFullDomainQuadraticIsRefusedByName(testCase)
        % Exact-or-named-refusal, never a silent fallback.
            testCase.verifyError(@() conjQ(QuaPol([1 1 1 0 0 0])), ...
                'PLQ:conjQ:notStrictlyConvex', 'a rank-one PSD form is not strictly convex');
            testCase.verifyError(@() conjQ(QuaPol([1 0 -1 0 0 0])), ...
                'PLQ:conjQ:notStrictlyConvex', 'an indefinite form is not strictly convex');
            testCase.verifyError(@() conjQ(QuaPol([-1 0 -1 0 0 0])), ...
                'PLQ:conjQ:notStrictlyConvex', 'a concave form is not strictly convex');
        end

        function aCaseNotYetImplementedIsRefusedByNameAndSaysWhy(testCase)
        % The bounded-domain cases are the next work item, and until they land the refusal must be
        % nameable rather than a fallback into the symbolic engine.
            f = QuaPol([0 0; 1 0; 1 1; 0 1], [1 2 1; 2 3 1; 3 4 1; 4 1 1], ...
                       [0 0 0 0 1 0 1 0 0 0], [1 0; 1 0; 1 0; 1 0]);
            testCase.verifyError(@() conjQ(f), 'PLQ:conjQ:notImplemented');
        end

        function anInputThatIsNotExactlyRationalIsRefusedRatherThanSnapped(testCase)
        % Bounding the vertex denominators does not bound the downstream ones -- DECISIONS.md's
        % attempt 3 carried 1e5 to 1e25 in a few squarings and the run hung. So an irrational
        % coefficient is a defect in the caller, not a case to round away.
            f = QuaPol([sqrt(2) 0 1 0 0 0]);
            testCase.verifyFalse(f.isExact(), ...
                'an irrational coefficient makes the whole object inexact, all or nothing');
            testCase.verifyError(@() conjQ(f), 'PLQ:QuaPol:notExact');
        end

        function theConjugateOfAConjugateReturnsTheOriginalFunction(testCase)
        % f** = f for a strictly convex quadratic on all of R^2. A round trip is a definition test
        % that no single-direction check can substitute for: it catches a consistent sign or
        % transpose error that both directions would otherwise share.
            rng(20260903);
            for k = 1:15
                H = randi([-4 4], 2, 2);  H = H + H';
                if ~(H(1,1) > 0 && det(H) > 0.5), continue, end
                L = randi([-5 5], 2, 1);  kappa = randi([-5 5]);
                fN0 = [0 0 0 0, H(1,1), H(1,2), H(2,2), L(1), L(2), kappa];

                g  = conjQ(QuaPol([H(1,1), H(1,2), H(2,2), L(1), L(2), kappa]));
                % conjQ takes a QuaPol, so read g's exact face back through the same 6-slot form.
                gq = QuaPol([g.fN(5)/g.fD, g.fN(6)/g.fD, g.fN(7)/g.fD, ...
                             g.fN(8)/g.fD, g.fN(9)/g.fD, g.fN(10)/g.fD]);
                h  = conjQ(gq);

                [n0, d0] = ratQ.canon(fN0, 1);
                testCase.verifyEqual(h.fN, n0, sprintf('case %d: f** numerator', k));
                testCase.verifyEqual(h.fD, d0, sprintf('case %d: f** denominator', k));
            end
        end
    end
end
