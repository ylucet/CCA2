classdef QuaPolExactTest < matlab.unittest.TestCase
% Unit tests for QuaPol's EXACT storage -- fN/fD/VN/VD, isExact and assertExact.
%
% BUCKET: fast (integer arithmetic only; no symbolic call anywhere).
%
% WHAT IS BEING PINNED. CCA2's design target is the conjugate of a QuaPol computed exactly over the
% rationals, so the INPUT has to carry exact coefficients and exact domain vertices. A caller writes
% doubles because that is what MATLAB source looks like, so the constructor reconstructs the
% rational they meant. The legacy FLOAT pipeline also builds QuaPol objects, out of COMPUTED
% doubles, and those are genuinely inexact -- so the contract is: exact when the data is a caller's,
% honestly not-exact when it is not, and a named refusal if anyone tries to compute exactly on the
% second kind.

    methods (Test)

        function aCallerWrittenQuadraticIsStoredExactly(testCase)
            f = QuaPol([3 1 3 0 0 0]);
            testCase.verifyTrue(f.isExact());
            [n, d] = f.faceQ(1);
            testCase.verifyEqual(n, [0 0 0 0, 3, 1, 3, 0, 0, 0]);
            testCase.verifyEqual(d, 1);
        end

        function fractionsAreRecoveredAsTheRationalTheCallerMeant(testCase)
        % 1/3 is not representable in binary at all -- the double is 6004799503160661/2^54 -- so a
        % bit-exact reading would return a 2^54 denominator, which is useless and then overflows.
        % What is wanted, and what is uniquely determined below maxDen, is the rational written.
            f = QuaPol([1/2 1/3 1/6 0 0 0]);
            testCase.verifyTrue(f.isExact());
            [n, d] = f.faceQ(1);
            testCase.verifyEqual(d, 6, 'the common denominator is the lcm of 2, 3 and 6');
            testCase.verifyEqual(n, [0 0 0 0, 3, 2, 1, 0, 0, 0]);
        end

        function theExactFieldsAreTheSameNumberAsTheDoubleOnes(testCase)
        % When both are present they are one value in two spellings, and the binary one is the
        % rounded view of the exact one -- not an independent claim that could disagree.
            rng(20260903);
            for k = 1:20
                num = randi([-9 9], 1, 6);  den = randi([1 12]);
                f = QuaPol(num / den);
                testCase.verifyTrue(f.isExact());
                [n, d] = f.faceQ(1);
                testCase.verifyEqual(n(5:10) / d, f.f(5:10), 'AbsTol', 0, ...
                    sprintf('case %d: the exact field and the double field disagree', k));
            end
        end

        function aPolyhedralDomainStoresItsVerticesExactlyToo(testCase)
        % Unlike a CONJUGATE's vertices, a QuaPol's are the caller's own data and are rational --
        % CONJ_FIELD_PROOF.md Theorem 1 is a statement about f*, not about f. So they are simply
        % stored, and no naming machinery is needed on the primal side.
            f = QuaPol([0 0; 1 0; 1 1; 0 1], [1 2 1; 2 3 1; 3 4 1; 4 1 1], ...
                       [0 0 0 0 1 0 1 0 0 0], [1 0; 1 0; 1 0; 1 0]);
            testCase.verifyTrue(f.isExact());
            for i = 1:f.nv
                [n, d] = f.vertexQ(i);
                testCase.verifyEqual(n / d, f.V(i,:), 'AbsTol', 0);
            end

            g = QuaPol([0 0; 1/2 0; 1/3 2/7], [1 2 1; 2 3 1; 3 1 1], ...
                       [0 0 0 0 1 0 1 0 0 0], [1 0; 1 0; 1 0]);
            testCase.verifyTrue(g.isExact(), 'fractional vertices are still exact');
            [n, d] = g.vertexQ(3);
            testCase.verifyEqual(d, 21);  testCase.verifyEqual(n, [7 6]);
        end

        % ---- the honest unknown ----------------------------------------------------------------

        function anIrrationalCoefficientMakesTheObjectInexactWithoutRaising(testCase)
        % The constructor must NOT raise: the legacy float pipeline builds QuaPol objects at nine
        % production call sites and is kept for unit testing, so raising here would stop it on day
        % one. The object is simply not exact, and says so.
            f = QuaPol([sqrt(2) 0 1 0 0 0]);
            testCase.verifyFalse(f.isExact());
            testCase.verifyEmpty(f.fN);
            testCase.verifyEqual(f.f(5), sqrt(2), 'AbsTol', 0, ...
                'the double field is untouched, so the legacy pipeline still works');
        end

        function aValueOneUlpOffASimpleRationalIsRecognisedAsInexact(testCase)
        % This is not hypothetical: DECISIONS.md 2026-09-03 measured today's conj emitting exactly
        % this value where the truth is 1, on x^2-y^2 over the unit square. Sixteen of sixty-four
        % face rows were one or two ULP out.
            f = QuaPol([0.99999999999999978 0 1 0 0 0]);
            testCase.verifyFalse(f.isExact(), ...
                'one ULP below 1 is not the rational 1, and must not be snapped to it');
        end

        function exactnessIsAllOrNothingAcrossTheWholeObject(testCase)
        % A partially exact object is the trap: every caller would have to ask per row, and the one
        % that forgot would compute exactly on rounded data and produce a plausible wrong answer.
            V = [0 0; 1 0; 1 1; 0 1];  E = [1 2 1; 2 3 1; 3 4 1; 4 1 1];
            fGood = [0 0 0 0 1 0 1 0 0 0];
            fMixed = [fGood; 0 0 0 0 pi 0 1 0 0 0];
            f = QuaPol(V, E, fMixed, [1 2; 1 2; 1 2; 1 2]);
            testCase.verifyFalse(f.isExact(), 'one inexact face makes the whole object inexact');
            testCase.verifyEmpty(f.fD);
        end

        function computingExactlyOnInexactDataIsRefusedByName(testCase)
        % The guard that keeps the two pipelines from contaminating each other. Exact arithmetic on
        % coefficients that are already one ULP wrong produces exactly the WRONG number, which
        % carries no warning -- pieceRecessionRays did that and it was a real defect (G10).
            f = QuaPol([sqrt(2) 0 1 0 0 0]);
            testCase.verifyError(@() f.assertExact(), 'PLQ:QuaPol:notExact');
            testCase.verifyError(@() f.faceQ(1),      'PLQ:QuaPol:notExact');
            testCase.verifyError(@() f.vertexQ(1),    'PLQ:QuaPol:notExact');
            testCase.verifyError(@() conjQ(f),        'PLQ:QuaPol:notExact');
        end

        % ---- the aliases must not have been left behind ------------------------------------------

        function theBackwardCompatibleAliasesCarryTheExactFieldsToo(testCase)
        % QuaPoly is what the SCIP bridge constructs by name, and PLQVC is the released name. Both
        % are constructor forwards, so both must inherit exactness rather than silently not.
            testCase.verifyTrue(QuaPoly([3 1 3 0 0 0]).isExact());
            testCase.verifyTrue(PLQVC([3 1 3 0 0 0]).isExact());
            testCase.verifyEqual(QuaPoly([1/2 0 1/2 0 0 0]).fD, 2);
        end

        function anEmptyQuaPolIsNotClaimedToBeExact(testCase)
        % The no-argument path writes nothing (RatCon's constructor protocol), so there is nothing
        % exact about it and it must not pretend otherwise.
            testCase.verifyFalse(QuaPol().isExact());
        end
    end
end
