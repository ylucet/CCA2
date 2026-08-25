classdef conicMeetTest < matlab.unittest.TestCase
% Unit tests for conicMeet, the vertex layer of the sym-free conjugate.
%
% BUCKET: fast (closed-form; the only iteration is a bounded Newton polish).
%
% WHAT IS BEING PINNED. Every assertion below is a PROPERTY of the intersection -- the returned
% points satisfy both conics, their count is right, the canonical order is order-independent, the
% quartic certificate is exact and integral -- rather than a golden coordinate list. That matters
% here more than usual, because the coordinates are the one part of the design that is deliberately
% floating point: a golden-value test on them would pin today's root-finder rather than the
% mathematics, and would go red the day the exact degree-<=4 kernel replaces `roots`.

    methods (Test)

        % ---- the cases the arrangement actually meets -------------------------------------------

        function twoLinesMeetInOnePoint(testCase)
            L1 = [0 0 0 1 0 -1];        % x = 1
            L2 = [0 0 0 0 1 -2];        % y = 2
            P = conicMeet(L1, L2);
            testCase.verifySize(P, [1 2]);
            testCase.verifyEqual(P, [1 2], 'AbsTol', 1e-12);
        end

        function aLineMeetsACircleInTwoPointsAndATangentInOne(testCase)
            C = [1 0 1 0 0 -1];                     % x^2 + y^2 = 1
            P = conicMeet([0 0 0 0 1 0], C);        % y = 0
            testCase.verifySize(P, [2 2]);
            testCase.verifyEqual(sort(P(:,1)), [-1; 1], 'AbsTol', 1e-9);
            T = conicMeet([0 0 0 1 0 -1], C);       % x = 1, tangent
            testCase.verifySize(T, [1 2], 'a tangency is ONE point, not two coincident ones');
            testCase.verifyEqual(T, [1 0], 'AbsTol', 1e-8);
        end

        function twoCirclesMeetInTwoPoints(testCase)
            A = [1 0 1 0 0 -4];                     % x^2 + y^2 = 4
            B = [1 0 1 -6 0 5];                     % (x-3)^2 + y^2 = 4
            [P, info] = conicMeet(A, B);
            testCase.verifySize(P, [2 2]);
            testCase.verifyEqual(P(:,1), [1.5; 1.5], 'AbsTol', 1e-9);
            testCase.verifyEqual(sort(P(:,2)), [-sqrt(7)/2; sqrt(7)/2], 'AbsTol', 1e-9);
            testCase.verifyLessThan(max(info.resid), 1e-12);
        end

        function anEllipseAndAHyperbolaCanMeetInFourPoints(testCase)
        % Four is the maximum and it is reachable, which is why the vertex degree bound is 4.
            A = [1 0 4 0 0 -4];                     % x^2 + 4y^2 = 4
            B = [1 0 -1 0 0 0];                     % x^2 = y^2
            P = conicMeet(A, B);
            testCase.verifySize(P, [4 2]);
            testCase.verifyLessThan(max(abs(conicMeetTest.evalAt(A, P))), 1e-9);
            testCase.verifyLessThan(max(abs(conicMeetTest.evalAt(B, P))), 1e-9);
        end

        function conicsThatDoNotMeetReturnNothing(testCase)
            A = [1 0 1 0 0 -1];                     % unit circle
            B = [1 0 1 -20 0 99];                   % (x-10)^2 + y^2 = 1
            testCase.verifyEmpty(conicMeet(A, B));
        end

        % ---- the certificate is exact, which is the whole point of the resultant route -----------

        function theQuarticCertificateIsIntegralAndVanishesAtEveryReturnedPoint(testCase)
            A = [3 1 2 -4 5 -6];
            B = [1 -2 5 7 -1 -9];
            [P, info] = conicMeet(A, B);
            testCase.verifyEqual(info.quartic, round(info.quartic), ...
                'the resultant must be an EXACT integer polynomial -- it is the vertex name');
            testCase.verifyLessThan(max(abs(info.quartic)), ratQ.LIMIT);
            for i = 1:size(P,1)
                xs = P(i,1) - info.shear * P(i,2);      % the certificate lives in sheared x
                testCase.verifyLessThan(abs(polyval(info.quartic, xs)) / ...
                    max(1, max(abs(info.quartic)) * max(1, abs(xs))^4), 1e-8);
            end
        end

        function theShearIsUsedExactlyWhenAConicHasNoYSquaredTerm(testCase)
        % b x y + ... has no y^2 term, so the 4x4 Sylvester determinant is not the resultant until
        % the shear supplies one. Without it this pair returns the wrong answer silently, which is
        % why the shear is not an optimisation.
            A = [0 1 0 0 0 -1];                     % xy = 1
            B = [0 1 0 0 0 -4];                     % xy = 4  (disjoint)
            [P, info] = conicMeet(A, B);
            testCase.verifyGreaterThan(info.shear, 0);
            testCase.verifyEmpty(P);
            C = [1 0 1 0 0 -5];                     % x^2 + y^2 = 5 meets xy = 1 in four points
            [P2, info2] = conicMeet(A, C);
            testCase.verifyGreaterThan(info2.shear, 0);
            testCase.verifySize(P2, [4 2]);
            testCase.verifyLessThan(max(abs(conicMeetTest.evalAt(A, P2))), 1e-8);
            testCase.verifyLessThan(max(abs(conicMeetTest.evalAt(C, P2))), 1e-8);
        end

        % ---- the canonical order is a function of the two conics and nothing else ---------------

        function theOrderDoesNotDependOnTheArgumentOrderOrOnTheConicScaling(testCase)
        % The order IS the vertex name's root index, so this is a correctness property, not a
        % cosmetic one: the same geometric vertex must get the same index however it is reached.
            A = [1 0 4 0 0 -4];
            B = [1 0 -1 0 0 0];
            P1 = conicMeet(A, B);
            P2 = conicMeet(B, A);
            P3 = conicMeet(A * 7, B * -3);
            testCase.verifyEqual(P2, P1, 'AbsTol', 1e-8, 'swapping the arguments must not reorder');
            testCase.verifyEqual(P3, P1, 'AbsTol', 1e-8, 'rescaling a conic must not reorder');
        end

        % ---- degeneracy is reported, never resolved ---------------------------------------------

        function twoConicsSharingAComponentAreReportedDegenerateRatherThanIntersected(testCase)
        % Infinitely many intersections is not a vertex. Returning four coincident points here is
        % the failure mode that would put a phantom vertex into the mesh.
            [P, info] = conicMeet([1 0 1 0 0 -1], [1 0 1 0 0 -1]);
            testCase.verifyTrue(info.degenerate);
            testCase.verifyEmpty(P);
            % xy = 0 and x^2 y = 0 are not both conics, so the reachable case is a repeated line
            % pair: {x^2 = 0} shares the line x = 0 with {x^2 + xy = 0}.
            [P2, info2] = conicMeet([1 0 0 0 0 0], [1 1 0 0 0 0]);
            testCase.verifyTrue(info2.degenerate || size(P2,1) <= 1, ...
                'a shared component must not be reported as a transversal vertex');
        end

        function aNearTangencyIsReportedAsTwoPointsWithASmallSep(testCase)
        % The filter needs to SEE a near-degeneracy: two intersections a hair apart must come back
        % as two points with a small .sep, not merged into one.
            A = [1 0 1 0 0 -1];                          % x^2 + y^2 = 1
            B = [10000 0 10000 -39800 0 29601];          % (x - 199/100)^2 + y^2 = 1, near-tangent
            [P, info] = conicMeet(A, B);
            testCase.verifySize(P, [2 2]);
            testCase.verifyLessThan(info.sep, 0.25);
            testCase.verifyGreaterThan(info.sep, 0);
            testCase.verifyLessThan(max(info.resid), 1e-9, ...
                'both points must still lie on both conics');
        end

        function aTighterNearTangencyOVERFLOWSLOUDLYRatherThanReturningAWrongAnswer(testCase)
        % THE PRECISION BUDGET, PINNED. Expressing a 1e-4 near-tangency EXACTLY needs conic
        % coefficients of order 1e8, and the resultant is degree 2 in each conic's coefficients, so
        % the quartic wants ~1e17 -- past 2^53. Measured: 2.3995600239996e17.
        %
        % This is the designed outcome, not a defect to route around. The alternative to raising is
        % a silently rounded resultant, whose roots are a plausible wrong answer, and that is the
        % one thing ratQ exists to prevent. What it bounds is how NEAR a degeneracy this layer can
        % express exactly -- and in the pipeline the conic coefficients come from differences of
        % face functions, whose size is set by the input data rather than by an arbitrarily fine
        % perturbation. If it ever fires in real use the answer is a wider integer type.
            A = [1 0 1 0 0 -1];
            B = [1e8 0 1e8 -2*19999*1e4 0 19999^2 - 1e8];   % (x - 19999/10000)^2 + y^2 = 1
            testCase.verifyError(@() conicMeet(A, B), 'ratQ:overflow');
        end

        % ---- a differential check against the definition -----------------------------------------

        function everyReturnedPointSatisfiesBothConicsOverARandomSweep(testCase)
        % The oracle is the DEFINITION -- a point is an intersection iff both conics vanish there --
        % rather than another implementation. Integer conics keep the certificate exact, so a
        % failure here is a real defect and not a conditioning artefact.
            rng(20260824);
            found = 0;
            for k = 1:300
                A = randi([-6 6], 1, 6);  B = randi([-6 6], 1, 6);
                if all(A == 0) || all(B == 0), continue, end
                try
                    [P, info] = conicMeet(A, B);
                catch ME
                    testCase.verifyEqual(ME.identifier, 'ratQ:overflow', ...
                        sprintf('case %d raised %s, which is not an accepted refusal', k, ME.identifier));
                    continue
                end
                if info.degenerate, continue, end
                for i = 1:size(P,1)
                    found = found + 1;
                    testCase.verifyLessThan(info.resid(i), 1e-8, ...
                        sprintf('case %d point %d does not lie on both conics', k, i));
                end
            end
            testCase.verifyGreaterThan(found, 200, 'the sweep must actually exercise intersections');
        end
    end

    methods (Static)
        function v = evalAt(C, P)
            v = C(1)*P(:,1).^2 + C(2)*P(:,1).*P(:,2) + C(3)*P(:,2).^2 ...
                + C(4)*P(:,1) + C(5)*P(:,2) + C(6);
        end
    end
end
