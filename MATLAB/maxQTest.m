classdef maxQTest < matlab.unittest.TestCase
% Unit tests for maxQ, the EXACT pointwise maximum of two QuaCons -- Step 3 of [COAP]/[JOGO].
%
% BUCKET: fast (integer arithmetic and closed-form conic intersection; no symbolic call).
%
% THE ORACLE IS THE DEFINITION AND IT IS UNUSUALLY GOOD HERE. max(g1,g2) evaluated at a point is
% just max(g1(s), g2(s)), which needs none of maxQ's machinery -- no overlay, no difference conic,
% no splitting. So the test compares the folded object against the two operands it was built from,
% and any error in the overlay, the sign of the split, or the cell assembly shows up immediately.
%
% MAXQ NEVER COMPUTES A NEW FUNCTION VALUE -- every cell of the result carries one of the two
% operands' face functions -- and the right way to assert that is on the STORED RATIONALS, not on
% the evaluation. `ratQ.canon` may reduce a face to a different spelling of the same number than
% the operand carried (dividing numerator and denominator by an odd gcd), and evaluating two equal
% spellings in double precision can differ in the last bit: evalPoly accumulates the rounding
% differently. Measured, not assumed -- 2 of 403 points at 4.4e-16, with the correct winner both
% times. So exactness is checked as exact rational EQUALITY of the face functions, and the
% evaluation comparison carries a few-ULP tolerance.

    methods (Test)

        function foldingTwoConjugatesReproducesTheirPointwiseMaximumExactly(testCase)
            [g1, g2] = maxQTest.twoTriangleConjugates();
            h = maxQ(g1, g2);
            rng(20260903);
            S = [randn(400,2)*2; 0 0; 1 1; -3 2];
            [got, idx] = h.eval(S);
            testCase.verifyTrue(all(idx > 0), 'the fold must cover everything the operands did');
            testCase.verifyEqual(got, max(g1.eval(S), g2.eval(S)), 'RelTol', 1e-14, ...
                'AbsTol', 1e-14, 'the fold must agree with the pointwise max of its operands');

            % The exactness claim itself, on the stored rationals rather than on the evaluation:
            % every face of the fold must BE one of the operands' faces, as an exact rational.
            for k = 1:h.nf
                found = false;
                for src = {g1, g2}
                    gg = src{1};
                    for i = 1:gg.nf
                        if ratQ.eqRat(h.fN(k,:), h.fD(k), gg.fN(i,:), gg.fD(i))
                            found = true;  break
                        end
                    end
                    if found, break, end
                end
                testCase.verifyTrue(found, ...
                    sprintf('face %d of the fold is not exactly a face of either operand', k));
            end
        end

        function theFoldIsCommutative(testCase)
        % max(a,b) = max(b,a) is a property no single evaluation pins, and the overlay loop is
        % asymmetric -- it iterates g1's faces outermost and attaches the split conic in one
        % order -- so this is a genuine check on the sign handling rather than a tautology.
            [g1, g2] = maxQTest.twoTriangleConjugates();
            rng(20260903);
            S = [randn(300,2)*2; 0 0];
            testCase.verifyEqual(maxQ(g1,g2).eval(S), maxQ(g2,g1).eval(S), 'AbsTol', 0);
        end

        function foldingAFunctionWithItselfChangesNothing(testCase)
        % max(g,g) = g. The difference is identically zero on every overlay cell, which is the one
        % branch that must NOT split -- splitting would produce two cells carrying one function
        % separated by a boundary that is not there.
            [g1, ~] = maxQTest.twoTriangleConjugates();
            h = maxQ(g1, g1);
            rng(20260903);
            S = [randn(200,2)*2; 0 0];
            testCase.verifyEqual(h.eval(S), g1.eval(S), 'AbsTol', 0);
            testCase.verifyEqual(h.nf, g1.nf, ...
                'the identical-function branch must not split a cell in two');
        end

        function theFoldOfTwoConvexFunctionsIsConvex(testCase)
        % A max of convex functions is convex -- a definition-level property that a wrong split
        % side breaks immediately, and that no cell-by-cell check implies.
            [g1, g2] = maxQTest.twoTriangleConjugates();
            h = maxQ(g1, g2);
            rng(20260903);
            A = randn(300,2)*2;  B = randn(300,2)*2;
            testCase.verifyLessThanOrEqual(h.eval((A+B)/2), (h.eval(A) + h.eval(B))/2 + 1e-12);
        end

        function theSplitBoundaryIsAConicAndItsSidesCarryDifferentFunctions(testCase)
        % WHY QuaCon EXISTS. Where two faces with different Hessians meet, the boundary lies on
        % {g_i = g_j}, whose quadratic part is (H_i - H_j)/2 -- a genuine conic, not a line, and
        % not storable as a QuaPar. This asserts that maxQ actually produces one.
            [g1, g2] = maxQTest.twoTriangleConjugates();
            h = maxQ(g1, g2);
            kinds = arrayfun(@(j) string(h.edgeKind(j)), 1:h.ne);
            testCase.verifyTrue(any(kinds ~= "line"), ...
                'folding faces with different Hessians must produce at least one curved edge');
            for j = 1:h.ne
                testCase.verifyEqual(h.EcQ(j,:), ratQ.conic(h.EcQ(j,:)), ...
                    sprintf('edge %d is not in canonical form', j));
            end
        end

        function everyStoredCurveAppearsExactlyOnceHoweverManyCellsUseIt(testCase)
        % The property the integer representation buys, asserted directly: a facet reached by two
        % different cells is ONE row, so the curve list has no duplicates. In doubles the two
        % spellings differ by an ULP and the shared facet becomes invisible to merge -- measured,
        % DECISIONS.md 2026-08-17.
            [g1, g2] = maxQTest.twoTriangleConjugates();
            h = maxQ(g1, g2);
            testCase.verifyEqual(size(unique(h.EcQ, 'rows'), 1), h.ne, ...
                'two rows of EcQ describe the same curve, so deduplication missed one');
        end

        function maxQRefusesAnythingThatIsNotAQuaCon(testCase)
            [g1, ~] = maxQTest.twoTriangleConjugates();
            testCase.verifyError(@() maxQ(g1, 42), 'PLQ:maxQ:input');
            testCase.verifyError(@() maxQ(QuaPol([1 0 1 0 0 0]), g1), 'PLQ:maxQ:input');
        end
    end

    methods (Static)
        function [g1, g2] = twoTriangleConjugates()
        % The conjugates of two DIFFERENT strictly convex quadratics on two different triangles.
        % Different Hessians is what makes the difference conic non-degenerate, which is the case
        % worth testing -- equal Hessians would make every split boundary a line.
            E = [1 2 1; 2 3 1; 3 1 1];  F = [1 0; 1 0; 1 0];
            f1 = QuaPol([0 0; 1 0; 0 1], E, [0 0 0 0 1 0 1 0 0 0], F);
            f2 = QuaPol([1 0; 1 1; 0 1], E, [0 0 0 0 4 1 3 -2 1 0], F);
            g1 = conjQ(f1);  g2 = conjQ(f2);
        end
    end
end
