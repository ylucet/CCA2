classdef ratQTest < matlab.unittest.TestCase
% Unit tests for ratQ, the exact rational coefficient layer the sym-free conjugate rests on.
%
% BUCKET: fast (closed-form integer arithmetic; no symbolic call on the compute path).
%
% These are deliberately heavier than the class looks like it needs, for the reason exactQTest
% gives: everything downstream inherits the guarantees. Two face functions merge exactly when
% their canonical forms are bitwise equal, and Step 3's cell count depends on that answer being
% exact -- DECISIONS.md 2026-08-17 measured what a nearly-right answer costs (a shared facet
% invisible to `merge`, and a cell count that then grows without bound).
%
% The oracle is the Symbolic Toolbox where an independent check is wanted, which is legitimate in a
% TEST and nowhere else: tests may use sym, the compute path may not.

    methods (Test)

        % ---- canonical form: the contract every downstream equality depends on ----------------

        function canonRemovesTheGcdAndForcesAPositiveDenominator(testCase)
            [n, d] = ratQ.canon([6 -9 12], 15);
            testCase.verifyEqual(n, [2 -3 4]);
            testCase.verifyEqual(d, 5);
            [n2, d2] = ratQ.canon([6 -9 12], -15);
            testCase.verifyEqual(n2, [-2 3 -4], 'a negative denominator moves to the numerator');
            testCase.verifyEqual(d2, 5);
        end

        function canonIsIdempotentAndFixesTheSpellingOfZero(testCase)
            [n, d] = ratQ.canon([0 0 0], 7);
            testCase.verifyEqual(n, [0 0 0]);
            testCase.verifyEqual(d, 1, 'zero over anything canonicalises to zero over one');
            [n2, d2] = ratQ.canon(n, d);
            testCase.verifyEqual([n2 d2], [n d], 'canon must be idempotent');
        end

        function equalValuesWithDifferentSpellingsCompareEqual(testCase)
        % THE POINT OF THE WHOLE CLASS. In the symbolic form these two can compare Unknown; in
        % doubles they can differ by an ULP. Here they are one object.
            testCase.verifyTrue(ratQ.eqRat([1 2 3], 2, [2 4 6], 4));
            testCase.verifyTrue(ratQ.eqRat([1 2 3], 2, [-3 -6 -9], -6));
            testCase.verifyFalse(ratQ.eqRat([1 2 3], 2, [1 2 3], 3));
        end

        % ---- conics have no scale, faces do ---------------------------------------------------

        function aConicIsPrimitiveAndSignNormalised(testCase)
            testCase.verifyEqual(ratQ.conic([2 -4 6 0 0 8]), [1 -2 3 0 0 4]);
            testCase.verifyEqual(ratQ.conic([-2 4 -6 0 0 -8]), [1 -2 3 0 0 4], ...
                'a conic and its negation are the same curve');
            testCase.verifyEqual(ratQ.conic([0 0 0 -3 6 -9]), [0 0 0 1 -2 3], ...
                'the sign is fixed by the FIRST NONZERO entry, not by entry 1');
        end

        function scalingAConicDoesNotChangeTheCurveButScalingAFaceChangesTheFunction(testCase)
        % The two canonical forms are different on purpose and this is the difference.
            testCase.verifyTrue(ratQ.sameConic([1 0 1 0 0 -1], [1000 0 1000 0 0 -1000]));
            testCase.verifyFalse(ratQ.eqRat([0 0 0 0 1 0 1 0 0 -1], 1, ...
                                            [0 0 0 0 1000 0 1000 0 0 -1000], 1));
        end

        function anAllZeroConicIsRefusedRatherThanAcceptedAsALine(testCase)
        % QuaPar spells a straight edge as an all-zero Ec row and recovers the line from the
        % edge's stored endpoint COORDINATES. H-form removes that dependence, so an edge that
        % names no curve is a defect here, not a shorthand.
            testCase.verifyError(@() ratQ.conic([0 0 0 0 0 0]), 'ratQ:zeroConic');
        end

        % ---- arithmetic ------------------------------------------------------------------------

        function addUsesTheLcmSoCoefficientsDoNotGrowForNoReason(testCase)
            [n, d] = ratQ.add([1 0], 6, [1 0], 10);
            testCase.verifyEqual(d, 15, 'lcm(6,10)=30, then the gcd of [5 0]+[3 0] with 30 is 2');
            testCase.verifyEqual(n, [4 0]);
            testCase.verifyEqual(n/d, [1/6 0] + [1/10 0], 'AbsTol', 0);
        end

        function subAndScaleAgreeWithExactRationalArithmetic(testCase)
            [n, d] = ratQ.sub([3 5], 4, [1 2], 8);
            testCase.verifyEqual(n/d, [3/4 - 1/8, 5/4 - 2/8], 'AbsTol', 0);
            [n2, d2] = ratQ.scale([3 5], 4, 2, 3);
            testCase.verifyEqual(n2/d2, [1/2, 5/6], 'AbsTol', 0);
        end

        function additionIsAssociativeOnValuesThatStressTheDenominators(testCase)
            a = {[7 -3], 12}; b = {[5 11], 18}; c = {[-1 4], 45};
            [n1, d1] = ratQ.add(a{:}, b{:});  [n1, d1] = ratQ.add(n1, d1, c{:});
            [n2, d2] = ratQ.add(b{:}, c{:});  [n2, d2] = ratQ.add(a{1}, a{2}, n2, d2);
            testCase.verifyTrue(ratQ.eqRat(n1, d1, n2, d2));
        end

        % ---- the difference conic, which the whole design turns on ------------------------------

        function theDifferenceOfTwoAdjacentQuadraticsIsDegenerate(testCase)
        % CONJ_FIELD_PROOF.md Theorem 3: two pieces of f that are ADJACENT give a difference conic
        % with b^2-4ac = 0. Here, two faces sharing a Hessian differ by an affine form, which is
        % the extreme case -- a LINE.
            g1 = ratQTest.face(2, 1, 3, 0, 0, 0);      % x^2 + xy + 3y^2/2
            g2 = ratQTest.face(2, 1, 3, 4, -5, 6);     % the same quadratic part
            c  = ratQ.diffConic(g1, 1, g2, 1);
            testCase.verifyEqual(ratQ.conicKind(c), 'line');
            testCase.verifyEqual(ratQ.discriminant(c), 0);
        end

        function anIndefiniteDifferenceIsHyperbolicAndADefiniteOneIsElliptic(testCase)
        % The sign of b^2-4ac is -det(H1-H2), so the classification is a statement about the two
        % Hessians and nothing else. Both branches are reachable from f*; the elliptic one is the
        % reason the Conic level of the lattice exists.
            g1 = ratQTest.face(2, 0, 2, 0, 0, 0);      % x^2 + y^2
            g2 = ratQTest.face(0, 0, 0, 0, 0, -1);     % the constant -1
            testCase.verifyEqual(ratQ.conicKind(ratQ.diffConic(g1, 1, g2, 1)), 'elliptic');
            g3 = ratQTest.face(2, 0, -2, 0, 0, 0);     % x^2 - y^2
            testCase.verifyEqual(ratQ.conicKind(ratQ.diffConic(g3, 1, g2, 1)), 'hyperbolic');
        end

        function theEllipseOfTheThreePieceCounterexampleIsReproducedExACTLY(testCase)
        % CONJ_FIELD_PROOF.md 7.5 exhibits a THREE-piece continuous PLQ whose f* has an edge on
        %     93 s1^2 + 38 s1 s2 - 6 s1 + 39 s2^2 - 482 s2 - 1003 = 0,   disc -13064
        % between two NON-adjacent pieces. That conic is what a QuaPar cannot store, so pinning it
        % here pins the reason the Conic level was added. Asserted as a property of the conic (its
        % canonical form and its exact discriminant), not as a golden output of any routine.
            c = ratQ.conic([93 38 39 -6 -482 -1003]);
            testCase.verifyEqual(c, [93 38 39 -6 -482 -1003], 'already primitive');
            testCase.verifyEqual(ratQ.discriminant(c), 38^2 - 4*93*39);
            testCase.verifyEqual(ratQ.discriminant(c), -13064, ...
                'the discriminant quoted in CONJ_FIELD_PROOF.md 7.5');
            testCase.verifyEqual(ratQ.conicKind(c), 'elliptic');
            testCase.verifyError(@() QuaPar.assertParabolicEdges(c), 'QuaPar:notParabola', ...
                'and this is exactly what QuaPar refuses to store');
        end

        function aCubicFaceIsRefusedByDiffConicRatherThanTruncated(testCase)
            cub = zeros(1,10); cub(1) = 1;             % x^3/6
            testCase.verifyError(@() ratQ.diffConic(cub, 1, zeros(1,10), 1), 'ratQ:cubicFace');
        end

        % ---- exactness is enforced, not hoped for ----------------------------------------------

        function chkRaisesOnANonIntegerAndOnOverflowRatherThanRounding(testCase)
            testCase.verifyError(@() ratQ.chk(1.5), 'ratQ:notInteger');
            testCase.verifyError(@() ratQ.chk(2^53), 'ratQ:overflow');
            testCase.verifyError(@() ratQ.chk(inf), 'ratQ:notFinite');
            testCase.verifyEqual(ratQ.chk(2^53 - 1), 2^53 - 1, 'the largest exact integer is fine');
        end

        function fromDoubleConvertsWhatIsExactAndRefusesWhatIsNot(testCase)
            [n, d] = ratQ.fromDouble([0.5 -0.25 3]);
            testCase.verifyEqual(d, 4);
            testCase.verifyEqual(n, [2 -1 12]);
            testCase.verifyEqual(n/d, [0.5 -0.25 3], 'AbsTol', 0);
            testCase.verifyError(@() ratQ.fromDouble(sqrt(2)), 'ratQ:notExact', ...
                'an irrational must be refused, never snapped -- DECISIONS.md on why snapping hung');
            testCase.verifyError(@() ratQ.fromDouble(1/3, 2), 'ratQ:notExact', ...
                'and so must a rational whose denominator exceeds the stated bound');
            % A NON-DYADIC rational must still be recovered: a double cannot hold 1/3 exactly, so
            % a bit-exact reading of the stored value gives a 2^54 denominator and is useless.
            % Reconstruction at the input boundary is the only thing that can work -- fromDouble's
            % header says why that is not the refuted snapping.
            [n3, d3] = ratQ.fromDouble([1/3 -3/10 50/11]);
            testCase.verifyEqual(d3, 330);
            testCase.verifyEqual(n3, [110 -99 1500]);
            testCase.verifyEqual(n3/d3, [1/3 -3/10 50/11], 'AbsTol', 0, ...
                'dividing back must reproduce the SAME doubles, bit for bit');
        end

        function aRoundTripThroughFromDoubleIsTheIdentityOnRationalData(testCase)
        % Property, not a golden value: for any vector of small-denominator rationals, converting
        % and dividing back must be bitwise exact.
            rng(20260824);
            for k = 1:200
                q = randi([-50 50], 1, 6) ./ randi([1 12], 1, 6);
                [n, d] = ratQ.fromDouble(q);
                testCase.verifyEqual(n/d, q, 'AbsTol', 0);
            end
        end

        % ---- an independent symbolic cross-check (legitimate in a test, not on the path) --------

        function theDifferenceConicAgreesWithTheSymbolicLevelSet(testCase)
        % Differential test against the Symbolic Toolbox: the canonical integer conic must define
        % the same curve as sym's own simplification of g1 - g2.
            rng(20260824);
            % NOT `syms x y real`. That declares an ASSUMPTION on the symbols named x and y, and
            % the assumption lives in the shared MuPAD session for the rest of the process -- not
            % in this function. `region.m` builds its own sym('x')/sym('y') and, with a stray
            % `real` assumption attached, four of regionTest's tests then failed with "Unable to
            % convert expression containing symbolic variables into double array". Measured: they
            % pass alone, fail when this suite runs first in the same MATLAB, and the fast bucket
            % runs the whole bucket in ONE process on purpose. Uniquely named symbols with no
            % assumptions cannot collide.
            x = sym('ratQTest_x');  y = sym('ratQTest_y');
            for k = 1:25
                a = randi([-6 6], 1, 6); b = randi([-6 6], 1, 6);
                g1 = ratQTest.face(a(1), a(2), a(3), a(4), a(5), a(6));
                g2 = ratQTest.face(b(1), b(2), b(3), b(4), b(5), b(6));
                if isequal(a, b), continue, end
                c = ratQ.diffConic(g1, 1, g2, 1);
                symDiff = 2 * (ratQTest.symOf(g1, x, y) - ratQTest.symOf(g2, x, y));
                symConic = c(1)*x^2 + c(2)*x*y + c(3)*y^2 + c(4)*x + c(5)*y + c(6);
                % Same CURVE means proportional, and the ratio must be a nonzero constant.
                r = simplify(symDiff / symConic);
                testCase.verifyTrue(isempty(symvar(r)) && double(r) ~= 0, ...
                    sprintf('case %d: the conic is not a nonzero constant multiple of g1-g2', k));
            end
        end
    end

    methods (Static)
        function c = face(c5, c6, c7, c8, c9, c10)
        % A quadratic face in CCA2's 10-wide weighted cubic basis (RatCon.m's `f`):
        %     value = c5*x^2/2 + c6*xy + c7*y^2/2 + c8*x + c9*y + c10
            c = [0 0 0 0, c5, c6, c7, c8, c9, c10];
        end

        function s = symOf(c, x, y)
            s = c(5)*x^2/2 + c(6)*x*y + c(7)*y^2/2 + c(8)*x + c(9)*y + c(10);
        end
    end
end
