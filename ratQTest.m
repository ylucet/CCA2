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

        % ---- the predicate leaf, and the arithmetic the per-piece closed forms need ------------

        function signQAgreesWithTheSignOfTheValueItself(testCase)
        % A definition test, not a golden one: signQ must agree with the sign of n/d wherever the
        % double division happens to be exact, and must be exact where the double is not.
            testCase.verifyEqual(ratQ.signQ([3 0 -3], 5), [1 0 -1]);
            testCase.verifyEqual(ratQ.signQ([3 0 -3], -5), [-1 0 1], ...
                'a negative denominator flips every sign');
            rng(20260903);
            for k = 1:50
                n = randi([-40 40], 1, 4);  d = randi([-40 40]);
                if d == 0, continue, end
                testCase.verifyEqual(ratQ.signQ(n, d), sign(n/d), ...
                    sprintf('case %d: signQ disagrees with the value it is the sign of', k));
            end
        end

        function divIsTheInverseOfScale(testCase)
        % The identity is the assertion: dividing by a value and multiplying it back is the
        % identity ON THE CANONICAL FORM, which is stronger than agreeing numerically.
            rng(20260903);
            for k = 1:40
                n1 = randi([-30 30], 1, 5);  d1 = randi([1 30]);
                n2 = randi([-30 30]);        d2 = randi([1 30]);
                if n2 == 0, continue, end
                [q, qd] = ratQ.div(n1, d1, n2, d2);
                [b, bd] = ratQ.scale(q, qd, n2, d2);
                [a, ad] = ratQ.canon(n1, d1);
                testCase.verifyEqual([b bd], [a ad], ...
                    sprintf('case %d: (v / w) * w is not v', k));
            end
        end

        function divRefusesZeroRatherThanReturningInf(testCase)
            testCase.verifyError(@() ratQ.div([1 2], 3, 0, 7), 'ratQ:divideByZero');
        end

        function combineDenSharesOneDenominatorWithoutChangingAnyValue(testCase)
        % What it must preserve is the VALUE of every row; what it must produce is ONE denominator.
            Ns = [1 2 3; 1 1 1; 5 0 -5];  ds = [2; 3; 10];
            [N, d] = ratQ.combineDen(Ns, ds);
            testCase.verifyEqual(d, 30, 'the shared denominator is the lcm, not the product');
            for i = 1:size(Ns,1)
                testCase.verifyTrue(ratQ.eqRat(N(i,:), d, Ns(i,:), ds(i)), ...
                    sprintf('row %d changed value', i));
            end
        end

        function combineDenIsIdempotentAndHandlesASingleRow(testCase)
            [N, d]   = ratQ.combineDen([4 6], 8);
            [N2, d2] = ratQ.combineDen(N, d);
            testCase.verifyEqual([N2 d2], [N d], 'combineDen must be idempotent');
            testCase.verifyTrue(ratQ.eqRat(N, d, [4 6], 8));
        end

        % ---- the exact linear solve: polyhedral vertices and the KKT systems -------------------

        function solve2SatisfiesItsOwnSystemExactly(testCase)
        % Recomputed, not pinned: A*x must equal b as an exact rational identity.
            rng(20260903);
            for k = 1:60
                n = 2 + (mod(k,2) == 0);                 % alternate 2x2 and 3x3
                A = randi([-9 9], n, n);  b = randi([-9 9], n, 1);
                % NOT `det(A) == 0`: det factorises in floating point, so a singular integer
                % matrix can come back as 4.4e-16 and a nonsingular one as 6.0000000000000009.
                % The true determinant IS an integer, so |det| < 0.5 decides singularity exactly.
                if abs(det(A)) < 0.5, continue, end
                [x, xd] = ratQ.solve2(A, b);
                testCase.verifyEqual(size(x), [n 1]);
                testCase.verifyTrue(all(A*x == b*xd), ...
                    sprintf('case %d: A*(x/xd) ~= b', k));
                testCase.verifyGreaterThan(xd, 0, 'the denominator is normalised positive');
            end
        end

        function solve2RefusesASingularSystemRatherThanReturningNaN(testCase)
            testCase.verifyError(@() ratQ.solve2([1 2; 2 4], [1; 2]), 'ratQ:singular');
        end

        % ---- exact Hessian and gradient, against the Symbolic Toolbox as an oracle -------------

        function hessQAndGradQAgreeWithSymbolicDifferentiation(testCase)
        % Differential test. The weighted basis exists so that the Hessian reads straight off
        % c5..c7 (RatCon.m says so); this is the check that it actually does.
            x = sym('ratQTest_hx');  y = sym('ratQTest_hy');
            rng(20260903);
            for k = 1:20
                c = randi([-8 8], 1, 6);  dn = randi([1 9]);
                g = ratQTest.face(c(1), c(2), c(3), c(4), c(5), c(6));
                s = ratQTest.symOf(g, x, y) / dn;

                [H, hd] = ratQ.hessQ(g, dn);
                testCase.verifyEqual(double(hessian(s, [x y])), H/hd, 'AbsTol', 0, ...
                    sprintf('case %d: exact Hessian disagrees with sym', k));

                [G, gd] = ratQ.gradQ(g, dn);
                pt = randi([-5 5], 2, 1);
                want = double(subs(gradient(s, [x y]), [x y], pt.'));
                testCase.verifyEqual(G*[pt; 1]/gd, want, 'AbsTol', 0, ...
                    sprintf('case %d: exact gradient disagrees with sym', k));
            end
        end

        function hessQRefusesACubicBecauseItsHessianIsNotConstant(testCase)
            testCase.verifyError(@() ratQ.hessQ([1 0 0 0 0 0 0 0 0 0], 1), 'ratQ:cubicFace');
        end

        % ---- the exact affine change of variables --------------------------------------------

        function substAffineAgreesWithSymbolicSubstitutionOnCubicsAndQuadratics(testCase)
        % Differential test against the Symbolic Toolbox, on the FULL cubic basis: h(u) = g(Mu+t)
        % must be the same polynomial sym gets by substituting. Cubic terms are included on
        % purpose -- RatCon stores f in the cubic basis and isConvex accepts one, so a routine that
        % silently truncated the cubic part would pass every quadratic test and be wrong.
            u = sym('ratQTest_su');  v = sym('ratQTest_sv');
            rng(20260903);
            for k = 1:20
                c  = randi([-5 5], 1, 10);
                d  = randi([1 6]);
                Mn = randi([-4 4], 2, 2);
                tn = randi([-4 4], 2, 1);
                md = randi([1 5]);
                % A random integer matrix is singular often enough to matter, and substAffine
                % REFUSES a singular change of variables by contract. Decide that exactly, with
                % the class's own determinant -- `det` here would be the floating-point test the
                % whole layer exists to avoid.
                if ratQ.detExact(Mn) == 0, continue, end

                [n2, d2] = ratQ.substAffine(c, d, Mn, tn, md);

                X = (Mn(1,1)*u + Mn(1,2)*v + tn(1)) / md;
                Y = (Mn(2,1)*u + Mn(2,2)*v + tn(2)) / md;
                want = ratQTest.symCubic(c, X, Y) / d;
                got  = ratQTest.symCubic(n2, u, v) / d2;
                testCase.verifyTrue(isAlways(expand(want - got) == 0, 'Unknown', 'false'), ...
                    sprintf('case %d: substAffine disagrees with symbolic substitution', k));
            end
        end

        function substAffineIsTheIdentityOnTheIdentityMap(testCase)
        % The map u -> u must return the value unchanged, in CANONICAL form -- a stronger
        % statement than returning something numerically equal.
            c = [1 -2 3 -4 5 -6 7 -8 9 -10];
            [n2, d2] = ratQ.substAffine(c, 3, eye(2), [0; 0], 1);
            [a, ad]  = ratQ.canon(c, 3);
            testCase.verifyEqual([n2 d2], [a ad]);
        end

        function substAffineComposesTheSameWayTheMapsDo(testCase)
        % g(M1(M2 u + t2) + t1) reached in one step or two must be the same object. This is the
        % property the x*y frame change and its read-back depend on, and it is what a sign or
        % transpose slip breaks -- neither of which a single round trip through the identity sees.
            rng(20260903);
            for k = 1:15
                c   = randi([-5 5], 1, 10);   d = randi([1 4]);
                M1  = randi([-3 3], 2, 2);   t1 = randi([-3 3], 2, 1);
                M2  = randi([-3 3], 2, 2);   t2 = randi([-3 3], 2, 1);
                % Both maps must be invertible, and so must their composite -- see the note in
                % substAffineAgreesWith... above for why this is decided exactly.
                if ratQ.detExact(M1) == 0 || ratQ.detExact(M2) == 0, continue, end

                [nA, dA] = ratQ.substAffine(c,  d,  M1, t1, 1);
                [nA, dA] = ratQ.substAffine(nA, dA, M2, t2, 1);

                [nB, dB] = ratQ.substAffine(c, d, M1*M2, M1*t2 + t1, 1);

                testCase.verifyEqual([nA dA], [nB dB], ...
                    sprintf('case %d: the two-step and one-step substitutions differ', k));
            end
        end

        function substAffineRefusesASingularChangeOfVariables(testCase)
        % A change of VARIABLES has to be invertible. A singular M is not a reparametrisation --
        % it collapses the plane onto a line, and the resulting face is a different function on a
        % lower-dimensional domain. Refusing names the caller's error where it happens.
            testCase.verifyError(@() ratQ.substAffine([0 0 0 0 1 0 1 0 0 0], 1, ...
                                                     [1 2; 2 4], [0; 0], 1), 'ratQ:singular');
        end

        % ---- exact 2-D feasibility, which decides whether a cell exists at all -----------------

        function feasible2DecidesTheStandardShapesCorrectly(testCase)
        % `strict` is the distinction between "nonempty" and "carries a two-dimensional face", and
        % the degenerate shapes are exactly where a cell subdivision goes wrong quietly.
            sq  = [1 0 0; -1 0 1; 0 1 0; 0 -1 1];    % the unit square
            emp = [1 0 -1; -1 0 0];                   % s1 >= 1 and s1 <= 0
            pt  = [1 0 0; -1 0 0; 0 1 0; 0 -1 0];     % the single point (0,0)
            seg = [1 0 0; -1 0 0; 0 1 0; 0 -1 1];     % a segment
            cone = [1 0 0; 0 1 0];                    % the first quadrant, unbounded

            testCase.verifyTrue (ratQ.feasible2(sq));   testCase.verifyTrue (ratQ.feasible2(sq, true));
            testCase.verifyFalse(ratQ.feasible2(emp));  testCase.verifyFalse(ratQ.feasible2(emp, true));
            testCase.verifyTrue (ratQ.feasible2(pt));   testCase.verifyFalse(ratQ.feasible2(pt, true), ...
                'a single point is nonempty but carries no 2-D face');
            testCase.verifyTrue (ratQ.feasible2(seg));  testCase.verifyFalse(ratQ.feasible2(seg, true), ...
                'nor does a segment');
            testCase.verifyTrue (ratQ.feasible2(cone)); testCase.verifyTrue (ratQ.feasible2(cone, true), ...
                'an unbounded cone has interior and no vertex enumeration would find it');
            testCase.verifyTrue (ratQ.feasible2(zeros(0,3)), 'no constraints is the whole plane');
        end

        function feasible2DecidesWhenItsOwnProductsPass2to53(testCase)
        % Rows well inside 2^53 whose DECISION needs a product outside it: s1 >= (B-1)/B against
        % s1 <= (C-1)/C with B, C ~ 1e9 compares (B-1)*C with (C-1)*B ~ 1e18. The answer is exact
        % arithmetic -- (C-1)/C > (B-1)/B iff C > B -- and must not be an overflow. Reached for
        % real by biconjQ of the poly2 kink at c = -2 and -0.5 (envelopeCellRegions, 2026-09-25),
        % whose own inputs are all below 32.
            B = 1e9 + 7;  C = 1e9 + 9;
            testCase.verifyTrue (ratQ.feasible2([B 0 -(B-1); -C 0 (C-1)]),        'C > B: an interval');
            testCase.verifyTrue (ratQ.feasible2([B 0 -(B-1); -C 0 (C-1)], true),  'with interior in s1 (s2 is free)');
            testCase.verifyFalse(ratQ.feasible2([C 0 -(C-1); -B 0 (B-1)]),        'swapped: empty');
            testCase.verifyTrue (ratQ.feasible2([B 0 -(B-1); -B 0 (B-1)]),        'equal: one line');
            testCase.verifyFalse(ratQ.feasible2([B 0 -(B-1); -B 0 (B-1)], true),  'equal: no interior');
            % the same comparison reached through the ELIMINATION of s2
            testCase.verifyFalse(ratQ.feasible2([C 1 -(C-1); -B 1 (B-1); 0 -1 0]), 'through FM: empty');
            testCase.verifyTrue (ratQ.feasible2([B 1 -(B-1); -C 1 (C-1); 0 -1 0]), 'through FM: nonempty');
        end

        function aCANCELLINGSumOfProductsIsExactOrRefusedNeverWrong(testCase)
        % The defect feasible2 had, asked of every ratQ routine that sums products: inputs whose
        % ANSWER is small but whose PRODUCTS pass 2^53, so each product is rounded and the rounded
        % difference is small enough to pass a check on the result. The contract is exact-or-refuse:
        % the true value, or ratQ:overflow -- never a plausible wrong number.
            B = 1e9;
            cases = { ...
                {'det 2x2', @() ratQ.detExact([B B+1; B-1 B]), 1}, ...
                {'det 3x3', @() ratQ.detExact([B B+1 0; B-1 B 0; 0 0 1]), 1}, ...
                {'det 3x3 off-diagonal', @() ratQ.detExact([0 0 1; B B+1 0; B-1 B 0]), 1}};
            for k = 1:numel(cases)
                c = cases{k};
                verifyExactOrRefused(testCase, c{2}, c{3}, c{1});
            end
            % add: 100000007/3 + n2/100000037, whose numerator over the lcm is exactly 1, while
            % 100000007 * 100000037 ~ 1e16 is past 2^53
            n2 = -3333334800000086;
            verifyExactOrRefused(testCase, @() addAsRow(100000007, 3, n2, 100000037), ...
                                 [1, 300000111], 'add');
        end

        function feasible2AgreesWithLinprogOnRandomSystems(testCase)
        % Differential test against an entirely different algorithm -- a simplex/interior-point LP
        % in floating point, from the Optimization Toolbox. Legitimate in a TEST and nowhere near
        % the compute path. 500 systems here; 1939 were run in the scratchpad with 0 disagreements.
            rng(20260903);
            opts = optimoptions('linprog', 'Display', 'none');
            for k = 1:500
                m = randi([2 6]);
                P = randi([-5 5], m, 3);
                if any(all(P(:,1:2) == 0, 2)), continue, end
                [~, ~, flag] = linprog(zeros(2,1), -P(:,1:2), P(:,3), [], [], [], [], opts);
                testCase.verifyEqual(ratQ.feasible2(P), flag == 1, ...
                    sprintf('case %d: %s', k, mat2str(P)));
            end
        end

        function strictFeasibilityImpliesFeasibility(testCase)
        % A property no single example pins: an interior point is a point.
            rng(20260903);
            for k = 1:300
                P = randi([-5 5], randi([2 6]), 3);
                if ratQ.feasible2(P, true)
                    testCase.verifyTrue(ratQ.feasible2(P), sprintf('case %d', k));
                end
            end
        end

        % ---- constant-sign conics, which is the only exact handle on CURVED emptiness ----------

        function conicSignRecognisesFormsThatNeverChangeSign(testCase)
            testCase.verifyEqual(ratQ.conicSign([1 0 1 0 0 1]),  1, ...
                'x^2 + y^2 + 1 is positive everywhere');
            testCase.verifyEqual(ratQ.conicSign([-1 0 -1 0 0 -1]), -1, ...
                'and its negation is negative everywhere');
            testCase.verifyEqual(ratQ.conicSign([1 0 1 0 0 0]),  1, ...
                'x^2 + y^2 touches zero at the origin and is still nonnegative');
            testCase.verifyEqual(ratQ.conicSign([1 0 1 0 0 -1]), 0, ...
                'the unit circle genuinely separates the plane');
            testCase.verifyEqual(ratQ.conicSign([0 0 0 1 0 0]),  0, ...
                'a line takes both signs');
            testCase.verifyEqual(ratQ.conicSign([1 0 -1 0 0 0]), 0, ...
                'a hyperbolic form takes both signs');
            testCase.verifyEqual(ratQ.conicSign([1 0 0 0 0 0]),  1, ...
                'x^2 is a degenerate PSD form, nonnegative everywhere');
        end

        function conicSignAgreesWithSamplingOnRandomConics(testCase)
        % Differential test: a form claimed one-signed must never be caught taking the other sign,
        % and a form claimed two-signed must actually be seen taking both. Sampling can only refute
        % the first claim and confirm the second, which is exactly how the two halves are checked.
            rng(20260904);
            [X, Y] = meshgrid(linspace(-8, 8, 61));
            for k = 1:300
                c = randi([-4 4], 1, 6);
                if all(c == 0), continue, end
                v = c(1)*X.^2 + c(2)*X.*Y + c(3)*Y.^2 + c(4)*X + c(5)*Y + c(6);
                s = ratQ.conicSign(c);
                switch s
                    case 1
                        testCase.verifyGreaterThanOrEqual(min(v(:)), -1e-9, ...
                            sprintf('case %d claimed nonnegative: %s', k, mat2str(c)));
                    case -1
                        testCase.verifyLessThanOrEqual(max(v(:)), 1e-9, ...
                            sprintf('case %d claimed nonpositive: %s', k, mat2str(c)));
                end
            end
        end

        function isPSD3NeedsEveryPrincipalMinorNotJustTheLeadingOnes(testCase)
        % Sylvester's LEADING-minor criterion characterises positive DEFINITE and is false in both
        % directions for semidefinite -- this matrix has all three leading minors zero and is not
        % PSD, so a leading-minor implementation would return true and conicSign would then declare
        % a two-signed conic one-signed, silently deleting a real constraint.
            M = diag([0 -1 0]);
            testCase.verifyEqual(det(M(1,1)), 0);
            testCase.verifyFalse(ratQ.isPSD3(M));
            testCase.verifyTrue(ratQ.isPSD3(diag([0 0 0])));
            testCase.verifyTrue(ratQ.isPSD3(diag([2 3 5])));
            testCase.verifyFalse(ratQ.isPSD3([1 0 0; 0 1 2; 0 2 1]));
        end

        % ---- the exact degree-<=4 real-algebraic kernel -----------------------------------------

        function sturmCountsRootsExactly(testCase)
        % Sturm's theorem: the number of DISTINCT real roots of p in (a,b] is V(a) - V(b), where V
        % counts sign changes in the Sturm chain. Checked against polynomials whose roots are known
        % by construction, so the test does not depend on any root finder.
            % (x-1)(x-2)(x-3) = x^3 - 6x^2 + 11x - 6
            p = [1 -6 11 -6];
            testCase.verifyEqual(ratQ.countRootsIn(p, -10, 10), 3);
            testCase.verifyEqual(ratQ.countRootsIn(p, 0, 2), 2, '1 and 2 lie in (0,2]');
            testCase.verifyEqual(ratQ.countRootsIn(p, 2, 10), 1, 'only 3 lies in (2,10]');
            testCase.verifyEqual(ratQ.countRootsIn(p, 4, 10), 0);

            % x^2 - 2: two irrational roots
            testCase.verifyEqual(ratQ.countRootsIn([1 0 -2], -2, 2), 2);
            testCase.verifyEqual(ratQ.countRootsIn([1 0 -2], 0, 2), 1);

            % x^2 + 1: none
            testCase.verifyEqual(ratQ.countRootsIn([1 0 1], -10, 10), 0);

            % a REPEATED root counts ONCE -- Sturm counts distinct roots, and the chain is built
            % from the squarefree part, so (x-1)^2 has exactly one.
            testCase.verifyEqual(ratQ.countRootsIn([1 -2 1], -5, 5), 1);
        end

        function isolateRootsSeparatesEveryRealRoot(testCase)
        % Each returned interval must contain exactly one root, and together they must account for
        % all of them. Both halves are checked with countRootsIn, which is independent of how the
        % intervals were produced.
            polys = {[1 -6 11 -6], [1 0 -2], [1 -2 1], [1 0 0 0 -1], [3 -24 10 160 -96]};
            for i = 1:numel(polys)
                p = polys{i};
                I = ratQ.isolateRoots(p);
                total = 0;
                for k = 1:size(I,1)
                    n = ratQ.countRootsIn(p, I(k,1), I(k,2));
                    testCase.verifyEqual(n, 1, ...
                        sprintf('poly %d interval %d holds %d roots, not 1', i, k, n));
                    total = total + n;
                end
                lo = min([I(:,1); -1]) - 1;  hi = max([I(:,2); 1]) + 1;
                testCase.verifyEqual(total, ratQ.countRootsIn(p, lo, hi), ...
                    sprintf('poly %d: the intervals miss a root', i));
            end
        end

        function theS4QuarticFromTheFieldProofIsIsolatedCorrectly(testCase)
        % CONJ_FIELD_PROOF.md section 6's quartic -- 3t^4 - 24t^3 + 10t^2 + 160t - 96 -- is
        % irreducible with Galois group S4, and its relevant root is recorded there to 18 digits.
        % This is the case the whole kernel exists for: a vertex that no tower of square roots can
        % reach.
            p = [3 -24 10 160 -96];
            I = ratQ.isolateRoots(p);
            want = 0.608050881512364091;
            hit = find(I(:,1) <= want & want <= I(:,2));
            testCase.verifyNumElements(hit, 1, 'the recorded root lies in exactly one interval');
        end

        function signAtDecidesTheSignOfAPolynomialAtAnAlgebraicNumber(testCase)
        % THE PREDICATE THE MESH NEEDS. Given a root alpha of p isolated in [a,b], decide sign(q(alpha))
        % EXACTLY -- no tolerance, and correct even when q(alpha) is zero.
            % alpha = sqrt(2), the positive root of x^2 - 2
            p = [1 0 -2];
            I = ratQ.isolateRoots(p);
            pos = I(I(:,1) >= 0, :);
            testCase.verifyNumElements(pos(:,1), 1);

            testCase.verifyEqual(ratQ.signAt([1 0 -2], p, pos), 0, ...
                'q = p vanishes at its own root');
            testCase.verifyEqual(ratQ.signAt([1 0], p, pos), 1, 'sqrt(2) > 0');
            testCase.verifyEqual(ratQ.signAt([1 0 -1], p, pos), 1, 'sqrt(2)^2 - 1 = 1 > 0');
            testCase.verifyEqual(ratQ.signAt([1 0 -3], p, pos), -1, 'sqrt(2)^2 - 3 = -1 < 0');
            testCase.verifyEqual(ratQ.signAt([1 -2], p, pos), -1, 'sqrt(2) - 2 < 0');
            testCase.verifyEqual(ratQ.signAt([1 -1], p, pos), 1, 'sqrt(2) - 1 > 0');
        end

        function signAtIsExactWhereFloatingPointIsAmbiguous(testCase)
        % The case a tolerance cannot decide: q vanishes at alpha to any precision one cares to
        % name, because it vanishes EXACTLY. (x^2-2) and (x^4-4) share the root sqrt(2).
            p = [1 0 -2];
            I = ratQ.isolateRoots(p);
            pos = I(I(:,1) >= 0, :);
            testCase.verifyEqual(ratQ.signAt([1 0 0 0 -4], p, pos), 0, ...
                'x^4 - 4 vanishes at sqrt(2), and the kernel must say 0 rather than a small number');
            % and a polynomial that is merely SMALL there is not zero
            testCase.verifyEqual(ratQ.signAt([1000000 0 -2000001], p, pos), -1, ...
                '10^6 x^2 - (2*10^6 + 1) is -1 at sqrt(2): small, but negative');
        end

        function signAtAgreesWithDoubleEvaluationWhereverThatIsSafe(testCase)
        % Differential test against MATLAB's `roots`, an entirely different algorithm: away from a
        % zero, the exact sign and the sign of q at the numerically computed root must agree.
        %
        % THE COMPARISON POINT IS THE ROOT, NOT THE INTERVAL'S MIDPOINT. Written the other way this
        % test failed 27 times with the code entirely correct: an isolating interval can be wide --
        % [0,3] for a root near 0.32 -- and q's sign at the midpoint has nothing to do with its sign
        % at the root. The signAt answers were right every time; the oracle was evaluating in the
        % wrong place.
            rng(20260905);
            checked = 0;
            for t = 1:60
                p = randi([-6 6], 1, 5);
                if p(1) == 0, continue, end
                I = ratQ.isolateRoots(p);
                if isempty(I), continue, end
                rr = roots(p);
                rr = real(rr(abs(imag(rr)) < 1e-9));
                for k = 1:size(I,1)
                    inIv = rr(rr >= I(k,1) - 1e-9 & rr <= I(k,2) + 1e-9);
                    if numel(inIv) ~= 1, continue, end      % ambiguous numerically: skip
                    alpha = inIv(1);
                    for j = 1:3
                        q = randi([-5 5], 1, 3);
                        v = polyval(q, alpha);
                        if abs(v) < 1e-6, continue, end     % too near a zero for the oracle
                        testCase.verifyEqual(ratQ.signAt(q, p, I(k,:)), sign(v), ...
                            sprintf('case %d root %d probe %d', t, k, j));
                        checked = checked + 1;
                    end
                end
            end
            testCase.verifyGreaterThan(checked, 50, 'the sweep must actually check something');
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
        % ---- the SYLVESTER RESULTANT -----------------------------------------------------------
        % The oracle is what a resultant MEANS -- R(x0) = 0 exactly when the two polynomials share
        % a root in y at x = x0 -- checked by exhibiting the shared root, not by comparing against
        % another implementation of the same formula.

        function theRESULTANTVanishesExactlyWhereTheyShareARootInY(testCase)
        % P = y^2 - x  and  Q = y - x. They share a root in y when x = x^2, i.e. at x = 0 and
        % x = 1, and nowhere else. So the resultant must be a nonzero multiple of x^2 - x.
        %   P(i+1,j+1) = coeff of x^i y^j:  y^2 -> (1,3);  -x -> (2,1)
            P = [0 0 1; -1 0 0];
            Q = [0 1; -1 0];                       % y - x
            R = ratQ.resultantY(P, Q);
            testCase.verifyEqual(numel(R), 3, 'a quadratic in x');
            testCase.verifyEqual(R/R(end), [0 -1 1], 'AbsTol', 0, 'proportional to x^2 - x');

            % and the meaning, checked directly: at x = 2 the two have no common y
            rootsP = roots([1 0 -2]);              % y^2 - 2
            testCase.verifyGreaterThan(min(abs(rootsP - 2)), 0.1, ...
                'y = x = 2 is not a root of y^2 - 2, so the resultant must not vanish at x = 2');
            val = R(1) + R(2)*2 + R(3)*4;
            testCase.verifyNotEqual(val, 0);
        end

        function theRESULTANTOfTwoCIRCLESIsTheirRadicalAxisCondition(testCase)
        % Two unit circles, centres (0,0) and (2,0): x^2 + y^2 - 1 and (x-2)^2 + y^2 - 1. They meet
        % only where the radical axis x = 1 crosses them, and there y^2 = 0 -- a double root. So the
        % resultant in y must vanish only at x = 1, to order 2, i.e. be proportional to (x-1)^2.
            P = [-1 0 1; 0 0 0; 1 0 0];            % x^2 + y^2 - 1
            Q = [3 0 1; -4 0 0; 1 0 0];            % (x-2)^2 + y^2 - 1 = x^2 - 4x + 3 + y^2
            R = ratQ.trim1(ratQ.resultantY(P, Q));
            testCase.verifyEqual(numel(R), 3, 'quadratic in x');
            testCase.verifyEqual(R/R(end), [1 -2 1], 'AbsTol', 0, 'proportional to (x-1)^2');
        end

        function theRESULTANTIsZeroWhenTheySHAREAFactor(testCase)
        % A common factor means a common root for EVERY x, so the resultant vanishes identically.
        % This is the degenerate case a caller has to handle, so it must be reported as the zero
        % polynomial rather than as an error.
            P = ratQ.polyMulBiv([0 1; 1 0], [0 1; 2 0]);   % (y+x)(y+2x)
            Q = ratQ.polyMulBiv([0 1; 1 0], [0 1; 5 0]);   % (y+x)(y+5x)
            R = ratQ.resultantY(P, Q);
            testCase.verifyTrue(all(R == 0), 'a shared factor gives the zero resultant');
        end

        % ---- the BIG-INTEGER fallback ----------------------------------------------------------
        % Added 2026-09-06. It is entered only where the double path would RAISE, so the tests have
        % to reach it deliberately: the whole point is that ordinary work never sees it.

        function bigIntegersROUNDTRIPAndRefuseWhatWillNotFit(testCase)
            A = [3 0 5; 0 7 0];
            testCase.verifyEqual(ratQ.fromBig(ratQ.toBig(A)), A, 'exact round trip');

            % a product past 2^53 is computed exactly and REFUSES to come back as a double
            D = ratQ.toBig([9007199254740991 0]);
            E = ratQ.bigMul(D, D);
            testCase.verifyEqual(char(E.c{1}.toString()), '81129638414606663681390495662081', ...
                'the big product is exact');
            testCase.verifyError(@() ratQ.fromBig(E), 'ratQ:overflow', ...
                'and coming back is refused rather than rounded');
        end

        function bigARITHMETICAgreesWithTheDoublePathWhereBothApply(testCase)
        % The differential check that matters for a fallback: on inputs the fast path handles, the
        % two must give the same answer, or the fallback is computing something else.
            rng(20260906);
            for t = 1:5
                A = randi([-6 6], 3, 2);   B = randi([-6 6], 2, 3);
                fast = ratQ.polyMulN(A, B);
                big  = ratQ.fromBig(ratQ.bigMul(ratQ.toBig(A), ratQ.toBig(B)));
                testCase.verifyEqual(big, fast, sprintf('product %d', t));
                testCase.verifyEqual(ratQ.fromBig(ratQ.bigAdd(ratQ.toBig(A), ratQ.toBig(A))), 2*A);
            end
        end

        function theRESULTANTFallsBackAndStillAgrees(testCase)
        % resultantVarBig returns the big form; on an input the double path handles, converting it
        % must reproduce resultantVar exactly.
            P = [0 0 1; -1 0 0];   Q = [0 1; -1 0];
            fast = ratQ.trimN(ratQ.resultantVar(P, Q, 2, 2));
            big  = ratQ.trimN(ratQ.fromBig(ratQ.resultantVarBig(P, Q, 2, 2)));
            testCase.verifyEqual(big, fast);
        end

        function bigDivExactDIVIDESOrRAISES(testCase)
            A = zeros(2,2);  A(2,1) = 1;  A(1,2) =  1;      % x + y
            B = zeros(2,2);  B(2,1) = 1;  B(1,2) = -2;      % x - 2y
            P = ratQ.toBig(ratQ.polyMulN(A, B));
            Q = ratQ.fromBig(ratQ.bigDivExact(P, ratQ.toBig(A), 2));
            testCase.verifyEqual(ratQ.trimN(Q), B);
            testCase.verifyError(@() ratQ.bigDivExact(ratQ.toBig(A), ratQ.toBig(B), 2), ...
                'ratQ:notDivisible');
        end

        function contentOfIsTheGCDAndComposesAcrossARRAYS(testCase)
        % contentOf exists because coefficients that share a scale must be reduced TOGETHER. Two
        % coefficients of one polynomial in an eliminated variable are exactly that case, and
        % reducing them independently moves the polynomial's roots -- measured, it made a
        % substituted root come out at half its true value.
            testCase.verifyEqual(ratQ.contentOf([20 -104 56 -40]), 4);
            testCase.verifyEqual(ratQ.contentOf([0 0 0]), 0, 'zero is gcd''s identity');
            testCase.verifyEqual(ratQ.contentOf([7 0 0]), 7);
            % the composing case: a common factor of two arrays that each have more of their own
            testCase.verifyEqual(gcd(ratQ.contentOf([6 12]), ratQ.contentOf([10 20])), 2);
        end

        function polyDivExactDIVIDESOrRAISES(testCase)
        % Exact multivariate division. The oracle is a product built on purpose: (x+y)(x-2y)
        % divided by (x+y) must give back (x-2y) exactly.
            A = zeros(2,2);  A(2,1) = 1;  A(1,2) =  1;      % x + y
            B = zeros(2,2);  B(2,1) = 1;  B(1,2) = -2;      % x - 2y
            P = ratQ.polyMulN(A, B);
            testCase.verifyEqual(ratQ.trimN(ratQ.polyDivExact(P, A, 2)), B);
            testCase.verifyEqual(ratQ.trimN(ratQ.polyDivExact(P, B, 2)), A);
        end

        function polyDivExactREFUSESRatherThanRounding(testCase)
        % A remainder or a non-integral quotient coefficient raises. Nothing here truncates, which
        % is the whole reason this layer exists.
            A = zeros(2,2);  A(2,1) = 1;  A(1,2) =  1;
            B = zeros(2,2);  B(2,1) = 1;  B(1,2) = -2;
            testCase.verifyError(@() ratQ.polyDivExact(A, B, 2), 'ratQ:notDivisible');
            testCase.verifyError(@() ratQ.polyDivExact(A, 0, 2), 'ratQ:divideByZero');
            % 2x + 2y over x + y is exact; 2x + 3y over x + y is not
            testCase.verifyEqual(ratQ.trimN(ratQ.polyDivExact(2*A, A, 2)), 2);
            C = zeros(2,2);  C(2,1) = 2;  C(1,2) = 3;
            testCase.verifyError(@() ratQ.polyDivExact(C, A, 2), 'ratQ:notDivisible');
        end

        function polyGcdBivFindsASHAREDFactor(testCase)
        % Two products built on a common factor: the gcd must be that factor, primitive.
            F  = zeros(3,2);  F(3,1) = 1;  F(1,2) =  1;     % x^2 + y
            G1 = zeros(2,2);  G1(2,1) = 1; G1(1,2) =  3;    % x + 3y
            G2 = zeros(2,2);  G2(2,1) = 1; G2(1,2) = -1;    % x - y
            G  = ratQ.polyGcdBiv(ratQ.polyMulN(F,G1), ratQ.polyMulN(F,G2));
            testCase.verifyEqual(ratQ.trimN(G), F, 'the shared factor, primitive');
        end

        function polyGcdBivRAISESRatherThanSwellingSilently(testCase)
        % The boundary of what a polynomial remainder sequence can do inside 2^53, pinned so it
        % cannot drift. These are the shapes the convex envelope's own bridge polynomial has --
        % large coefficients on a modest degree -- and a PRS multiplies by a leading coefficient at
        % every elimination. Raising is the CORRECT behaviour for an exact library; the wrong one
        % would be a rounded answer.
        %
        % The fix, if this is ever needed at scale, is a MODULAR gcd, not a better PRS. See
        % polyGcdBiv's header and .claude/gcdSwellProbe.m.
            A = zeros(3,2);  A(3,1) =  2664000000000;  A(1,2) = -1332000000000;
            B = zeros(3,2);  B(3,1) =  1332000000000;  B(1,2) =   666000000000;
            % after the integer content comes off these are tiny, so this pair SUCCEEDS
            testCase.verifyEqual(ratQ.trimN(ratQ.polyGcdBiv(A, B)), 1, ...
                'a large common factor is content, and content comes off first');

            % genuinely large coefficients on a genuine degree do not survive the swell
            rng(20260906);
            P = round(1e11 * (2*rand(7,5) - 1));
            Q = round(1e11 * (2*rand(6,4) - 1));
            f = @() ratQ.polyGcdBiv(P, Q);
            try
                f();
                ok = true;
            catch e
                ok = false;
                testCase.verifyEqual(e.identifier, 'ratQ:overflow', ...
                    'and when it cannot, it says so exactly');
            end
            if ok
                testCase.verifyTrue(true, 'this pair happened to stay inside the ceiling');
            end
        end

        function polyGcdBivReturnsACONSTANTWhenCoprime(testCase)
            G1 = zeros(2,2);  G1(2,1) = 1; G1(1,2) =  3;
            G2 = zeros(2,2);  G2(2,1) = 1; G2(1,2) = -1;
            G = ratQ.trimN(ratQ.polyGcdBiv(G1, G2));
            testCase.verifyTrue(isscalar(G), 'coprime inputs give a constant gcd');
            testCase.verifyEqual(abs(G), 1);
        end

        function contentFreeDividesOutTheCOMMONFactorAndNothingElse(testCase)
        % A polynomial EQUATION has no scale, and in a resultant the scale is not free: an NxN
        % Sylvester determinant multiplies the content roughly N times over. Measured on the bridge
        % elimination -- leaving a factor of 4 in one conic took it from 2.7e15 to an overflow at
        % 1.3e16, and reducing every intermediate brought it down to 2.7e12.
            testCase.verifyEqual(ratQ.contentFree([20 -104 56 -40]), [5 -26 14 -10]);
            testCase.verifyEqual(ratQ.contentFree([5 -26 14 -10]), [5 -26 14 -10], ...
                'already primitive, so unchanged');
            testCase.verifyEqual(ratQ.contentFree([0 0 0]), [0 0 0], ...
                'an all-zero array has no content to remove');
            testCase.verifyEqual(ratQ.contentFree([-4 -8]), [-4 -8]/4, ...
                'the SIGN is left alone: this is content, not orientation');
            A = zeros(2,2,2);  A(1,1,1) = 6;  A(2,2,2) = -9;
            R = ratQ.contentFree(A);
            testCase.verifyEqual(R(1,1,1), 2);
            testCase.verifyEqual(R(2,2,2), -3);
        end

        function theGENERALResultantAgreesWithTheBivariateOne(testCase)
        % resultantY is the nv = 2 case of resultantVar, so on a bivariate input they must agree.
        % Written as a test rather than a comment because the general one PERMUTES its arrays, and
        % a permutation read backwards would eliminate the wrong variable and still return
        % something plausible.
            P = [0 0 1; -1 0 0];              % y^2 - x
            Q = [0 1; -1 0];                  % y - x
            R1 = ratQ.trim1(ratQ.resultantY(P, Q));
            R2 = ratQ.trimN(ratQ.resultantVar(P, Q, 2, 2));
            testCase.verifyEqual(R2(:).', R1, 'the two routines must give the same polynomial');
        end

        function theGENERALResultantEliminatesTheNAMEDVariable(testCase)
        % Three variables (x,y,z), eliminating the MIDDLE one: y - x and y - z share a root in y
        % exactly when x = z, so the resultant must be x - z up to sign. Eliminating the wrong
        % variable would give x - y or y - z instead, all three of which look equally reasonable --
        % which is why the answer is checked and not the shape.
            A = zeros(2,2,1);  A(1,2,1) =  1;  A(2,1,1) = -1;     % y - x
            B = zeros(1,2,2);  B(1,2,1) =  1;  B(1,1,2) = -1;     % y - z
            R = ratQ.trimN(ratQ.resultantVar(A, B, 2, 3));
            testCase.verifyEqual(size(R,2), 1, 'y is gone');
            M = reshape(R, size(R,1), size(R,3));                 % rows x, cols z
            testCase.verifyEqual(M/M(2,1), [0 -1; 1 0], 'AbsTol', 0, 'proportional to x - z');
        end

        function theGENERALResultantRefusesAVariableThatIsNotThere(testCase)
            A = zeros(2,2);  A(1,2) = 1;  A(2,1) = -1;
            testCase.verifyError(@() ratQ.resultantVar(A, A, 3, 2), 'ratQ:badVariable');
            testCase.verifyError(@() ratQ.resultantVar(A, A, 0, 2), 'ratQ:badVariable');
        end

        function theRESULTANTRefusesWhatItCannotEliminate(testCase)
            testCase.verifyError(@() ratQ.resultantY([1 2 3].', [0 1]), 'ratQ:resultantDegree');
            big = zeros(1, 6);  big(end) = 1;             % y^5
            testCase.verifyError(@() ratQ.resultantY(big, big), 'ratQ:resultantTooLarge');
        end

        % ---- the ARGUMENT VALIDATORS, which had no test ------------------------------------
        % Added 2026-09-05. ratQ is the layer everything exact rests on, and its guards are what
        % stand between an exact answer and a plausible wrong one -- but an untested validator is
        % the worst state a check can be in: it passes forever whatever it does. The same
        % measurement that took QuaCon from 70.7% to 100% pointed here next.
        %
        % Only guards a CALLER can reach are tested. The rest are internal invariants and are
        % classified in SUPPORT_MATRIX.md rather than chased -- notably signAtNoConverge,
        % isolateNoConverge and pseudoRemNotIntegral, which the file's own comments argue cannot
        % fire for a squarefree integer polynomial.

        function aDenominatorMustBeANONZEROSCALAR(testCase)
            testCase.verifyError(@() ratQ.canon([1 2], [1 2]), 'ratQ:badDenominator');
            testCase.verifyError(@() ratQ.canon([1 2], 0),     'ratQ:badDenominator');
            testCase.verifyError(@() ratQ.signQ([1 2], [1 2]), 'ratQ:badDenominator');
            testCase.verifyError(@() ratQ.signQ([1 2], 0),     'ratQ:badDenominator');
        end

        function divideByAVECTORIsRefusedBecauseItIsNotDefined(testCase)
            testCase.verifyError(@() ratQ.div([1 2], 1, [1 2], 1), 'ratQ:badDivisor');
            testCase.verifyError(@() ratQ.div([1 2], 1, 0, 1),     'ratQ:divideByZero');
        end

        function combineDenNeedsONEDenominatorPerROW(testCase)
            testCase.verifyError(@() ratQ.combineDen([1 2; 3 4], [1 1 1]), 'ratQ:sizeMismatch');
            testCase.verifyError(@() ratQ.combineDen([1 2; 3 4], [1; 0]),  'ratQ:badDenominator');
        end

        function solve2TakesA2x2OrA3x3AndNothingElse(testCase)
            testCase.verifyError(@() ratQ.solve2([1 2; 3 4], [1; 2; 3]), 'ratQ:badSystem');
            testCase.verifyError(@() ratQ.solve2([1 2 3], [1; 2]),       'ratQ:badSystem');
            testCase.verifyError(@() ratQ.solve2(eye(4), ones(4,1)),     'ratQ:badSystem');
        end

        function detExactCoversOneByOneAndRefusesAboveThree(testCase)
        % The 1x1 case is a real branch, not a guard, and was the one uncovered line here that
        % returns a value rather than raising.
            testCase.verifyEqual(ratQ.detExact(7), 7);
            testCase.verifyError(@() ratQ.detExact(eye(4)), 'ratQ:badSystem');
        end

        function aCUBICHasNoConstantHessianAndNoAffineGradient(testCase)
        % The documented contract of hessQ/gradQ: a cubic's Hessian varies with the point, so
        % returning one would be a wrong ANSWER rather than a missing feature.
            cubic = [1 0 0 0, 0, 0, 0, 0, 0, 0];      % x^3 term present
            testCase.verifyError(@() ratQ.hessQ(cubic, 1), 'ratQ:cubicFace');
            testCase.verifyError(@() ratQ.gradQ(cubic, 1), 'ratQ:cubicFace');
        end

        function substAffineChecksTheSHAPEOfItsMap(testCase)
            q = [0 0 0 0, 2, 0, 2, 0, 0, 0];
            testCase.verifyError(@() ratQ.substAffine(q, 1, eye(3), [0;0]), 'ratQ:badMap');
            testCase.verifyError(@() ratQ.substAffine(q, 1, eye(2), [0;0;0]), 'ratQ:badMap');
            testCase.verifyError(@() ratQ.substAffine(q, 1, eye(2), [0;0], 0), ...
                'ratQ:badDenominator');
            testCase.verifyError(@() ratQ.substAffine(q, 1, zeros(2), [0;0]), 'ratQ:singular');
        end

        function aPRODUCTAboveTheCUBICBasisIsRefusedRatherThanTruncated(testCase)
        % Truncating would silently change the function, which is the one failure this layer is
        % built to prevent. x^2 * x^2 is degree 4 and has nowhere to go.
            A = zeros(4); A(3,1) = 1;                  % x^2 in the (i,j) degree grid
            testCase.verifyError(@() ratQ.polyMul(A, A), 'ratQ:degreeOverflow');
        end

        function fromDoubleREFUSESWhatItCannotRepresentExactly(testCase)
        % The counting argument in fromDouble's header, exercised: an irrational is rejected
        % because the best rational with denominator <= maxDen is off by ~1/q^2, four orders
        % coarser than the window it accepts.
            testCase.verifyError(@() ratQ.fromDouble(sqrt(2)), 'ratQ:notExact');
            [n, d] = ratQ.fromDouble(50/11);
            testCase.verifyEqual([n d], [50 11], 'and a genuine rational still converts');
        end

        function theBigPseudoRemainderAGREESWithTheDoubleOneWhereverBothCanRun(testCase)
        % THE ORACLE FOR THE BIG PATH, and the only reason to trust it. The thing being replaced is
        % the thing that checks the replacement: wherever the double path can run at all, the
        % BigInteger path must return the SAME vector, bitwise. Random small integer polynomials,
        % because the interesting behaviour is the elimination bookkeeping and small coefficients
        % are what let the double path run for comparison.
            rng(2026);
            n = 0;
            for it = 1:300
                da = randi([1 5]);  db = randi([1 da]);
                a = randi([-6 6], 1, da + 1);
                b = randi([-6 6], 1, db + 1);
                if a(1) == 0 || b(1) == 0, continue, end
                ap = ratQ.polyPrimitive(a);  bp = ratQ.polyPrimitive(b);
                try
                    rd = ratQ.polyPseudoRemDouble(ap, bp, numel(ap)-1, numel(bp)-1);
                catch
                    continue                        % the double path cannot run: nothing to compare
                end
                rb = ratQ.polyPseudoRemBig(a, b);
                testCase.verifyEqual(ratQ.polyTrim(rd), ratQ.polyTrim(rb), ...
                    sprintf('a = [%s], b = [%s]', num2str(a), num2str(b)));
                n = n + 1;
            end
            testCase.verifyGreaterThan(n, 100, ...
                'the generator must actually produce comparable cases, or this test asserts nothing');
        end

        function theSturmChainSURVIVESAnInputThatOverflowsTheDoublePath(testCase)
        % THE CASE THAT MOTIVATED THE BIG PATH. ratQ.sturmChain squares its leading coefficient at
        % every pseudo-remainder step, so a quartic with large coefficients overflowed 2^53 and the
        % Sturm count could not be taken at all -- which is what blocked the bridge region's exact
        % membership test on the anisotropic instance (DECISIONS.md 2026-09-06).
        %
        % The root count is asserted against a property, not a stored value: this quartic is
        % (x^2 - 2)*(x^2 - 3) scaled, whose roots are +-sqrt(2) and +-sqrt(3), so exactly two lie
        % in (0, 2] and exactly one in (0, 1.5].
            p = 10^6 * [1 0 -5 0 6];
            testCase.verifyEqual(ratQ.countRootsIn(p, 0, 2), 2, ...
                'sqrt(2) and sqrt(3) both lie in (0,2]');
            testCase.verifyEqual(ratQ.countRootsIn(p, 0, 1.5), 1, ...
                'only sqrt(2) lies in (0,1.5]');
            testCase.verifyEqual(ratQ.countRootsIn(p, -2, 2), 4, 'and all four are in (-2,2]');
        end

        function chkREFUSESAnythingThatIsNotAnExactInteger(testCase)
        % The single thing standing between an exact answer and a plausible wrong one.
            testCase.verifyError(@() ratQ.chk(1.5), 'ratQ:notInteger');
            testCase.verifyError(@() ratQ.chk(Inf), 'ratQ:notFinite');
            testCase.verifyError(@() ratQ.chk(NaN), 'ratQ:notFinite');
            testCase.verifyError(@() ratQ.chk(2^53), 'ratQ:overflow');
            testCase.verifyEqual(ratQ.chk(2^53 - 1), 2^53 - 1, 'and one below it is fine');
        end

        % ---- quadSignOnPoly: filter 2 of the bridge assembly ---------------------------------

        function quadSignOnPolyAnswersTheThreeOutcomesOnKnownGeometry(testCase)
        % The contract, on cases whose answer is arithmetic rather than a stored value.
            sq = [1 0 0; 0 1 0; -1 0 1; 0 -1 1];               % the unit square
            testCase.verifyEqual(ratQ.quadSignOnPoly([1 0 1 0 0 1], sq), 1, ...
                'x^2 + y^2 + 1 is positive everywhere, so certainly on the square');
            testCase.verifyEqual(ratQ.quadSignOnPoly([-1 0 -1 0 0 -1], sq), -1, ...
                'and its negative is negative everywhere');
            testCase.verifyEqual(ratQ.quadSignOnPoly([1 0 1 0 0 -1], sq), 0, ...
                ['x^2 + y^2 - 1 is negative at the origin and positive at (1,1), both corners ' ...
                 'of the square, so it must be reported as changing sign']);
            testCase.verifyEqual(ratQ.quadSignOnPoly([0 0 0 1 -1 0], sq), 0, ...
                'x - y changes sign across the square''s diagonal');
        end

        function quadSignOnPolyIsNotFooledByAMinimumOFF_theVertices(testCase)
        % THE CASE A VERTEX-ONLY TEST GETS WRONG, and the reason the facet stationary points are
        % enumerated at all. On the square [0,1]^2 the quadratic (x - 1/2)^2 - 1/16 is positive at
        % all four corners (1/4 - 1/16 = 3/16) and negative along the whole line x = 1/2. Anything
        % that inspected only the vertices would answer +1 and discard a genuine sign change.
        %
        % Written as 16x^2 - 16x + 3 so the coefficients stay integral.
            sq = [1 0 0; 0 1 0; -1 0 1; 0 -1 1];
            q  = [16 0 0 -16 0 3];
            testCase.verifyEqual(ratQ.quadSignOnPoly(q, sq), 0);
            testCase.verifyEqual(sign(16*(1/2)^2 - 16*(1/2) + 3), -1, ...
                'the interior value really is negative, so 0 is the right answer and not a decline');
        end

        function quadSignOnPolySeesUNBOUNDEDNESSRatherThanOnlyTheNearbyPoints(testCase)
        % The recession analysis, which is what separates this from sampling. On the first quadrant
        % the affine form 1 - x is positive near the origin and negative far along the x axis; a
        % routine that enumerated only the vertices of an unbounded cell would answer +1.
            quad = [1 0 0; 0 1 0];                              % x >= 0, y >= 0
            testCase.verifyEqual(ratQ.quadSignOnPoly([0 0 0 -1 0 1], quad), 0, ...
                '1 - x goes negative far along a recession direction');
            testCase.verifyEqual(ratQ.quadSignOnPoly([0 0 0 1 0 1], quad), 1, ...
                'while 1 + x stays positive on the whole quadrant, including at infinity');
        end

        function quadSignOnPolyNEVERClaimsOneSignWhenTheQuadraticChangesSign(testCase)
        % THE PROPERTY THAT MATTERS, asserted rather than sampled for a value.
        %
        % The two errors this routine can make are not symmetric. Answering 0 on a quadratic that
        % is really one-signed costs one elimination in the bridge assembly. Answering +1 or -1 on
        % one that changes sign DISCARDS A BRIDGE FACE, and the envelope is then wrong with no
        % symptom. So the test asserts only the unsafe direction, over random integer cases, and
        % says nothing about how often the safe direction is taken -- that is a rate, not a
        % contract, and .claude/quadSignOracle.m reports it.
        %
        % Small coefficients on purpose: they are what make the degenerate configurations the
        % guards exist for -- an exactly flat recession direction, exactly parallel constraints --
        % actually occur. A wide random range would test only the generic case.
            rng(2026);
            for it = 1:400
                c = randi([-4 4], 1, 6);
                k = randi([2 5]);
                P = randi([-3 3], k, 3);
                P(all(P(:,1:2) == 0, 2), 1) = 1;
                if ~ratQ.feasible2(P, true), continue, end
                s = ratQ.quadSignOnPoly(c, P);
                if s == 0, continue, end
                X = ratQTest.pointsIn(P);
                if isempty(X), continue, end
                v = c(1)*X(:,1).^2 + c(2)*X(:,1).*X(:,2) + c(3)*X(:,2).^2 + ...
                    c(4)*X(:,1)      + c(5)*X(:,2)        + c(6);
                tol = 1e-9 * max(1, max(abs(v)));
                if s > 0
                    testCase.verifyGreaterThanOrEqual(min(v), -tol, ...
                        sprintf('claimed q >= 0 on P for c = [%s]', num2str(c)));
                else
                    testCase.verifyLessThanOrEqual(max(v), tol, ...
                        sprintf('claimed q <= 0 on P for c = [%s]', num2str(c)));
                end
            end
        end

        function quadSignOnPolyIsNeverWeakerThanTheWholePlaneTest(testCase)
        % A REGRESSION GUARD WITH A HISTORY. One sign over the whole plane implies one sign on any
        % subset, so conicSign is a sound answer wherever it fires -- but the polyhedral route
        % DECLINES on degenerate configurations, and before the fallback was added it lost pairs
        % that the free test settled (DECISIONS 2026-09-06). A filter that is sometimes weaker than
        % the test it replaces is a regression however much stronger it is on average.
            rng(7);
            for it = 1:300
                c = randi([-4 4], 1, 6);
                if all(c == 0), continue, end
                P = randi([-3 3], randi([2 5]), 3);
                P(all(P(:,1:2) == 0, 2), 1) = 1;
                if ~ratQ.feasible2(P, true), continue, end
                sp = ratQ.conicSign(c);
                if sp ~= 0
                    testCase.verifyEqual(ratQ.quadSignOnPoly(c, P), sp, ...
                        sprintf('the plane says %d for c = [%s]; the cell must not say less', ...
                                sp, num2str(c)));
                end
            end
        end

        function quadSignOnPolySaysWhetherTheSignIsSTRICTBecauseATANGENTIALTieIsStillATie(testCase)
        % THE DEFECT THE END-TO-END RUN CAUGHT, pinned so it cannot come back.
        %
        % `+1` from this routine means `q >= 0`, which ALLOWS q = 0 somewhere. Two conjugate cells
        % whose difference touches zero without crossing are TANGENTIALLY tied, and a tangential
        % tie carries a bridge face exactly as a transversal one does. The bridge assembly used
        % `s ~= 0` as its prune test and therefore dropped a real bridge -- cell pair (5,5) of the
        % diagonal-kink triangle pair, where the definition's own maximiser ties. Only the second
        % output licenses a prune.
            sq = [1 0 0; 0 1 0; -1 0 1; 0 -1 1];               % the unit square
            [s, strict] = ratQ.quadSignOnPoly([1 0 1 0 0 1], sq);
            testCase.verifyEqual(s, 1);
            testCase.verifyTrue(strict, 'x^2 + y^2 + 1 never reaches zero');

            [s, strict] = ratQ.quadSignOnPoly([1 0 1 0 0 0], sq);
            testCase.verifyEqual(s, 1, 'x^2 + y^2 is nonnegative on the square');
            testCase.verifyFalse(strict, ...
                'but it VANISHES at the corner (0,0), so the sign is not strict and a tie is possible');

            [s, strict] = ratQ.quadSignOnPoly([1 0 1 0 0 -1], sq);
            testCase.verifyEqual(s, 0);
            testCase.verifyFalse(strict, 'a sign change is never strict');
        end

        function conicDefiniteIsTheLEADINGMinorTestNotTheSemidefiniteOne(testCase)
        % Two different questions with two different criteria, and using one for the other is the
        % classic error. Semidefinite needs EVERY principal minor (ratQ.conicSign); definite needs
        % only the LEADING ones (Sylvester). A parabola touching the axis is PSD and does vanish,
        % so it must not be called definite.
            testCase.verifyTrue(ratQ.conicDefinite([1 0 1 0 0 1]), 'x^2 + y^2 + 1 is definite');
            testCase.verifyTrue(ratQ.conicDefinite([-1 0 -1 0 0 -1]), 'and so is its negative');
            testCase.verifyFalse(ratQ.conicDefinite([1 0 0 0 -1 0]), ...
                'y = x^2 is PSD as a form and vanishes along the parabola, so it is NOT definite');
            testCase.verifyFalse(ratQ.conicDefinite([1 0 1 0 0 -1]), ...
                'and a form that changes sign is not definite either');
        end

        function conicSignMustNotBeFedThroughConicBecauseConicNormalisesTheSIGN(testCase)
        % THE DEFECT THE ORACLE CAUGHT, pinned so it cannot come back.
        %
        % ratQ.conic returns the canonical spelling of a CURVE, and a curve has no scale, so it
        % makes the leading nonzero coefficient positive -- negating half of all inputs. For a
        % curve that is the identity; for a SIGN it is the opposite answer. -(y-1)^2 is <= 0
        % everywhere, and conic() turns it into (y-1)^2.
            c = [0 0 -1 0 2 -1];                                % -(y - 1)^2
            testCase.verifyEqual(ratQ.conicSign(c), -1, 'the form itself is nonpositive');
            testCase.verifyEqual(ratQ.conicSign(ratQ.conic(c)), 1, ...
                'and conic() flips it, which is correct for a curve and wrong for a sign');
            halfPlane = [-1 -1 2; -3 -3 0];
            testCase.verifyEqual(ratQ.quadSignOnPoly(c, halfPlane), -1, ...
                'so quadSignOnPoly must report the nonpositive sign on this half-plane');
        end
    end

    methods (Static)
        function c = face(c5, c6, c7, c8, c9, c10)
        % A quadratic face in CCA2's 10-wide weighted cubic basis (RatCon.m's `f`):
        %     value = c5*x^2/2 + c6*xy + c7*y^2/2 + c8*x + c9*y + c10
            c = [0 0 0 0, c5, c6, c7, c8, c9, c10];
        end

        function s = symCubic(c, x, y)
        % The FULL 10-wide weighted basis of RatCon.m's `f`, as a symbolic expression:
        %   c1 x^3/6 + c2 x^2 y/2 + c3 x y^2/2 + c4 y^3/6
        % + c5 x^2/2 + c6 x y + c7 y^2/2 + c8 x + c9 y + c10
        % symOf below is the quadratic special case, kept because most tests only need it.
            s = c(1)*x^3/6 + c(2)*x^2*y/2 + c(3)*x*y^2/2 + c(4)*y^3/6 ...
              + c(5)*x^2/2 + c(6)*x*y + c(7)*y^2/2 + c(8)*x + c(9)*y + c(10);
        end

        function X = pointsIn(P)
        % Points of the polyhedron P for a TEST oracle: a rejection sample from a box, plus points
        % pushed far along each recession direction. The far ones are the point of it -- a
        % quadratic that only dips negative at distance is exactly what a box sample calls
        % one-signed, and exactly what the recession analysis in quadSignOnPoly exists to catch.
            X = (rand(4000, 2) * 2 - 1) * 6;
            X = X(all(X * P(:,1:2).' + P(:,3).' >= 0, 2), :);
            D = zeros(0,2);
            for i = 1:size(P,1)
                for u = [[P(i,2), -P(i,1)]; [-P(i,2), P(i,1)]].'
                    if any(u ~= 0) && all(P(:,1:2) * u >= 0), D = [D; u.']; end %#ok<AGROW>
                end
            end
            if isempty(X) || isempty(D), return, end
            base = X(1:min(end,10), :);
            for t = [1e2 1e4 1e6]
                for i = 1:size(D,1), X = [X; base + t * D(i,:)]; end %#ok<AGROW>
            end
        end

        function s = symOf(c, x, y)
            s = c(5)*x^2/2 + c(6)*x*y + c(7)*y^2/2 + c(8)*x + c(9)*y + c(10);
        end



        function aDEGENERATEConicNamesItsTwoLinesAndTheSIGNIsExact(testCase)
        % THE BRANCH ON AN EDGE, which is the edge's counterpart of AlgCon.FB on a face. Neither
        % line of a degenerate conic is rational; the PAIR is, and the branch NAMES one of them by
        % the canonical order of the slope quadratic's roots. The sign at a rational point is then
        % the sign of a linear polynomial AT A ROOT, which signAt decides -- no square root is ever
        % formed and nothing is compared in floating point.
        %
        % It exists because the convex envelope's bridge face is cut by a segment whose endpoints
        % are quadratic irrationals (DECISIONS.md 2026-09-17).
            % (x - y)(x + y - 2) = x^2 - y^2 - 2x + 2y: lines y = x and y = -x + 2, apex (1,1)
            C = ratQ.conic([1 0 -1 -2 2 0]);
            info = ratQ.pairLineData(C);
            testCase.verifyEqual(info.apex, [1; 1; 1], 'the apex is RATIONAL, though the lines are');
            testCase.verifyFalse(info.vertical);
            % branch 1 is the SMALLER slope: -1, i.e. y = -x + 2
            testCase.verifyEqual(ratQ.pairLineEval(C, 1, [3 0; 0 3]), [1; 1], 'AbsTol', 1e-12);
            testCase.verifyEqual(ratQ.pairLineEval(C, 2, [3 0; 0 3]), [-3; 3], 'AbsTol', 1e-12);

            % THE EXACT SIGN AGREES WITH THE CACHE, everywhere off the lines. Two routes that share
            % no arithmetic: signAt on the slope quadratic, and a double evaluation.
            bad = 0;  onLine = 0;
            for a = -4:4
                for b = -4:4
                    for br = 1:2
                        s = ratQ.pairLineSign(C, br, [a; b], 1);
                        v = ratQ.pairLineEval(C, br, [a b]);
                        if s == 0
                            onLine = onLine + 1;
                            if abs(v) > 1e-12, bad = bad + 1; end
                        elseif s ~= sign(v)
                            bad = bad + 1;
                        end
                    end
                end
            end
            testCase.verifyEqual(bad, 0, 'the exact sign and the cache must agree');
            testCase.verifyGreaterThan(onLine, 0, 'and points ON a line must be seen, as sign 0');

            % A VERTICAL line has no finite slope and is branch 2 by convention:
            % x(x - y) = x^2 - xy, lines x = 0 and y = x
            Cv = ratQ.conic([1 -1 0 0 0 0]);
            iv = ratQ.pairLineData(Cv);
            testCase.verifyTrue(iv.vertical);
            testCase.verifyEqual(ratQ.pairLineSign(Cv, 2, [3; 7], 1), 1, 'x = 0, to its right');
            testCase.verifyEqual(ratQ.pairLineSign(Cv, 2, [-3; 7], 1), -1);
            testCase.verifyEqual(ratQ.pairLineSign(Cv, 1, [0; 5], 1), 1, 'y = x, above it');
        end

        function aPARALLELPairNamesItsTwoLinesByTHEIRLEVELSAscending(testCase)
        % [fast] A corner face's wedge has one edge whose NORMAL is rational and whose OFFSET is not,
        % so its conjugate pair is a pair of PARALLEL lines -- no apex, no slopes to order by. The
        % branches are then the two LEVELS of the common linear form, ascending, which is the same
        % convention one dimension lower. Asserted against the definition, not a stored value:
        % branch b is the line { t = t_b }, so the exact sign must be the sign of t(x) - t_b.
            C = ratQ.conic([1 2 1 0 0 -2]);                 % (x + y)^2 = 2: x + y = -sqrt2, +sqrt2
            info = ratQ.pairLineData(C);
            testCase.verifyEqual(info.kind, 'parallel');
            lv = sort(roots(info.quad(:).'));               % the levels of the form, ascending
            testCase.verifyEqual(numel(lv), 2);

            bad = 0;  onLine = 0;
            for a = -3:3
                for b = -3:3
                    for br = [1 2]
                        t = info.form(1)*a + info.form(2)*b;
                        s = ratQ.pairLineSign(C, br, [a; b], 1);
                        testCase.verifyEqual(s, ratQ.pairLineSign(C, br, [2*a; 2*b], 2), ...
                                             'the sign must not depend on how the point is spelled');
                        if abs(t - lv(br)) < 1e-9
                            onLine = onLine + 1;
                            if s ~= 0, bad = bad + 1; end
                        elseif s ~= sign(t - lv(br))
                            bad = bad + 1;
                        end
                        v = ratQ.pairLineEval(C, br, [a b]);
                        if sign(v) ~= s && abs(v) > 1e-9, bad = bad + 1; end
                    end
                end
            end
            testCase.verifyEqual(bad, 0, 'the exact sign must be the sign of t(x) - t_b');
            testCase.verifyGreaterThan(onLine, 0, 'and the levels here are irrational, so none');

            % The form may be y alone, when the quadratic part has no x at all: y^2 - 3.
            Cy = ratQ.conic([0 0 1 0 0 -3]);
            iy = ratQ.pairLineData(Cy);
            testCase.verifyEqual(iy.kind, 'parallel');
            testCase.verifyEqual(ratQ.pairLineSign(Cy, 2, [0; 5], 1), 1, 'above y = +sqrt3');
            testCase.verifyEqual(ratQ.pairLineSign(Cy, 2, [0; 1], 1), -1, 'below it');
            testCase.verifyEqual(ratQ.pairLineSign(Cy, 1, [0; 1], 1), 1, 'but above y = -sqrt3');

            % A PARABOLA has the same square as its quadratic part and is NOT a pair of lines, so it
            % is refused rather than answered with a branch of something that is not there.
            testCase.verifyError(@() ratQ.pairLineData(ratQ.conic([1 0 0 0 -1 0])), ...
                'ratQ:notALinePair');                       % x^2 - y
            testCase.verifyError(@() ratQ.pairLineData(ratQ.conic([1 2 1 0 0 0])), ...
                'ratQ:notALinePair');                       % (x + y)^2: ONE line, twice
        end

        function aCONICThatIsNOTTwoCrossingLinesNamesNoBranch(testCase)
        % The refusals, because a branch on the wrong kind of conic would be read as the whole curve
        % by any consumer that forgot to look -- so it is refused where it is asked for.
            testCase.verifyError(@() ratQ.pairLineData(ratQ.conic([1 0 1 0 0 -1])), ...
                'ratQ:notALinePair');                       % a circle: definite quadratic part
            testCase.verifyError(@() ratQ.pairLineData(ratQ.conic([1 0 -1 0 0 -1])), ...
                'ratQ:notDegenerate');                      % a hyperbola: indefinite but NOT a pair
            testCase.verifyError(@() ratQ.pairLineSign(ratQ.conic([1 0 -1 -2 2 0]), 3, [0;0], 1), ...
                'ratQ:branch');
        end
    end
end

function r = addAsRow(n1, d1, n2, d2)
    [n, d] = ratQ.add(n1, d1, n2, d2);
    r = [n, d];
end

function verifyExactOrRefused(testCase, f, want, name)
    try
        got = f();
    catch ME
        testCase.verifyEqual(ME.identifier, 'ratQ:overflow', ...
            sprintf('%s: the only acceptable refusal is ratQ:overflow; got %s', name, ME.identifier));
        return
    end
    testCase.verifyEqual(got, want, sprintf('%s: a value came back, so it must be the EXACT one', name));
end
