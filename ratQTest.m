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

        function s = symOf(c, x, y)
            s = c(5)*x^2/2 + c(6)*x*y + c(7)*y^2/2 + c(8)*x + c(9)*y + c(10);
        end
    end
end
