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
        % CHANGED 2026-09-03: this used to hand in the unit square, which Case B now covers. The
        % assertion is unaltered -- an uncovered case must be a NAMED refusal and never a silent
        % fallback into the symbolic engine -- so the fixture moves to one that is genuinely
        % uncovered: a MULTI-FACE input, which needs Step 3 (the pointwise max across pieces).
            V = [0 0; 1 0; 1 1; 0 1];
            E = [1 2 1; 2 3 1; 3 4 1; 4 1 1; 1 3 1];
            f = QuaPol(V, E, [0 0 0 0 1 0 1 0 0 0; 0 0 0 0 2 0 2 0 0 0], ...
                       [1 0; 1 0; 2 0; 2 0; 1 2]);
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
        % ---- Case B: a strictly convex quadratic on a bounded convex polygon --------------------

        function theSquareConjugateMatchesTheSeparableClosedForm(testCase)
        % (x^2+y^2)/2 over [0,1]^2 SEPARATES, so f*(s) = g(s1) + g(s2) with g the 1-D conjugate
        %       g(t) = 0 for t <= 0,  t^2/2 for 0 <= t <= 1,  t - 1/2 for t >= 1.
        % That closed form is written down here independently of anything conjQ computes, which is
        % what makes it an oracle rather than a restatement.
            g = conjQ(conjQTest.unitSquareEnergy());
            gg = @(t) (t<=0).*0 + (t>0 & t<1).*(t.^2/2) + (t>=1).*(t-0.5);

            rng(20260903);
            S = [randn(300,2)*1.5; 0 0; 1 1; 0.5 0.5; -2 3; 2 2; 1 0; 0 1];
            [got, idx] = g.eval(S);
            testCase.verifyEqual(got, gg(S(:,1)) + gg(S(:,2)), 'RelTol', 1e-12, 'AbsTol', 1e-12);
            testCase.verifyTrue(all(idx > 0), ...
                'the cells must COVER the plane: dom f* is R^2 for a bounded domain');
        end

        function theSquareConjugateHasTheSubdivisionTheMathematicsPredicts(testCase)
        % Nine cells -- one per vertex (4), one per edge (4), one interior -- and the cell
        % boundaries are exactly the four lines s1 = 0, s1 = 1, s2 = 0, s2 = 1.
            g = conjQ(conjQTest.unitSquareEnergy());
            testCase.verifyEqual(g.nf, 9, '4 vertex cells + 4 edge cells + 1 interior');
            testCase.verifyEqual(g.ne, 4, 'and only four distinct bounding lines');
            for j = 1:g.ne
                testCase.verifyEqual(g.edgeKind(j), 'line');
            end
            want = sortrows([ratQ.conic([0 0 0 1 0  0]);      % s1 = 0
                             ratQ.conic([0 0 0 1 0 -1]);      % s1 = 1
                             ratQ.conic([0 0 0 0 1  0]);      % s2 = 0
                             ratQ.conic([0 0 0 0 1 -1])]);    % s2 = 1
            testCase.verifyEqual(sortrows(g.EcQ), want);
            testCase.verifyEqual(g.nv, 4, 'the four dual corners where those lines meet');
        end

        function aFacetSHAREDByTwoCellsIsONEEdgeIndexNotTwoSpellings(testCase)
        % THE property the exact representation exists for. DECISIONS.md 2026-08-17 measured the
        % alternative: two cells reached one facet by different routes, the doubles differed by an
        % ULP, region.merge could not see that they matched, and Step 3's cell count grew without
        % bound -- 57 cells carrying 10 distinct functions, 4 merges out of 612 attempts.
        %
        % Here the constraint rows are canonical integers, so a shared facet is the SAME INDEX.
        % Asserted structurally: the cells reference far more constraints than there are distinct
        % lines, which can only happen through sharing.
            g = conjQ(conjQTest.unitSquareEnergy());
            total = 0;
            for k = 1:g.nf, total = total + size(g.FC{k},1); end
            testCase.verifyGreaterThan(total, 3 * g.ne, ...
                'the cells must reuse the same line indices, not each carry their own copy');
            for k = 1:g.nf
                testCase.verifyTrue(all(g.FC{k}(:,1) >= 1 & g.FC{k}(:,1) <= g.ne));
            end
        end

        function conjugatesOverRandomPolygonsMatchAnIndependentMaximisation(testCase)
        % Differential sweep. The oracle enumerates the candidate maximisers of a concave objective
        % over a polytope -- the unconstrained stationary point when interior, each edge's clamped
        % 1-D stationary point, and the vertices -- which is a DIFFERENT computation from conjQ's
        % cell subdivision: no normal cones, no multipliers, no arrangement.
            rng(20260903);
            for c = 1:12
                [V, E, F] = conjQTest.randomConvexPolygon();
                [H, L, k0] = conjQTest.randomStrictlyConvexQ();
                f = QuaPol(V, E, [0 0 0 0, H(1,1), H(1,2), H(2,2), L(1), L(2), k0], F);
                testCase.assumeTrue(f.isExact());

                g = conjQ(f);
                S = [randn(60,2)*3; 0 0];
                [got, idx] = g.eval(S);
                testCase.verifyTrue(all(idx > 0), sprintf('case %d: a point fell in no cell', c));

                want = zeros(size(S,1),1);
                for i = 1:size(S,1)
                    want(i) = conjQTest.supOverPolygon(S(i,:).', V, H, L, k0);
                end
                testCase.verifyEqual(got, want, 'RelTol', 1e-9, 'AbsTol', 1e-9, ...
                    sprintf('case %d: conjQ disagrees with the independent maximisation', c));
            end
        end

        function theConjugateOverAPolygonIsConvex(testCase)
        % f* is convex whatever f is -- a definition-level property no cell-by-cell check implies,
        % and one a wrong cell boundary breaks immediately.
            g = conjQ(conjQTest.unitSquareEnergy());
            rng(20260903);
            A = randn(200,2)*2;  B = randn(200,2)*2;  M = (A+B)/2;
            testCase.verifyLessThanOrEqual(g.eval(M), (g.eval(A) + g.eval(B))/2 + 1e-12, ...
                'the midpoint value must not exceed the average of the endpoints');
        end

        function theCellFunctionsAreExactRationalsWithSmallDenominators(testCase)
        % Not a golden value: the assertion is that the coefficients ARE rationals stored as
        % integers, and that the denominators have not blown up -- which is what would happen if
        % the routine multiplied denominators instead of taking the lcm.
            g = conjQ(conjQTest.unitSquareEnergy());
            testCase.verifyTrue(all(g.fD > 0));
            testCase.verifyEqual(g.fN, round(g.fN), 'AbsTol', 0, 'numerators are integers');
            testCase.verifyLessThan(max(g.fD), 100, ...
                'the unit square carries halves, not products of every intermediate denominator');
        end

        % ---- Case B refuses what it does not cover, by name ---------------------------------

        function anUnboundedPolygonIsRefusedByName(testCase)
            V = [0 0; 1 0; 0 1];
            f = QuaPol(V, [1 2 1; 2 3 0; 3 1 0], [0 0 0 0 1 0 1 0 0 0], [1 0; 1 0; 1 0]);
            testCase.verifyError(@() conjQ(f), 'PLQ:conjQ:notImplemented');
        end

        function aSemidefiniteSingularQuadraticOnAPolygonIsRefusedByName(testCase)
        % IDENTIFIER UPDATED 2026-09-03 when Case C landed: the refusal moved UP to the dispatch,
        % which now classifies H exactly and names all three outcomes (strictly convex, concave or
        % affine, neither), so it raises notImplemented rather than notStrictlyConvex. The
        % assertion this test exists for is unchanged -- an uncovered case is a NAMED refusal and
        % never a silent fallback.
        %
        % A rank-one PSD Q is genuinely a third case: the interior cell degenerates, because the
        % unconstrained sup is finite only on a measure-zero set.
            V = [0 0; 1 0; 1 1; 0 1];  E = [1 2 1; 2 3 1; 3 4 1; 4 1 1];  F = [1 0;1 0;1 0;1 0];
            f = QuaPol(V, E, [0 0 0 0 1 1 1 0 0 0], F);      % rank-one PSD
            testCase.verifyError(@() conjQ(f), 'PLQ:conjQ:notImplemented');
        end
        % ---- Case C: a concave or affine quadratic on a bounded convex polygon ------------------

        function theConcaveConjugateIsTheMaxOverTheVerticesAndNothingElse(testCase)
        % If q is concave then <s,x> - q(x) is CONVEX in x, so its max over a polytope sits at an
        % extreme point and f* is the max of the affine functions <s,v> - q(v). The oracle is that
        % max, written directly -- no cells, no normal fan.
            V = [0 0; 1 0; 1 1; 0 1];  E = [1 2 1; 2 3 1; 3 4 1; 4 1 1];  F = [1 0;1 0;1 0;1 0];
            H = [-2 -1; -1 -2];  L = [1; -2];  k0 = 3;
            f = QuaPol(V, E, [0 0 0 0, H(1,1), H(1,2), H(2,2), L(1), L(2), k0], F);
            g = conjQ(f);

            rng(20260903);
            S = [randn(300,2)*2.5; 0 0];
            want = conjQTest.maxOverVertices(S, V, H, L, k0);
            [got, idx] = g.eval(S);
            testCase.verifyEqual(got, want, 'RelTol', 1e-12, 'AbsTol', 1e-12);
            testCase.verifyTrue(all(idx > 0), 'the normal fan covers the plane');
            testCase.verifyTrue(all(g.fN(:,5) == 0 & g.fN(:,6) == 0 & g.fN(:,7) == 0), ...
                'every cell of a concave conjugate is AFFINE -- no curvature can appear');
        end

        function anAffineFunctionConjugatesToTheSupportFunctionOfItsDomain(testCase)
        % q affine is the same statement with a zero Hessian, and then f* is the polygon's support
        % function shifted -- the elementary case, and the one a sign error shows up in first.
            V = [0 0; 2 0; 2 1; 0 1];  E = [1 2 1; 2 3 1; 3 4 1; 4 1 1];  F = [1 0;1 0;1 0;1 0];
            L = [1; 2];  k0 = -1;
            f = QuaPol(V, E, [0 0 0 0, 0, 0, 0, L(1), L(2), k0], F);
            g = conjQ(f);

            rng(20260903);
            S = [randn(200,2)*2; 0 0];
            want = conjQTest.maxOverVertices(S, V, zeros(2), L, k0);
            testCase.verifyEqual(g.eval(S), want, 'RelTol', 1e-12, 'AbsTol', 1e-12);
            testCase.verifyEqual(g.nf, 4, 'each of the four corners wins on its own normal cone');
        end

        function aDominatedVertexContributesNoCellAtAll(testCase)
        % WHERE THE EXACT FEASIBILITY TEST EARNS ITS KEEP. Add a vertex at the midpoint of an edge
        % and give q strict concavity: q(midpoint) then EXCEEDS the average of the two endpoints,
        % so the lifted point sits strictly above the chord and that vertex never attains the max.
        % Its cell is EMPTY -- not small, empty -- and must not appear in the mesh.
            V = [0 0; 1/2 0; 1 0; 1 1; 0 1];
            E = [1 2 1; 2 3 1; 3 4 1; 4 5 1; 5 1 1];
            F = [ones(5,1), zeros(5,1)];
            f = QuaPol(V, E, [0 0 0 0, -2, 0, -2, 0, 0, 0], F);   % q = -(x^2+y^2)
            g = conjQ(f);

            testCase.verifyEqual(g.nf, 4, ...
                'the midpoint vertex is dominated, so four cells remain, not five');
            rng(20260903);
            S = [randn(200,2)*2; 0 0];
            want = conjQTest.maxOverVertices(S, V, [-2 0; 0 -2], [0;0], 0);
            testCase.verifyEqual(g.eval(S), want, 'RelTol', 1e-12, 'AbsTol', 1e-12, ...
                'and dropping it changes no value -- it never attained the max anywhere');
        end

        function concaveConjugatesOverRandomPolygonsMatchTheVertexMax(testCase)
            rng(20260903);
            for c = 1:10
                [V, E, F] = conjQTest.randomConvexPolygon();
                H = randi([0 4], 2, 2);  H = -(H + H.');
                if ~(H(1,1) <= 0 && H(2,2) <= 0 && (H(1,1)*H(2,2)-H(1,2)^2) >= 0), continue, end
                L = randi([-4 4], 2, 1);  k0 = randi([-4 4]);
                f = QuaPol(V, E, [0 0 0 0, H(1,1), H(1,2), H(2,2), L(1), L(2), k0], F);
                testCase.assumeTrue(f.isExact());
                g = conjQ(f);
                S = [randn(60,2)*3; 0 0];
                [got, idx] = g.eval(S);
                testCase.verifyTrue(all(idx > 0), sprintf('case %d: uncovered dual point', c));
                testCase.verifyEqual(got, conjQTest.maxOverVertices(S, V, H, L, k0), ...
                    'RelTol', 1e-9, 'AbsTol', 1e-9, sprintf('case %d', c));
            end
        end

        function anIndefiniteQuadraticOnAPolygonIsRefusedByName(testCase)
        % Neither branch applies, and the refusal must say which case is missing rather than
        % falling through to something that happens to run.
            V = [0 0; 1 0; 1 1; 0 1];  E = [1 2 1; 2 3 1; 3 4 1; 4 1 1];  F = [1 0;1 0;1 0;1 0];
            f = QuaPol(V, E, [0 0 0 0, 1, 0, -1, 0, 0, 0], F);      % x^2/2 - y^2/2
            testCase.verifyError(@() conjQ(f), 'PLQ:conjQ:notImplemented');
        end
    end

    methods (Static)

        function w = maxOverVertices(S, V, H, L, k0)
        % THE ORACLE for Case C: max over the polygon's VERTICES of <s,v> - q(v), written directly.
        % For a concave or affine q that IS the conjugate, because <s,x> - q(x) is then convex in x
        % and a convex function attains its maximum over a polytope at an extreme point.
            w = -inf(size(S,1), 1);
            for j = 1:size(V,1)
                v = V(j,:).';
                qv = 0.5 * v.' * H * v + L.' * v + k0;
                w = max(w, S * v - qv);
            end
        end

        function f = unitSquareEnergy()
        % (x^2 + y^2)/2 over [0,1]^2
            f = QuaPol([0 0; 1 0; 1 1; 0 1], [1 2 1; 2 3 1; 3 4 1; 4 1 1], ...
                       [0 0 0 0 1 0 1 0 0 0], [1 0; 1 0; 1 0; 1 0]);
        end

        function v = supOverPolygon(s, V, H, L, k0)
        % THE ORACLE. max over the polygon of <s,x> - q(x), by enumerating the candidate
        % maximisers of a concave objective over a polytope: the unconstrained stationary point
        % when it lies inside, each edge's clamped 1-D stationary point, and the vertices. One of
        % those is always the maximiser, so this is exact up to double rounding -- and it shares no
        % machinery with conjQ, which is the point.
            q   = @(x) 0.5 * x.' * H * x + L.' * x + k0;
            obj = @(x) s.' * x - q(x);
            cand = V.';
            xs = H \ (s - L);
            if conjQTest.inPolygon(xs.', V), cand = [cand, xs]; end
            m = size(V,1);
            for j = 1:m
                a = V(j,:).';  b = V(mod(j,m)+1,:).';  d = b - a;
                t = min(1, max(0, (s - (H*a + L)).' * d / (d.' * H * d)));
                cand = [cand, a + t*d]; %#ok<AGROW>
            end
            v = -inf;
            for i = 1:size(cand,2), v = max(v, obj(cand(:,i))); end
        end

        function tf = inPolygon(x, V)
            m = size(V,1);  ctr = mean(V,1);  tf = true;
            for j = 1:m
                a = V(j,:);  b = V(mod(j,m)+1,:);  d = b - a;  n = [d(2), -d(1)];
                if n * (ctr - a).' > 0, n = -n; end
                if n * (x - a).' > 1e-12, tf = false; return, end
            end
        end

        function [V, E, F] = randomConvexPolygon()
        % Vertices on a coarse half-integer grid, so they are exactly rational and QuaPol reads
        % them exactly -- an irrational vertex would make the whole object inexact and there would
        % be nothing to test.
            m  = randi([3 6]);
            th = sort(rand(1,m) * 2*pi);
            r  = randi([1 6], 1, m);
            V  = round([r'.*cos(th'), r'.*sin(th')] * 2) / 2;
            V  = uniquetol(V, 1e-9, 'ByRows', true);
            k  = convhull(V(:,1), V(:,2));
            V  = V(k(1:end-1), :);
            m  = size(V,1);
            E  = [(1:m).', [2:m, 1].', ones(m,1)];
            F  = [ones(m,1), zeros(m,1)];
        end

        function [H, L, k0] = randomStrictlyConvexQ()
            while true
                H = randi([-4 4], 2, 2);  H = H + H.';
                if H(1,1) > 0 && (H(1,1)*H(2,2) - H(1,2)^2) > 0, break, end
            end
            L = randi([-4 4], 2, 1);  k0 = randi([-4 4]);
        end
    end
end
