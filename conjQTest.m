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


        function aMultiFaceQuaPolConjugatesByFoldingItsPieces(testCase)
        % CHANGED TWICE, and the trail is the point. It first pinned the unit square as
        % unimplemented (Case B covered it), then a MULTI-FACE input as unimplemented -- and the
        % per-piece fold now covers that too. Rather than hunt for a third uncovered fixture, it
        % becomes the POSITIVE test of what replaced it: f is q_k on P_k, so
        %       f*(s) = max_k sup_{x in P_k} <s,x> - q_k(x),
        % and conjugation turns a union into a max. The remaining refusals are pinned by
        % anIndefinite... and aSemidefiniteSingular... below.
            V  = [0 0; 1 0; 1 1; 0 1];
            E  = [1 2 1; 2 3 1; 3 4 1; 4 1 1; 1 3 1];
            % F(5,:) = [2 1], NOT [1 2]. F(j,:) is [left, right] of edge j, and face 1 -- the
            % triangle BELOW the diagonal -- is on its RIGHT. Written the other way the mesh is
            % MALFORMED: P{1}'s signs then describe {y >= 0, y >= x, x <= 1}, an unbounded
            % region whose corners are not the endpoints of face 1's own edges, and QuaPol.eval
            % reports the wrong piece at interior points. This test passed anyway, because the
            % code and the oracle both read the face as the polygon its edges bound -- which is
            % what a user means and what a WELL-FORMED mesh says. Measured 2026-09-04.
            F  = [1 0; 1 0; 2 0; 2 0; 2 1];       % face 1 = edges 1,2,5; face 2 = edges 3,4,5
            H1 = [1 0; 0 1];  L1 = [0;0];  k1 = 0;
            H2 = [4 1; 1 3];  L2 = [-2;1]; k2 = 0;
            f  = QuaPol(V, E, [0 0 0 0, H1(1,1), H1(1,2), H1(2,2), L1(1), L1(2), k1; ...
                               0 0 0 0, H2(1,1), H2(1,2), H2(2,2), L2(1), L2(2), k2], F);
            testCase.verifyEqual(f.nf, 2);
            g = conjQ(f);

            T1 = [0 0; 1 0; 1 1];  T2 = [1 1; 0 1; 0 0];
            rng(20260903);
            S = [randn(150,2)*2; 0 0];
            want = zeros(size(S,1),1);
            for i = 1:size(S,1)
                want(i) = max(conjQTest.supOverPolygon(S(i,:).', T1, H1, L1, k1), ...
                              conjQTest.supOverPolygon(S(i,:).', T2, H2, L2, k2));
            end
            [got, idx] = g.eval(S);
            testCase.verifyTrue(all(idx > 0), 'the folded cells must cover the plane');
            testCase.verifyEqual(got, want, 'RelTol', 1e-9, 'AbsTol', 1e-9);
        end


        function aCONVEXQuadraticOnAnUNBOUNDEDConeIsComputedNotRefused(testCase)
        % REPLACES three refusal tests deleted 2026-09-04, and it is the gap the coverage table
        % named: q = (x^2+y^2)/2 on the first quadrant has a FINITE conjugate at every s, and the
        % routine used to refuse it. The quadrant makes the problem separable, so the answer is
        % written down independently: sup_{x>=0} t*x - x^2/2 is t^2/2 for t >= 0 and 0 below.
            f = conjQTest.wedge([0 0], [1 0], [0 1], [0 0 0 0 1 0 1 0 0 0]);
            g = conjQ(f);
            gg = @(t) (t >= 0) .* (t.^2/2);
            rng(20260904);
            S = [randn(300,2)*2; 0 0; 1 1; -1 -1; 2 -3];
            [got, idx] = g.eval(S);
            testCase.verifyTrue(all(idx > 0), 'dom f* is the whole plane here');
            testCase.verifyEqual(got, gg(S(:,1)) + gg(S(:,2)), 'RelTol', 1e-12, 'AbsTol', 1e-12);
        end

        function aNULLRECESSIONDIRECTIONGivesAHalfPlaneNotARefusal(testCase)
        % The second half of the same gap, and the reason it is a computation rather than a
        % refusal: along a recession direction of ZERO curvature the slope is
        % <s-L,d> - <Hd,x>, so the condition is <s-L,d> <= inf_P <Hd,x> -- a LINEAR PROGRAM over
        % the piece, hence the minimum over its vertices. For q = xy on the first quadrant and
        % d = (1,0) that infimum is 0, giving s1 <= 0, and by symmetry s2 <= 0; on that cone the
        % sup is 0, attained at the origin.
            f = conjQTest.wedge([0 0], [1 0], [0 1], [0 0 0 0 0 1 0 0 0 0]);
            g = conjQ(f);
            inCone  = [-1 -1; -3 -0.5; 0 0; -2 -7];
            outCone = [1 -1; -1 1; 2 2; 0.5 -3];
            testCase.verifyEqual(g.eval(inCone), zeros(4,1), 'AbsTol', 1e-12);
            testCase.verifyTrue(all(isinf(g.eval(outCone))));
        end


        function aDivergentSupReturnsThePlusInfinityFunctionRatherThanRaising(testCase)
        % REPLACES three refusal tests, 2026-09-04. They asserted that an empty dual domain was
        % REFUSED, on the belief that a QuaCon could not hold such a function. It can: an object
        % with ZERO faces matches nothing in eval and so answers +infinity everywhere, which IS
        % the function. The routines were declining to return an answer the type already stores.
        %
        % A strictly concave q on a cone runs to +infinity along every recession direction, so
        % f* is +infinity for every s.
            g = conjQ(conjQTest.wedge([0 0], [1 0], [0 1], [0 0 0 0, -2, 0, -2, 1, 1, 0]));
            testCase.verifyEqual(g.nf, 0, 'no faces is how "+infinity everywhere" is spelled');
            testCase.verifyTrue(all(isinf(g.eval([0 0; 1 2; -3 4; 1e6 -1e6]))));
            testCase.verifyTrue(g.isMeshed(), 'it is still a well-formed object');
        end

        function aFullPlaneQuadraticWithNegativeCurvatureAlsoConjugatesToPlusInfinity(testCase)
        % Three of the five non-PD full-plane cases are this, not a thin domain: if H has any
        % direction of negative curvature then <s,x> - q(x) rises without bound along it for EVERY
        % s. Negative definite, negative semidefinite-singular and indefinite all qualify.
            for H = {[-2 0; 0 -2], [1 0; 0 -1], [-1 0; 0 0]}
                g = conjQ(QuaPol([H{1}(1,1), H{1}(1,2), H{1}(2,2), 1, 2, 3]));
                testCase.verifyEqual(g.nf, 0, sprintf('H = %s', mat2str(H{1})));
                testCase.verifyTrue(all(isinf(g.eval([0 0; 3 4]))));
            end
        end

        function anAFFINEFunctionOnThePlaneConjugatesToAPOINTSupportedFunction(testCase)
        % THE THIN DUAL DOMAIN, first kind. q(x) = <a,x> + b on all of R^2 gives
        %       f*(s) = -b at s = a,  +infinity elsewhere,
        % so dom f* is a single POINT. It is stored with the H-form's EQUALITY side -- a side of 0
        % on a curve means ON it -- which is one value in a column that already existed rather
        % than a new class.
            a = [1; 2];  b = 3;
            g = conjQ(QuaPol([0 0 0 a(1) a(2) b]));
            testCase.verifyEqual(g.nf, 1);
            testCase.verifyEqual(g.ne, 2, 'two equalities cut the plane down to one point');
            testCase.verifyEqual(g.eval(a.'), -b, 'AbsTol', 1e-12);
            testCase.verifyTrue(all(isinf(g.eval([1.5 2; 1 2.5; 0 0; -3 7]))), ...
                'and +infinity at every other point');
            for k = 1:g.nf
                testCase.verifyTrue(all(g.FC{k}(:,2) == 0), 'both sides are EQUALITIES');
            end
        end

        function aPSDSINGULARQuadraticOnThePlaneConjugatesToALINESupportedFunction(testCase)
        % The second kind. With H = lambda*w*w' the objective splits: along w it is a concave
        % parabola, and along w-perp it is LINEAR, so the sup is +infinity unless s - L is parallel
        % to w. Hence dom f* is the LINE {u parallel to w} and there f* = (u.w)^2/(2 lambda) - kappa.
        % Written down independently of conjQ's own algebra below.
            rng(20260904);
            for t = 1:12
                w = randi([-3 3], 2, 1);
                if all(w == 0), continue, end
                lam = randi([1 4]);
                H = lam * (w * w.');
                L = randi([-3 3], 2, 1);  k0 = randi([-3 3]);
                g = conjQ(QuaPol([H(1,1), H(1,2), H(2,2), L(1), L(2), k0]));
                testCase.verifyEqual(g.nf, 1, sprintf('case %d', t));
                testCase.verifyEqual(g.ne, 1, 'one equality: a line');

                % on the line: s = L + r*w for any r
                for r = [-2 -0.5 0 1 3]
                    sOn = (L + r*w).';
                    want = (sOn*w - L.'*w)^2 / (2 * lam * (w.'*w)^2) * (w.'*w) - k0;
                    % (u.w)^2/(2 lambda |w|^4) * |w|^2 simplifies to (u.w)^2/(2 lambda |w|^2)
                    want = (r*(w.'*w))^2 / (2*lam*(w.'*w)^2) - k0;
                    testCase.verifyEqual(g.eval(sOn), want, 'RelTol', 1e-9, 'AbsTol', 1e-9, ...
                        sprintf('case %d on the line at r=%g', t, r));
                end
                % off the line
                off = (L + [w(2); -w(1)]).';
                testCase.verifyTrue(isinf(g.eval(off)), 'and +infinity off it');
            end
        end

        function everyRefusalLeftIsByDESIGNRatherThanAGap(testCase)
        % After the thin-domain work the coverage table has no gaps: what still raises is what the
        % toolbox deliberately does not accept, and both are contract rather than absence.
            V = [0 0; 1 0; 1 1; 0 1];  E = [1 2 1; 2 3 1; 3 4 1; 4 1 1];  F = [1 0;1 0;1 0;1 0];
            cubic = QuaPol(V, E, [1 0 0 0 1 0 1 0 0 0], F);
            testCase.verifyError(@() conjQ(cubic), ?MException, ...
                'a cubic numerator is out of scope by assertOperable');
            testCase.verifyError(@() conjQ(QuaPol([sqrt(2) 0 1 0 0 0])), 'PLQ:QuaPol:notExact', ...
                'and an inexact input is refused deliberately, not for want of an algorithm');
        end


        function aNEEDLEConjugatesToAnAffineFunctionOnTheWHOLEPlane(testCase)
        % A low-dimensional INPUT has a full-dimensional OUTPUT, which is why this needed no new
        % storage: the conjugate of the function that is c at p and +infinity elsewhere is
        % <s,p> - c, affine and finite everywhere.
            p = [2 -1];  c = 5;
            f = QuaPol(p, zeros(0,3), [0 0 0 0 0 0 0 0 0 c], zeros(0,2));
            testCase.verifyEqual(f.dom.dim, 0);
            g = conjQ(f);
            rng(20260904);
            S = [randn(200,2)*3; 0 0; 1 1];
            testCase.verifyEqual(g.eval(S), S*p.' - c, 'AbsTol', 1e-12);
            testCase.verifyEqual(g.nf, 1, 'one affine piece, and no boundary at all');
            testCase.verifyEqual(g.ne, 0);
        end

        function aSEGMENTConjugatesToItsOneDimensionalClampedMaximum(testCase)
        % The dimension-1 input. The sup over a segment of <s,x> - q(x) is a ONE-dimensional
        % problem in the segment parameter, so the candidates are the two endpoints and, where the
        % restriction is concave, its interior stationary point -- which is exactly caseD's
        % candidate set, computed here with no new machinery.
            f = QuaPol([0 0; 2 0], [1 2 1], [0 0 0 0 1 0 0 0 0 0], [0 0]);   % q = x^2/2 on [0,2]x{0}
            testCase.verifyEqual(f.dom.dim, 1);
            g = conjQ(f);
            rng(20260904);
            S = [randn(200,2)*2; -1 0; 0.5 0; 1 5; 3 -2; 2 0];
            want = zeros(size(S,1),1);
            for i = 1:size(S,1)
                t = min(2, max(0, S(i,1)));            % the clamped maximiser along the segment
                want(i) = S(i,1)*t - t^2/2;
            end
            testCase.verifyEqual(g.eval(S), want, 'RelTol', 1e-12, 'AbsTol', 1e-12);
        end

        function aNONCONVEXFaceIsSPLITAndConjugatedCorrectly(testCase)
        % WAS a refusal, now an answer. A piece is described here as an intersection of half-planes,
        % which cannot express a REFLEX corner -- conjugating an L-shaped face through that
        % description silently took the sup over the smaller set the half-planes cut out, too low
        % at 54 of 303 dual points. The fix is not a better description but a SPLIT: conj turns a
        % union into a MAX, so f* over P = union of P_i is max_i (q + I_{P_i})*, and triangulating
        % the face is exact rather than approximate. The triangulation is exact too -- ear clipping,
        % whose every decision is the sign of an integer cross product.
            V = [0 0; 2 0; 2 1; 1 1; 1 2; 0 2];          % an L, reflex at (1,1)
            m = size(V,1);
            E = [(1:m).', [2:m,1].', ones(m,1)];  F = [ones(m,1), zeros(m,1)];
            f = QuaPol(V, E, [0 0 0 0 1 0 1 0 0 0], F);
            [~, isConv] = f.orderEdges(1);
            testCase.verifyFalse(logical(isConv), 'the fixture really is non-convex');

            g = conjQ(f);

            % The oracle is a sup sampled over the TRUE L -- a LOWER bound, so f* below it is a
            % definite defect, and on a grid this fine the gap above it must be tiny as well.
            [X, Y] = meshgrid(linspace(0, 2, 501));
            P = [X(:), Y(:)];
            P = P((P(:,1) <= 1 & P(:,2) <= 2) | (P(:,2) <= 1 & P(:,1) <= 2), :);
            rng(20260904);
            S = [randn(120,2)*2.5; 0 0; 3 3; -3 -3];
            lo = zeros(size(S,1),1);
            for i = 1:size(S,1)
                lo(i) = max(P*S(i,:).' - 0.5*sum(P.^2, 2));
            end
            got = g.eval(S);
            testCase.verifyGreaterThanOrEqual(got, lo - 1e-9, ...
                'f* may never fall below a sup actually attained on the domain');
            testCase.verifyLessThan(max(got - lo), 1e-4, ...
                'nor exceed it by more than the sampling grid can explain');

            A = randn(200,2)*2;  B = randn(200,2)*2;
            testCase.verifyLessThanOrEqual(g.eval((A+B)/2), (g.eval(A) + g.eval(B))/2 + 1e-9, ...
                'and the fold of the pieces is still convex');
        end

        function aNonConvexUNBOUNDEDPieceIsStillRefusedByName(testCase)
        % The split needs a closed boundary to clip ears from, so a non-convex piece that also
        % recedes is a separate case -- cutting it would mean splitting along its recession
        % directions first. Refused by name rather than silently mis-described.
            V = [0 0; 2 0; 1 1; 0 2; 3 0];
            E = [1 2 1; 2 3 1; 3 4 1; 4 1 1; 2 5 0];
            testCase.verifyError(@() conjQ(QuaPol(V, E, [0 0 0 0 1 0 1 0 0 0], ...
                [1 0; 1 0; 1 0; 1 0; 1 0])), ?MException);
        end


        function aFOURPieceFanConjugatesCorrectly(testCase)
        % Also found by widening the axes: the corner-naming loop used the PRE-merge cell count, so
        % once mergeAdjacentCells actually started firing -- which it only did after Case D was
        % restructured -- it ran off the end of the list and crashed with MATLAB:badsubscript, an
        % unnamed crash rather than a named refusal. Four pieces is the smallest input that
        % triggered a merge.
        %
        % The oracle is closed form: each quadrant carries a DIAGONAL H, so its own sup separates
        % into sum over axes of (t>0) t^2/(2h), and f* is the max of the four.
            V = [0 0; 1 0; 0 1; -1 0; 0 -1];
            E = [1 2 0; 1 3 0; 1 4 0; 1 5 0];
            F = [1 4; 2 1; 3 2; 4 3];
            Hs = {[1 0;0 1], [2 0;0 1], [1 0;0 2], [3 0;0 3]};
            fc = zeros(4,10);
            for k = 1:4
                fc(k,:) = [0 0 0 0, Hs{k}(1,1), Hs{k}(1,2), Hs{k}(2,2), 0, 0, 0];
            end
            g = conjQ(QuaPol(V, E, fc, F));

            rng(20260904);
            S = [randn(300,2)*2; 0 0];
            sgn = {[1 1], [-1 1], [-1 -1], [1 -1]};
            want = -inf(size(S,1),1);
            for k = 1:4
                u = S(:,1)*sgn{k}(1);  v = S(:,2)*sgn{k}(2);
                want = max(want, (u>0).*(u.^2/(2*Hs{k}(1,1))) + (v>0).*(v.^2/(2*Hs{k}(2,2))));
            end
            testCase.verifyEqual(g.eval(S), want, 'RelTol', 1e-12, 'AbsTol', 1e-12);
        end

        function aHALFPLANEPieceConjugatesCorrectly(testCase)
        % UN-QUARANTINED 2026-09-04. This was an assumeFail: QuaPol.examples{11} and {12} -- the
        % plane cut by a line into two HALF-PLANES -- produced cells that OVERLAPPED with different
        % values, and eval's first-match returned the smaller. TWO defects, both specific to a
        % piece bounded BY A LINE, where every vertex and every recession direction lies ON that
        % line so the piece's own geometry says nothing:
        %
        %   1. The edge's SIDE could not be derived, so no half-plane was recorded and the domain
        %      was silently the whole plane. Only F knows; it is now consulted there.
        %   2. The OUTWARD normal was re-derived by orienting against those same vertices and
        %      directions, so one of the two opposite rays got the wrong side and its KKT
        %      multiplier condition was flipped. It is now read from the half-plane already decided.
        %   3. The RECESSION CONE of a half-plane is two-dimensional and no two rays generate it,
        %      so a null direction pointing into it was missed and the sup came out FINITE where it
        %      is +infinity.
        %
        % Verified against the legacy conjugate, an independent implementation: 0 value and 0
        % domain disagreements on examples(12), where before there were 28 and 91.
            P = QuaPol.examples();
            rng(20260904);
            S = [randn(150,2)*2; 0 0; 1 1; -1 -1; -2.670 -0.895; 2.503 0.013];
            for k = [11 12]
                f = P{k};
                testCase.assumeTrue(f.isExact());
                g = conjQ(f);

                % the invariant the defect broke, asserted for BOTH fixtures
                testCase.verifyTrue(checkQuaConConsistent(g, 400), ...
                    sprintf('examples(%d): cells overlap carrying different functions', k));

                % NO LEGACY COMPARISON HERE, deliberately. Calling f.conj() on these fixtures runs
                % the SYMBOLIC engine and took this suite from 10 s to 128 s -- CLAUDE.md section 7
                % budgets the fast bucket in seconds, and section 9 puts symbolic work in the slow
                % bucket. The legacy cross-check lives in .claude/legacy-diff.m, which is where a
                % deliberate, occasional run belongs; what is asserted HERE is the two properties
                % that need no oracle at all.

                % and, oracle or not, f* may never fall below a sup actually attained on the
                % domain. A dozen dual points is enough: this is a guard against a whole region
                % being wrong, which is what both defects were, not a hunt for a lone bad point --
                % and it keeps the suite inside the fast bucket's seconds budget.
                lo = min(f.V,[],1) - 3;  hi = max(f.V,[],1) + 3;
                X = lo + rand(1500,2) .* (hi - lo);
                q = f.eval(X);  X = X(isfinite(q),:);  q = q(isfinite(q));
                Sfew = S([1:6, end-1, end], :);
                vf = g.eval(Sfew);
                for i = 1:size(Sfew,1)
                    testCase.verifyGreaterThanOrEqual(vf(i), max(X*Sfew(i,:).' - q) - 1e-6, ...
                        sprintf('examples(%d) at s=(%g,%g)', k, Sfew(i,1), Sfew(i,2)));
                end
            end
        end

        function noConjugateInTheFixtureCorpusHasOverlappingCells(testCase)
        % THE INVARIANT ITSELF, applied broadly. A mesh whose cells overlap while carrying
        % different functions does not define a function at all -- eval resolves it by first match,
        % so the answer depends on cell ORDER, and the failure is silent. Measured across the
        % corpus: 2 of 29 conjugates were inconsistent before the half-plane fixes, 0 of 29 after.
        %
        % The sample counts are small on purpose. checkQuaConConsistent's EXACT half -- pairs whose
        % constraints are all straight, decided by ratQ.feasible2 -- is what caught both defects,
        % and it does not depend on nSample at all; the sampled half is a backstop for curved pairs.
        % Keeping the counts low keeps this inside the fast bucket's seconds budget.
            names = {'energy', 'oneNorm', 'oneNormConjugate'};
            for i = 1:numel(names)
                g = conjQ(QuaPol.(names{i})());
                testCase.verifyTrue(checkQuaConConsistent(g, 250), names{i});
            end
            P = QuaPol.examples();
            for k = 1:numel(P)
                if ~isa(P{k}, 'QuaPol') || ~P{k}.isExact(), continue, end
                try
                    g = conjQ(P{k});
                catch
                    continue                     % a refusal is not this test's business
                end
                testCase.verifyTrue(checkQuaConConsistent(g, 250), sprintf('examples(%d)', k));
            end
        end


        function aDualDomainThatIsASINGLEPOINTIsRecoveredNotDiscarded(testCase)
        % QuaPol.examples{19} is nine AFFINE pieces on a fan whose dual domain is the single point
        % s = (0,0), where f*(0,0) = -inf f = 0. assembleQuaConCells drops cells with no
        % two-dimensional interior, so the whole answer came back with ZERO faces -- +infinity
        % everywhere -- losing a correct finite value.
        %
        % The point is now recovered exactly: a thin cell's extreme points are intersections of
        % pairs of its own bounding LINES, hence rational and computable with ratQ.solve2, and if
        % every candidate that satisfies its cell is the SAME point, the domain is that point. It
        % is then emitted with EQUALITY sides, the machinery caseAFullDomain already uses.
        %
        % Only the single-point case is recovered, deliberately. Keeping every cell that is merely
        % nonempty AS A SET was tried and measured WORSE: 982 degenerate cells survived and eval's
        % tolerance then admitted points genuinely outside the domain, turning one wrong answer into
        % two. A point is exact, checkable, and cannot over-cover.
            P = QuaPol.examples();
            f = P{19};
            testCase.assumeTrue(f.isExact());
            g = conjQ(f);

            testCase.verifyEqual(g.nf, 1, 'one cell: the point');
            testCase.verifyEqual(g.eval([0 0]), 0, 'AbsTol', 1e-12, ...
                'f*(0,0) = -inf f = 0, which the legacy conjugate also gives');
            testCase.verifyTrue(all(isinf(g.eval([1 1; 0.01 0; 0 0.01; -1 2; 5 -3]))), ...
                'and +infinity everywhere else -- the domain really is that one point');
            for k = 1:g.nf
                testCase.verifyTrue(all(g.FC{k}(:,2) == 0), ...
                    'the cell is carried by EQUALITY sides');
            end
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

        % ---- Case D: every quadratic that is not positive definite, on a bounded polygon --------

        function anIndefiniteQuadraticMatchesItsSeparableClosedForm(testCase)
        % x^2/2 - y^2/2 over [0,1]^2 SEPARATES, and each axis has an elementary 1-D conjugate:
        %       sup_{0<=x<=1} s1 x - x^2/2  = s1*xc - xc^2/2 with xc = clamp(s1, 0, 1)   (concave)
        %       sup_{0<=y<=1} s2 y + y^2/2  = max(0, s2 + 1/2)                            (convex,
        %                                     so the max is at an endpoint)
        % Written down independently of anything conjQ computes.
            V = [0 0; 1 0; 1 1; 0 1];  E = [1 2 1; 2 3 1; 3 4 1; 4 1 1];  F = [1 0;1 0;1 0;1 0];
            g = conjQ(QuaPol(V, E, [0 0 0 0 1 0 -1 0 0 0], F));

            rng(20260904);
            S = [randn(400,2)*2; 0 0; 1 1; -1 -1; 0.5 -0.5];
            xc = min(1, max(0, S(:,1)));
            want = S(:,1).*xc - xc.^2/2 + max(0, S(:,2) + 0.5);
            [got, idx] = g.eval(S);
            testCase.verifyTrue(all(idx > 0), 'the cells must cover the plane');
            testCase.verifyEqual(got, want, 'RelTol', 1e-12, 'AbsTol', 1e-12);
        end

        function everyHessianClassOnAPolygonAgreesWithAnIndependentMaximisation(testCase)
        % The dichotomy conjQ dispatches on is "H positive definite or not", so the sweep has to
        % cover BOTH sides of it and the boundary between them: positive definite, indefinite,
        % negative definite, PSD-singular, NSD-singular and identically zero.
            rng(20260904);
            classes = {'PD', 'indefinite', 'ND', 'PSD-singular', 'NSD-singular', 'affine'};
            for c = 1:numel(classes)
                for rep = 1:2
                    [V, E, F] = conjQTest.randomConvexPolygon();
                    H = conjQTest.hessianOfClass(classes{c});
                    L = randi([-4 4], 2, 1);  k0 = randi([-4 4]);
                    f = QuaPol(V, E, [0 0 0 0, H(1,1), H(1,2), H(2,2), L(1), L(2), k0], F);
                    testCase.assumeTrue(f.isExact());
                    g = conjQ(f);
                    S = [randn(40,2)*3; 0 0];
                    [got, idx] = g.eval(S);
                    testCase.verifyTrue(all(idx > 0), ...
                        sprintf('%s rep %d: a dual point fell in no cell', classes{c}, rep));
                    want = zeros(size(S,1),1);
                    for i = 1:size(S,1)
                        want(i) = conjQTest.supAnyQ(S(i,:).', V, H, L, k0);
                    end
                    testCase.verifyEqual(got, want, 'RelTol', 1e-9, 'AbsTol', 1e-9, ...
                        sprintf('%s rep %d', classes{c}, rep));
                end
            end
        end

        function theConcaveCaseIsUNCHANGEDByTheGeneralBoundaryBranch(testCase)
        % Case D subsumes the concave case rather than replacing it: for a negative semidefinite H
        % no edge has positive curvature, so no edge object is built, the parts list has one entry,
        % no fold happens, and the result IS the vertex max. This pins that the subsumption is
        % exact and not merely close -- the object must still be all-affine with one cell per
        % winning vertex.
            V = [0 0; 1 0; 1 1; 0 1];  E = [1 2 1; 2 3 1; 3 4 1; 4 1 1];  F = [1 0;1 0;1 0;1 0];
            g = conjQ(QuaPol(V, E, [0 0 0 0, -2, -1, -2, 1, -2, 3], F));
            testCase.verifyEqual(g.nf, 4, 'one cell per corner, and no edge cells at all');
            testCase.verifyTrue(all(g.fN(:,5) == 0 & g.fN(:,6) == 0 & g.fN(:,7) == 0), ...
                'a concave conjugate is piecewise AFFINE -- no curvature can appear');
            for j = 1:g.ne
                testCase.verifyEqual(g.edgeKind(j), 'line');
            end
        end

        function anIndefiniteConjugateIsConvexAndCarriesACurvedEdge(testCase)
        % f* is convex whatever f is -- broken immediately by a wrong split side. And the fold of
        % pieces with different Hessians must produce a genuine conic boundary, which is the whole
        % reason QuaCon exists rather than QuaPar.
            V = [0 0; 2 0; 2 1; 0 2];  m = size(V,1);
            E = [(1:m).', [2:m,1].', ones(m,1)];  F = [ones(m,1), zeros(m,1)];
            g = conjQ(QuaPol(V, E, [0 0 0 0, 2, 3, -1, 1, 1, 0], F));

            rng(20260904);
            A = randn(200,2)*2;  B = randn(200,2)*2;
            testCase.verifyLessThanOrEqual(g.eval((A+B)/2), (g.eval(A) + g.eval(B))/2 + 1e-9);

            kinds = arrayfun(@(j) string(g.edgeKind(j)), 1:g.ne);
            testCase.verifyTrue(any(kinds ~= "line"), ...
                'an indefinite piece produces edge cells whose differences are genuine conics');
        end

        function theReportedFaceCountIsAnUpperBoundAndTheValuesAreStillRight(testCase)
        % HONEST PIN OF A KNOWN GAP, and it asserts the gap's BOUNDARY rather than its size. A cell
        % that is empty only because of a CURVED constraint is not detected -- assembleQuaConCells
        % says so -- so nf over-counts. What must NOT happen is a wrong value or an uncovered
        % point, and that is what is asserted here. Measured 2026-09-04 on this fixture: 274 faces
        % reported, 75 ever occupied over 200k samples, carrying 10 distinct functions. The size is
        % deliberately not pinned: it will move as the filters improve, and a number pinned here
        % would turn an improvement into a failure.
            u = [2;1];  H = u*u.';                         % PSD-singular: the worst case for this
            V = [0 0; 3 0; 3 2; 1 3; 0 2];  m = size(V,1);
            E = [(1:m).', [2:m,1].', ones(m,1)];  F = [ones(m,1), zeros(m,1)];
            g = conjQ(QuaPol(V, E, [0 0 0 0, H(1,1), H(1,2), H(2,2), 1, -1, 2], F));

            rng(20260904);
            S = randn(500,2)*6;
            [got, idx] = g.eval(S);
            testCase.verifyTrue(all(idx > 0), 'every dual point still lands in a face');
            want = zeros(size(S,1),1);
            for i = 1:size(S,1)
                want(i) = conjQTest.supAnyQ(S(i,:).', V, H, [1;-1], 2);
            end
            testCase.verifyEqual(got, want, 'RelTol', 1e-9, 'AbsTol', 1e-9, ...
                'the over-counting must not disturb a single value');

            occupied = numel(unique(idx(idx > 0)));
            testCase.verifyLessThanOrEqual(occupied, g.nf, ...
                'occupancy cannot exceed the reported face count');
        end
        % ---- Case E: a concave or affine quadratic on an UNBOUNDED piece ------------------------

        function theOneNormConjugatesToTheLibrarysOwnExpectedAnswer(testCase)
        % THE BEST ORACLE IN THE REPOSITORY FOR THIS CASE: QuaPol.oneNormConjugate is an
        % independently written fixture, in the library since before any of this work, stating what
        % |x|_1 conjugates to -- the indicator of the unit infinity-ball. Nothing about it came from
        % conjQ, so agreement is a genuine cross-check rather than a restatement.
        %
        % This is also the first case where dom f* is a PROPER SUBSET, which is the mathematics and
        % not a gap: an unbounded domain makes the sup infinite outside a cone.
            g = conjQ(QuaPol.oneNorm());
            want = QuaPol.oneNormConjugate();

            rng(20260904);
            S = [randn(400,2)*1.5; 0 0; 1 1; -1 -1; 0.999 0.999; 1.001 0; 2 0; -3 1];
            testCase.verifyEqual(g.eval(S), want.eval(S), 'AbsTol', 1e-12, ...
                'the exact conjugate must agree with the library''s own expected answer');

            inside  = all(abs(S) <= 1 - 1e-9, 2);
            outside = any(abs(S) >= 1 + 1e-9, 2);
            testCase.verifyTrue(all(g.eval(S(inside,:))  == 0),   'zero on the infinity-ball');
            testCase.verifyTrue(all(isinf(g.eval(S(outside,:)))), '+inf off it');
            testCase.verifyEqual(g.nf, 1, 'the answer is ONE cell -- an indicator of a square');
            testCase.verifyEqual(g.ne, 4, 'bounded by four lines');
        end

        function anIndicatorOfAConeConjugatesToTheIndicatorOfItsPolar(testCase)
        % q = 0 on the first quadrant is the indicator of that cone, and the conjugate of an
        % indicator is the support function -- which for a cone is the indicator of its POLAR,
        % here the third quadrant. Elementary, and it is where a sign slip would surface first.
            f = conjQTest.wedge([0 0], [1 0], [0 1], zeros(1,10));
            g = conjQ(f);
            testCase.verifyEqual(g.nf, 1);
            testCase.verifyEqual(sortrows(g.EcQ), sortrows([0 0 0 1 0 0; 0 0 0 0 1 0]), ...
                'the two bounding lines are the axes');

            S = [-1 -1; -2 -3; 0 0; -0.5 -0.5];
            testCase.verifyEqual(g.eval(S), zeros(4,1), 'zero on the polar cone');
            T = [1 -1; -1 1; 1 1; 3 0];
            testCase.verifyTrue(all(isinf(g.eval(T))), '+inf off it');
        end

        function anAffinePieceOnAWedgeIsTheSupportFunctionShifted(testCase)
        % q(x) = <a,x> + b on a cone C gives f*(s) = sup_C <s-a,x> - b, which is -b on the polar of
        % C shifted by a, and +inf elsewhere. Written down independently below.
            a = [-1; -2];  b = 3;
            f = conjQTest.wedge([0 0], [1 0], [0 1], [0 0 0 0, 0, 0, 0, a(1), a(2), b]);
            g = conjQ(f);
            rng(20260904);
            S = [randn(300,2)*2; a.'; a.' + [-1 -1]; a.' + [1 0]];
            want = inf(size(S,1),1);
            in = (S(:,1) - a(1) <= 1e-12) & (S(:,2) - a(2) <= 1e-12);
            want(in) = -b;
            testCase.verifyEqual(g.eval(S), want, 'AbsTol', 1e-12);
        end



        function aRaysDirectionMarkerIsNotTreatedAsAVertex(testCase)
        % A ray edge [v1 v2 0] stores the ray [V(v1), V(v2)) -- v1 is the base and V(v2) only fixes
        % the DIRECTION. Counting V(v2) as a vertex cannot change a VALUE, since it is a point of
        % the domain and so bounded by the sup, but it manufactures cells that are empty. The wedge
        % has exactly one extreme point, so its conjugate has exactly one cell.
            g = conjQ(conjQTest.wedge([0 0], [1 0], [0 1], zeros(1,10)));
            testCase.verifyEqual(g.nf, 1, ...
                'one extreme point, so one cell -- not one per direction marker as well');
        end
    end

    methods (Static)

        function f = wedge(apex, d1, d2, fc)
        % A two-ray cone: the apex, plus one direction marker along each ray. The F columns are
        % OPPOSITE because the two rays bound the wedge from opposite sides -- F(j,:) is
        % [left, right] of edge j with 0 meaning +inf, so putting the face on the left of both
        % describes no convex wedge at all, and QuaPol.eval then cannot locate its own interior
        % (measured while building this: q(1,1) came back +Inf inside the first quadrant).
            V = [apex; apex + d1; apex + d2];
            f = QuaPol(V, [1 2 0; 1 3 0], fc, [1 0; 0 1]);
        end

        function H = hessianOfClass(cls)
        % A random integer symmetric H of a named definiteness class. The classes ARE the dispatch:
        % conjQ branches on "positive definite or not", so a sweep has to cover both sides and the
        % singular boundary between them.
            switch cls
                case 'PD'
                    while true, H = randi([-4 4],2,2); H = H+H.';
                        if H(1,1) > 0 && det(H) > 0.5, return, end, end
                case 'indefinite'
                    while true, H = randi([-4 4],2,2); H = H+H.';
                        if det(H) < -0.5, return, end, end
                case 'ND'
                    while true, H = randi([-4 4],2,2); H = H+H.';
                        if H(1,1) < 0 && det(H) > 0.5, return, end, end
                case 'PSD-singular'
                    u = randi([-3 3],2,1);  if all(u == 0), u = [1;0]; end
                    H = u*u.';
                case 'NSD-singular'
                    u = randi([-3 3],2,1);  if all(u == 0), u = [1;0]; end
                    H = -u*u.';
                case 'affine'
                    H = zeros(2);
                otherwise
                    error('conjQTest:badClass', 'unknown Hessian class %s', cls);
            end
        end

        function v = supAnyQ(s, V, H, L, k0)
        % THE GENERAL ORACLE: sup over the polygon of <s,x> - q(x), for ANY H. Enumerates the
        % candidate maximisers -- every vertex, every edge's clamped 1-D stationary point whatever
        % that edge's curvature, and the unconstrained stationary point when it lies inside. For a
        % concave objective one of those IS the maximiser; for an indefinite one the maximiser is
        % on the boundary, hence a vertex or an edge stationary point, so the list is complete
        % either way. Every candidate is a genuine point of P, so a spurious one can only ever be
        % dominated -- which is why adding the interior point unconditionally is safe.
            q   = @(x) 0.5 * x.' * H * x + L.' * x + k0;
            obj = @(x) s.' * x - q(x);
            cand = V.';
            if abs(det(H)) > 1e-12
                xs = H \ (s - L);
                if conjQTest.inPolygon(xs.', V), cand = [cand, xs]; end
            end
            m = size(V,1);
            for j = 1:m
                a = V(j,:).';  b = V(mod(j,m)+1,:).';  d = b - a;
                al = d.' * H * d;
                if abs(al) < 1e-14, continue, end     % affine along this edge: endpoints suffice
                t = min(1, max(0, (s - (H*a + L)).' * d / al));
                cand = [cand, a + t*d]; %#ok<AGROW>
            end
            v = -inf;
            for i = 1:size(cand,2), v = max(v, obj(cand(:,i))); end
        end

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
