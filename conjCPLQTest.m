classdef conjCPLQTest < matlab.unittest.TestCase
    % Tests for the 'cplq' conjugate engine (conjCPLQ) and the QuaPol operator interface.
    % Covers the currently implemented case (full-domain quadratics) plus the
    % validation/rejection behaviour. See DESIGN.md II.5.1.

    methods (Test)
        function energyIsSelfConjugate(testCase)
            % f = 1/2(x^2+y^2) is self-conjugate.
            p = QuaPol.energy();
            q = p.conj('cplq');
            testCase.verifyEqual(q.f, p.f, 'AbsTol', 1e-12);
        end

        function defaultEngineIsCPLQ(testCase)
            % conj() with no engine argument defaults to 'cplq'.
            p = QuaPol.energy();
            testCase.verifyEqual(p.conj().f, p.conj('cplq').f, 'AbsTol', 1e-12);
        end

        function biconjOfConvexQuadraticIsItself(testCase)
            % f** = f for a closed convex function.
            p = QuaPol.energy();
            b = p.biconj('cplq');
            testCase.verifyEqual(b.f, p.f, 'AbsTol', 1e-12);
        end

        function biconjCoverageByInputCase(testCase)
            % Pins which inputs biconj handles, and by which route, case by case.
            %
            % HISTORY (closed 2026-07-28): Case B used to ERROR here. biconj was literally
            % conj-of-conj, and the conjugate of a bounded-domain function is finite everywhere --
            % an UNBOUNDED multi-face domain, which conjCPLQ still rejects (Case C requires
            % isDomBounded), so the second conjugation had nowhere to go. biconj now routes
            % through biconjCPLQ, which returns Step 1's convex envelope for a bounded triangle:
            % conv(q + I_T) IS f** when T is compact, so no second conjugation is needed. The
            % underlying conjCPLQ gap (unbounded multi-face conj) is untouched and still open --
            % SUPPORT_MATRIX.md section 1.2. Correctness of the new route is checked against a
            % pipeline-free ground truth in biconjCPLQTest.
            E3 = [1 2 1; 2 3 1; 3 1 1]; F3 = [1 0; 1 0; 1 0];

            % Case A -- full-domain strictly convex quadratic: conjugate is again a full-domain
            % quadratic, so the second conjugation is Case A too, and biconj closes the loop.
            caseA = QuaPol([1 0 1 0 0 0]);
            testCase.verifyEqual(caseA.conj().kind(), 'QuaPol');
            testCase.verifyEqual(caseA.biconj().kind(), 'QuaPol');

            % Case B -- single bounded triangle: conj gives a genuine mesh QuaPar; biconj is the
            % convex envelope, a RatPol on the same triangle. f = xy vanishes at all three
            % vertices and is >= 0 on this triangle, so f** is identically zero on it.
            caseB = QuaPol([0 0; 1 0; 0 1], E3, [0 1 0 0 0 0], F3);
            testCase.verifyEqual(caseB.conj().kind(), 'QuaPar');
            bB = caseB.biconj();
            testCase.verifyEqual(bB.kind(), 'RatPol');
            testCase.verifyEqual(bB.eval([0.2 0.2; 1/3 1/3]), [0; 0], 'AbsTol', 1e-12);

            % ...and a convex triangle piece is its own biconjugate -- now WITHOUT taking Case B's
            % route at all. `biconj` IS the closed-convex-envelope operator, so a convex f is
            % returned unchanged by the short-circuit at the top of biconjCPLQ (commit 0027900,
            % 436 s -> 0.05 s), which means the answer keeps the input's own type: QuaPol, not the
            % RatPol Step 1 would have widened it to. This expectation said RatPol until
            % 2026-08-19 and was stale from that commit onward -- the slow bucket had not been run
            % since, which is exactly what the gap in the handoff was about. The VALUE check below
            % is what this line is really for and is unchanged.
            caseBconvex = QuaPol([0 0; 1 0; 0 1], E3, [1 0 1 0 0 0], F3);
            bBc = caseBconvex.biconj();
            S = [0.2 0.2; 0.6 0.2; 1/3 1/3];
            testCase.verifyEqual(bBc.kind(), 'QuaPol');
            testCase.verifyEqual(bBc.eval(S), caseBconvex.eval(S), 'AbsTol', 1e-12);

            % Case C -- general bounded multi-face domain. Its CONVEX member no longer takes the
            % symbolic route at all: conjCPLQ's Case B2 fan-triangulates a bounded polygon and
            % conjugates each triangle in closed form, so this comes back as a mesh QuaPar. The
            % values are unchanged by that (measured across five quadratics on this exact domain,
            % worst |new - old| = 0); what changed is the route, and hence the type.
            V = [0 0; 1 0; 1 1; 0 1];
            E = [1 2 1; 2 3 1; 1 3 1; 3 4 1; 4 1 1];
            F = [1 0; 1 0; 2 1; 2 0; 2 0];
            caseC = QuaPol(V, E, [1 0 1 0 0 0; 1 0 1 0 0 0], F);
            % ASSERTION CHANGED 2026-08-24, and the reason is what this line is FOR. It read
            % `verifyEqual(caseC.conj().kind(), 'QuaPar')` -- the EXACT class. Its own comment
            % above says what it is about: "what changed is the route, and hence the type", i.e.
            % closed-form numerics rather than the symbolic fallback. A convex face is now
            % conjugated WHOLE by conjConvexPolygon instead of being fan-triangulated, and the
            % result is a QuaPol: a QuaPar (the lattice's own subsumption) that additionally
            % guarantees no edge is curved. It satisfies the property this test exists to pin
            % strictly better, so pinning the exact class was pinning a representation detail the
            % test itself says it is not about.
            gC = caseC.conj();
            testCase.verifyInstanceOf(gC, 'QuaPar', 'must be a mesh, not the symbolic form');
            testCase.verifyFalse(isa(gC, 'QuaParCPLQ'), 'must not fall back to the symbolic path');

            % An INDEFINITE quadratic on the same domain used to go to Case C, because every
            % triangle's conjugate is then curved and Step 3 took at most one curved operand.
            % Arc-vs-arc assembly (2026-08-13) lifted that limit, so this now completes on the
            % numeric route and returns a MESHED QuaPar -- pinned here the other way round, with
            % the values checked, since the route is the thing this test is about.
            caseCind = QuaPol(V, E, [0 1 0 0 0 0; 0 1 0 0 0 0], F);
            gInd = caseCind.conj();
            testCase.verifyEqual(gInd.kind(), 'QuaPar');
            for sPt = {[3 -1], [-2 3], [1 1], [4 4], [-3 -3], [0.5 0.5]}
                sv = sPt{1};
                testCase.verifyEqual(gInd.eval(sv), max([0, sv(1), sv(2), sv(1)+sv(2)-1]), ...
                    'AbsTol', 1e-9, sprintf('indefinite Case C at (%g,%g)', sv(1), sv(2)));
            end

            % CASE C's BICONJUGATE NOW WORKS, and this assertion used to pin it as broken.
            %
            % HISTORY. It first read `verifyEqual(caseC.biconj().kind(), 'QuaParCPLQ')`, which
            % passes on an EMPTY piece list -- QuaParCPLQ(functionNDomain.empty()).kind() is still
            % 'QuaParCPLQ' -- so it hid the defect entirely. Measured on pristine HEAD
            % (2026-07-31): caseC.conj() gave 9 pieces and caseC.biconj() gave ZERO, i.e.
            % f** = +inf everywhere, for an f that is CONVEX and hence its own biconjugate. It was
            % then changed to verifyError, pinning the loud failure instead of the silent one,
            % with the failing chain traced into functionNDomain.getInterior.
            %
            % None of that chain was fixed. What changed is that the FIRST conjugation no longer
            % goes through it: Case B2 answers this domain in closed form, so the second
            % conjugation is handed a clean mesh QuaPar instead of the symbolic pieces whose
            % latent getInterior bug it used to trip over. f = (x^2+y^2)/2 is convex, so f** = f
            % exactly, and that is what is asserted now -- the strongest available statement, and
            % one no return-type check can fake.
            %
            % The getInterior chain is NOT closed: it is still on the path of any INDEFINITE
            % domain, which is why caseCind above stays on the symbolic route. See
            % SUPPORT_MATRIX.md section 7.
            bC = caseC.biconj();
            SC = [0.2 0.2; 0.5 0.5; 0.8 0.3; 0.1 0.9];
            for i = 1:size(SC,1)
                testCase.verifyEqual(conjCPLQTest.evalConjResult(bC, SC(i,:)), ...
                    caseC.eval(SC(i,:)), 'AbsTol', 1e-9, sprintf( ...
                    'f** must equal f for a convex f, at (%g,%g)', SC(i,1), SC(i,2)));
            end
        end

        function generalPositiveDefiniteQuadratic(testCase)
            % f(x) = 1/2 x'Q x + L'x + k,  Q=[2 0;0 4], L=[1;-1], k=3.
            % f*(s) = 1/2 (s-L)' inv(Q) (s-L) - k.
            f6 = [2 0 4 1 -1 3];          % [x^2 xy y^2 x y const]
            p  = QuaPol(f6);
            q  = p.conj('cplq');
            Q = [2 0; 0 4]; L = [1; -1]; k = 3;
            M = inv(Q); grad = -M*L; d = 0.5*(L'*M*L) - k;
            expf6 = [M(1,1) M(1,2) M(2,2) grad(1) grad(2) d];
            expq  = QuaPol(expf6);
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
            p  = QuaPol(f6);
            q  = p.conj('cplq');
            s  = [3; 2];
            % direct: sup_x <s,x> - f(x), maximizer x* = Q\(s-L)
            Q = [2 0; 0 4]; L = [1; -1]; k = 3;
            xstar = Q \ (s - L);
            fval  = 0.5*xstar'*Q*xstar + L'*xstar + k;
            direct = s'*xstar - fval;
            testCase.verifyEqual(q.eval(s'), direct, 'AbsTol', 1e-10);
        end

        function fullDomainNonStrictlyConvexIsANSWEREDInAllTHREECases(testCase)
        % A full-domain quadratic that is not strictly convex gives three different objects, and as
        % of 2026-08-27 `conj` RETURNS all three. They were one lumped `notImplemented` until
        % 2026-08-25, then three classified refusals, and each fell in turn -- every time because
        % the representation already existed and one routine elsewhere did not know about it.
        %
        %   Q indefinite / a negative eigenvalue  -> conv f = -inf, so f* = +inf EVERYWHERE.
        %                                            dom f* is empty; the FUNCTION is a full-domain
        %                                            mesh whose constant is +inf.
        %   Q = 0 (affine)                        -> f* is the indicator of the single point L:
        %                                            a NEEDLE (nv=1, ne=0).
        %   Q PSD of rank 1                       -> f* is finite only on a LINE through L:
        %                                            two opposite RAYS from one point.
        %
        % This test covers the first; the other two have their own tests below, which assert their
        % values against the closed form.
            for c = {{[0 1 0 0 0 0], 'f = xy, indefinite'}, ...
                     {[-2 0 -2 0 0 0], 'concave'}, ...
                     {[1 0 -1 0 0 0], 'x^2 - y^2'}}
                f6 = c{1}{1}; nm = c{1}{2};
                g = QuaPol(f6).conj('cplq');
                testCase.verifyEqual(g.kind(), 'QuaPar', nm);
                for s = {[0 0], [3 -2], [-1 5], [1e4 1e4]}
                    testCase.verifyTrue(isinf(g.eval(s{1})) && g.eval(s{1}) > 0, sprintf( ...
                        '%s: f* must be +inf at (%g,%g)', nm, s{1}(1), s{1}(2)));
                end
            end
        end

        function theEMPTYDomainConjugateRoundTripsToMinusInfinity(testCase)
        % REPLACES theFullDomainRefusalsCarryTheClosedForm, which checked that a refusal message
        % carried its closed form. There is no refusal left in this family to check -- all three
        % cases are answered as of 2026-08-27 -- so what is worth pinning instead is the one
        % consequence a caller has to know about.
        %
        % f* for the empty-domain case is +inf everywhere. Conjugating THAT gives -inf everywhere,
        % which is mathematically correct ((+inf)* = -inf) and is NOT a PLQ function: it is outside
        % the class this library represents. So a caller computing f** on a function with a
        % negative curvature direction gets a -inf mesh rather than an error.
        %
        % Pinned so the behaviour is deliberate rather than discovered. If a future return type
        % gains a way to say "not a PLQ function", this is the site that should use it.
            g = QuaPol([0 1 0 0 0 0]).conj('cplq');
            testCase.verifyTrue(isinf(g.eval([0 0])) && g.eval([0 0]) > 0, ...
                'f* of an indefinite full-domain quadratic is +inf everywhere');
            h = g.conj('cplq');
            for s = {[0 0], [2 -3]}
                v = h.eval(s{1});
                testCase.verifyTrue(isinf(v) && v < 0, sprintf( ...
                    'f** must be -inf at (%g,%g) -- correct, and outside the PLQ class', ...
                    s{1}(1), s{1}(2)));
            end
        end

        function theAffineFullDomainConjugateIsAPointAndIsRETURNED(testCase)
        % B3's point case, closed 2026-08-25 (overnight). f = 3x - 2y + 5 is affine on all of R^2,
        % so f*(s) = sup_x <s - L, x> - kappa is +inf unless s = L = (3,-2), and -kappa = -5 there.
        % dom f* is a single POINT.
        %
        % It used to raise `PLQ:conjCPLQ:conjugateIsAPoint`, and the reason was never the
        % mathematics -- it was that `QuaPar.eval` returned +inf on a NEEDLE domain (nv=1, ne=0),
        % including at the needle's own vertex, so the answer had nowhere to live. eval now has
        % that branch, so the answer is simply returned.
            g = QuaPol([0 0 0 3 -2 5]).conj('cplq');
            testCase.verifyEqual(g.kind(), 'QuaPar');
            testCase.verifyEqual(g.nv, 1, 'dom f* is a single point');
            testCase.verifyEqual(g.ne, 0, 'a needle has no edges');
            testCase.verifyEqual(g.eval([3 -2]), -5, 'AbsTol', 1e-12, ...
                'f* at the point it is supported on');
            for s = {[3 -1.9], [0 0], [-3 2], [3.001 -2]}
                v = g.eval(s{1});
                testCase.verifyTrue(isinf(v), sprintf( ...
                    'f*(%g,%g) must be +inf off the point, got %g', s{1}(1), s{1}(2), v));
            end
        end

        function theRank1PSDConjugateIsALineAndIsRETURNED(testCase)
        % B3's LINE case, closed 2026-08-26. For Q positive semidefinite of rank 1, f is affine
        % along null(Q), so f* is finite only on the LINE {s : <s - L, n> = 0} with n spanning
        % null(Q), and equals 1/2 (s-L)' pinv(Q) (s-L) - kappa there.
        %
        % It used to be refused because "QuaPar's dim<2 domain is a SEGMENT or ray between two
        % vertices, not a line". A line IS two opposite RAYS from one point, which the constructor
        % has always accepted; what was missing was in `eval` (its edge-chain branch treated every
        % edge as a segment, ignoring the ray marker) and in `belongToEdge` (its ray range test was
        % coordinate-wise, so it only worked for rays pointing into the positive quadrant). Both
        % are fixed and pinned in QuaParTest.
            for c = {{[1 0 0 0 0 0], 'x^2/2'}, {[1 1 1 2 0 0], '(x+y)^2/2 + 2x'}}
                f6 = c{1}{1}; nm = c{1}{2};
                g = QuaPol(f6).conj('cplq');
                testCase.verifyEqual(g.kind(), 'QuaPar', nm);
                testCase.verifyEqual(g.nv, 3, [nm ': a line is an apex plus two ray directions']);
                testCase.verifyEqual(g.ne, 2, [nm ': two opposite rays']);
                Q = [f6(1) f6(2); f6(2) f6(3)]; L = [f6(4); f6(5)]; kappa = f6(6);
                M = pinv(Q);
                n = null(Q); n = n(:,1).' / norm(n(:,1));
                t = [-n(2), n(1)];
                for a = [0 1 -2.5 4]
                    s = L(:).' + a*t;
                    want = 0.5*(s(:)-L).' * M * (s(:)-L) - kappa;
                    testCase.verifyEqual(g.eval(s), want, 'AbsTol', 1e-9, sprintf( ...
                        '%s: f* on the line at (%g,%g)', nm, s(1), s(2)));
                end
                for a = [0.7 -1.3]
                    s = L(:).' + a*n;
                    testCase.verifyTrue(isinf(g.eval(s)), sprintf( ...
                        '%s: f* must be +inf off the line at (%g,%g)', nm, s(1), s(2)));
                end
            end
        end

        function affineTriangleViaOrchestrator(testCase)
            % conjCPLQ dispatches a single bounded-triangle piece straight to conjPieceCPLQ (no
            % Step 1 needed for an affine piece). ell = -x over (0,0),(1,0),(0,1).
            V = [0 0; 1 0; 0 1]; E = [1 2 1; 2 3 1; 3 1 1]; F = [1 0; 1 0; 1 0];
            p = QuaPol(V, E, [0 0 0 -1 0 0], F);
            g = p.conj('cplq');
            testCase.verifyClass(g, 'QuaPar');
            testCase.verifyEqual(g.nf, 3);
            S = [-2 -1; 3 -1; -1 2; 0.5 0.5];
            for i = 1:size(S,1)
                expected = max([0; S(i,1)+1; S(i,2)]);
                testCase.verifyEqual(g.eval(S(i,:)), expected, 'AbsTol', 1e-10);
            end
        end

        function convexQuadraticTriangleViaOrchestrator(testCase)
            % PD quadratic over a triangle: conjCPLQ should match conjPieceCPLQ directly (no
            % envelope needed, since a PD piece is already its own convex envelope).
            A = [2 1; 1 3]; b = [1; -2]; cc = 0.5;
            V = [0 0; 2 0; 1 2]; E = [1 2 1; 2 3 1; 3 1 1]; F = [1 0; 1 0; 1 0];
            f6 = [A(1,1) A(1,2) A(2,2) b(1) b(2) cc];
            p = QuaPol(V, E, f6, F);
            g = p.conj('cplq');
            testCase.verifyEqual(g.f, conjPieceCPLQ(p).f, 'AbsTol', 1e-12);
            qf = @(x) 0.5*x'*A*x + b'*x + cc;
            xin = [0.8; 0.5]; s = A*xin + b;
            testCase.verifyEqual(g.eval(s'), s'*xin - qf(xin), 'AbsTol', 1e-8);
        end

        function concaveTriangleViaOrchestratorSidestepsToEnvelope(testCase)
            % A concave piece cannot be conjugated directly (conjPieceCPLQ rejects it); conjCPLQ
            % must fall back to Step 1's affine envelope automatically. q = -(x^2+y^2) over
            % (0,0),(2,0),(0,2): f*(s) = max_i(<s,v_i> - q(v_i)).
            V = [0 0; 2 0; 0 2]; E = [1 2 1; 2 3 1; 3 1 1]; F = [1 0; 1 0; 1 0];
            q = QuaPol(V, E, [-2 0 -2 0 0 0], F);
            g = q.conj('cplq');
            testCase.verifyClass(g, 'QuaPar');
            testCase.verifyEqual(g.nf, 3);
            qf = @(x) -(x(1)^2 + x(2)^2);
            v1=[0;0]; v2=[2;0]; v3=[0;2];
            S = [1 1; 5 -1; -1 5; 3 3];
            for i = 1:size(S,1)
                s = S(i,:)';
                expected = max([s'*v1-qf(v1); s'*v2-qf(v2); s'*v3-qf(v3)]);
                testCase.verifyEqual(g.eval(S(i,:)), expected, 'AbsTol', 1e-10);
            end
        end

        function indefiniteTriangleZeroOrOneConvexEdgeViaOrchestrator(testCase)
            % Genuinely indefinite pieces with 0 or 1 convex edge are conjugated directly (no
            % envelope needed) -- exercise both through the orchestrator.
            V0 = [0 0; 1 0; 0 1]; E0 = [1 2 1; 2 3 1; 3 1 1]; F0 = [1 0; 1 0; 1 0];
            g0 = QuaPol(V0, E0, [0 1 0 0 0 0], F0).conj('cplq');
            testCase.verifyEqual(g0.nf, 3);   % zero convex edges -> 3-cone piecewise-linear

            V1 = [0 0; 2 0; 1 1]; E1 = [1 2 1; 2 3 1; 3 1 1]; F1 = [1 0; 1 0; 1 0];
            g1 = QuaPol(V1, E1, [0 1 0 0 0 0], F1).conj('cplq');
            testCase.verifyEqual(g1.nf, 6);   % one convex edge -> 6-face parabolic QuaPar
        end

        function indefiniteTriangleTwoConvexEdgesSidestepsToEnvelope(testCase)
            % Two convex edges: conjPieceCPLQ rejects the raw piece directly, so conjCPLQ must
            % fall back to Step 1's rank-1-PSD envelope (COAP Appendix A.4) automatically, matching
            % what conjPieceCPLQTest/psdRank1QuadraticEndToEnd does by hand.
            V = [0 0; 2 1; 1 2]; E = [1 2 1; 2 3 1; 3 1 1]; F = [1 0; 1 0; 1 0];
            q = QuaPol(V, E, [0 1 0 0 0 0], F);
            g = q.conj('cplq');
            env = convEnvCPLQ(q);
            testCase.verifyEqual(env.nf, 1);   % single-face envelope: no Step 3 needed here
            testCase.verifyEqual(g.f, conjPieceCPLQ(env).f, 'AbsTol', 1e-12);
        end

        function indefiniteTriangleThreeConvexEdgesUsesStep3(testCase)
            % Three convex edges: Step 1 splits the triangle into sub-triangle pieces, so the
            % orchestrator needs Step 3 (max of conjugates, maxQuaPar) to combine their Step-2
            % conjugates. A strictly increasing vertex chain (both x and y increasing) makes all 3
            % pairwise slopes positive, i.e. all 3 edges convex for f=xy -- the same T=conv{(0,0),
            % (3,3),(1,2)} example maxQuaParTest.m validates by hand.
            %
            % UPDATE (Part 2c, 2026-07-17/18 session): Step 1 used to always split a 3-convex-edge
            % triangle into exactly 2 plain-quadratic sub-triangles; it now correctly recurses
            % (convEnvCPLQ.m's solveTriangleBF/splitTwoConvexEdges) into each of those 2
            % sub-triangles, since each is itself an ordinary 2-convex-edge triangle subject to the
            % SAME tightness criterion (tangentCevian) already proven for the standalone
            % nCE==2 case -- there is nothing special about a 3CE sub-triangle that would exempt
            % it. For THIS triangle both sub-triangles need the further split, so Step 1 now
            % produces 4 pieces (2 plain quadratic + 2 genuinely rational), not 2 -- confirmed
            % correct (not a regression) per DESIGN.md's Part 2c. `conjPieceCPLQ` cannot yet
            % conjugate a genuinely rational piece (a pre-existing, separately-tracked gap -- see
            % conjPieceCPLQ.m's own TODO and DESIGN.md), so `conj('cplq')` on this triangle must now
            % correctly error rather than silently combining only the 2 formerly-assumed pieces.
            % This test now pins BOTH the corrected split count and the new, loud, correct failure
            % mode (same "pin the loud failure instead of the old silent success" pattern already
            % used elsewhere in this file's history -- see
            % maxQuaParTest.matchHalfEdgesRejectsSameSideRayPairingAndDropsSubsumedPieces).
            % UPDATE (2026-07-28): Step 2 now falls back to cPLQ's own symbolic Step 2/3 for a
            % rational envelope face (conjCPLQ's conjEnvelopeViaCPLQ). That fixes the
            % 2-convex-edge split outright -- see
            % indefiniteTriangleTwoConvexEdgesSplitViaCPLQStep2 below -- and for THIS triangle it
            % gets Step 2 to complete as well. The remaining failure is in Step 3, cPLQ's
            % cross-piece maximum (plq.maximum:185 -> maximumConjugate -> functionNDomain.maximumP
            % -> region.maximum), not in the conjugate of any piece.
            %
            % UPDATE (2026-07-29): Step 3 now RUNS to completion on this triangle -- four
            % separate blockers were fixed (region.maxArray's degenerate probe directions,
            % functionNDomain.mergeL's stale removeTangent vertices, region.linear3pt's
            % self-overwriting index, and plq_1p.conjugateFunction's hard-coded nCE==2 grad
            % half-planes). Step 1 and Step 2 are now EXACT here: the pointwise max of the four
            % per-piece conjugates matches sup_{x in T} <s,x> - xy at all 289 points of a 17x17
            % dual grid.
            %
            % What is still wrong is only the ASSEMBLY. region.merge unions two regions by
            % deleting their shared facet and intersecting what remains -- an over-approximation
            % unless the two carry identical other constraints -- and simplifyUnboundedRegion
            % deletes any constraint not passing through a finite vertex. Both leave regions
            % covering territory that was never theirs, with the wrong value on it (~12% of that
            % grid). conjCPLQ's assertStep3MatchesPieces cross-checks the assembled maximum
            % against the per-piece max and raises cplqFailed on disagreement, so this stays a
            % LOUD failure rather than a silently wrong answer -- same convention as before, but
            % now triggered by the real defect instead of by a crash.
            %
            % Nothing here needs a "3 convex edges" case: [COAP] Appendix A.5's split reduces such
            % a triangle to 2-convex-edge sub-triangles, and Step 1 already does it -- that is why
            % the envelope below has 4 faces. The edge count is how the input is DESCRIBED; it is
            % not what fails.
            %
            % An earlier version of this comment blamed the numeric/symbolic boundary and called
            % for exact arithmetic in Step 1. That was wrong and is worth recording: the envelope
            % coefficients are accurate to ~1e-16, and at the vertex where the denominator
            % vanishes the numerator vanishes too with grad N parallel to grad den (residual 0 to
            % 1.3e-16), so the 0/0 is cleanly REMOVABLE. What was missing was resolving it by a
            % limit -- see symbolicFunction.limitDirectional.
            % UPDATE (2026-07-29, later session): the assembly is FIXED and this case now closes
            % end to end. region.simplifyUnboundedRegion deletes a constraint only when
            % region.redundantSubset certifies it redundant by LP, and region.merge deletes a
            % shared facet only when region.unionIsExact certifies A u B = A' n B' (equivalently,
            % that the union is convex). Measured on the 17x17 dual grid this test's own gate
            % samples: the assembled partition went from 57 of 289 points wrong to 0 -- every one
            % of the fold's seven rounds is now exact, and the whole fold got FASTER (1645 s vs
            % 1768 s), because a correct partition needs fewer regions than a damaged one.
            % Unit-level coverage of the two certificates is in regionTest.m.
            %
            % This test therefore no longer pins a loud failure; it pins the answer. The gate
            % (conjCPLQ's assertStep3MatchesPieces) stays -- it is a real invariant, and cheap
            % next to the pipeline it follows.
            V = [0 0; 3 3; 1 2]; E = [1 2 1; 2 3 1; 3 1 1]; F = [1 0; 1 0; 1 0];
            q = QuaPol(V, E, [0 1 0 0 0 0], F);
            testCase.verifyEqual(convEnvCPLQ(q).nf, 4);   % confirms the (now recursive) split
            g = q.conj('cplq');
            % ROUTE CHANGED 2026-08-13: the four-face envelope's conjugate now assembles
            % numerically (arc-vs-arc), so this returns a meshed QuaPar instead of the symbolic
            % QuaParCPLQ. The ground-truth checks below are unchanged and are what matters.
            testCase.verifyEqual(g.kind(), 'QuaPar');

            % Ground truth is the closed-form sup over T, not anything the pipeline produced.
            S = [1 1; -3 -3; -7 5.25; 0 0; 3 3; 2 -1; -1 2; 5 5];
            for i = 1:size(S,1)
                s = S(i,:);
                expected = convEnvCPLQTest.supBilinearOverPoly(s, V);
                got = g.eval(s);
                testCase.verifyTrue(isfinite(got), ...
                    sprintf('no region covers s=(%.6f,%.6f)', s(1), s(2)));
                testCase.verifyEqual(got, expected, 'AbsTol', 1e-9, ...
                    sprintf('s=(%.6f,%.6f)', s(1), s(2)));
            end
        end

        function indefiniteTriangleTwoConvexEdgesSplitIsFullyNumeric(testCase)
            % The 2-convex-edge tightness split (convEnvCPLQTest.bilinearTwoConvexEdgesSplitIsTight)
            % leaves Step 1 with a RATIONAL face, and Step 2 has no rational branch. This case
            % therefore used to fall back to cPLQ's SYMBOLIC Step 2/3 and return a QuaParCPLQ.
            %
            % It no longer does, and no rational-conjugate formula was needed to close it. Step 1
            % puts on each sub-triangle T_k the envelope of the ORIGINAL q over THAT sub-triangle,
            % so (R_k + I_{T_k})* = (conv(q + I_{T_k}))* = (q + I_{T_k})* exactly -- the rational
            % face is a re-expression of a piece whose conjugate Step 2 already has in closed form
            % (an indefinite quadratic over a 1-convex-edge triangle). conjCPLQ's
            % conjFaceOrOriginal conjugates the original over that face's domain instead; see its
            % header for the derivation and the measurement.
            %
            % What this test pins is therefore the RETURN TYPE as much as the values: a QuaPar
            % means the whole triangle went through closed-form numerics with no call into the
            % Symbolic Math Toolbox. If a future change sends it back to the symbolic fallback the
            % values would still be right and only this line would notice.
            V = [2 1; 0 0; 1 0]; E = [1 2 1; 2 3 1; 3 1 1]; F = [1 0; 1 0; 1 0];
            q = QuaPol(V, E, [0 1 0 0 0 0], F);           % f = xy
            testCase.verifyEqual(convEnvCPLQ(q).nf, 2);   % the tightness split
            g = q.conj('cplq');
            testCase.verifyClass(g, 'QuaPar');

            sBad = [-0.008727 -0.999962];                 % the paper's own flagged dual point
            S = [sBad; 1.90 2.50; -1 -1; 0.5 0.5; 3 -2; -3 2; 2 2];
            for i = 1:size(S,1)
                s = S(i,:);
                expected = convEnvCPLQTest.supBilinearOverPoly(s, V);
                got = g.eval(s);
                testCase.verifyTrue(isfinite(got), ...
                    sprintf('no region covers s=(%.6f,%.6f)', s(1), s(2)));
                testCase.verifyEqual(got, expected, 'AbsTol', 1e-9, ...
                    sprintf('s=(%.6f,%.6f)', s(1), s(2)));
            end
        end

        function caseCValuesAreCorrectForAGeneralQuadratic(testCase)
        % THE TEST THAT WAS MISSING. Case C was only ever asserted by RESULT TYPE, never by a
        % value, and it was returning the conjugate of co(x*y) whatever the face carried --
        % measured on pristine HEAD, q=(x^2+y^2)/2 over the unit square gave f*(0.3,0.4) = 0.4
        % where the truth is 0.125, f*(-0.5,0.2) = 0.2 for 0.02, f*(2,-1) = 2 for 1.5.
        %
        % Step 1 now classifies by the SIGNS OF THE EIGENVALUES of Q rather than by nCE (which
        % tests edge slopes and so only classifies x*y): convex and affine keep q as their own
        % envelope, concave get the affine interpolant through the actual values of q, and
        % indefinite are moved into the frame where q IS x*y (xyFrame.m) so that cPLQ's own
        % closed forms apply to the function they were written for.
        %
        % Evaluated through evalConjResult rather than g.fnd, because these five inputs no longer
        % all take the same ROUTE: since Case B2 (conjCPLQ's bounded-polygon branch) the convex
        % and concave ones are answered by closed-form numerics and come back as a QuaPar, while
        % the three indefinite ones still fall back to the symbolic pipeline and come back as a
        % QuaParCPLQ, which is the only kind that has an .fnd. What this test is FOR is the
        % values, and those are unchanged: measured across all five, worst |new - old| = 0.
            V = [0 0; 1 0; 1 1; 0 1];
            E = [1 2 1; 2 3 1; 1 3 1; 3 4 1; 4 1 1];
            F = [1 0; 1 0; 2 1; 2 0; 2 0];
            S = [0.3 0.4; 1 1; -0.5 0.2; 2 -1; 0 0; -2 -2; 1.5 0.5];
            % f6 = [x^2 xy y^2 x y const]; matrixForm reads Q = [c5 c6; c6 c7].
            fs = { [1 0 1 0 0 0], ...      % (x^2+y^2)/2   convex
                   [0 1 0 0 0 0], ...      % x*y           indefinite, cPLQ's own canonical case
                   [2 0 -2 0 0 0], ...     % x^2 - y^2     indefinite, needs a frame change
                   [0 3 0 7 -2 5], ...     % 3xy+7x-2y+5   indefinite with an affine part
                   [-1 0 -1 0 0 0] };      % -(x^2+y^2)/2  concave
            [gx, gy] = meshgrid(linspace(0,1,801));
            G = [gx(:), gy(:)];
            for k = 1:numel(fs)
                f6 = fs{k};
                g = QuaPol(V, E, [f6; f6], F).conj('cplq');
                Q = [f6(1) f6(2); f6(2) f6(3)]; L = [f6(4); f6(5)]; c = f6(6);
                qg = 0.5*sum((G*Q).*G,2) + G*L + c;
                for t = 1:size(S,1)
                    ref = max(G*S(t,:)' - qg);           % sup over the square, to grid resolution
                    got = conjCPLQTest.evalConjResult(g, S(t,:));
                    testCase.verifyFalse(isnan(got), sprintf( ...
                        'f = %s: s=(%g,%g) is covered by no dual region', ...
                        mat2str(f6), S(t,1), S(t,2)));
                    testCase.verifyEqual(got, ref, 'AbsTol', 1e-6, sprintf( ...
                        'f = %s at s=(%g,%g)', mat2str(f6), S(t,1), S(t,2)));
                end
            end
        end

        function multiFaceUnboundedConvexFacesConjugateExactly(testCase)
        % An UNBOUNDED multi-face domain whose faces carry CONVEX quadratics, end to end through
        % the public entry. Each face's co q = q, so Step 1 does nothing and Step 2 conjugates a
        % curved function via conjConvexOverPiece's active-set cells; Step 3 then assembles them.
        % This used to be refused outright.
        %
        % Truth is available in closed form because each face is SEPARABLE on its own quadrant:
        % sup over a half-line of s*t - a*t^2/2 is s^2/(2a) when s points into the half-line and
        % 0 otherwise, and f* is the max over the four faces.
            V = [0 0;-1 0; 0 1;1 0;0 -1];
            E = [1 2 0;1 3 0;1 4 0;1 5 0];
            f = [1 0 1 0 0 0;1 0 2 0 0 0;2 0 2 0 0 0;2 0 1 0 0 0];
            F = [1 2;2 3;3 4;4 1];
            p = QuaPol(V,E,f,F);
            testCase.verifyFalse(p.isDomBounded);
            g = p.conj('cplq');
            % Which quadrant each face occupies is DERIVED from the face's own half-planes,
            % not hardcoded: F fixes the face-to-cone assignment and getting it wrong by hand
            % silently compares against the wrong closed form.
            pp = quaPolToPlq(p);
            ax = zeros(4,2); sg = zeros(4,2);
            for k = 1:4
                ax(k,:) = [f(k,1), f(k,3)];                  % Q = [c5 c6; c6 c7], c6 = 0 here
                [A, ~, lin] = pp.pieces(k).d.polygon.linearForm;
                A = A(lin,:);
                for t = [1 -1]
                    for u = [1 -1]
                        if all(A*[t;u]*0.5 <= 1e-9), sg(k,:) = [t u]; end
                    end
                end
            end
            half = @(s,a,dir) (sign(s)==dir) * s^2/(2*a);
            S = [0 0; 1 1; -1 -1; 2 0.5; -0.5 2; 3 -1; -2 3; 0.25 -0.75];
            for t = 1:size(S,1)
                best = -inf;
                for k = 1:4
                    best = max(best, half(S(t,1),ax(k,1),sg(k,1)) + half(S(t,2),ax(k,2),sg(k,2)));
                end
                % READ CHANGED 2026-08-24: was evalFunctionNDomain(g.fnd, ...), which only a
                % QuaParCPLQ has. This shape now takes the NUMERIC route -- conjCPLQ's polygon path
                % is no longer gated on isDomBounded, and conjConvexPolygon handles a convex face
                % with recession directions -- so g is a meshed QuaPar here. evalConjResult is this
                % file's own route-agnostic reader and its header says it exists for exactly this.
                got = conjCPLQTest.evalConjResult(g, S(t,:));
                testCase.verifyEqual(got, best, 'AbsTol', 1e-9, ...
                    sprintf('at s=(%g,%g)', S(t,1), S(t,2)));
            end
        end

        function step3UnboundedAssemblyAgreesWithItsOwnPieces(testCase)
        % RENAMED 2026-08-15 from step3DropsCellsOnSomeUnboundedAssemblies, whose assertion was
        % `verifyError(..., 'PLQ:conjCPLQ:cplqFailed')` -- it pinned the GATE firing on this
        % input. The defect the gate was catching is FIXED, so the gate no longer fires and the
        % old assertion is now backwards. What this pins instead is the property the gate exists
        % to protect: the assembled cross-piece maximum agrees with the pointwise max of the
        % per-piece conjugates -- the same f*, computed the other way, which is exactly what
        % conjCPLQ's own assertStep3MatchesPieces checks before returning.
        %
        % THE TWO DEFECTS THIS INPUT FOUND, both now closed, both worth keeping:
        %
        % (1) A DROP, fixed 2026-08-02. The assembled maximum kept only 4 of the 16 cells, losing
        % face 1's s_2^2/2 cell on {s1<=0, s2>=0}, so f*(-0.5,2) came back 1.125 for a truth of 2.
        % region.simplifyUnboundedRegion declared a genuine 2-D cone EMPTY, because it decides
        % that from probe directions built out of constraint SLOPES at a vertex and the split
        % conic's gradient VANISHES at exactly that vertex. Fixed by region.witnessAwayFrom.
        %
        % (2) An OVER-CLAIM, fixed 2026-08-15, and it is the SAME trap one level over. At
        % s = (-3,-2.4) the assembly gave 5.130 where the per-piece max gives 4.500 (right: the
        % four cone suprema are 0, 4.5, 3.69, 2.88). The cell carrying s1^2/4 + s2^2/2 had lost
        % its -s1 <= 0 and kept only {s2 <= 0, s2^2/2 - s1^2 <= 0, s1^2 - 2*s2^2 <= 0} -- two
        % constraints BLIND TO THE SIGN of s1, so the region was symmetric under s1 -> -s1 and
        % claimed the mirror wedge. region.removeTangent deleted that constraint after building a
        % "tangent" to a quadratic AT ITS OWN APEX, where the gradient vanishes and no tangent
        % exists. Fixed by refusing to conclude anything from a vanishing gradient.
        %
        % Measure both halves here, since one input exhibited both.
            V = [0 0;-1 0; 0 1;1 0;0 -1];
            E = [1 2 0;1 3 0;1 4 0;1 5 0];
            f = [1 0 1 0 0 0;1 0 2 0 0 0;2 0 2 0 0 0;2 0 1 0 0 0];
            F = [3 2;2 1;1 4;4 3];
            p = QuaPol(V,E,f,F);
            g = p.conj('cplq');     % assertStep3MatchesPieces runs inside, and must not fire
            testCase.verifyTrue(isa(g, 'RatPar'));
            for sPt = {[-0.5 2], [-3 -2.4], [3 3], [-2 -3], [1 1]}
                s = sPt{1};
                testCase.verifyEqual(conjCPLQTest.evalConjResult(g, s), ...
                    conjCPLQTest.fourConeTruth(s), 'AbsTol', 1e-9, ...
                    sprintf('4-cone fan at (%g,%g)', s(1), s(2)));
            end
        end

        function step3UnboundedAssemblyMatchesTheTruth(testCase)
        % The SAME 4-cone fan as step3UnboundedAssemblyAgreesWithItsOwnPieces, pinned against a
        % CLOSED-FORM truth rather than against the pipeline's own per-piece max. The two are
        % complementary and were resolved together on 2026-08-15: that one says the assembly
        % agrees with its own pieces, this one says both are right.
        %
        % Ground truth needs no pipeline. Each face carries q = a*x^2 + b*y^2 on one quadrant, and
        % for a cone C = {sigma1*x >= 0, sigma2*y >= 0} the sup separates:
        %       (q + I_C)*(s) = [s1^2/(4a) if sigma1*s1 >= 0 else 0]
        %                     + [s2^2/(4b) if sigma2*s2 >= 0 else 0],
        % and f* is the max of the four. The face-to-quadrant assignment below is read off F's
        % left/right convention: edge 3 runs along +x with face 1 on its left, so face 1 is the
        % first quadrant, and so on around the fan.
        %
        % s = (-3,-2.4) is the point the gate currently catches: the four cone suprema there are
        % 0, 4.5, 3.69 and 2.88, so f* = 4.5, while the assembly returns 5.130 = s1^2/4 + s2^2/2 --
        % face 4's cell, which belongs on {s1>=0, s2<=0}. Some region has grown across s1 = 0.
            V = [0 0;-1 0; 0 1;1 0;0 -1];
            E = [1 2 0;1 3 0;1 4 0;1 5 0];
            f = [1 0 1 0 0 0;1 0 2 0 0 0;2 0 2 0 0 0;2 0 1 0 0 0];
            F = [3 2;2 1;1 4;4 3];
            p = QuaPol(V,E,f,F);
            g = p.conj('cplq');
            S = [-3 -2.4; -0.5 2; 1 1; -1 -1; 2 -0.5; 0 0; 3 3; -2 -3];
            for t = 1:size(S,1)
                s = S(t,:);
                got = conjCPLQTest.evalConjResult(g, s);   % route-agnostic: see the note above
                testCase.verifyEqual(got, conjCPLQTest.fourConeTruth(s), 'AbsTol', 1e-9, ...
                    sprintf('at s=(%g,%g)', s(1), s(2)));
            end
        end

        function multiFaceBoundedDomainViaCPLQIntegration(testCase)
            % Case C (Phase 1 cPLQ integration, this session): a genuinely multi-triangle BOUNDED
            % nonconvex PLQ now conjugates end to end through the wired-in conjCPLQ dispatch --
            % exactly the case Case B's own numeric path (conjPieceCPLQ+maxQuaPar) cannot do
            % (maxQuaPar refuses curved-edge QuaPar inputs from independent triangles). Reuses
            % cplqAdapterTest's f=xy-over-a-diagonally-split-square example, but now going through
            % the PUBLIC conj('cplq') entry point rather than calling quaPolToPlq directly.
            V = [0 0; 1 0; 1 1; 0 1];
            E = [1 2 1; 2 3 1; 3 1 1; 3 4 1; 4 1 1];
            F = [1 0; 1 0; 1 2; 2 0; 2 0];
            f = [0 1 0 0 0 0; 0 1 0 0 0 0];   % xy on both faces
            q = QuaPol(V, E, f, F);
            testCase.verifyEqual(q.nf, 2);
            testCase.verifyTrue(q.isDomBounded);

            g = q.conj('cplq');
            % ROUTING CHANGED 2026-08-13: the numeric route now completes here, so the result is a MESHED QuaPar rather than the symbolic QuaParCPLQ this used to pin. That is the arc-vs-arc work of 2026-08-13: maxQuaPar assembles two curved operands, so conjCPLQ no longer has to fall back to cPLQ's Step 2/3. The VALUES are what this test is really for, and they are checked below (exact, error 0, against the closed-form sup).
            testCase.verifyTrue(isa(g, 'RatPar'));

            nt = 220; [uu,vv] = meshgrid(linspace(0,1,nt));
            Xg = uu(:); Yg = vv(:); xyg = Xg.*Yg;
            % s=(0.5,0.5): exact tie point between the two triangles' vertex cones -- see
            % cplqAdapterTest.m / functionNDomain.maximumP's HISTORY comment / SESSION_HANDOFF.md.
            S = [3 -1; -2 3; 1 1; 0 -3; 4 4; -3 -3; 6 2; -1 6; 2 2; 0.5 0.5];
            for i = 1:size(S,1)
                sup = max(S(i,1)*Xg + S(i,2)*Yg - xyg);
                gv = g.eval(S(i,:));
                testCase.verifyEqual(gv, sup, 'AbsTol', 2e-3, sprintf('s=%d', i));
            end
        end

        function cubicRejectedByOperators(testCase)
            % Cubic numerator is storable but rejected by operators (allowed for isConvex only).
            p = QuaPol([1 0 0 0 0 0 0 0 0 0]);   % x^3/6 term present -> degree 3
            testCase.verifyEqual(p.degree, 3);
            testCase.verifyError(@() p.conj('cplq'), 'PLQ:op:unsupportedType');
        end

        function unimplementedEnginesError(testCase)
            p = QuaPol.energy();
            testCase.verifyError(@() p.conj('pqp'),   'PLQ:conj:engine');
            testCase.verifyError(@() p.conj('graph'), 'PLQ:conj:engine');
        end

        function partialConjEngineRestriction(testCase)
            p = QuaPol.energy();
            testCase.verifyError(@() p.partialConj(1,'graph'), 'PLQ:partialConj:engine');
        end

        function plqvcAliasStillWorks(testCase)
            % Backward compatibility: PLQVC is an alias of QuaPol.
            p = PLQVC.energy();                 % inherited static factory
            testCase.verifyTrue(isa(p, 'QuaPol'));
            V = [0 0;-1 0; 0 1;1 0;0 -1];
            E = [1 2 0;1 3 0;1 4 0;1 5 0];
            f = [1 0 1 0 0 0;1 0 2 0 0 0;2 0 2 0 0 0;2 0 1 0 0 0];
            F = [1 2;2 3;3 4;4 1];
            p2 = PLQVC(V,E,f,F);                % construct via the alias
            testCase.verifyTrue(isa(p2, 'PLQVC'));
            testCase.verifyEqual(p2.nf, 4);
        end

        function boundedPolygonsTakeTheClosedFormPath(testCase)
        % Case B2. A bounded polygon is a union of triangles and Case B conjugates a triangle in
        % closed form, so f* = max_k (q + I_T_k)* needs no symbolic engine. Before this branch
        % existed EVERY polygon went to Case C's symbolic pipeline, which is the single largest
        % source of runtime in this suite.
        %
        % The RETURN TYPE is the assertion that matters: a QuaPar means closed-form numerics
        % throughout, a QuaParCPLQ means it fell back to the symbolic path. Values are checked
        % against an EXACT reference -- the max of <s,x> - q(x) over the vertices, the edge
        % stationary points and the interior stationary point of each triangle -- never a sampled
        % one. A sampled reference reported ~3e-5 here and that was the reference's own error.
            polys = { [0 0; 1 0; 1 1; 0 1]                    % unit square
                      [0 0; 2 0.3; 2.4 1.7; 0.6 2.1]          % general quadrilateral
                      [0 0; 2 1; 3 3; 1 2]                    % parallelogram
                      [0 0; 2 0; 2.6 1.5; 1 2.5; -0.4 1.3] }; % pentagon
            quads = { [1 0 1 0 0 0]        % convex   (x^2+y^2)/2
                      [0 0 0 3 -2 1]       % affine   3x-2y+1
                      [-1 0 -1 0 0 0] };   % concave -(x^2+y^2)/2
            S = [0.3 0.4; -1 -1; 1 1; 2 -0.5; 0 0; 3 -2; -3 2; 1.5 0.2];
            for pi = 1:numel(polys)
                W = polys{pi};
                q = conjCPLQTest.polygonQuaPol(W);
                for qi = 1:numel(quads)
                    f6 = quads{qi};
                    q.f = repmat([0 0 0 0, f6], q.nf, 1);
                    g = q.conj('cplq');
                    % ASSERTION CHANGED 2026-08-24: was verifyClass(g,'QuaPar'), which demands the
                    % EXACT class. This test's header states the property it wants -- "a QuaPar
                    % means closed-form numerics throughout, a QuaParCPLQ means it fell back to the
                    % symbolic path" -- and since conjConvexPolygon a CONVEX polygon comes back as
                    % a QuaPol, which is a QuaPar with every edge conic pinned to zero. That is the
                    % same property, more strongly. The VALUE checks below are untouched.
                    testCase.verifyInstanceOf(g, 'QuaPar', sprintf( ...
                        'polygon %d quadratic %d is not a mesh', pi, qi));
                    testCase.verifyFalse(isa(g, 'QuaParCPLQ'), sprintf( ...
                        'polygon %d quadratic %d fell back to the symbolic path', pi, qi));
                    for t = 1:size(S,1)
                        testCase.verifyEqual(g.eval(S(t,:)), ...
                            conjCPLQTest.supQuadOverPoly(S(t,:), f6, W), 'AbsTol', 1e-9, ...
                            sprintf('polygon %d quadratic %d at s=(%g,%g)', ...
                                    pi, qi, S(t,1), S(t,2)));
                    end
                end
            end
        end

        function indefiniteOverAPolygonStillFallsBack(testCase)
        % The BOUNDARY of what Case B2 buys, pinned so it is not mistaken for a bug. An
        % indefinite triangle conjugates to a PARABOLIC QuaPar, and maxQuaPar takes at most one
        % curved operand, so a polygon carrying an indefinite quadratic cannot be assembled
        % numerically and must fall through to Case C. Closing this is the arc-vs-arc face
        % clipping gap in maxQuaPar, not anything in conjCPLQ.
            V = [0 0; 1 0; 1 1; 0 1];
            E = [1 2 1; 2 3 1; 1 3 1; 3 4 1; 4 1 1];
            F = [1 0; 1 0; 2 1; 2 0; 2 0];
            q = QuaPol(V, E, [0 1 0 0 0 0; 0 1 0 0 0 0], F);   % f = xy on both triangles
            g = q.conj('cplq');
            % NO LONGER FALLS BACK. maxQuaPar took at most one curved operand when this was
            % written, which is why an indefinite quadratic over a polygon had to go to Case C.
            % Arc-vs-arc assembly (2026-08-13) removed that limit, so the numeric route completes
            % and returns a meshed QuaPar. Kept as a pin on the ROUTE, now the other way round,
            % with the values checked against the closed-form sup over the unit box.
            testCase.verifyEqual(g.kind(), 'QuaPar');
            S = [3 -1; -2 3; 1 1; 0 -3; 4 4; -3 -3; 2 2; 0.5 0.5];
            for i = 1:size(S,1)
                s = S(i,:);
                truth = max([0, s(1), s(2), s(1)+s(2)-1]);
                testCase.verifyEqual(g.eval(s), truth, 'AbsTol', 1e-9, sprintf('s=%d', i));
            end
        end

        function unboundedNonConvexFaceWithAFiniteEnvelopeIsANSWERED(testCase)
        % `TODO.md` G3 / B4 says an unbounded face carrying a non-convex quadratic is a GAP, and
        % `SUPPORT_MATRIX.md` had it as one. MEASURED 2026-08-25: that is the NUMERIC route's
        % decline, and it is not the whole story -- the symbolic route answers this family, and
        % answers it EXACTLY. Pinned here so the classification cannot drift back.
        %
        % f = x*y on the first quadrant. On the recession cone {x,y >= 0} f is bounded below by 0,
        % so conv f is finite and f* is a genuine function:
        %       f*(s) = sup_{x,y >= 0} s1 x + s2 y - x y
        % Taking y = 0 gives s1 x, which runs to +inf when s1 > 0, and symmetrically for s2; for
        % s <= 0 the interior critical point (x,y) = (s2,s1) is infeasible and the sup is attained
        % at the origin. So f* is the INDICATOR of the third quadrant: 0 there, +inf elsewhere.
        %
        % Cost ~14 s, symbolic -- this suite's bucket.
            q = QuaPol([0 0; 1 0; 0 1], [1 2 0; 1 3 0], [0 1 0 0 0 0], [1 0; 0 1]);
            g = q.conj('cplq');
            inside  = [-1 -1; -3 -2; -0.5 -0.25; 0 0; -5 0; 0 -5];
            outside = [1 1; 2 -1; -1 2; 0.5 0.5];
            for i = 1:size(inside,1)
                s = inside(i,:);
                testCase.verifyEqual(conjCPLQTest.evalConjResult(g, s), 0, 'AbsTol', 1e-9, ...
                    sprintf('f*(%g,%g) must be 0 on the third quadrant', s(1), s(2)));
            end
            for i = 1:size(outside,1)
                s = outside(i,:);
                v = conjCPLQTest.evalConjResult(g, s);
                testCase.verifyTrue(isnan(v) || isinf(v), sprintf( ...
                    'f*(%g,%g) must be +inf off the third quadrant, got %g', s(1), s(2), v));
            end
        end

        function unboundedNonConvexFaceUNBOUNDEDBelowIsRefusedByName(testCase)
        % The other half of B4, and it is not the same gap. f = x^2 - y^2 on the first quadrant
        % runs to -inf along the y-axis, which is IN the recession cone, so conv f = -inf and
        % f* = +inf everywhere: dom f* is EMPTY.
        %
        % That is the right answer and it is refused rather than returned, for the same reason as
        % `PLQ:conjCPLQ:conjugateHasEmptyDomain` in Case A -- no mesh encodes an empty domain
        % (nf = 0 means dim < 2, not empty). What this test pins is that the refusal is BY NAME and
        % identifies the cause, rather than a generic notImplemented that would read as "we have
        % not got to it yet".
            q = QuaPol([0 0; 1 0; 0 1], [1 2 0; 1 3 0], [1 0 -1 0 0 0], [1 0; 0 1]);
            testCase.verifyError(@() q.conj('cplq'), 'convEnvUnbounded:minusInfinity');
        end

        function aGENUINETwoPieceDomainFoldsCorrectlyAndCHEAPLY(testCase)
        % THE CROSS-PIECE MAXIMUM, on input that actually reaches it. Every two-piece fixture in
        % this repository carries the SAME quadratic on both faces, so Step 0's mergeSameQuadFaces
        % deletes the shared edge and hands the rest of the pipeline ONE face -- the fold is never
        % performed. Measured 2026-08-25: that is exactly why `testcPLQ`'s hole at (-0.5,2) does
        % not reach `conj`.
        %
        % These two carry DIFFERENT functions per face, so they survive Step 0 (asserted below --
        % if a future Step 0 merged them this test would silently stop testing anything) and the
        % fold runs for real. The truth is the max of the two per-piece sups, each computed in
        % closed form: both faces are triangles carrying a convex or affine function, so the sup is
        % attained at a vertex, an edge stationary point, or the interior stationary point.
        %
        % CHEAP ON PURPOSE. The obvious fixture -- PRect's two polygons with different quadratics
        % -- did not finish in 40 minutes. The unit square split by its diagonal does the same job
        % in 1.3 s, which is what lets this live in a fast bucket rather than being a thing nobody
        % runs.
            V = [0 0; 1 0; 1 1; 0 1];
            E = [1 2 1; 2 3 1; 1 3 1; 3 4 1; 4 1 1];
            F = [1 0; 1 0; 2 1; 2 0; 2 0];
            T1 = V([1 2 3],:); T2 = V([1 3 4],:);
            cases = { {[0 0 0 1 0 0; 0 0 0 0 1 0], 'affine x | affine y'}, ...
                      {[2 0 2 0 0 0; 2 0 2 1 0 0], 'convex | convex shifted'} };
            S = [0 0; 1 0; 0 1; 1 1; -1 1; 1 -1; -1 -1; 2 0.5; -0.5 2; 3 -2; 0.4 0.6];
            for c = cases
                f6 = c{1}{1}; nm = c{1}{2};
                q = QuaPol(V, E, f6, F);
                testCase.verifyEqual(q.nf, 2, ...
                    [nm ': Step 0 must NOT merge these, or the fold is never exercised']);
                g = q.conj('cplq');
                for i = 1:size(S,1)
                    sPt = S(i,:);
                    want = max(conjCPLQTest.supQuadOverPoly(sPt, f6(1,:), T1), ...
                               conjCPLQTest.supQuadOverPoly(sPt, f6(2,:), T2));
                    got = conjCPLQTest.evalConjResult(g, sPt);
                    testCase.verifyEqual(got, want, 'AbsTol', 1e-9, sprintf( ...
                        '%s: f*(%g,%g)', nm, sPt(1), sPt(2)));
                end
            end
        end

        function frameChangedPieceKeepsItsEnvelopeBLOCKS(testCase)
        % G5, the `MATLAB:badsubscript` crash. Found by checkConjAgainstDefinition's random sweep
        % (seed 20260824) as case 29, a 4-gon carrying x*y plus an affine part; case 17, a 5-gon,
        % crashes the same way. NOT 5-gon-specific, which is what `TODO.md` recorded.
        %
        % THE DEFECT, measured. `plq_1p.conjugate` has a FRAME-CHANGE branch: an indefinite
        % quadratic that is not literally x*y is redone in the z-frame where it is, and the answer
        % read back. That branch copied the z-frame object's ENVELOPE (2 faces here) but replaced
        % its `conjfia` -- the per-face block boundaries into `conjugates` -- with the single block
        % [1 nConj+1]. `maximumConjugate` then loops over `size(envelope,2)` blocks and indexes
        % `conjfia(i+1)`, so the second face asks for `conjfia(3)` of a 2-element array.
        %
        % Observed before the fix, both triangles of case 29:
        %       envelope faces = 2, numel(conjfia) = 2, nConjugates = 11
        %
        % WHY THIS ASSERTS THE INVARIANT AND NOT THE VALUE. The end-to-end value test was written
        % first and cannot be run: with the index repaired, `maximumConjugate` goes on to take the
        % cross-face symbolic max this input always needed, and `maximumP` on this fixture's two
        % envelope faces did not finish in 25 minutes -- measured both in the z-frame (rational
        % coefficients) and after the read-back (surds), with no difference. So G5's crash was
        % standing in front of a symbolic max that does not terminate in practical time, and a
        % value assertion here would be a test nobody can run rather than a check on the defect.
        %
        % What IS the defect is the broken invariant, and it is exact, cheap and one stage early:
        % after `conjugate`, `conjfia` must declare one block per envelope face plus a terminator.
        % Before the fix this read 2 against 2 faces; it is the mismatch, not the crash, that the
        % out-of-bounds index is a symptom of. Stage cost: ~40 s, all of it `conjugate`.
            W = [ 0.21821484235235822 -0.82638242938157225
                  0.38208834489342713 -1.0607216447629879
                  1.2127969933300677  -0.6725798784826531
                  0.20468128742742464 -0.69679526375087664];
            f6 = [0 1 0 -0.73790021007093154 0.0021521450703039899 1.0552416486837408];
            n = size(W,1);
            E = [(1:n)', mod((1:n),n)'+1, ones(n,1)];
            F = [ones(n,1), zeros(n,1)];
            q = QuaPol(W, E, f6, F);

            p = quaPolToPlq(q);
            p = p.triangulate;
            sawFrame = false;
            for i = 1:p.nPieces
                pc = p.pieces(i).convexEnvelope;
                pc = pc.conjugate;
                nFaces = size(pc.envelope, 2);
                if isempty(pc.frame), continue, end
                sawFrame = true;
                testCase.verifyEqual(numel(pc.conjfia), nFaces + 1, sprintf( ...
                    ['piece %d took the frame-change branch and its envelope has %d faces, so ' ...
                     'conjfia must declare %d block boundaries; it declares %d. maximumConjugate ' ...
                     'indexes conjfia(i+1) for i up to the face count, which is the G5 ' ...
                     'MATLAB:badsubscript.'], i, nFaces, nFaces + 1, numel(pc.conjfia)));
                testCase.verifyEqual(pc.conjfia(end), size(pc.conjugates,2) + 1, sprintf( ...
                    'piece %d: the last block boundary must be one past the conjugate cells', i));
            end
            testCase.verifyTrue(sawFrame, ...
                'no piece took the frame-change branch: the fixture no longer exercises G5');
        end

        function sweepCase21HitsTheKnownStep3LegacyGapNotANewFailure(testCase)
        % REVISED 2026-08-28. This test used to pin case 21's PRE-diff behaviour (a fast, named
        % foldDroppedACell refusal, 2.4 s) and fail if the parked `assemblePieces` diff (G1/G10)
        % ever landed, since applying it changed case 21's outcome to a 292.5 s crash. That
        % premise is gone: the diff has now LANDED, deliberately, because this session traced
        % the "crash" to `SUPPORT_MATRIX.md` section 1.2's already-tracked, pre-existing Step 3
        % legacy unreliability (`conjEnvelopeViaCPLQ` -> `plq.maximumConjugate` ->
        % `functionNDomain.maximumP` -> `region.maximum`), not a new defect the diff introduced
        % -- confirmed by reading the error's own message text, which already names that section
        % (DECISIONS.md 2026-08-28, "the parked assembly diff's MuPAD crash is the KNOWN Step 3
        % gap, not new"). The trade was made deliberately: the diff fixes G1/G4/G10 on the
        % scaling fixture that actually matters for SCIP (a numeric-path defect), at the cost of
        % this ONE symbolic Case-C fallback input reaching an already-open, already-scoped
        % problem instead of refusing early. What this test protects now is narrower and still
        % real: that case 21 does not develop a THIRD, different failure mode nobody has seen.
        %
        % Fixture: `checkConjAgainstDefinition`'s own sweep case 21 (seed 20260824, the G4/G10
        % triangle: 3-gon, indefinite xy), reconstructed by hand since `randomCase` is a private
        % local function -- the coefficients below are that seed's case 21, verified byte-for-
        % byte against a fresh sweep run before this test was written.
            W = [ 0.605704715097184  0.930075181116088
                 -0.335394747201423  0.525152429255813
                 -1.08249961702189   0.0844860974378195];
            f6 = [0 1 0 -0.717791341344706 -0.607564534720844 -0.678183523294723];
            E = [1 2 1; 2 3 1; 3 1 1];
            F = [1 0; 1 0; 1 0];
            q = QuaPol(W, E, f6, F);

            testCase.verifyError(@() q.conj('cplq'), 'PLQ:conjCPLQ:cplqFailed');
        end

        function threePieceEllipticalEdgeWitnessNowCompletes(testCase)
        % doc/QuaConExample.md's minimal 3-piece counterexample to [JOGO] Theorem 6 / [COAP]
        % section 4: a triangle (0,0),(60,10),(-5,10) cut by two cevians into pieces 1,2,3, all
        % positive definite, with pieces 1 and 3 NON-adjacent -- so f*'s edge between them is a
        % genuine ELLIPSE (b^2-4ac=-71/2254 there), not a parabola. Historically this input DIED
        % in conj('cplq'): first a PLQ:conjCPLQ:cplqFailed coverage gap traced (2026-08-30) to a
        % single piece's own vertex-cone collapsing to a point (conjConvexOverPiece/edgeDirsAt),
        % then -- after that fix -- a MATLAB:badsubscript in symbolicFunction.tangent (a
        % degenerate conic missing one ambient variable). Both fixed the same day; this pins the
        % result correct, not just "does not crash".
        %
        % f6 rows are QuaPol's raw-Hessian convention [Q11 Q12 Q22 beta1 beta2 gamma], matching
        % CONJ_FIELD_PROOF.md 4.1 / doc/QuaConExample.md section 2 exactly (pieces 4,5 of the
        % five-piece example deleted).
            V = [0 0; 60 10; 15 10; 5 10; -5 10];
            E = [1 2 1; 1 3 1; 1 4 1; 1 5 1; 2 3 1; 3 4 1; 4 5 1];
            F = [1 0; 2 1; 3 2; 0 3; 1 0; 2 0; 3 0];
            f6 = [ 3 -1  5  0  1 0;
                   3 -5 17  6 -8 0;
                  15 -2 11  2 -6 0];
            q = QuaPol(V, E, f6, F);
            g = q.conj('cplq');

            % doc/QuaConExample.md section 2: the exact triple-value point on {g1=g3}, where
            % pieces 1 and 3 tie and piece 2 is strictly below -- f*(s*) = 2.9278190688.
            sstar = double([-17/62 + sqrt(sym(612030))/186, sym(3)/2]);
            got = conjCPLQTest.evalConjResult(g, sstar);
            testCase.verifyEqual(got, 2.9278190688, 'AbsTol', 1e-9);

            % A random sweep against the SAME per-piece oracle every other test in this file
            % uses (supQuadOverPoly, max over the 3 pieces) -- no dependency on the pipeline.
            W = {V([1 2 3],:), V([1 3 4],:), V([1 4 5],:)};
            rng(7);
            for t = 1:30
                s = (rand(1,2)-0.5)*300;
                ref = -inf;
                for k = 1:3
                    ref = max(ref, conjCPLQTest.supQuadOverPoly(s, f6(k,:), W{k}));
                end
                got = conjCPLQTest.evalConjResult(g, s);
                testCase.verifyEqual(got, ref, 'RelTol', 1e-6, 'AbsTol', 1e-9, ...
                    sprintf('at s=(%g,%g)', s(1), s(2)));
            end
        end
    end

    methods (Static)
        function v = evalConjResult(g, s)
        % Evaluate a conj('cplq') result whatever route produced it. A QuaParCPLQ wraps a cPLQ
        % functionNDomain array and is read with evalFunctionNDomain; a QuaPar/QuaPol is a mesh
        % and has its own eval. Since Case B2 the same input family can produce either, so a test
        % that is about VALUES must not hard-code one of the two.
            if isa(g, 'QuaParCPLQ')
                v = evalFunctionNDomain(g.fnd, s);
            else
                v = g.eval(s);
                if ~isfinite(v), v = NaN; end   % uncovered reads as NaN either way
            end
        end

        function q = polygonQuaPol(W)
        % A single convex polygonal face from its CCW vertex list.
            n = size(W,1);
            E = zeros(n,3); F = zeros(n,2);
            for i = 1:n
                E(i,:) = [i, mod(i,n)+1, 1];
                F(i,:) = [1, 0];
            end
            q = QuaPol(W, E, [0 0 0 0 0 0], F);
        end

        function v = supQuadOverPoly(s, f6, W)
        % EXACT sup of <s,x> - q(x) over the polygon, for a CONCAVE-or-affine integrand (i.e. any
        % q that is convex, affine, or concave). Candidates are the vertices, the per-edge
        % stationary points and the interior stationary point -- no sampling anywhere.
            Q = [f6(1) f6(2); f6(2) f6(3)]; L = [f6(4); f6(5)];
            n = size(W,1); v = -inf;
            cand = W;
            if abs(det(Q)) > 1e-12
                cand = [cand; (Q \ (s' - L))'];
            end
            for e = 1:n
                a = W(e,:); b = W(mod(e,n)+1,:); d = b - a;
                den = d*Q*d';
                if abs(den) > 1e-14
                    tt = (s*d' - (Q*a' + L)'*d') / den;
                    if tt > 0 && tt < 1, cand = [cand; a + tt*d]; end %#ok<AGROW>
                end
            end
            for k = 1:size(cand,1)
                p = cand(k,:);
                if ~conjCPLQTest.inPolygon(p, W), continue, end
                z = 0.5*f6(1)*p(1)^2 + f6(2)*p(1)*p(2) + 0.5*f6(3)*p(2)^2 ...
                    + f6(4)*p(1) + f6(5)*p(2) + f6(6);
                v = max(v, s*p' - z);
            end
        end

        function tf = inPolygon(p, W)
            n = size(W,1); tf = true;
            for i = 1:n
                a = W(i,:); b = W(mod(i,n)+1,:);
                if (b(1)-a(1))*(p(2)-a(2)) - (b(2)-a(2))*(p(1)-a(1)) < -1e-9
                    tf = false; return
                end
            end
        end

        function v = fourConeTruth(s)
        % f* of the 4-cone fan of step3UnboundedAssemblyMatchesTheTruth, in closed form. Each row
        % is {quadrant signs, a, b} for q = a*x^2 + b*y^2 on that quadrant; see that test's header.
            C = { [ 1  1], 0.5, 0.5     % face 1, first quadrant,  q = (x^2+y^2)/2
                  [-1  1], 0.5, 1.0     % face 2, second,          q = x^2/2 + y^2
                  [-1 -1], 1.0, 1.0     % face 3, third,           q = x^2 + y^2
                  [ 1 -1], 1.0, 0.5 };  % face 4, fourth,          q = x^2 + y^2/2
            v = -inf;
            for k = 1:size(C,1)
                sg = C{k,1}; a = C{k,2}; b = C{k,3};
                t1 = 0; if sg(1)*s(1) >= 0, t1 = s(1)^2/(4*a); end
                t2 = 0; if sg(2)*s(2) >= 0, t2 = s(2)^2/(4*b); end
                v = max(v, t1 + t2);
            end
        end

        function h = supBilinearOverPoly(s, T)
            % Exact sup_{(x,y) in T} [s1 x + s2 y - x y] over a triangle T: the Hessian of the
            % objective is indefinite (eigenvalues +-1), so no interior point can be a local max,
            % and the sup is attained on T's boundary -- checked in closed form (quadratic-in-t
            % along each edge). Same construction as maxQuaParTest.m's own ground-truth helper.
            s1 = s(1); s2 = s(2);
            best = -inf;
            n = size(T,1);
            for i = 1:n
                va = T(i,:); vb = T(mod(i,n)+1,:);
                dx = vb(1)-va(1); dy = vb(2)-va(2);
                A = -dx*dy;
                B = s1*dx + s2*dy - va(1)*dy - va(2)*dx;
                C = s1*va(1) + s2*va(2) - va(1)*va(2);
                cand = [0 1];
                if abs(A) > 1e-14
                    tstar = -B/(2*A);
                    if tstar > 0 && tstar < 1, cand(end+1) = tstar; end %#ok<AGROW>
                end
                for t = cand
                    val = A*t^2 + B*t + C;
                    if val > best, best = val; end
                end
            end
            h = best;
        end
    end
end
