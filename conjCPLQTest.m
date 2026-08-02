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

            % ...and a convex triangle piece is its own biconjugate, on the same route.
            caseBconvex = QuaPol([0 0; 1 0; 0 1], E3, [1 0 1 0 0 0], F3);
            bBc = caseBconvex.biconj();
            S = [0.2 0.2; 0.6 0.2; 1/3 1/3];
            testCase.verifyEqual(bBc.kind(), 'RatPol');
            testCase.verifyEqual(bBc.eval(S), caseBconvex.eval(S), 'AbsTol', 1e-12);

            % Case C -- general bounded multi-face domain: still the literal double conjugation,
            % through cPLQ's own symbolic machinery (QuaParCPLQ.conj).
            V = [0 0; 1 0; 1 1; 0 1];
            E = [1 2 1; 2 3 1; 1 3 1; 3 4 1; 4 1 1];
            F = [1 0; 1 0; 2 1; 2 0; 2 0];
            caseC = QuaPol(V, E, [1 0 1 0 0 0; 1 0 1 0 0 0], F);
            testCase.verifyEqual(caseC.conj().kind(), 'QuaParCPLQ');

            % Case C's BICONJUGATE does not work, and this assertion used to hide that. It read
            % `verifyEqual(caseC.biconj().kind(), 'QuaParCPLQ')`, which passes on an EMPTY piece
            % list -- QuaParCPLQ(functionNDomain.empty()).kind() is still 'QuaParCPLQ'. Measured
            % on pristine HEAD (2026-07-31): caseC.conj() gives 9 pieces, caseC.biconj() gives
            % ZERO, i.e. f** = +inf everywhere, for an f that is convex and hence its own
            % biconjugate.
            %
            % It fails in conjugateOfPiecePoly, behind a CHAIN of latent bugs each of which only
            % becomes reachable once the one before it is fixed. Two are now fixed --
            % region.getNormalConeVertexQ indexed py(1) before its own isempty(py) guard (dead
            % code), and region.splitmax3 left its output unassigned when f1 < f2 at every vertex
            % -- and the next one down is functionNDomain.getInterior, which indexes c2(2) under
            % a guard that only tests size(c1,2). None of this is caused by the unbounded or
            % general-quadratic work; the first conjugation is now RICHER (11 pieces rather than
            % 9), which is what carries the second one far enough to reach these.
            %
            % Pinned as "errors" rather than by identifier, because the identifier moves as the
            % chain is peeled. The invariant that matters is that it does not silently return a
            % wrong f**. See SUPPORT_MATRIX.md section 7.
            testCase.verifyError(@() caseC.biconj(), ?MException);
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

        function indefiniteQuadraticNotImplemented(testCase)
            % f = xy is indefinite; its conjugate is not a full-domain quadratic (QuaPar, TODO).
            p = QuaPol([0 1 0 0 0 0]);
            testCase.verifyError(@() p.conj('cplq'), 'PLQ:conjCPLQ:notImplemented');
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
            testCase.verifyEqual(g.kind(), 'QuaParCPLQ');

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

        function indefiniteTriangleTwoConvexEdgesSplitViaCPLQStep2(testCase)
            % The 2-convex-edge tightness split (convEnvCPLQTest.bilinearTwoConvexEdgesSplitIsTight)
            % leaves Step 1 with a RATIONAL face, which CCA2's own Step 2 cannot conjugate. It now
            % falls back to cPLQ's symbolic Step 2/3 on that envelope, so conj closes end to end.
            %
            % Note WHICH half of cPLQ is reused: its Step 2/3, on CCA2's Step 1 output. Running
            % cPLQ end to end on this triangle instead gives the WRONG answer, because cPLQ's own
            % envelope for a 2-convex-edge triangle is the single, untight Appendix A.4 formula --
            % it leaves sBad below covered by no region at all. CCA2's Step 1 is ahead of cPLQ's
            % here, so the working combination is CCA2 Step 1 + cPLQ Step 2/3.
            V = [2 1; 0 0; 1 0]; E = [1 2 1; 2 3 1; 3 1 1]; F = [1 0; 1 0; 1 0];
            q = QuaPol(V, E, [0 1 0 0 0 0], F);           % f = xy
            testCase.verifyEqual(convEnvCPLQ(q).nf, 2);   % the tightness split
            g = q.conj('cplq');
            testCase.verifyEqual(g.kind(), 'QuaParCPLQ');

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
                    got = evalFunctionNDomain(g.fnd, S(t,:));
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
                got = evalFunctionNDomain(g.fnd, S(t,:));
                testCase.verifyEqual(got, best, 'AbsTol', 1e-9, ...
                    sprintf('at s=(%g,%g)', S(t,1), S(t,2)));
            end
        end

        function step3DropsCellsOnSomeUnboundedAssemblies(testCase)
        % THE REMAINING BLOCKER, pinned. The SAME 4-cone geometry as above, differing only in
        % which quadratic sits on which cone, makes cPLQ's cross-piece maximum disagree with its
        % own per-piece conjugates. Still true, but the FAILURE HAS MOVED (2026-08-02) and the
        % name of this test is now half right -- read this before working on it.
        %
        % WAS: the assembled maximum kept only 4 of the 16 cells, LOSING face 1's s_2^2/2 cell on
        % {s1<=0, s2>=0}, so f*(-0.5,2) came back 1.125 for a truth of 2. Cause: splitting that
        % quadrant on s2^2/2 = s1^2/2 + s2^2/4 produced the half {s1<=0, s2>=0, s1^2/2-s2^2/4<=0}
        % -- a genuine 2-D cone containing (-0.5,2), (-0.1,3), (-1,4) -- and
        % region.simplifyUnboundedRegion declared it EMPTY, because it decides that from probe
        % directions built out of constraint SLOPES at a vertex and the split conic's gradient
        % vanishes at exactly that vertex. Fixed by region.witnessAwayFrom, which refutes an
        % emptiness verdict with an actual feasible point (sound: an empty region has none).
        %
        % IS: 8 cells assemble and f*(-0.5,2) is 2, correct, along with 7 other probes. What
        % assertStep3MatchesPieces now catches is the OPPOSITE error at a different point --
        % s = (-3,-2.4), where the assembly gives 5.130 and the per-piece max gives 4.500. The
        % per-piece value is the right one (by hand: the four cone suprema are 0, 4.5, 3.69 and
        % 2.88). 5.13 = s1^2/4 + s2^2/2 there, which is face 4's cell -- and face 4's cell should
        % live on {s1>=0, s2<=0}, so a region is claiming territory across s1 = 0. That is an
        % OVER-claim, not a drop; the next session should start from which region grew.
        %
        % The defect is DATA-DEPENDENT, not universal, which is exactly why the cross-check
        % matters: without it Step 3 returns plausible numbers on the cases it gets wrong.
        % assertStep3MatchesPieces compares the assembled maximum against the pointwise max of
        % the per-piece conjugates -- the same f*, computed the other way.
            V = [0 0;-1 0; 0 1;1 0;0 -1];
            E = [1 2 0;1 3 0;1 4 0;1 5 0];
            f = [1 0 1 0 0 0;1 0 2 0 0 0;2 0 2 0 0 0;2 0 1 0 0 0];
            F = [3 2;2 1;1 4;4 3];
            p = QuaPol(V,E,f,F);
            testCase.verifyError(@() p.conj('cplq'), 'PLQ:conjCPLQ:cplqFailed');
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
            testCase.verifyClass(g, 'QuaParCPLQ');

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
    end

    methods (Static)
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
