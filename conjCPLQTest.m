classdef conjCPLQTest < matlab.unittest.TestCase
    % Tests for the 'cplq' conjugate engine (conjCPLQ) and the QuaPoly operator interface.
    % Covers the currently implemented case (full-domain quadratics) plus the
    % validation/rejection behaviour. See DESIGN.md II.5.1.

    methods (Test)
        function energyIsSelfConjugate(testCase)
            % f = 1/2(x^2+y^2) is self-conjugate.
            p = QuaPoly.energy();
            q = p.conj('cplq');
            testCase.verifyEqual(q.f, p.f, 'AbsTol', 1e-12);
        end

        function defaultEngineIsCPLQ(testCase)
            % conj() with no engine argument defaults to 'cplq'.
            p = QuaPoly.energy();
            testCase.verifyEqual(p.conj().f, p.conj('cplq').f, 'AbsTol', 1e-12);
        end

        function biconjOfConvexQuadraticIsItself(testCase)
            % f** = f for a closed convex function.
            p = QuaPoly.energy();
            b = p.biconj('cplq');
            testCase.verifyEqual(b.f, p.f, 'AbsTol', 1e-12);
        end

        function generalPositiveDefiniteQuadratic(testCase)
            % f(x) = 1/2 x'Q x + L'x + k,  Q=[2 0;0 4], L=[1;-1], k=3.
            % f*(s) = 1/2 (s-L)' inv(Q) (s-L) - k.
            f6 = [2 0 4 1 -1 3];          % [x^2 xy y^2 x y const]
            p  = QuaPoly(f6);
            q  = p.conj('cplq');
            Q = [2 0; 0 4]; L = [1; -1]; k = 3;
            M = inv(Q); grad = -M*L; d = 0.5*(L'*M*L) - k;
            expf6 = [M(1,1) M(1,2) M(2,2) grad(1) grad(2) d];
            expq  = QuaPoly(expf6);
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
            p  = QuaPoly(f6);
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
            p = QuaPoly([0 1 0 0 0 0]);
            testCase.verifyError(@() p.conj('cplq'), 'PLQ:conjCPLQ:notImplemented');
        end

        function affineTriangleViaOrchestrator(testCase)
            % conjCPLQ dispatches a single bounded-triangle piece straight to conjPieceCPLQ (no
            % Step 1 needed for an affine piece). ell = -x over (0,0),(1,0),(0,1).
            V = [0 0; 1 0; 0 1]; E = [1 2 1; 2 3 1; 3 1 1]; F = [1 0; 1 0; 1 0];
            p = QuaPoly(V, E, [0 0 0 -1 0 0], F);
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
            p = QuaPoly(V, E, f6, F);
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
            q = QuaPoly(V, E, [-2 0 -2 0 0 0], F);
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
            g0 = QuaPoly(V0, E0, [0 1 0 0 0 0], F0).conj('cplq');
            testCase.verifyEqual(g0.nf, 3);   % zero convex edges -> 3-cone piecewise-linear

            V1 = [0 0; 2 0; 1 1]; E1 = [1 2 1; 2 3 1; 3 1 1]; F1 = [1 0; 1 0; 1 0];
            g1 = QuaPoly(V1, E1, [0 1 0 0 0 0], F1).conj('cplq');
            testCase.verifyEqual(g1.nf, 6);   % one convex edge -> 6-face parabolic QuaPar
        end

        function indefiniteTriangleTwoConvexEdgesSidestepsToEnvelope(testCase)
            % Two convex edges: conjPieceCPLQ rejects the raw piece directly, so conjCPLQ must
            % fall back to Step 1's rank-1-PSD envelope (COAP Appendix A.4) automatically, matching
            % what conjPieceCPLQTest/psdRank1QuadraticEndToEnd does by hand.
            V = [0 0; 2 1; 1 2]; E = [1 2 1; 2 3 1; 3 1 1]; F = [1 0; 1 0; 1 0];
            q = QuaPoly(V, E, [0 1 0 0 0 0], F);
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
            V = [0 0; 3 3; 1 2]; E = [1 2 1; 2 3 1; 3 1 1]; F = [1 0; 1 0; 1 0];
            q = QuaPoly(V, E, [0 1 0 0 0 0], F);
            testCase.verifyEqual(convEnvCPLQ(q).nf, 4);   % confirms the (now recursive) split
            testCase.verifyError(@() q.conj('cplq'), 'conjPieceCPLQ:notImplemented');
        end

        function multiFacePieceStillNotImplemented(testCase)
            % A multi-face domain (nf>1) still needs Step 3 (not implemented), even though the
            % single-triangle case is now handled. Reuses the known-good 4-face V/E/F geometry
            % from plqvcAliasStillWorks below (a fan of 4 unbounded cones around the origin).
            V = [0 0;-1 0; 0 1;1 0;0 -1];
            E = [1 2 0;1 3 0;1 4 0;1 5 0];
            f = [1 0 1 0 0 0;1 0 2 0 0 0;2 0 2 0 0 0;2 0 1 0 0 0];
            F = [1 2;2 3;3 4;4 1];
            p = QuaPoly(V,E,f,F);
            testCase.verifyEqual(p.nf, 4);
            testCase.verifyError(@() p.conj('cplq'), 'PLQ:conjCPLQ:notImplemented');
        end

        function cubicRejectedByOperators(testCase)
            % Cubic numerator is storable but rejected by operators (allowed for isConvex only).
            p = QuaPoly([1 0 0 0 0 0 0 0 0 0]);   % x^3/6 term present -> degree 3
            testCase.verifyEqual(p.degree, 3);
            testCase.verifyError(@() p.conj('cplq'), 'PLQ:op:unsupportedType');
        end

        function unimplementedEnginesError(testCase)
            p = QuaPoly.energy();
            testCase.verifyError(@() p.conj('pqp'),   'PLQ:conj:engine');
            testCase.verifyError(@() p.conj('graph'), 'PLQ:conj:engine');
        end

        function partialConjEngineRestriction(testCase)
            p = QuaPoly.energy();
            testCase.verifyError(@() p.partialConj(1,'graph'), 'PLQ:partialConj:engine');
        end

        function plqvcAliasStillWorks(testCase)
            % Backward compatibility: PLQVC is an alias of QuaPoly.
            p = PLQVC.energy();                 % inherited static factory
            testCase.verifyTrue(isa(p, 'QuaPoly'));
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
