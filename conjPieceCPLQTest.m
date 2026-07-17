classdef conjPieceCPLQTest < matlab.unittest.TestCase
    % Tests for conjPieceCPLQ (Step 2): conjugate of a convex quadratic over a triangle -> QuaPar.
    % For each dual region the maximizer x* is known, so f*(s) = <s,x*> - q(x*) gives an exact check.

    methods (Test)
        function sevenRegionConjugateValues(testCase)
            % q = x^2 + y^2 over the triangle (0,0),(1,0),(0,1).
            A = [2 0; 0 2]; b = [0; 0]; cc = 0;
            V = [0 0; 1 0; 0 1]; E = [1 2 1; 2 3 1; 3 1 1]; F = [1 0; 1 0; 1 0];
            p = QuaPoly(V, E, [2 0 2 0 0 0], F);          % x^2 + y^2
            g = conjPieceCPLQ(p);
            testCase.verifyClass(g, 'QuaPar');
            testCase.verifyEqual(g.nf, 7);

            qf = @(x) 0.5*x'*A*x + b'*x + cc;
            gq = @(x) A*x + b;
            un = @(e) [e(2); -e(1)] / norm([e(2); -e(1)]);
            v1 = [0;0]; v2 = [1;0]; v3 = [0;1];
            n1 = un(v2-v1); n2 = un(v3-v2); n3 = un(v1-v3);

            X = {}; S = {};
            X{1} = [0.2;0.2];      S{1} = gq(X{1});                 % central region F0
            X{2} = v1;             S{2} = gq(v1) + 5*(n1+n3);       % cone C1
            X{3} = v2;             S{3} = gq(v2) + 5*(n1+n2);       % cone C2
            X{4} = v3;             S{4} = gq(v3) + 5*(n2+n3);       % cone C3
            X{5} = (v1+v2)/2;      S{5} = gq(X{5}) + 5*n1;          % strip S1
            X{6} = (v2+v3)/2;      S{6} = gq(X{6}) + 5*n2;          % strip S2
            X{7} = (v3+v1)/2;      S{7} = gq(X{7}) + 5*n3;          % strip S3
            for i = 1:7
                expected = S{i}'*X{i} - qf(X{i});
                testCase.verifyEqual(g.eval(S{i}'), expected, 'AbsTol', 1e-10, ...
                    sprintf('region %d', i));
            end
        end

        function conjugateFiniteEverywhere(testCase)
            V = [0 0; 1 0; 0 1]; E = [1 2 1; 2 3 1; 3 1 1]; F = [1 0; 1 0; 1 0];
            g = conjPieceCPLQ(QuaPoly(V, E, [2 0 2 0 0 0], F));
            vals = g.eval([100 100; -100 50; 0 0; -30 -70]);
            testCase.verifyTrue(all(isfinite(vals)));   % conjugate of a piece over a bounded set
        end

        function generalConvexQuadraticTriangle(testCase)
            % q = 1/2 x'A x + b'x + c with A=[2 1;1 3] (PD), b=[1;-2], c=0.5, over a triangle.
            A = [2 1; 1 3]; b = [1; -2]; cc = 0.5;
            V = [0 0; 2 0; 1 2]; E = [1 2 1; 2 3 1; 3 1 1]; F = [1 0; 1 0; 1 0];
            f6 = [A(1,1) A(1,2) A(2,2) b(1) b(2) cc];     % stored: matrixForm reads Q=[..],L=b,const
            g = conjPieceCPLQ(QuaPoly(V, E, f6, F));
            testCase.verifyEqual(g.nf, 7);
            qf = @(x) 0.5*x'*A*x + b'*x + cc;
            gq = @(x) A*x + b;
            un = @(e) [e(2); -e(1)] / norm([e(2); -e(1)]);
            v1 = [0;0]; v2 = [2;0]; v3 = [1;2];
            n1 = un(v2-v1); n2 = un(v3-v2); n3 = un(v1-v3);
            X = {[0.8;0.5], v1, v2, v3, (v1+v2)/2, (v2+v3)/2, (v3+v1)/2};
            S = {gq(X{1}), gq(v1)+8*(n1+n3), gq(v2)+8*(n1+n2), gq(v3)+8*(n2+n3), ...
                 gq(X{5})+8*n1, gq(X{6})+8*n2, gq(X{7})+8*n3};
            for i = 1:7
                testCase.verifyEqual(g.eval(S{i}'), S{i}'*X{i} - qf(X{i}), 'AbsTol', 1e-9, ...
                    sprintf('region %d', i));
            end
        end

        function linearOverTriangleConjugate(testCase)
            % ell(x) = -x over the triangle (0,0),(1,0),(0,1). Vertex values ell = 0,-1,0, so
            % f*(s) = max_i(<s,v_i> - ell(v_i)) = max(0, s1+1, s2), a 3-cone piecewise-linear QuaPar.
            V = [0 0; 1 0; 0 1]; E = [1 2 1; 2 3 1; 3 1 1]; F = [1 0; 1 0; 1 0];
            g = conjPieceCPLQ(QuaPoly(V, E, [0 0 0 -1 0 0], F));   % ell = -x
            testCase.verifyClass(g, 'QuaPar');
            testCase.verifyEqual(g.nf, 3);
            S = [-2 -1; 3 -1; -1 2; 0.5 0.5; -1 0];
            for i = 1:size(S,1)
                expected = max([0; S(i,1)+1; S(i,2)]);        % max(L1,L2,L3)
                testCase.verifyEqual(g.eval(S(i,:)), expected, 'AbsTol', 1e-10, sprintf('s=%d', i));
            end
        end

        function concaveQuadraticConjugateEndToEnd(testCase)
            % conj(concave quadratic over T) = conj(its affine envelope). For q = -(x^2+y^2) over
            % the triangle (0,0),(2,0),(0,2): sup_{x in T} <s,x> - q(x) is attained at a vertex, so
            % f*(s) = max_i(<s,v_i> - q(v_i)) with q-values 0,-4,-4.
            V = [0 0; 2 0; 0 2]; E = [1 2 1; 2 3 1; 3 1 1]; F = [1 0; 1 0; 1 0];
            q = QuaPoly(V, E, [-2 0 -2 0 0 0], F);            % -(x^2+y^2)
            r = convEnvCPLQ(q);                                % Step 1 -> affine envelope (RatPol)
            g = conjPieceCPLQ(r);                              % Step 2 -> 3-cone conjugate
            testCase.verifyEqual(g.nf, 3);
            qf = @(x) -(x(1)^2 + x(2)^2);
            v1=[0;0]; v2=[2;0]; v3=[0;2];
            S = [1 1; 5 -1; -1 5; 3 3];
            for i = 1:size(S,1)
                s = S(i,:)';
                expected = max([s'*v1-qf(v1); s'*v2-qf(v2); s'*v3-qf(v3)]);
                testCase.verifyEqual(g.eval(S(i,:)), expected, 'AbsTol', 1e-10, sprintf('s=%d', i));
            end
        end

        function bilinearOneConvexEdgeConjugate(testCase)
            % f = xy over a triangle with ONE convex edge -> parabolic QuaPar (6 faces). Validate
            % against the numeric sup f*(s) = max_{x in T} <s,x> - xy over a fine triangle grid.
            V = [0 0; 2 0; 1 1];                          % CCW; one convex edge (1,1)-(0,0), m=1
            E = [1 2 1; 2 3 1; 3 1 1]; F = [1 0; 1 0; 1 0];
            g = conjPieceCPLQ(QuaPoly(V, E, [0 1 0 0 0 0], F));
            testCase.verifyClass(g, 'QuaPar');
            testCase.verifyEqual(g.nf, 6);
            nt = 200; [uu,vv] = meshgrid(linspace(0,1,nt)); uu = uu(:); vv = vv(:);
            kk = (uu+vv <= 1); uu = uu(kk); vv = vv(kk);
            Xg = V(1,1)+uu*(V(2,1)-V(1,1))+vv*(V(3,1)-V(1,1));
            Yg = V(1,2)+uu*(V(2,2)-V(1,2))+vv*(V(3,2)-V(1,2));
            xyg = Xg.*Yg;
            S = [0.5 0.5; 3 -1; -2 3; 1 1; 0 -3; 4 4; -3 -3; 0.5 1.5];
            for i = 1:size(S,1)
                sup = max(S(i,1)*Xg + S(i,2)*Yg - xyg);
                testCase.verifyEqual(g.eval(S(i,:)), sup, 'AbsTol', 1e-3, sprintf('s=%d', i));
            end
        end

        function bilinearZeroConvexEdgeConjugate(testCase)
            % f = xy over (0,0),(1,0),(0,1): no convex edge (slopes 0,-1,vertical) -> vertex values
            % xy = 0,0,0, so f*(s) = max(0, s1, s2), a 3-cone piecewise-linear QuaPar. Cross-check
            % against the numeric sup over the triangle too.
            V = [0 0; 1 0; 0 1]; E = [1 2 1; 2 3 1; 3 1 1]; F = [1 0; 1 0; 1 0];
            g = conjPieceCPLQ(QuaPoly(V, E, [0 1 0 0 0 0], F));
            testCase.verifyClass(g, 'QuaPar');
            testCase.verifyEqual(g.nf, 3);
            S = [-2 -1; 3 -1; -1 2; 0.5 0.5; -1 -1];
            for i = 1:size(S,1)
                testCase.verifyEqual(g.eval(S(i,:)), max([0, S(i,1), S(i,2)]), ...
                    'AbsTol', 1e-10, sprintf('s=%d', i));
            end
            nt = 200; [uu,vv] = meshgrid(linspace(0,1,nt)); uu = uu(:); vv = vv(:);
            kk = (uu+vv <= 1); uu = uu(kk); vv = vv(kk);
            Xg = V(1,1)+uu*(V(2,1)-V(1,1))+vv*(V(3,1)-V(1,1));
            Yg = V(1,2)+uu*(V(2,2)-V(1,2))+vv*(V(3,2)-V(1,2));
            xyg = Xg.*Yg;
            for i = 1:size(S,1)
                sup = max(S(i,1)*Xg + S(i,2)*Yg - xyg);
                testCase.verifyEqual(g.eval(S(i,:)), sup, 'AbsTol', 1e-3, sprintf('numeric s=%d', i));
            end
        end

        function bilinearTwoConvexEdgesDocumentedLimitation(testCase)
            % xy over (0,0),(2,1),(1,2): two convex edges (m=1/2 and m=2, sharing vertex (0,0)).
            % This is a DELIBERATE, DOCUMENTED non-implementation, not a TODO: the conjugate's
            % dual-space arrangement needs a boundary between the two edges' own quadratic
            % conjugate formulas (COAP B.2 applied to each edge separately), and that boundary is
            % a genuine HYPERBOLA whenever the two slopes differ -- verified both numerically
            % (0/2e6 mismatches across two test triangles when the true 2D sup was compared
            % against max(vertex linears, gated per-edge quadratics)) and symbolically: the
            % quadratic parts of the two edges' formulas share b=1/2 always, so their difference's
            % discriminant is -4*(1/(4m1)-1/(4m2))*(m1/4-m2/4) = (m1-m2)^2/(4*m1*m2) > 0 whenever
            % m1~=m2. QuaPar only supports parabolic/linear (degenerate) conics per edge, so this
            % case cannot be represented as constructed. It also does not arise from the wired
            % pipeline: Step 1 (convEnvCPLQ) always convexifies a 2-convex-edge piece into a
            % rank-1 PSD quadratic (see conjPieceCPLQTest/psdRank1QuadraticEndToEnd) before Step 2
            % ever sees it, so this raw-indefinite-input case is never actually hit.
            V = [0 0; 2 1; 1 2]; E = [1 2 1; 2 3 1; 3 1 1]; F = [1 0; 1 0; 1 0];
            q = QuaPoly(V, E, [0 1 0 0 0 0], F);
            testCase.verifyError(@() conjPieceCPLQ(q), 'conjPieceCPLQ:notImplemented');
        end

        function psdRank1QuadraticConjugate(testCase)
            % q = 1/2 x'Ax + b'x + c, A = lam*u*u' rank-1 PSD (u a unit vector, not axis-aligned),
            % over a generic triangle -> six-face QuaPar (three parabolic strips + three linear
            % vertex cones). Validated against the numeric sup over a fine triangle grid.
            u = [cos(0.7); sin(0.7)]; lam = 1.3;
            A = lam*(u*u'); b = [0.4; -0.6]; cc = 0.2;
            V = [1 1; 4 3; 3 5]; E = [1 2 1; 2 3 1; 3 1 1]; F = [1 0; 1 0; 1 0];
            f6 = [A(1,1) A(1,2) A(2,2) b(1) b(2) cc];
            g = conjPieceCPLQ(QuaPoly(V, E, f6, F));
            testCase.verifyClass(g, 'QuaPar');
            testCase.verifyEqual(g.nf, 6);
            nt = 220; [uu,vv] = meshgrid(linspace(0,1,nt)); uu = uu(:); vv = vv(:);
            kk = (uu+vv <= 1); uu = uu(kk); vv = vv(kk);
            Xg = V(1,1)+uu*(V(2,1)-V(1,1))+vv*(V(3,1)-V(1,1));
            Yg = V(1,2)+uu*(V(2,2)-V(1,2))+vv*(V(3,2)-V(1,2));
            qg = A(1,1)*Xg.^2 + 2*A(1,2)*Xg.*Yg + A(2,2)*Yg.^2; % 1/2*x'Ax with A symmetric plain form
            qg = 0.5*qg + b(1)*Xg + b(2)*Yg + cc;
            S = [0.5 0.5; 3 -1; -2 3; 1 1; 0 -3; 4 4; -3 -3; 6 2; -1 6];
            for i = 1:size(S,1)
                sup = max(S(i,1)*Xg + S(i,2)*Yg - qg);
                testCase.verifyEqual(g.eval(S(i,:)), sup, 'AbsTol', 2e-3, sprintf('s=%d', i));
            end
        end

        function psdRank1QuadraticEndToEnd(testCase)
            % conj(Step 1 envelope of xy over a 2-convex-edge triangle) end to end. Step 1
            % (convEnvCPLQ) turns the indefinite xy into TWO rank-1 PSD quadratic sub-triangle
            % pieces (COAP A.4, corrected -- the single-quadratic formula alone is not always the
            % tightest envelope over the whole triangle; see convEnvCPLQ.m's splitTwoConvexEdges
            % HISTORY / DESIGN.md). Extracts face 1 (the sub-triangle containing the two convex
            % edges' shared vertex, where the ORIGINAL single-quadratic formula is unchanged) as a
            % standalone triangle and feeds it to Step 2, which must dispatch to
            % conjPSDRank1QuadTriangle (generic) or its tied-vertex sub-case
            % (conjPSDRank1QuadTriangleTie), not error out. Checked against the numeric sup of that
            % sub-triangle's quadratic itself (not the original xy).
            % Uses (1,1),(4,3),(3,5), also used (whole, unsplit) by psdRank1QuadraticConjugate to
            % exercise the generic 6-face path directly. With the corrected split-cevian location
            % (this session's Part 2 fix -- the OLD, buggy cevian placed the seam at roughly half
            % the correct distance from the shared vertex, see DESIGN.md), face 1's sub-triangle
            % shape changed enough that it now hits the TIED-vertex sub-case here (5 faces, not 6)
            % -- confirmed still numerically exact against ground truth below; the 5-face path
            % itself is exercised independently (from a different, mirror-symmetric starting
            % triangle) by psdRank1QuadraticTieEndToEnd.
            V = [1 1; 4 3; 3 5]; E = [1 2 1; 2 3 1; 3 1 1]; F = [1 0; 1 0; 1 0];
            q = QuaPoly(V, E, [0 1 0 0 0 0], F);          % xy, two convex edges
            r = convEnvCPLQ(q);                            % Step 1 -> 2 rank-1 PSD sub-triangles
            testCase.verifyEqual(r.nf, 2);
            iVs = unique(r.E(r.F(:,1)==1 | r.F(:,2)==1, 1:2));
            V1 = r.V(iVs, :);
            [L, Q, C] = QuaPoly.matrixForm(r.f(1,:));
            testCase.verifyEmpty(C);
            testCase.verifyEqual(min(eig(Q)), 0, 'AbsTol', 1e-8);   % rank-1 PSD, as proven
            p1 = QuaPoly(V1, E, r.f(1,5:10), F);
            g = conjPieceCPLQ(p1);                          % Step 2
            testCase.verifyClass(g, 'QuaPar');
            testCase.verifyEqual(g.nf, 5);
            nt = 220; [uu,vv] = meshgrid(linspace(0,1,nt)); uu = uu(:); vv = vv(:);
            kk = (uu+vv <= 1); uu = uu(kk); vv = vv(kk);
            Xg = V1(1,1)+uu*(V1(2,1)-V1(1,1))+vv*(V1(3,1)-V1(1,1));
            Yg = V1(1,2)+uu*(V1(2,2)-V1(1,2))+vv*(V1(3,2)-V1(1,2));
            qg = Q(1,1)*Xg.^2 + 2*Q(1,2)*Xg.*Yg + Q(2,2)*Yg.^2;
            qg = 0.5*qg + L(1)*Xg + L(2)*Yg + r.f(1,end);
            S = [0.5 0.5; 3 -1; -2 3; 1 1; 0 -3; 4 4; -3 -3];
            for i = 1:size(S,1)
                sup = max(S(i,1)*Xg + S(i,2)*Yg - qg);
                testCase.verifyEqual(g.eval(S(i,:)), sup, 'AbsTol', 2e-3, sprintf('s=%d', i));
            end
        end

        function pickRepFindsThinEdgeStripFace(testCase)
            % Regression test for a bug found via a randomized f(x,y)=xy triangle stress test:
            % conjPSDRank1QuadTriangle's pickRep helper searches base+/-dir*mag over a fixed list
            % of magnitudes (scale*[1,0.3,...,0.001,3,10]) to find a point inside each candidate
            % face; Step 1's triangle split produces a sub-triangle so thin that its E23
            % edge-strip face needed a magnitude one order of magnitude below the smallest value
            % the old list tried (scale*0.001) -- so conjPieceCPLQ threw conjPieceCPLQ:internal
            % ("could not locate a representative point for face 3"). Fixed by replacing the fixed
            % list with a dense geometric sweep (see pickRep's HISTORY).
            %
            % T was originally (0,0),(4.46,1.83),(5.81,2.38): after splitThreeConvex's split-point
            % formula was itself corrected (see convEnvCPLQ.m's HISTORY -- the old horizontal-
            % through-vmid split was mathematically wrong in general), that triangle's sub-triangle
            % sliver became thinner still (needing mag ~1e-9, not just ~1e-4) and started tripping
            % a SEPARATE, deeper, still-open bug several layers downstream (QuaPar.orderEdges
            % rejecting the resulting face topology -- see session handoff notes on near-degenerate
            % sliver triangles). Swapped in this triangle, which still produces a thin (~0.027,
            % same order of magnitude as the original 0.0074) edge-strip sub-triangle exercising
            % pickRep's search exactly as before, without hitting that separate unresolved issue.
            T = [3.1436 2.4929; 5.0857 4.1038; 9.0757 7.5555];
            E = [1 2 1; 2 3 1; 3 1 1]; F = [1 0; 1 0; 1 0];
            q = QuaPoly(T, E, [0 1 0 0 0 0], F);           % f(x,y) = x*y over T
            r = convEnvCPLQ(q);                             % Step 1: 2 sub-triangles
            testCase.verifyEqual(r.nf, 2);
            V1 = r.V(unique(r.E(r.F(:,1)==1 | r.F(:,2)==1, 1:2)), :);
            f1 = r.f(1,:);
            testCase.verifyLessThan(min(abs(diff(sort(V1(:,1))))), 0.05, ...
                'expected the thin sliver sub-triangle this bug depends on');

            p1 = QuaPoly(V1, E, f1(5:10), F);
            g = conjPieceCPLQ(p1);   % used to throw conjPieceCPLQ:internal here
            testCase.verifyClass(g, 'QuaPar');
            testCase.verifyEqual(g.nf, 6);

            nt = 220; [uu,vv] = meshgrid(linspace(0,1,nt)); uu = uu(:); vv = vv(:);
            kk = (uu+vv <= 1); uu = uu(kk); vv = vv(kk);
            Xg = V1(1,1)+uu*(V1(2,1)-V1(1,1))+vv*(V1(3,1)-V1(1,1));
            Yg = V1(1,2)+uu*(V1(2,2)-V1(1,2))+vv*(V1(3,2)-V1(1,2));
            qg = QuaPoly.evalPoly(f1, [Xg Yg]);
            S = [1 1; 0.5 0.2; 2 1; -1 1];
            for i = 1:size(S,1)
                sup = max(S(i,1)*Xg + S(i,2)*Yg - qg);
                testCase.verifyEqual(g.eval(S(i,:)), sup, 'AbsTol', 1e-4, sprintf('s=%d', i));
            end
        end

        function psdRank1QuadraticTieConjugate(testCase)
            % Same construction as psdRank1QuadraticConjugate but with va,vb manufactured to tie
            % in the rotated-t coordinate (the triangle edge va-vb is exactly perpendicular to A's
            % nonzero eigenvector u) -- dispatches to conjPSDRank1QuadTriangleTie, five faces (two
            % parabolic edge strips + three linear vertex cones) instead of six.
            u = [cos(0.7); sin(0.7)]; lam = 1.3; uperp = [-u(2); u(1)];
            A = lam*(u*u'); b = [0.4; -0.6]; cc = 0.2;
            v1 = [0;0]; t0 = 2;
            v2 = t0*u + 1.0*uperp; v3 = t0*u - 1.5*uperp;   % v2,v3 tie at t=t0
            V = [v1'; v2'; v3'];
            if 0.5*((V(2,1)-V(1,1))*(V(3,2)-V(1,2)) - (V(3,1)-V(1,1))*(V(2,2)-V(1,2))) < 0
                V = V([1 3 2],:);
            end
            f6 = [A(1,1) A(1,2) A(2,2) b(1) b(2) cc];
            g = conjPieceCPLQ(QuaPoly(V, [1 2 1; 2 3 1; 3 1 1], f6, [1 0; 1 0; 1 0]));
            testCase.verifyClass(g, 'QuaPar');
            testCase.verifyEqual(g.nf, 5);
            nt = 220; [uu2,vv2] = meshgrid(linspace(0,1,nt)); uu2 = uu2(:); vv2 = vv2(:);
            kk = (uu2+vv2 <= 1); uu2 = uu2(kk); vv2 = vv2(kk);
            Xg = V(1,1)+uu2*(V(2,1)-V(1,1))+vv2*(V(3,1)-V(1,1));
            Yg = V(1,2)+uu2*(V(2,2)-V(1,2))+vv2*(V(3,2)-V(1,2));
            qg = A(1,1)*Xg.^2 + 2*A(1,2)*Xg.*Yg + A(2,2)*Yg.^2;
            qg = 0.5*qg + b(1)*Xg + b(2)*Yg + cc;
            S = [0.5 0.5; 3 -1; -2 3; 1 1; 0 -3; 4 4; -3 -3; 6 2; -1 6];
            for i = 1:size(S,1)
                sup = max(S(i,1)*Xg + S(i,2)*Yg - qg);
                testCase.verifyEqual(g.eval(S(i,:)), sup, 'AbsTol', 2e-3, sprintf('s=%d', i));
            end
        end

        function psdRank1QuadraticTieEndToEnd(testCase)
            % conj(Step 1 envelope of xy over (0,0),(2,1),(1,2)) end to end -- previously a
            % documented UNHANDLED case (see psdRank1QuadraticEndToEnd): the harmonic-mean
            % envelope's Hessian eigenvector is exactly (1,1)/sqrt(2), which ties the two
            % non-origin vertices in the rotated t-coordinate. Now dispatches to
            % conjPSDRank1QuadTriangleTie (5 faces).
            V = [0 0; 2 1; 1 2]; E = [1 2 1; 2 3 1; 3 1 1]; F = [1 0; 1 0; 1 0];
            q = QuaPoly(V, E, [0 1 0 0 0 0], F);          % xy, two convex edges
            r = convEnvCPLQ(q);                            % Step 1 -> rank-1 PSD quadratic
            [L, Q, C] = QuaPoly.matrixForm(r.f(1,:));
            testCase.verifyEmpty(C);
            testCase.verifyEqual(min(eig(Q)), 0, 'AbsTol', 1e-8);
            g = conjPieceCPLQ(r);                          % Step 2 -> conjPSDRank1QuadTriangleTie
            testCase.verifyClass(g, 'QuaPar');
            testCase.verifyEqual(g.nf, 5);
            nt = 220; [uu,vv] = meshgrid(linspace(0,1,nt)); uu = uu(:); vv = vv(:);
            kk = (uu+vv <= 1); uu = uu(kk); vv = vv(kk);
            Xg = V(1,1)+uu*(V(2,1)-V(1,1))+vv*(V(3,1)-V(1,1));
            Yg = V(1,2)+uu*(V(2,2)-V(1,2))+vv*(V(3,2)-V(1,2));
            qg = Q(1,1)*Xg.^2 + 2*Q(1,2)*Xg.*Yg + Q(2,2)*Yg.^2;
            qg = 0.5*qg + L(1)*Xg + L(2)*Yg + r.f(1,end);
            S = [0.5 0.5; 3 -1; -2 3; 1 1; 0 -3; 4 4; -3 -3];
            for i = 1:size(S,1)
                sup = max(S(i,1)*Xg + S(i,2)*Yg - qg);
                testCase.verifyEqual(g.eval(S(i,:)), sup, 'AbsTol', 2e-3, sprintf('s=%d', i));
            end
        end

        function indefiniteQuadraticZeroConvexEdgeConjugate(testCase)
            % q = 3*x*y + 2x + y + 0.5 (indefinite, beta=3 not 1, plus a linear shift) over
            % (0,0),(1,0),(0,1): generalizes bilinearZeroConvexEdgeConjugate (which only covered
            % the exact pure-xy, no-shift case) via conjIndefiniteQuadTriangle. In the reduced
            % bilinear frame this triangle still has zero convex edges, so the conjugate is again
            % the 3-cone piecewise-linear QuaPar of vertex values (COAP B.1), just with the actual
            % vertex values of q (not xy). Cross-checked against the numeric sup over the triangle.
            A = [0 3; 3 0]; b = [2; 1]; cc = 0.5;
            V = [0 0; 1 0; 0 1]; E = [1 2 1; 2 3 1; 3 1 1]; F = [1 0; 1 0; 1 0];
            f6 = [A(1,1) A(1,2) A(2,2) b(1) b(2) cc];
            g = conjPieceCPLQ(QuaPoly(V, E, f6, F));
            testCase.verifyClass(g, 'QuaPar');
            testCase.verifyEqual(g.nf, 3);
            qf = @(x) 0.5*x'*A*x + b'*x + cc;
            v1=[0;0]; v2=[1;0]; v3=[0;1];
            S = [-2 -1; 3 -1; -1 2; 0.5 0.5; -1 -1];
            for i = 1:size(S,1)
                s = S(i,:)';
                expected = max([s'*v1-qf(v1); s'*v2-qf(v2); s'*v3-qf(v3)]);
                testCase.verifyEqual(g.eval(S(i,:)), expected, 'AbsTol', 1e-9, sprintf('s=%d', i));
            end
            nt = 200; [uu,vv] = meshgrid(linspace(0,1,nt)); uu = uu(:); vv = vv(:);
            kk = (uu+vv <= 1); uu = uu(kk); vv = vv(kk);
            Xg = V(1,1)+uu*(V(2,1)-V(1,1))+vv*(V(3,1)-V(1,1));
            Yg = V(1,2)+uu*(V(2,2)-V(1,2))+vv*(V(3,2)-V(1,2));
            qg = A(1,1)*Xg.^2 + 2*A(1,2)*Xg.*Yg + A(2,2)*Yg.^2;
            qg = 0.5*qg + b(1)*Xg + b(2)*Yg + cc;
            for i = 1:size(S,1)
                sup = max(S(i,1)*Xg + S(i,2)*Yg - qg);
                testCase.verifyEqual(g.eval(S(i,:)), sup, 'AbsTol', 2e-3, sprintf('numeric s=%d', i));
            end
        end

        function indefiniteQuadraticOneConvexEdgeConjugate(testCase)
            % q = 1/2 x'Ax+b'x+c with A=[2 1;1 -2] indefinite (rotated, not axis- or diagonal-
            % aligned) and a nonzero shift b, over a triangle -- generalizes
            % bilinearOneConvexEdgeConjugate to a genuinely rotated/shifted indefinite quadratic
            % via conjIndefiniteQuadTriangle (bilinearFrame reduction + shift correction +
            % pushforwardQuaParDual). In the reduced frame this hits exactly one convex edge, so
            % the conjugate is a six-face parabolic QuaPar, validated against the numeric sup.
            A = [2 1; 1 -2]; b = [0.5; 0.5]; cc = 0.2;
            V = [0 0; 2 0; 1 2]; E = [1 2 1; 2 3 1; 3 1 1]; F = [1 0; 1 0; 1 0];
            f6 = [A(1,1) A(1,2) A(2,2) b(1) b(2) cc];
            g = conjPieceCPLQ(QuaPoly(V, E, f6, F));
            testCase.verifyClass(g, 'QuaPar');
            testCase.verifyEqual(g.nf, 6);
            nt = 220; [uu,vv] = meshgrid(linspace(0,1,nt)); uu = uu(:); vv = vv(:);
            kk = (uu+vv <= 1); uu = uu(kk); vv = vv(kk);
            Xg = V(1,1)+uu*(V(2,1)-V(1,1))+vv*(V(3,1)-V(1,1));
            Yg = V(1,2)+uu*(V(2,2)-V(1,2))+vv*(V(3,2)-V(1,2));
            qg = A(1,1)*Xg.^2 + 2*A(1,2)*Xg.*Yg + A(2,2)*Yg.^2;
            qg = 0.5*qg + b(1)*Xg + b(2)*Yg + cc;
            S = [0.5 0.5; 3 -1; -2 3; 1 1; 0 -3; 4 4; -3 -3; 6 2; -1 6];
            for i = 1:size(S,1)
                sup = max(S(i,1)*Xg + S(i,2)*Yg - qg);
                testCase.verifyEqual(g.eval(S(i,:)), sup, 'AbsTol', 2e-3, sprintf('s=%d', i));
            end
        end

        function rationalRejected(testCase)
            V = [0 0; 1 0; 0 1]; E = [1 2 1; 2 3 1; 3 1 1]; F = [1 0; 1 0; 1 0];
            r = RatPol(V, E, [2 0 2 0 0 0], F, [1 1 3]);  % rational denominator
            testCase.verifyError(@() conjPieceCPLQ(r), 'conjPieceCPLQ:notImplemented');
        end

        function nonTriangleRejected(testCase)
            V = [0 0; 1 0; 1 1; 0 1]; E = [1 2 1; 2 3 1; 3 4 1; 4 1 1]; F = [1 0;1 0;1 0;1 0];
            q = QuaPoly(V, E, [2 0 2 0 0 0], F);
            testCase.verifyError(@() conjPieceCPLQ(q), 'conjPieceCPLQ:notImplemented');
        end
    end
end
