classdef convEnvCPLQTest < matlab.unittest.TestCase
    % Tests for convEnvCPLQ (Step 1 of the 'cplq' pipeline: convex envelope of one quadratic piece).
    % Implemented cases: convex (envelope = self) and concave-over-triangle (affine interpolation),
    % both returning a RatPol. Indefinite/multi-piece/cubic are guarded as not-implemented/rejected.

    methods (Test)
        function convexQuadraticOverTriangleIsItself(testCase)
            % q = 1/2(x^2+y^2) (convex) over the triangle {x>=0,y>=0,x+y<=1}: envelope = q.
            V = [0 0; 1 0; 0 1]; E = [1 2 1; 2 3 1; 3 1 1]; F = [1 0; 1 0; 1 0];
            q = QuaPoly(V,E,[1 0 1 0 0 0],F);
            r = convEnvCPLQ(q);
            testCase.verifyClass(r, 'RatPol');
            S = [0.25 0.25; 0.1 0.2; 0.4 0.4];   % interior points
            testCase.verifyEqual(r.eval(S), q.eval(S), 'AbsTol', 1e-12);
            testCase.verifyEqual(r.eval([0.25 0.25]), 0.0625, 'AbsTol', 1e-12);
        end

        function convexQuadraticFullDomainIsItself(testCase)
            q = QuaPoly.energy();                % 1/2(x^2+y^2) on all of R^2
            r = convEnvCPLQ(q);
            testCase.verifyClass(r, 'RatPol');
            S = [0 0; 1 1; -2 3];
            testCase.verifyEqual(r.eval(S), 0.5*(S(:,1).^2+S(:,2).^2), 'AbsTol', 1e-12);
        end

        function convexQuadraticOverSquareIsItself(testCase)
            % Convex case is valid over any convex face, not just triangles.
            V = [0 0; 1 0; 1 1; 0 1]; E = [1 2 1; 2 3 1; 3 4 1; 4 1 1]; F = [1 0;1 0;1 0;1 0];
            q = QuaPoly(V,E,[2 0 2 0 0 0],F);    % x^2 + y^2
            r = convEnvCPLQ(q);
            testCase.verifyEqual(r.eval([0.5 0.5]), 0.5, 'AbsTol', 1e-12);  % 0.25+0.25
        end

        function concaveOverTriangleIsAffineInterpolation(testCase)
            % q = -1/2(x^2+y^2) over triangle (0,0),(2,0),(0,2). Vertex values 0,-2,-2 =>
            % affine envelope z = -x - y on the triangle.
            V = [0 0; 2 0; 0 2]; E = [1 2 1; 2 3 1; 3 1 1]; F = [1 0; 1 0; 1 0];
            q = QuaPoly(V,E,[-1 0 -1 0 0 0],F);
            r = convEnvCPLQ(q);
            testCase.verifyClass(r, 'RatPol');
            % matches q at the vertices
            testCase.verifyEqual(r.eval(V), [0; -2; -2], 'AbsTol', 1e-12);
            % interior: affine value -x-y, and underestimates q
            testCase.verifyEqual(r.eval([0.5 0.5]), -1, 'AbsTol', 1e-12);
            testCase.verifyLessThanOrEqual(r.eval([0.5 0.5]) - q.eval([0.5 0.5]), 1e-12);
        end

        function bilinearZeroConvexEdgesIsAffine(testCase)
            % xy over the triangle (0,0),(1,0),(0,1): no convex edges (slopes 0, -1, vertical),
            % vertex values xy = 0,0,0 => convex envelope is the zero function.
            V = [0 0; 1 0; 0 1]; E = [1 2 1; 2 3 1; 3 1 1]; F = [1 0; 1 0; 1 0];
            q = convEnvCPLQ(QuaPoly(V,E,[0 1 0 0 0 0],F));
            testCase.verifyClass(q, 'RatPol');
            testCase.verifyEqual(q.eval([0.25 0.25; 0.1 0.6; 0.5 0.2]), [0;0;0], 'AbsTol', 1e-12);
        end

        function bilinearOneConvexEdgeRational(testCase)
            % [COAP] Appendix A.3.3 Example 2: conv(xy) over conv{(1,1),(0,0),(2,0)} = 2y^2/(y-x+2).
            V = [1 1; 0 0; 2 0]; E = [1 2 1; 2 3 1; 3 1 1]; F = [1 0; 1 0; 1 0];
            r = convEnvCPLQ(QuaPoly(V,E,[0 1 0 0 0 0],F));
            testCase.verifyClass(r, 'RatPol');
            S = [1 0.3; 1.2 0.3; 1.5 0.2; 0.6 0.3];               % interior points
            expected = 2*S(:,2).^2 ./ (S(:,2) - S(:,1) + 2);
            testCase.verifyEqual(r.eval(S), expected, 'AbsTol', 1e-12);
            % envelope touches f = xy along the convex edge (0,0)-(1,1)
            testCase.verifyEqual(r.eval([0.5 0.5]), 0.25, 'AbsTol', 1e-12);
        end

        function bilinearOneConvexEdgePlusLinear(testCase)
            % conv(xy + 2x - y) over the same triangle = 2y^2/(y-x+2) + 2x - y.
            V = [1 1; 0 0; 2 0]; E = [1 2 1; 2 3 1; 3 1 1]; F = [1 0; 1 0; 1 0];
            r = convEnvCPLQ(QuaPoly(V,E,[0 1 0 2 -1 0],F));
            S = [1 0.3; 1.2 0.3; 0.6 0.3];
            expected = 2*S(:,2).^2 ./ (S(:,2) - S(:,1) + 2) + 2*S(:,1) - S(:,2);
            testCase.verifyEqual(r.eval(S), expected, 'AbsTol', 1e-12);
        end

        function bilinearTwoConvexEdgesQuadratic(testCase)
            % [COAP] Appendix A.4.3's single-quadratic formula, conv(xy) over triangle (2,1),(0,0),
            % (1,0) = (x^2 + 2 sqrt(2) xy + 2 y^2 - x + 2 y) / (3 + 2 sqrt(2)) -- the paper states
            % "the domain is the entire triangle", but this is FALSE in general (confirmed against
            % ground truth; see bilinearTwoConvexEdgesSplitIsTight and DESIGN.md's Part 2 diagnosis
            % / fix): the formula is only tight in the sub-region containing the shared vertex
            % (2,1), bounded by the line where the formula's own implied dual point reaches the far
            % endpoint of a classified convex edge. (0.5,0.2) used to be asserted here too, but it
            % lies on the OTHER side of that boundary, where convEnvCPLQ now correctly returns a
            % different (rational) piece instead -- see bilinearTwoConvexEdgesSplitIsTight, which
            % checks tightness comprehensively via ground truth rather than this one closed form.
            V = [2 1; 0 0; 1 0]; E = [1 2 1; 2 3 1; 3 1 1]; F = [1 0; 1 0; 1 0];
            r = convEnvCPLQ(QuaPoly(V,E,[0 1 0 0 0 0],F));
            S = [1 0.3; 1.5 0.6; 1.8 0.85];                        % interior points, P-side only
            x = S(:,1); y = S(:,2);
            expected = (x.^2 + 2*sqrt(2)*x.*y + 2*y.^2 - x + 2*y) / (3 + 2*sqrt(2));
            testCase.verifyEqual(r.eval(S), expected, 'AbsTol', 1e-12);
        end

        function bilinearTwoConvexEdgesSplitIsTight(testCase)
            % Regression test for the tightness bug diagnosed in DESIGN.md / the session handoff:
            % the single-quadratic Appendix A.4 formula checked by bilinearTwoConvexEdgesQuadratic
            % above is a valid minorant but NOT always the tightest envelope over the WHOLE
            % triangle -- it is only tight on the sub-triangle containing the two convex edges'
            % shared vertex. Reproduced on the paper's OWN Appendix A.4.3 example, V=(2,1),(0,0),
            % (1,0): the weak (non-convex) edge is (0,0)-(1,0) (f=xy is identically 0 there), and
            % the old single-quadratic formula dipped to q1(0.474343,0)=-0.042780 instead of the
            % true envelope value ~0 (matching f exactly, since f=0 along the whole weak edge here
            % and the true envelope must equal the affine chord -- here 0 -- along it).
            %
            % Tightness is checked via each sub-triangle's OWN conjugate (Step 2, conjPieceCPLQ)
            % rather than the full q.conj('cplq') pipeline: this split's internal seam is a genuine
            % kink (unlike splitThreeConvex's C1-smooth seam for the 3-convex-edge case), which
            % Step 3 (maxQuaPar) does not yet support combining -- see DESIGN.md / session handoff
            % for that separate, newly-found gap. max(g1,g2) at a dual point s is exactly f*(s)
            % whenever the underlying 2-piece split is the TRUE envelope, regardless of Step 3.
            V = [2 1; 0 0; 1 0]; E = [1 2 1; 2 3 1; 3 1 1]; F = [1 0; 1 0; 1 0];
            q = QuaPoly(V, E, [0 1 0 0 0 0], F);
            r = convEnvCPLQ(q);
            testCase.verifyClass(r, 'RatPol');
            testCase.verifyEqual(r.nf, 2);           % genuine split (unlike the mirror-symmetric,
                                                      % no-split-needed case exercised elsewhere)

            x0 = [0.474343 0];                        % the weak-edge "dip" point from the diagnosis
            testCase.verifyEqual(r.eval(x0), 0, 'AbsTol', 1e-6);

            % underestimates f on interior sample points
            S = [1 0.3; 0.5 0.2; 1.5 0.6; 0.9 0.05; 1.9 0.95];
            testCase.verifyLessThanOrEqual(r.eval(S) - S(:,1).*S(:,2), 1e-9);

            g1 = conjPieceCPLQ(convEnvCPLQTest.extractTriFace(r, 1));
            g2 = conjPieceCPLQ(convEnvCPLQTest.extractTriFace(r, 2));
            sBad = [-0.008727 -0.999962];             % the paper's own flagged bad dual point:
            expected = convEnvCPLQTest.supBilinearOverPoly(sBad, V);   % old code gave 0.03864091
            testCase.verifyEqual(max(g1.eval(sBad), g2.eval(sBad)), expected, 'AbsTol', 1e-6);
            S2 = [1.90 2.50; -1 -1; 0.5 0.5; 3 -2; -3 2];
            for i = 1:size(S2,1)
                s = S2(i,:);
                expected = convEnvCPLQTest.supBilinearOverPoly(s, V);
                testCase.verifyEqual(max(g1.eval(s), g2.eval(s)), expected, ...
                    'AbsTol', 1e-8, sprintf('s=(%.4f,%.4f)', s(1), s(2)));
            end
        end

        function generalIndefiniteViaRotation(testCase)
            % conv(x^2 - y^2) over the triangle {(0,0),(1,0),(1,1)} = (x-y)^2/(1-y).
            % (It rotates to conv(xy) over conv{(1,1),(0,0),(2,0)} = 2y^2/(y-x+2).)
            % Vertices listed counter-clockwise so F=[1 0;...] (face on the left) is consistent.
            V = [0 0; 1 0; 1 1]; E = [1 2 1; 2 3 1; 3 1 1]; F = [1 0; 1 0; 1 0];
            r = convEnvCPLQ(QuaPoly(V,E,[2 0 -2 0 0 0],F));   % x^2 - y^2
            testCase.verifyClass(r, 'RatPol');
            S = [0.7 0.3; 0.8 0.2; 0.9 0.5];                  % interior, away from y=1
            x = S(:,1); y = S(:,2);
            testCase.verifyEqual(r.eval(S), (x-y).^2 ./ (1-y), 'AbsTol', 1e-12);
        end

        function negativeBilinearViaRotation(testCase)
            % conv(-xy) over the triangle {(0,0),(1,-1),(2,0)} = 2y^2/(2-x-y).
            % Vertices listed counter-clockwise so F=[1 0;...] (face on the left) is consistent.
            V = [0 0; 1 -1; 2 0]; E = [1 2 1; 2 3 1; 3 1 1]; F = [1 0; 1 0; 1 0];
            r = convEnvCPLQ(QuaPoly(V,E,[0 -1 0 0 0 0],F));   % -xy
            S = [1 -0.3; 1.2 -0.3; 0.6 -0.3];                 % interior
            x = S(:,1); y = S(:,2);
            testCase.verifyEqual(r.eval(S), 2*y.^2 ./ (2 - x - y), 'AbsTol', 1e-12);
        end

        function threeConvexEdgesSplit(testCase)
            % triangle (0,0),(1,1),(3,2): all three edge slopes (1, 1/2, 2/3) positive => split
            % into two sub-triangles (a 2-face RatPol). No closed form; check key properties.
            V = [0 0; 1 1; 3 2]; E = [1 2 1; 2 3 1; 3 1 1]; F = [1 0; 1 0; 1 0];
            r = convEnvCPLQ(QuaPoly(V,E,[0 1 0 0 0 0],F));    % xy
            testCase.verifyClass(r, 'RatPol');
            testCase.verifyEqual(r.nf, 2);
            % envelope equals xy at the three original vertices
            testCase.verifyEqual(r.eval([0 0; 1 1; 3 2]), [0; 1; 6], 'AbsTol', 1e-10);
            % finite and underestimates xy at interior points
            S = [1 0.8; 1.5 1.1; 2 1.4];
            vals = r.eval(S);
            testCase.verifyTrue(all(isfinite(vals)));
            testCase.verifyLessThanOrEqual(vals - S(:,1).*S(:,2), 1e-9);
        end

        function multiFaceConvexSquare(testCase)
            % Unit square split by the diagonal (0,0)-(1,1) into two triangles, convex
            % 1/2(x^2+y^2) on both: per-triangle envelope = the function, assembled over the square.
            V = [0 0; 1 0; 1 1; 0 1];
            E = [1 2 1; 2 3 1; 1 3 1; 3 4 1; 4 1 1];   % bottom, right, diagonal, top, left
            F = [1 0; 1 0; 2 1; 2 0; 2 0];             % face1 = lower-right, face2 = upper-left
            f = [1 0 1 0 0 0; 1 0 1 0 0 0];
            q = QuaPoly(V,E,f,F);
            % validate the input itself
            testCase.verifyEqual(q.eval([0.5 0.3; 0.3 0.5]), 0.5*[0.34; 0.34], 'AbsTol', 1e-12);
            r = convEnvCPLQ(q);
            testCase.verifyClass(r, 'RatPol');
            S = [0.5 0.3; 0.3 0.5; 0.5 0.5];
            testCase.verifyEqual(r.eval(S), 0.5*(S(:,1).^2 + S(:,2).^2), 'AbsTol', 1e-12);
        end

        function multiFaceBilinearSquare(testCase)
            % Same 2-triangle square, f = xy on both faces. Per-triangle convex envelopes:
            %   lower-right triangle (y<x): y^2/(y-x+1);  upper-left (y>x): x^2/(x-y+1).
            V = [0 0; 1 0; 1 1; 0 1];
            E = [1 2 1; 2 3 1; 1 3 1; 3 4 1; 4 1 1];
            F = [1 0; 1 0; 2 1; 2 0; 2 0];
            f = [0 1 0 0 0 0; 0 1 0 0 0 0];
            q = QuaPoly(V,E,f,F);
            testCase.verifyEqual(q.eval([0.5 0.3; 0.3 0.5]), [0.15; 0.15], 'AbsTol', 1e-12); % xy
            r = convEnvCPLQ(q);
            testCase.verifyEqual(r.nf, 2);
            testCase.verifyEqual(r.eval([0.5 0.3]), 0.3^2/(0.3-0.5+1), 'AbsTol', 1e-12); % T1: 0.1125
            testCase.verifyEqual(r.eval([0.3 0.5]), 0.3^2/(0.3-0.5+1), 'AbsTol', 1e-12); % T2: 0.1125
            testCase.verifyEqual(r.eval([0.5 0.5]), 0.25, 'AbsTol', 1e-12);               % touches xy
            % underestimates xy
            S = [0.6 0.2; 0.2 0.6];
            testCase.verifyLessThanOrEqual(r.eval(S) - S(:,1).*S(:,2), 1e-12);
        end

        function concaveOverSquareNotImplemented(testCase)
            % Non-convex over a single non-triangular face (nf==1) is not implemented yet.
            V = [0 0; 1 0; 1 1; 0 1]; E = [1 2 1; 2 3 1; 3 4 1; 4 1 1]; F = [1 0;1 0;1 0;1 0];
            q = QuaPoly(V,E,[-1 0 -1 0 0 0],F);
            testCase.verifyError(@() convEnvCPLQ(q), 'convEnvCPLQ:notImplemented');
        end

        function unboundedMultiFaceRejected(testCase)
            % oneNorm has 4 unbounded (ray) faces -> multi-face path requires bounded faces.
            testCase.verifyError(@() convEnvCPLQ(QuaPoly.oneNorm()), 'convEnvCPLQ:notImplemented');
        end

        function cubicRejected(testCase)
            q = QuaPoly([1 0 0 0 0 0 0 0 0 0]);  % cubic
            testCase.verifyError(@() convEnvCPLQ(q), 'PLQ:op:unsupportedType');
        end
    end

    methods (Static)
        function p = extractTriFace(r, k)
            % r: a 2-face RatPol (as produced by splitTwoConvexEdges). Returns face k's 3 vertices
            % (CCW) and quadratic coefficients as a standalone QuaPoly, the shape conjPieceCPLQ
            % requires. Same construction as maxQuaParTest.extractTriFace, duplicated here since
            % it is file-local there.
            edgeIdx = find(r.F(:,1)==k | r.F(:,2)==k);
            vids = unique(r.E(edgeIdx,1:2));
            Vt = r.V(vids,:);
            area2 = (Vt(2,1)-Vt(1,1))*(Vt(3,2)-Vt(1,2)) - (Vt(2,2)-Vt(1,2))*(Vt(3,1)-Vt(1,1));
            if area2 < 0, Vt = Vt([1 3 2],:); end
            E = [1 2 1; 2 3 1; 3 1 1]; F = [1 0; 1 0; 1 0];
            p = QuaPoly(Vt, E, r.f(k,5:10), F);
        end

        function h = supBilinearOverPoly(s, T)
            % Exact sup_{(x,y) in T} [s1 x + s2 y - x y]: the Hessian of the objective is
            % indefinite (eigenvalues +-1), so no interior point can be a local max, and the sup
            % is attained on T's boundary -- checked in closed form (quadratic-in-t along each
            % edge). Same construction as maxQuaParTest.supBilinearOverPoly, duplicated here since
            % it is file-local there.
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
