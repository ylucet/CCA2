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
            % [COAP] Appendix A.4.3 Example: conv(xy) over triangle (2,1),(0,0),(1,0) =
            %   (x^2 + 2 sqrt(2) xy + 2 y^2 - x + 2 y) / (3 + 2 sqrt(2)).
            V = [2 1; 0 0; 1 0]; E = [1 2 1; 2 3 1; 3 1 1]; F = [1 0; 1 0; 1 0];
            r = convEnvCPLQ(QuaPoly(V,E,[0 1 0 0 0 0],F));
            S = [1 0.3; 0.5 0.2; 1.5 0.6];                         % interior points
            x = S(:,1); y = S(:,2);
            expected = (x.^2 + 2*sqrt(2)*x.*y + 2*y.^2 - x + 2*y) / (3 + 2*sqrt(2));
            testCase.verifyEqual(r.eval(S), expected, 'AbsTol', 1e-12);
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
end
