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

        function indefiniteOverTriangleNotImplemented(testCase)
            V = [0 0; 1 0; 0 1]; E = [1 2 1; 2 3 1; 3 1 1]; F = [1 0; 1 0; 1 0];
            q = QuaPoly(V,E,[0 1 0 0 0 0],F);    % xy (indefinite)
            testCase.verifyError(@() convEnvCPLQ(q), 'convEnvCPLQ:notImplemented');
        end

        function concaveOverSquareNotImplemented(testCase)
            % Non-convex over a non-triangle is not implemented yet.
            V = [0 0; 1 0; 1 1; 0 1]; E = [1 2 1; 2 3 1; 3 4 1; 4 1 1]; F = [1 0;1 0;1 0;1 0];
            q = QuaPoly(V,E,[-1 0 -1 0 0 0],F);
            testCase.verifyError(@() convEnvCPLQ(q), 'convEnvCPLQ:notImplemented');
        end

        function multiPieceNotImplemented(testCase)
            testCase.verifyError(@() convEnvCPLQ(QuaPoly.oneNorm()), 'convEnvCPLQ:notImplemented');
        end

        function cubicRejected(testCase)
            q = QuaPoly([1 0 0 0 0 0 0 0 0 0]);  % cubic
            testCase.verifyError(@() convEnvCPLQ(q), 'PLQ:op:unsupportedType');
        end
    end
end
