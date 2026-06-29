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

        function indefiniteRejected(testCase)
            V = [0 0; 1 0; 0 1]; E = [1 2 1; 2 3 1; 3 1 1]; F = [1 0; 1 0; 1 0];
            q = QuaPoly(V, E, [0 1 0 0 0 0], F);          % xy (not PD)
            testCase.verifyError(@() conjPieceCPLQ(q), 'conjPieceCPLQ:notImplemented');
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
