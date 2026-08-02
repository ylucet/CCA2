classdef biconjugateTest < matlab.unittest.TestCase
    % BICONJUGATE tests: f** = conv f, checked against closed forms.
    %
    % WHY THIS FILE EXISTS. Before it, NOTHING checked a biconjugate value anywhere in the
    % toolbox. testcPLQ/testRectBiconj runs the pipeline and asserts nothing at all;
    % conjCPLQTest.biconjCoverageByInputCase asserted only the RESULT CLASS, which passes on an
    % EMPTY piece list because QuaParCPLQ(functionNDomain.empty()).kind() is still 'QuaParCPLQ'.
    % So `f** = +inf everywhere` -- the actual behaviour for every input tried -- looked green.
    %
    % The algorithm ([JOGO]/[COAP]) is, in both passes:
    %     triangulate -> convex envelope per triangle (splitting as needed)
    %                 -> conjugate per piece -> MAX of those conjugates
    % applied to f to get f*, then to f* to get f**. The max is one shared implementation,
    % functionNDomain.maxOfList, used by plq.maximumConjugate, plq.biconjugateF and
    % QuaParCPLQ.conj alike -- it must never be duplicated per caller.
    %
    % Every case below has a closed form that owes nothing to the pipeline.

    methods (Static)
        function v = evalAny(g, p)
            if isa(g, 'QuaParCPLQ'), v = evalFunctionNDomain(g.fnd, p); else, v = g.eval(p); end
            if isnan(v), v = inf; end
        end
        function q = box(V)
            % QuaPol carrying f = 0 on the polygon V (i.e. the indicator of V), CCW.
            n = size(V,1);
            E = [(1:n)', [2:n 1]', ones(n,1)];
            F = [ones(n,1), zeros(n,1)];
            if biconjugateTest.signedArea(V) < 0, V = flipud(V); end
            q = QuaPol(V, E, zeros(1,6), F);
        end
        function a = signedArea(V)
            n = size(V,1); a = 0;
            for i = 1:n
                j = mod(i,n)+1;
                a = a + V(i,1)*V(j,2) - V(j,1)*V(i,2);
            end
            a = a/2;
        end
    end

    methods (Test)

        function oneNormBiconjugateIsItself(testCase)
        % THE elementary case, and convexdb's norm1_2d. f = ||x||_1 is convex, so f** = f, and
        % f* = I_{||s||_inf <= 1}. Both fixtures already shipped in QuaPol (oneNorm,
        % oneNormConjugate) and were used by a dozen eval/structure tests but NEVER for conj.
        % Note f's domain is four UNBOUNDED cones, so this also covers the unbounded case.
            p = QuaPol.oneNorm();
            testCase.verifyFalse(p.isDomBounded);

            g = p.conj('cplq');
            for s = [0 0; 0.5 0.5; 1 1; -1 -1; 0.9 -0.3]'
                testCase.verifyEqual(biconjugateTest.evalAny(g, s'), 0, 'AbsTol', 1e-9, ...
                    sprintf('f* must be 0 inside the l-inf ball at (%g,%g)', s(1), s(2)));
            end
            for s = [1.5 0; 0 2; -1.5 0.5; 2 2]'
                testCase.verifyEqual(biconjugateTest.evalAny(g, s'), inf, ...
                    sprintf('f* must be +inf outside the l-inf ball at (%g,%g)', s(1), s(2)));
            end

            h = p.biconj('cplq');
            X = [0 0; 1 0; 0.5 0.5; -1 2; 2 -3; 0.25 0.75; -2 -2];
            for i = 1:size(X,1)
                testCase.verifyEqual(biconjugateTest.evalAny(h, X(i,:)), ...
                    abs(X(i,1)) + abs(X(i,2)), 'AbsTol', 1e-9, ...
                    sprintf('f** must equal |x|+|y| at (%g,%g)', X(i,1), X(i,2)));
            end
        end

        function indicatorOfUnitSquareBiconjugateIsItself(testCase)
        % f = I_{[0,1]^2}: convex, so f** = f. Its conjugate is the support function of the
        % square, max(0,s1) + max(0,s2). A NON-TRIANGULAR bounded domain, and one whose f is
        % constant -- the shape whose conjugate cells are pure normal cones.
            p = biconjugateTest.box([0 0; 1 0; 1 1; 0 1]);

            g = p.conj('cplq');
            S = [0 0; 1 1; -1 -1; 2 -0.5; -0.5 2; 0.3 0.4];
            for i = 1:size(S,1)
                want = max(S(i,1),0) + max(S(i,2),0);
                testCase.verifyEqual(biconjugateTest.evalAny(g, S(i,:)), want, 'AbsTol', 1e-9, ...
                    sprintf('support function of [0,1]^2 at (%g,%g)', S(i,1), S(i,2)));
            end

            h = p.biconj('cplq');
            for x = [0.5 0.5; 0 0; 1 1; 0.25 0.75]'
                testCase.verifyEqual(biconjugateTest.evalAny(h, x'), 0, 'AbsTol', 1e-9, ...
                    sprintf('f** = 0 inside the square at (%g,%g)', x(1), x(2)));
            end
            for x = [1.5 0.5; -0.5 0.5; 2 2]'
                testCase.verifyEqual(biconjugateTest.evalAny(h, x'), inf, ...
                    sprintf('f** = +inf outside the square at (%g,%g)', x(1), x(2)));
            end
        end

        function bilinearOverABoxGivesTheMcCormickEnvelope(testCase)
        % THE case that motivated all of this: f = x*y on [0,1]^2 is NONCONVEX, so f** is a
        % genuine convex envelope, and for a bilinear term over a box it is exactly the McCormick
        % envelope max(0, x+y-1) (Al-Khayyal-Falk). A non-triangular domain whose envelope is NOT
        % obtainable from Step 1 alone -- Step 1 returns per-triangle envelopes, and the box's
        % true envelope only appears after conjugating twice.
            x = sym('x'); y = sym('y');
            d = domain([0 0; 0 1; 1 1; 1 0], x, y);
            P = plq([plq_1p(d, symbolicFunction(x*y))]);
            P = P.triangulate;
            P = P.maximum;
            P = P.biconjugateF;
            B = P.biconjugate;
            testCase.verifyNotEmpty(B, 'the biconjugate must not be empty');

            X = [0.75 0.25; 0.5 0.5; 0.25 0.25; 0.9 0.6; 0.2 0.8; 0.5 0.9; 0.1 0.1; ...
                 0.95 0.95; 0.4 0.3; 0.05 0.9; 0.6 0.6];
            for i = 1:size(X,1)
                got = evalFunctionNDomain(B, X(i,:));
                want = max(0, X(i,1) + X(i,2) - 1);
                testCase.verifyFalse(isnan(got), sprintf( ...
                    'f** uncovered at (%g,%g) -- inside the box', X(i,1), X(i,2)));
                testCase.verifyEqual(got, want, 'AbsTol', 1e-9, sprintf( ...
                    'McCormick envelope at (%g,%g)', X(i,1), X(i,2)));
            end
        end

        function biconjugateOfAConvexQuadraticOverATriangleIsItself(testCase)
        % A bounded TRIANGLE with a convex q: f is already convex, so f** = f on the triangle
        % and +inf off it. Guards the triangle path against the non-triangular ones above.
            V = [0 0; 1 0; 0 1]; E = [1 2 1; 2 3 1; 3 1 1]; F = [1 0; 1 0; 1 0];
            p = QuaPol(V, E, [1 0 1 0 0 0], F);           % q = (x^2+y^2)/2
            h = p.biconj('cplq');
            X = [0.2 0.2; 0.5 0.25; 1/3 1/3; 0 0];
            for i = 1:size(X,1)
                want = 0.5*(X(i,1)^2 + X(i,2)^2);
                testCase.verifyEqual(biconjugateTest.evalAny(h, X(i,:)), want, 'AbsTol', 1e-9, ...
                    sprintf('f** = q at (%g,%g)', X(i,1), X(i,2)));
            end
        end

        % KNOWN FAILING, left in deliberately -- see the note below. Remove the tag when fixed.
        function biconjugateOverATwoFaceSubdivisionIsTheEnvelope(testCase)
        % A genuine multi-FACE subdivision (nf = 2), not one face that triangulate happens to
        % split: the unit square given as two triangles sharing the diagonal, f = x*y on both.
        % f is nonconvex, so f** is a real envelope, and over the square it is again McCormick,
        % max(0, x+y-1) -- so the answer must not depend on whether the square arrives as one
        % face or two.
        %
        % A PENTAGON was tried here first and removed: correct-looking but far too slow for a
        % unit test. Measured -- triangulate 3 pieces in 1 s, then f* took 885 s and came out
        % with 41 REGIONS, after which the second pass has to conjugate and max all 41. The cost
        % is in the number of f* regions, which grows fast with the vertex count; that is worth
        % knowing before anyone reaches for a bigger polygon as a test case.
            V = [0 0; 1 0; 1 1; 0 1];
            E = [1 2 1; 2 3 1; 1 3 1; 3 4 1; 4 1 1];
            F = [1 0; 1 0; 2 1; 2 0; 2 0];
            p = QuaPol(V, E, [0 1 0 0 0 0; 0 1 0 0 0 0], F);
            % OPEN DEFECT this test pins. The SAME function on the SAME domain gives the right
            % answer as ONE face (bilinearOverABoxGivesTheMcCormickEnvelope) and the wrong one
            % as two. Both first conjugates are correct -- verified against a brute-force sup at
            % 7 dual points, one-face f* (6 regions) and two-face f* (9 regions) agree exactly.
            % The difference is that the two-face f* carries PARABOLIC pieces (e.g. the conic
            % 2*s1*s2 - 4*s1 + s1^2 + s2^2 <= 0), and on those the second pass returns the
            % PER-FACE envelope: at (0.5,0.5), (0.25,0.25), (0.1,0.1), (0.6,0.6) -- all on the
            % shared diagonal y = x, where x*y is convex and a diagonal split relaxes nothing --
            % it returns x*y itself instead of max(0,x+y-1).
            % Blocks 1 and 4 of that f* differ ONLY in the sign of that conic (opposite sides of
            % one parabola) yet produce IDENTICAL conjugate cells, so at least one is wrong.
            % That is the parabolic branch (getNormalConeVertexQ / getSubdiffVertexT2Q / the
            % isQuad chord rewrite), NOT the polyhedral one fixed alongside these tests.
            % Replacing the isQuad chord's hardcoded vx(1),vx(2) endpoints with the vertices the
            % conic actually touches was tried and changed nothing here -- do not re-try it
            % without new evidence.
            testCase.verifyEqual(p.nf, 2);
            testCase.verifyTrue(p.isDomBounded);

            h = p.biconj('cplq');
            X = [0.75 0.25; 0.5 0.5; 0.25 0.25; 0.9 0.6; 0.2 0.8; 0.1 0.1; 0.6 0.6];
            for i = 1:size(X,1)
                got = biconjugateTest.evalAny(h, X(i,:));
                want = max(0, X(i,1) + X(i,2) - 1);
                testCase.verifyEqual(got, want, 'AbsTol', 1e-9, sprintf( ...
                    'two-face square: McCormick at (%g,%g)', X(i,1), X(i,2)));
            end
        end

        function biconjugateOfAnUnboundedPiecewiseAffineIsItself(testCase)
        % UNBOUNDED, non-triangular, and multi-face: f = max(0, x, y) as three unbounded wedges.
        % f is convex, so f** = f. Its conjugate is the indicator of the simplex
        % {s >= 0, s1+s2 <= 1}, which is where the second pass has to conjugate an INDICATOR on a
        % bounded polygon -- the shape that used to come back empty.
            V = [0 0; -1 0; 1 1; 0 -1];
            E = [1 2 0; 1 3 0; 1 4 0];
            F = [1 2; 2 3; 3 1];
            f = [0 0 0 0 0 0; 0 0 0 0 1 0; 0 0 0 1 0 0];   % 0, y, x
            p = QuaPol(V, E, f, F);
            testCase.verifyFalse(p.isDomBounded);

            g = p.conj('cplq');
            for s = [0.2 0.2; 0 0; 1 0; 0 1; 0.5 0.5]'
                testCase.verifyEqual(biconjugateTest.evalAny(g, s'), 0, 'AbsTol', 1e-9, ...
                    sprintf('f* = 0 on the simplex at (%g,%g)', s(1), s(2)));
            end
            for s = [-0.2 0.3; 1 1; 0.6 0.6]'
                testCase.verifyEqual(biconjugateTest.evalAny(g, s'), inf, ...
                    sprintf('f* = +inf off the simplex at (%g,%g)', s(1), s(2)));
            end

            h = p.biconj('cplq');
            X = [0 0; 2 1; -1 3; -2 -3; 0.5 0.25];
            for i = 1:size(X,1)
                want = max([0, X(i,1), X(i,2)]);
                testCase.verifyEqual(biconjugateTest.evalAny(h, X(i,:)), want, 'AbsTol', 1e-9, ...
                    sprintf('f** = max(0,x,y) at (%g,%g)', X(i,1), X(i,2)));
            end
        end
    end
end
