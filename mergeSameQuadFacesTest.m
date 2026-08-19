classdef mergeSameQuadFacesTest < matlab.unittest.TestCase
    % Step 0 of `conj`/`biconj`: an interior edge whose two sides carry the same quadratic is a
    % line the CALLER drew, not a property of f, and every dispatch downstream reads the mesh.
    %
    % The case this exists for is the unit square handed in as two triangles sharing its diagonal:
    % the same function as the one-face square, which answers by McCormick in 0 s, while the
    % two-face spelling used to go triangulate -> per-piece conjugate -> Step 3 -> second
    % conjugation and come out WRONG (biconjugateTest.biconjugateOverATwoFaceSubdivisionIsTheEnvelope
    % was a known-failing test for exactly this). The answer must not depend on how the domain was
    % drawn.

    methods (Static)
        function p = twoTriangleSquare(f)
        % The unit square as two triangles sharing the diagonal (1,1)-(0,0), f per face.
            V = [0 0; 1 0; 1 1; 0 1];
            E = [1 2 1; 2 3 1; 1 3 1; 3 4 1; 4 1 1];
            F = [1 0; 1 0; 2 1; 2 0; 2 0];
            p = QuaPol(V, E, f, F);
        end
    end

    methods (Test)
        function twoTrianglesSharingADiagonalBecomeOneSquare(testCase)
            p = mergeSameQuadFacesTest.twoTriangleSquare([0 1 0 0 0 0; 0 1 0 0 0 0]);
            [q, n] = mergeSameQuadFaces(p);

            testCase.verifyEqual(n, 1, 'one face must disappear');
            testCase.verifyEqual([q.nv q.ne q.nf], [4 4 1], ...
                'the merged square is four corners, four edges, one face');
            testCase.verifyEqual(q.f(1,:), p.f(1,:), 'the quadratic is carried over unchanged');

            % SHAPE IS THE POINT: mccormickEnvelope tests nv == 4, ne == 4, nf == 1 before it
            % will do anything, so anything else here leaves biconj on the long route.
            X = [0.75 0.25; 0.5 0.5; 0.1 0.9; 0 0; 1 1; 0.3 0.3];
            for i = 1:size(X,1)
                testCase.verifyEqual(q.eval(X(i,:)), p.eval(X(i,:)), 'AbsTol', 1e-12, ...
                    sprintf('the function changed at (%g,%g)', X(i,1), X(i,2)));
            end
            testCase.verifyEqual(q.eval([2 2]), inf, 'the domain is unchanged too');
        end

        function differentQuadraticsAreLeftAlone(testCase)
            p = mergeSameQuadFacesTest.twoTriangleSquare([0 1 0 0 0 0; 0 2 0 0 0 0]);
            [q, n] = mergeSameQuadFaces(p);
            testCase.verifyEqual(n, 0);
            testCase.verifyEqual([q.nv q.ne q.nf], [p.nv p.ne p.nf]);
            testCase.verifyEqual(q.E, p.E);
            testCase.verifyEqual(q.F, p.F);
        end

        function coefficientsThatMerelyAgreeCloselyAreTwoFaces(testCase)
        % NO TOLERANCE, deliberately: a normaliser may not change the function, and merging two
        % faces that differ in the twelfth digit does. Callers wanting fuzzy behaviour round first.
            p = mergeSameQuadFacesTest.twoTriangleSquare([0 1 0 0 0 0; 0 1+1e-12 0 0 0 0]);
            [~, n] = mergeSameQuadFaces(p);
            testCase.verifyEqual(n, 0);
        end

        function everyEdgeInteriorGivesTheWholePlane(testCase)
        % Two half-planes carrying the same quadratic, split by a line drawn as two opposite rays.
        % Nothing bounds anything, so the normal form is the one-face full-domain QuaPol.
            V = [0 0; 1 0; -1 0];
            E = [1 2 0; 1 3 0];
            F = [1 2; 2 1];
            p = QuaPol(V, E, [0 1 0 0 0 0; 0 1 0 0 0 0], F);
            [q, n] = mergeSameQuadFaces(p);
            testCase.verifyEqual(n, 1);
            testCase.verifyEqual([q.nv q.ne q.nf], [0 0 1]);
            testCase.verifyEqual(q.eval([2 3]), p.eval([2 3]), 'AbsTol', 1e-12);
        end

        function aReflexUnionIsRefused(testCase)
        % MEASURED, not assumed. Three wedges round the origin, the outer two carrying the same
        % affine piece: their union is a 225-degree REFLEX wedge, and this representation cannot
        % express one -- the merged mesh builds without complaint and then eval returns +inf at
        % (2,3), a point inside the reflex face, because a face is read as an intersection of
        % half-planes. So the merge is refused and the edge kept.
            V = [0 0; -1 0; 1 1; 0 -1];
            E = [1 2 0; 1 3 0; 1 4 0];
            F = [1 2; 2 3; 3 1];
            f = [0 0 0 0 1 0; 0 0 0 0 1 0; 0 0 0 1 0 0];   % y, y, x
            p = QuaPol(V, E, f, F);
            [q, n] = mergeSameQuadFaces(p);
            testCase.verifyEqual(n, 0, 'a reflex union must not be merged');
            testCase.verifyEqual([q.ne q.nf], [3 3]);
            for X = [2 3; -2 3; 0.5 -2; -1 -1]'
                testCase.verifyEqual(q.eval(X'), p.eval(X'), 'AbsTol', 1e-12, ...
                    sprintf('the function changed at (%g,%g)', X(1), X(2)));
            end
        end

        function aConvexUnionOfWedgesStillMerges(testCase)
        % The other side of the guard, and note how narrow it is for a FAN: three wedges around
        % one vertex sum to 360 degrees, so a pair is convex only when the third wedge is itself
        % at least 180. Rays east, north and west give two quadrants and a lower half-plane; the
        % two quadrants carry the same piece and fuse into the upper half-plane.
            V = [0 0; 1 0; 0 1; -1 0];
            E = [1 2 0; 1 3 0; 1 4 0];
            F = [1 3; 2 1; 3 2];
            f = [0 0 0 0 1 0; 0 0 0 0 1 0; 0 0 0 1 0 0];   % quadrants carry y, lower half x
            p = QuaPol(V, E, f, F);
            [q, n] = mergeSameQuadFaces(p);
            testCase.verifyEqual(n, 1);
            testCase.verifyEqual([q.ne q.nf], [2 2]);
            for X = [2 3; -2 3; 0.5 -2; -1 -1]'
                testCase.verifyEqual(q.eval(X'), p.eval(X'), 'AbsTol', 1e-12, ...
                    sprintf('the function changed at (%g,%g)', X(1), X(2)));
            end
        end

        function theOptOutSuppressesIt(testCase)
        % CCA2_NO_STEP0 exists for biconjugateTest's lens regression, which must keep exercising
        % the two-face route; if the flag ever stopped working that test would silently stop
        % testing anything.
            testCase.applyFixture(matlab.unittest.fixtures.EnvironmentVariableFixture( ...
                'CCA2_NO_STEP0', '1'));
            p = mergeSameQuadFacesTest.twoTriangleSquare([0 1 0 0 0 0; 0 1 0 0 0 0]);
            [q, n] = mergeSameQuadFaces(p);
            testCase.verifyEqual(n, 0);
            testCase.verifyEqual(q.nf, 2);
        end

        function biconjOfTheTwoFaceSquareIsMcCormickAndImmediate(testCase)
        % END TO END, and it is the reason Step 0 is in the dispatch rather than in a utility
        % nobody calls: with the mesh normalised, biconj of the two-face square takes the
        % McCormick short-circuit and returns a MESHED QuaPol -- max(0, x+y-1) over the square.
            p = mergeSameQuadFacesTest.twoTriangleSquare([0 1 0 0 0 0; 0 1 0 0 0 0]);
            h = p.biconj('cplq');
            testCase.verifyEqual(kind(h), 'QuaPol', ...
                'the long route returns a mesh-less QuaParCPLQ; the short-circuit returns a mesh');
            X = [0.75 0.25; 0.5 0.5; 0.25 0.25; 0.9 0.6; 0.2 0.8; 0.1 0.1; 0.6 0.6];
            for i = 1:size(X,1)
                testCase.verifyEqual(h.eval(X(i,:)), max(0, X(i,1) + X(i,2) - 1), ...
                    'AbsTol', 1e-9, sprintf('McCormick at (%g,%g)', X(i,1), X(i,2)));
            end
        end
    end
end
