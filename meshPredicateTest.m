classdef meshPredicateTest < matlab.unittest.TestCase
% Unit tests for the mesh PREDICATES the three sibling mesh classes each implement:
% `isFaceConvex`, `isEdgeContinuous`, `isEdgeConvex`, `isConvex`, `isDomBounded`,
% `belongToEdge`, `evalHessian`, `matrixForm`, `oneNormConjugate`.
%
% WHY ONE SUITE FOR THREE CLASSES. `QuaPol`, `QuaPar` and `RatPol` are siblings, not a chain, and
% each carries its own copy of these methods. Only QuaPol's were reached by any test -- through
% `PLQVCTest`, and then only in passing, as part of larger pipeline checks. Measured 2026-08-31:
% QuaPar 36.7% and RatPol 28.5% under fast+normal, with fourteen identical predicates dark in both.
%
% Running the SAME assertions against all three is what makes this a unit suite rather than three
% of them: the properties are defined on the mathematics, not on the representation, so the three
% implementations must agree, and a divergence between siblings is itself a defect worth failing
% on. All three constructors accept (V,E,f,F), which is what makes the sharing possible.
%
% BUCKET: fast. Everything here is closed-form numerics on 2-face meshes; nothing calls the
% Symbolic Math Toolbox.

    methods (Static)

        % ---- fixtures ---------------------------------------------------------------------------
        % One mesh shape throughout: the unit square split by the diagonal (0,0)-(1,1) into a lower
        % and an upper triangle. Two faces, one INTERIOR edge (the diagonal) and four boundary
        % edges, which is the smallest mesh on which "across an edge" means anything.
        function [V, E, F] = splitSquare()
            V = [0 0; 1 0; 1 1; 0 1];
            E = [1 2 1;      % boundary, lower face only
                 2 3 1;      % boundary, lower face only
                 3 4 1;      % boundary, upper face only
                 4 1 1;      % boundary, upper face only
                 1 3 1];     % the DIAGONAL: interior, both faces
            F = [1 0; 1 0; 2 0; 2 0; 1 2];
        end

        function p = build(cls, f)
        % The split square carrying the per-face coefficient rows f (2 x 6, [x^2 xy y^2 x y 1]).
            [V, E, F] = meshPredicateTest.splitSquare();
            switch cls
                case 'QuaPol', p = QuaPol(V, E, f, F);
                case 'QuaPar', p = QuaPar(V, E, f, F);
                case 'RatPol', p = RatPol(V, E, f, F);
                otherwise, error('meshPredicateTest:cls', 'unknown class %s', cls);
            end
        end

        function cs = classes()
            cs = {'QuaPol', 'QuaPar', 'RatPol'};
        end

        function j = diagonalEdge()
            j = 5;                                   % the interior edge of splitSquare
        end

        % ---- an independent evaluator, so nothing is asserted against the class under test ------
        function v = polyVal(c, x, y)
        % The stored 6-vector evaluated directly, from the WEIGHTED basis the classes actually use.
        %
        % THE WEIGHTS ARE THE CONVENTION AND THEY ARE NOT OBVIOUS. `evalPoly` multiplies the
        % monomials by [1/6 1/2 1/2 1/6 1/2 1 1/2 1 1 1] before contracting with the coefficient
        % vector, "to avoid constant coefficients in computation of Hessian". So a quadratic row
        % [a b c d e f] means
        %       a/2 x^2  +  b xy  +  c/2 y^2  +  d x  +  e y  +  f
        % -- the pure squares carry a half, the cross term does not. That is exactly what makes
        % Q = [a b; b c] the Hessian, since f = 1/2 z'Qz + L'z + const.
        %
        % Written out here rather than by calling evalPoly, so that the checks below test the
        % classes against the convention rather than against themselves.
            v = 0.5*c(1)*x.^2 + c(2)*x.*y + 0.5*c(3)*y.^2 + c(4)*x + c(5)*y + c(6);
        end
    end

    methods (Test)

        % =========================================================================================
        % matrixForm -- the (L, Q, C) split must reconstruct the polynomial it came from
        % =========================================================================================
        function matrixFormReconstructsTheQuadratic(testCase)
        % f(x) = 1/2 x'Qx + L'x + c is the documented meaning of the split, so the check is to
        % rebuild f from the pieces and compare to the coefficients as written -- at random points,
        % not at a golden value. Run for all three classes, which must agree.
            rng(20260831);
            for cls = meshPredicateTest.classes()
                for trial = 1:8
                    c = round((rand(1,6)-0.5)*8, 3);
                    [L, Q, C] = feval([cls{1} '.matrixForm'], c);
                    testCase.verifyLessThanOrEqual(sum(abs(C),'all'), sqrt(eps), sprintf( ...
                        '%s.matrixForm: a QUADRATIC must give a zero cubic tensor', cls{1}));
                    P = (rand(6,2)-0.5)*4;
                    for i = 1:size(P,1)
                        z = P(i,:).';
                        got  = 0.5*(z.'*Q*z) + L.'*z + c(6);
                        want = meshPredicateTest.polyVal(c, z(1), z(2));
                        testCase.verifyEqual(got, want, 'AbsTol', 1e-12, sprintf( ...
                            '%s.matrixForm trial %d at (%g,%g)', cls{1}, trial, z(1), z(2)));
                    end
                end
            end
        end

        % =========================================================================================
        % evalHessian -- against a central finite difference, which knows nothing of the storage
        % =========================================================================================
        function evalHessianMatchesFiniteDifferences(testCase)
            rng(20260831);
            h = 1e-4;
            for cls = meshPredicateTest.classes()
                for trial = 1:6
                    c = round((rand(1,6)-0.5)*6, 3);
                    X = (rand(3,2)-0.5)*3;
                    [H, g] = feval([cls{1} '.evalHessian'], c, X);
                    for i = 1:size(X,1)
                        x = X(i,1); y = X(i,2);
                        f  = @(a,b) meshPredicateTest.polyVal(c, a, b);
                        fdH = [ (f(x+h,y) - 2*f(x,y) + f(x-h,y))/h^2, ...
                                (f(x+h,y+h) - f(x+h,y-h) - f(x-h,y+h) + f(x-h,y-h))/(4*h^2); 0 0];
                        fdH(2,1) = fdH(1,2);
                        fdH(2,2) = (f(x,y+h) - 2*f(x,y) + f(x,y-h))/h^2;
                        fdg = [ (f(x+h,y) - f(x-h,y))/(2*h); (f(x,y+h) - f(x,y-h))/(2*h) ];
                        testCase.verifyEqual(H(:,:,i), fdH, 'AbsTol', 1e-4, sprintf( ...
                            '%s.evalHessian: Hessian at (%g,%g), trial %d', cls{1}, x, y, trial));
                        testCase.verifyEqual(g(:,i), fdg, 'AbsTol', 1e-6, sprintf( ...
                            '%s.evalHessian: gradient at (%g,%g), trial %d', cls{1}, x, y, trial));
                    end
                end
            end
        end

        function evalHessianOfAQuadraticNeedsNoPoint(testCase)
        % Documented shortcut: for a quadratic, H is constant, so the one-output call may omit x.
        % Asserted because it is a separate branch, and a wrong constant there is invisible to the
        % test above, which always passes points.
            c = [2 -1 3 5 -4 7];                     % x^2 - xy + 1.5y^2 + 5x - 4y + 7
            want = [2 -1; -1 3];                     % = Q, which IS the Hessian in this basis
            for cls = meshPredicateTest.classes()
                H = feval([cls{1} '.evalHessian'], c);
                testCase.verifyEqual(H, want, 'AbsTol', 1e-12, sprintf( ...
                    '%s.evalHessian(c) with no point', cls{1}));
            end
        end

        % =========================================================================================
        % isFaceConvex / isEdgeContinuous / isEdgeConvex / isConvex
        %
        % Four fixtures on the same mesh, chosen so that each predicate is the ONLY thing that
        % separates two of them -- otherwise a predicate that always returned the same answer
        % would pass.
        %
        %   A  x^2+y^2 on both faces      convex faces, continuous, convex overall
        %   B  0 / 1                      convex faces, DISCONTINUOUS
        %   C  -(x^2+y^2) on both         continuous, faces NOT convex
        %   D  0 / -x + y                 convex (affine) faces, continuous, but a CONCAVE kink
        %                                 across the diagonal: convex per face, not convex overall
        % =========================================================================================
        function convexFixtureIsConvexEverywhere(testCase)
            f = [1 0 1 0 0 0; 1 0 1 0 0 0];          % x^2 + y^2 on both faces
            j = meshPredicateTest.diagonalEdge();
            for cls = meshPredicateTest.classes()
                p = meshPredicateTest.build(cls{1}, f);
                testCase.verifyTrue(p.isFaceConvex(1), sprintf('%s: face 1 of x^2+y^2', cls{1}));
                testCase.verifyTrue(p.isFaceConvex(2), sprintf('%s: face 2 of x^2+y^2', cls{1}));
                testCase.verifyTrue(p.isEdgeContinuous(j), sprintf('%s: same function both sides', cls{1}));
                [cvx, cont] = p.isEdgeConvex(j);
                testCase.verifyTrue(cont, sprintf('%s: isEdgeConvex must also report continuity', cls{1}));
                testCase.verifyTrue(cvx, sprintf('%s: no kink, so the edge is convex', cls{1}));
                testCase.verifyTrue(p.isConvex(), sprintf('%s: x^2+y^2 is convex', cls{1}));
            end
        end

        function discontinuousFixtureIsCaughtAtTheEdge(testCase)
        % The two faces carry 0 and 1, so the function jumps by exactly 1 across the diagonal.
        % isEdgeContinuous must say so, isEdgeConvex must report the discontinuity in BOTH of its
        % outputs, and isConvex must be false -- a discontinuous function on a connected domain
        % cannot be convex.
            f = [0 0 0 0 0 0; 0 0 0 0 0 1];
            j = meshPredicateTest.diagonalEdge();
            for cls = meshPredicateTest.classes()
                p = meshPredicateTest.build(cls{1}, f);
                testCase.verifyFalse(p.isEdgeContinuous(j), sprintf( ...
                    '%s: a jump of 1 across the diagonal was not detected', cls{1}));
                [cvx, cont] = p.isEdgeConvex(j);
                testCase.verifyFalse(cont, sprintf('%s: isEdgeConvex continuity output', cls{1}));
                testCase.verifyFalse(cvx, sprintf( ...
                    '%s: a discontinuous edge cannot be convex', cls{1}));
                testCase.verifyFalse(p.isConvex(), sprintf('%s: overall convexity', cls{1}));
                % ...but each FACE on its own is a constant, hence convex. This is the pair that
                % stops isFaceConvex and isEdgeConvex being the same test.
                testCase.verifyTrue(p.isFaceConvex(1), sprintf('%s: a constant face is convex', cls{1}));
                testCase.verifyTrue(p.isFaceConvex(2), sprintf('%s: a constant face is convex', cls{1}));
            end
        end

        function concaveFacesAreNotFaceConvex(testCase)
            f = [-1 0 -1 0 0 0; -1 0 -1 0 0 0];      % -(x^2 + y^2)
            j = meshPredicateTest.diagonalEdge();
            for cls = meshPredicateTest.classes()
                p = meshPredicateTest.build(cls{1}, f);
                testCase.verifyFalse(p.isFaceConvex(1), sprintf( ...
                    '%s: -(x^2+y^2) is concave, so the face is not convex', cls{1}));
                testCase.verifyFalse(p.isConvex(), sprintf('%s: overall convexity', cls{1}));
                % Continuity is INDEPENDENT of convexity and must still be reported true here.
                testCase.verifyTrue(p.isEdgeContinuous(j), sprintf( ...
                    '%s: the same concave function on both sides is still continuous', cls{1}));
            end
        end

        function aConcaveKinkAcrossTheEdgeIsNotEdgeConvex(testCase)
        % The discriminating fixture. Lower face 0, upper face -x + y. On the diagonal y = x the
        % two agree, so the function is continuous; both faces are AFFINE, so both are convex; and
        % yet the function has a downward kink across the diagonal, so its subdifferential there is
        % empty and it is not convex. Only isEdgeConvex can tell -- which is what it is for.
            f = [0 0 0 0 0 0; 0 0 0 -1 1 0];
            j = meshPredicateTest.diagonalEdge();
            for cls = meshPredicateTest.classes()
                p = meshPredicateTest.build(cls{1}, f);
                testCase.verifyTrue(p.isFaceConvex(1), sprintf('%s: an affine face is convex', cls{1}));
                testCase.verifyTrue(p.isFaceConvex(2), sprintf('%s: an affine face is convex', cls{1}));
                testCase.verifyTrue(p.isEdgeContinuous(j), sprintf( ...
                    '%s: the two affine pieces agree on y = x', cls{1}));
                [cvx, cont] = p.isEdgeConvex(j);
                testCase.verifyTrue(cont, sprintf('%s: continuity output', cls{1}));
                testCase.verifyFalse(cvx, sprintf( ...
                    '%s: a concave kink across the edge must not read as convex', cls{1}));
                testCase.verifyFalse(p.isConvex(), sprintf( ...
                    '%s: convex faces plus a concave kink is not a convex function', cls{1}));
            end
        end

        function boundaryEdgesAreVacuouslyContinuous(testCase)
        % An edge with only one face beside it has nothing to be discontinuous across, and the
        % documented answer is true. Asserted because it is the early-return branch and a predicate
        % that got it backwards would make every bounded mesh look discontinuous.
            f = [1 0 1 0 0 0; 1 0 1 0 0 0];
            for cls = meshPredicateTest.classes()
                p = meshPredicateTest.build(cls{1}, f);
                for j = 1:4                                   % the four boundary edges
                    testCase.verifyTrue(p.isEdgeContinuous(j), sprintf( ...
                        '%s: boundary edge %d has one face and must read continuous', cls{1}, j));
                end
            end
        end

        % =========================================================================================
        % isDomBounded
        % =========================================================================================
        function isDomBoundedSeparatesSegmentsFromRays(testCase)
        % The documented rule is "bounded iff every edge is a segment".
        %
        % The unbounded fixture is a genuine WEDGE -- the first quadrant, one vertex at the origin
        % and two rays -- not the square with a segment flag flipped to 0. Flipping a flag makes a
        % mesh whose face boundary no longer closes, and the constructor rejects it
        % (QuaPol/orderEdges: "Face 1 has ..."), which tests the constructor's validation rather
        % than isDomBounded. A ray's second vertex is a DIRECTION MARKER, so an unbounded fixture
        % has to be built as one from the start.
            f  = [1 0 1 0 0 0; 1 0 1 0 0 0];
            fw = [1 0 1 0 0 0];
            Vw = [0 0; 1 0; 0 1];
            Ew = [1 2 0; 1 3 0];                              % two rays from the origin
            Fw = [1 0; 0 1];                                  % the quadrant is left of e1, right of e2
            for cls = meshPredicateTest.classes()
                pB = meshPredicateTest.build(cls{1}, f);
                testCase.verifyTrue(pB.isDomBounded(), sprintf( ...
                    '%s: every edge is a segment, so the domain is bounded', cls{1}));
                switch cls{1}
                    case 'QuaPol', pR = QuaPol(Vw, Ew, fw, Fw);
                    case 'QuaPar', pR = QuaPar(Vw, Ew, fw, Fw);
                    case 'RatPol', pR = RatPol(Vw, Ew, fw, Fw);
                end
                testCase.verifyFalse(pR.isDomBounded(), sprintf( ...
                    '%s: a wedge has two rays, so its domain is unbounded', cls{1}));
            end
        end

        % =========================================================================================
        % belongToEdge -- the segment/ray distinction is the whole content
        % =========================================================================================
        function belongToEdgeSeparatesSegmentFromRay(testCase)
        % [V1,V2] as a SEGMENT contains only 0 <= t <= 1; as a RAY it contains every t >= 0. The
        % same three probes therefore have to give different answers depending on the flag, and a
        % point off the line must be rejected by both.
            V1 = [0 0]; V2 = [1 0];
            X  = [0.5 0;      % inside the segment
                  2.0 0;      % beyond V2: on the ray, off the segment
                 -1.0 0;      % behind V1: off both
                  0.5 0.3];   % off the line entirely
            wantSeg = [true; false; false; false];
            wantRay = [true; true;  false; false];
            for cls = meshPredicateTest.classes()
                fn = [cls{1} '.belongToEdge'];
                if ~ismember('belongToEdge', methods(cls{1})), continue, end
                gotSeg = logical(feval(fn, V1, V2, X, true));
                testCase.verifyEqual(gotSeg(:), wantSeg, sprintf( ...
                    '%s.belongToEdge as a SEGMENT', cls{1}));
                gotRay = logical(feval(fn, V1, V2, X, false));
                testCase.verifyEqual(gotRay(:), wantRay, sprintf( ...
                    '%s.belongToEdge as a RAY', cls{1}));
            end
        end

        % =========================================================================================
        % oneNormConjugate -- a named fixture that is only worth keeping if it is the right one
        % =========================================================================================
        function oneNormConjugateIsTheIndicatorOfTheUnitInfinityBall(testCase)
        % (||.||_1)* is the indicator of the unit ball of the dual norm, ||.||_inf. So the mesh
        % must evaluate to 0 strictly inside [-1,1]^2 and to +inf outside it. Asserted rather than
        % assumed: the method is a hand-written vertex list, and nothing else in the toolbox has
        % ever checked that the list says what its name claims.
            inside  = [0 0; 0.5 0.5; -0.9 0.2; 0.2 -0.9];
            outside = [1.5 0; 0 -2; 3 3; -1.2 -1.2];
            for cls = {'QuaPol', 'QuaPar'}                    % RatPol has no such static
                if ~ismember('oneNormConjugate', methods(cls{1})), continue, end
                p = feval([cls{1} '.oneNormConjugate']);
                for i = 1:size(inside,1)
                    testCase.verifyEqual(p.eval(inside(i,:)), 0, 'AbsTol', 1e-12, sprintf( ...
                        '%s.oneNormConjugate at (%g,%g), inside the unit inf-ball', ...
                        cls{1}, inside(i,1), inside(i,2)));
                end
                for i = 1:size(outside,1)
                    testCase.verifyEqual(p.eval(outside(i,:)), Inf, sprintf( ...
                        '%s.oneNormConjugate at (%g,%g), outside the unit inf-ball', ...
                        cls{1}, outside(i,1), outside(i,2)));
                end
            end
        end

        % =========================================================================================
        % The siblings must agree with each other
        % =========================================================================================
        function theThreeClassesAgreeOnEveryPredicate(testCase)
        % Three independent implementations of the same mathematics. Any disagreement is a defect
        % in at least one of them, and this catches it without knowing which is right -- the kind
        % of check that only exists once the three are driven from one place.
            fixtures = { [1 0 1 0 0 0; 1 0 1 0 0 0], 'convex'
                         [0 0 0 0 0 0; 0 0 0 0 0 1], 'discontinuous'
                         [-1 0 -1 0 0 0; -1 0 -1 0 0 0], 'concave'
                         [0 0 0 0 0 0; 0 0 0 -1 1 0], 'kinked' };
            j = meshPredicateTest.diagonalEdge();
            for k = 1:size(fixtures,1)
                f = fixtures{k,1}; nm = fixtures{k,2};
                got = struct('faceConvex',[], 'edgeCont',[], 'edgeConvex',[], 'convex',[], 'bounded',[]);
                for cls = meshPredicateTest.classes()
                    p = meshPredicateTest.build(cls{1}, f);
                    [cvx, cont] = p.isEdgeConvex(j);
                    got.faceConvex(end+1) = p.isFaceConvex(1);
                    got.edgeCont(end+1)   = p.isEdgeContinuous(j);
                    got.edgeConvex(end+1) = cvx;
                    got.convex(end+1)     = p.isConvex();
                    got.bounded(end+1)    = p.isDomBounded();
                end
                for fld = fieldnames(got).'
                    v = got.(fld{1});
                    testCase.verifyEqual(numel(unique(v)), 1, sprintf( ...
                        ['fixture "%s": QuaPol/QuaPar/RatPol disagree on %s -- got [%s]. ' ...
                         'These are three copies of one predicate; a divergence is a defect ' ...
                         'in at least one of them.'], nm, fld{1}, num2str(v)));
                end
            end
        end
    end
end
