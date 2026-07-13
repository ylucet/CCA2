classdef addQuaPolyTest < matlab.unittest.TestCase
    % Tests for QuaPoly.add / addQuaPoly.m: pointwise sum of two QuaPoly functions.

    methods (Test)
        function bothFullDomainSumsCoefficients(testCase)
            f = QuaPoly.energy();   % 0.5*(x^2+y^2), full domain
            g = QuaPoly.energy();
            h = f.add(g);
            testCase.verifyEqual(h.nv, 0);
            testCase.verifyEqual(h.nf, 1);
            testCase.verifyEqual(h.eval([1 2]), 1^2 + 2^2, 'AbsTol', 1e-10);   % x^2+y^2 at (1,2)
        end

        function fullDomainPlusBoundedKeepsBoundedDomain(testCase)
            f = QuaPoly.energy();              % 0.5*(x^2+y^2), full domain
            g = addQuaPolyTest.squareLinear([0 0], [0 1 0]);   % g(x,y)=x on [0,2]x[0,2]
            h = f.add(g);
            testCase.verifyEqual(h.nv, g.nv);
            testCase.verifyEqual(h.nf, g.nf);
            testCase.verifyEqual(h.eval([1 1]), 0.5*(1+1) + 1, 'AbsTol', 1e-10);   % inside g's domain
            testCase.verifyEqual(h.eval([5 5]), Inf);                              % outside g's domain
        end

        function twoOverlappingSquaresSumOnTheirIntersection(testCase)
            % f(x,y)=x on [0,2]x[0,2]; g(x,y)=y on [1,3]x[1,3]. Overlap = [1,2]x[1,2],
            % h(x,y)=x+y there, +Inf everywhere else (independently computable ground truth).
            f = addQuaPolyTest.squareLinear([0 0], [0 1 0]);
            g = addQuaPolyTest.squareLinear([1 1], [0 1 0]);
            h = f.add(g);

            testCase.verifyClass(h, 'QuaPoly');
            testCase.verifyEqual(h.eval([1.5 1.5]), 1.5+1.5, 'AbsTol', 1e-10);   % inside the overlap
            testCase.verifyEqual(h.eval([0.5 0.5]), Inf);                        % inside f only
            testCase.verifyEqual(h.eval([2.5 2.5]), Inf);                        % inside g only
            testCase.verifyEqual(h.eval([-1 -1]), Inf);                          % inside neither
            testCase.verifyTrue(h.dom.isConvex);
        end

        function missingClosingEdgeConstraintIsEnforced(testCase)
            % Regression test for a bug found while implementing addQuaPar.m: polyConstraints
            % used to loop only 1:nv-1 for a bounded poly, silently dropping its closing edge
            % (nv,1) -- e.g. clipping f=[0,2]x[0,2] by g=[1,3]x[1,3] never enforced one of g's
            % four sides, so a point outside g along exactly that side (but inside f) was wrongly
            % treated as inside the overlap. twoOverlappingSquaresSumOnTheirIntersection's own
            % sample points didn't happen to probe that missing side; this one does.
            f = addQuaPolyTest.squareLinear([0 0], [0 1 0]);   % [0,2]x[0,2]
            g = addQuaPolyTest.squareLinear([1 1], [0 1 0]);   % [1,3]x[1,3]
            h = f.add(g);
            testCase.verifyEqual(h.eval([1.5 0.5]), Inf);   % inside f, outside g (y<1): must be Inf
            testCase.verifyEqual(h.eval([0.5 1.5]), Inf);   % inside f, outside g (x<1): must be Inf
        end

        function noOverlapErrors(testCase)
            f = addQuaPolyTest.squareLinear([0 0], [0 1 0]);
            g = addQuaPolyTest.squareLinear([10 10], [0 0 1]);
            testCase.verifyError(@() f.add(g), 'addQuaPoly:noOverlap');
        end

        function nestedUnboundedQuadrantsExerciseRayClipping(testCase)
            % f(x,y)=x on the quadrant {x>=0,y>=0} (apex (0,0)); g(x,y)=y on the quadrant
            % {x>=-1,y>=-1} (apex (-1,-1)). The first quadrant is a strict SUBSET of the second,
            % so the overlap is exactly the first quadrant, and h(x,y)=x+y there -- ground truth
            % independently computable since the domains nest.
            f = addQuaPolyTest.quadrant([0 0], [0 1 0]);
            g = addQuaPolyTest.quadrant([-1 -1], [0 1 0]);
            h = f.add(g);

            testCase.verifyClass(h, 'QuaPoly');
            testCase.verifyEqual(h.eval([1 1]), 1+1, 'AbsTol', 1e-10);      % inside both
            testCase.verifyEqual(h.eval([-0.5 -0.5]), Inf);                 % inside g only, not f
            testCase.verifyEqual(h.eval([-2 -2]), Inf);                     % inside neither
        end
    end

    methods (Static)
        function p = squareLinear(corner, coeffXYconst)
            % A single unit^2*2 axis-aligned square [corner, corner+(2,2)], CCW, one face,
            % f(x,y) = coeffXYconst(1)*x + coeffXYconst(2)*y + coeffXYconst(3).
            c = corner;
            V = [c; c+[2 0]; c+[2 2]; c+[0 2]];
            E = [1 2 1; 2 3 1; 3 4 1; 4 1 1];
            F = [1 0; 1 0; 1 0; 1 0];
            f = zeros(1,10);
            f(8) = coeffXYconst(1); f(9) = coeffXYconst(2); f(10) = coeffXYconst(3);
            p = QuaPoly(V, E, f, F);
        end

        function p = quadrant(apex, coeffXYconst)
            % The unbounded quadrant {x>=apex(1), y>=apex(2)} (2 rays from apex, +x and +y
            % direction), f(x,y) = coeffXYconst(1)*x + coeffXYconst(2)*y + coeffXYconst(3).
            V = [apex; apex+[1 0]; apex+[0 1]];
            E = [1 2 0; 1 3 0];
            F = [1 0; 0 1];
            f = zeros(1,10);
            f(8) = coeffXYconst(1); f(9) = coeffXYconst(2); f(10) = coeffXYconst(3);
            p = QuaPoly(V, E, f, F);
        end
    end
end
