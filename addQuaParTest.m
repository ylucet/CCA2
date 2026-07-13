classdef addQuaParTest < matlab.unittest.TestCase
    % Tests for QuaPar.add / addQuaPar.m: pointwise sum of two QuaPar functions, generalizing
    % addQuaPoly.m to allow one curved (parabolic-arc) edge per face. See addQuaPar.m's header for
    % the algorithm (uniform 6-coefficient constraint rows, straight edges parametrized in t,
    % the one curved edge parametrized in u via an axis-rotation of its parabola).

    methods (Test)
        function bothFullDomainSumsCoefficients(testCase)
            f = QuaPar([1 0 1 0 0 0]);   % 0.5*(x^2+y^2), full domain
            g = QuaPar([1 0 1 0 0 0]);
            h = f.add(g);
            testCase.verifyEqual(h.nv, 0);
            testCase.verifyEqual(h.nf, 1);
            testCase.verifyEqual(h.eval([1 2]), 1^2 + 2^2, 'AbsTol', 1e-10);
        end

        function straightOnlyMatchesQuaPolyBehaviour(testCase)
            % Both operands purely polyhedral (all-zero Ec): addQuaPar must behave exactly like
            % addQuaPoly's own overlapping-squares case (same regression the QuaPoly test covers).
            f = addQuaParTest.squarePar([0 0], [0 1 0]);   % f(x,y)=x on [0,2]x[0,2]
            g = addQuaParTest.squarePar([1 1], [0 1 0]);   % g(x,y)=y on [1,3]x[1,3]
            h = f.add(g);
            testCase.verifyClass(h, 'QuaPar');
            testCase.verifyEqual(h.eval([1.5 1.5]), 1.5+1.5, 'AbsTol', 1e-10);   % inside the overlap
            testCase.verifyEqual(h.eval([0.5 0.5]), Inf);                        % inside f only
            testCase.verifyEqual(h.eval([1.5 0.5]), Inf);   % inside f, outside g (y<1): regression
            testCase.verifyEqual(h.eval([0.5 1.5]), Inf);   % inside f, outside g (x<1): regression
        end

        function curvedFacePlusFullDomainPreservesTheArc(testCase)
            % Adding a full-domain quadratic never clips anything (addQuaPar's trivial-case
            % shortcut), so the curved edge must survive completely unaffected.
            trough = addQuaParTest.trough();
            quad = QuaPar([1 0 1 0 0 0]);   % 0.5*(x^2+y^2)
            h = quad.add(trough);
            testCase.verifyClass(h, 'QuaPar');
            testCase.verifyTrue(any(h.Ec(:) ~= 0));
            testCase.verifyEqual(h.eval([0 2]), 0.5*(0+4), 'AbsTol', 1e-10);    % inside trough
            testCase.verifyEqual(h.eval([0 -1]), Inf);                          % below the parabola
            testCase.verifyEqual(h.eval([1 3]), 0.5*(1+9), 'AbsTol', 1e-10);
            testCase.verifyEqual(h.eval([1 0.5]), Inf);                         % below parabola at x=1
        end

        function curvedFaceClippedByStraightSquare(testCase)
            % The real geometry test: trough={-2<=x<=2, y>=x^2} (one curved edge, two rays) added
            % to square=[-1,1]x[0,5] (straight, f(x,y)=x). Ground truth: for x in [-1,1], trough's
            % own y>=x^2 already implies y in [0,5] (since x^2<=1<=5), so the overlap is EXACTLY
            % {-1<=x<=1, x^2<=y<=5} -- independently computable since the domains nest that way.
            trough = addQuaParTest.trough();
            square = addQuaParTest.squarePar([-1 0], [1 0 0], 2, 5);   % f(x,y)=x on [-1,1]x[0,5]
            h = trough.add(square);
            testCase.verifyClass(h, 'QuaPar');
            testCase.verifyTrue(any(h.Ec(:) ~= 0));
            testCase.verifyEqual(h.eval([0 2]), 0+0, 'AbsTol', 1e-10);          % inside both
            testCase.verifyEqual(h.eval([0 -0.5]), Inf);                       % below the parabola
            testCase.verifyEqual(h.eval([0.5 0.3]), 0+0.5, 'AbsTol', 1e-10);   % y=0.3>=0.25=x^2
            testCase.verifyEqual(h.eval([1.5 3]), Inf);                        % x outside square
            testCase.verifyEqual(h.eval([-1 2]), 0+(-1), 'AbsTol', 1e-10);
            testCase.verifyEqual(h.eval([1 4.9]), 0+1, 'AbsTol', 1e-10);
            testCase.verifyEqual(h.eval([1 -0.5]), Inf);                       % y<x^2=1 at x=1
        end
    end

    methods (Static)
        function p = squarePar(corner, coeffXYconst, side, heightArg)
            % Axis-aligned rectangle [corner, corner+(side,height)] (default 2x2), CCW, one face,
            % f(x,y) = coeffXYconst(1)*x + coeffXYconst(2)*y + coeffXYconst(3). Straight edges only
            % (all-zero Ec), as a QuaPar.
            if nargin < 3, side = 2; end
            if nargin < 4, heightArg = side; end
            c = corner;
            V = [c; c+[side 0]; c+[side heightArg]; c+[0 heightArg]];
            E = [1 2 1; 2 3 1; 3 4 1; 4 1 1];
            F = [1 0; 1 0; 1 0; 1 0];
            f = zeros(1,10);
            f(8) = coeffXYconst(1); f(9) = coeffXYconst(2); f(10) = coeffXYconst(3);
            p = QuaPar(V, E, zeros(4,6), f, F);
        end

        function p = trough()
            % {-2<=x<=2, y>=x^2}: unbounded upward, convex, one face (f=0). Boundary: the parabola
            % arc A(-2,4)->B(2,4) (curved, y=x^2) plus a ray straight up from each of A and B.
            V = [-2 4; 2 4; -2 5; 2 5];      % A, B, upA(dir ref), upB(dir ref)
            E = [1 2 1; 1 3 0; 2 4 0];        % A-B curved segment; ray A->upA; ray B->upB
            Ec = [-1 0 0 0 1 0; 0 0 0 0 0 0; 0 0 0 0 0 0];   % edge1: y-x^2, >0 on left of A->B
            F = [1 0; 0 1; 1 0];
            f = zeros(1,10);
            p = QuaPar(V, E, Ec, f, F);
        end
    end
end
