classdef QuaParTest < matlab.unittest.TestCase
    % Tests for QuaPar (quadratic on a parabolic subdivision).
    % Covers: equivalence to QuaPoly on linear-edge inputs, full-domain evaluation,
    % parabola validation, evaluation across a genuine parabolic edge, and the conic helpers.

    methods (Test)
        function linearEdgesMatchQuaPoly(testCase)
            % A QuaPar built from linear edges must evaluate identically to the QuaPoly
            % built from the same data (polyhedral subdivision = special parabolic subdivision).
            pPoly = QuaPoly.oneNorm();      % |x| + |y|, finite everywhere
            pPar  = QuaPar.oneNorm();       % same data via the (V,E,f,F) constructor (Ec=0)
            [X,Y] = meshgrid(linspace(-1.7,1.7,23), linspace(-1.3,1.3,19));
            P = [X(:)+0.013, Y(:)+0.007];   % small offset to avoid exact edge points
            testCase.verifyEqual(pPar.eval(P), pPoly.eval(P), 'AbsTol', 1e-10);
        end

        function fullDomainQuadraticEval(testCase)
            p = QuaPar.energy();            % 1/2(x^2 + y^2)
            S = [0 0; 1 1; -2 3; 0.5 -1.5];
            expected = 0.5*(S(:,1).^2 + S(:,2).^2);
            testCase.verifyEqual(p.eval(S), expected, 'AbsTol', 1e-12);
        end

        function rejectsNonParabolicConic(testCase)
            % A circle x^2 + y^2 - 1 = 0 has b^2-4ac = -4 ~= 0 and must be rejected.
            V  = [-1 0; 1 0];
            E  = [1 2 1; 1 2 1];
            Ec = [1 0 1 0 0 -1;     % circle: NOT a parabola
                  0 0 0 0 0  0];
            f  = [1 0 1 0 0 0];
            F  = [1 0; 1 0];
            P  = {[1 2]};
            testCase.verifyError(@() QuaPar(V,E,Ec,f,F,P), 'QuaPar:notParabola');
        end

        function acceptsParabolaAndLineConics(testCase)
            % y = x^2 (parabola) and x = 0 (line, degenerate parabola) are both accepted.
            testCase.verifyTrue(QuaPar.isParabola([1 0 0 0 -1 0]));   % x^2 - y = 0
            testCase.verifyTrue(QuaPar.isParabola([0 0 0 1 0 0]));    % x = 0 (a line)
            testCase.verifyFalse(QuaPar.isParabola([1 0 1 0 0 -1]));  % circle
            testCase.verifyFalse(QuaPar.isParabola([0 0 0 0 0 0]));   % all zero (not a curve)
        end

        function evalConicValues(testCase)
            % a x^2 + b xy + c y^2 + d x + e y + f at sample points.
            c = [1 0 0 0 -1 0];                     % x^2 - y
            X = [2 1; 0 0; -3 9];
            testCase.verifyEqual(QuaPar.evalConic(c, X), [3; 0; 0], 'AbsTol', 1e-12);
        end

        function parabolicSliceEval(testCase)
            % f = 1/2(x^2+y^2) on { y >= x^2 } ∩ { y <= 1 }, +inf elsewhere.
            p = QuaPar.parabolicSlice();
            % interior points
            testCase.verifyEqual(p.eval([0 0.5]),  0.125, 'AbsTol', 1e-12);  % 1/2(0+0.25)
            testCase.verifyEqual(p.eval([0.5 0.8]), 0.445, 'AbsTol', 1e-12); % 1/2(0.25+0.64)
            % boundary point (inclusive)
            testCase.verifyEqual(p.eval([0 1]),    0.5,   'AbsTol', 1e-12);  % 1/2(0+1)
            % exterior: above the top line, and below the parabola
            testCase.verifyEqual(p.eval([0 2]),    Inf);
            testCase.verifyEqual(p.eval([0.95 0.5]), Inf);
        end

        function parabolicSliceVectorized(testCase)
            p = QuaPar.parabolicSlice();
            S = [0 0.5; 0 2; 0.95 0.5; 0.5 0.8];
            expected = [0.125; Inf; Inf; 0.445];
            testCase.verifyEqual(p.eval(S), expected, 'AbsTol', 1e-12);
        end

        function parabolicFaceAutoBuildsP(testCase)
            % A parabolic "triangle" built with the 5-arg constructor (no explicit P): edge 1 is
            % the arc of y=x^2 from (-1,1) to (1,1); edges 2,3 are the straight edges to (0,2).
            % The conic y-x^2 is oriented >0 above the arc = left of the directed edge (1->2), so
            % createP/orderEdges build P automatically and eval locates points correctly.
            V  = [-1 1; 1 1; 0 2];
            E  = [1 2 1; 2 3 1; 3 1 1];
            Ec = [-1 0 0 0 1 0; 0 0 0 0 0 0; 0 0 0 0 0 0];   % edge1: y - x^2 = 0
            f  = [1 0 1 0 0 0];                               % 1/2(x^2+y^2)
            p  = QuaPar(V, E, Ec, f, [1 0; 1 0; 1 0]);        % no P supplied
            testCase.verifyEqual(p.nf, 1);
            % interior points (above parabola y>=x^2, below both lines to (0,2))
            Sin = [0 1.5; 0 1.2; 0 0.5; 0.3 0.8];
            testCase.verifyEqual(p.eval(Sin), 0.5*(Sin(:,1).^2 + Sin(:,2).^2), 'AbsTol', 1e-12);
            % exterior: below the parabola, above the apex, and outside a top edge
            testCase.verifyEqual(p.eval([0 -1; 0 3; -0.5 1.7]), [Inf; Inf; Inf]);
        end

        function fullDomainConjugateViaCplq(testCase)
            % conj of a full-domain strictly convex quadratic QuaPar works through conjCPLQ.
            p = QuaPar.energy();
            q = p.conj('cplq');
            S = [0 0; 1 -1; 2 3];
            testCase.verifyEqual(q.eval(S), 0.5*(S(:,1).^2+S(:,2).^2), 'AbsTol', 1e-12);
        end
    end
end
