classdef QuaParTest < matlab.unittest.TestCase
    % Tests for QuaPar (quadratic on a parabolic subdivision).
    % Covers: equivalence to QuaPol on linear-edge inputs, full-domain evaluation,
    % parabola validation, evaluation across a genuine parabolic edge, and the conic helpers.

    methods (Test)
        function linearEdgesMatchQuaPol(testCase)
            % A QuaPar built from linear edges must evaluate identically to the QuaPol
            % built from the same data (polyhedral subdivision = special parabolic subdivision).
            pPoly = QuaPol.oneNorm();      % |x| + |y|, finite everywhere
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

        function evalLocatesAPointExactlyAtItsOwnVertex(testCase)
        % THE "QuaPar.eval exactly AT a vertex" defect, reproduced deterministically.
        %
        % A vertex is a point where several edge conics ought to evaluate to exactly zero. In
        % floating point they come out at +-1e-16 instead, and the point location eval used to
        % do -- `all(vals <= 0, 2)`, an exact comparison with no tolerance -- then admitted the
        % point to NO face at all, so eval returned its Inf initialization. The fix is a
        % conic-magnitude-relative tolerance; this test pins it.
        %
        % The mesh below is not hand-made: it is case 2 of
        % sweepQuaParEvalAtVertices(20260802, 200), the first subdivision that sweep produces
        % with a vertex the exact test cannot locate. That sweep is committed alongside this
        % test and reports 225 of 1205 vertices (18.7%) unlocatable under the exact test and 0
        % under the current one, with every ring of radius 1e-8 around them evaluating
        % correctly -- the signature of an exact test applied to a quantity that is only zero
        % in exact arithmetic. See SUPPORT_MATRIX.md sections 0.1 and 7.
        %
        % ONE quadratic on all five faces, so the function is globally smooth and any face that
        % claims the point gives the same value: what this asserts is purely point location.
            V = [ 0.37841944881914641   1.7270572762519276
                 -0.30045265946634558   1.023881684721043
                 -1.0491536570933224   -0.89771602108316162
                 -1.0550773423633844   -1.2557220832312725
                 -0.75039498409826855  -1.6784901757118644
                  1.1867556445418206   -0.99215711362818615
                  1.8670090437313969    1.0619937756582127];
            coef = [0.28440589539080674 -0.45957766659064853 -0.17312711301216793 ...
                   -0.72098049529602004 -0.86330613273515178  1.3392036513028369];
            n = size(V,1);
            E = zeros(0,3); F = zeros(0,2);
            for j = 1:n
                j2 = mod(j,n) + 1;
                if j == 1,      fk = 1;
                elseif j == n,  fk = n - 2;
                else,           fk = j - 1;
                end
                E(end+1,:) = [j j2 1]; %#ok<AGROW>
                F(end+1,:) = [fk 0];   %#ok<AGROW>
            end
            for j = 3:n-1
                E(end+1,:) = [1 j 1];   %#ok<AGROW>
                F(end+1,:) = [j-1 j-2]; %#ok<AGROW>
            end
            q = QuaPar(V, E, repmat(coef, n-2, 1), F);

            for iv = 1:q.nv
                v = q.V(iv,:);
                want = QuaPar.evalPoly(coef, v);
                got = q.eval(v);
                testCase.verifyTrue(isfinite(got), sprintf( ...
                    'eval returned %g at its OWN vertex %d (%.17g,%.17g)', got, iv, v(1), v(2)));
                testCase.verifyEqual(got, want, 'AbsTol', 1e-9, sprintf( ...
                    'eval at vertex %d', iv));
            end
        end

        function scalarMulAndNegateScaleCoefficientsOnly(testCase)
            p = QuaPar.energy();   % 0.5*(x^2+y^2), full domain
            S = [1 2; -3 0.5];
            q = p.scalarMul(3);
            testCase.verifyEqual(q.nv, p.nv);   % domain untouched
            testCase.verifyEqual(q.eval(S), 3*p.eval(S), 'AbsTol', 1e-12);
            r = p.negate();
            testCase.verifyEqual(r.eval(S), -p.eval(S), 'AbsTol', 1e-12);
        end
    end
end
