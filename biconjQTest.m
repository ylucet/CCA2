classdef biconjQTest < matlab.unittest.TestCase
% Unit tests for biconjQ, the EXACT biconjugate f** = cl co f.
%
% BUCKET: fast (integer arithmetic only; no symbolic call anywhere).
%
% THE ORACLE IS THE DEFINITION OF A CONVEX ENVELOPE, in the three parts plqCheck.m also uses:
%   co f <= f with EQUALITY at the extreme points of the domain; co f is CONVEX; and co f is the
%   LARGEST convex minorant, so no affine function lying below f may exceed it anywhere. None of
%   the three mentions a lower hull, which is how biconjQ computes -- so they check the answer
%   rather than restate the method.

    methods (Test)

        function aConvexFunctionIsItsOwnEnvelope(testCase)
        % co f = f when f is convex, so there is nothing to compute -- and getting this wrong is
        % the silent failure ALGORITHM.md warns about: biconj returning a non-convex f unchanged.
            f = biconjQTest.square([0 0 0 0, 2, 0, 2, -1, 3, 1]);   % strictly convex
            h = biconjQ(f);
            X = biconjQTest.insideSquare(200);
            testCase.verifyEqual(h.eval(X), f.eval(X), 'AbsTol', 1e-12);
        end

        function anAffineFunctionIsItsOwnEnvelopeToo(testCase)
        % The H = 0 edge of both branches: affine is convex AND concave, so whichever route runs
        % must return f unchanged.
            f = biconjQTest.square([0 0 0 0, 0, 0, 0, -1, 2, 3]);
            h = biconjQ(f);
            X = biconjQTest.insideSquare(200);
            testCase.verifyEqual(h.eval(X), f.eval(X), 'AbsTol', 1e-12);
            testCase.verifyEqual(h.nf, 1, 'an affine function needs exactly one piece');
        end

        function theEnvelopeOfAConcaveQuadraticInterpolatesItsVertexValues(testCase)
        % The defining property of this case: a concave function on a polytope has its whole
        % envelope determined by the VERTEX values, so co f must match f at every corner and lie
        % strictly below it in between.
            f = biconjQTest.square([0 0 0 0, -1, 0, -1, 0, 0, 0]);   % -(x^2+y^2)/2
            h = biconjQ(f);
            V = [0 0; 1 0; 1 1; 0 1];
            testCase.verifyEqual(h.eval(V), f.eval(V), 'AbsTol', 1e-12, ...
                'equality at every extreme point');

            X = biconjQTest.insideSquare(300);
            testCase.verifyLessThanOrEqual(h.eval(X), f.eval(X) + 1e-12, 'co f <= f');
            testCase.verifyLessThan(h.eval([0.5 0.5]), f.eval([0.5 0.5]) - 1e-9, ...
                'and strictly below in the interior, since q is strictly concave there');
        end

        function theEnvelopeIsConvexAndIsTheLargestConvexMinorant(testCase)
        % Two definition-level properties that no cell-by-cell check implies. The second is what
        % separates the true envelope from any convex function that merely sits below f.
            f = biconjQTest.polygon([0 0; 3 0; 4 3; 1 4], [0 0 0 0, -1, 0, -1, 0, 0, 0]);
            h = biconjQ(f);
            V = [0 0; 3 0; 4 3; 1 4];

            rng(20260904);
            A = biconjQTest.inHull(V, 200);  B = biconjQTest.inHull(V, 200);
            testCase.verifyLessThanOrEqual(h.eval((A+B)/2), (h.eval(A) + h.eval(B))/2 + 1e-9, ...
                'the envelope must be convex');

            qv = f.eval(V);
            X  = biconjQTest.inHull(V, 300);
            hx = h.eval(X);
            for t = 1:25
                g = randn(1,2);
                b = min(qv - V*g.');                 % the highest affine minorant with slope g
                testCase.verifyLessThanOrEqual(X*g.' + b, hx + 1e-9, ...
                    'an affine function below f may never exceed co f');
            end
        end

        function aCollinearVertexIsDominatedAndTheEnvelopeSaysSo(testCase)
        % NOT a defect, and it cost a diagnosis to establish: a vertex lying ON the segment between
        % its neighbours is not an extreme point, so for a strictly concave q its lifted point sits
        % strictly ABOVE the chord and is not on the lower hull. co f is then strictly BELOW q
        % there. The dual of conjQ's aDominatedVertexContributesNoCellAtAll.
            V = [0 0; 1 0; 2 0; 2 2; 0 2];           % (1,0) is collinear with (0,0) and (2,0)
            f = biconjQTest.polygon(V, [0 0 0 0, -2, 0, -2, 0, 0, 0]);
            h = biconjQ(f);
            ext = [1 3 4 5];
            testCase.verifyEqual(h.eval(V(ext,:)), f.eval(V(ext,:)), 'AbsTol', 1e-12, ...
                'equality at the four EXTREME points');
            testCase.verifyLessThan(h.eval(V(2,:)), f.eval(V(2,:)) - 1e-9, ...
                'and strictly below at the collinear one, which the hull does not use');
        end

        function theEnvelopeCoefficientsAreExactRationals(testCase)
        % The point of the exact route: the stored coefficients are integers over integers, not
        % doubles that happen to be close.
            f = biconjQTest.square([0 0 0 0, -1, 0, -1, 0, 0, 0]);
            h = biconjQ(f);
            testCase.verifyEqual(h.fN, round(h.fN), 'AbsTol', 0);
            testCase.verifyTrue(all(h.fD > 0));
            testCase.verifyTrue(all(h.fN(:,5) == 0 & h.fN(:,6) == 0 & h.fN(:,7) == 0), ...
                'the envelope of a concave q is piecewise AFFINE -- no curvature may appear');
        end

        % ---- what is refused, by name ------------------------------------------------------

        function anIndefinitePieceIsRefusedByName(testCase)
        % The trigger for AlgAlg: this is where the first face that cannot be written rationally
        % appears, since an affine cell of a general f** names a dual vertex of degree up to 4.
            f = biconjQTest.square([0 0 0 0, 1, 0, -1, 0, 0, 0]);
            testCase.verifyError(@() biconjQ(f), 'PLQ:biconjQ:notImplemented');
        end

        function aMultiPieceNonConvexInputIsRefusedByName(testCase)
        % The envelope COUPLES pieces -- the convex hull of a union is not determined piecewise --
        % so unlike the conjugate it is not a fold, and refusing is the honest answer.
            V = [0 0; 1 0; 1 1; 0 1];
            E = [1 2 1; 2 3 1; 3 4 1; 4 1 1; 1 3 1];
            F = [1 0; 1 0; 2 0; 2 0; 1 2];
            f = QuaPol(V, E, [0 0 0 0, -1, 0, -1, 0, 0, 0; 0 0 0 0, -2, 0, -2, 0, 0, 0], F);
            testCase.verifyError(@() biconjQ(f), 'PLQ:biconjQ:notImplemented');
        end

        function anInexactInputIsRefusedBeforeAnythingElse(testCase)
            testCase.verifyError(@() biconjQ(QuaPol([sqrt(2) 0 1 0 0 0])), 'PLQ:QuaPol:notExact');
        end

        function aFalseConvexityAssertionIsRefusedLoudly(testCase)
        % ALGORITHM.md's rule: the free necessary condition is still checked, because the failure
        % it prevents is silent -- biconj returning a non-convex f as its own convex envelope.
            f = biconjQTest.square([0 0 0 0, 1, 0, -1, 0, 0, 0]);
            f.fIsConvex = true;
            testCase.verifyError(@() biconjQ(f), 'PLQ:biconjQ:notConvex');
        end
    end

    methods (Static)
        function f = square(fc)
            f = QuaPol([0 0; 1 0; 1 1; 0 1], [1 2 1; 2 3 1; 3 4 1; 4 1 1], fc, ...
                       [1 0; 1 0; 1 0; 1 0]);
        end

        function f = polygon(V, fc)
            m = size(V,1);
            f = QuaPol(V, [(1:m).', [2:m,1].', ones(m,1)], fc, [ones(m,1), zeros(m,1)]);
        end

        function X = insideSquare(n)
            rng(20260904);
            X = rand(n,2);
        end

        function X = inHull(V, n)
            m = size(V,1);  X = zeros(n,2);
            for i = 1:n
                w = rand(m,1);  w = w/sum(w);
                X(i,:) = w.' * V;
            end
        end
    end
end
