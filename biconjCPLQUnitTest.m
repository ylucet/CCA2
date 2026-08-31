classdef biconjCPLQUnitTest < matlab.unittest.TestCase
% Branch-level tests for `biconjCPLQ`, at 33.7% under fast+normal.
%
% WHY THE ENTRY POINT AND NOT THE HELPERS. The routines that are dark here -- `oned`, `onedEnv`,
% `diamondEnvelope`, `convexEnough`, `separableEnvelopeCoefs`, `hasCurvedEdge`, `asQuaPol` -- are
% file-private subfunctions of biconjCPLQ.m and cannot be called from anywhere else. The smallest
% callable unit that contains them IS `biconjCPLQ`, so each test below drives it with the ONE
% input shape that selects one branch, and asserts that branch's own answer. That is a different
% thing from a pipeline test: nothing here runs triangulate, conjugate or Step 3, and when one
% goes red it names the branch.
%
% THE ASSERTION IS THE DEFINITION OF f**, in every case. For a proper function bounded below,
% f** = conv f, so the answer must
%   * UNDERESTIMATE f on the domain -- the direction that is a definite defect,
%   * be CONVEX along sampled segments,
%   * equal f wherever f is already convex there -- in particular at the vertices of a domain on
%     which the envelope is the affine interpolant.
% Closed forms are asserted against an independent formula, never against the routine's own
% output.
%
% BUCKET: fast. Every fixture is a single bounded face and every branch below is closed-form; the
% general Case C path, which is symbolic, is deliberately not exercised here.

    methods (Static)
        function [V, E, F] = boxMesh(lo, hi)
        % One face, four segments, counter-clockwise -- the axis-aligned box shape the separable
        % branch requires.
            V = [lo(1) lo(2); hi(1) lo(2); hi(1) hi(2); lo(1) hi(2)];
            E = [1 2 1; 2 3 1; 3 4 1; 4 1 1];
            F = [1 0; 1 0; 1 0; 1 0];
        end

        function v = evalCoefRow(c, x, y)
        % The stored WEIGHTED basis, written out: a row [a b c d e f] means
        % a/2 x^2 + b xy + c/2 y^2 + d x + e y + f. See meshPredicateTest for why the halves are
        % on the pure squares only.
            c = c(end-5:end);
            v = 0.5*c(1)*x.^2 + c(2)*x.*y + 0.5*c(3)*y.^2 + c(4)*x + c(5)*y + c(6);
        end

        function P = boxSample(lo, hi, n)
            [gx, gy] = meshgrid(linspace(lo(1), hi(1), n), linspace(lo(2), hi(2), n));
            P = [gx(:), gy(:)];
        end
    end

    methods (Test)

        % =========================================================================================
        % convexEnough -- f already convex, so f** = f and nothing else may happen
        % =========================================================================================
        function anAlreadyConvexQuadraticIsItsOwnBiconjugate(testCase)
        % The cheapest branch, and the one whose failure would be silent: if the short-cut ever
        % returned something OTHER than f for a convex f, every value would still look plausible.
        % Asserted pointwise against f itself.
            f = [1 0 1 0 0 0];                             % x^2/2 + y^2/2 in the weighted basis
            p = QuaPol(f);                                 % full domain, strictly convex
            h = biconjCPLQ(p);
            P = [0 0; 1 2; -3 0.5; 4 -4; 0.25 0.25];
            for i = 1:size(P,1)
                testCase.verifyEqual(h.eval(P(i,:)), p.eval(P(i,:)), 'AbsTol', 1e-9, sprintf( ...
                    'f** must equal f for an already-convex f, at (%g,%g)', P(i,1), P(i,2)));
            end
        end

        % =========================================================================================
        % separableEnvelopeCoefs -- a SEPARABLE quadratic over an axis-aligned box
        % =========================================================================================
        function theSeparableBoxEnvelopeIsTheChordInEachVariable(testCase)
        % A separable f(x,y) = phi(x) + psi(y) on a box has conv f = conv(phi) + conv(psi), and
        % for a CONCAVE quadratic on an interval the convex envelope is the CHORD. So the answer
        % is available in closed form independently of the routine, and that is what it is checked
        % against -- not against a coefficient row.
            lo = [-1 -2]; hi = [2 3];
            [V, E, F] = biconjCPLQUnitTest.boxMesh(lo, hi);
            % -x^2 - 2y^2 in the weighted basis: c1/2 = -1, c3/2 = -2.
            a = -1; b = -2;
            p = QuaPol(V, E, [2*a 0 2*b 0 0 0], F);
            h = biconjCPLQ(p);

            chord = @(c, l, hh, t) ((c*hh^2 - c*l^2)/(hh - l))*(t - l) + c*l^2;
            P = biconjCPLQUnitTest.boxSample(lo, hi, 11);
            for i = 1:size(P,1)
                x = P(i,1); y = P(i,2);
                want = chord(a, lo(1), hi(1), x) + chord(b, lo(2), hi(2), y);
                got  = h.eval([x y]);
                testCase.verifyEqual(got, want, 'AbsTol', 1e-7, sprintf( ...
                    'the separable box envelope at (%g,%g)', x, y));
            end
            % ...and it is a genuine underestimator that touches at the four corners, which is
            % where a chord in each variable must meet the function.
            for i = 1:size(V,1)
                fv = a*V(i,1)^2 + b*V(i,2)^2;
                testCase.verifyEqual(h.eval(V(i,:)), fv, 'AbsTol', 1e-7, sprintf( ...
                    'the envelope must touch f at the corner (%g,%g)', V(i,1), V(i,2)));
            end
        end

        function theSeparableBoxEnvelopeUnderestimatesAndIsConvex(testCase)
        % The two properties that hold for ANY f**, asserted on the same fixture so that a wrong
        % closed form is caught even where it happens to match the chord at the corners.
            lo = [0 0]; hi = [3 2];
            [V, E, F] = biconjCPLQUnitTest.boxMesh(lo, hi);
            c = [-2 0 -1 1 -3 4];                          % concave separable plus a linear part
            p = QuaPol(V, E, c, F);
            h = biconjCPLQ(p);
            P = biconjCPLQUnitTest.boxSample(lo, hi, 13);
            worst = inf;
            for i = 1:size(P,1)
                fv = biconjCPLQUnitTest.evalCoefRow(c, P(i,1), P(i,2));
                worst = min(worst, fv - h.eval(P(i,:)));
            end
            testCase.verifyGreaterThanOrEqual(worst, -1e-7, sprintf( ...
                'f** exceeds f by %.3g somewhere on the box -- it is not an underestimator', -worst));

            % Convex along sampled chords: the midpoint may not sit above the chord.
            rng(20260831);
            idx = randi(size(P,1), 60, 2);
            over = -inf;
            for k = 1:size(idx,1)
                A = P(idx(k,1),:); B = P(idx(k,2),:);
                over = max(over, h.eval((A+B)/2) - (h.eval(A) + h.eval(B))/2);
            end
            testCase.verifyLessThanOrEqual(over, 1e-7, sprintf( ...
                'f** is not convex: a midpoint sits %.3g above its chord', over));
        end

        % =========================================================================================
        % diamondEnvelope -- a BILINEAR function over a rotated box
        % =========================================================================================
        function theDiamondEnvelopeUnderestimatesAndTouchesAtTheVertices(testCase)
        % x*y over the diamond conv{(1,0),(0,1),(-1,0),(0,-1)}. The branch exists because the
        % rotation u = (x+y)/sqrt2, v = (x-y)/sqrt2 turns x*y into (u^2 - v^2)/2, which is
        % separable on a box -- so the answer must still be an underestimator that meets f at the
        % four vertices, whatever the rotation does internally.
            V = [1 0; 0 1; -1 0; 0 -1];
            E = [1 2 1; 2 3 1; 3 4 1; 4 1 1];
            F = [1 0; 1 0; 1 0; 1 0];
            p = QuaPol(V, E, [0 1 0 0 0 0], F);            % x*y
            h = biconjCPLQ(p);
            testCase.verifyNotEmpty(h, 'the diamond branch returned nothing');

            % inside the diamond: |x| + |y| <= 1
            [gx, gy] = meshgrid(linspace(-1, 1, 21), linspace(-1, 1, 21));
            P = [gx(:), gy(:)];
            P = P(abs(P(:,1)) + abs(P(:,2)) <= 1 - 1e-9, :);
            testCase.verifyGreaterThan(size(P,1), 0, 'nothing sampled inside the diamond');
            worst = inf;
            for i = 1:size(P,1)
                worst = min(worst, P(i,1)*P(i,2) - h.eval(P(i,:)));
            end
            testCase.verifyGreaterThanOrEqual(worst, -1e-6, sprintf( ...
                'the diamond envelope exceeds x*y by %.3g -- not an underestimator', -worst));
            for i = 1:size(V,1)
                testCase.verifyEqual(h.eval(V(i,:)), V(i,1)*V(i,2), 'AbsTol', 1e-6, sprintf( ...
                    'the envelope must touch x*y at the vertex (%g,%g)', V(i,1), V(i,2)));
            end
        end

        % =========================================================================================
        % asQuaPol -- a polyhedral QuaPar must be accepted and give the same answer as its QuaPol
        % =========================================================================================
        function aPolyhedralQuaParGivesTheSameAnswerAsTheEquivalentQuaPol(testCase)
        % `asQuaPol` rebuilds a polyhedral QuaPar as a QuaPol because convEnvCPLQ documents a
        % QuaPol input. The two carry the identical mesh, so the biconjugate must not depend on
        % which class it arrived in -- and that is exactly what the rebuild is for.
            lo = [0 0]; hi = [2 2];
            [V, E, F] = biconjCPLQUnitTest.boxMesh(lo, hi);
            c = [-2 0 -2 0 0 0];
            hPol = biconjCPLQ(QuaPol(V, E, c, F));
            hPar = biconjCPLQ(QuaPar(V, E, c, F));         % all-zero Ec: polyhedral
            P = biconjCPLQUnitTest.boxSample(lo, hi, 9);
            for i = 1:size(P,1)
                testCase.verifyEqual(hPar.eval(P(i,:)), hPol.eval(P(i,:)), 'AbsTol', 1e-7, ...
                    sprintf('QuaPar and QuaPol inputs disagree at (%g,%g)', P(i,1), P(i,2)));
            end
        end
    end
end
