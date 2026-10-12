classdef QuaPolEvalTest < matlab.unittest.TestCase
% QuaPol.eval against its DEFINITION: a point of the domain gets the polynomial of the piece it
% lies in, a point outside gets +inf. The domain is read here independently of eval -- each case
% names its points by hand, inside / on the boundary / outside.
%
% BUCKET: fast (a few dozen double evaluations).
%
% Found 2026-10-11 checking conjQ against a numeric sup: eval returned +inf at points IN the domain
% on a ray (read as the segment between its two vertices), on a non-convex piece (the intersection
% of the edge half-planes is only the polygon's kernel), and on a needle.

    methods (Test)

        function aNEEDLEIsFiniteAtItsVertexOnly(t)
            P = QuaPol.examples2();  p = P{1};                         % V = (0,0), f = x + 2y + 3
            t.verifyEqual(p.eval([0 0]), QuaPol.evalPoly(p.f, [0 0]));
            t.verifyEqual(p.eval([1 0; 0 -1]), [inf; inf]);
        end

        function aRAYRunsThroughItsSecondVertexAndBeyond(t)
            P = QuaPol.examples2();  p = P{3};                         % ray from (0,0) along (1,0)
            in = [0 0; 0.5 0; 1 0; 2 0; 3 0; 1e3 0];
            t.verifyEqual(p.eval(in), QuaPol.evalPoly(p.f, in), 'AbsTol', 1e-9);
            t.verifyEqual(p.eval([-1 0; 2 1]), [inf; inf]);
        end

        function aRAYPointingTowardNegativeCoordinatesIsReadByItsDirection(t)
        % a ray from (1,1) through (0,2): direction (-1,1)
            p = QuaPol([1 1; 0 2], [1 2 0], [1 0 0 0 0 0], [0 0]);
            in = [1 1; 0 2; -3 5];
            t.verifyEqual(p.eval(in), QuaPol.evalPoly(p.f, in), 'AbsTol', 1e-9);
            t.verifyEqual(p.eval([2 0; 0 0]), [inf; inf]);
        end

        function aNOTCHEDPolygonIsTheWholePolygonNotItsKernel(t)
        % (0,0) (1,1) (2,1) (2,2) (0,2): reflex at (1,1)
            P = QuaPol.examplesNonconvex();  p = P{1};
            in = [0 0; 2 1; 0.5 0.5; 1 1; 0.2 1; 1 1.5; 2 2; 0 1];
            t.verifyEqual(p.eval(in), QuaPol.evalPoly(p.f, in), 'AbsTol', 1e-9);
            t.verifyEqual(p.eval([1 0.5; 1.5 0.9; 3 1.5; -0.1 1]), inf(4,1));
        end

        function aNONCONVEXUnboundedPieceIsEverythingAboveItsBoundary(t)
        % the W epigraph: rays from (-1,0) through (-2,1) and from (1,0) through (2,1), segments
        % to the peak (0,1); the piece lies above
            P = QuaPol.examplesNonconvex();  p = P{2};
            in = [0 2; -3 2.5; -1 0.5; 1 0.5; -5 4; 0 1; 50 60];
            t.verifyEqual(p.eval(in), QuaPol.evalPoly(p.f, in), 'AbsTol', 1e-9);
            t.verifyEqual(p.eval([0 0.5; -3 1; 5 3; 0 -1]), inf(4,1));
        end

        function aCONVEXPieceIsUnchanged(t)
        % regression guard: every example with faces, on a grid, agrees with the half-plane reading
        % wherever the domain is convex (the old code was right there)
            P = QuaPol.examples();
            [X, Y] = meshgrid(-3:0.5:3);  x = [X(:) Y(:)];
            for k = 1:numel(P)
                p = P{k};
                if p.nv == 0 || p.nf == 0 || ~p.createDom().isConvex, continue; end
                for i = 1:p.nf
                    t.verifyTrue(QuaPolEvalTest.halfPlaneAgree(p, i, x), sprintf('examples{%d} face %d', k, i));
                end
            end
        end
    end

    methods (Static)
        function ok = halfPlaneAgree(p, i, x)
        % points strictly inside face i by its half-planes get that face's value (or another
        % face's, on a shared boundary -- excluded by the strict test)
            [a, b] = deal(p.V(p.E(:,1),:), p.V(p.E(:,2),:));
            c = QuaPol.lineEquation(a, b);
            s = ([x ones(size(x,1),1)] * c(abs(p.P{i}),:)') .* sign(p.P{i});
            strict = all(s < -1e-9, 2);
            v = p.eval(x(strict,:));
            ok = all(abs(v - QuaPol.evalPoly(p.f(i,:), x(strict,:))) <= 1e-9 * max(1, abs(v)));
        end
    end
end
