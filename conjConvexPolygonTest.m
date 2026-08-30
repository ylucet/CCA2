classdef conjConvexPolygonTest < matlab.unittest.TestCase
% Unit tests for conjConvexPolygon -- the conjugate of a CONVEX quadratic over any convex polygon,
% bounded or not, in closed form.
%
% BUCKET: fast (closed form; the only loops are the brute-force references, which are small).
%
% THE ORACLE IS THE DEFINITION. Every value assertion compares against
%     f*(s) = sup_{x in P} <s,x> - q(x)
% computed either in closed form (the unbounded fixtures separate, so their sup is a one-line
% formula per axis) or numerically, by a coarse scan followed by a pattern search that refines the
% maximiser to full double accuracy. The two-stage reference matters: a SCAN ALONE is only as good
% as its step, and on the wedge fixture that was 3.2e-3 -- absorbing it with a loose tolerance
% would have loosened the tolerance on the code under test by the same amount. Following
% DECISIONS.md 2026-08-18,
% the polygon's own VERTICES are always in the sample -- for a sup that a vertex attains, a grid
% is not a reference, and reading its shortfall as a code defect has already cost this project a
% session.

    methods (Test)

        % ---- the bounded cases, against the definition -------------------------------------------

        function aTriangleAgreesWithTheClosedFormAlreadyInConjPieceCPLQ(testCase)
        % DIFFERENTIAL TEST. conjPieceCPLQ.conjConvexQuadTriangle has conjugated a convex quadratic
        % over a triangle since long before this routine existed, so it is the natural oracle for
        % the n = 3 case: a general construction that disagrees with the special one it generalises
        % is wrong, whatever the definition says about either.
            A = [2 0.5; 0.5 3]; L = [-1; 2]; c = 0.5;
            W = [0 0; 2 0; 0.5 1.5];
            g = conjConvexPolygon(W, [], [], A, L, c);
            q = QuaPol(W, [1 2 1; 2 3 1; 3 1 1], ...
                       [A(1,1) A(1,2) A(2,2) L(1) L(2) c], [1 0; 1 0; 1 0]);
            ref = q.conj('cplq');
            testCase.verifyEqual(g.nf, 7, 'a triangle gives 1 interior + 3 strips + 3 cones');
            S = conjConvexPolygonTest.probes();
            for i = 1:size(S,1)
                testCase.verifyEqual(g.eval(S(i,:)), ref.eval(S(i,:)), 'AbsTol', 1e-9, ...
                    sprintf('at s = (%g,%g)', S(i,1), S(i,2)));
            end
        end

        function boundedPolygonsMatchTheDefinitionSup(testCase)
        % A square, a general quadrilateral and a pentagon. The face count is asserted too, because
        % a missing cell is the failure mode that leaves most probe points right.
            A = [2 0.5; 0.5 3]; L = [-1; 2]; c = 0.5;
            polys = { [0 0; 1 0; 1 1; 0 1], ...
                      [0 0; 2 0; 2.5 1; 0.5 1], ...
                      [0 0; 2 0; 3 1.5; 1 3; -1 1] };
            for k = 1:numel(polys)
                W = polys{k};
                g = conjConvexPolygon(W, [], [], A, L, c);
                testCase.verifyEqual(g.nf, 2*size(W,1) + 1, ...
                    sprintf('polygon %d: n cones + n strips + 1 interior', k));
                S = conjConvexPolygonTest.probes();
                for i = 1:size(S,1)
                    brute = conjConvexPolygonTest.supOverPolygon(S(i,:), W, A, L, c);
                    got = g.eval(S(i,:));
                    testCase.verifyEqual(got, brute, 'AbsTol', 1e-8, ...
                        sprintf('poly %d at s=(%g,%g)', k, S(i,1), S(i,2)));
                end
            end
        end

        function theSubdivisionIsPOLYHEDRALWithNoParabolaAnywhere(testCase)
        % CONJ_FIELD_PROOF.md 7.3 as an assertion, and it is the reason this routine exists: a
        % convex piece's conjugate needs no curved edge, so nothing curved enters Step 3 for it and
        % the arc-vs-arc gap in maxQuaPar is never reached. A QuaPol return type IS that statement.
            A = [1 0; 0 1];
            g = conjConvexPolygon([0 0; 2 0; 3 1.5; 1 3; -1 1], [], [], A, [0;0], 0);
            testCase.verifyEqual(g.kind(), 'QuaPol');
            testCase.verifyTrue(isa(g, 'Pol'));
            testCase.verifyTrue(all(g.Ec(:) == 0), 'a QuaPol pins every edge conic to zero');
        end

        % ---- the unbounded cases, which had no numeric route at all before -----------------------

        function theFirstQuadrantMatchesItsSeparableClosedForm(testCase)
        % P = {x >= 0, y >= 0} written as ONE vertex and two rays. q = (x^2+y^2)/2 separates, so
        %     f*(s) = max(s1,0)^2/2 + max(s2,0)^2/2
        % exactly, with no sampling in the reference at all.
            g = conjConvexPolygon([0 0], [0 1], [1 0], eye(2), [0;0], 0);
            testCase.verifyEqual(g.nf, 4, '1 interior + 2 strips + 1 cone');
            testCase.verifyFalse(g.isDomBounded);
            S = conjConvexPolygonTest.probes();
            for i = 1:size(S,1)
                s = S(i,:);
                truth = 0.5*max(s(1),0)^2 + 0.5*max(s(2),0)^2;
                testCase.verifyEqual(g.eval(s), truth, 'AbsTol', 1e-9, ...
                    sprintf('at s = (%g,%g)', s(1), s(2)));
            end
        end

        function aHalfStripMatchesItsSeparableClosedForm(testCase)
        % P = {0 <= x <= 1, y >= 0}: two vertices and two parallel rays, so the recession cone is a
        % ray rather than a full cone -- the case where the two strip cells are unbounded in the
        % SAME direction.
            g = conjConvexPolygon([0 0; 1 0], [0 1], [0 1], eye(2), [0;0], 0);
            testCase.verifyFalse(g.isDomBounded);
            S = conjConvexPolygonTest.probes();
            for i = 1:size(S,1)
                s = S(i,:);
                if     s(1) <= 0, tx = 0;
                elseif s(1) >= 1, tx = s(1) - 0.5;
                else,             tx = s(1)^2/2;
                end
                truth = tx + 0.5*max(s(2),0)^2;
                testCase.verifyEqual(g.eval(s), truth, 'AbsTol', 1e-9, ...
                    sprintf('at s = (%g,%g)', s(1), s(2)));
            end
        end

        function aHalfPlaneHasNoCornerAtAllAndIsStillHandled(testCase)
        % The extreme of the "vertex that is not a corner" case: {y >= 0} is written as one marker
        % point with two OPPOSITE rays, so the normal cone there has zero width and must be dropped
        % rather than built as an empty face. The two strips then carry the same function, which is
        % what the arithmetic of stripFun says (it is even in the edge direction).
            g = conjConvexPolygon([0 0], [-1 0], [1 0], eye(2), [0;0], 0);
            testCase.verifyEqual(g.nf, 3, 'interior + two strips, and NO cone');
            S = conjConvexPolygonTest.probes();
            for i = 1:size(S,1)
                s = S(i,:);
                truth = 0.5*s(1)^2 + 0.5*max(s(2),0)^2;
                testCase.verifyEqual(g.eval(s), truth, 'AbsTol', 1e-9, ...
                    sprintf('at s = (%g,%g)', s(1), s(2)));
            end
        end

        function aWedgeWithAGeneralQuadraticMatchesTheDefinitionSup(testCase)
        % The general unbounded case: a non-axis-aligned wedge, a quadratic with cross and linear
        % terms, and a nonzero constant. The reference samples the wedge out to a radius where the
        % quadratic already dominates, which for a PD A is where the maximiser must be.
            A = [3 -1; -1 2]; L = [0.5; -2]; c = 1.25;
            W = [1 1]; dF = [1 2]; dL = [2 -1];
            g = conjConvexPolygon(W, dF, dL, A, L, c);
            S = conjConvexPolygonTest.probes();
            for i = 1:size(S,1)
                brute = conjConvexPolygonTest.supOverWedge(S(i,:), W, dF, dL, A, L, c);
                got = g.eval(S(i,:));
                testCase.verifyEqual(got, brute, 'AbsTol', 1e-8, ...
                    sprintf('at s=(%g,%g)', S(i,1), S(i,2)));
            end
        end

        % ---- what it refuses, and refuses by name ------------------------------------------------

        function anIndefiniteOrSemidefiniteQuadraticIsRefusedByName(testCase)
            W = [0 0; 1 0; 0 1];
            testCase.verifyError(@() conjConvexPolygon(W, [], [], [0 1; 1 0], [0;0], 0), ...
                'conjConvexPolygon:notPositiveDefinite');
            testCase.verifyError(@() conjConvexPolygon(W, [], [], [1 0; 0 0], [0;0], 0), ...
                'conjConvexPolygon:notPositiveDefinite', ...
                'a semidefinite A can have an infinite sup in a recession direction');
        end

        function aClockwiseOrNonConvexPolygonIsRefusedRatherThanMirrored(testCase)
        % An orientation slip produces a perfectly well-formed mesh built on the wrong side, whose
        % values are wrong on a whole region and right at the dual points. Refusing is the only
        % safe answer.
            testCase.verifyError(@() conjConvexPolygon([0 0; 0 1; 1 0], [], [], eye(2), [0;0], 0), ...
                'conjConvexPolygon:orientation', 'clockwise');
            testCase.verifyError(@() conjConvexPolygon([0 0; 2 0; 1 0.1; 2 2; 0 2], [], [], ...
                                                        eye(2), [0;0], 0), ...
                'conjConvexPolygon:orientation', 'a reflex vertex');
        end

        function oneRayIsRefused(testCase)
            testCase.verifyError(@() conjConvexPolygon([0 0; 1 0], [0 1], [], eye(2), [0;0], 0), ...
                'conjConvexPolygon:oneRay');
        end

        function agreesAcrossRandomPolygonsAndCoordinateScales(testCase)
        % This is the SCIP-relevant fast path (the box envelope calls it directly), and it is a
        % DIFFERENT implementation from conjConvexOverPiece.m -- the legacy Case-C routine whose
        % vertex-cone construction had a real coverage-gap bug found and fixed 2026-08-30 (a
        % numeric feasibility-probe step too small relative to its own tolerance, on a
        % coordinate-scale this file's own fixtures never reached). That bug does not transfer
        % here mechanically (this file builds cones analytically via outwardNormal, not by a
        % numeric probe), but nothing in the existing fixtures tests scales bigger than a handful
        % of units either. This sweeps random convex polygons (via convhull, so genuinely convex)
        % and random SPD Hessians across scales from 0.1 to 1000 -- the range real coefficient/box
        % data can plausibly span -- against exactPolygonOracle: for a BOUNDED convex polygon and
        % a genuinely convex q the maximiser is exactly a vertex, an edge's clamped 1-D stationary
        % point, or the interior unconstrained one, so there is a CLOSED FORM and no need for
        % supOverPolygon's coarser sample-then-search reference. Measured directly: at one
        % elongated, ill-conditioned trial this test's first version (using supOverPolygon)
        % reported a spurious ~1e-3 disagreement that conjConvexPolygon's OWN answer disproved --
        % it matched the exact oracle to ~1e-12 at every one of those points, so the mismatch was
        % supOverPolygon's own search stalling early on that geometry, not a defect here.
            rng(20260830);
            nTrials = 40;
            for trial = 1:nTrials
                pts = (rand(18, 2) - 0.5) * (10^(rand*4-1)) + (rand(1,2)-0.5)*10;
                k = convhull(pts(:,1), pts(:,2));
                W = pts(k(1:end-1), :);
                x = W(:,1); y = W(:,2); n = size(W,1); a = 0;
                for i = 1:n, j = mod(i,n)+1; a = a + (x(i)*y(j) - x(j)*y(i)); end
                if a < 0, W = flipud(W); end
                scale = max(1, max(abs(W(:))));
                M = randn(2) * scale; A = M*M'/scale + eye(2)*1e-3*scale;
                L = (rand(2,1)-0.5) * scale;
                c = (rand-0.5) * scale;

                g = conjConvexPolygon(W, [], [], A, L, c);
                for s = [conjConvexPolygonTest.probes(); (rand(4,2)-0.5)*scale*3]'
                    got = g.eval(s');
                    ref = conjConvexPolygonTest.exactPolygonOracle(s', W, A, L, c);
                    testCase.verifyEqual(got, ref, 'RelTol', 1e-9, 'AbsTol', 1e-6, ...
                        sprintf('trial %d, scale %.3g, s=(%g,%g)', trial, scale, s(1), s(2)));
                end
            end
        end
    end

    methods (Static)

        function S = probes()
        % Probe points spread over all four quadrants and across several magnitudes, so that every
        % cell of the dual subdivision is visited by something.
            S = [ 0 0; 2 2; -1 3; 0.3 0.4; 5 -2; -4 -4; 1.5 0.5; -3 1; 1 -6; 8 8; -0.2 -0.9; 0 4 ];
        end

        function v = supOverPolygon(s, W, A, L, c)
        % objective: sup_{x in P} <s,x> - q(x) over the convex hull of W, to full double accuracy.
        %   A coarse scan locates the maximiser's cell, then a pattern search refines it. The two
        %   stages are separated on purpose: a scan alone is only as good as its step, and a
        %   tolerance chosen to absorb that error is a tolerance on the CODE UNDER TEST as well.
        %   Measured on the wedge fixture: the 300x300 scan alone was off by up to 3.2e-3, which
        %   was the whole of the apparent disagreement.
            n = size(W,1);
            cand = W;                                        % the vertices, always -- a sup that a
            for i = 1:n                                      % vertex attains is invisible to a grid
                j = mod(i, n) + 1;
                t = linspace(0, 1, 200).';
                cand = [cand; (1-t)*W(i,:) + t*W(j,:)]; %#ok<AGROW>
            end
            for a = linspace(0, 1, 25)
                for i = 2:n-1
                    b = linspace(0, 1, 25).';
                    cand = [cand; W(1,:) + a*(W(i,:) - W(1,:)) + a*b*(W(i+1,:) - W(i,:))]; %#ok<AGROW>
                end
            end
            vals = zeros(size(cand,1),1);
            for i = 1:numel(vals), vals(i) = conjConvexPolygonTest.obj(s, cand(i,:), A, L, c); end
            [v, bi] = max(vals);
            % The direction set must include the polygon's own EDGE directions. A maximiser that
            % lies ON an edge can only be improved by SLIDING along it, and a fixed axis/diagonal
            % set has no such direction in general -- the search then stalls at the scan's own
            % accuracy while reporting convergence. Measured on the wedge fixture below: the
            % maximiser sat exactly on a boundary ray and the refined value came out WORSE than the
            % coarse scan's.
            D = conjConvexPolygonTest.baseDirections();
            n2 = size(W,1);
            for i = 1:n2
                e = W(mod(i,n2)+1,:) - W(i,:);
                D = [D; e/norm(e); -e/norm(e)]; %#ok<AGROW>
            end
            inside = @(p) conjConvexPolygonTest.inHull(p, W);
            v = conjConvexPolygonTest.patternSearch(s, cand(bi,:), inside, A, L, c, v, D);
        end

        function v = supOverWedge(s, w, dF, dL, A, L, c)
        % objective: the same, over the wedge w + cone(dF, dL). A is PD so the objective is
        %   strictly concave and its maximiser is within a bounded radius; the scan runs past it
        %   and the pattern search finishes the job.
            dF = dF / norm(dF); dL = dL / norm(dL);
            [Ag, Bg] = ndgrid(linspace(0, 40, 120), linspace(0, 40, 120));
            cand = w + Ag(:)*dF + Bg(:)*dL;
            vals = zeros(size(cand,1),1);
            for i = 1:numel(vals), vals(i) = conjConvexPolygonTest.obj(s, cand(i,:), A, L, c); end
            [v, bi] = max(vals);
            % Same reason as supOverPolygon: the maximiser can sit exactly on one of the two
            % boundary rays, so those two directions must be in the set.
            D = [conjConvexPolygonTest.baseDirections(); dF; -dF; dL; -dL];
            M = [dF(:), dL(:)];
            inside = @(p) all(M \ (p(:) - w(:)) >= -1e-12);
            v = conjConvexPolygonTest.patternSearch(s, cand(bi,:), inside, A, L, c, v, D);
        end

        function D = baseDirections()
            D = [1 0; -1 0; 0 1; 0 -1; 1 1; 1 -1; -1 1; -1 -1];
            D = D ./ vecnorm(D, 2, 2);
        end

        function v = patternSearch(s, p, inside, A, L, c, v, D)
        % objective: refine a maximiser by a shrinking pattern search over the direction set D,
        %   staying feasible. Deliberately a GENERIC optimiser rather than the closed form: a
        %   reference that reuses the construction under test checks nothing.
            step = 1;
            for it = 1:80
                improved = false;
                for k = 1:size(D,1)
                    qp = p + step * D(k,:);
                    if inside(qp)
                        vq = conjConvexPolygonTest.obj(s, qp, A, L, c);
                        if vq > v, v = vq; p = qp; improved = true; end
                    end
                end
                if ~improved
                    step = step / 2;
                    if step < 1e-14, break, end
                end
            end
        end

        function tf = inHull(p, W)
        % objective: is p in the convex hull of the CCW polygon W. Half-plane test, one per edge.
            n = size(W,1); tf = true;
            tol = 1e-12 * (1 + max(abs(W(:))));
            for i = 1:n
                j = mod(i, n) + 1;
                e = W(j,:) - W(i,:);
                if (p(:).' - W(i,:)) * [e(2); -e(1)] > tol, tf = false; return, end
            end
        end

        function v = obj(s, x, A, L, c)
            x = x(:);
            v = s(:).'*x - (0.5*(x.'*A*x) + L(:).'*x + c);
        end

        function best = exactPolygonOracle(s0, W, Q, L, c)
        % EXACT sup_{x in P} <s,x>-q(x) for q convex on a BOUNDED convex polygon P: since q is
        % strictly convex on a compact set, the maximiser is exactly a vertex, an edge's own
        % clamped 1-D stationary point, or the unconstrained interior stationary point -- no
        % sampling, no search, full double precision regardless of scale or conditioning.
            s0 = s0(:);
            n = size(W,1);
            cand = W;
            if rcond(Q) > 1e-14
                xs = Q \ (s0 - L(:));
                if conjConvexPolygonTest.inHull(xs.', W)
                    cand = [cand; xs.'];
                end
            end
            for e = 1:n
                a = W(e,:).'; b = W(mod(e,n)+1,:).';
                d = b - a;
                den = d.'*Q*d;
                if abs(den) > 1e-14
                    t = min(1, max(0, (d.'*s0 - d.'*Q*a - L(:).'*d) / den));
                    cand = [cand; (a + t*d).']; %#ok<AGROW>
                end
            end
            best = -inf;
            for k = 1:size(cand,1)
                best = max(best, conjConvexPolygonTest.obj(s0, cand(k,:), Q, L, c));
            end
        end
    end
end
