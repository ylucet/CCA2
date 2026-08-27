classdef unboundedFaceTest < matlab.unittest.TestCase
    % Step 1 and Step 2 over an UNBOUNDED face -- next-step 1 of the session handoff.
    %
    % The pipeline used to reject such a face outright (quaPolToPlq:unboundedFace) because both
    % failure modes behind the gate were silent: quaPolToPlq threw the ray DIRECTION away, and
    % plq_1p read region's +/-intmax direction markers as ordinary coordinates, returning
    % envelopes and conjugates carrying 2147483647 and intmax^2 = 4611686014132420609 (max
    % error 1.15e18 as measured in the handoff). These tests pin the three new pieces --
    % fanUnboundedFace, convEnvUnbounded, and conjugateFunction's envelope-keyed dispatch --
    % and, importantly, pin the REFUSALS too, since the affine-envelope class is not all of them.

    methods (Static)
        function v = evalAnyConj(g, x)
        % Evaluate a conj result whatever representation it came back as -- the same idiom as
        % conjCPLQTest.evalConjResult and biconjCPLQTest.evalAnyResult. Since the numeric routes
        % were widened, one input family can legitimately produce a QuaParCPLQ or a mesh, and a
        % test about VALUES must not hard-code either. Note the two disagree on how "uncovered"
        % reads: NaN from the symbolic form, +inf from a mesh.
            if isa(g, 'QuaParCPLQ')
                v = evalFunctionNDomain(g.fnd, x);
            else
                v = g.eval(x);
            end
        end

        function d = faceOf(ineqs, x, y)
            d = domain();
            d = d.domainEdge(ineqs, [x, y]);
        end
        function m = inFace(r, P)
            [A, b, lin] = r.linearForm;
            m = all(A(lin,:)*P' <= b(lin) + 1e-9, 1)';
        end
    end

    methods (Test)

        function fanCoversTheFaceExactlyAndOverclaimsNothing(testCase)
        % The fan must be a COVER of the face (a sup over a union is the max of the sups, so a
        % cover is enough) made of SUBSETS of it (a superset would inflate its own sup and hence
        % the max). Both halves are checked on a grid, since either one failing is silent.
            x = sym('x'); y = sym('y');
            faces = { [x, -y], [y, -x-1, x-1], [-x, -y, 1-x-y], ...
                      [-y, -x-y, x-y-2], [-y, -x-y, x-y-2, x-3] };
            [gx, gy] = meshgrid(linspace(-9, 9, 61), linspace(-9, 9, 61));
            P = [gx(:), gy(:)];
            for k = 1:numel(faces)
                r  = region(faces{k}, [x y]);
                ds = fanUnboundedFace(r, x, y);
                inR = unboundedFaceTest.inFace(r, P);
                covered = false(size(P,1),1);
                for i = 1:numel(ds)
                    s = unboundedFaceTest.inFace(ds{i}.polygon, P);
                    testCase.verifyEqual(sum(s & ~inR), 0, ...
                        sprintf('face %d piece %d claims points outside the face', k, i));
                    covered = covered | s;
                end
                testCase.verifyEqual(sum(inR & ~covered), 0, ...
                    sprintf('face %d is not fully covered by its fan', k));
            end
        end

        function affineEnvelopeOverAWedgeIsTheTangentPlaneAtTheApex(testCase)
            x = sym('x'); y = sym('y');
            % co(x*y) over the first quadrant is 0: the tangent plane at the apex.
            r = region([-x, -y], [x y]);
            [e, why] = convEnvUnbounded(r, x*y, [x y]);
            testCase.verifyEqual(why, 'wedge');
            testCase.verifyTrue(isAlways(e == 0));

            % Shifting the apex and scaling: q = 2xy has grad (2y,2x) = (4,2) at (1,2), q = 4.
            r = region([-(x-1), -(y-2)], [x y]);
            e = convEnvUnbounded(r, 2*x*y, [x y]);
            testCase.verifyTrue(isAlways(simplify(e - (4*x + 2*y - 4)) == 0));
        end

        function affineEnvelopeOverAHalfStripUsesTheSmallerRaySlope(testCase)
        % <a,d> = min over the two base vertices of the directional derivative along d. Pinned
        % on a case where the two DISAGREE (w'Qd ~= 0), which is what the min is there for.
            x = sym('x'); y = sym('y');
            % {0<=x<=1, y>=0}, q = x*y: grad q = (y,x), so along d = (0,1) the slope is x, i.e.
            % 0 at (0,0) and 1 at (1,0). The min is 0, and the envelope is 0.
            r = region([-y, -x, x-1], [x y]);
            [e, why] = convEnvUnbounded(r, x*y, [x y]);
            testCase.verifyEqual(why, 'half-strip');
            testCase.verifyTrue(isAlways(e == 0));

            % Same strip, q = -x^2: concave along the base, flat along the rays -> the chord -x.
            e = convEnvUnbounded(r, -x^2, [x y]);
            testCase.verifyTrue(isAlways(simplify(e + x) == 0));
        end

        function theResultFollowsTheFaceFunctionNotAlwaysThatOfXY(testCase)
        % REGRESSION for the defect that motivated the whole classify-by-eigenvalues change.
        % plq_1p.convexEnvelope1 was a closed form in the vertex COORDINATES alone and never
        % referenced obj.f, so it computed co(x*y) whatever the piece carried -- pinned by
        % measurement, it returned a byte-identical envelope for x*y, x^2-y^2, (x^2+y^2)/2 and
        % 3xy+7x-2y+5.
        %
        % The invariant checked here is the one that must hold whatever the internals do:
        % adding an affine l to f shifts the conjugate by the exact rule
        %       (f + l)*(s) = f*(s - a) - c,     l(x) = <a,x> + c.
        % Testing the ENVELOPE expression directly is no longer meaningful, because an indefinite
        % q with an affine part is now handled by a frame change (xyFrame) that deliberately
        % leaves obj.envelope empty in the original frame.
            x = sym('x'); y = sym('y');
            t = [0 0; 2 0; 0 2];
            a = [3; -2]; c = 7;
            base  = plq(plq_1p(domain(t,x,y), symbolicFunction(x*y)));
            shift = plq(plq_1p(domain(t,x,y), symbolicFunction(x*y + a(1)*x + a(2)*y + c)));
            base  = base.triangulate;  base  = base.maximum;
            shift = shift.triangulate; shift = shift.maximum;
            S = [0.5 0.5; 1 -1; -1 2; 2 2; 0 0];
            for k = 1:size(S,1)
                sh = S(k,:);
                got  = evalFunctionNDomain(shift.maxConjugate, sh);
                want = evalFunctionNDomain(base.maxConjugate,  sh - a') - c;
                testCase.verifyFalse(isnan(got), sprintf('(f+l)* uncovered at (%g,%g)', sh));
                testCase.verifyEqual(got, want, 'AbsTol', 1e-9, sprintf( ...
                    '(f+l)*(s) must equal f*(s-a)-c at s=(%g,%g)', sh(1), sh(2)));
            end
        end

        function conjugateOfAnAffineEnvelopeOverAWedgeIsExact(testCase)
        % End to end: an unbounded piece through triangulate -> maximum, against the closed form.
        % f = x*y on the first quadrant has convex envelope 0, so f* is the indicator of the
        % third quadrant -- 0 where every component of s is <= 0, +inf elsewhere.
        %
        % The (-10,5) and (-1,0.5) probes are the point of this test. An earlier version of the
        % dual construction used region.getNormalConeVertex, which reads the VERTEX LIST, and
        % getVertices fabricates corner vertices at (+/-intmax,+/-intmax) for an unbounded
        % region. The corner (intmax,intmax) supplied the direction (1,1), so the cone came out
        % {s1+s2 <= 0, s1 <= 0} instead of {s1 <= 0, s2 <= 0} -- a different set, agreeing with
        % the truth on the obvious probes and reporting 0 instead of +inf on these two.
            x = sym('x'); y = sym('y');
            d = unboundedFaceTest.faceOf([-x, -y], x, y);
            P = plq(plq_1p(d, symbolicFunction(x*y)));
            P = P.triangulate;
            P = P.maximum;
            for s = [-1 -1; -4 -2; 0 0]'
                testCase.verifyEqual(evalFunctionNDomain(P.maxConjugate, s'), 0, 'AbsTol', 1e-12);
            end
            for s = [0.5 -2; 1.5 1.5; -0.5 3; -10 5; -1 0.5; 5 -10]'
                testCase.verifyTrue(isnan(evalFunctionNDomain(P.maxConjugate, s')), ...
                    sprintf('f* must be +inf at (%g,%g), off the third quadrant', s(1), s(2)));
            end
        end

        function unboundedMultiFaceQuaPolConjugateIsExact(testCase)
        % The whole point of next-step 1: an UNBOUNDED multi-face QuaPol through the public
        % conj('cplq') entry. The 4-cone fan with f = |x|+|y|, whose conjugate is the indicator
        % of [-1,1]^2. This used to be a hard rejection (quaPolToPlq:unboundedFace).
            V = [0 0; -1 0; 0 1; 1 0; 0 -1];
            E = [1 2 0; 1 3 0; 1 4 0; 1 5 0];
            F = [3 2; 2 1; 1 4; 4 3];
            f = [0 0 0 1 1 0; 0 0 0 -1 1 0; 0 0 0 -1 -1 0; 0 0 0 1 -1 0];
            g = QuaPol(V, E, f, F).conj('cplq');
            % ACCESSOR CHANGED 2026-08-26, not the assertion. This read `g.fnd` and tested the
            % outside with `isnan`, both of which assume a QuaParCPLQ. f = |x|+|y| is ALL-AFFINE
            % over four cones, so it now takes the `conjAffinePLQ` route (TODO.md G2) and comes
            % back as a QuaPol MESH in 0.95 s instead of going through the symbolic pipeline --
            % where "uncovered" reads as +inf rather than NaN.
            %
            % The values are unchanged and were re-measured before this edit: 0 at all six probe
            % points inside (and at the corners (-1,1) and (1,-1) too), +inf at all four outside
            % plus (-3,-3). This test's subject is the VALUE -- "the conjugate is the indicator of
            % [-1,1]^2" -- so it must read whichever representation it is handed.
            ev = @(x) unboundedFaceTest.evalAnyConj(g, x);
            inBox  = [0 0; 0.5 0.5; 1 1; -1 -1; 0.9 -0.3; -0.25 0.75];
            outBox = [1.5 0; 0 2; 2 2; -1.5 0.5];
            for i = 1:size(inBox,1)
                testCase.verifyEqual(ev(inBox(i,:)), 0, 'AbsTol', 1e-12, ...
                    sprintf('f* must be 0 inside [-1,1]^2 at (%g,%g)', inBox(i,1), inBox(i,2)));
            end
            for i = 1:size(outBox,1)
                v = ev(outBox(i,:));
                testCase.verifyTrue(isnan(v) || isinf(v), sprintf( ...
                    'f* must be +inf outside [-1,1]^2 at (%g,%g), got %g', ...
                    outBox(i,1), outBox(i,2), v));
            end
        end

        function conjugateOfAnAffineEnvelopeOverAHalfStripIsExact(testCase)
        % q = -x^2 on {0<=x<=1, y>=0} has envelope -x, so
        %   f*(s) = sup{ s1 x + s2 y + x : 0<=x<=1, y>=0 } = max(0, s1+1) when s2 <= 0, +inf else.
            x = sym('x'); y = sym('y');
            d = unboundedFaceTest.faceOf([-y, -x, x-1], x, y);
            P = plq(plq_1p(d, symbolicFunction(-x^2)));
            P = P.triangulate;
            P = P.maximum;
            probes = [-1 -1, 0; 0.5 -2, 1.5; 0 0, 1; 2 -0.5, 3; -4 -4, 0];
            for i = 1:size(probes,1)
                got = evalFunctionNDomain(P.maxConjugate, probes(i,1:2));
                testCase.verifyEqual(got, probes(i,3), 'AbsTol', 1e-12);
            end
        end

        function minusInfinityIsReportedNotApproximated(testCase)
        % co q = -inf makes f* identically +inf, which is a different kind of answer than a
        % function. It must never come back as a finite-looking affine envelope.
            x = sym('x'); y = sym('y');
            % x*y on the SECOND quadrant: (-1,1) is a recession direction with d'Qd = -2 < 0.
            r = region([x, -y], [x y]);
            testCase.verifyError(@() convEnvUnbounded(r, x*y, [x y]), ...
                'convEnvUnbounded:minusInfinity');
            % A d'Qd == 0 direction with a negative linear slope is just as fatal: on the first
            % quadrant, q = xy + 3x - 2y + 7 decreases along +y at x = 0.
            r = region([-x, -y], [x y]);
            testCase.verifyError(@() convEnvUnbounded(r, x*y + 3*x - 2*y + 7, [x y]), ...
                'convEnvUnbounded:minusInfinity');
        end

        function convexQuadraticOverAWedgeIsExact(testCase)
        % The handoff's headline case. q = (x^2+y^2)/2 over the wedge {x<=0, y>=0} is CONVEX, so
        % co q = q and Step 1 has nothing to do -- but Step 2 must then conjugate a CURVED
        % function, which cPLQ's support-function branch cannot. conjConvexOverPiece does it by
        % the KKT active set: apex cell, one cell per ray, and the interior cell. The closed form
        % is min(s1,0)^2/2 + max(s2,0)^2/2.
            x = sym('x'); y = sym('y');
            s1 = sym('s_1'); s2 = sym('s_2');
            r = region([x, -y], [x y]);
            cj = conjConvexOverPiece(r, eye(2), [0;0], 0, [s1 s2]);
            testCase.verifyEqual(numel(cj), 4, 'apex + 2 rays + interior');
            S = [0 0; 1 1; -1 -1; 2 -0.5; -0.5 2; 3 1; -2 -3; 0.5 0.5; -3 2; 1.5 -2.5];
            for t = 1:size(S,1)
                ref = min(S(t,1),0)^2/2 + max(S(t,2),0)^2/2;
                got = evalFunctionNDomain(cj, S(t,:));
                testCase.verifyEqual(got, ref, 'AbsTol', 1e-12, ...
                    sprintf('at s=(%g,%g)', S(t,1), S(t,2)));
            end
        end

        function convexQuadraticOverATriangleCoversTheInterior(testCase)
        % cPLQ's conjugateOfPiecePoly returns only the 3 VERTEX cells for a strictly convex q on
        % a bounded triangle, leaving every interior dual point uncovered -- by design there,
        % since its Step 1 never hands Step 2 a curved envelope. The active-set construction
        % supplies the edge and interior cells too, so the dual plane is actually covered.
            x = sym('x'); y = sym('y');
            s1 = sym('s_1'); s2 = sym('s_2');
            r = region([-x, -y, x+y-2], [x y]);
            cj = conjConvexOverPiece(r, eye(2), [0;0], 0, [s1 s2]);
            testCase.verifyGreaterThan(numel(cj), 3, 'vertex cells alone do not cover the dual');
            % (0.5,0.5) is interior: the unconstrained maximizer x* = s lies inside the triangle,
            % so f*(s) = |s|^2/2 -- a point no vertex cell can be right at.
            testCase.verifyEqual(evalFunctionNDomain(cj, [0.5 0.5]), 0.25, 'AbsTol', 1e-12);
            testCase.verifyEqual(evalFunctionNDomain(cj, [-1 -1]),   0,    'AbsTol', 1e-12);
            testCase.verifyEqual(evalFunctionNDomain(cj, [3 1]),     4,    'AbsTol', 1e-12);
        end

        function convexFaceOverAWedgeGoesThroughTheWholePipeline(testCase)
        % End to end through triangulate -> maximum for a CONVEX face on an unbounded domain.
        % This used to raise plq_1p:conjugateFunction:unboundedNonAffine: Step 1 correctly
        % returned q as its own envelope, and Step 2 had no branch for a curved one. It now
        % goes through conjConvexOverPiece, and the answer is the closed form
        % min(s1,0)^2/2 + max(s2,0)^2/2.
            x = sym('x'); y = sym('y');
            d = unboundedFaceTest.faceOf([x, -y], x, y);
            P = plq(plq_1p(d, symbolicFunction((x^2 + y^2)/2)));
            P = P.triangulate;
            P = P.maximum;
            S = [0 0; 1 1; -1 -1; 2 -0.5; -0.5 2; -3 2; 1.5 -2.5];
            for t = 1:size(S,1)
                ref = min(S(t,1),0)^2/2 + max(S(t,2),0)^2/2;
                got = evalFunctionNDomain(P.maxConjugate, S(t,:));
                testCase.verifyEqual(got, ref, 'AbsTol', 1e-10, ...
                    sprintf('at s=(%g,%g)', S(t,1), S(t,2)));
            end
        end

        function quaPolToPlqKeepsTheRayDirection(testCase)
        % Defect (1) of quaPolToPlq's header: faceVertexIndices reads E(j,1) for every edge, but
        % for a RAY that is the base point and E(j,2) is a DIRECTION point. Two rays off a shared
        % apex therefore both reported the apex, and the bounded domain() constructor turned that
        % into a degenerate NaN <= 0. The half-plane route keeps both.
            V = [0 0; -1 0; 0 1; 1 0; 0 -1];
            E = [1 2 0; 1 3 0; 1 4 0; 1 5 0];
            F = [3 2; 2 1; 1 4; 4 3];
            f = [0 0 0 1 1 0; 0 0 0 -1 1 0; 0 0 0 -1 -1 0; 0 0 0 1 -1 0];   % |x| + |y|
            q = QuaPol(V, E, f, F);
            testCase.verifyFalse(q.isDomBounded);
            p = quaPolToPlq(q);
            testCase.verifyEqual(p.nPieces, 4);
            for k = 1:4
                r = p.pieces(k).d.polygon;
                testCase.verifyFalse(any(isnan(double(r.vx))) || any(isnan(double(r.vy))), ...
                    sprintf('face %d produced NaN vertices', k));
                [~, kind] = r.recessionRays;
                testCase.verifyEqual(kind, 'wedge', ...
                    sprintf('face %d should be a pointed unbounded cone', k));
            end
        end

        % ==========================================================================================
        % NON-POINTED FACES. A face whose recession cone contains a LINE -- a half-plane or a slab
        % -- has no apex to fan from, so fanUnboundedFace refuses it outright. It does not need a
        % new envelope theory: cutting the face along a line transverse to its lineality direction
        % leaves two pieces whose recession cones are pointed, and the fan already handles those.
        % A cut is a cover, so the conjugate decomposes over it exactly as it does over the fan.
        % ==========================================================================================

        function fanCoversANonPointedFace(testCase)
        % Same contract as fanCoversTheFaceExactlyAndOverclaimsNothing (subset + cover, on a grid),
        % now for the two non-pointed shapes: a HALF-PLANE and a SLAB.
            x = sym('x'); y = sym('y');
            faces = { [-y], [-y, y-1], [x-2], [-x-y, x+y-3] };   %#ok<NBRAK>
            [gx, gy] = meshgrid(linspace(-9, 9, 61), linspace(-9, 9, 61));
            P = [gx(:), gy(:)];
            for k = 1:numel(faces)
                r  = region(faces{k}, [x y]);
                ds = fanUnboundedFace(r, x, y);
                inR = unboundedFaceTest.inFace(r, P);
                covered = false(size(P,1),1);
                for i = 1:numel(ds)
                    s = unboundedFaceTest.inFace(ds{i}.polygon, P);
                    testCase.verifyEqual(sum(s & ~inR), 0, ...
                        sprintf('face %d piece %d claims points outside the face', k, i));
                    covered = covered | s;
                end
                testCase.verifyEqual(sum(inR & ~covered), 0, ...
                    sprintf('face %d is not fully covered by its fan', k));
            end
        end

        function everyPieceOfANonPointedFanIsPointed(testCase)
        % The reduction's whole point: after the cut, every emitted piece is a shape the existing
        % machinery already handles. Pinned separately from the cover check because a fan could
        % cover correctly and still hand convEnvUnbounded something it must refuse.
            x = sym('x'); y = sym('y');
            faces = { [-y], [-y, y-1], [x-2], [-x-y, x+y-3] };   %#ok<NBRAK>
            for k = 1:numel(faces)
                ds = fanUnboundedFace(region(faces{k}, [x y]), x, y);
                testCase.verifyGreaterThanOrEqual(numel(ds), 2, ...
                    sprintf('face %d: a non-pointed face must be cut, not passed through', k));
                for i = 1:numel(ds)
                    [~, kind] = ds{i}.polygon.recessionRays;
                    testCase.verifyTrue(ismember(kind, {'bounded','ray','wedge'}), ...
                        sprintf('face %d piece %d is still %s', k, i, kind));
                end
            end
        end

        function convexQuadraticOverAHalfPlaneIsExact(testCase)
        % End to end on a HALF-PLANE. q = (x^2+y^2)/2 over {y >= 0} is convex and separable, so
        %       f*(s) = sup_x (s1 x - x^2/2) + sup_{y>=0} (s2 y - y^2/2)
        %             = s1^2/2 + max(s2,0)^2/2.
            x = sym('x'); y = sym('y');
            d = unboundedFaceTest.faceOf([-y], x, y);   %#ok<NBRAK>
            P = plq(plq_1p(d, symbolicFunction((x^2 + y^2)/2)));
            P = P.triangulate;
            P = P.maximum;
            S = [0 0; 1 1; -1 -1; 2 -0.5; -0.5 2; -3 2; 1.5 -2.5];
            for t = 1:size(S,1)
                ref = S(t,1)^2/2 + max(S(t,2),0)^2/2;
                got = evalFunctionNDomain(P.maxConjugate, S(t,:));
                testCase.verifyEqual(got, ref, 'AbsTol', 1e-10, ...
                    sprintf('at s=(%g,%g)', S(t,1), S(t,2)));
            end
        end

        function convexQuadraticOverASlabIsExact(testCase)
        % End to end on a SLAB {0 <= y <= 1}, whose recession cone is a LINE (both +y and -y are
        % excluded but +/-x recede), the second non-pointed shape. Same q, so
        %       f*(s) = s1^2/2 + h(s2),  h = 0 (s2<=0) | s2^2/2 (0<=s2<=1) | s2-1/2 (s2>=1).
            x = sym('x'); y = sym('y');
            d = unboundedFaceTest.faceOf([-y, y-1], x, y);
            P = plq(plq_1p(d, symbolicFunction((x^2 + y^2)/2)));
            P = P.triangulate;
            P = P.maximum;
            S = [0 0; 1 1; -1 -1; 2 -0.5; -0.5 2; 0.5 0.5; 3 4];
            for t = 1:size(S,1)
                s2 = S(t,2);
                if     s2 <= 0, h = 0;
                elseif s2 <= 1, h = s2^2/2;
                else,           h = s2 - 0.5;
                end
                ref = S(t,1)^2/2 + h;
                got = evalFunctionNDomain(P.maxConjugate, S(t,:));
                testCase.verifyEqual(got, ref, 'AbsTol', 1e-10, ...
                    sprintf('at s=(%g,%g)', S(t,1), S(t,2)));
            end
        end

        % ==========================================================================================
        % NONCONVEX q WHOSE ENVELOPE OVER AN UNBOUNDED PIECE IS CURVED. convEnvUnbounded only ever
        % returns an AFFINE envelope, and refuses (convEnvUnbounded:convexAlongRay) as soon as some
        % ray direction has positive curvature -- because then the envelope bends along that ray and
        % no affine function can be it. That refusal is correct about ITS OWN formula and wrong as a
        % verdict on the case: the envelope exists, it is simply not affine.
        % ==========================================================================================

        function nonconvexQuadraticWithACurvedEnvelopeOverAHalfStripIsExact(testCase)
        % q = -x^2 + y^2 over the half-strip {0 <= x <= 1, y >= 0}. Indefinite, so Step 1 must
        % produce an envelope; the domain is a product and q is separable, so
        %       co q = co_x(-x^2) + co_y(y^2) = -x + y^2,
        % convex but NOT affine -- exactly the shape convEnvUnbounded refuses. Ground truth:
        %       f*(s) = sup_{0<=x<=1}(s1 x + x^2) + sup_{y>=0}(s2 y - y^2)
        %             = max(0, s1+1) + max(s2,0)^2/4.
            x = sym('x'); y = sym('y');
            d = unboundedFaceTest.faceOf([-x, x-1, -y], x, y);
            P = plq(plq_1p(d, symbolicFunction(-x^2 + y^2)));
            P = P.triangulate;
            P = P.maximum;
            S = [0 0; 1 1; -1 -1; -2 0.5; 0.5 2; -3 -2; 2 3];
            for t = 1:size(S,1)
                ref = max(0, S(t,1)+1) + max(S(t,2),0)^2/4;
                got = evalFunctionNDomain(P.maxConjugate, S(t,:));
                testCase.verifyEqual(got, ref, 'AbsTol', 1e-10, ...
                    sprintf('at s=(%g,%g)', S(t,1), S(t,2)));
            end
        end

        function curvedEnvelopeOverAWedgeIsExact(testCase)
        % The same gap on a WEDGE rather than a half-strip. q = x*y over K = {0 <= y <= x}, whose
        % extreme rays are (1,0) with d'Qd = 0 and (1,1) with d'Qd = 2 > 0 -- so convEnvUnbounded
        % refuses it as convexAlongRay, and that refusal is right about its OWN formula: the apex
        % tangent plane is not the envelope here. Checked by hand at h = (1,1), where the best
        % affine minorant is worth 1 (take l = 2*sqrt(e)*y - e at e = 1) and the tangent plane is
        % worth 0, so co is strictly above the tangent plane and cannot be affine.
        %
        % Ground truth, parametrizing y = mu*x on mu in [0,1], x >= 0:
        %       sup_x [ x*(s1 + s2*mu) - mu*x^2 ] = (s1+s2*mu)^2/(4*mu)  when s1 + s2*mu > 0,
        % which for s1 <= 0 increases in mu, so the max sits at mu = 1; and mu = 0 makes the sup
        % +inf as soon as s1 > 0. Hence
        %       f*(s) = +inf                 if s1 > 0
        %             = (s1+s2)^2/4          if s1 <= 0 and s1 + s2 > 0
        %             = 0                    otherwise.
        %
        % NOTE for whoever implements this: it is a DERIVATION, not a wiring gap. Writing the
        % wedge as v + alpha*d1 + beta*d2 turns the conjugate into a QP over {alpha,beta >= 0}
        % whose quadratic form has determinant -(d11*d22 - d12*d21)^2/4 <= 0 -- never positive
        % definite, but nonnegative ON THE CONE by AM-GM whenever both rays lie in one closed
        % quadrant of the bilinear frame. So it is a convex QP on the cone with a SINGULAR form,
        % which is a different active-set analysis from conjConvexOverPiece's.
            x = sym('x'); y = sym('y');
            d = unboundedFaceTest.faceOf([-y, y-x], x, y);
            P = plq(plq_1p(d, symbolicFunction(x*y)));
            P = P.triangulate;
            P = P.maximum;
            finite = [-1 3; -2 0.5; -1 -1; -3 1; 0 -2];
            for t = 1:size(finite,1)
                s = finite(t,:);
                if s(1) + s(2) > 0, ref = (s(1)+s(2))^2/4; else, ref = 0; end
                got = evalFunctionNDomain(P.maxConjugate, s);
                testCase.verifyEqual(got, ref, 'AbsTol', 1e-10, ...
                    sprintf('at s=(%g,%g)', s(1), s(2)));
            end
            for s = [1 0; 2 -3; 0.5 0.5]'
                testCase.verifyTrue(isnan(evalFunctionNDomain(P.maxConjugate, s')), ...
                    sprintf('f* must be +inf at (%g,%g)', s(1), s(2)));
            end
        end
    end
end
