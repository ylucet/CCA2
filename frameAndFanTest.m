classdef frameAndFanTest < matlab.unittest.TestCase
% Unit tests for the four free functions Case C's unbounded route is built from, each of which is
% exercised today only from inside a whole symbolic pipeline run:
%
%   xyFrame            the congruence that turns an indefinite quadratic into exactly x*y
%   transformDomain    the same change of variables applied to a domain
%   fanUnboundedFace   the cover of an unbounded face by triangles, half-strips and wedges
%   convEnvUnbounded   the affine convex envelope on one such piece
%
% WHY SEPARATELY. All four are reached by the slow bucket, so they are not untested in the sense
% of "never run" -- but they are only ever run in composition, where a defect in one shows up as a
% wrong number many stages later, if at all. Measured 2026-08-31 on fast+normal:
% convEnvUnbounded 0%, fanUnboundedFace 0%, xyFrame 0%, transformDomain 0%. Each has a contract
% stated in its own header in one or two sentences; this asserts those sentences directly.
%
% BUCKET: fast. `convEnvUnbounded` and `fanUnboundedFace` build symbolic regions, so this is the
% slowest member of the fast bucket -- everything is kept to small fixtures for that reason.

    properties (Constant)
        X = sym('x')
        Y = sym('y')
    end

    methods (Static)
        function v = vars()
            v = [frameAndFanTest.X, frameAndFanTest.Y];
        end

        function q = quadSym(Q, L, c, v)
            z = [v(1); v(2)];
            q = expand(0.5*(z.'*Q*z) + L(:).'*z + c);
        end

        function val = quadNum(Q, L, c, p)
            z = p(:);
            val = 0.5*(z.'*Q*z) + L(:).'*z + c;
        end
    end

    methods (Test)

        % =========================================================================================
        % xyFrame -- the congruence is the whole contract, and it is checkable in one line
        % =========================================================================================
        function xyFrameProducesTheXYCongruence(testCase)
        % The header promises M'QM = [0 1; 1 0] and q(Mz) = z1*z2 + a'z + c0. Both are asserted,
        % on a family of indefinite Q rather than one -- including a nearly-degenerate one, since
        % the construction divides by sqrt of each eigenvalue and that is where it would fail.
            v = frameAndFanTest.vars();
            Qs = { [0 1; 1 0], [1 0; 0 -1], [2 3; 3 -1], [1 4; 4 1], [1e-3 1; 1 -1e-3], [5 0; 0 -0.01] };
            rng(20260831);
            for k = 1:numel(Qs)
                Q = Qs{k}; L = round((rand(2,1)-0.5)*4, 3); c = round((rand-0.5)*4, 3);
                testCase.assumeLessThan(det(Q), 0, 'the fixture must be indefinite');
                [M, a, c0] = xyFrame(Q, L, c);

                % xyFrame answers in SYM, deliberately -- its output feeds expressions, not
                % arithmetic (see conjConvexOverPiece's exact-values/numeric-decisions rule). The
                % assertions are about the numbers, so they are made on double(...).
                Md = double(M);
                testCase.verifyEqual(Md.'*Q*Md, [0 1; 1 0], 'AbsTol', 1e-9, sprintf( ...
                    'xyFrame case %d: M''QM is not the x*y form', k));
                testCase.verifyEqual(double(a(:)), Md.'*L(:), 'AbsTol', 1e-9, sprintf( ...
                    'xyFrame case %d: a must be M''L', k));
                testCase.verifyEqual(double(c0), c, 'AbsTol', 1e-12, sprintf( ...
                    'xyFrame case %d: the constant is unchanged', k));

                % The point of the frame: q(Mz) really is z1*z2 + a'z + c0, at arbitrary z.
                Z = (rand(6,2)-0.5)*6;
                for i = 1:size(Z,1)
                    z = Z(i,:).';
                    lhs = frameAndFanTest.quadNum(Q, L, c, Md*z);
                    rhs = z(1)*z(2) + double(a(:)).'*z + double(c0);
                    testCase.verifyEqual(lhs, rhs, 'AbsTol', 1e-8*max(1,abs(rhs)), sprintf( ...
                        'xyFrame case %d: q(Mz) at z = (%g,%g)', k, z(1), z(2)));
                end
                testCase.verifyGreaterThan(abs(det(Md)), 1e-12, sprintf( ...
                    'xyFrame case %d: M must be invertible or the frame cannot be undone', k));
                % `quadSym` exists so the symbolic and numeric readings of the same q cannot drift
                % apart unnoticed in the tests below.
                qs = frameAndFanTest.quadSym(Q, L, c, v);
                testCase.verifyEqual(double(subs(qs, v, [1 -2])), ...
                    frameAndFanTest.quadNum(Q, L, c, [1 -2]), 'AbsTol', 1e-12, ...
                    'the symbolic and numeric fixtures disagree');
            end
        end

        % =========================================================================================
        % transformDomain -- membership must be carried across exactly
        % =========================================================================================
        function transformDomainPreservesMembership(testCase)
        % z = Minv*x, so a point x lies in the original domain exactly when Minv*x lies in the
        % transformed one. Checked in BOTH directions on a sampled grid, for a bounded triangle
        % and for an unbounded wedge -- the header's claim that a ray needs no special handling is
        % the half worth testing.
            v = frameAndFanTest.vars();
            Minv = [2 1; -1 1];                        % invertible, not orthogonal, not diagonal
            M = inv(Minv);

            % `domain` is constructed from a VERTEX LIST, not from a region -- domain(V,x,y).
            % An UNBOUNDED domain cannot be written that way at all; fanUnboundedFace builds those
            % field by field, so the unbounded fixture here is taken from its output, which is
            % also the only shape transformDomain is ever handed in production.
            wedge = plqCheck.wedgeRegion([1 1], [1 0], [0 1], v);
            fanned = fanUnboundedFace(wedge, v(1), v(2));
            cases = { domain([0 0; 2 0; 0 3], v(1), v(2)), plqCheck.triRegion([0 0; 2 0; 0 3], v), 'triangle'
                      fanned{1},                           fanned{1}.polygon,                     'fanned wedge piece' };
            for k = 1:size(cases,1)
                d0 = cases{k,1}; r0 = cases{k,2}; nm = cases{k,3};
                d  = transformDomain(d0, Minv, v);
                box = plqCheck.regionBox(r0);
                P = plqCheck.regionSample(r0, box, 41, 1e-7);
                testCase.verifyGreaterThan(size(P,1), 0, sprintf('%s: nothing sampled', nm));
                Z = (Minv * P.').';
                inNew = plqCheck.inRegion(d.polygon, Z, 1e-7);
                testCase.verifyTrue(all(inNew), sprintf( ...
                    '%s: %d of %d points of the original domain left the transformed one', ...
                    nm, sum(~inNew), numel(inNew)));

                % ...and nothing was gained: points of the new domain map back inside the old.
                boxZ = plqCheck.regionBox(d.polygon);
                Q = plqCheck.regionSample(d.polygon, boxZ, 41, 1e-7);
                if ~isempty(Q)
                    Xback = (M * Q.').';
                    inOld = plqCheck.inRegion(r0, Xback, 1e-7);
                    testCase.verifyTrue(all(inOld), sprintf( ...
                        '%s: the transformed domain gained %d points the original does not have', ...
                        nm, sum(~inOld)));
                end
            end
        end

        % =========================================================================================
        % fanUnboundedFace -- a COVER by SUBSETS, which is exactly what its header insists on
        % =========================================================================================
        function fanCoversTheFaceWithoutLeavingIt(testCase)
        % Two assertions, and the header says why the second is the dangerous one: the convex
        % envelope does not decompose over a cover, only the conjugate does, and that argument
        % needs every piece to be a SUBSET. A superset inflates its own sup and hence the max.
            v = frameAndFanTest.vars();
        % The third fixture is the one with more than one finite vertex, which is the only case
        % where the fan has anything to do: conv{(1,0),(0,1)} + the first quadrant, i.e. the
        % quadrant with a corner cut off. (A quadrant cut by x + y <= 4 is BOUNDED and
        % fanUnboundedFace refuses it by design -- see the guard test below.)
            cases = { plqCheck.wedgeRegion([0 0], [1 0], [0 1], v),            'quadrant wedge'
                      plqCheck.halfStripRegion([0 0], [2 0], [0 1], v),        'half-strip'
                      plqCheck.halfPlanes([-1 0; 0 -1; -1 -1], [0; 0; -1], v), 'quadrant less a corner' };
            for k = 1:size(cases,1)
                r = cases{k,1}; nm = cases{k,2};
                ds = fanUnboundedFace(r, v(1), v(2));
                testCase.verifyGreaterThanOrEqual(numel(ds), 1, sprintf('%s: no pieces', nm));

                box = plqCheck.regionBox(r);
                % SUBSET: nothing any piece contains may lie outside r.
                for i = 1:numel(ds)
                    Pi = plqCheck.regionSample(ds{i}.polygon, box, 31, 1e-7);
                    if isempty(Pi), continue, end
                    inR = plqCheck.inRegion(r, Pi, 1e-6);
                    bad = find(~inR, 1);
                    testCase.verifyTrue(all(inR), sprintf( ...
                        ['%s: piece %d contains %d points OUTSIDE the face (first (%g,%g)) -- a ' ...
                         'superset inflates its own sup and corrupts the max across pieces'], ...
                        nm, i, sum(~inR), Pi(max(bad,1),1), Pi(max(bad,1),2)));
                end
                % COVER: every point of r must land in some piece.
                P = plqCheck.regionSample(r, box, 31, 1e-6);
                testCase.verifyGreaterThan(size(P,1), 0, sprintf('%s: nothing sampled', nm));
                covered = false(size(P,1),1);
                for i = 1:numel(ds)
                    covered = covered | plqCheck.inRegion(ds{i}.polygon, P, 1e-6);
                end
                miss = find(~covered, 1);
                testCase.verifyTrue(all(covered), sprintf( ...
                    '%s: %d of %d points of the face are in no piece (first (%g,%g))', ...
                    nm, sum(~covered), numel(covered), P(max(miss,1),1), P(max(miss,1),2)));
            end
        end

        function fanRefusesAFaceItIsNotFor(testCase)
        % The guard is part of the contract: the routine covers an UNBOUNDED, POINTED face, and a
        % bounded one belongs on plq_1p.triangulate. Asserted by identifier rather than by message
        % so a reworded error does not make this red -- and asserted at all because a guard that
        % silently passed a bounded face would hand Step 1 a cover it cannot use.
            v = frameAndFanTest.vars();
            bounded = plqCheck.triRegion([0 0; 4 0; 0 4], v);
            testCase.verifyError(@() fanUnboundedFace(bounded, v(1), v(2)), ...
                'fanUnboundedFace:notPointed', ...
                'a bounded face must be refused, not fanned');
        end

        % =========================================================================================
        % convEnvUnbounded -- an affine MINORANT that TOUCHES, which is what "envelope" means here
        % =========================================================================================
        function convEnvUnboundedReturnsATouchingAffineMinorant(testCase)
        % The routine returns the envelope only in the case where the envelope is affine, so the
        % two halves of "co q" are both checkable directly:
        %   MINORANT   expr <= q at every sampled point of the piece. A violation is a definite
        %              defect -- an envelope above the function is not an underestimator.
        %   TIGHT      the gap q - expr must reach ~0 somewhere on the piece, or `expr` is merely
        %              some minorant and not the greatest one.
        % Plus: the result must be AFFINE, which is the routine's own precondition and the thing a
        % future change would break silently.
            v = frameAndFanTest.vars();
            % Indefinite q on a wedge whose two rays are both q-null directions -- the header's
            % first shape. q = x*y on the first quadrant: both axes are null and co q = 0 there.
            cases = { plqCheck.wedgeRegion([0 0], [1 0], [0 1], v), [0 1; 1 0], [0;0], 0, 'x*y on the quadrant'
                      plqCheck.wedgeRegion([1 2], [1 0], [0 1], v), [0 1; 1 0], [1;-1], 3, 'x*y shifted' };
            for k = 1:size(cases,1)
                r = cases{k,1}; Q = cases{k,2}; L = cases{k,3}; c = cases{k,4}; nm = cases{k,5};
                q = frameAndFanTest.quadSym(Q, L, c, v);
                [expr, why] = convEnvUnbounded(r, q, v);
                testCase.verifyNotEmpty(why, sprintf('%s: the shape verdict must be reported', nm));

                testCase.verifyEqual(double(diff(expr, v(1), 2)), 0, 'AbsTol', 1e-12, sprintf( ...
                    '%s: the envelope must be AFFINE in x', nm));
                testCase.verifyEqual(double(diff(expr, v(2), 2)), 0, 'AbsTol', 1e-12, sprintf( ...
                    '%s: the envelope must be AFFINE in y', nm));
                testCase.verifyEqual(double(diff(diff(expr, v(1)), v(2))), 0, 'AbsTol', 1e-12, ...
                    sprintf('%s: the envelope must have no cross term', nm));

                box = plqCheck.regionBox(r);
                P = plqCheck.regionSample(r, box, 35, 1e-7);
                testCase.verifyGreaterThan(size(P,1), 0, sprintf('%s: nothing sampled', nm));
                h = matlabFunction(expr, 'Vars', {v(1), v(2)});
                gap = zeros(size(P,1),1);
                for i = 1:size(P,1)
                    gap(i) = frameAndFanTest.quadNum(Q, L, c, P(i,:)) - double(h(P(i,1), P(i,2)));
                end
                sc = max(1, max(abs(gap)));
                [worst, iw] = min(gap);
                testCase.verifyGreaterThanOrEqual(worst, -1e-7*sc, sprintf( ...
                    '%s: the "envelope" exceeds q by %.3g at (%g,%g) -- it is not a minorant', ...
                    nm, -worst, P(iw,1), P(iw,2)));

                % TIGHTNESS is checked at the piece's own VERTICES, not on the interior grid. For
                % q = x*y on a quadrant the envelope is 0 and the gap is x*y, which vanishes only
                % ON the boundary -- the grid samples strict interior points, so its smallest gap
                % is the grid step squared and says nothing. The vertices are where an affine
                % envelope of a concave-along-every-line q has to touch.
                Vf = plqCheck.finiteVertexSet(r);
                testCase.verifyGreaterThan(size(Vf,1), 0, sprintf('%s: no finite vertex', nm));
                for i = 1:size(Vf,1)
                    gv = frameAndFanTest.quadNum(Q, L, c, Vf(i,:)) - double(h(Vf(i,1), Vf(i,2)));
                    testCase.verifyEqual(gv, 0, 'AbsTol', 1e-7*sc, sprintf( ...
                        ['%s: the minorant does not touch q at the vertex (%g,%g) -- gap %.3g, ' ...
                         'so it is not the greatest affine minorant'], nm, Vf(i,1), Vf(i,2), gv));
                end
            end
        end

        function convEnvUnboundedOnAHalfStripIsAlsoATouchingAffineMinorant(testCase)
        % The 'ray' branch -- the largest of the three and the one the two tests above do not
        % reach, since a wedge dispatches on 'wedge' and a triangle on 'bounded'.
        %
        % ITS HYPOTHESIS, from the header: conv{v1,v2} + cone(d) with d'Qd = 0 and w'Qw <= 0 for
        % w = v2 - v1. For q = x*y the null directions of Q are the two axes, so d = (0,1); and
        % w'Qw = 2*w1*w2 <= 0 needs w with components of opposite sign, so w = (1,-1). The
        % fixture is built to satisfy exactly that, because a fixture that misses the hypothesis
        % tests the guard rather than the branch.
            v = frameAndFanTest.vars();
            r = plqCheck.halfStripRegion([0 0], [1 -1], [0 1], v);
            Q = [0 1; 1 0]; L = [0;0]; c = 0;
            q = frameAndFanTest.quadSym(Q, L, c, v);
            [expr, why] = convEnvUnbounded(r, q, v);
            testCase.verifyNotEmpty(why, 'the shape verdict must be reported');

            testCase.verifyEqual(double(diff(expr, v(1), 2)), 0, 'AbsTol', 1e-12, ...
                'the envelope must be affine in x');
            testCase.verifyEqual(double(diff(expr, v(2), 2)), 0, 'AbsTol', 1e-12, ...
                'the envelope must be affine in y');
            testCase.verifyEqual(double(diff(diff(expr, v(1)), v(2))), 0, 'AbsTol', 1e-12, ...
                'the envelope must have no cross term');

            box = plqCheck.regionBox(r);
            P = plqCheck.regionSample(r, box, 35, 1e-7);
            testCase.verifyGreaterThan(size(P,1), 0, 'nothing sampled on the half-strip');
            h = matlabFunction(expr, 'Vars', {v(1), v(2)});
            gap = zeros(size(P,1),1);
            for i = 1:size(P,1)
                gap(i) = frameAndFanTest.quadNum(Q, L, c, P(i,:)) - double(h(P(i,1), P(i,2)));
            end
            sc = max(1, max(abs(gap)));
            [worst, iw] = min(gap);
            testCase.verifyGreaterThanOrEqual(worst, -1e-7*sc, sprintf( ...
                'the half-strip envelope exceeds q by %.3g at (%g,%g)', -worst, P(iw,1), P(iw,2)));

            Vf = plqCheck.finiteVertexSet(r);
            testCase.verifyGreaterThan(size(Vf,1), 0, 'the half-strip must report finite vertices');
            for i = 1:size(Vf,1)
                gv = frameAndFanTest.quadNum(Q, L, c, Vf(i,:)) - double(h(Vf(i,1), Vf(i,2)));
                testCase.verifyEqual(gv, 0, 'AbsTol', 1e-7*sc, sprintf( ...
                    'the envelope must touch q at the finite vertex (%g,%g) -- gap %.3g', ...
                    Vf(i,1), Vf(i,2), gv));
            end
        end

        function theThreeShapesDispatchToThreeDifferentBranches(testCase)
        % `why` is the routine's own report of which of the three hypotheses decided the answer.
        % If two shapes came back with the same verdict, one of them would be going down a branch
        % written for the other -- which the value checks above could still pass by accident on a
        % symmetric fixture. Cheap, and it is the only thing that pins the dispatch itself.
            v = frameAndFanTest.vars();
            Q = [0 1; 1 0]; L = [0;0]; c = 0;
            q = frameAndFanTest.quadSym(Q, L, c, v);
            shapes = { plqCheck.wedgeRegion([0 0], [1 0], [0 1], v),      'wedge'
                       plqCheck.halfStripRegion([0 0], [1 -1], [0 1], v), 'half-strip'
                       plqCheck.triRegion([0 0; 3 0; 0 2], v),            'bounded triangle' };
            whys = cell(size(shapes,1),1);
            for k = 1:size(shapes,1)
                [~, whys{k}] = convEnvUnbounded(shapes{k,1}, q, v);
                testCase.verifyNotEmpty(whys{k}, sprintf('%s: no verdict', shapes{k,2}));
            end
            testCase.verifyEqual(numel(unique(whys)), numel(whys), sprintf( ...
                ['the three shapes must take three different branches, but the verdicts were ' ...
                 '%s -- a repeat means one shape is being handled by another''s formula'], ...
                strjoin(whys.', ', ')));
        end

        function convEnvUnboundedOnABoundedTriangleAgreesWithTheVertexHull(testCase)
        % The header says a bounded triangle is also accepted. There the affine envelope of a
        % quadratic that is concave along every edge is the affine interpolant of the three VERTEX
        % values -- an independent formula, so it is worth asserting against rather than resampling.
            v = frameAndFanTest.vars();
            V = [0 0; 3 0; 0 2];
            r = plqCheck.triRegion(V, v);
            Q = [0 1; 1 0]; L = [0;0]; c = 0;                 % x*y: concave along every straight line
            q = frameAndFanTest.quadSym(Q, L, c, v);
            [expr, ~] = convEnvUnbounded(r, q, v);
            h = matlabFunction(expr, 'Vars', {v(1), v(2)});
            for i = 1:3
                testCase.verifyEqual(double(h(V(i,1), V(i,2))), ...
                    frameAndFanTest.quadNum(Q, L, c, V(i,:)), 'AbsTol', 1e-8, sprintf( ...
                    'the affine envelope must equal q at vertex %d of the triangle', i));
            end
        end
    end
end
