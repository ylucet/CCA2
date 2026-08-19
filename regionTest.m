classdef regionTest < matlab.unittest.TestCase
    % Tests for region.m's LP certificates -- the real redundancy test that replaced
    % simplifyUnboundedRegion's "does this constraint pass through a finite vertex" proxy, and
    % the real union-exactness test that region.merge now requires before deleting a shared
    % facet. Both proxies used to delete constraints that carry weight, leaving a region that
    % claims territory belonging to neither operand and carries the wrong function value on it
    % (SUPPORT_MATRIX.md section 1.2).
    %
    % These are unit tests on region alone: the end-to-end consequence for Step 3's assembly is
    % pinned by conjCPLQTest.indefiniteTriangleThreeConvexEdgesUsesStep3.

    methods (Static)
        function [i1, i2] = sharedFacet(r1, r2)
        % The (i,j) with r1.ineqs(i) == -r2.ineqs(j), i.e. the facet the two regions meet on.
            i1 = 0; i2 = 0;
            for i = 1:size(r1.ineqs,2)
                for j = 1:size(r2.ineqs,2)
                    if r1.ineqs(i) == -r2.ineqs(j)
                        i1 = i; i2 = j; return
                    end
                end
            end
        end

        function [bad, per] = coneVsDefinition(r, NC, s1, s2)
        % THE SPECIFICATION OF A VERTEX CONE, as an oracle built from the definition rather than
        % from any implementation: `u` is in the normal cone at `v` iff `v` maximises `<u,.>`
        % over the region NEAR `v` (the local form -- these regions are not all convex).
        %
        % What the rows of NC mean is fixed by the CONSUMER, functionNDomain.getSubdiffVertexT1:
        % it throws the constant term away and re-anchors each line at `grad f(v_j)`, and
        % region's convention is that `expr <= 0` is feasible. So the object being specified is
        % the LINEAR PART under `<= 0`, and nothing else about the row is observable.
        %
        % Returns the number of directions (of 72) where the two disagree CLEARLY -- directions
        % whose supporting value is within sampling noise of zero are the cone's own boundary and
        % are not counted.
            x = r.vars(1); y = r.vars(2);
            gf = matlabFunction([r.ineqs.f], 'Vars', {[x y]});
            nd = 72; ang = (0:nd-1)*(2*pi/nd); U = [cos(ang)', sin(ang)']; epsr = 0.05;
            per = zeros(1, r.nv);
            for j = 1:r.nv
                v = [double(r.vx(j)), double(r.vy(j))];
                W = [];
                for rr = epsr*(0.2:0.2:1.0)
                    th = (0:359)*pi/180;
                    W = [W; v + rr*[cos(th)', sin(th)']];   %#ok<AGROW>
                end
                keep = false(size(W,1),1);
                for i = 1:size(W,1), keep(i) = all(gf(W(i,:)) <= 1e-12); end
                D = (W(keep,:) - v)/epsr;

                inCode = true(nd,1);
                for k = 1:size(NC,2)
                    e = NC(j,k);
                    if isAlways(e == 0, 'Unknown', 'false'), continue, end
                    c0 = double(subs(e, [s1 s2], [0 0]));
                    L = zeros(nd,1);
                    for i = 1:nd, L(i) = double(subs(e, [s1 s2], U(i,:))) - c0; end
                    inCode = inCode & (L <= 1e-9);
                end
                m = zeros(nd,1);
                for i = 1:nd, m(i) = max(D*U(i,:)'); end
                inTrue = m <= 1e-9; amb = (m > 1e-9) & (m < 1e-3);
                per(j) = sum((inCode ~= inTrue) & ~amb);
            end
            bad = sum(per);
        end
    end

    methods (Test)
        function normalConesOnCurvedEdgesAreUnchanged(testCase)
        % A CHARACTERIZATION test, and it exists to make one refactor safe.
        %
        % getNormalConeVertexQ finds its probe points with solve(), takes the FIRST root, tests
        % it, and flips the probe to the other side if it fails -- it never tries the second
        % root. So anything that changes root ORDER can silently change the cone, and the suites
        % that exercise curved cones (conjCPLQTest, testMaxMultiRegion) are in the SLOW bucket,
        % ~100 minutes. This pins the same behaviour in about a second, so the ordering can be
        % fixed and checked immediately.
        %
        % The values below are the output of the suite-green implementation on 2026-08-18, not a
        % hand derivation: the point is that a refactor must not MOVE them. They are NOT normal cones,
        % and are not claimed to be: called without eIdx this routine falls back to the slot
        % convention of the UNBOUNDED layout, and all three fixtures are BOUNDED, so it reads the
        % wrong pair of constraints at some vertices (measured: 32, 43 and 5 wrong directions of
        % 72 per region). The cones themselves are specified and checked in
        % vertexConesMatchTheDefinition, which supplies the edge list. Three regions, each
        % with a genuine conic facet -- piece 9 of the parallelogram's f*, a half-lens, and a
        % parabola capped by two lines.
            x = sym('x'); y = sym('y'); s1 = sym('s_1'); s2 = sym('s_2');

            cases = { ...
                region([16*x - 4*x*y - x^2 - 4*y^2, -x-2*y, (2*y)/3 - x - sym(1)/3, ...
                        x + 2*y - 2], [x y]), ...
                    {'s_2 - 1/8', '2*s_1 - s_2 + 5/8'; ...
                     's_2 - 2*s_1', '(2*s_1)/3 + s_2'; ...
                     '25/24 - s_2 - (2*s_1)/3', '2*s_1 - s_2 + 3/8'}, 'piece 9'; ...
                region([(x+y)^2 - 4*x, -x, -y], [x y]), ...
                    {'-s_2', 's_2'; '-s_2', '4 - s_1'}, 'half lens'; ...
                region([x^2 - y, y - 1, -x - 1], [x y]), ...
                    {'s_1/2 - s_2 + 3/2', 's_1 + 1'; '1 - s_1', '1 - s_2'}, 'parabola' };

            for k = 1:size(cases,1)
                r = cases{k,1}; want = cases{k,2}; nm = cases{k,3};
                NC = r.getNormalConeVertexQ(s1, s2);
                testCase.verifyEqual(size(NC), size(want), ...
                    sprintf('%s: cone matrix changed shape', nm));
                for a = 1:size(want,1)
                    for b = 1:size(want,2)
                        testCase.verifyTrue( ...
                            isAlways(simplify(NC(a,b) - str2sym(want{a,b})) == 0, ...
                                     'Unknown', 'false'), ...
                            sprintf('%s: NC(%d,%d) moved -- got %s, pinned %s', ...
                                    nm, a, b, char(NC(a,b)), want{a,b}));
                    end
                end
            end
        end

        function vertexConesMatchTheDefinition(testCase)
        % THE SPECIFICATION, and it is executable. `getNormalConeVertexQ` must return, for each
        % vertex j, two linear forms in the dual variables whose LINEAR PARTS, read as `<= 0`,
        % cut out the normal cone of the region at that vertex. Two things fix that reading, and
        % both come from the consumer `functionNDomain.getSubdiffVertexT1`:
        %
        %   * it DISCARDS the constant term and re-anchors each line at `grad f(v_j)` -- so the
        %     conjugate cell it builds is `grad f(v_j) + N_D(v_j) = subdiff f(v_j)`, which is
        %     the identity the whole vertex branch of the conjugate rests on;
        %   * the region built from those rows is `region(subdV(j,:), ...)`, and region's
        %     convention is `expr <= 0` feasible (`ptFeasible`).
        %
        % MEASURED 2026-08-18: given the correct edge list, the committed routine satisfies this
        % on every vertex of all three curved cases -- including the CUSP (half-lens vertex 1)
        % and the vertex where three constraints are active (piece 9 vertex 3). The earlier
        % finding that it "does not compute a normal cone" was measured on the eIdx-less SLOT
        % fallback, which picks the wrong pair of constraints; see DECISIONS.md.
            x = sym('x'); y = sym('y'); s1 = sym('s_1'); s2 = sym('s_2');

            % {region, edge list, name}. eIdx(j) is the constraint bounding the edge from vertex
            % j to vertex j+1; for the lens it must be the TRAVERSAL order, not just the right
            % pair -- reversing it moves the cone at vertex 1 (the routine deduces orientation
            % from a probe on the ARRIVING edge).
            cases = { ...
                region([16*x - 4*x*y - x^2 - 4*y^2, -x-2*y, (2*y)/3 - x - sym(1)/3, ...
                        x + 2*y - 2], [x y]), [2 1 3], 'piece 9'; ...
                region([(x+y)^2 - 4*x, -x, -y], [x y]), [1 3], 'half lens'; ...
                region([x^2 - y, y - 1, -x - 1], [x y]), [2 1], 'parabola' };

            for k = 1:size(cases,1)
                r = cases{k,1}; eIdx = cases{k,2}; nm = cases{k,3};
                NC = r.getNormalConeVertexQ(s1, s2, eIdx);

                % NORMALISATION is part of the contract, not cosmetics: getSubdiffVertexT1 uses
                % `m = d(row)/d(s1)` as the slope of the re-anchored line, which is the row's own
                % slope only when the coefficient of s2 is +-1.
                for a = 1:size(NC,1)
                    for b = 1:size(NC,2)
                        if isAlways(NC(a,b) == 0, 'Unknown', 'false'), continue, end
                        c = symbolicFunction(NC(a,b)).getLinearCoeffs([s1 s2]);
                        if c(2) == 0
                            ok = abs(double(c(1))) == 1;
                        else
                            ok = abs(double(c(2))) == 1;
                        end
                        testCase.verifyTrue(logical(ok), sprintf( ...
                            '%s: NC(%d,%d) = %s is not normalised for getSubdiffVertexT1', ...
                            nm, a, b, char(NC(a,b))));
                    end
                end

                [bad, per] = regionTest.coneVsDefinition(r, NC, s1, s2);
                testCase.verifyEqual(bad, 0, sprintf( ...
                    '%s: cone disagrees with the definition on %d of 72 directions [%s]', ...
                    nm, bad, num2str(per)));
            end
        end

        % ---- simplifyUnboundedRegion ----------------------------------------------------
        function aHalfPlaneIsNotEmpty(testCase)
        % simplifyUnboundedRegion decided a region had no interior from the count of its FINITE
        % vertices, and a half-plane has none -- nor does a slab, nor the whole plane. Every one
        % of them came back region.empty(), which deletes an answer silently.
        %
        % A HALF-PLANE IS EXACTLY WHAT A TANGENT VERTEX PRODUCES, so this is not a corner case.
        % getNormalConeVertexQ builds the cone at a vertex from the two edges meeting there;
        % when those edges are TANGENT -- an arc and its chord touching, the way a curvilinear
        % piece ends -- both half-planes are the SAME one and the cone is a half-plane. On piece
        % 9 of f* for x*y over the parallelogram conv{(0,0),(2,0),(2.5,1),(0.5,1)} that cone is
        % {2x/3 + y >= 4/3}; dropping it left the conjugate uncovered on exactly that half-plane
        % and the biconjugate then collapsed to nothing.
        %
        % All four spellings are checked because the caller produces the duplicated one: two
        % identical half-planes plus a trivial 0 is what a tangent vertex's two constraints and
        % the unused third slot come to.
            x = sym('x'); y = sym('y');
            a = sym(4)/3 - y - (2*x)/3;              % the cone at that vertex, {2x/3 + y >= 4/3}
            forms = { a, [a, a], [a, sym(0)], [a, a, sym(0)] };
            for k = 1:numel(forms)
                r = region(forms{k}, [x,y]);
                testCase.verifyTrue(r.ptFeasible([x,y], [3,3]), ...
                    'precondition: (3,3) is in this half-plane');
                r2 = r.simplifyUnboundedRegion;
                testCase.verifyFalse(isempty(r2), sprintf( ...
                    'a half-plane is not empty (form %d)', k));
                testCase.verifyTrue(r2.ptFeasible([x,y], [3,3]), sprintf( ...
                    'and it still contains its own points (form %d)', k));
            end
        end

        function aGenuinelyEmptyRegionIsStillReportedEmpty(testCase)
        % The other side of the test above: refuting an emptiness verdict must need a WITNESS, so
        % a region that really is empty has to survive the change. Two opposite half-planes with
        % no common point have none.
            x = sym('x'); y = sym('y');
            r = region([x - 1, 2 - x], [x,y]);       % x >= 2 and x <= 1
            r2 = r.simplifyUnboundedRegion;
            testCase.verifyTrue(isempty(r2), 'x >= 2 and x <= 1 is empty');
        end

        % ---- maxLinear ------------------------------------------------------------------
        function maxLinearReportsOptimalUnboundedAndInfeasible(testCase)
            % The three answers callers branch on. Unboundedness and infeasibility are
            % first-class results here, not failures -- these regions are routinely unbounded,
            % so a primitive that could only answer "optimal" would be useless.
            [v, st] = region.maxLinear([1 0; 0 1; -1 -1], [1;1;0], [1 1]);
            testCase.verifyEqual(st, 0);
            testCase.verifyEqual(v, 2, 'AbsTol', 1e-9);

            [v, st] = region.maxLinear([-1 0], 0, [1 0]);
            testCase.verifyEqual(st, 1);
            testCase.verifyEqual(v, inf);

            [~, st] = region.maxLinear([1 0; -1 0], [-1; -1], [1 0]);
            testCase.verifyEqual(st, -1);
        end

        % ---- linearForm -----------------------------------------------------------------
        function linearFormFlagsOnlyTheAffineConstraints(testCase)
            % A quadratic facet gets no row, and lin says so -- every caller must consult it
            % rather than reading a NaN row as a constraint.
            syms x y
            r = region([-x, -y, x+y-2, x^2+y^2-9], [x y]);
            [A, b, lin] = r.linearForm;
            testCase.verifyEqual(nnz(lin), 3);
            testCase.verifyEqual(sort(b(lin))', [0 0 2], 'AbsTol', 1e-12);
            testCase.verifyTrue(all(isfinite(A(lin,:)), 'all'));
        end

        % ---- redundancy -----------------------------------------------------------------
        function redundantSubsetFindsTheImpliedConstraint(testCase)
            % x<=3 is implied by x<=1, so it -- and only it -- may be deleted.
            syms x y
            r = region([-x, -y, x-1, x-3, y-1], [x y]);
            testCase.verifyEqual(r.redundantSubset(1:size(r.ineqs,2)), 4);
        end

        function redundantSubsetKeepsEveryFacetOfASimplex(testCase)
            % No constraint of a triangle is implied by the other two. This is the case the old
            % proxy got wrong whenever the region was unbounded: it asked whether a constraint
            % passes through a finite vertex, which is not a redundancy test at all.
            syms x y
            r = region([-x, -y, x+y-2], [x y]);
            testCase.verifyEmpty(r.redundantSubset(1:size(r.ineqs,2)));
        end

        function redundantSubsetKeepsAnIrredundantConstraintOfAnUnboundedRegion(testCase)
            % The measured instance, reduced: the halfplane -x-y<=0 bounds this unbounded region
            % but touches no other constraint at a finite vertex, so the old rule deleted it.
            % Dropping it lets the region reach (-3,-3), which is exactly how a Step 3 region
            % carrying (s1+s2)^2/4 came to report f*(-3,-3)=9 instead of 0.
            syms x y
            r = region([-x-y, x-y-1, y-x-1], [x y]);
            testCase.verifyEmpty(r.redundantSubset(1:size(r.ineqs,2)));
            testCase.verifyFalse(r.ptFeasible([x y], [-3 -3]));
        end

        function redundantSubsetNeverDeletesBothCopiesOfAConstraint(testCase)
            % Each of two identical constraints is redundant GIVEN the other, so a test that
            % judged all candidates against the original set would delete both and open the
            % region up. Candidates are judged against the constraints still standing.
            syms x y
            r = region([-x, -y, x+y-2], [x y]);
            r.ineqs(end+1) = r.ineqs(3);                  % a literal duplicate of x+y-2<=0
            del = r.redundantSubset(1:size(r.ineqs,2));
            testCase.verifyEqual(numel(del), 1);
            r2 = r.deleteIfRedundant(1:size(r.ineqs,2));
            testCase.verifyTrue(r2.ptFeasible([x y], [0.5 0.5]));
            testCase.verifyFalse(r2.ptFeasible([x y], [2 2]));
        end

        % ---- union exactness / merge ----------------------------------------------------
        function mergeUnionsTwoRegionsWhoseUnionIsConvex(testCase)
            % [0,1]x[0,1] and [1,2]x[0,1] meet on x=1 and their union IS [0,2]x[0,1], so the
            % delete-the-facet-and-intersect recipe is exact here and the merge must happen.
            syms x y
            r1 = region([-x, -y, x-1, y-1], [x y]);
            r2 = region([1-x, -y, x-2, y-1], [x y]);
            [i1, i2] = regionTest.sharedFacet(r1, r2);
            testCase.verifyGreaterThan(i1, 0);
            testCase.verifyTrue(r1.unionIsExact(r2, i1, i2));

            [l, m] = r1.merge(r2);
            testCase.verifyTrue(l);
            testCase.verifyTrue(m.ptFeasible([x y], [1.5 0.5]));   % r2's half
            testCase.verifyTrue(m.ptFeasible([x y], [0.5 0.5]));   % r1's half
            testCase.verifyFalse(m.ptFeasible([x y], [1.5 1.5]));  % in neither
        end

        function mergeRefusesWhenTheUnionWouldOverClaim(testCase)
            % Same shared facet x=1, but the second box is taller: the union is an L, not
            % convex, and A' n B' is the bounding box [0,2]x[0,2]. The old recipe returned that
            % box and carried r1's function value onto (0.5,1.5), a point in NEITHER operand --
            % the defect that made three Step 3 regions merge into one covering s=(1,1).
            syms x y
            r1 = region([-x, -y, x-1, y-1], [x y]);
            r2 = region([1-x, -y, x-2, y-2], [x y]);
            [i1, i2] = regionTest.sharedFacet(r1, r2);
            testCase.verifyGreaterThan(i1, 0);
            testCase.verifyFalse(r1.unionIsExact(r2, i1, i2));

            [l, m] = r1.merge(r2);
            testCase.verifyFalse(l);
            testCase.verifyFalse(m.ptFeasible([x y], [0.5 1.5]));  % r1 returned untouched
        end

        function mergeRefusesWhenTwoFacetsAreShared(testCase)
            % With two shared facets the "a merged point is in one operand or the other"
            % argument fails outright -- a point on the <=0 side of one and the >=0 side of the
            % other lies in neither yet survives the intersection. A vertical strip and its
            % complement-side strip share BOTH walls.
            syms x y
            r1 = region([x-1, -x-1, y-1, -y-1], [x y]);        % [-1,1]x[-1,1]
            r2 = region([1-x, x+1, y-1, -y-1], [x y]);         % 1<=x<=-1: infeasible/degenerate
            [l, ~] = r1.merge(r2);
            testCase.verifyFalse(l);
        end

        function quadUnboundedBelowDecidesTheEnvelopeGate(testCase)
        % conv q over an UNBOUNDED face is -inf exactly when q is unbounded below on it, so this
        % is the gate Step 1 has to pass before there is any finite envelope to compute. The test
        % is "does the recession cone meet {d : d'Qd < 0}", with the d'Qd == 0 directions decided
        % by the linear slope instead.
            x = sym('x'); y = sym('y');
            bounded = region([-x, -y, x+y-2], [x y]);      % triangle (0,0),(2,0),(0,2)
            wedge   = region([x, -y], [x y]);              % {x<=0, y>=0}: recession cone itself
            halfpl  = region([-y], [x y]);                 % {y>=0}

            % A bounded region recedes along nothing, so no quadratic is unbounded below on it.
            testCase.verifyFalse(bounded.quadUnboundedBelow([1 0; 0 -1], [0; 0]));
            testCase.verifyFalse(bounded.quadUnboundedBelow(-eye(2), [3; -7]));

            % On the wedge, x^2+y^2 grows along every recession direction ...
            testCase.verifyFalse(wedge.quadUnboundedBelow(eye(2), [0; 0]));
            % ... while x^2-y^2 decays along (0,1), which the wedge never leaves.
            testCase.verifyTrue(wedge.quadUnboundedBelow(diag([1 -1]), [0; 0]));
            testCase.verifyTrue(wedge.quadUnboundedBelow(diag([-1 1]), [0; 0]));

            % f = x*y has Q = [0 1;1 0] and d'Qd = 2*d1*d2, negative throughout the wedge's own
            % recession cone {d1<=0, d2>=0}. This is what makes "the envelope over an unbounded
            % face can be -inf" concrete. (Here Q's negative-eigenvalue direction (-1,1) is itself
            % a recession direction, so the eigenvector form of the test already catches it.)
            testCase.verifyTrue(wedge.quadUnboundedBelow([0 1; 1 0], [0; 0]));

            % But the eigenvector form is only SUFFICIENT, which is why the implementation
            % minimizes d'Qd over the whole recession cone instead. Take a region whose recession
            % cone is the single ray through (1,2):
            %     0 <= 2x - y + 1,  2x - y <= 0,  x >= 0   =>   {d : d2 = 2*d1, d1 >= 0}
            % With Q = diag(1,-1) neither eigenvector (1,0) nor (0,1) is a recession direction,
            % yet d = (1,2)/sqrt(5) gives d'Qd = (1-4)/5 < 0, so the value IS -inf.
            ray = region([2*x - y, -2*x + y - 1, -x], [x y]);
            testCase.verifyTrue(ray.quadUnboundedBelow(diag([1 -1]), [0; 0]));
            % ... and the same region is bounded below for a Q that grows along that ray.
            testCase.verifyFalse(ray.quadUnboundedBelow(eye(2), [0; 0]));

            % Q == 0: purely affine, so the LINEAR slope decides. On {y>=0} the recession cone
            % contains (1,0) and (-1,0), so any nonzero x-component is unbounded below, while
            % f = y is bounded below (slope 0 along those, +1 along (0,1)).
            testCase.verifyTrue(halfpl.quadUnboundedBelow(zeros(2), [1; 0]));
            testCase.verifyFalse(halfpl.quadUnboundedBelow(zeros(2), [0; 1]));
        end

        function quadUnboundedBelowRefusesACurvedFacet(testCase)
        % Dropping a non-affine facet enlarges the region and so enlarges the recession cone,
        % which is unsound in the one direction that matters: it could certify -inf for a
        % direction the true region never recedes along. Refuse rather than approximate.
            x = sym('x'); y = sym('y');
            r = region([x^2 - y, y - 4], [x y]);           % x^2 <= y <= 4 : a parabolic facet
            testCase.verifyError(@() r.quadUnboundedBelow(eye(2), [0; 0]), ...
                'region:quadUnboundedBelow:nonAffineFacet');
        end

        function recessionRaysClassifiesTheConeAndReturnsItsExtremeRays(testCase)
        % The directions must come from the INEQUALITIES. An unbounded region stores its
        % directions as +/-intmax vertices, and Step 1 reads vertices as coordinates -- that is
        % how an envelope ends up carrying 2147483647 and intmax^2 = 4611686014132420609.
            x = sym('x'); y = sym('y');

            % A pointed 2-dimensional cone: the second quadrant recedes along -x and +y.
            [D, kind] = region([x, -y], [x y]).recessionRays;
            testCase.verifyEqual(kind, 'wedge');
            testCase.verifyEqual(size(D,1), 2);
            testCase.verifyTrue(any(all(abs(D - [-1 0]) < 1e-12, 2)));
            testCase.verifyTrue(any(all(abs(D - [ 0 1]) < 1e-12, 2)));

            % A half-strip recedes along one ray only.
            [D, kind] = region([y, -x-1, x-1], [x y]).recessionRays;   % y<=0, -1<=x<=1
            testCase.verifyEqual(kind, 'ray');
            testCase.verifyEqual(size(D,1), 1);
            testCase.verifyLessThan(norm(D(1,:) - [0 -1]), 1e-12);

            % A bounded region recedes along nothing at all.
            [D, kind] = region([-x, -y, x+y-2], [x y]).recessionRays;
            testCase.verifyEqual(kind, 'bounded');
            testCase.verifyEmpty(D);

            % A half-plane's recession cone contains a line, so there is no apex to fan from.
            [~, kind] = region([-y], [x y]).recessionRays;
            testCase.verifyEqual(kind, 'nonpointed');
        end

        function recessionRaysReturnsExACTdirectionsNotTrigRoundTrips(testCase)
        % REGRESSION. The directions were once rebuilt as (cos t, sin t) from the candidate
        % angle. That put 6.123e-17 where a 0 belonged, and since these directions go on to
        % BUILD the sub-face half-planes in fanUnboundedFace, `x <= 0` came back as
        % `x - 4967757600021511/81129638414606681695789005144064*y <= 0` -- a different, and no
        % longer pointed, half-plane. A direction of a rational half-plane is rational; keep it.
            x = sym('x'); y = sym('y');
            D = region([x, -y], [x y]).recessionRays;
            testCase.verifyEqual(D(abs(D) < 0.5), zeros(2,1), ...
                'a zero component must be exactly zero, not 6e-17');
        end
    end
end
