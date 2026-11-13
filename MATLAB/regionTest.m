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
            % VECTORISED, and the resolution matters: the oracle can only see a sliver of the
            % cone as wide as the angular step it samples the boundary with. At 1 degree it
            % reported a direction 0.96 degrees OUTSIDE a cone as inside (measured, on the
            % unbounded parabola strip), because the nearest sampled boundary direction landed
            % exactly on the degenerate value. One handle per constraint, evaluated on the whole
            % point cloud at once, is what makes a 0.2 degree step affordable.
            gh = cell(1, size(r.ineqs,2));
            for k = 1:size(r.ineqs,2)
                gh{k} = matlabFunction(r.ineqs(k).f, 'Vars', {x, y});
            end
            nd = 72; ang = (0:nd-1)*(2*pi/nd); U = [cos(ang)', sin(ang)']; epsr = 0.05;
            per = zeros(1, r.nv);
            for j = 1:r.nv
                v = [double(r.vx(j)), double(r.vy(j))];
                th = (0:0.2:359.8)'*pi/180;
                W = [];
                for rr = epsr*(0.2:0.2:1.0)
                    W = [W; v + rr*[cos(th), sin(th)]];   %#ok<AGROW>
                end
                keep = true(size(W,1),1);
                for k = 1:numel(gh)
                    keep = keep & (gh{k}(W(:,1), W(:,2)) <= 1e-12);
                end
                D = (W(keep,:) - v)/epsr;

                inCode = true(nd,1); onBdry = false(nd,1);
                for k = 1:size(NC,2)
                    e = NC(j,k);
                    if isAlways(e == 0, 'Unknown', 'false'), continue, end
                    c0 = double(subs(e, [s1 s2], [0 0]));
                    L = zeros(nd,1);
                    for i = 1:nd, L(i) = double(subs(e, [s1 s2], U(i,:))) - c0; end
                    inCode = inCode & (L <= 1e-9);
                    onBdry = onBdry | (abs(L) <= 1e-9);
                end
                m = zeros(nd,1);
                for i = 1:nd, m(i) = max(D*U(i,:)'); end
                % A direction ON THE CONE'S OWN BOUNDARY is never a clear disagreement, and this
                % is a mathematical exclusion, not a numerical fudge: where the region sits on the
                % CONCAVE side of a conic the exact normal cone is NOT CLOSED -- piece 9's vertex
                % 2 is locally {x <= y^2/4}, and (1,0), the perpendicular to the conic's tangent,
                % has points of the region strictly ahead of it (x = y^2/8 > 0) even though every
                % direction just inside does not. The routine returns the closed cone; the two
                % differ by one ray, which is exactly what the cell decomposition does not
                % distinguish (adjacent conjugate cells share their boundaries anyway).
                inTrue = m <= 1e-9; amb = ((m > 1e-9) & (m < 1e-3)) | onBdry;
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
        % RE-PINNED 2026-08-18 when the probe sites were rewritten to pick a root by FEASIBILITY
        % instead of by position (region.probeOnConstraint). Three rows changed ORIENTATION --
        % piece 9's (3,2) and both of the half-lens's first column -- because a different, and
        % feasible, probe now decides the sign. That is the improvement, not a regression: scored
        % against the definition the same three regions went 32 -> 29, 43 -> 29 and 5 -> 5 wrong
        % directions, and WITH the edge list they are still exactly 0
        % (vertexConesMatchTheDefinition, which is the contract; this test is only a pin).
        %
        % The values below are the output of the suite-green implementation, not a
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
                     '25/24 - s_2 - (2*s_1)/3', 's_2 - 2*s_1 - 3/8'}, 'piece 9'; ...
                region([(x+y)^2 - 4*x, -x, -y], [x y]), ...
                    {'s_2', 's_2'; 's_2', '4 - s_1'}, 'half lens'; ...
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

        function theSlotFallbackIsRightOnTheUnboundedLayout(testCase)
        % THE OTHER HALF OF THE SPECIFICATION. Called without an edge list, getNormalConeVertexQ
        % pairs vertex j with constraints j and j+1. That is the layout `getEdgeNosInf` builds for
        % an UNBOUNDED region -- slot 1 holds the constraint carrying the ray at vertex 1, and
        % edge j (vertex j to vertex j+1) lands in slot j+1 -- so the arriving element at vertex j
        % is slot j and the leaving one is slot j+1, exactly the pair the fallback takes. On a
        % BOUNDED region edge j is slot j and the same pair is off by one, which is why the
        % characterization test's three bounded fixtures come out wrong without eIdx.
        %
        % So the fallback is built here the way the pipeline builds it -- removeInfV,
        % poly2orderUnbounded, then the getEdgeNosInf scatter, as conjugateOfPiecePoly does --
        % rather than by hand, because it is the LAYOUT that is under test.
        %
        % MEASURED 2026-08-18: exact, 0 of 72 directions wrong at both vertices. The path is live
        % (the caller reaches it whenever edgeIndexList refuses, which it does for every unbounded
        % region), and this is what says it is sound there.
            x = sym('x'); y = sym('y'); s1 = sym('s_1'); s2 = sym('s_2');

            % {y >= x^2, -2 <= x <= 2}: two finite vertices, two upward rays, a curved edge
            % between them, and nIneq = nv + 1 -- the unbounded layout's own shape.
            d = region([x^2 - y, x - 2, -x - 2], [x y]);
            d = d.removeInfV;
            d = d.poly2orderUnbounded;
            edgeNo = d.getEdgeNosInf(d.vars);
            old = d.ineqs;
            for k = 1:numel(edgeNo)
                if edgeNo(k) > 0, d.ineqs(edgeNo(k)) = old(k); end
            end
            testCase.verifyEqual(char(d.ineqs(1).f), '- x - 2', ...
                'precondition: slot 1 carries the ray at vertex 1');

            NC = d.getNormalConeVertexQ(s1, s2);
            [bad, per] = regionTest.coneVsDefinition(d, NC, s1, s2);
            testCase.verifyEqual(bad, 0, sprintf( ...
                'the slot fallback disagrees with the definition on %d directions [%s]', ...
                bad, num2str(per)));
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

        function signOnRegionDecidesAffineFactorsAndRefusesASignChange(testCase)
        % FAST. The certificate the rational cross-face max rests on: a denominator may be
        % cleared only where its sign is CONSTANT on the cell, and "unknown" must never be read
        % as "positive". Products are decided factor by factor, which is the shape that actually
        % arrives -- the difference of two rational conjugates has denominator q1*q2.
            x = sym('x'); y = sym('y');
            r = region([-x, -y, x+y-1], [x y]);              % the unit simplex

            testCase.verifyEqual(r.signOnRegion(x + y + 1),  1);   % >= 1 on the simplex
            testCase.verifyEqual(r.signOnRegion(2 - x - y), 1);    % >= 1
            testCase.verifyEqual(r.signOnRegion(x + y - 2), -1);
            testCase.verifyEqual(r.signOnRegion((x+y+1)*(x+y-2)), -1, ...
                'a product is the product of its factors'' signs');
            testCase.verifyEqual(r.signOnRegion(sym(-3)), -1);

            % REFUSALS. x - y changes sign on the simplex; x + y - 1 VANISHES on a whole facet,
            % and touching zero is exactly what clearing a denominator may not do.
            testCase.verifyEqual(r.signOnRegion(x - y), 0);
            testCase.verifyEqual(r.signOnRegion(x + y - 1), 0);
            testCase.verifyEqual(r.signOnRegion(x*y - y*x + 0), 0, 'the zero expression');
        end

        function maxArrayRefusesATiedRationalPairItCanSplit(testCase)
        % THE DEFECT (DECISIONS 2026-08-19 night, later; biconjugateTest lines 246-251). On a
        % cell where every VERTEX ties, maxArray refuses to guess -- but the refusal was
        % restricted to a POLYNOMIAL pair, because splitting on f1 = f2 needs a polynomial
        % constraint. Every second-pass conjugate is rational, so on the biconjugate path the
        % refusal never applied and the vertex comparison picked f2 for no better reason than
        % operand order. Here f1 wins on half the simplex and f2 on the other half, so ANY
        % single winner is wrong.
            x = sym('x'); y = sym('y');
            r = region([-x, -y, x+y-1], [x y]);
            f2 = symbolicFunction(x^2, x+y+1);
            f1 = symbolicFunction(x^2 + x*y*(x-y), x+y+1);     % f1 - f2 = x*y*(x-y)/(x+y+1)

            testCase.verifyFalse(f1.isPolynomial);
            testCase.verifyFalse(f2.isPolynomial);
            for v = [0 0; 1 0; 0 1]'
                testCase.verifyEqual(double(subs(f1.f - f2.f, [x y], v')), 0, 'AbsTol', 1e-14, ...
                    'the fixture must TIE at every vertex');
            end

            [l, ~, ~] = r.maxArray(f1, f2);
            testCase.verifyFalse(l, ...
                'neither dominates on the simplex, so maxArray must refuse and let the caller split');
        end

        function splitmax3ClearsDenominatorsSoTheSplitIsRepresentable(testCase)
        % The other half of the same defect: once maxArray refuses, the caller splits on
        % f1 = f2, and region() takes only POLYNOMIAL constraints -- normalize1 raises
        % symbolic:coeffs:NotAPolynomial on a rational one. Clearing a denominator whose sign is
        % certified turns the split into the polynomial {x*y*(x-y) >= 0} / {<= 0}, and the
        % ORIENTATION must still be right: the first half is where f1 wins.
            x = sym('x'); y = sym('y');
            r = region([-x, -y, x+y-1], [x y]);
            f2 = symbolicFunction(x^2, x+y+1);
            f1 = symbolicFunction(x^2 + x*y*(x-y), x+y+1);

            ineqs = r.splitmax3(f1, f2);
            testCase.verifyTrue(ineqs(1).isPolynomial);
            testCase.verifyTrue(ineqs(2).isPolynomial);

            % region() must accept them -- this is the raise the defect note names.
            keep = sym.empty;
            for k = 1:size(r.ineqs,2), keep(k) = r.ineqs(k).f; end
            d1 = region([keep, ineqs(1).f], [x y]);
            d2 = region([keep, ineqs(2).f], [x y]);
            testCase.verifyNotEmpty(d1);
            testCase.verifyNotEmpty(d2);

            % The halves carry f1 and f2 in that order (maximumP's convention), so the first
            % half must be exactly where f1 - f2 >= 0.
            pF1 = [0.5 0.1];        % x > y: f1 wins
            pF2 = [0.1 0.5];        % x < y: f2 wins
            testCase.verifyGreaterThan(double(subs(f1.f - f2.f, [x y], pF1)), 0);
            testCase.verifyLessThan   (double(subs(f1.f - f2.f, [x y], pF2)), 0);
            testCase.verifyLessThanOrEqual(double(subs(ineqs(1).f, [x y], pF1)), 0, ...
                'the f1 half must admit a point where f1 wins');
            testCase.verifyGreaterThan(double(subs(ineqs(1).f, [x y], pF2)), 0, ...
                'and must exclude a point where f2 wins');
            testCase.verifyLessThanOrEqual(double(subs(ineqs(2).f, [x y], pF2)), 0);
            testCase.verifyGreaterThan(double(subs(ineqs(2).f, [x y], pF1)), 0);
        end

        function clearedDifferenceRefusesWhenTheDenominatorChangesSign(testCase)
        % SOUNDNESS. Clearing p1/q1 - p2/q2 to p1*q2 - p2*q1 flips the inequality wherever the
        % denominator does, so a denominator that is not sign-definite on the cell must be a
        % REFUSAL, not a silent multiplication. Refusing leaves today's behaviour in place; it
        % never produces a wrong split.
            x = sym('x'); y = sym('y');
            r = region([-x, -y, x+y-1], [x y]);
            f1 = symbolicFunction(x^2, x-y);                 % pole ON the simplex
            f2 = symbolicFunction(sym(0));
            [ok, ~, why] = r.clearedDifference(f1, f2);
            testCase.verifyFalse(ok);
            testCase.verifyEqual(why, 'denominatorSignUnknown');

            % ... and a polynomial pair needs no clearing at all.
            [ok2, ineq2, why2] = r.clearedDifference(symbolicFunction(x^2), symbolicFunction(y));
            testCase.verifyTrue(ok2);
            testCase.verifyEqual(why2, 'alreadyPolynomial');
            testCase.verifyEqual(simplify(ineq2.f - (x^2 - y)), sym(0));
        end


        function signEverywhereProvesAPerfectSquareInExpandedForm(testCase)
        % THE DEFECT, measured on the A.4 split. Two conjugate cells that meet TANGENTIALLY have
        % a perfect-square difference -- here f1 - f2 = ((s1+s2) - 1)^2/4, an edge cell's
        % quadratic against the affine cell it touches along s1 + s2 = 1. maxArray asked
        % isAlways about that difference as simplifyFraction leaves it, which is EXPANDED, and
        % MuPAD does not complete the square:
        %     isAlways((s1*s2)/2 - s2/2 - s1/2 + s1^2/4 + s2^2/4 + 1/4 >= 0)   UNKNOWN
        %     isAlways(((s1+s2) - 1)^2/4 >= 0)                                  TRUE
        % Same function, opposite verdicts. The undecided answer makes the caller SPLIT the cell
        % on f1 = f2 -- and that split boundary is `-(s1+s2-1)^2 <= 0`, which is VACUOUS. That
        % is how the constraint in aDegenerateConicPairStillYieldsItsVertex below came to exist.
        %
        % A cell that must not be split, split anyway, is also the Step 3 cost problem: every
        % undecided comparison is a cell that never merges back.
            s1 = sym('s_1'); s2 = sym('s_2');
            sq = (s1*s2)/2 - s2/2 - s1/2 + s1^2/4 + s2^2/4 + 1/4;   % = ((s1+s2)-1)^2/4
            testCase.verifyEqual(double(subs(sq, [s1 s2], [0.4 0.2])), 0.04, 'AbsTol', 1e-12, ...
                'the fixture must be the square, and strictly positive off the line');

            testCase.verifyEqual(region.signEverywhere(sq),  1);
            testCase.verifyEqual(region.signEverywhere(-sq), -1);

            % REFUSALS. An indefinite difference has no sign -- this one is maxArray's own
            % recorded example, the 4-cone fan's s2^2/4 - s1^2/2, which changes sign inside the
            % quadrant and must NOT be decided.
            testCase.verifyEqual(region.signEverywhere(s2^2/4 - s1^2/2), 0);
            testCase.verifyEqual(region.signEverywhere(s1 + s2), 0);

            % And the variables are treated as REAL: over the complex numbers -s1^2 <= 0 is not
            % true, and a version of this test that omitted the substitution decided nothing.
            testCase.verifyEqual(region.signEverywhere(-s1^2), -1);
        end

        function aDegenerateConicPairStillYieldsItsVertex(testCase)
        % REGRESSION, from the A.4 split. A line meets the degenerate conic -(s1+s2-1)^2 = 0 in a
        % DOUBLE root, so region.lineMeetsConicSym's radical is an unevaluated form of sqrt(0):
        %     t1 = 3 - 2*sqrt(2) - 2*((sqrt(2)-3/2)^2 - (sqrt(2)-2)^2 - sqrt(2) + 7/4)^(1/2)
        % Its VALUE is 3 - 2*sqrt(2) = 0.1715728752538099, and `double` reads it without
        % complaint -- but `simplify` raises symbolic:kernel:DivisionByZero on it, and
        % getVertices called simplify on the way to STORING the coordinate. The whole region
        % constructor then failed, from inside maximumP, with nothing pointing at the cause.
        %
        % simplify NORMALISES a coordinate here; it does not decide whether the vertex exists.
            s1 = sym('s_1'); s2 = sym('s_2');
            r = region([s1 + s2/2 + sqrt(sym(2)) - 2, -(s1+s2-1)^2], [s1 s2]);
            testCase.verifyNotEmpty(r);
            testCase.verifyGreaterThanOrEqual(r.nv, 1, 'the intersection point is a vertex');
            hit = false;
            for k = 1:r.nv
                if abs(double(r.vx(k)) - (3 - 2*sqrt(2))) < 1e-12 && ...
                   abs(double(r.vy(k)) - (2*sqrt(2) - 2)) < 1e-12
                    hit = true;
                end
            end
            testCase.verifyTrue(hit, 'the vertex at (3-2*sqrt(2), 2*sqrt(2)-2) must survive');
        end

        function tightenUnboundedFacetHandlesTWOCurvedFacets(testCase)
        % `tightenUnboundedFacet` used to abort unconditionally whenever a region had more than
        % one curved facet (`numel(qidx) ~= 1`), so `maxAffineOverRegion` fell back to whatever
        % the LINEAR-only relaxation said, with no chance to confirm or correct it using the
        % curved facets at all. Mechanism 1 (a straight recession ray receding every linear facet
        % AND every conic's own recession condition) never actually used that there is only one
        % conic, so it now runs for TWO as well (SCIP_READINESS.md task list, 2026-08-28).
        %
        % Region: {y >= x^2} intersect {y >= (x-1)^2} intersect {-5 <= x <= 5} -- two parabolas
        % open the same way, TWO curved facets, genuinely unbounded upward (no linear facet
        % bounds y at all). direction (0,1) recedes both parabolas (grad = (2x,-1) or
        % (2(x-1),-1), grad.(0,1) = -1 <= 0 for every x) and both linear facets (normals (1,0)
        % and (-1,0), dot 0). Cross-checked against a brute-force sample: feasible points exist
        % arbitrarily far out along y (built as its own oracle, CLAUDE.md sec 6).
            x = sym('x'); y = sym('y');
            r = region([x^2 - y, (x-1)^2 - y, x - 5, -x - 5], [x y]);
            [A, b, lin] = r.linearForm;
            testCase.verifyEqual(sum(~lin), 2, 'fixture must have exactly two curved facets');

            [val, st, ok] = region.tightenUnboundedFacet(r, [0 1], A, b, lin);
            testCase.verifyTrue(ok, 'mechanism 1 must fire: a straight ray recedes both conics');
            testCase.verifyEqual(st, 1);
            testCase.verifyEqual(val, inf);

            [valFull, stFull] = r.maxAffineOverRegion([0 1]);
            testCase.verifyEqual(stFull, 1);
            testCase.verifyEqual(valFull, inf);

            % The oracle: sample the true (non-relaxed) constraint set directly and confirm it
            % stays feasible at y far beyond any candidate finite bound.
            rng(7);
            hit = false;
            for k = 1:2000
                xx = -5 + 10*rand();
                yy = 1e8 * rand();
                if xx^2 <= yy && (xx-1)^2 <= yy
                    hit = true; break
                end
            end
            testCase.verifyTrue(hit, 'the true region must contain points at very large y');
        end

        function twoDifferentAxisConvexConicsAreProvablyBounded(testCase)
        % REGRESSION for a real bug, not a hypothetical: before this fix, `maxAffineOverRegion`
        % answered `Inf` on a region with NO linear facets at all (so the linear-only relaxation
        % is trivially "unbounded", region.maxLinear's own isempty(A) branch) bounded by two
        % CONVEX (PSD) parabolas on DIFFERENT axes -- {x>=y^2} n {y>=x^2}, a small lens between
        % (0,0) and (1,1). That answer is not conservative, it is WRONG: the true region is
        % bounded (a rank<=1 PSD conic can only be receded along its own null/axis direction --
        % see `tightenUnboundedFacet`'s mechanism-3 comment for the closed-form argument -- and
        % two DIFFERENT axes share no common receding direction, so neither facet's curvature
        % admits an escape).
        %
        % Cross-checked against a brute-force sample (built as its own oracle, CLAUDE.md sec 6)
        % on 6 directions, including one (`[-1 1]`) whose true maximiser is a SMOOTH tangency
        % point on one arc, `(0.25, 0.5)`, not a vertex -- confirming the fix reuses
        % `tightenBoundedFacet`'s candidate machinery correctly, not just the two known corners.
            x = sym('x'); y = sym('y');
            r = region([y^2 - x, x^2 - y], [x y]);
            [A, b, lin] = r.linearForm;
            testCase.verifyEqual(sum(lin), 0, 'fixture must have NO linear facets');
            testCase.verifyEqual(sum(~lin), 2, 'fixture must have exactly two curved facets');

            cases = {[-1 1], 0.25; [1 1], 2; [1 0], 1; [0 1], 1; [2 -1], 1; [1 -3], 1/12};
            rng(11);
            for k = 1:size(cases,1)
                c = cases{k,1}; want = cases{k,2};
                [val, st] = r.maxAffineOverRegion(c);
                testCase.verifyEqual(st, 0, sprintf('direction [%g %g] must be decided bounded', c));
                testCase.verifyEqual(double(val), want, 'AbsTol', 1e-9, sprintf( ...
                    'direction [%g %g]: got %s, want %g', c(1), c(2), char(val), want));

                best = -inf;
                for it = 1:200000
                    xx = rand(); yy = rand();
                    if xx >= yy^2 && yy >= xx^2
                        best = max(best, c(1)*xx + c(2)*yy);
                    end
                end
                testCase.verifyLessThanOrEqual(best, double(val) + 1e-6, sprintf( ...
                    'brute-force sample %.6f exceeds the computed bound %.6f for [%g %g]', ...
                    best, double(val), c(1), c(2)));
            end

            % A GENUINELY unbounded same-axis pair (both facets recede the SAME way) must still
            % report Inf -- mechanism 3 must not turn a real unboundedness into a false bound.
            rSameAxis = region([x^2-y, (x-1)^2-y, x-5, -x-5], [x y]);
            [valSA, stSA] = rSameAxis.maxAffineOverRegion([0 1]);
            testCase.verifyEqual(stSA, 1);
            testCase.verifyEqual(valSA, inf);

            % Mixed convex/concave must still abstain (mechanism 3 only fires when BOTH are
            % convex; curvature alone does not rule out an escape when one facet is concave).
            rMixed = region([y^2-x, y-x^2], [x y]);   % x>=y^2 (convex) AND y<=x^2 (concave)
            [valM, stM] = rMixed.maxAffineOverRegion([1 0]);
            testCase.verifyEqual(stM, 1, 'mixed convex/concave must still report unbounded');
            testCase.verifyEqual(valM, inf);
        end

        function sameAxisOppositeSenseConvexConicsAreProvablyBounded(testCase)
        % A SECOND, separate previously-latent bug, found while generalizing the different-axis
        % mechanism above: two CONVEX facets sharing the SAME axis, but receding in OPPOSITE
        % senses along it (one needs y->+inf, the other y->-inf), also bound the region --
        % mechanism 1's own candidate list only tested ONE sign of each conic's null direction
        % (`region.quadNullDirsNumeric`'s own convention, not a geometric necessity), so a facet
        % receding only in the UNTESTED sign was invisible to it. Fixed alongside this test by
        % adding the negated null direction to mechanism 1's candidate list too.
        %
        % `y >= x^2-1` recedes only at (0,1); `y <= 1-x^2` recedes only at (0,-1); no common
        % direction, so the true region is bounded: `x^2-1 <= y <= 1-x^2` forces `x^2<=1`, a
        % lens with true max(x)=1 and max(y)=1, both at smooth points, neither a vertex of the
        % OTHER facet. Before this fix, both directions answered `Inf`.
            x = sym('x'); y = sym('y');
            r = region([x^2-y-1, x^2+y-1], [x y]);
            [A, b, lin] = r.linearForm;
            testCase.verifyEqual(sum(lin), 0, 'fixture must have NO linear facets');

            [valX, stX] = r.maxAffineOverRegion([1 0]);
            testCase.verifyEqual(stX, 0);
            testCase.verifyEqual(double(valX), 1, 'AbsTol', 1e-9);

            [valY, stY] = r.maxAffineOverRegion([0 1]);
            testCase.verifyEqual(stY, 0);
            testCase.verifyEqual(double(valY), 1, 'AbsTol', 1e-9);

            rng(13);
            for c = {[1 0], [0 1], [1 1], [1 -1]}
                cc = c{1};
                [val, st] = r.maxAffineOverRegion(cc);
                testCase.verifyEqual(st, 0);
                best = -inf;
                for it = 1:200000
                    xx = -1 + 2*rand();
                    if xx^2 > 1, continue, end
                    lo = xx^2 - 1; hi = 1 - xx^2;
                    yy = lo + (hi-lo)*rand();
                    best = max(best, cc(1)*xx + cc(2)*yy);
                end
                testCase.verifyLessThanOrEqual(best, double(val) + 1e-6, sprintf( ...
                    'brute-force sample %.6f exceeds the computed bound %.6f for [%g %g]', ...
                    best, double(val), cc(1), cc(2)));
            end
        end
    end
end
