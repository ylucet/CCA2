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
    end

    methods (Test)
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
    end
end
