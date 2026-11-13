classdef regionUnitTest < matlab.unittest.TestCase
% Method-level unit tests for `region`, the largest file in the toolbox (5936 lines) and the one
% whose coverage depends most on the symbolic pipeline: 37.7% from the fast bucket alone, 56%
% once the normal bucket runs, and the rest only when the slow bucket does.
%
% HOW THIS DIFFERS FROM `regionTest` AND `testRegion`. Those exercise region through the
% operations the pipeline performs -- merge chains, simplification, envelope gates. This calls one
% method at a time, on the smallest fixture that method is defined for, and asserts that method's
% own contract. A pipeline test that happens to reach a helper is not coverage of it: when it goes
% red it names the pipeline, and when the helper is subtly wrong it usually stays green.
%
% EVERY ASSERTION IS THE DEFINITION. A region is the point set {p : g_k(p) <= 0}, so the checks
% sample it (plqCheck.inRegion and friends); the interval routines are checked against dense
% membership; the predicates are checked against what they promise, in the direction they promise
% it. Nothing pins a constraint list or a vertex order.
%
% BUCKET: fast. Affine fixtures only, so no `solve` on a conic and no symbolic blow-up.

    properties (Constant)
        X = sym('x')
        Y = sym('y')
    end

    methods (Static)
        function v = vars()
            v = [regionUnitTest.X, regionUnitTest.Y];
        end

        function tf = inIntervals(ivs, u)
        % Membership in a set of [lo hi] rows -- the independent reading of what the interval
        % routines return.
            tf = false(size(u));
            for k = 1:size(ivs,1)
                tf = tf | (u >= ivs(k,1) & u <= ivs(k,2));
            end
        end
    end

    methods (Test)

        % =========================================================================================
        % quadIneqIntervals -- the admissible set of a scalar quadratic inequality
        % =========================================================================================
        function quadIneqIntervalsMatchesDirectEvaluation(testCase)
        % `a u^2 + b u + c <= 0` has five shapes, and the routine's header names the one that is
        % easy to get wrong: a DOWNWARD parabola excludes a middle interval, so the answer is two
        % rows, not one. All five are asserted the same way -- against evaluating the quadratic on
        % a dense grid, which knows nothing of the case split.
            cases = { [1 0 -4],   'upward, two roots: [-2,2]'
                      [1 0 1],    'upward, no real root: empty'
                      [1 -2 1],   'upward, double root: the single point u = 1'
                      [-1 0 4],   'downward: two unbounded rows, |u| >= 2'
                      [-1 0 -1],  'downward, no real root: the whole line'
                      [0 2 -6],   'linear, positive slope: u <= 3'
                      [0 -2 -6],  'linear, negative slope: u >= -3'
                      [0 0 -1],   'constant <= 0: the whole line'
                      [0 0 1],    'constant > 0: empty' };
            u = linspace(-8, 8, 3001);
            for k = 1:size(cases,1)
                coef = cases{k,1}; nm = cases{k,2};
                ivs = region.quadIneqIntervals(coef);
                want = (coef(1)*u.^2 + coef(2)*u + coef(3)) <= 1e-12;
                got  = regionUnitTest.inIntervals(ivs, u);
                % Points within one grid step of a root are genuinely ambiguous at this
                % resolution, so they are excluded rather than given a tolerance that would also
                % hide a wrong interval.
                val = abs(coef(1)*u.^2 + coef(2)*u + coef(3));
                decisive = val > 1e-3;
                testCase.verifyEqual(got(decisive), want(decisive), sprintf( ...
                    'quadIneqIntervals(%s) -- %s -- disagrees with direct evaluation at %d of %d points', ...
                    mat2str(coef), nm, sum(got(decisive) ~= want(decisive)), sum(decisive)));
                % The rows must be genuinely disjoint and ordered, which the callers rely on.
                if size(ivs,1) > 1
                    testCase.verifyTrue(all(ivs(1:end-1,2) < ivs(2:end,1)), sprintf( ...
                        'quadIneqIntervals(%s): rows must be disjoint and sorted, got %s', ...
                        mat2str(coef), mat2str(ivs)));
                end
                testCase.verifyTrue(all(ivs(:,1) <= ivs(:,2)), sprintf( ...
                    'quadIneqIntervals(%s): every row must have lo <= hi', mat2str(coef)));
            end
        end

        function intersectIntervalUnionsIsTheSetIntersection(testCase)
        % Asserted as a set identity on a dense grid, over pairs chosen to cover disjoint,
        % nested, overlapping, touching and empty inputs -- and the unbounded rows the routine
        % above produces, since those are what it is actually fed.
            pairs = { [-2 2],              [0 5]
                      [-2 2],              [3 5]
                      [-inf -2; 2 inf],    [-3 3]
                      [-inf inf],          [1 4]
                      zeros(0,2),          [0 1]
                      [0 1; 2 3],          [0.5 2.5]
                      [0 1],               [1 2] };
            u = linspace(-8, 8, 3001);
            for k = 1:size(pairs,1)
                A = pairs{k,1}; B = pairs{k,2};
                out = region.intersectIntervalUnions(A, B);
                want = regionUnitTest.inIntervals(A, u) & regionUnitTest.inIntervals(B, u);
                got  = regionUnitTest.inIntervals(out, u);
                testCase.verifyEqual(got, want, sprintf( ...
                    'intersectIntervalUnions(%s, %s) = %s is not the set intersection', ...
                    mat2str(A), mat2str(B), mat2str(out)));
            end
        end

        % =========================================================================================
        % isFeasible -- a NECESSARY condition, so only one direction may be asserted
        % =========================================================================================
        function isFeasibleIsSoundInTheDirectionItClaims(testCase)
        % Its own header: "false means provably empty; true means nothing was found to contradict
        % it". So `false` on a NON-empty region is a defect and is asserted; `true` on an empty one
        % is permitted and is not. The infeasible fixture is the one the header names, because
        % pairwise checking cannot see it: x >= 1, y >= 1, x + y <= 1 is jointly infeasible with
        % every PAIR perfectly satisfiable.
            v = regionUnitTest.vars();
            feasible = { plqCheck.triRegion([0 0; 1 0; 0 1], v),           'unit triangle'
                         plqCheck.boxRegion([-1 -1], [2 3], v),            'box'
                         plqCheck.wedgeRegion([0 0], [1 0], [0 1], v),     'quadrant'
                         plqCheck.halfStripRegion([0 0], [2 0], [0 1], v), 'half-strip' };
            for k = 1:size(feasible,1)
                r = feasible{k,1}; nm = feasible{k,2};
                P = plqCheck.regionSample(r, plqCheck.regionBox(r), 21, 1e-7);
                testCase.assumeGreaterThan(size(P,1), 0, sprintf('%s: fixture has no interior', nm));
                testCase.verifyTrue(r.isFeasible(), sprintf( ...
                    ['%s: isFeasible said FALSE for a region with %d sampled interior points. ' ...
                     'False is documented to mean provably empty.'], nm, size(P,1)));
            end
            % The three-constraint trap from the header.
            trap = plqCheck.halfPlanes([-1 0; 0 -1; 1 1], [-1; -1; 1], v);
            testCase.verifyFalse(trap.isFeasible(), ...
                ['x >= 1, y >= 1, x + y <= 1 is jointly infeasible and must be reported so -- ' ...
                 'every PAIR of those three is satisfiable, which is why pairwise checking missed it']);
        end

        % =========================================================================================
        % intersection -- the point-set intersection, both directions
        % =========================================================================================
        function intersectionIsThePointSetIntersection(testCase)
            v = regionUnitTest.vars();
            A = plqCheck.boxRegion([0 0], [3 3], v);
            B = plqCheck.boxRegion([1 -1], [5 2], v);
            C = A.intersection(B);
            box = plqCheck.regionBox([A B]);
            P = plqCheck.regionSample(C, box, 41, 1e-7);
            testCase.verifyGreaterThan(size(P,1), 0, 'the intersection has no interior to sample');
            testCase.verifyTrue(all(plqCheck.inRegion(A, P, 1e-6) & plqCheck.inRegion(B, P, 1e-6)), ...
                'a point of the intersection lies outside one of the operands');
            % ...and nothing was lost: a point in both operands must be in the intersection.
            Q = plqCheck.regionSample(A, box, 41, 1e-7);
            both = Q(plqCheck.inRegion(B, Q, -1e-7), :);
            if ~isempty(both)
                inC = plqCheck.inRegion(C, both, 1e-6);
                testCase.verifyTrue(all(inC), sprintf( ...
                    'intersection lost %d of %d points that lie in both operands', sum(~inC), numel(inC)));
            end
        end

        % =========================================================================================
        % minus -- the cheap half. The full set-difference check is `regionMinusTest`, in the
        % NORMAL bucket: measured 143 s for that one test against 18 s for the other nine here,
        % because region.minus does symbolic work per facet pair. Splitting keeps this suite
        % inside the fast bucket's budget without weakening either assertion.
        % =========================================================================================
        function minusOfARegionWithItselfIsEmpty(testCase)
        % The early return in `minus`, and the only case where the answer is exactly nothing.
            v = regionUnitTest.vars();
            A = plqCheck.triRegion([0 0; 2 0; 0 2], v);
            testCase.verifyEmpty(A - A, 'A \ A must be empty');
        end

        % =========================================================================================
        % deleteIneq -- dropping redundant facets must not move the boundary
        % =========================================================================================
        function deleteIneqDropsOnlyRedundantFacets(testCase)
        % The unit triangle plus two constraints it already implies. Whatever the routine chooses
        % to remove, the point set must not change and the list must not grow -- the two halves of
        % "redundant".
            v = regionUnitTest.vars();
            base = plqCheck.triRegion([0 0; 1 0; 0 1], v);
            g = sym.empty(1,0);
            for k = 1:size(base.ineqs,2), g(k) = base.ineqs(k).f; end
            g(end+1) = v(1) - 5;                       % x <= 5, implied
            g(end+1) = v(2) - 5;                       % y <= 5, implied
            padded = region(g, v);

            [lelim, trimmed] = padded.deleteIneq(v);
            testCase.verifyClass(lelim, 'logical', 'deleteIneq reports whether it removed anything');
            testCase.verifyLessThanOrEqual(size(trimmed.ineqs,2), size(padded.ineqs,2), ...
                'deleteIneq must never add a constraint');
            plqCheck.verifyRegionsAgree(testCase, padded, trimmed, ...
                'deleteIneq changed the point set');
        end

        % =========================================================================================
        % eq2 -- an equality predicate is only useful if it separates
        % =========================================================================================
        function eq2SeparatesEqualFromUnequalRegions(testCase)
            v = regionUnitTest.vars();
            A = plqCheck.triRegion([0 0; 1 0; 0 1], v);
            A2 = plqCheck.triRegion([0 0; 1 0; 0 1], v);
            B = plqCheck.triRegion([0 0; 2 0; 0 1], v);
            testCase.verifyTrue(A.eq2(A), 'eq2 must be reflexive');
            testCase.verifyTrue(A.eq2(A2), 'two identically built regions must compare equal');
            testCase.verifyFalse(A.eq2(B), 'a different triangle must not compare equal');
        end

        % =========================================================================================
        % getFeasiblePtNearV -- the contract is in the name: FEASIBLE, and NEAR that vertex
        % =========================================================================================
        function getFeasiblePtNearVReturnsAFeasiblePoint(testCase)
        % It bisects the two edge directions at vertex i, so the answer must satisfy every
        % constraint and must be closer to vertex i than to any other. Both are asserted; neither
        % pins a distance, which is an implementation choice.
            v = regionUnitTest.vars();
            r = plqCheck.triRegion([0 0; 4 0; 0 3], v);
            V = plqCheck.finiteVertexSet(r);
            testCase.assumeEqual(size(V,1), 3, 'the fixture must have three vertices');
            for i = 1:r.nv
                vi = [double(r.vx(i)), double(r.vy(i))];
                if ~all(isfinite(vi)), continue, end
                [tx, ty] = r.getFeasiblePtNearV(i);
                p = [double(tx), double(ty)];
                testCase.verifyTrue(all(isfinite(p)), sprintf( ...
                    'getFeasiblePtNearV(%d) returned %s', i, mat2str(p)));
                testCase.verifyTrue(plqCheck.inRegion(r, p, 1e-7), sprintf( ...
                    'getFeasiblePtNearV(%d) returned (%g,%g), which is not in the region', ...
                    i, p(1), p(2)));
                d = vecnorm(V - p, 2, 2);
                [~, nearest] = min(d);
                testCase.verifyEqual(norm(V(nearest,:) - vi), 0, 'AbsTol', 1e-9, sprintf( ...
                    'getFeasiblePtNearV(%d) returned a point nearer another vertex', i));
            end
        end

        % =========================================================================================
        % certifiesNonPositive -- sound one way only, which is exactly how it must be tested
        % =========================================================================================
        function certifiesNonPositiveIsSoundWhenItSaysYes(testCase)
        % A `true` verdict is a certificate that the form is <= 0 everywhere on the region, so a
        % true verdict with a positive sampled value is a defect. A `false` verdict asserts
        % NOTHING -- the routine is documented to refuse whatever it cannot prove -- so this test
        % deliberately makes no claim about it, and instead asserts that at least one fixture is
        % certified, or the check would be vacuous.
            v = regionUnitTest.vars();
            r = plqCheck.boxRegion([0 0], [1 1], v);
            forms = { v(1) - 2, v(2) - 2, v(1) + v(2) - 5, -v(1), v(1) - 0.5 };
            box = plqCheck.regionBox(r);
            P = plqCheck.regionSample(r, box, 31, 1e-7);
            testCase.assumeGreaterThan(size(P,1), 0, 'nothing sampled');
            nCertified = 0;
            for k = 1:numel(forms)
                f = symbolicFunction(forms{k});
                [ok, ~] = r.certifiesNonPositive(f);
                if ~ok, continue, end
                nCertified = nCertified + 1;
                h = plqCheck.handleFor(forms{k}, v);
                val = arrayfun(@(i) double(h(P(i,1), P(i,2))), (1:size(P,1))');
                testCase.verifyLessThanOrEqual(max(val), 1e-7, sprintf( ...
                    ['certifiesNonPositive said YES for %s on the unit box, but it reaches %.3g ' ...
                     'there -- a certificate that is not true'], char(forms{k}), max(val)));
            end
            testCase.verifyGreaterThan(nCertified, 0, ...
                'no form was certified, so the soundness check above never ran');
        end
    end
end
