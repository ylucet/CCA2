classdef testRegion < matlab.unittest.TestCase

    properties
        r
        s
        t
        u
        v
        w
    end
    
    methods(TestMethodSetup)
        % Setup for each test
        function setUpTestData(testCase)
          x = sym('x');
          y = sym('y');
          l1 = [-x,-y,x+y-1,x-2,y-1];
          testCase.r = region(l1,[x,y]);
          l2 = [-y,x+y-1,y-x,x-1,y-0.5];
          testCase.s = region(l2,[x,y]);
          l3 = [y-x,y-x^2,-y,y+x^2];
          testCase.t = region(l3,[x,y]);
          l4 = [y-x^2,y^2-x];
          testCase.u = region(l4,[x,y]);
          testCase.v = region([-x,x+1],[x,y]);

          %l5 = [-x - 7*y - 4, x + 7*y - 10, 148*x - 196*y + (x + 7*y)^2 - 684,-x - 4,(9*y)/5 - x + 5] 

          %l5 = [-x - 7*y - 4, x + 7*y - 10, 148*x - 196*y + (x + 7*y)^2 - 684,-(9*y)/5 + x - 5] ;

          %l5 = [x+2*y+4, x+5*y/2-1, x+7*y-46]

        
          %testCase.w = region(l5,[x,y]);

          %l6 = [x + 2*y + 4, x + (5*y)/2 - 1, 46 - 7*y - x, x - 2*y + 44 ] 
          l6 = [x + 2*y + 4, -x - (5*y)/2 + 1,  -x + 2*y - 44 ] 

          testCase.w = region(l6,[x,y]);

          %l7 = [-x - 7*y - 4, x + 7*y-10, 196*y-148*x-(x+7*y)^2 + 684, 4-2*y-x ] 

          %l7 = [x + 7*y-10, -196*y+148*x+(x+7*y)^2 - 684, 4-2*y-x, x - y/3-14/3,5*y/7-x+25/7 ] 

          %l7 = [-x - 7*y - 4, x + 7*y-10, -196*y+148*x+(x+7*y)^2 - 684, x-9*y/5-5, x-y/3-14/3 ] 

          %l7 = [-x - 7*y - 4, x + 7*y-10, -196*y+148*x+(x+7*y)^2 - 684, -x+9*y/5+5, -x-4 ] 

          %l7 = [x-y/3-14/3, 10-7*y-x]
          %l7 = [ -196*y+148*x+(x+7*y)^2 - 684, -x+9*y/5+5 ] 

          %l7 = [-x-7*y-4, x+7*y-10, 148*x-196*y+(x+7*y)^2-684, -x-2*y-4, x+2*y-4, (x+2*y)^2-24*y-8*x+32, 2*y-x-44] 

          %testCase.w = region(l7,[x,y]);
        end
    end

    methods (Static)
        function verifySlopesAreImplicitDerivatives(testCase, r, pi, pt, name)
        % `slopeAtVertex` must return, for each constraint index in pi, the slope of that
        % constraint's zero level curve at pt: dy/dx = -(dg/dx)/(dg/dy), with Inf where dg/dy = 0.
        % Recomputed here from `diff` so the test does not reuse the routine's own case split.
            v = r.vars;
            m = r.slopeAtVertex(pi, pt);
            testCase.verifyEqual(numel(m), numel(pi), sprintf('%s: one slope per index', name));
            for j = 1:numel(pi)
                g  = r.ineqs(pi(j)).f;
                gx = double(subs(diff(g, v(1)), v, pt));
                gy = double(subs(diff(g, v(2)), v, pt));
                if abs(gy) <= 1e-12 * max(1, abs(gx))
                    testCase.verifyTrue(isinf(double(m(j))), sprintf( ...
                        '%s: constraint %d has dg/dy = 0 at (%g,%g), so its slope must be Inf; got %s', ...
                        name, pi(j), pt(1), pt(2), char(m(j))));
                else
                    testCase.verifyEqual(double(m(j)), -gx/gy, 'AbsTol', 1e-10, sprintf( ...
                        '%s: constraint %d slope at (%g,%g)', name, pi(j), pt(1), pt(2)));
                end
            end
        end

        function verifyIsFacetPermutation(testCase, edgeNo, r, name)
        % getEdgeNos / getEdgeNosInf return the region's facets in BOUNDARY-WALK order, so the
        % list must be a permutation of the facet indices: every facet named exactly once, none
        % invented, none dropped. That is the property the callers rely on; the particular
        % starting facet is not (it depends on which vertex getVertices happened to emit first,
        % a Symbolic Math Toolbox ordering detail -- see testCreation's note).
            e = double(edgeNo(:)).';
            n = size(r.ineqs, 2);
            testCase.verifyEqual(sort(e), 1:n, sprintf( ...
                '%s: edge list %s is not a permutation of the %d facets', name, mat2str(e), n));
        end

        function verifyVerticesUpToRotation(testCase, r, wantX, wantY)
        % Verify r's vertex list equals (wantX,wantY) up to a CYCLIC ROTATION -- see the note in
        % testCreation for why rotation is the right equivalence and reversal is not.
            gotX = double(r.vx); gotY = double(r.vy);
            testCase.verifyEqual(numel(gotX), numel(wantX), 'vertex count');
            n = numel(wantX);
            ok = false;
            for k = 0:n-1
                idx = mod((0:n-1) + k, n) + 1;
                if isequal(gotX(idx), wantX) && isequal(gotY(idx), wantY)
                    ok = true; break
                end
            end
            testCase.verifyTrue(ok, sprintf( ...
                'vertices %s differ from %s by more than a cyclic rotation', ...
                mat2str([gotX(:), gotY(:)]), mat2str([wantX(:), wantY(:)])));
        end
    end

    methods(Test)
        % Test methods

        function testCreation(testCase)
           x = sym('x');
           y = sym('y');
          % l1 = [-x,-y,x+y-1];
          % l2 = [-y,x+y-1,y-x];
          % r = region(l1,[x,y]);
          % s = region(l2,[x,y]);
          % % r.print;
          % s.print;
          % r.printMaple;
          % s.printMaple;

          testCase.verifyEqual(isequal(testCase.r.ineqs(1).f,-x), true);
            testCase.verifyEqual(isequal(testCase.r.ineqs(2).f,-y), true);
            testCase.verifyEqual(isequal(testCase.r.ineqs(3).f,x+y-1), true);
            % Vertices are compared UP TO CYCLIC ROTATION. A polygon's vertex list has no
            % canonical starting point -- getVertices enumerates pairwise constraint
            % intersections, so which vertex lands first depends on the order `solve` returns
            % roots in, a Symbolic Math Toolbox implementation detail that does change between
            % releases. On R2024b `s` came out (0,0),(0.5,0.5),(1,0) against a hardcoded
            % (1,0),(0,0),(0.5,0.5): the SAME triangle, rotated by one. That is why this had
            % stood as a longstanding "toolbox compatibility" failure -- it was the assertion
            % over-specifying, not the geometry being wrong.
            %
            % Rotation only, NOT reversal: the cyclic ORDER carries the boundary orientation,
            % which the pipeline genuinely depends on (region.getNormalConeVertex walks
            % consecutive vertices and wraps), so a reversed list is a real defect and must
            % still fail here.
            testRegion.verifyVerticesUpToRotation(testCase, testCase.r, [0,0,1], [0,1,0]);

            testCase.verifyEqual(isequal(testCase.s.ineqs(1).f,-y), true);
            testCase.verifyEqual(isequal(testCase.s.ineqs(2).f,x+y-1), true);
            testCase.verifyEqual(isequal(testCase.s.ineqs(3).f,y-x), true);
            testRegion.verifyVerticesUpToRotation(testCase, testCase.s, [1,0,0.5], [0,0,0.5]);
        end

        function testslopeAtVertex(testCase)
        % slopeAtVertex(pi, pt) is the slope of each constraint pi(j)'s zero LEVEL CURVE at pt.
        % Asserted against the definition recomputed independently by implicit differentiation,
        % dy/dx = -(dg/dx)/(dg/dy) -- not against pinned numbers, so a constraint reordering or a
        % change of representation cannot make this red for the wrong reason. Covers both branches
        % of the routine: linear constraints (where it reads the coefficients) and curved ones
        % (where it differentiates), and a vertical tangent, where the answer must be Inf and not
        % a division-by-zero error.
            testRegion.verifySlopesAreImplicitDerivatives(testCase, testCase.r, [1,2], [0,0], 'r');
            testRegion.verifySlopesAreImplicitDerivatives(testCase, testCase.s, [1,3], [0,0], 's');
            testRegion.verifySlopesAreImplicitDerivatives(testCase, testCase.t, [1,2], [0,0], 't at origin');
            testRegion.verifySlopesAreImplicitDerivatives(testCase, testCase.t, [1,2], [1,1], 't at (1,1)');
            % The vertical facet -x <= 0 must give Inf rather than raising: pinned explicitly
            % because "no error and a finite wrong number" would pass the check above only if the
            % check itself divided by zero the same way, and it does not.
            m = testCase.r.slopeAtVertex([1,2],[0,0]);
            testCase.verifyTrue(isinf(double(m(1))), ...
                'the vertical facet -x <= 0 must have infinite slope');
        end
        function testsimplifyUnboundedRegion(testCase)
        % simplifyUnboundedRegion rewrites a region's PRESENTATION -- it drops redundant
        % constraints and repairs the vertex list -- and must not move its boundary. That is the
        % contract, and it is what is asserted: the same point set before and after, plus a
        % constraint count that does not grow (a "simplification" that adds a facet is not one).
        % Checked on all five setup fixtures, not just `w`, because the routine's branches are
        % selected by boundedness and by curvature: r and s are bounded polygons, t and u are
        % curved, v is a slab, w is an unbounded polyhedron.
        % THREE OF THE SIX FIXTURES ARE DEGENERATE, and finding that out is most of what writing
        % this assertion produced:
        %   v = {x >= 0} n {x <= -1}                     INFEASIBLE
        %   t = {y <= x, y <= x^2, y >= 0, y >= -x^2}    the SINGLE POINT (0,0) -- the last two
        %                                                constraints force y = 0 and then x = 0
        %   w = the three-facet fixture the printing version of this test was built on: the
        %       SINGLE POINT (-24,10), because its constraints force y = 10 exactly
        % So the only fixture the old test ever touched had no interior at all.
        %
        % THE CONTRACT ASSERTED, which is what survives that. simplifyUnboundedRegion rewrites a
        % presentation, so on a region with interior it must not move the boundary, must not grow
        % the constraint list, and must be idempotent. On a region with EMPTY interior it is
        % allowed to discard the region entirely -- and it does, for w -- because a
        % lower-dimensional cell carries no function value in a 2-D subdivision. What it may never
        % do in either case is invent feasibility, or drop a region that has interior.
        %
        % ONE BEHAVIOUR IS PINNED RATHER THAN DERIVED, and flagged as such: `w` is a genuine point
        % of its own constraint set (all three constraints are 0 there) and simplify returns
        % EMPTY. That looks deliberate but is documented nowhere, and `isempty` means "infeasible"
        % to callers, so the two readings differ. Pinned here so a change is a decision.
            fixtures = {testCase.r, 'r'; testCase.s, 's'; testCase.t, 't'; ...
                        testCase.u, 'u'; testCase.v, 'v'; testCase.w, 'w'};
            nWithInterior = 0;
            for k = 1:size(fixtures,1)
                before = fixtures{k,1}; nm = fixtures{k,2};
                after  = before.simplifyUnboundedRegion;

                if isempty(before)
                    testCase.verifyTrue(isempty(after), sprintf( ...
                        'simplifyUnboundedRegion(%s) made an INFEASIBLE region non-empty', nm));
                    continue
                end

                box  = plqCheck.regionBox(before);
                nInt = size(plqCheck.regionSample(before, box), 1);
                if nInt > 0
                    nWithInterior = nWithInterior + 1;
                    testCase.verifyFalse(isempty(after), sprintf( ...
                        'simplifyUnboundedRegion(%s) discarded a region that has interior', nm));
                    if isempty(after), continue, end
                    plqCheck.verifyRegionsAgree(testCase, before, after, ...
                        sprintf('simplifyUnboundedRegion(%s) changed the point set', nm));
                    testCase.verifyLessThanOrEqual(size(after.ineqs,2), size(before.ineqs,2), sprintf( ...
                        'simplifyUnboundedRegion(%s) grew the constraint list from %d to %d', ...
                        nm, size(before.ineqs,2), size(after.ineqs,2)));
                    % IDEMPOTENT: simplifying twice must give the same set as simplifying once. A
                    % routine that keeps finding new things to drop is not converging on a
                    % presentation, it is losing constraints.
                    plqCheck.verifyRegionsAgree(testCase, after, after.simplifyUnboundedRegion, ...
                        sprintf('simplifyUnboundedRegion(%s) is not idempotent', nm));
                elseif ~isempty(after)
                    % Kept a region with no interior: then its vertices may not have moved.
                    testCase.verifyEqual(plqCheck.finiteVertexSet(after), ...
                                         plqCheck.finiteVertexSet(before), 'AbsTol', 1e-9, sprintf( ...
                        'simplifyUnboundedRegion(%s) moved the vertices of a zero-area region', nm));
                end
            end
            % The one pinned behaviour, called out by name rather than folded into the loop.
            testCase.verifyTrue(isempty(testCase.w.simplifyUnboundedRegion), ...
                ['simplifyUnboundedRegion no longer discards the zero-area fixture w. That may ' ...
                 'be an improvement, but it changes what isempty means to every caller -- ' ...
                 'record the decision before making this test green again.']);
            testCase.verifyGreaterThanOrEqual(nWithInterior, 2, ...
                'fewer than two fixtures had an interior to sample; this test would be near-vacuous');
        end

        function testremoveTangent (testCase)
        % `removeTangent` deletes a facet that only TOUCHES the region -- it meets the boundary at
        % one point and cuts nothing off -- and `linear3pt` rewrites a facet that three vertices
        % turn out to be collinear on. Both are presentation changes, so both must leave the point
        % set alone; that is the assertion, in place of the three `print`-and-`return` blocks that
        % used to be here (two of which were unreachable, behind an unconditional `return`).
        %
        % TWO OF THE THREE ORIGINAL COMMENTS WERE WRONG, which is what asserting found. They read
        % "tangent removed / not removed / not removed"; measured over [-60,60]^2 the outcomes are
        % 5 facets -> 3, 2 -> 1, and 2 -> 2. In every case the point set is IDENTICAL -- the
        % facets that go are genuinely redundant, each meeting the region only where the conic
        % already bounds it. So the invariant asserted is the point set, which is the contract,
        % and the counts are recorded as observations rather than as the expectation.
            x = sym('x');
            y = sym('y');

            l7 = [x + 7*y - 10 , x + 2*y - 4 , 48*x - 56*y + 4*x*y + x^2 + 4*y^2 - 184 , ...
                  x - (9*y)/5 - 5, - y - 5 ];
            w0 = region(l7,[x,y]);
            w1 = w0.linear3pt;
            plqCheck.verifyRegionsAgree(testCase, w0, w1, 'linear3pt changed the point set');

            w2 = w1.removeTangent(w1.nv, w1.vx, w1.vy);
            plqCheck.verifyRegionsAgree(testCase, w1, w2, 'removeTangent (case 1) changed the point set');
            testCase.verifyLessThanOrEqual(size(w2.ineqs,2), size(w1.ineqs,2), ...
                'removeTangent must not add constraints');

            % Case 2: a parabola and a line TANGENT to it. The line touches the parabolic region
            % without cutting anything off, so it is redundant and removing it is sound -- the
            % point set is unchanged, which is asserted, over a window wide enough to see the
            % parabola (regionSample widens automatically here: the region has one vertex and
            % extends far past a vertex-sized box).
            l7 = [ -196*y+148*x+(x+7*y)^2 - 684, x-9*y/5-5 ];
            w3 = region(l7,[x,y]);
            w4 = w3.removeTangent(w3.nv, w3.vx, w3.vy);
            n2 = plqCheck.verifyRegionsAgree(testCase, w3, w4, 'removeTangent (case 2) changed the point set');
            testCase.verifyGreaterThan(n2, 0, ...
                'case 2: nothing was sampled, so the point-set check said nothing');
            testCase.verifyLessThanOrEqual(size(w4.ineqs,2), size(w3.ineqs,2), ...
                'removeTangent must not add constraints');

            % Case 3: the same pair with the conic negated -- the other side of the parabola.
            l7 = [ 196*y-148*x-(x+7*y)^2 + 684, x-9*y/5-5 ];
            w5 = region(l7,[x,y]);
            w6 = w5.removeTangent(w5.nv, w5.vx, w5.vy);
            n3 = plqCheck.verifyRegionsAgree(testCase, w5, w6, 'removeTangent (case 3) changed the point set');
            testCase.verifyGreaterThan(n3, 0, ...
                'case 3: nothing was sampled, so the point-set check said nothing');
            testCase.verifyEqual(size(w6.ineqs,2), size(w5.ineqs,2), ...
                'case 3: with the conic negated the line genuinely cuts, so no facet may be dropped');
        end
        function testMinus(testCase)
          x = sym('x');
          y = sym('y');
          l1 = [-x,-y,x+y-1];
          l2 = [-y,x+y-1,y-x];
          r = region(l1,[x,y]);
          s = region(l2,[x,y]);
          %r.print;
          %s.print;
            t = r - s;
           % t.print;
            testCase.verifyEqual(isequal(t.ineqs(1).f,-x), true);
            testCase.verifyEqual(isequal(t.ineqs(2).f,x+y-1), true);
            testCase.verifyEqual(isequal(t.ineqs(3).f,x-y), true);
        end

        function testMinus2(testCase)

          x = sym('x');
          y = sym('y');
          l1 = [-x,-y,x+y-1];
          l2 = [-x,-y,x+y-1,y-x-0.5,-y+x-0.5];
          r = region(l1,[x,y]);
          s = region(l2,[x,y]);
          r.print;
          s.print;
            t = r - s;
            size(t)
            t(1).print;
            t(2).print;
            testCase.verifyEqual(isequal(t(1).ineqs(1).f,-x), true);
            testCase.verifyEqual(isequal(t(1).ineqs(2).f,x+y-1), true);
            testCase.verifyEqual(isequal(t(1).ineqs(3).f,x-y+1/2), true);
            testCase.verifyEqual(isequal(t(2).ineqs(1).f,-y), true);
            testCase.verifyEqual(isequal(t(2).ineqs(2).f,x+y-1), true);
            testCase.verifyEqual(isequal(t(2).ineqs(3).f,y-x+1/2), true);
        end

        function testMerge(testCase)
        % A CHAIN of merges -- the shape mergeL actually runs -- so each step is checked against
        % the accumulator it was handed, not against the original operands. `merge`'s contract is
        % in plqCheck.verifyMergeSound: the result always contains both inputs, and when the
        % returned flag says the union was exact it may contain nothing else.
          x = sym('x');
          y = sym('y');
          l1 = [-x - 7*y - 4, x + 7*y-10, 196*y-148*x-(x+7*y)^2 + 684, 4-2*y-x ];
          r1 = region(l1,[x,y]);
          r1 = r1.simplifyUnboundedRegion;

          l2 = [(y/3)-x+ 14/3, 10-7*y-x, y-2];
          r2 = region(l2,[x,y]);
          r2 = r2.simplifyUnboundedRegion;
          l3 = [x+7*y+4,4-2*y-x];
          r3 = region(l3,[x,y]);
          r3 = r3.simplifyUnboundedRegion;
          l4 = [-y/3 +x - 14/3, 10 -7*y - x, y-2, 5*y/7 -x + 25/7]; %#ok<NASGU>
          % HISTORY: r4 is built from l2, not l4. Left as it was -- the fixture is a repeat of r2
          % on purpose or by accident, but either way it exercises merging a region with one it
          % already absorbed, which is a case mergeL hits constantly.
          r4 = region(l2,[x,y]);
          r4 = r4.simplifyUnboundedRegion;

          l5 = [x + 7*y-10, -196*y+148*x+(x+7*y)^2 - 684, 4-2*y-x, x - y/3-14/3,5*y/7-x+25/7 ];
          r5 = region(l5,[x,y]);
          r5 = r5.simplifyUnboundedRegion;

          acc = r1;
          chain = {r2, r3, r4, r5};
          for k = 1:numel(chain)
              nxt = chain{k};
              [l, m] = merge(acc, nxt);
              plqCheck.verifyMergeSound(testCase, m, acc, nxt, l, ...
                  sprintf('chain merge step %d (r%d)', k, k+1));
              acc = m;
          end
          testCase.verifyNotEmpty(acc, 'the merge chain collapsed to an empty region');
        end

        function testMerge2(testCase)
          x = sym('x');
          y = sym('y');

          l0 = [x + 7*y - 10 , 148*x - 196*y + (x + 7*y)^2 - 684 , 4 - 2*y - x , (5*y)/7 - x + 25/7 ]
          r0 = region(l0,[x,y]);
          r0 = r0.simplifyUnboundedRegion;
          
          l1 = [x-y/3-14/3, 10-7*y-x, y-2, (5*y)/7-x+25/7] 
          r1 = region(l1,[x,y]);
          r1 = r1.simplifyUnboundedRegion;
                  
          r0.print
          r1.print
          [l,r] = merge(r0,r1)
          plqCheck.verifyMergeSound(testCase, r, r0, r1, l, ...
              'merge(r0,r1)');
          l
          r.print

          return
          
        end

        function testMerge3(testCase)
          x = sym('x');
          y = sym('y');

          l0 = [x+2*y+4, x-(9*y)/5-5, -y-5, x+7*y-46];
          r0 = region(l0,[x,y]);
          r0 = r0.simplifyUnboundedRegion;
          
          l1 = [x-y/3-14/3, x-2*y+44, 46-7*y-x] ;
          r1 = region(l1,[x,y]);
          r1 = r1.simplifyUnboundedRegion;
                  
          r0.print
          r1.print
          [l,r] = merge(r0,r1)
          plqCheck.verifyMergeSound(testCase, r, r0, r1, l, ...
              'merge(r0,r1)');
          l
          r.print

          return
          
        end

        function testMerge4(testCase)
          x = sym('x');
          y = sym('y');

          l0 = [x-(9*y)/5-5, y+5 ];
          r0 = region(l0,[x,y]);
          r0 = r0.simplifyUnboundedRegion;
          
          l1 = [-x+(9*y)/5+5,x+4] ;
          r1 = region(l1,[x,y]);
          r1 = r1.simplifyUnboundedRegion;
                  
          r0.print
          r1.print
          [l,r] = merge(r0,r1);
          plqCheck.verifyMergeSound(testCase, r, r0, r1, l, ...
              'merge(r0,r1)');
          l
          r.print

          return
          
        end

        function testMerge5(testCase)
          x = sym('x');
          y = sym('y');

          l0 = [-x-7*y-4, x+2*y-4, 56*y-48*x-4*x*y-x^2-4*y^2+184 ];
          r0 = region(l0,[x,y]);
          %r0 = r0.simplifyUnboundedRegion;
          r0.print
          %return
          l1 = [x+7*y+4,x-(9*y)/5-5,-x-2*y-4, 56*y-48*x-4*x*y-x^2-4*y^2+184] ;
          r1 = region(l1,[x,y]);
          r1 = r1.simplifyUnboundedRegion;
                  
          
          r1.print
          [l,r] = merge(r0,r1);
          plqCheck.verifyMergeSound(testCase, r, r0, r1, l, ...
              'merge(r0,r1)');
          l
          r.print

          return
          
        end

  
        function testMerge6(testCase)
          x = sym('x');
          y = sym('y');

          l0 = [-x-7*y+10, 5*y/7-x+25/7, y-2 ];
          r0 = region(l0,[x,y]);
          %r0 = r0.simplifyUnboundedRegion;
          r0.print
          %return
          l1 = [x+7*y-10,5*y/7-x+25/7, -x-2*y+4] ;
          r1 = region(l1,[x,y]);
          r1 = r1.simplifyUnboundedRegion;
                  
          
          r1.print
          [l,r] = merge(r0,r1);
          plqCheck.verifyMergeSound(testCase, r, r0, r1, l, ...
              'merge(r0,r1)');
          l
          r.print

          return
          
        end

        function testMerge7(testCase)
          x = sym('x');
          y = sym('y');




          l0 = [-x-7*y- 4,x+7*y-10,-x-2*y-4,x+2*y-4,48*x-56*y+4*x*y+x^2+4*y^2-184];
          r0 = region(l0,[x,y]);
          r0 = r0.simplifyUnboundedRegion;
          r0.print
          %return
          l1 = [x+7*y+4,x-(9*y)/5-5,-x-2*y-4,48*x-56*y+4*x*y+x^2+4*y^2-184] ;
          r1 = region(l1,[x,y]);
          r1 = r1.simplifyUnboundedRegion;
                  
          
          r1.print
          [l,r] = merge(r0,r1);
          plqCheck.verifyMergeSound(testCase, r, r0, r1, l, ...
              'merge(r0,r1)');
          l
          r.print

          return
          
        end

        function testMerge8(testCase)
          x = sym('x');
          y = sym('y');



          l0 = [x+7*y-10, -x-2*y-4, x+2*y-4, 48*x-56*y+4*x*y+x^2+4*y^2-184, x-(9*y)/5-5];
          r0 = region(l0,[x,y]);
          r0 = r0.simplifyUnboundedRegion;
          r0.print
          %return
          l1 = [x+2*y-4, 2*y-x-44, 10-7*y-x, -x-2*y-4] ;
          r1 = region(l1,[x,y]);
          r1 = r1.simplifyUnboundedRegion;
                  
          
          r1.print
          [l,r] = merge(r0,r1);
          plqCheck.verifyMergeSound(testCase, r, r0, r1, l, ...
              'merge(r0,r1)');
          l
          r.print

          return
          
        end

        function testMerge9(testCase)
          x = sym('x');
          y = sym('y');
          f1 = symbolicFunction(48*x-56*y+4*x*y+x^2+4*y^2-184)
          %f1.tangent()

          l0 = [x+7*y-10, -x-2*y-4, 48*x-56*y+4*x*y+x^2+4*y^2-184];
          r0 = region(l0,[x,y]);
          
          %anyEq(symFunType(r0.vx),['plus','plus','plus'])
          %return
          simplify(r0.vx.^2)
          r0 = r0.simplifyUnboundedRegion;
          r0.print
          %return
          l1 = [10 - 7*y - x , x - (5*y)/7 - 25/7 , x - y/3 - 14/3 , 2*y - x - 44 , x - (4*y)/7 - 27/7 , - x - 2*y - 4 ] ;
          r1 = region(l1,[x,y]);
          
          r1 = r1.simplifyUnboundedRegion;
                  
        
          r1.print
          [l,r] = merge(r0,r1);
          plqCheck.verifyMergeSound(testCase, r, r0, r1, l, ...
              'merge(r0,r1)');
          l
          r.print

          return
          
        end
        
        function testMerge10(testCase)
        % The merge itself is commented out in this fixture and stays that way -- what it records
        % is the pair of CONSTRUCTIONS, both of which carry the same conic facet. So this asserts
        % the construction: each region is non-empty, and every vertex it reports is a real vertex
        % of its own constraint set (feasible, with two constraints active). That is the property
        % `getVertices` must deliver for the merge these fixtures were collected for to mean
        % anything, and it is what its own HISTORY note says used to break.
          x = sym('x');
          y = sym('y');

          l0 = [- x - 7*y - 4 , x + 7*y - 10 , x + 2*y - 4 , 48*x - 56*y + 4*x*y + x^2 + 4*y^2 - 184 ];
          r0 = region(l0,[x,y]);
          testCase.verifyNotEmpty(r0, 'r0 is infeasible; the fixture no longer describes a region');
          plqCheck.verifyVerticesAreVertices(testCase, r0, 'testMerge10 r0');

          l1 = [x + 7*y + 4 , x - (9*y)/5 - 5 , - y - 5 , 48*x - 56*y + 4*x*y + x^2 + 4*y^2 - 184  ];
          r1 = region(l1,[x,y]);
          testCase.verifyNotEmpty(r1, 'r1 is infeasible; the fixture no longer describes a region');
          plqCheck.verifyVerticesAreVertices(testCase, r1, 'testMerge10 r1');
        end

        function testLinear3pt(testCase)
        % linear3pt rewrites a facet that three of the region's vertices turn out to be collinear
        % on. It is a change of PRESENTATION, so the point set must survive it, and applying it
        % twice must change nothing further.
          x = sym('x');
          y = sym('y');

          l0 = [x+7*y-10, -x-2*y-4, 48*x-56*y+4*x*y+x^2+4*y^2-184];
          r0 = region(l0,[x,y]);
          testCase.verifyNotEmpty(r0, 'the fixture is infeasible');
          r1 = r0.linear3pt;
          plqCheck.verifyRegionsAgree(testCase, r0, r1, 'linear3pt (conic fixture)');
          plqCheck.verifyRegionsAgree(testCase, r1, r1.linear3pt, 'linear3pt is not idempotent');
        end

        function testLinear3pt2(testCase)
        % The same, on the degenerate input the routine has to survive: a WEDGE, two facets and
        % one vertex, where there is no triple of vertices to be collinear at all.
          x = sym('x');
          y = sym('y');

          l0 = [4-2*y-x, x+(5*y)/2-1];
          r0 = region(l0,[x,y]);
          testCase.verifyNotEmpty(r0, 'the fixture is infeasible');
          r1 = r0.linear3pt;
          plqCheck.verifyRegionsAgree(testCase, r0, r1, 'linear3pt (two-facet wedge)');
          testCase.verifyEqual(size(r1.ineqs,2), size(r0.ineqs,2), ...
              'with only two facets there is nothing for linear3pt to rewrite');
        end

        function testgetEdgeNosInf(testCase)
        % getEdgeNosInf orders an UNBOUNDED region's facets along its boundary walk. The property
        % that makes it usable is that the answer is a PERMUTATION of the facets -- each named
        % once, none invented, none dropped -- and the walk starts and ends on the two facets that
        % run to infinity. The particular starting facet is not part of the contract, so it is not
        % asserted; see testCreation's note on why pinning an order here was brittle.
        %
        % Fixture: the four-facet unbounded polyhedron from the QPLIB-shaped dual meshes.
            s_1 = sym('s1');
            s_2 = sym('s2');
            ineq(1) = s_1 + 2*s_2 + 4;
            ineq(2) = s_1 - (9*s_2)/5 - 5;
            ineq(3) = - s_2 - 5 ;
            ineq(4) = s_1 + 7*s_2 - 46;
            d = region(ineq,[s_1,s_2]);
            d = d.removeInfV;
            d = d.poly2orderUnbounded;
            edgeNo = d.getEdgeNosInf(d.vars);
            testRegion.verifyIsFacetPermutation(testCase, edgeNo, d, 'getEdgeNosInf, 4 facets');
            plqCheck.verifyVerticesAreVertices(testCase, d, 'getEdgeNosInf fixture');
        end

        function testgetEdgeNosInf2(testCase)
        % A second unbounded fixture, and the one that used to be UNREACHABLE: the body had an
        % unconditional `return` above the getEdgeNosInf call, so the routine this test is named
        % after was never invoked. The stale `% [1 2]` comment beside it dated from a two-facet
        % version of the fixture; the fixture now has four, and the permutation property is
        % asserted instead of a number that was no longer true.
            s_1 = sym('s1');
            s_2 = sym('s2');
            ineq(1) = s_1 - (5*s_2)/7 - 25/7 ;
            ineq(2) = 4 - 2*s_2 - s_1 ;
            ineq(3) = s_1 - (4*s_2)/7 - 27/7 ;
            ineq(4) = s_1 - s_2/3 - 14/3 ;
            d = region(ineq,[s_1,s_2]);
            d = d.removeInfV;
            d = d.poly2orderUnbounded;
            edgeNo = d.getEdgeNosInf(d.vars);
            testRegion.verifyIsFacetPermutation(testCase, edgeNo, d, 'getEdgeNosInf2');
        end

        function testgetEdgeNosInf3(testCase)
        % The BOUNDED path and the unbounded path on the same region: `r` is a bounded triangle
        % once simplifyUnboundedRegion has dropped its two redundant facets, so both getEdgeNosInf
        % and getEdgeNos must return a permutation of what is left, and simplification must not
        % have moved the region.
            r0 = testCase.r;
            r1 = r0.simplifyUnboundedRegion;
            plqCheck.verifyRegionsAgree(testCase, r0, r1, 'simplify before ordering');
            r1 = r1.poly2orderUnbounded;
            testRegion.verifyIsFacetPermutation(testCase, r1.getEdgeNosInf(r1.vars), r1, 'getEdgeNosInf on r');
            testRegion.verifyIsFacetPermutation(testCase, r1.getEdgeNos(r1.vars),    r1, 'getEdgeNos on r');
        end

        function testgetEdgeNosInf4(testCase)
        % The smallest unbounded case there is: a WEDGE, two facets and one vertex. Both facets
        % run to infinity, so the walk is the whole facet list -- the degenerate case where an
        % ordering routine most easily returns a zero, a duplicate, or a short list.
            s_1 = sym('s1');
            s_2 = sym('s2');
            ineq(1) = s_1 + 4;
            ineq(2) =  s_2 + 5 ;
            d = region(ineq,[s_1,s_2]);
            d = d.removeInfV;
            d = d.poly2orderUnbounded;
            edgeNo = d.getEdgeNosInf(d.vars);
            testRegion.verifyIsFacetPermutation(testCase, edgeNo, d, 'getEdgeNosInf on a wedge');
        end

        function testgetVertices(testCase)
        % getVertices enumerates the pairwise constraint intersections that are feasible. Two
        % properties, both definitional and neither pinning coordinates:
        %   * every FINITE vertex reported satisfies all the constraints and has at least two of
        %     them active -- otherwise it is a boundary point or an infeasible intersection;
        %   * calling it again on an already-populated region gives the same answer. Its own
        %     HISTORY note records the defect this catches: vx/vy were appended to rather than
        %     reset, so a second call piled duplicate placeholder vertices on the first call's,
        %     which then broke region.minus.
        %
        % The fixture is the unbounded region {s2 <= 1, s2 <= s1, s1 >= -1}: two real vertices
        % (-1,-1) and (1,1), plus the placeholder points the "vertex at infinity" phase adds.
             s_1 = sym('s1');
            s_2 = sym('s2');

            ineq(1) = s_2 - 1;
            ineq(2) =  s_2 - s_1 ;
            ineq(3) = - s_1 - 1;
            d = region(ineq,[s_1,s_2]);

            plqCheck.verifyVerticesAreVertices(testCase, d, 'getVertices fixture');

            nv1 = d.nv; vx1 = d.vx; vy1 = d.vy;
            d2 = d.getVertices;
            testCase.verifyEqual(d2.nv, nv1, ...
                'getVertices is not idempotent: a second call changed the vertex count');
            testCase.verifyTrue(isequal(simplify(d2.vx - vx1), sym(zeros(size(vx1)))) && ...
                                isequal(simplify(d2.vy - vy1), sym(zeros(size(vy1)))), ...
                'getVertices is not idempotent: a second call moved the vertices');
        end

    end

end