classdef plqCheck
% plqCheck  Numeric verification for the cPLQ pipeline's stages, so a "crash test" becomes a
%   test.
%
% WHY THIS EXISTS. `testcPLQ` and `testMaxMultiRegion` between them held 32 tests and ZERO
% assertions: each ran `triangulate -> convexEnvelope -> conjugate -> maximum -> biconjugateF`,
% printed the answer, and returned. They passed if nothing threw. That is worth something -- most
% of this pipeline's defects this month WERE crashes -- but it is not verification, and it cost
% 90% of a two-hour bucket to learn that a function returns. Measured 2026-08-19:
% `testcPLQ/testRectBiconj` alone was 3198 s.
%
% WHAT IS CHECKED, and each is a DEFINITION rather than a golden value, so nothing here needs
% re-pinning when a representation changes:
%
%   convex envelope   co f <= f on the domain, and co f = f at every vertex of it (an envelope
%                     that is not an underestimator is wrong; one that dips below the vertex
%                     values is wrong; and on these fixtures the vertices are where the envelope
%                     touches).
%   conjugate         f*(s) = sup_{x in D} <s,x> - f(x), against a numeric sup over the domain.
%                     This is the strongest of the three and needs no convexity assumption.
%   biconjugate       f** <= f on the domain, f** convex along segments, and f** = f wherever f
%                     is already convex there.
%
% NUMERIC SUP, and why sampling is sound here. The domains are convex polygons and f is a
% quadratic, so `<s,x> - f(x)` is a quadratic whose sup over the polygon is attained on the
% boundary or at an interior critical point. Both are sampled: the whole boundary at a fine step,
% every vertex exactly, plus an interior grid. That gives a LOWER bound on the true sup, so
% `f*(s) < sup_sampled - tol` is a definite failure, while `f*(s) > sup_sampled` within tolerance
% is expected (the sampler cannot exceed the truth). Both directions are checked, with the
% one-sided one given the tolerance it deserves.
%
% TOLERANCES are absolute-and-relative and deliberately loose (1e-6 relative): the point is to
% catch a wrong branch, a dropped cell or a mis-oriented constraint, all of which are O(1) errors
% -- not to police the last ulp, which the exact-arithmetic tests elsewhere already do.

    methods (Static)

        function pts = domainSample(d, nB, nI)
        % Points of the polygon `d`: every vertex, nB points along each edge, and an interior
        % grid of about nI points. `d` is a `domain` (its polygon carries vx/vy).
            if nargin < 2, nB = 40; end
            if nargin < 3, nI = 200; end
            V = [double(d.polygon.vx(:)), double(d.polygon.vy(:))];
            n = size(V,1);
            pts = V;
            for k = 1:n
                a = V(k,:); b = V(mod(k,n)+1,:);
                t = linspace(0, 1, nB+2)'; t = t(2:end-1);
                pts = [pts; a + t.*(b - a)];                          %#ok<AGROW>
            end
            % interior: rejection-sample the bounding box against the polygon
            lo = min(V,[],1); hi = max(V,[],1);
            m = ceil(sqrt(nI)) * 4;
            [gx, gy] = meshgrid(linspace(lo(1), hi(1), m), linspace(lo(2), hi(2), m));
            C = [gx(:), gy(:)];
            inside = plqCheck.inPolygon(C, V);
            C = C(inside,:);
            if size(C,1) > nI
                C = C(round(linspace(1, size(C,1), nI)), :);
            end
            pts = [pts; C];
        end

        function tf = inPolygon(P, V)
        % Convex-polygon membership by half-planes, orientation-agnostic. `inpolygon` needs the
        % Mapping/Image toolboxes on some installs; this needs nothing and these domains are
        % convex by construction.
            n = size(V,1);
            s = zeros(size(P,1), n);
            for k = 1:n
                a = V(k,:); b = V(mod(k,n)+1,:);
                s(:,k) = (b(1)-a(1))*(P(:,2)-a(2)) - (b(2)-a(2))*(P(:,1)-a(1));
            end
            tol = 1e-9 * max(1, max(abs(s), [], 'all'));
            tf = all(s >= -tol, 2) | all(s <= tol, 2);
        end

        % ==========================================================================================
        % REGION-LEVEL DEFINITION CHECKS
        %
        % WHY THESE ARE HERE. `testRegion` and `testfunctionNDomain` between them held 22 tests
        % that printed a region and returned -- the same shape this class already replaced for the
        % cPLQ pipeline. A region is a POINT SET, {p : g_k(p) <= 0 for every k}, and every
        % operation on one (merge, simplifyUnboundedRegion, removeTangent, linear3pt) has a
        % contract stated in terms of that point set. So the check is the definition: SAMPLE the
        % set and compare memberships. Nothing here pins a constraint list or a vertex order,
        % which is what made the previous assertions in this suite brittle.
        %
        % THE MARGIN, and why the sample is one-sided. A sampled point sitting exactly on a facet
        % belongs to both sides and decides nothing, so `regionSample` returns only points that
        % are STRICTLY inside by `margin`, and a containment test then accepts the other region
        % up to `+margin`. A defect that moves a facet by less than the margin is invisible to
        % this check -- but every defect these operations can have (a dropped facet, a facet on
        % the wrong side, a merge that swallows a non-convex gap) is O(1), not O(margin).

        function tf = inRegion(r, P, margin)
        % Membership in a `region` by its own definition: every constraint <= margin at P.
        % [input]  r      : a scalar region
        %          P      : k x 2 points
        %          margin : slack added to the right-hand side (default 0)
        % [output] tf     : k x 1 logical
        %
        % MEMOIZED HANDLES, and it is not an optimisation. `matlabFunction` calls the symbolic
        % engine, so it costs tens of milliseconds per constraint per call -- and the callers below
        % ask the same question of the same region dozens of times (both directions of a
        % containment, then again after the operation, then again for idempotence). Measured
        % 2026-08-31: without this cache the assertions added to `testRegion` took that suite from
        % ~120 s to ~13 min, which is a fast-bucket budget spent on re-compiling the same three
        % inequalities. The key is the constraint's own `char`, which is exactly what determines
        % the handle, so the cache cannot answer a different question than the caller asked.
            if nargin < 3, margin = 0; end
            if isempty(r), tf = false(size(P,1),1); return, end
            tf = true(size(P,1), 1);
            v = r.vars;
            for k = 1:size(r.ineqs, 2)
                h = plqCheck.handleFor(r.ineqs(k).f, v);
                g = arrayfun(@(i) double(h(P(i,1), P(i,2))), (1:size(P,1))');
                tf = tf & (g <= margin);
            end
        end

        function h = handleFor(f, v)
        % A numeric function handle for the symbolic expression f in the variables v, cached on
        % (char(f), char(v)). See inRegion for why.
            persistent cache
            if isempty(cache), cache = containers.Map('KeyType','char','ValueType','any'); end
            key = [char(f) '##' char(v(1)) ',' char(v(2))];
            if isKey(cache, key), h = cache(key); return, end
            h = matlabFunction(f, 'Vars', {v(1), v(2)});
            cache(key) = h;
        end

        % ---- FIXTURE BUILDERS -------------------------------------------------------------------
        % Every unit suite needs the same four shapes -- a bounded triangle, a box, a wedge (one
        % vertex, two rays) and a half-strip (two vertices, two parallel rays) -- because those are
        % the only shapes `plq_1p` and `fanUnboundedFace` emit and therefore the only ones the
        % routines downstream are specified for. Each suite used to rebuild them by hand, which is
        % why no two fixtures agreed. Built from EXPLICIT half-planes so the sign convention is
        % written once: `region` reads every constraint as g(x,y) <= 0.

        function r = halfPlanes(A, b, vars)
        % {x : A x <= b} as a `region`, from numeric rows. The rows are kept exactly as given --
        % no normalisation -- so a caller can hand in the badly scaled rows a scale test needs.
            g = sym.empty(1,0);
            for i = 1:size(A,1)
                g(i) = A(i,1)*vars(1) + A(i,2)*vars(2) - b(i);
            end
            r = region(g, vars);
        end

        function r = triRegion(V, vars)
        % The closed triangle on the three rows of V, given COUNTER-CLOCKWISE. Each edge a->b
        % contributes "the interior is on the left": cross(b-a, p-a) >= 0.
            if size(V,1) ~= 3, error('plqCheck:triRegion', 'V must be 3 x 2.'); end
            A = zeros(3,2); b = zeros(3,1);
            for i = 1:3
                a = V(i,:); c = V(mod(i,3)+1,:); e = c - a;
                A(i,:) = [e(2), -e(1)];          % -cross(e, p-a) <= 0
                b(i)   = e(2)*a(1) - e(1)*a(2);
            end
            r = plqCheck.halfPlanes(A, b, vars);
        end

        function r = boxRegion(lo, hi, vars)
        % The axis-aligned box [lo(1),hi(1)] x [lo(2),hi(2)] -- the SCIP/QPLIB domain shape.
            A = [-1 0; 1 0; 0 -1; 0 1];
            b = [-lo(1); hi(1); -lo(2); hi(2)];
            r = plqCheck.halfPlanes(A, b, vars);
        end

        function r = wedgeRegion(v, d1, d2, vars)
        % The pointed cone v + cone(d1,d2), with d1 -> d2 turning COUNTER-CLOCKWISE. Refuses a
        % non-pointed pair rather than returning a half-plane the caller did not ask for.
            c = d1(1)*d2(2) - d1(2)*d2(1);
            if c <= 0
                error('plqCheck:wedgeRegion', ...
                    ['d1 -> d2 must turn counter-clockwise through less than pi (cross = %g). ' ...
                     'A flat or reflex pair is a half-plane or a non-convex set, not a wedge.'], c);
            end
            A = [ d1(2), -d1(1);          % -cross(d1, p-v) <= 0
                 -d2(2),  d2(1)];         % -cross(p-v, d2) <= 0
            b = [ d1(2)*v(1) - d1(1)*v(2);
                 -d2(2)*v(1) + d2(1)*v(2)];
            r = plqCheck.halfPlanes(A, b, vars);
        end

        function r = halfStripRegion(v1, v2, d, vars)
        % conv{v1,v2} + cone(d): the two-vertex, two-parallel-ray shape. Three facets -- the two
        % rays and the segment -- so it is the smallest fixture with both a bounded and an
        % unbounded facet, which is where clipping routines usually part company with their spec.
        % Written in the (s,t) COORDINATES the set is defined by rather than by picking signs for
        % three normals: p = v1 + s e + t d with e = v2 - v1, so the set is exactly
        % {s >= 0, s <= 1, t >= 0}, and (s,t) is an affine function of p through inv([e d]).
        % Choosing normals by hand got two of the three facets backwards and produced an empty
        % region, which then failed several stages downstream instead of here.
            e = (v2(:) - v1(:));
            C = [e, d(:)];
            if abs(det(C)) <= eps * max(1, norm(C))
                error('plqCheck:halfStripRegion', 'd must not be parallel to v2 - v1.');
            end
            W = inv(C);                                   %#ok<MINV>  2x2, exactness not at stake
            A = [-W(1,:); W(1,:); -W(2,:)];
            b = [-W(1,:)*v1(:); 1 + W(1,:)*v1(:); -W(2,:)*v1(:)];
            r = plqCheck.halfPlanes(A, b, vars);
        end

        function b = regionBox(rs, pad)
        % A sampling window big enough to see the interesting part of a family of regions: the
        % bounding box of every FINITE vertex any of them has, padded. Unbounded regions run off
        % this window, which is fine -- the check is over the window, not over the plane.
        %
        % The pad is RELATIVE to the vertices' own spread, not a constant. A fixed pad of 8 around
        % the unit triangle left a 45 x 45 grid with three interior points -- technically not
        % vacuous, and useless.
            V = zeros(0,2);
            for i = 1:numel(rs)
                r = rs(i);
                if isempty(r), continue, end
                for j = 1:r.nv
                    xv = double(r.vx(j)); yv = double(r.vy(j));
                    if isfinite(xv) && isfinite(yv) && abs(xv) < 1e6 && abs(yv) < 1e6
                        V(end+1,:) = [xv yv];                                  %#ok<AGROW>
                    end
                end
            end
            if isempty(V), b = [-10 10 -10 10]; return, end
            if nargin < 2 || isempty(pad)
                pad = max(1, 0.6 * max(max(V,[],1) - min(V,[],1)));
            end
            lo = min(V,[],1) - pad; hi = max(V,[],1) + pad;
            b = [lo(1) hi(1) lo(2) hi(2)];
        end

        function P = regionSample(r, box, n, margin)
        % Points STRICTLY inside r (by `margin`) on an n x n grid over `box = [x0 x1 y0 y1]`.
        % Returns 0 x 2 when the region misses the window entirely -- callers must treat that as
        % "nothing checked" rather than "check passed", which `verifyRegionSubset` does.
        %
        % ONE ADAPTIVE RETRY, and it is not a nicety. `regionBox` sizes the window from the finite
        % vertices, so a region with a single vertex gets a window of pad-size around it -- and a
        % PARABOLIC region with one vertex extends far beyond that. Measured on testRegion's
        % removeTangent fixture: the default window found 0 interior points and the check went
        % silently vacuous, while [-60,60]^2 found 3712 and the check was decisive. Widening
        % unconditionally instead would cost density on the small fixtures, so widen only when
        % there is nothing to look at.
            if nargin < 3, n = 45; end
            if nargin < 4, margin = 1e-6; end
            P = gridInside(box);
            if isempty(P)
                cx = (box(1)+box(2))/2; cy = (box(3)+box(4))/2;
                hx = max(box(2)-box(1), 1) * 10; hy = max(box(4)-box(3), 1) * 10;
                P = gridInside([cx-hx, cx+hx, cy-hy, cy+hy]);
            end

            function Q = gridInside(b)
                [gx, gy] = meshgrid(linspace(b(1), b(2), n), linspace(b(3), b(4), n));
                C = [gx(:), gy(:)];
                Q = C(plqCheck.inRegion(r, C, -margin), :);
            end
        end

        function nChecked = verifyRegionSubset(tc, rSub, rSup, box, name, margin)
        % Every point strictly inside rSub must lie in rSup. Returns how many points were
        % actually tested so the caller can assert the check was not vacuous.
            if nargin < 6, margin = 1e-6; end
            P = plqCheck.regionSample(rSub, box, 45, margin);
            nChecked = size(P,1);
            if nChecked == 0, return, end
            in = plqCheck.inRegion(rSup, P, margin);
            bad = find(~in, 1);
            if isempty(bad), bad = 1; end
            tc.verifyTrue(all(in), sprintf( ...
                '%s: %d of %d sampled points of the subset are outside the superset (first at (%g,%g))', ...
                name, sum(~in), nChecked, P(bad,1), P(bad,2)));
        end

        function nChecked = verifyRegionsAgree(tc, rA, rB, name, box)
        % rA and rB describe the SAME point set over the window. Used for every operation that is
        % documented to rewrite a region's presentation without moving its boundary --
        % simplifyUnboundedRegion, removeTangent, linear3pt.
        %
        % RETURNS how many points were actually compared, and does NOT itself assert that the
        % number is positive. A region with empty interior -- and testRegion's own fixtures
        % include an infeasible one and one that is a single point -- has nothing to sample, and
        % that is a property of the fixture, not a failure. The caller knows which of its inputs
        % are supposed to be two-dimensional and asserts non-vacuity there; `finiteVertexSet`
        % below is the check for the degenerate ones.
            if nargin < 5, box = plqCheck.regionBox([rA, rB]); end
            nA = plqCheck.verifyRegionSubset(tc, rA, rB, box, [name ' (A subset B)']);
            nB = plqCheck.verifyRegionSubset(tc, rB, rA, box, [name ' (B subset A)']);
            nChecked = nA + nB;
        end

        function V = finiteVertexSet(r)
        % The region's finite vertices as a sorted, de-duplicated numeric list -- the only thing
        % left to compare when a region has no interior to sample. Sorted, because the vertex
        % ORDER is a Symbolic Math Toolbox detail (see testRegion/testCreation).
            V = zeros(0,2);
            if isempty(r), return, end
            for j = 1:r.nv
                xv = double(r.vx(j)); yv = double(r.vy(j));
                if isfinite(xv) && isfinite(yv) && abs(xv) < 1e6 && abs(yv) < 1e6
                    V(end+1,:) = [xv yv];                                      %#ok<AGROW>
                end
            end
            V = unique(round(V, 9), 'rows');
        end

        function verifyMergeSound(tc, rM, rA, rB, isExact, name, box)
        % `region.merge`'s own contract, and it has TWO branches, not one.
        %
        % The signature is `[l, obj] = merge(obj, obj2)`: the result overwrites the FIRST
        % argument. So when merge declines -- no shared facet, more than one shared facet, or
        % `unionIsExact` withholding its certificate -- it returns `l = false` and the first
        % operand UNCHANGED. It has not produced a union and must not be asserted as one.
        %
        %   l TRUE   rM is exactly A u B. Both halves are asserted: rM contains A and contains B
        %            (merge deletes the shared facet from each and intersects the rest, which
        %            `unionIsExactCompute`'s header proves can never lose a point), and rM adds
        %            nothing (every point of rM lies in A or in B). The second half is the one
        %            that catches the over-claiming defect merge's own header describes -- three
        %            Step 3 cells merging into one that covered a point none of them did.
        %   l FALSE  no merge happened, so the contract is that A came back untouched. Asserted
        %            as a point-set identity rather than as a constraint-list identity, because
        %            merge may return the pre-restore copy `obj3` on the empty-union path.
            if nargin < 7, box = plqCheck.regionBox([rA, rB, rM]); end
            if isempty(rM)
                tc.verifyTrue(isempty(rA) && isempty(rB), sprintf( ...
                    '%s: merge returned empty for non-empty inputs', name));
                return
            end
            if ~isExact
                plqCheck.verifyRegionsAgree(tc, rA, rM, ...
                    [name ' declined, so the first operand must come back unchanged'], box);
                tc.verifyEqual(size(rM.ineqs,2), size(rA.ineqs,2), sprintf( ...
                    ['%s: merge declined but returned %d constraints against the first ' ...
                     'operand''s %d -- a declined merge must change nothing'], ...
                    name, size(rM.ineqs,2), size(rA.ineqs,2)));
                return
            end
            nA = plqCheck.verifyRegionSubset(tc, rA, rM, box, [name ' (A subset merge)']);
            nB = plqCheck.verifyRegionSubset(tc, rB, rM, box, [name ' (B subset merge)']);
            tc.verifyGreaterThan(nA + nB, 0, sprintf( ...
                '%s: neither input has a sampled interior point -- the check was vacuous', name));
            P = plqCheck.regionSample(rM, box, 45, 1e-6);
            if isempty(P), return, end
            in = plqCheck.inRegion(rA, P, 1e-6) | plqCheck.inRegion(rB, P, 1e-6);
            bad = find(~in, 1);
            if isempty(bad), bad = 1; end
            tc.verifyTrue(all(in), sprintf( ...
                ['%s: merge was reported EXACT but %d of %d of its own points lie in neither ' ...
                 'input (first at (%g,%g)) -- the union is not convex there'], ...
                name, sum(~in), numel(in), P(bad,1), P(bad,2)));
        end

        function verifyVerticesAreVertices(tc, r, name)
        % Every finite vertex a region reports must (a) satisfy all of its constraints and
        % (b) make at least TWO of them active -- that is what makes it a vertex rather than a
        % boundary point. Catches the duplicate-vertex and stale-vertex failures getVertices'
        % own HISTORY note describes, without pinning coordinates.
            nFinite = 0;
            for j = 1:r.nv
                xv = double(r.vx(j)); yv = double(r.vy(j));
                if ~isfinite(xv) || ~isfinite(yv) || abs(xv) >= 1e6 || abs(yv) >= 1e6, continue, end
                nFinite = nFinite + 1;
                g = zeros(1, size(r.ineqs,2));
                for k = 1:size(r.ineqs,2)
                    g(k) = double(subs(r.ineqs(k).f, r.vars, [xv yv]));
                end
                sc = max(1, max(abs(g)));
                tc.verifyLessThanOrEqual(max(g), 1e-7*sc, sprintf( ...
                    '%s: reported vertex %d = (%g,%g) violates its own constraints by %.3g', ...
                    name, j, xv, yv, max(g)));
                tc.verifyGreaterThanOrEqual(sum(abs(g) <= 1e-7*sc), 2, sprintf( ...
                    '%s: reported vertex %d = (%g,%g) has %d active constraints; a vertex needs 2', ...
                    name, j, xv, yv, sum(abs(g) <= 1e-7*sc)));
            end
            tc.verifyGreaterThan(nFinite, 0, sprintf( ...
                '%s: no finite vertex was reported, so nothing was checked', name));
        end

        function [v, ok] = safeEval(fnd, pt)
        % evalFunctionNDomain at one point, with a RATIONAL face's pole treated as "no value
        % here" rather than an error. Step 1's envelope is quadratic-over-LINEAR, and its
        % denominator vanishes exactly at some domain vertices -- which is where an underestimator
        % check most wants to sample. MATLAB raises symbolic:kernel:DivisionByZero there, so the
        % sampler must skip the point instead of the test dying on the fixture's own geometry.
            v = NaN; ok = false;
            try
                v = evalFunctionNDomain(fnd, pt);
                ok = isfinite(v);
            catch
                v = NaN; ok = false;
            end
        end

        function v = evalSym(f, vars, X)
        % f evaluated at the rows of X, numerically.
            h = matlabFunction(f, 'Vars', {vars(1), vars(2)});
            v = arrayfun(@(i) double(h(X(i,1), X(i,2))), (1:size(X,1))');
        end

        function s = supOverDomain(fExpr, vars, d, sPt, nB, nI)
        % sup_{x in D} <sPt,x> - f(x). A LOWER bound on the true sup; see the header.
        %
        % THE GRID ALONE IS NOT ACCURATE ENOUGH FOR THE TOLERANCE THE CALLERS ASSERT, and that was
        % costing a false red. `conjugateMatchesSup` checks BOTH directions at 1e-5 relative; the
        % upper one asks the reference to be accurate to 1e-5, and a 200-per-edge / 900-interior
        % grid is not. MEASURED on `testcPLQ`'s PRect piece 3 -- the triangle {(0,-4),(1,3),
        % (1.8708,-0.2583)} carrying x*y -- at s = (0,0): the true sup is EXACTLY 2, `plq_1p`
        % returns 2, and the grid returns 1.9999666056. Short by 3.34e-05, asserted at 1e-05, and
        % the test had been reading that as the conjugate exceeding the sup.
        %
        % So refine the grid with the CLOSED-FORM candidates. For a quadratic f the objective
        % `<s,x> - f(x)` is quadratic, and its sup over a convex polygon is attained at a vertex,
        % at a per-edge stationary point, or at the interior stationary point -- a finite list,
        % computed exactly. Taking the max of the grid and that list is SAFE in the strong sense:
        % both are lower bounds, so the result is still a lower bound and can only get closer to
        % the truth. For a quadratic on a bounded polygon it IS the truth.
        %
        % One consequence to expect rather than be surprised by: a better lower bound makes the
        % `f* < sup - tol` direction fire more readily. That direction is the definite-defect one,
        % so anything it newly catches is a real minorant the coarse grid was hiding.
            if nargin < 5, nB = 200; end
            if nargin < 6, nI = 900; end
            X = plqCheck.domainSample(d, nB, nI);
            fv = plqCheck.evalSym(fExpr, vars, X);
            s = max(X * sPt(:) - fv);
            s = max(s, plqCheck.supQuadExact(fExpr, vars, d, sPt));
        end

        function v = supQuadExact(fExpr, vars, d, sPt)
        % The closed-form candidate maximisers of `<s,x> - f(x)` over the polygon of `d`, for a
        % QUADRATIC f: the vertices, the stationary point along each edge, and the interior
        % stationary point. Returns -inf when f is not quadratic or the polygon is unusable, so the
        % caller's grid stands alone -- this only ever RAISES a lower bound.
            v = -inf;
            try
                H = double(hessian(fExpr, vars));
                g0 = double(subs(gradient(fExpr, vars), vars, [0 0]));
                c0 = double(subs(fExpr, vars, [0 0]));
                if any(~isfinite(H(:))) || any(~isfinite(g0)) || ~isfinite(c0), return, end
                W = [double(d.polygon.vx(:)), double(d.polygon.vy(:))];
            catch
                return                      % not a quadratic, or no polygon: leave it to the grid
            end
            W = W(all(isfinite(W), 2) & all(abs(W) < 1e12, 2), :);   % drop ray markers
            n = size(W,1);
            if n < 2, return, end
            obj = @(p) sPt(:).' * p(:) - (0.5*p(:).'*H*p(:) + g0(:).'*p(:) + c0);
            cand = W;
            % interior stationary point of the objective: H x = s - g0
            if abs(det(H)) > 1e-12 * max(1, norm(H))
                cand = [cand; (H \ (sPt(:) - g0(:))).'];
            end
            for e = 1:n
                a = W(e,:); b = W(mod(e,n)+1,:); dir = b - a;
                den = dir * H * dir.';
                if abs(den) > 1e-14
                    t = (sPt(:).'*dir(:) - (H*a(:) + g0(:)).'*dir(:)) / den;
                    if t > 0 && t < 1, cand = [cand; a + t*dir]; end %#ok<AGROW>
                end
            end
            for k = 1:size(cand,1)
                p = cand(k,:);
                if ~plqCheck.inConvexPolygon(p, W), continue, end
                v = max(v, obj(p));
            end
        end

        function tf = inConvexPolygon(p, W)
        % p inside the CCW-or-CW convex polygon W, with a slack tolerance. Used only to reject
        % candidate maximisers that fall outside, so a false accept would merely make the lower
        % bound wrong -- hence the tolerance is tight and the orientation is handled by sign vote.
            n = size(W,1); s = zeros(n,1);
            for i = 1:n
                a = W(i,:); b = W(mod(i,n)+1,:);
                s(i) = (b(1)-a(1))*(p(2)-a(2)) - (b(2)-a(2))*(p(1)-a(1));
            end
            sc = max(1, max(abs(s)));
            tf = all(s >= -1e-9*sc) || all(s <= 1e-9*sc);
        end

        % ---- the three checks -----------------------------------------------------------

        function envelopeUnderestimates(tc, piece, name)
        % co f <= f on the domain, and co f = f at every vertex. `piece` is a plq_1p AFTER
        % convexEnvelope, so piece.envelope is the functionNDomain array of the envelope.
            d = piece.d; vars = d.polygon.vars;
            X = plqCheck.domainSample(d, 12, 60);
            fv = plqCheck.evalSym(piece.f.f, vars, X);
            ev = arrayfun(@(i) plqCheck.safeEval(piece.envelope, X(i,:)), (1:size(X,1))');
            got = ~isnan(ev) & ~isnan(fv);
            tc.verifyGreaterThan(nnz(got), 0, ...
                sprintf('%s: the envelope covers no sampled point of its own domain', name));
            tol = 1e-6 * max(1, max(abs(fv)));
            tc.verifyLessThanOrEqual(max(ev(got) - fv(got)), tol, ...
                sprintf('%s: the convex envelope EXCEEDS f somewhere on the domain', name));

            V = [double(d.polygon.vx(:)), double(d.polygon.vy(:))];
            fvV = plqCheck.evalSym(piece.f.f, vars, V);
            evV = arrayfun(@(i) plqCheck.safeEval(piece.envelope, V(i,:)), (1:size(V,1))');
            ok = ~isnan(evV) & ~isnan(fvV);
            tc.verifyLessThanOrEqual(max(abs(evV(ok) - fvV(ok))), 1e-6 * max(1, max(abs(fvV))), ...
                sprintf('%s: the envelope does not touch f at the domain vertices', name));
        end

        function conjugateMatchesSup(tc, fnd, fExpr, vars, d, S, name)
        % f*(s) = sup_{x in D} <s,x> - f(x) at each dual point in the rows of S.
            for i = 1:size(S,1)
                s = S(i,:);
                got = plqCheck.safeEval(fnd, s);
                want = plqCheck.supOverDomain(fExpr, vars, d, s);
                tc.verifyFalse(isnan(got), sprintf( ...
                    '%s: f* is uncovered at (%g,%g) -- a cell is missing', name, s(1), s(2)));
                if isnan(got), continue, end
                tol = 1e-5 * max(1, abs(want));
                % The sampler is a LOWER bound, so got < want is a definite defect; got > want
                % beyond tolerance means f* claims more than the definition allows.
                tc.verifyGreaterThanOrEqual(got, want - tol, sprintf( ...
                    '%s: f*(%g,%g) = %.9g is BELOW the sampled sup %.9g', ...
                    name, s(1), s(2), got, want));
                tc.verifyLessThanOrEqual(got, want + tol, sprintf( ...
                    '%s: f*(%g,%g) = %.9g EXCEEDS the sup %.9g over the domain', ...
                    name, s(1), s(2), got, want));
            end
        end

        function biconjugateIsAConvexUnderestimator(tc, bnd, fExpr, vars, d, name)
        % f** <= f on the domain, and f** convex along sampled segments.
            X = plqCheck.domainSample(d, 10, 50);
            fv = plqCheck.evalSym(fExpr, vars, X);
            bv = arrayfun(@(i) plqCheck.safeEval(bnd, X(i,:)), (1:size(X,1))');
            ok = ~isnan(bv) & ~isnan(fv);
            tc.verifyGreaterThan(nnz(ok), 0, ...
                sprintf('%s: f** covers no sampled point of the domain', name));
            tol = 1e-6 * max(1, max(abs(fv)));
            tc.verifyLessThanOrEqual(max(bv(ok) - fv(ok)), tol, ...
                sprintf('%s: f** EXCEEDS f somewhere on the domain', name));

            % convexity along segments between sampled points
            rng(20260819);
            idx = randi(size(X,1), 60, 2);
            worst = -inf;
            for k = 1:size(idx,1)
                a = X(idx(k,1),:); b = X(idx(k,2),:);
                m = (a + b)/2;
                va = plqCheck.safeEval(bnd, a);
                vb = plqCheck.safeEval(bnd, b);
                vm = plqCheck.safeEval(bnd, m);
                if any(isnan([va vb vm])), continue, end
                worst = max(worst, vm - (va + vb)/2);
            end
            if isfinite(worst)
                tc.verifyLessThanOrEqual(worst, 1e-6 * max(1, max(abs(fv))), sprintf( ...
                    '%s: f** is not convex -- a midpoint sits above the chord by %.3g', ...
                    name, worst));
            end
        end
    end
end
