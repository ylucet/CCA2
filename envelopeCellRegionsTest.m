classdef envelopeCellRegionsTest < matlab.unittest.TestCase
% Unit tests for envelopeCellRegions, which CONSTRUCTS the primal subdivision of one cell's
% contribution to the convex envelope -- each candidate's region as exact integer half-planes in x,
% and the quadratic it carries there.
%
% BUCKET: fast (conjQ of one triangle, then integer arithmetic; no symbolic call).
%
% THE ORACLES, in the order of how much they are trusted.
%   1. NO OVERSHOOT, which needs no oracle at all. h is a SUP, so for any s feasible in the cell,
%      h(x) >= <s,x> - g(s). A face value that some feasible s beats is a definite defect, and this
%      is the check that found the one real bug in this routine (see the transpose test below).
%   2. THE COMPLEMENT, decided by DIRECTIONS rather than by sampling. Along a recession direction r
%      the objective behaves as t*<r,Y> - t^2 * r'Nr/2, so the sup is +infinity exactly when some r
%      in rec(P) has r'Nr = 0 and r'Y > 0. The regions must be EXACTLY the set where that fails.
%   3. An independent enumeration of the sup, written here in doubles. Weakest of the three,
%      because it shares the mathematics -- it is a cross-check on the algebra, not on the idea.

    methods (Test)

        function NOFeasiblePointMayBeatAFaceValue(testCase)
        % THE ORACLE-FREE CHECK, and the strongest statement available. h(x) is a supremum over the
        % cell, so any feasible s gives a lower bound. If a face value can be beaten by a feasible
        % s then that face UNDERSTATES the sup, which is a defect whatever the region says.
            rng(20260913);
            g = conjQ(envelopeCellRegionsTest.tri([0 0; 1 0; 1 1], 0));
            worst = -inf;  nTested = 0;
            for k = 1:g.nf
                [N, m, q, dd, P] = envelopeCellRegionsTest.cellData(g, k);
                if isempty(P), continue, end
                faces = envelopeCellRegions(N, m, q, dd, P);
                X = (rand(60,2) - 0.5) * 8;
                for t = 1:size(X,1)
                    x = X(t,:);
                    f = envelopeCellRegionsTest.locate(faces, x);
                    if isempty(f), continue, end
                    got = ratQ.evalFace(f.fN, f.fD, x);
                    S = envelopeCellRegionsTest.feasiblePts(P, 200);
                    if isempty(S), continue, end
                    lo = max(S*x.' - (0.5*sum((S*N).*S,2) + S*m + q)/dd);
                    worst = max(worst, lo - got);
                    nTested = nTested + 1;
                end
            end
            testCase.verifyGreaterThan(nTested, 100, 'the check was actually exercised');
            testCase.verifyLessThan(worst, 1e-9, ...
                sprintf('a feasible s beat its own face by %.3e', worst));
        end

        function theREGIONSAreEXACTLYWhereTheSupIsFinite(testCase)
        % Coverage is NOT 100% and that is the answer, not a gap: a cell of a conjugate is generally
        % an unbounded cone, so for many x the sup diverges, there is no maximiser and correctly no
        % region. The union of the regions must be EXACTLY the set where the sup is finite; a hole
        % and a false face are both defects, and they are opposite ones.
        %
        % DECIDED EXACTLY, NOT SAMPLED, and the two earlier attempts at this check are why. A
        % rejection sampler starves in a narrow cone at large radius and reported false "not
        % divergent" with NEGATIVE growth -- the tell that the measurement had run out of samples
        % rather than that the sup had settled. Replacing it with a grid of DIRECTIONS only moved
        % the problem into the grid's own resolution: measured, 501 directions manufacture 61 false
        % holes on this fixture while 2001 and 10001 both give none. That is the defect class
        % DECISIONS.md 2026-09-07 records -- a threshold read below the measurement's own
        % resolution -- so the sampling is gone and `supDiverges` decides it from the rank of N and
        % the recession cone, in integers.
            rng(20260913);
            holes = 0;  falseFaces = 0;  nCov = 0;  nUnc = 0;
            for c = [0 -1 -2]
                g = conjQ(envelopeCellRegionsTest.tri([0 0; 1 1; 0 1], c));
                for k = 1:g.nf
                    [N, m, q, dd, P] = envelopeCellRegionsTest.cellData(g, k);
                    if isempty(P), continue, end
                    faces = envelopeCellRegions(N, m, q, dd, P);
                    X = (rand(80,2) - 0.5) * 8;
                    for t = 1:size(X,1)
                        x = X(t,:);
                        diverges = envelopeCellRegionsTest.supDiverges(N, m, dd, P, x);
                        covered  = ~isempty(envelopeCellRegionsTest.locate(faces, x));
                        if covered
                            nCov = nCov + 1;
                            if diverges, falseFaces = falseFaces + 1; end
                        else
                            nUnc = nUnc + 1;
                            if ~diverges, holes = holes + 1; end
                        end
                    end
                end
            end
            testCase.verifyGreaterThan(nCov, 100, 'some x was covered');
            testCase.verifyGreaterThan(nUnc, 100, 'and some x was not');
            testCase.verifyEqual(holes, 0, 'an x where the sup is FINITE but no region covers it');
            testCase.verifyEqual(falseFaces, 0, 'a face claiming an x where the sup is INFINITE');
        end

        function theVERTEXSolveStacksROWSAndNotColumns(testCase)
        % RED-FIRST REGRESSION for the one real bug this routine had. Two matrices appear in the
        % vertex case and they are TRANSPOSES of each other: v solves the ROW system
        % [a_i'; a_j'] v = [-c_i; -c_j], while the normal cone combines the COLUMNS [a_i a_j].
        % Using one for the other still type-checks and still produces a plausible vertex.
        %
        % THE WITNESS, shrunk from a random failure by the no-overshoot check above. The cell is
        %       -s1 + 2 >= 0 ,   s2 >= 0 ,   s1 - s2 >= 0
        % with g(s) = |s|^2/4. Rows 1 and 3 are active at the apex (2,2) -- the column form puts it
        % at (2,0) instead, which silently drops that face's whole x2 term, so the face came out as
        % 2*x1 - 1 where it must be 2*x1 + 2*x2 - 2.
            P  = [-1 0 2; 0 1 0; 1 -1 0];
            N  = [1 0; 0 1];  m = [0; 0];  q = 0;  dd = 2;   % g(s) = (s'Ns/2)/dd = |s|^2/4
            faces = envelopeCellRegions(N, m, q, dd, P);

            v = [2; 2];
            isOurs = arrayfun(@(f) strcmp(f.who.kind,'vertex') && f.who.facet == 1 ...
                                   && f.who.facet2 == 3, faces);
            testCase.verifyEqual(nnz(isOurs), 1, 'the vertex of rows 1 and 3 owns exactly one face');
            f = faces(isOurs);

            % the face IS x |-> <v,x> - g(v), everywhere -- an identity, not a spot value
            X = [3 4; 2.9 3.7; 10 -2; 0 0; -5 7];
            want = X*v - (0.5*(v.'*N*v) + m.'*v + q)/dd;
            testCase.verifyEqual(ratQ.evalFace(f.fN, f.fD, X), want, 'AbsTol', 1e-12, ...
                'the vertex face must be <v,x> - g(v) with v = (2,2)');

            % and the x2 term is present, which is precisely what the bug removed
            testCase.verifyNotEqual(f.fN(9), 0, 'the face must depend on x2');

            % the apex really is where the two rows are active, and it is in the cell
            testCase.verifyEqual(P(1,1:2)*v + P(1,3), 0, 'row 1 active at v');
            testCase.verifyEqual(P(3,1:2)*v + P(3,3), 0, 'row 3 active at v');
            testCase.verifyGreaterThanOrEqual(P*[v;1], 0, 'and v lies in the cell');
        end

        function theTHREECandidateKindsAllOCCURAndTheCountsArePinned(testCase)
        % A routine with three branches needs every branch to fire on some fixture, or two of them
        % are untested code. Pinned per cell, because the shape of the answer is the measurement:
        % the three VERTEX cells of a conjugate carry AFFINE faces (N = 0), so they have neither a
        % stationary point nor a facet-interior maximiser and contribute exactly one face each --
        % which is why a routine that only ever saw the interior cell would look complete.
            g = conjQ(envelopeCellRegionsTest.tri([0 0; 1 0; 1 1], 0));
            testCase.verifyEqual(g.nf, 7, 'one interior, three edge and three vertex cells');
            kinds = cell(1, g.nf);
            for k = 1:g.nf
                [N, m, q, dd, P] = envelopeCellRegionsTest.cellData(g, k);
                faces = envelopeCellRegions(N, m, q, dd, P);
                kinds{k} = sort(arrayfun(@(f) f.who.kind(1), faces));
            end
            n = cellfun(@numel, kinds);
            testCase.verifyEqual(sort(n), [1 1 1 3 3 3 7], 'faces per cell');
            testCase.verifyEqual(sum(n), 19, 'nineteen faces over the seven cells');
            all1 = [kinds{:}];
            testCase.verifyEqual(nnz(all1 == 's'), 1, 'one STATIONARY candidate, in the interior cell');
            testCase.verifyEqual(nnz(all1 == 'f'), 6, 'six FACET candidates');
            testCase.verifyEqual(nnz(all1 == 'v'), 12, 'twelve VERTEX candidates');
        end

        function hIsCONTINUOUSAcrossEveryRegionBoundary(testCase)
        % h(x) = sup_s [<s,x> - g(s)] is a supremum of AFFINE functions of x, hence convex and
        % continuous wherever it is finite. So two regions meeting along a boundary must carry
        % faces that AGREE there -- that is what makes the subdivision define ONE function rather
        % than a collection of unrelated formulas, and it is the property a wrong region silently
        % breaks.
        %
        % ON CONSTRUCTED POINTS, because a random one is never on a boundary. Measured: 1050 random
        % draws found ZERO points lying in two regions at once -- the regions are disjoint and the
        % ties live on a measure-zero set, which is why the first version of this test, demanding a
        % sampled overlap, could only ever fail. The boundary is reached instead by BISECTING
        % between two points that landed in different regions.
            rng(20260913);
            worst = 0;  nBoundary = 0;
            for c = [0 -1]
                g = conjQ(envelopeCellRegionsTest.tri([0 0; 1 0; 1 1], c));
                for k = 1:g.nf
                    [N, m, q, dd, P] = envelopeCellRegionsTest.cellData(g, k);
                    if isempty(P), continue, end
                    faces = envelopeCellRegions(N, m, q, dd, P);
                    if numel(faces) < 2, continue, end
                    X = (rand(40,2) - 0.5) * 6;
                    idx = arrayfun(@(t) envelopeCellRegionsTest.locateIdx(faces, X(t,:)), ...
                                   (1:size(X,1)).');
                    for t = 1:size(X,1)
                        for u = t+1:size(X,1)
                            if idx(t) == 0 || idx(u) == 0 || idx(t) == idx(u), continue, end
                            lo = X(t,:);  hi = X(u,:);  a = idx(t);
                            for it = 1:60
                                mid = (lo + hi)/2;
                                if envelopeCellRegionsTest.locateIdx(faces, mid) == a
                                    lo = mid;
                                else
                                    hi = mid;
                                end
                            end
                            ia = envelopeCellRegionsTest.locateIdx(faces, lo);
                            ib = envelopeCellRegionsTest.locateIdx(faces, hi);
                            if ia == 0 || ib == 0 || ia == ib, continue, end
                            va = ratQ.evalFace(faces(ia).fN, faces(ia).fD, lo);
                            vb = ratQ.evalFace(faces(ib).fN, faces(ib).fD, hi);
                            worst = max(worst, abs(va - vb));
                            nBoundary = nBoundary + 1;
                        end
                    end
                end
            end
            testCase.verifyGreaterThan(nBoundary, 20, 'some boundary was actually constructed');
            testCase.verifyLessThan(worst, 1e-9, ...
                sprintf('two faces disagree across their shared boundary by %.3e', worst));
        end

        function RANDOMPointsNeverLieInTwoRegionsAtOnce(testCase)
        % The complement of the test above, as a MEASUREMENT rather than an assumption: the regions
        % are disjoint, so a point in two at once would mean the subdivision overlaps and does not
        % define a function. Zero is the expected count -- which is precisely why the continuity
        % above has to be checked on constructed points instead.
            rng(20260913);
            g = conjQ(envelopeCellRegionsTest.tri([0 0; 1 0; 1 1], 0));
            nOverlap = 0;  nPts = 0;
            for k = 1:g.nf
                [N, m, q, dd, P] = envelopeCellRegionsTest.cellData(g, k);
                if isempty(P), continue, end
                faces = envelopeCellRegions(N, m, q, dd, P);
                X = (rand(150,2) - 0.5) * 8;
                for t = 1:size(X,1)
                    nPts = nPts + 1;
                    in = arrayfun(@(f) all(f.R(:,1:2)*X(t,:).' + f.R(:,3) >= -1e-12), faces);
                    if nnz(in) < 2, continue, end
                    nOverlap = nOverlap + 1;
                    vals = arrayfun(@(f) ratQ.evalFace(f.fN, f.fD, X(t,:)), faces(in));
                    testCase.verifyLessThan(max(vals) - min(vals), 1e-9, ...
                        'two regions overlap AND disagree');
                end
            end
            testCase.verifyGreaterThan(nPts, 500, 'enough points were drawn to mean something');
            testCase.verifyEqual(nOverlap, 0, 'the regions are disjoint off their boundaries');
        end

        function theINPUTContractIsRefusedByName(testCase)
            P = [-1 0 2; 0 1 0; 1 -1 0];
            testCase.verifyError(@() envelopeCellRegions([1 0; 0 1], [0;0], 0, 2), ...
                'envelopeCellRegions:nargin');
            testCase.verifyError(@() envelopeCellRegions([1 2; 3 4], [0;0], 0, 2, P), ...
                'envelopeCellRegions:shape', 'a NON-symmetric N is not a quadratic form');
            testCase.verifyError(@() envelopeCellRegions([1 0; 0 1], [0;0], 0, 0, P), ...
                'envelopeCellRegions:denominator');
            testCase.verifyError(@() envelopeCellRegions([1 0; 0 1], [0;0], 0, -2, P), ...
                'envelopeCellRegions:denominator', 'a negative denominator flips every inequality');
        end

    end

    methods (Static)

        function f = locate(faces, x)
        % The first face whose region contains x, or empty.
            f = [];
            for k = 1:numel(faces)
                if all(faces(k).R(:,1:2)*x(:) + faces(k).R(:,3) >= -1e-12)
                    f = faces(k);  return
                end
            end
        end

        function [N, m, q, dd, P] = cellData(g, k)
        % Cell k of a QuaCon as the (N, m, q, dd, P) envelopeCellRegions takes.
            n = g.fN(k,:);  dd = g.fD(k);
            N = [n(5) n(6); n(6) n(7)];  m = [n(8); n(9)];  q = n(10);
            FC = g.FC{k};  P = zeros(0,3);
            for t = 1:size(FC,1)
                cc = g.EcQ(FC(t,1),:);
                if any(cc(1:3) ~= 0), P = zeros(0,3); return, end
                if FC(t,2) == 0, continue, end
                P(end+1,:) = FC(t,2) * cc(4:6); %#ok<AGROW>
            end
        end

        function tf = supDiverges(N, m, dd, P, x)
        % EXACT: is sup { <s,Y> - s'Ns/2 : s in P } equal to +infinity, with Y = dd*x - m?
        %
        % Along a recession direction r of P the objective behaves as t*<r,Y> - t^2*(r'Nr)/2, so it
        % diverges exactly when some r in rec(P) has r'Nr = 0 AND r'Y > 0. N is positive
        % SEMIdefinite, so r'Nr = 0 is equivalent to N r = 0 -- which turns the question into one
        % about the NULL SPACE of N, decided by its rank in integers:
        %
        %   rank 2 : the null space is {0}, so nothing diverges and the sup is always finite.
        %   rank 1 : the null space is one line, spanned by u = [n12; -n11] (or [1;0] when the
        %            first column vanishes). Both of its directions are tested.
        %   rank 0 : N is zero and the objective is LINEAR, so by Farkas the sup is finite exactly
        %            when Y lies in the cone spanned by the NEGATED rows of P.
            Y = dd*x(:) - m;
            A = P(:,1:2);
            if N(1,1)*N(2,2) - N(1,2)*N(2,1) ~= 0, tf = false; return, end
            if any(N(:) ~= 0)
                if N(1,1) ~= 0 || N(1,2) ~= 0, u = [N(1,2); -N(1,1)]; else, u = [1; 0]; end
                tf = (all(A*u >= 0) && u.'*Y > 0) || (all(A*(-u) >= 0) && (-u).'*Y > 0);
                return
            end
            tf = ~envelopeCellRegionsTest.inCone(-A.', Y);
        end

        function tf = inCone(G, Y)
        % Is Y in the cone generated by the COLUMNS of G? In the plane a cone is spanned by at most
        % two of its generators, so checking every single and every pair is complete.
            tf = all(Y == 0);
            if tf, return, end
            n = size(G,2);
            for i = 1:n
                gi = G(:,i);
                if all(gi == 0), continue, end
                if abs(gi(1)*Y(2) - gi(2)*Y(1)) < 1e-12 && gi.'*Y > 0, tf = true; return, end
            end
            for i = 1:n
                for j = i+1:n
                    M = [G(:,i), G(:,j)];
                    d = M(1,1)*M(2,2) - M(1,2)*M(2,1);
                    if d == 0, continue, end
                    mu = [M(2,2)*Y(1) - M(1,2)*Y(2); -M(2,1)*Y(1) + M(1,1)*Y(2)] / d;
                    if all(mu >= -1e-12), tf = true; return, end
                end
            end
        end

        function i = locateIdx(faces, x)
            i = 0;
            for k = 1:numel(faces)
                if all(faces(k).R(:,1:2)*x(:) + faces(k).R(:,3) >= -1e-12), i = k; return, end
            end
        end

        function S = feasiblePts(P, n)
            S = zeros(0,2);
            for sc = [1 4 16 64]
                X = (rand(40*n,2) - 0.5) * 2 * sc;
                keep = all(X*P(:,1:2).' + P(:,3).' >= 0, 2);
                S = [S; X(keep,:)]; %#ok<AGROW>
                if size(S,1) >= n, S = S(1:n,:); return, end
            end
        end

        function f = tri(V, c)
            q = [0 0 0 0, 2, 0, 2, 0, 0, 0];
            q(8) = q(8) - c;  q(9) = q(9) + c;
            f = QuaPol(V, [1 2 1; 2 3 1; 3 1 1], q, [1 0; 1 0; 1 0]);
        end

    end
end
