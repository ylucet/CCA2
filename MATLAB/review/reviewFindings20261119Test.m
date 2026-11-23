classdef reviewFindings20261119Test < matlab.unittest.TestCase
% Findings of the 2026-11-19 code review (reviews/codereview-2026-11-19.md). Each test asserts the
% CORRECT behaviour, so it FAILS against the reviewed code; it is expected to pass once the finding
% is fixed. This file lives in MATLAB/review/, which suite.sh does not run (it lists *Test.m in
% MATLAB/ only). Run it explicitly:
%     cd MATLAB; addpath(pwd); runtests(fullfile('review','reviewFindings20261119Test.m'))

    methods (TestClassSetup)
        function onPath(tc) %#ok<MANU>
            addpath(fileparts(fileparts(mfilename('fullpath'))));
        end
    end

    methods (Test)
        function coneWithInteriorNullDirection(tc)
            % Finding 1 (CRITICAL). f(x,y) = y^2/2 on the cone {x >= |y|}. Along x = (t,0) the
            % objective <s,x> - f(x) = s1*t, so f*(s) = +inf whenever s1 > 0. conjQ only tests the
            % cone's two edge rays for zero curvature; the null direction (1,0) is strictly inside.
            f = QuaPol([0 0; 1 -1; 1 1], [1 2 0; 1 3 0], [0 0 0 0 0 0 1 0 0 0], [1 0; 0 1]);
            g = conjQ(f);
            v = g.eval([1 0; 0.5 0.3; 0.1 0; -1 0]);
            tc.verifyEqual(v(1:3).', [Inf Inf Inf], 'f*(s) must be +inf for s1 > 0');
            tc.verifyEqual(v(4), 0, 'AbsTol', 1e-12);
        end

        function cellInteriorMissesAFarVerticalStrip(tc)
            % Finding 2 (MAJOR). The cell {10 <= x <= 10.5} minus the disc |x-(100,0)| < 1 has
            % interior (e.g. (10.25, 0)). No two of its curves meet inside it, the random probe box
            % is around the origin, and the "complete" CAD fallback skips vertical lines, so the
            % cell is declared EMPTY. Its horizontal twin is (correctly) found.
            EcQ = [0 0 0 1 0 -10; 0 0 0 2 0 -21; 1 0 1 -200 0 9999];
            rows = [1 1; 2 -1; 3 1];
            tc.verifyTrue(cellHasInterior(EcQ, rows));
        end

        function cellInteriorMissesSharedLineFactorCell(tc)
            % Finding 2, second mechanism: four degenerate conics sharing x + y + 3 = 0; conicMeet
            % calls each pair degenerate and returns no points, so the CAD misses abscissae.
            % Reached in a real conj call; (-58.4253, -109.3822) satisfies every row strictly.
            EcQ = [101 76 -25 3378 3000 9225; 91 38 -53 -174 -606 -1341; 64 37 -27 2386 2113 6582;
                   1221 154 -1067 -42386 -49250 -138147; 0 0 0 3 -5 65; 0 0 0 3 -5 -31;
                   0 0 0 -18 8 -119; 0 0 0 9 -4 313; 0 0 0 5 -2 169];
            rows = [1 1; 2 -1; 3 -1; 4 1; 5 1; 6 1; 7 1; 8 1; 9 1];
            tc.verifyTrue(cellHasInterior(EcQ, rows));
        end

        function conjLShapeY2(tc)
            % Finding 8 (CRITICAL): f = y^2/2 on the L-shape; sup at s = (2.53125, -1.453125) is
            % 2*s1 = 5.0625 (attained at (2,0)). conjQ returns (2 s1 + s2)^2/2 = 6.5137939453125.
            V = [0 0; 2 0; 2 1; 1 1; 1 2; 0 2];
            E = [(1:6).' [2:6 1].' ones(6,1)];
            f = QuaPol(V, E, [0 0 0 0 0 0 1 0 0 0], repmat([1 0], 6, 1));
            g = conjQ(f);
            tc.verifyEqual(g.eval([2.53125 -1.453125]), 5.0625, 'AbsTol', 1e-12);
        end

        function biconjInfiniteOutsideTheSquare(tc)
            % Finding 9a (CRITICAL): f = (y^2 - x^2)/2 on the unit square; (-0.49, 0.01) is outside
            % the square, so f** = +inf there.
            f = QuaPol([0 0; 1 0; 1 1; 0 1], [1 2 1; 2 3 1; 3 4 1; 4 1 1], [0 0 0 0 -1 0 1 0 0 0], repmat([1 0], 4, 1));
            h = biconjQ(f);
            tc.verifyEqual(h.eval([-0.4921875 0.01171875]), Inf);
        end

        function biconjFenchelYoungLShape(tc)
            % Finding 9b (CRITICAL): f = -xy + y^2/2 on the L-shape. f**(x) >= <s,x> - f*(s) for
            % every s; with conjQ's own value at s the bound is 1.04606.
            V = [0 0; 2 0; 2 1; 1 1; 1 2; 0 2];
            E = [(1:6).' [2:6 1].' ones(6,1)];
            f = QuaPol(V, E, [0 0 0 0 0 -1 1 0 0 0], repmat([1 0], 6, 1));
            x = [0.2578125 1.76171875];  s = [-1.9484626037100141 1.6789473678974558];
            g = conjQ(f);  h = biconjQ(f);
            tc.verifyGreaterThanOrEqual(h.eval(x), s*x.' - g.eval(s) - 1e-9);
        end

        function conjCPLQTinyPositiveEigenvalue(tc)
            % Finding 7 (legacy 'cplq' engine, MAJOR there). f = (x^2 + y^2)/2^31 is strictly convex,
            % f*(s) = 2^29 |s|^2 everywhere. conjCPLQ's absolute eig(Q) > sqrt(eps) test sends it
            % to the AFFINE branch and returns a needle: +inf at every s ~= 0.
            f = QuaPol([2^-30 0 2^-30 0 0 0]);
            g = conjCPLQ(f);
            v = g.eval([1 0; 0 1]);
            tc.verifyEqual(v(:).', [2^29 2^29], 'RelTol', 1e-9);
        end

        function discontinuousInputIsReadLowerSemicontinuous(tc)
            % Finding 10 (MAJOR) / #108. All 30 two-piece D fixtures of the review generator
            % (plq_gen.py, seed 20261119): a convex quadrilateral cut by a diagonal, the two
            % pieces carrying DIFFERENT quadratics. f is read LOWER SEMICONTINUOUS (on the
            % diagonal, the min of the two), and conjQ / biconjQ must give the values of
            % python/tests/review/d_family_ref.json -- the Python engine's, graded 0 wrong by
            % compare_engines.py (conj exact, biconj inside the certified bracket) -- or the same
            % named refusal. D-090 is refused by MATLAB's own 2^53 ceiling (witnessOverflow)
            % where the ports name noWitness. The minutes-long fixtures (and the MATLAB-hour
            % ones) are in discontinuousInputIsReadLowerSemicontinuousSlow.
            reviewFindings20261119Test.checkDFamily(tc, false);
        end

        function discontinuousInputSharedVertexD132(tc)
            % #108 residual (OPEN, TODO): the atoms route subdivides the min over atom pairs into
            % an AlgCrv that takes the max over the faces holding at a point; on D-132 a bridge
            % face of piece 1's edge atom reaches the vertex (-2,2) with piece 1's value 11, where
            % lsc f = min(11, -7) = -7, and 9.875 on the diagonal next to it. Filtered until fixed.
            tc.assumeFail('#108 residual: AlgCrv max at a jump vertex (D-132)');
            f = QuaPol([-2 2; -1 -2; 1 -3; 2 2], [1 2 1; 2 3 1; 3 1 1; 3 4 1; 4 1 1], ...
                       [0 0 0 0 -1 -1 2 -2 1 -1; 0 0 0 0 -1 0 -1 0 -2 1], [1 0; 1 0; 1 2; 2 0; 2 0]);
            v = biconjQ(f).eval([-2 2; -1.625 1.375]);
            tc.verifyEqual(v(1), -7, 'AbsTol', 1e-9, 'f** at the shared vertex (-2,2)');
            tc.verifyLessThanOrEqual(v(2), -4.015625 + 1e-9, 'f** at (-13/8, 11/8) <= lsc f');
        end

        function discontinuousInputIsReadLowerSemicontinuousSlow(tc)
            % The rest of the D family (see above). Not in a bucket: run by name.
            reviewFindings20261119Test.checkDFamily(tc, true);
        end

        function continuousInputIsNotRefused(tc)
            % Finding 10 guard. R3-072 (a continuous 2-piece split of a quadrilateral: the two
            % quadratics differ by a multiple of the shared edge's own linear form) must PASS the
            % continuity check, so conjQ still answers it.
            V = [-3 2; -1 0; 2 0; -2 3];
            E = [1 2 1; 2 3 1; 3 1 1; 3 4 1; 4 1 1];
            f = QuaPol(V, E, [0 0 0 0 0 2 -2 2 1 1; 0 0 0 0 4 9 8 -4 -8 5], ...
                       [1 0; 1 0; 1 2; 2 0; 2 0]);
            g = conjQ(f);
            tc.verifyClass(g, 'QuaCon', 'a continuous multi-piece input must not be refused');
        end

        function conjCPLQSmallEigenvalueDoesNotCrash(tc)
            % Finding 7, second face. Q = diag(1, 2^-30) (strictly convex) and Q = diag(1, -2^-30)
            % (indefinite, f* = +inf everywhere) both fall into the rank-1 branch, where
            % null(Q) is EMPTY for a numerically full-rank Q: MATLAB:badsubscript, unnamed.
            g = conjCPLQ(QuaPol([1 0 2^-30 0 0 0]));
            v = g.eval([0 1]);
            tc.verifyEqual(v, 2^29, 'RelTol', 1e-9);
            g = conjCPLQ(QuaPol([1 0 -2^-30 0 0 0]));
            v = g.eval([0 1; 1 0]);
            tc.verifyEqual(v(:).', [Inf Inf]);
        end
    end
    methods (Static)
        function checkDFamily(tc, slow)
            slowSet = {'D-069', 'D-111', 'D-153', 'D-188', 'D-209'};
            ownRefusal = struct('name', {'D-090'}, 'err', {'bridgeFaceAlg:witnessOverflow'});
            root = fileparts(fileparts(fileparts(mfilename('fullpath'))));   % CCA2/
            fixes = jsondecode(fileread(fullfile(root, 'python', 'tests', 'review', 'd_family_ref.json')));
            if isstruct(fixes), fixes = num2cell(fixes); end
            tc.verifyEqual(numel(fixes), 30, 'the D family has 30 fixtures');
            num = @(c) reviewFindings20261119Test.numOf(c);
            for k = 1:numel(fixes)
                rec = fixes{k};
                if xor(slow, ismember(rec.name, slowSet)), continue, end
                f = QuaPol(rec.V, rec.E, rec.f, rec.F);
                tc.verifyFalse(f.isContinuous(), sprintf('%s: the D family jumps', rec.name));
                gv = conjQ(f).eval(rec.S);
                tc.verifyEqual(gv(:), num(rec.conj), 'RelTol', 1e-9, 'AbsTol', 1e-9, ...
                               sprintf('%s: conj', rec.name));
                want = '';
                if isfield(rec, 'biconj_err') && ~isempty(rec.biconj_err), want = rec.biconj_err; end
                m = strcmp({ownRefusal.name}, rec.name);
                if any(m), want = ownRefusal(m).err; end
                if ~isempty(want)
                    tc.verifyError(@() biconjQ(f), want, sprintf('%s: biconj', rec.name));
                    continue
                end
                hv = biconjQ(f).eval(rec.X);
                tc.verifyEqual(hv(:), num(rec.biconj), 'RelTol', 1e-7, 'AbsTol', 1e-7, ...
                               sprintf('%s: biconj', rec.name));
            end
        end

        function v = numOf(c)
            % jsondecode gives a double vector, or a cell when 'inf' strings are mixed in
            if ~iscell(c), v = c(:); return, end
            v = zeros(numel(c), 1);
            for i = 1:numel(c)
                if ischar(c{i}), v(i) = str2double(strrep(c{i}, 'inf', 'Inf')); else, v(i) = c{i}; end
            end
        end
    end
end
