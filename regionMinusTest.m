classdef regionMinusTest < matlab.unittest.TestCase
% The set-difference half of region's unit tests, split out from `regionUnitTest` on cost alone.
%
% BUCKET: normal. MEASURED 2026-08-31: this one test is 143 s, against 18 s for the nine tests
% left in regionUnitTest, because `region.minus` does symbolic work per facet pair and the fixture
% has eight of them. The fast bucket's whole budget is five minutes, so one test taking half of it
% belongs in the next bucket up -- not in a weaker form here.
%
% The assertion is the definition of A \ B and nothing else: no piece may contain a point of B,
% every piece must stay inside A, and every point of A outside B must be in some piece. The first
% of those is the one that would corrupt a Step 3 partition silently.

    properties (Constant)
        X = sym('x')
        Y = sym('y')
    end

    methods (Static)
        function v = vars()
            v = [regionMinusTest.X, regionMinusTest.Y];
        end
    end

    methods (Test)

        function minusIsTheSetDifference(testCase)
        % A \ B, as pieces. Two directions, and the first is the one that would corrupt a Step 3
        % partition silently: no piece may contain a point of B.
            v = regionMinusTest.vars();
            A = plqCheck.boxRegion([0 0], [4 4], v);
            B = plqCheck.boxRegion([1 1], [2 2], v);
            pieces = A - B;
            testCase.verifyGreaterThanOrEqual(numel(pieces), 1, 'A \ B should not be empty here');
            box = plqCheck.regionBox(A);

            for i = 1:numel(pieces)
                Pi = plqCheck.regionSample(pieces(i), box, 31, 1e-6);
                if isempty(Pi), continue, end
                testCase.verifyTrue(all(plqCheck.inRegion(A, Pi, 1e-5)), sprintf( ...
                    'piece %d of A \\ B leaves A', i));
                inB = plqCheck.inRegion(B, Pi, -1e-5);
                bad = find(inB, 1);
                testCase.verifyFalse(any(inB), sprintf( ...
                    'piece %d of A \\ B contains %d points strictly inside B (first (%g,%g))', ...
                    i, sum(inB), Pi(max(bad,1),1), Pi(max(bad,1),2)));
            end

            % COVER: a point of A strictly outside B must be in some piece.
            Q = plqCheck.regionSample(A, box, 31, 1e-6);
            outside = Q(~plqCheck.inRegion(B, Q, 1e-5), :);
            testCase.verifyGreaterThan(size(outside,1), 0, 'nothing to cover -- the fixture is wrong');
            covered = false(size(outside,1),1);
            for i = 1:numel(pieces)
                covered = covered | plqCheck.inRegion(pieces(i), outside, 1e-5);
            end
            miss = find(~covered, 1);
            testCase.verifyTrue(all(covered), sprintf( ...
                'A \\ B misses %d of %d points of A that are outside B (first (%g,%g))', ...
                sum(~covered), numel(covered), outside(max(miss,1),1), outside(max(miss,1),2)));
        end

    end
end
