classdef envelopeCellPairsTest < matlab.unittest.TestCase
% Unit tests for envelopeCellPairs, the exact cell-pair enumeration and classification that opens
% the convex envelope's assembly.
%
% BUCKET: fast (conjQ of two triangles per fixture, then integer predicates; no symbolic call).
%
% TWO ORACLES, and they answer different questions.
%   * THE DEFINITION. The tie conic must really be the curve where the two faces agree, and the
%     classification must really say which face dominates. Both are checked by evaluating the
%     faces -- not by re-deriving the same expression a second time.
%   * THE PROBE IT WAS PROMOTED FROM. `.claude/assemblyProbe.m` measured these counts on these
%     fixtures before any of this was product code, so a count that moves is either a real change
%     or a transcription defect. It is the differential test CLAUDE.md section 6 asks for, with the
%     probe as the thing being replaced and therefore as the oracle.

    methods (Test)

        function theCOUNTSAgreeWithTheProbeItWasPromotedFrom(testCase)
        % Measured by `.claude/assemblyProbe.m` before this was a routine, on exactly these four
        % fixtures. The probe reports non-bridge and identical-face pairs separately and then adds
        % them into one "polyhedral" total; this routine counts the same way, so `nPoly` is what
        % the probe prints as polyhedral.
        %
        % THE FOURTH FIXTURE IS NOT THE UNIT SQUARE ON PURPOSE: three runs on one geometry answer
        % one geometry, and the 2x1 rectangle is anisotropic, so its dual mesh is not symmetric
        % under swapping the coordinates.
            %            name                        c    meet same poly bridge
            want = { 'unit square, c = 0 (CONVEX)',  0,    18,   2,   8,   10; ...
                     'unit square, c = -1',        -1,    23,   2,  12,   11; ...
                     'unit square, c = -2',        -2,    16,   2,   8,    8 };
            for k = 1:size(want,1)
                [gA, gB] = envelopeCellPairsTest.squarePieces(want{k,2});
                [pairs, st] = envelopeCellPairs(gA, gB);
                nm = want{k,1};
                testCase.verifyEqual(st.nPairs, 49, nm);
                testCase.verifyEqual(st.nMeet,  want{k,3}, ['meet in 2-D: ' nm]);
                testCase.verifyEqual(st.nSame,  want{k,4}, ['identical faces: ' nm]);
                testCase.verifyEqual(st.nPoly,  want{k,5}, ['polyhedral: ' nm]);
                testCase.verifyEqual(st.nBridge,want{k,6}, ['bridge: ' nm]);
                testCase.verifyEqual(numel(pairs), st.nMeet, ...
                    'every pair that meets in 2-D is returned, and nothing else is');
            end
            % the anisotropic fixture, whose numbers differ from all three above
            [gA, gB] = envelopeCellPairsTest.rectanglePieces(-2);
            [~, st] = envelopeCellPairs(gA, gB);
            testCase.verifyEqual([st.nMeet st.nSame st.nPoly st.nBridge], [20 1 12 8], ...
                'rectangle 2x1, c = -2');
        end

        function theTIECONICIsTheCurveWhereTheTwoFacesAGREE(testCase)
        % THE WEIGHTED-BASIS TRAP, pinned by an identity rather than by a comment. A face is stored
        % in the 10-wide cubic basis whose squared terms carry halves, so the conic of a difference
        % is [n5, 2*n6, n7, 2*n8, 2*n9, 2*n10] and NOT the raw entries. Reading it wrong names a
        % different curve and moves real numbers (DECISIONS.md 2026-09-06).
        %
        % The identity that says it is right: with dA, dB the two faces' denominators,
        %       q(s) = 2 * dA * dB * ( gA_i(s) - gB_j(s) )   for EVERY s.
        % Checked at points that have nothing to do with the cells, because the identity is about
        % the algebra and not about where the pair lives.
            [gA, gB] = envelopeCellPairsTest.squarePieces(-1);
            pairs = envelopeCellPairs(gA, gB);
            S = [0 0; 1 0; 0 1; 2 -3; -1.5 0.25; 7 11];
            nChecked = 0;
            for p = pairs
                i = p.ij(1);  j = p.ij(2);
                dA = gA.fD(i);  dB = gB.fD(j);
                want = 2 * dA * dB * ( ratQ.evalFace(gA.fN(i,:), dA, S) ...
                                     - ratQ.evalFace(gB.fN(j,:), dB, S) );
                got = QuaPar.evalConic(p.q, S);
                testCase.verifyEqual(got, want, 'AbsTol', 1e-9, 'RelTol', 1e-12, ...
                    sprintf('pair (%d,%d)', i, j));
                nChecked = nChecked + 1;
            end
            testCase.verifyGreaterThan(nChecked, 10, 'the identity was actually exercised');
        end

        function theWINNERReallyDominatesOnTheWHOLEPair(testCase)
        % What 'poly' CLAIMS is that one face is at least the other everywhere on the pair -- that
        % is what makes the pair carry a single quadratic. Sampled inside each such pair, which is
        % the definition and not a second derivation of the classification.
            [gA, gB] = envelopeCellPairsTest.squarePieces(-1);
            pairs = envelopeCellPairs(gA, gB);
            rng(20260913);
            nTested = 0;
            for p = pairs
                if ~strcmp(p.kind, 'poly'), continue, end
                if all(p.q == 0), continue, end          % identical faces: nothing to dominate
                S = envelopeCellPairsTest.pointsInside(p.P, 40);
                if isempty(S), continue, end
                dA = gA.fD(p.ij(1));  dB = gB.fD(p.ij(2));
                gAv = ratQ.evalFace(gA.fN(p.ij(1),:), dA, S);
                gBv = ratQ.evalFace(gB.fN(p.ij(2),:), dB, S);
                if strcmp(p.win, 'A')
                    testCase.verifyGreaterThanOrEqual(gAv - gBv, -1e-9, ...
                        sprintf('A was declared the winner on pair (%d,%d)', p.ij(1), p.ij(2)));
                else
                    testCase.verifyGreaterThanOrEqual(gBv - gAv, -1e-9, ...
                        sprintf('B was declared the winner on pair (%d,%d)', p.ij(1), p.ij(2)));
                end
                nTested = nTested + 1;
            end
            testCase.verifyGreaterThan(nTested, 3, 'some poly pair was actually sampled');
        end

        function aSAMPLEDSignChangeForcesTheBridgeClassification(testCase)
        % The sound direction, and the only one sampling can supply. Finding both signs on a pair
        % PROVES the tie changes sign there, so such a pair must be classified 'bridge'; failing to
        % find them proves nothing, so the converse is not asserted here (see the next test, which
        % measures why it would be wrong).
            [gA, gB] = envelopeCellPairsTest.squarePieces(-1);
            pairs = envelopeCellPairs(gA, gB);
            nSeen = 0;
            for p = pairs
                S = envelopeCellPairsTest.pointsInside(p.P, 200);
                if isempty(S), continue, end
                v = QuaPar.evalConic(p.q, S);
                if any(v > 1e-9) && any(v < -1e-9)
                    testCase.verifyEqual(p.kind, 'bridge', ...
                        sprintf('pair (%d,%d) changes sign, so it carries a bridge', ...
                                p.ij(1), p.ij(2)));
                    nSeen = nSeen + 1;
                end
            end
            testCase.verifyGreaterThan(nSeen, 0, 'some pair actually changed sign in the sample');
        end

        function MOSTBridgesAreTANGENTIALAndAPruneOnTheSignAloneWouldDropThem(testCase)
        % THE MEASUREMENT THAT PINS WHY FILTER 2 READS THE SECOND OUTPUT. `quadSignOnPoly` returns
        % `+1` for "c >= 0 on P", which ALLOWS c = 0 somewhere: the two faces then TOUCH without
        % crossing -- a TANGENTIAL tie, which carries a bridge face exactly as a transversal one
        % does. Only `strict` says the tie never reaches zero, and only that licenses a prune.
        % Using `s ~= 0` as the prune test dropped a real bridge once (DECISIONS.md 2026-09-06).
        %
        % This test exists because the FIRST version of it asserted that every bridge takes both
        % signs, and that is false: measured below, the tangential case is the MAJORITY, and on the
        % CONVEX fixture it is every single one. A prune on the sign alone would drop all ten there
        % -- which is the whole envelope -- so the stake is not a corner case.
            want = { 0, 0, 10;      % c,  transversal, tangential
                    -1, 3,  8; ...
                    -2, 2,  6 };
            for k = 1:size(want,1)
                [gA, gB] = envelopeCellPairsTest.squarePieces(want{k,1});
                pairs = envelopeCellPairs(gA, gB);
                nX = 0;  nT = 0;
                for p = pairs
                    if ~strcmp(p.kind, 'bridge'), continue, end
                    [s, strict] = ratQ.quadSignOnPoly(p.q, p.P);
                    testCase.verifyFalse(s ~= 0 && strict, ...
                        'a STRICTLY one-signed pair must never have been called a bridge');
                    if s == 0, nX = nX + 1; else, nT = nT + 1; end
                end
                nm = sprintf('c = %d', want{k,1});
                testCase.verifyEqual(nX, want{k,2}, ['transversal bridges, ' nm]);
                testCase.verifyEqual(nT, want{k,3}, ['TANGENTIAL bridges, ' nm]);
            end
        end

        function aPOLYPairNeverShowsBothSigns(testCase)
        % The other half, asserted at the strength sampling supports: a pair declared 'poly' claims
        % one face dominates throughout, so no sample may contradict it.
            [gA, gB] = envelopeCellPairsTest.squarePieces(-1);
            pairs = envelopeCellPairs(gA, gB);
            nTested = 0;
            for p = pairs
                if ~strcmp(p.kind, 'poly') || all(p.q == 0), continue, end
                S = envelopeCellPairsTest.pointsInside(p.P, 200);
                if isempty(S), continue, end
                v = QuaPar.evalConic(p.q, S);
                testCase.verifyFalse(any(v > 1e-9) && any(v < -1e-9), ...
                    sprintf('poly pair (%d,%d) must NOT change sign', p.ij(1), p.ij(2)));
                nTested = nTested + 1;
            end
            testCase.verifyGreaterThan(nTested, 3, 'some poly pair was actually sampled');
        end

        function IDENTICALFacesAreNotABridge(testCase)
        % When both primal pieces share the vertex a pair's cells are built from, the two
        % conjugates agree on the whole pair: q is identically zero, there is no tie curve, and
        % there is no bridge. Calling it one is the mistake `.claude/bridgeEnvelopeProbe.m` made on
        % its first run, and an all-zero q reaching a sign test would be asked which side of
        % nothing a point is on.
            [gA, gB] = envelopeCellPairsTest.squarePieces(-1);
            [pairs, st] = envelopeCellPairs(gA, gB);
            same = pairs(arrayfun(@(p) all(p.q == 0), pairs));
            testCase.verifyEqual(numel(same), st.nSame, 'the count and the pairs agree');
            testCase.verifyEqual(st.nSame, 2, 'this fixture has two');
            for p = same
                testCase.verifyEqual(p.kind, 'poly', 'an identical-face pair is never a bridge');
            end
        end

        function EVERYReturnedPairMeetsInTwoDimensions(testCase)
        % Filter 1 as an invariant of the OUTPUT rather than as a line in the routine. A pair
        % meeting only along an edge or at a point carries no region, so a face built on one would
        % have empty interior -- and the count says the filter does real work here: most pairs are
        % dropped.
            [gA, gB] = envelopeCellPairsTest.squarePieces(-1);
            [pairs, st] = envelopeCellPairs(gA, gB);
            for p = pairs
                testCase.verifyTrue(ratQ.feasible2(p.P, true), ...
                    sprintf('pair (%d,%d) was returned but does not meet in 2-D', ...
                            p.ij(1), p.ij(2)));
            end
            testCase.verifyLessThan(st.nMeet, st.nPairs, 'the filter drops something');
        end

        function aWHOLEPLANECellIsPairedNotSkipped(testCase)
        % The conjugate of a single point v with value c is <s,v> - c on the WHOLE plane: one cell
        % with no boundary rows. Until 2026-10-10 an empty row list was read as "no region" and the
        % cell was skipped, so a point paired with anything produced NO pairs and the envelope was
        % +inf everywhere. Paired with a triangle's conjugate, every cell of the triangle's must
        % meet it.
            gv = conjQ(QuaPol([1 -3], zeros(0,3), [0 0 0 0 0 0 0 0 0 -11], zeros(0,2)));
            gT = conjQ(envelopeCellPairsTest.tri([-2 0; 2 0; 0 2], 0));
            [pairs, stats] = envelopeCellPairs(gT, gv);
            testCase.verifyEqual(stats.nMeet, gT.nf, ...
                'every cell of the triangle meets the point''s whole-plane cell');
            testCase.verifyEqual(numel(pairs), gT.nf);
            [pairs2, ~] = envelopeCellPairs(gv, gv);
            testCase.verifyEqual(numel(pairs2), 1, 'a whole-plane cell meets itself');
            testCase.verifyEqual(pairs2(1).kind, 'poly');
        end

        function aCURVEDCellIsRefusedByNameRatherThanDroppedInSilence(testCase)
        % The conjugate of a CONVEX piece is polyhedral -- every cell boundary is the affine
        % condition that a KKT multiplier changes sign (CONJ_FIELD_PROOF.md 7.3) -- so a genuinely
        % conic cell edge means the caller passed something that is not one. Dropping it would
        % surface much later as a missing or wrong region; naming it surfaces it here.
            [V, E, Ec, F, sq] = QuaConTest.overlapParts();
            fA = [0 0 0 0 0 0 0 0 0 1];
            cur = [5 -1];                                    % edge 5 is the unit circle
            g = QuaCon(V, E, Ec, [fA; fA], [1;1], F, {[sq; cur]; [sq; cur]});
            [gA, ~] = envelopeCellPairsTest.squarePieces(-1);
            testCase.verifyError(@() envelopeCellPairs(g, gA), 'envelopeCellPairs:curvedCell');
            testCase.verifyError(@() envelopeCellPairs(gA, g), 'envelopeCellPairs:curvedCell');
        end

        function theINPUTContractIsRefusedByName(testCase)
        % The classification is integer arithmetic on fN/fD and EcQ, which only an exact conjugate
        % carries -- so a legacy double type is refused rather than silently read.
            [gA, ~] = envelopeCellPairsTest.squarePieces(-1);
            p = QuaPol([3 1 3 0 0 0]);
            testCase.verifyError(@() envelopeCellPairs(p, gA), 'envelopeCellPairs:type');
            testCase.verifyError(@() envelopeCellPairs(gA, p), 'envelopeCellPairs:type');
            testCase.verifyError(@() envelopeCellPairs(gA), 'envelopeCellPairs:nargin');
        end

    end

    methods (Static)

        function [gA, gB] = squarePieces(c)
        % The unit square cut along its diagonal into two triangles, with `c` the coefficient of
        % the linear form vanishing ON that diagonal -- so the two pieces agree there for every c
        % and the input is always a FUNCTION. c = 0 is convex; c < 0 is a kink of increasing depth,
        % which is the case a bridge exists for. DECISIONS.md 2026-09-05.
            gA = conjQ(envelopeCellPairsTest.tri([0 0; 1 0; 1 1], 0));
            gB = conjQ(envelopeCellPairsTest.tri([0 0; 1 1; 0 1], c));
        end

        function [gA, gB] = rectanglePieces(c)
            gA = conjQ(envelopeCellPairsTest.tri([0 0; 2 0; 2 1], 0));
            gB = conjQ(envelopeCellPairsTest.tri([0 0; 2 1; 0 1], c));
        end

        function f = tri(V, c)
            q = [0 0 0 0, 2, 0, 2, 0, 0, 0];
            q(8) = q(8) - c;  q(9) = q(9) + c;
            f = QuaPol(V, [1 2 1; 2 3 1; 3 1 1], q, [1 0; 1 0; 1 0]);
        end

        function S = pointsInside(P, n)
        % Points satisfying every row of P, found by rejection from a box that grows until it
        % catches some -- the cells are unbounded cones, so a fixed box can miss a pair entirely
        % and return "no counterexample" for the wrong reason.
            S = zeros(0,2);
            for scale = [1 4 16 64 256]
                X = (rand(20*n, 2) - 0.5) * 2 * scale;
                keep = all(X * P(:,1:2).' + P(:,3).' >= 0, 2);
                S = [S; X(keep,:)]; %#ok<AGROW>
                if size(S,1) >= n, S = S(1:n,:); return, end
            end
        end

    end
end
