classdef maxQuaParTest < matlab.unittest.TestCase
    % Tests for maxQuaPar (Step 3 of the 'cplq' pipeline): pointwise maximum of two full-domain
    % QuaPar objects. Centred on the f(x,y)=xy, three-convex-edge-triangle example from
    % /home/ylucet/CCA2/3-edge.tex: T=conv{(0,0),(3,3),(1,2)}, split by [COAP, Appendix A.5] into
    % T1=conv{(0,0),(1,2),(2,2)} and T2=conv{(1,2),(2,2),(3,3)}, each convexified to a rank-1 PSD
    % envelope (q1, q2 below) whose conjugates g1=conjPieceCPLQ(q1), g2=conjPieceCPLQ(q2) are the
    % inputs to maxQuaPar.
    %
    % An earlier version of this note (and of maxQuaPar.m's splitCell) claimed that two of the
    % comparison boundaries between g1 and g2 are genuine hyperbolas -- checking only the
    % quadratic-part discriminant delta=b^2-4ac (equivalently rank of the 2x2 submatrix
    % [a b/2;b/2 c]) of their difference, which is positive for both. That is necessary but not
    % sufficient for a genuine (irreducible) hyperbola: the FULL 3x3 conic discriminant Delta is
    % exactly zero for both, so each difference is in fact a DEGENERATE conic -- a pair of two
    % real straight lines -- fully representable by QuaPar, not a hyperbola. See the corrected
    % Proposition 1 in 3-edge.tex for the exact factorizations.

    methods (Test)
        function stepTwoProducesSixFaceEnvelopeConjugatesForBothSubTriangles(testCase)
            [g1, g2] = maxQuaParTest.buildG1G2();
            testCase.verifyClass(g1, 'QuaPar');
            testCase.verifyClass(g2, 'QuaPar');
            testCase.verifyEqual(g1.nf, 6);   % rank-1 PSD triangle conjugate: 3 edge strips + 3
            testCase.verifyEqual(g2.nf, 6);   % vertex cones (see conjPieceCPLQ), no tied-vertex
        end                                   % degeneracy for this q1/q2 (would give nf=5).

        function boundaryConicsAreDegenerateLinePairsNotHyperbolas(testCase)
            % The mathematical fact the whole fix depends on: find g1's shared-edge piece (the
            % edge-strip for the primal edge (1,2)-(2,2), shared by T1 and T2 but a convex edge of
            % NEITHER -- eq (eq:E1shared) in 3-edge.tex) and g2's "outer" and "far" pieces (eq
            % (eq:E2outer)/(eq:E2far)) by their known closed-form coefficients, then verify each
            % pairwise difference has delta=b^2-4ac>0 (matching 3-edge.tex's original, correct
            % discriminant computation) but the full conic discriminant Delta~=0 (the correction).
            [g1, g2] = maxQuaParTest.buildG1G2();

            s2v2 = sqrt(2);
            e1Shared = [ (0.75+s2v2/2), 0, 0, -s2v2, 2, 0 ];                 % eq (eq:E1shared), 3rd line
            e2Outer  = [ 1, 0.5, 0.25, -1.5, 0.75, 1.125 ];                  % eq (eq:E2outer)
            e2Far    = [ 0.5, 0.5, 0.5, 0, 0, 0 ];                          % eq (eq:E2far)

            kShared = maxQuaParTest.findRow(g1.f, e1Shared, testCase);
            lOuter  = maxQuaParTest.findRow(g2.f, e2Outer, testCase);
            lFar    = maxQuaParTest.findRow(g2.f, e2Far, testCase);

            [deltaA, DeltaA] = maxQuaParTest.conicInvariants(g1.f(kShared,:) - g2.f(lOuter,:));
            [deltaB, DeltaB] = maxQuaParTest.conicInvariants(g1.f(kShared,:) - g2.f(lFar,:));

            testCase.verifyGreaterThan(deltaA, 0);              % Boundary A: hyperbolic TYPE...
            testCase.verifyEqual(DeltaA, 0, 'AbsTol', 1e-9);    % ...but degenerate (a line pair).
            testCase.verifyGreaterThan(deltaB, 0);              % Boundary B: same story.
            testCase.verifyEqual(DeltaB, 0, 'AbsTol', 1e-9);
        end

        function splitCellAcceptsGenuineNonDegenerateParabola(testCase)
            % Regression test for a bug found via a randomized stress test across many f(x,y)=xy
            % triangles (not just the one hard-coded 3-edge.tex example used elsewhere in this
            % file): splitCell's degeneracy guard used to reject ANY g1-face/g2-face boundary
            % whose full 3x3 discriminant Delta was nonzero, mislabelling it "a genuine
            % ellipse/hyperbola" -- but delta=b^2-4ac==0 (parabolic TYPE) with Delta~=0 is a
            % genuine, representable PARABOLA, which QuaPar's curved Ec edges exist to handle.
            %
            % T=(0,0),(9.31,7.63),(3.80,7.40) is a plain triangle with all 3 edges convex for
            % u1*u2 (so f(x,y)=xy is genuinely nonconvex on it -- the only case that matters, per
            % the reduction argument: remove the affine part, rotate any indefinite quadratic to
            % xy, triangulate). One of its g1-face/g2-face boundaries has delta~0 (parabolic type)
            % but Delta significantly nonzero (genuinely non-degenerate) -- maxQuaPar(g1,g2) used
            % to crash on this legitimate input with maxQuaPar:notDegenerate. (A different random
            % triangle, T=(0,0),(7.02,0.67),(8.43,7.63), hits the same guard bug but ALSO exposes a
            % separate, deeper Step-1 correctness issue unrelated to this guard -- see session notes
            % -- so it is deliberately not used here to keep this test isolated to the one bug it
            % documents.)
            %
            % g1,g2 are a frozen fixture (see frozenG1G2's header): Step 1 now correctly splits
            % this T into 4 pieces (2 rational, unconjugatable yet -- see DESIGN.md's Part 2c), so
            % buildG1G2ForTriangle can no longer build g1,g2 for it live; this test is about
            % maxQuaPar's own splitCell logic, not Step 1's tightness.
            T = [0 0; 9.31 7.63; 3.80 7.40];
            [g1, g2] = maxQuaParTest.frozenG1G2('splitCell');

            foundGenuineParabola = false;
            for k = 1:g1.nf
                for l = 1:g2.nf
                    [delta, Delta] = maxQuaParTest.conicInvariants(g1.f(k,:) - g2.f(l,:));
                    if abs(delta) < 1e-6 && abs(Delta) > 1e-3
                        foundGenuineParabola = true;
                    end
                end
            end
            testCase.verifyTrue(foundGenuineParabola, ...
                'expected at least one genuinely non-degenerate (parabolic-type) boundary pair');

            g = maxQuaPar(g1, g2);   % used to throw maxQuaPar:notDegenerate here
            testCase.verifyClass(g, 'QuaPar');

            testPts = [1 1; 3 2; 5 3; 2 6; -1 3];
            for i = 1:size(testPts,1)
                s = testPts(i,:);
                testCase.verifyEqual(g.eval(s), maxQuaParTest.supBilinearOverPoly(s, T), ...
                    'AbsTol', 1e-6, sprintf('s=(%.4f,%.4f)', s(1), s(2)));
            end
        end

        function dedupHitsMergesCrossingsAtACellCorner(testCase)
            % Regression test for a bug found via the same randomized triangle stress test:
            % splitCell's dedupHits used tol=sqrt(eps) (~1.5e-8) to merge boundary crossings that
            % are the SAME physical point computed via two different cell edges' independent
            % quadratic-root arithmetic. That arithmetic can disagree by ~1e-7 for a genuinely
            % coincident hit (the curve crossing exactly at a cell corner) -- the same
            % cross-arithmetic noise floor already documented (with a 1e-6 absolute tolerance) in
            % assemblePieces' global vertex merge -- so two hits ~7e-8 apart at a corner were
            % wrongly kept distinct, inflating a genuine 2-crossing split into 3 and tripping
            % splitCell's "expected exactly 2 boundary crossings" assertion on legitimate input.
            %
            % T=(0,0),(2.11,1.43),(8.84,4.50) is a plain triangle with all 3 edges convex for
            % u1*u2; maxQuaPar(g1,g2) used to crash on it with maxQuaPar:internal.
            %
            % g1,g2 are a frozen fixture (see frozenG1G2's header): Step 1 now correctly splits
            % this T into 4 pieces (2 rational, unconjugatable yet), so buildG1G2ForTriangle can
            % no longer build g1,g2 for it live; this test is about maxQuaPar's own dedupHits
            % logic, not Step 1's tightness.
            T = [0 0; 2.11 1.43; 8.84 4.50];
            [g1, g2] = maxQuaParTest.frozenG1G2('dedup');

            g = maxQuaPar(g1, g2);   % used to throw maxQuaPar:internal here
            testCase.verifyClass(g, 'QuaPar');

            testPts = [1 1; 3 2; 0.5 0.3; 5 3; -2 1];
            for i = 1:size(testPts,1)
                s = testPts(i,:);
                testCase.verifyEqual(g.eval(s), maxQuaParTest.supBilinearOverPoly(s, T), ...
                    'AbsTol', 1e-6, sprintf('s=(%.4f,%.4f)', s(1), s(2)));
            end
        end

        function assemblePiecesResolvesNearDuplicateApexCluster(testCase)
            % Regression test for the "assemblePieces: boundary edge (r,c) has no matching
            % neighbour" crash reported in the session handoff (see maxQuaPar.m's assemblePieces
            % HISTORY for the full fix). Root cause: a fan of many pieces meeting near one apex
            % produced several near-duplicate vertex computations, some agreeing only to ~1e-5
            % (cross-arithmetic noise across two different formulas for the "same" point) while
            % OTHER, genuinely distinct vertices in the same tiling were only ~5e-3 apart -- no
            % single coordinate-clustering tolerance can separate these two cases. Fixed by
            % matching half-edges directly by geometry (never by pre-clustering vertices into a
            % shared index) and deriving vertex identity afterwards, via union-find, purely from
            % those confirmed matches.
            %
            % T=(6.0365,4.9504),(9.8960,6.3015),(1.4908,3.3753): the exact reproduction case from
            % the session handoff; used to throw maxQuaPar:internal here.
            %
            % g1,g2 are a frozen fixture (see frozenG1G2's header): Step 1 now correctly splits
            % this T into 4 pieces (2 rational, unconjugatable yet), so buildG1G2ForTriangle can
            % no longer build g1,g2 for it live; this test is about maxQuaPar's own assemblePieces
            % logic, not Step 1's tightness.
            T = [6.0365 4.9504; 9.8960 6.3015; 1.4908 3.3753];
            [g1, g2] = maxQuaParTest.frozenG1G2('apexCluster');

            g = maxQuaPar(g1, g2);   % used to throw maxQuaPar:internal here
            testCase.verifyClass(g, 'QuaPar');

            testPts = [T(1,:); T(2,:); T(3,:); mean(T,1); mean(T,1)+[1 -2]];
            for i = 1:size(testPts,1)
                s = testPts(i,:);
                testCase.verifyEqual(g.eval(s), maxQuaParTest.supBilinearOverPoly(s, T), ...
                    'AbsTol', 1e-6, sprintf('s=(%.4f,%.4f)', s(1), s(2)));
            end
        end

        function checkOrphanHalfEdgesDropsProvablyDegenerateOrphanEdges(testCase)
            % Regression test for the residual "assemblePieces: ... has no matching neighbour"
            % crash left open by assemblePiecesResolvesNearDuplicateApexCluster's fix (see this
            % session's entry in maxQuaPar.m's assemblePieces HISTORY, and checkOrphanHalfEdges'
            % own header). Root cause: a genuinely AMBIGUOUS 3-way vertex cluster -- three
            % different pieces' edges all mutually within matchHalfEdges' tolPos, at a scale
            % (~1e-4) too coarse to be cross-arithmetic noise but too fine to be a distinct
            % feature -- so the best-first greedy matcher can pair up only 2 of the 3, always
            % orphaning the third. Diagnosis showed the orphaned edge's own two endpoints ALWAYS
            % resolve to the very same global vertex via the OTHER, confirmed matches on its own
            % piece's boundary, i.e. it is provably zero-length once the rest of the topology is
            % resolved -- safe to drop rather than error. Fixed by checkOrphanHalfEdges.
            %
            % 5 randomly-found near-degenerate (nearly-collinear) triangles that all used to throw
            % maxQuaPar:internal here; verified end-to-end against ground truth at several points.
            %
            % g1,g2 for each T are a frozen fixture (see frozenG1G2's header): Step 1 now correctly
            % splits each of these into 4 pieces (2 rational, unconjugatable yet), so
            % buildG1G2ForTriangle can no longer build g1,g2 for them live; this test is about
            % maxQuaPar's own checkOrphanHalfEdges logic, not Step 1's tightness.
            Ts = { [8.5697 2.6142; 5.0151 1.8051; 1.3296 0.9185]
                   [7.4090 3.9129; 4.3476 2.4669; 1.6019 1.1161]
                   [8.1673 7.7777; 5.0209 4.8767; 0.6393 0.7022]
                   [6.6140 9.9669; 6.0687 8.8885; 4.8459 6.3676]
                   [2.1767 3.3847; 0.5909 2.7773; 7.0977 5.1590] };
            for ti = 1:numel(Ts)
                T = Ts{ti};
                [g1, g2] = maxQuaParTest.frozenG1G2(sprintf('orphan%d', ti));
                g = maxQuaPar(g1, g2);   % used to throw maxQuaPar:internal here
                testCase.verifyClass(g, 'QuaPar');

                c = mean(T,1);
                testPts = [T(1,:); T(2,:); T(3,:); c; c+[1 -2]; c+[-1 3]];
                for i = 1:size(testPts,1)
                    s = testPts(i,:);
                    testCase.verifyEqual(g.eval(s), maxQuaParTest.supBilinearOverPoly(s, T), ...
                        'AbsTol', 1e-6, sprintf('T#%d s=(%.4f,%.4f)', ti, s(1), s(2)));
                end
            end
        end

        function insertPassthroughVerticesDropsNearDuplicateCrossingPoint(testCase)
            % Regression test for the still-open silent-wrong-answer issue from the prior session's
            % handoff: g.eval used to return Inf over a real, comfortably-covered region of the
            % plane instead of the correct finite value. Root cause: insertPassthroughVertices'
            % "already a vertex" check shared its tolerance (1e-7) with onOpenSegment/onOpenRay's
            % own edge-matching tolerance, but a genuine cross-arithmetic-noise gap between the
            % ORIGINAL face vertex used as a passthrough candidate and the ALREADY-PRESENT cell
            % vertex it geometrically coincides with can be ~3e-5 -- too coarse for that shared
            % 1e-7. This face (an unbounded piece whose two rays point in the SAME direction, a
            % degenerate "parallel-strip" shape rather than a normal wedge) ended up with an extra,
            % near-zero-length sliver edge whose line equation was dominated by floating-point noise
            % in the tiny direction vector, wrongly excluding a real chunk of the plane from
            % QuaPar.eval's exact (no-tolerance) membership test. Fixed by decoupling the
            % "already a vertex" pre-check into its own, wider tolerance (tolSnap=1e-4) -- see
            % insertPassthroughVertices' own header for why the shared-tolerance fix was tried first
            % and reverted (it broke two OTHER regression tests in this file with a different wrong
            % topology, each at a different noise scale).
            %
            % T=(7.8665,4.6784),(2.6908,1.9477),(0.3892,0.7130): a randomly-found repro triangle;
            % s=(3.380265,-0.644943) used to give g.eval==Inf here (true value ~2.6946).
            %
            % g1,g2 are a frozen fixture (see frozenG1G2's header): Step 1 now correctly splits
            % this T into 4 pieces (2 rational, unconjugatable yet), so buildG1G2ForTriangle can
            % no longer build g1,g2 for it live; this test is about maxQuaPar's own
            % insertPassthroughVertices logic, not Step 1's tightness.
            T = [7.8665 4.6784; 2.6908 1.9477; 0.3892 0.7130];
            [g1, g2] = maxQuaParTest.frozenG1G2('passthrough');
            g = maxQuaPar(g1, g2);
            testCase.verifyClass(g, 'QuaPar');

            testPts = [3.380265 -0.644943; T(1,:); T(2,:); T(3,:); mean(T,1)];
            for i = 1:size(testPts,1)
                s = testPts(i,:);
                testCase.verifyEqual(g.eval(s), maxQuaParTest.supBilinearOverPoly(s, T), ...
                    'AbsTol', 1e-6, sprintf('s=(%.6f,%.6f)', s(1), s(2)));
            end
        end

        function matchHalfEdgesRejectsSameSideRayPairingAndDropsSubsumedPieces(testCase)
            % Regression test for three bugs found while diagnosing the 2-convex-edge
            % (splitTwoConvexEdges) "vertex fan" gap left open by the prior session (see DESIGN.md /
            % session handoff), all confirmed on COAP Appendix A.4.3's own example below:
            %
            % 1) matchHalfEdges used to pair two rays purely by (apex, direction) equality, with no
            %    check that the two candidate pieces are on OPPOSITE sides of that ray. When a g1
            %    (or g2) face is cut into several (k,l) sub-pieces by different opposing faces,
            %    every sub-piece independently inherits the parent face's own boundary rays -- so
            %    two sub-pieces descending from the SAME side could be wrongly paired as if mutually
            %    adjacent. Confirmed via a coverage stress test (this session's scratch scripts, not
            %    committed): two of the pieces this triangle produces geometrically OVERLAPPED over
            %    a positive-area region, not just a shared boundary. Fixed by oppositeSides, which
            %    requires each candidate's own adjacent geometry to fall on opposite sides of the
            %    shared ray before accepting the pairing.
            % 2) Separately, two DIFFERENT (k,l) pairs sharing the same winning row `f` produced
            %    pieces that are not exact geometric duplicates (dedupPieces only collapses those)
            %    yet still overlapped, one wholly containing the other -- pure redundant territory.
            %    Fixed by dropSubsumedPieces (called right after dedupPieces), confirmed to cut this
            %    triangle's piece count from 12 to 10 with zero change in any resolved value.
            % 3) THE remaining "third, undiagnosed pattern" this test used to only pin via
            %    verifyError: facePoly/polyConstraints computed a ray edge's outward half-plane
            %    normal from a FIXED role-based formula (rot90cw(-dirIn) for the "incoming" ray,
            %    rot90cw(dirOut) for the "outgoing" one), implicitly assuming that ray's own sign
            %    in its original P{k} entry is always +1 for "incoming" and always -1 for
            %    "outgoing". That assumption holds for a typical unbounded face (real vertices
            %    between the two rays constrain orderEdges' walk enough to guarantee it), but NOT
            %    for a face whose ENTIRE boundary is just two rays sharing one apex (no real
            %    vertices at all) -- there, orderEdges' own ray-selection rule can legitimately give
            %    BOTH rays the SAME sign. g2's face 3 here is exactly this case (P{3}=[3 2], both
            %    positive): the old formula computed a ~12-degree sliver as "inside" when the true
            %    region (confirmed against g2.eval, i.e. QuaPar's own ground-truth face-membership
            %    test) is the ~170-degree complementary wedge -- so g1 face 3 (the piece that
            %    should tile exactly that wedge from the other side) was left with a genuinely
            %    unmatched boundary edge, not a T-junction sliver like checkOrphanHalfEdges' other
            %    known-safe case, hence the loud crash. Fixed by having facePoly capture each ray's
            %    true sign(Pk(t)) (dirInSign/dirOutSign, swapped alongside dirIn/dirOut through the
            %    CW->CCW reversal, then threaded through clipPolyHalfPlane/splitCell so every piece
            %    -- not just the original g1/g2 faces -- carries a correct sign) and having
            %    polyConstraints use it: outward normal = sign*rot90ccw(direction), which reduces
            %    to EXACTLY the old formula whenever dirInSign==+1/dirOutSign==-1 (the previously
            %    assumed, still-correctly-handled case), so no previously-passing case is affected.
            %
            % All three fixes were independently correct AND, together, closed the then-open gap:
            % this triangle fully ASSEMBLED and matched ground truth (supBilinearOverPoly) at both
            % the sample points spot-checked and, in a separate randomized sweep, at 200 random
            % points to machine precision -- when g1/g2 came from the PRE-2026-07-17-session
            % splitTwoConvexEdges, which returned two plain quadratics (q1 and a bespoke
            % buildEdgeAffinePiece quadratic) for this triangle's split.
            %
            % UPDATE (2026-07-17, later session): that buildEdgeAffinePiece quadratic was found to
            % be a genuine tightness bug (see convEnvCPLQ.m's header and the session handoff) and
            % replaced with the correct Appendix A.3 rational (quadratic/linear) envelope for the
            % "other" sub-triangle. conjPieceCPLQ cannot yet conjugate a rational piece (a
            % pre-existing, separately-tracked gap -- see DESIGN.md), so buildG1G2ForTriangle can
            % no longer build g1/g2 for THIS triangle at all; the three maxQuaPar assembly fixes
            % above are no longer exercisable via this specific repro now that Step 1 refuses to
            % (wrongly) hand it two plain quadratics. They remain covered by this file's OTHER
            % assembly regression tests (e.g. maxQuaParResolvesBothHyperbolaCellsWithoutMisclassifying,
            % assemblePiecesResolvesNearDuplicateApexCluster), which use different repro triangles
            % untouched by the Step 1 fix. This test now just pins the new, correct, LOUD failure
            % mode in place of the old silent-wrong-answer one.
            V = [2 1; 0 0; 1 0];
            testCase.verifyError(@() maxQuaParTest.buildG1G2ForTriangle(V), ...
                'extractTriFace:rationalFaceNotSupported');
        end

        function maxQuaParResolvesBothHyperbolaCellsWithoutMisclassifying(testCase)
            % Calling the public entry point end-to-end. maxQuaPar(g1,g2) now fully ASSEMBLES (the
            % face-clipping topology gap this test used to pin -- a plain decided cell, g1 face 1
            % vs g2 face 4, whose boundary edge had no matching neighbour -- is fixed; see
            % maxQuaPar.m's header HISTORY for the four distinct bugs that were involved: a
            % wraparound-index bug and an endpoint-ordering bug in clipPolyHalfPlane, a ray
            % half-edge orientation bug and a missing collinear-vertex split in assemblePieces/
            % clipByFace) AND matches ground truth at all 7 sample points (see
            % verifyAgainstGroundTruth), including s=(-3,2), which used to be wrong due to two
            % further, separate bugs (also now fixed, see maxQuaPar.m's header HISTORY): a
            % pivot-vertex bug in QuaPar.m's orderEdges, and a ray left/right assignment bug in
            % this file's assemblePieces.
            [g1, g2] = maxQuaParTest.buildG1G2();
            g = maxQuaPar(g1, g2);
            testCase.verifyClass(g, 'QuaPar');
            maxQuaParTest.verifyAgainstGroundTruth(testCase, g);
        end
    end

    methods (Static)
        function [g1, g2] = buildG1G2()
            % q1 over T1=(0,0),(1,2),(2,2) and q2 over T2=(1,2),(2,2),(3,3): the rank-1 PSD convex
            % envelopes of f(x,y)=xy given in closed form by 3-edge.tex eq (eq:q1)/(eq:q2).
            s2v2 = sqrt(2);
            A1 = [12-8*s2v2, 6*s2v2-8; 6*s2v2-8, 6-4*s2v2];
            V1 = [0 0; 1 2; 2 2]; E1 = [1 2 1; 2 3 1; 3 1 1]; F1 = [1 0; 1 0; 1 0];
            p1 = QuaPoly(V1, E1, [A1(1,1) A1(1,2) A1(2,2) 0 0 0], F1);
            g1 = conjPieceCPLQ(p1);

            A2 = [6-4*s2v2, 6*s2v2-8; 6*s2v2-8, 12-8*s2v2];
            b2 = [9-6*s2v2; 6*s2v2-9];
            V2 = [1 2; 2 2; 3 3]; E2 = [1 2 1; 2 3 1; 3 1 1]; F2 = [1 0; 1 0; 1 0];
            p2 = QuaPoly(V2, E2, [A2(1,1) A2(1,2) A2(2,2) b2(1) b2(2) 0], F2);
            g2 = conjPieceCPLQ(p2);
        end

        function [g1, g2] = buildG1G2ForTriangle(T)
            % Generalization of buildG1G2 to an ARBITRARY triangle T (CCW), derived from the real
            % Step-1/Step-2 pipeline (convEnvCPLQ then conjPieceCPLQ) instead of hard-coded
            % closed-form coefficients, so it works for any T with all 3 edges convex for u1*u2
            % (f(x,y)=xy genuinely nonconvex on it), not just the one 3-edge.tex example.
            E = [1 2 1; 2 3 1; 3 1 1]; F = [1 0; 1 0; 1 0];
            q = QuaPoly(T, E, [0 1 0 0 0 0], F);   % f(x,y) = x*y over T (c=[A11 A12 A22 d e f])
            r = convEnvCPLQ(q);
            if r.nf ~= 2
                error('buildG1G2ForTriangle:unexpectedSplit', ...
                    'expected 2 sub-triangles from convEnvCPLQ, got %d -- T does not have all 3 edges convex', r.nf);
            end
            [V1, f1] = maxQuaParTest.extractTriFace(r, 1);
            [V2, f2] = maxQuaParTest.extractTriFace(r, 2);
            p1 = QuaPoly(V1, E, f1(5:10), F);
            p2 = QuaPoly(V2, E, f2(5:10), F);
            g1 = conjPieceCPLQ(p1);
            g2 = conjPieceCPLQ(p2);
        end

        function [g1, g2] = frozenG1G2(key)
            % Hardcoded g1,g2 QuaPar fixtures for the 5 tests below that used to call
            % buildG1G2ForTriangle on a hand-picked 3-convex-edge triangle T. Captured ONCE, from
            % the PRE-Part-2c pipeline (commit f19fe16, before convEnvCPLQ.m's solveTriangleBF/
            % assemblePiecesBF recursion into splitThreeConvex's own sub-triangles), via
            % buildG1G2ForTriangle(T) itself run against that old code.
            %
            % Why freeze rather than call the live pipeline: Step 1 now correctly recurses into
            % each 3-convex-edge sub-triangle (see DESIGN.md's Part 2c -- the tangentCevian
            % tightness criterion depends only on a triangle's own mh,mw, not on how it arose, so
            % there is no basis for exempting a splitThreeConvex sub-triangle from it), which for
            % every T used below now correctly produces 2 rational sub-pieces (alongside 2 plain
            % quadratic ones) that conjPieceCPLQ cannot yet conjugate (a pre-existing, separately-
            % tracked gap -- see conjPieceCPLQ.m's own TODO). These 5 tests exercise maxQuaPar's
            % OWN assembly logic (dedup, orphan half-edges, near-duplicate vertices) given a KNOWN
            % g1,g2 pair -- not Step 1's tightness -- so freezing the fixture keeps that regression
            % coverage intact and decouples it from Step 1's now-stricter (and correct) behaviour.
            % Each fixture was verified (at the point this freeze was made) to reconstruct into a
            % valid QuaPar and to match maxQuaPar(g1,g2)'s ground truth (supBilinearOverPoly over T)
            % exactly at that test's own check points.
            switch key
                case 'splitCell'
                    V1 = [6.9213742009476302 5.4787418770627756;5.8233887614520903 4.6096111765516294;0 0;7.5553364060137298 4.705191663541755;8.1080407608331768 7.6934444929145327;7.0100553213376369 6.8243137924033865;4.8461515396335137 5.1114356958638716;0.63396220506609946 -0.7735502135210206;-0.97723722181857697 0.50182451931224215];
                    E1 = [1 2 1;2 3 1;1 4 0;1 5 0;2 6 0;2 7 0;3 8 0;3 9 0];
                    Ec1 = zeros(8,6);
                    f1 = [0 0 0 0 0.61009174311926606 0.49999999999999989 0.40977443609022535 0 0 0;0 0 0 0 5.0009091382765476 -2.6795523701210318 1.4357391237657371 -12.970543519417404 16.385876046362146 0;0 0 0 0 0.25675675675675669 0.49999999999999994 0.97368421052631582 0 0 0;0 0 0 0 0 0 0 6.9620441895682443 5.7057354636311173 -39.723582431786497;0 0 0 0 0 0 0 3.7999999999999998 7.4000000000000004 -28.119999999999997;0 0 0 0 0 0 0 0 0 0];
                    F1 = [1 2;1 3;4 1;2 4;5 2;3 5;1 6;6 3];
                    V2 = [6.9213742009476285 5.4787418770627756;7.254699338735203 7.2808984503001479;7.6299999999999999 9.3100000000000005;7.7510832593688797 4.4663445855946424;5.3683037084251772 2.5802113220908738;5.7016288462127518 4.3823678953282457;7.0674473143075041 11.766805644198524;8.459709058421252 8.2976027085318673;7.4427479755723009 13.795907193898376];
                    E2 = [1 2 1;2 3 1;1 4 0;1 5 0;2 6 0;2 7 0;3 8 0;3 9 0];
                    Ec2 = zeros(8,6);
                    f2 = [0 0 0 0 0.61009174311926662 0.50000000000000011 0.40977443609022551 -2.8172001158457836e-15 -5.5422868159664301e-16 1.0319738051329496e-14;0 0 0 0 5.0009091382765902 -2.6795523701210531 1.4357391237657482 -12.970543519417546 16.38587604636222 2.5579538487363607e-13;0 0 0 0 11.978260869565259 0.49999999999999906 0.020871143375680429 -86.739130434782908 3.6206896551724204 314.05547226386915;0 0 0 0 0 0 0 6.9620441895682443 5.7057354636311173 -39.72358243178649;0 0 0 0 0 0 0 3.7999999999999998 7.4000000000000004 -28.119999999999997;0 0 0 0 0 0 0 9.3100000000000005 7.6299999999999999 -71.035300000000007];
                    F2 = [2 1;3 1;1 4;4 2;2 5;5 3;6 1;3 6];
                case 'dedup'
                    V1 = [0 0;1.3278591739894494 2.2607112887288547;1.3372819655304191 2.27675381162984;0.53845132786679284 -1.0577577196316552;-0.62128821971814063 0.91672597454914462;0.70657095427130878 3.1774372632779992;12.46277728354762 23.192653396003031;1.875733293397212 1.2189960919981848;12.47220007508859 23.20869591890402];
                    E1 = [1 2 1;2 3 1;1 4 0;1 5 0;2 6 0;2 7 0;3 8 0;3 9 0];
                    Ec1 = zeros(8,6);
                    f1 = [0 0 0 0 0.98222222222222177 0.49999999999999989 0.25452488687782798 0 0 0;0 0 0 0 0.73776223776223759 0.49999999999999989 0.33886255924170605 0 0 0;0 0 0 0 384.6417291168562 -204.6133193606845 108.84573173047006 -46.068407804000543 27.058898780491809 0;0 0 0 0 0 0 0 0 0 0;0 0 0 0 0 0 0 2.1099999999999999 1.4299999999999999 -3.0173000000000001;0 0 0 0 0 0 0 2.4518849697359095 1.2481314891189583 -3.0602748384248728];
                    F1 = [2 1;3 1;1 4;4 2;2 5;5 3;6 1;3 6];
                    V2 = [1.3372819655304191 2.27675381162984;1.3458354768443477 2.2945039872434982;4.5 8.8399999999999999;1.907523156292924 1.1565466724430531;-10.455035432949291 -19.89099835923669;-10.446481921635362 -19.873248183623033;0.80602648297702628 3.4778637685878211;5.0702411907625047 7.7197928608132127;3.9601910061326784 10.023359781344322];
                    E2 = [1 2 1;2 3 1;1 4 0;1 5 0;2 6 0;2 7 0;3 8 0;3 9 0];
                    Ec2 = zeros(8,6);
                    f2 = [0 0 0 0 0.98222222222222189 0.5 0.25452488687782815 2.243678392160306e-16 -8.4003054283321348e-17 -4.8037992065489674e-17;0 0 0 0 384.6417291169351 -204.61331936072656 108.8457317304925 -46.068407804010235 27.058898780496975 5.9507954119908391e-13;0 0 0 0 1.0960912052117258 0.49999999999999989 0.22808320950965835 -0.51241042345276844 0.233744427934621 0.11977308129770434;0 0 0 0 0 0 0 2.4518849697359095 1.2481314891189583 -3.0602748384248728;0 0 0 0 0 0 0 2.1099999999999999 1.4299999999999999 -3.0172999999999996;0 0 0 0 0 0 0 8.8399999999999999 4.5 -39.780000000000001];
                    F2 = [2 1;3 1;1 4;4 2;2 5;5 3;6 1;3 6];
                case 'apexCluster'
                    V1 = [3.3753000000000002 1.4908000000000001;4.9522574593895978 6.0311394177212279;4.9522622221441015 6.0311531305337596;2.8453762270784715 3.0129499884355244;3.9039754045753527 -0.034944261683816258;5.4809328639649504 4.5053951560374115;418.81941580694604 1194.5799883499049;4.4223384492225728 7.5533031189692839;418.81942056970053 1194.5800020627173];
                    E1 = [1 2 1;2 3 1;1 4 0;1 5 0;2 6 0;2 7 0;3 8 0;3 9 0];
                    Ec1 = zeros(8,6);
                    f1 = [0 0 0 0 1.4361971157132118 0.5 0.17407081330604859 -4.1021961246668042 1.4281452315233427 5.8585318342164312;0 0 0 0 1.4429877468097261 0.49999999999999994 0.17325164441120181 -4.1251165418068672 1.4293664485117799 5.8963031810596762;0 0 0 0 875655.3881096642 -304913.85145361535 106174.70990383932 -2497486.9426096831 869662.3676662168 3561571.67010098;0 0 0 0 0 0 0 1.4907999999999999 3.3753000000000002 -5.0318972400000002;0 0 0 0 0 0 0 6.0365000000000002 4.9504000000000001 -29.883089600000005;0 0 0 0 0 0 0 6.025805160298936 4.9541240732007266 -29.852586405054122];
                    F1 = [1 2;1 3;4 1;2 4;5 2;3 5;1 6;6 3];
                    V2 = [4.9522622221441015 6.0311531305337605;4.9522669971108817 6.0311668082677476;6.3015000000000008 9.8960000000000008;4.4234011753527067 7.5502505505635042;-408.08491477095936 -1180.1341450325901;-408.08490999599258 -1180.1341313548562;5.4825916636689191 4.5162618785210622;5.772638953208606 11.415097420029745;6.8318246665580382 8.3810950702533162];
                    E2 = [1 2 1;2 3 1;1 4 0;1 5 0;2 6 0;2 7 0;3 8 0;3 9 0];
                    Ec2 = zeros(8,6);
                    f2 = [0 0 0 0 1.4361971157132138 0.50000000000000011 0.17407081330604848 -4.102196124666813 1.428145231523344 5.8585318342164499;0 0 0 0 875655.39285227843 -304913.85310511192 106174.71047893159 -2497486.9561359221 869662.37237639062 3561571.68938983;0 0 0 0 1.4282806602027993 0.50000000000000011 0.17503562637647349 -4.0523105802679362 1.4185974413784179 5.748597420838788;0 0 0 0 0 0 0 6.025805160298936 4.9541240732007266 -29.852586405054122;0 0 0 0 0 0 0 6.0365000000000002 4.9504000000000001 -29.883089600000005;0 0 0 0 0 0 0 9.8960000000000008 6.3014999999999999 -62.359644000000003];
                    F2 = [1 2;1 3;4 1;2 4;5 2;3 5;1 6;6 3];
                case 'orphan1'
                    V1 = [0.91849999999999998 1.3296000000000001;1.7991661135465991 5.0397665221340056;1.7992476142271701 5.040109877114892;1.428953199752891 -0.84987290884643985;0.40116798140642318 3.4800930684938272;1.2818340949530223 7.1902595906278322;75.618676931316003 320.36464135706939;2.3097008139800614 2.860636968268452;75.618758431996568 320.36498471205027];
                    E1 = [1 2 1;2 3 1;1 4 0;1 5 0;2 6 0;2 7 0;3 8 0;3 9 0];
                    Ec1 = zeros(8,6);
                    f1 = [0 0 0 0 2.1348410685852461 0.50000000000000011 0.11710473612242922 -1.2960515214955481 0.30354754285161817 0.39341325475907479;0 0 0 0 2.0784457478005862 0.50000000000000011 0.12028218694885366 -1.2442524193548388 0.29932280423280427 0.3724331233347416;0 0 0 0 44686.78010489064 -10461.452649719033 2449.0910127206894 -27670.646551283513 6480.8493102777647 8561.0684962831838;0 0 0 0 0 0 0 1.3295999999999999 0.91849999999999998 -1.2212376;0 0 0 0 0 0 0 5.0151000000000003 1.8050999999999999 -9.0527570100000005;0 0 0 0 0 0 0 5.0651111164680822 1.793392087152792 -9.0837301968235025];
                    F1 = [2 1;3 1;1 4;4 2;2 5;5 3;6 1;3 6];
                    V2 = [1.7992476142271694 5.040109877114892;1.7993279955130483 5.0404580115552058;2.6141999999999994 8.5696999999999992;2.3160630633986856 2.8334721855735241;-72.940343918258279 -314.21494998603282;-72.940263536972395 -314.21460185159253;1.289834091504743 7.2788055967584224;3.1310154491715156 6.3630623084586313;2.1047060959916939 10.808047585203216];
                    E2 = [1 2 1;2 3 1;1 4 0;1 5 0;2 6 0;2 7 0;3 8 0;3 9 0];
                    Ec2 = zeros(8,6);
                    f2 = [0 0 0 0 2.1348410685852475 0.50000000000000022 0.11710473612242918 -1.2960515214955515 0.30354754285161883 0.39341325475907662;0 0 0 0 44686.78011677786 -10461.452652501868 2449.09101337216 -27670.646558645756 6480.849312001289 8561.0684985630578;0 0 0 0 2.1966382400197766 0.50000000000000011 0.11381027401114049 -1.4576016870596964 0.33178009480672943 0.48360322592311478;0 0 0 0 0 0 0 5.0651111164680822 1.793392087152792 -9.0837301968235007;0 0 0 0 0 0 0 5.0151000000000003 1.8050999999999999 -9.0527570099999988;0 0 0 0 0 0 0 8.5696999999999992 2.6141999999999999 -22.402909739999998];
                    F2 = [2 1;3 1;1 4;4 2;2 5;5 3;6 1;3 6];
                case 'orphan2'
                    V1 = [1.1161000000000001 1.6018999999999999;2.4597184184576664 4.3621976224761525;2.459791161416752 4.3623470638757791;1.6692327585982112 0.45340985324096361;0.55705430071112916 2.738242742476646;1.9006727191687953 5.4985403649527989;111.12314928289825 229.88259035877093;3.0129239200149631 3.2138569171167428;111.12322202585733 229.88273980017055];
                    E1 = [1 2 1;2 3 1;1 4 0;1 5 0;2 6 0;2 7 0;3 8 0;3 9 0];
                    Ec1 = zeros(8,6);
                    f1 = [0 0 0 0 1.0381686212814651 0.50000000000000022 0.24080866525460212 -0.35774999821224274 0.17229859912865289 0.061639823530247516;0 0 0 0 1.016323660053303 0.50000000000000022 0.24598463051316602 -0.33336883698549069 0.1640072203809595 0.054674896315623563;0 0 0 0 40029.944622205629 -19287.795073364916 9293.5186971443218 -14320.8706652399 6904.846528272069 2552.5045113402507;0 0 0 0 0 0 0 1.6019000000000001 1.1161000000000001 -1.7878805900000001;0 0 0 0 0 0 0 4.3475999999999999 2.4668999999999999 -10.725094440000001;0 0 0 0 0 0 0 4.3771015324140086 2.4526851536662875 -10.735651944641797];
                    F1 = [2 1;3 1;1 4;4 2;2 5;5 3;6 1;3 6];
                    V2 = [2.459791161416752 4.3623470638757782;2.459863167555123 4.3624980351637248;3.9128999999999996 7.4089999999999998;3.016446277547228 3.2065433149908786;-106.89561001706399 -222.59416282554301;-106.89553801092561 -222.59401185425506;1.908599640535575 5.5296060307499237;4.4695551161304756 6.2531962511151002;3.3616364729804515 8.5761079955861987];
                    E2 = [1 2 1;2 3 1;1 4 0;1 5 0;2 6 0;2 7 0;3 8 0;3 9 0];
                    Ec2 = zeros(8,6);
                    f2 = [0 0 0 0 1.038168621281464 0.5 0.24080866525460218 -0.35774999821224118 0.17229859912865217 0.061639823530246968;0 0 0 0 40029.944614899483 -19287.795069846488 9293.5186954499513 -14320.870662616895 6904.8465270089064 2552.5045108694098;0 0 0 0 1.0585753803596125 0.50000000000000011 0.23616645978963882 -0.43759960580912782 0.20669269941856627 0.090448643789189181;0 0 0 0 0 0 0 4.3771015324140086 2.4526851536662875 -10.735651944641795;0 0 0 0 0 0 0 4.3475999999999999 2.4668999999999999 -10.725094440000001;0 0 0 0 0 0 0 7.4089999999999998 3.9129 -28.990676099999995];
                    F2 = [2 1;3 1;1 4;4 2;2 5;5 3;6 1;3 6];
                case 'orphan3'
                    V1 = [0.70220000000000027 0.63929999999999976;4.8625360530425601 5.0357666307315174;4.8626517477549003 5.0358888919949516;1.3882440290507412 -0.090618656023458088;0.011484669826439475 1.3642822231856446;4.1718207228689996 5.7607488539171623;172.15463106567802 183.27785583582906;5.5486957768056406 4.3059702359714933;172.15474676039037 183.27797809709247];
                    E1 = [1 2 1;2 3 1;1 4 0;1 5 0;2 6 0;2 7 0;3 8 0;3 9 0];
                    Ec1 = zeros(8,6);
                    f1 = [0 0 0 0 0.53197653876051154 0.5 0.46994553666312444 -0.053903925517631228 0.050663818411264605 0.002730978694079603;0 0 0 0 0.52480536591208549 0.50000000000000011 0.47636708051853199 -0.048868327943466398 0.046558525424502582 0.0022752372890088093;0 0 0 0 31722.257015954303 -29773.455071107288 27944.374400138408 -4313.4259952939256 4058.0272860169807 269.45620253110803;0 0 0 0 0 0 0 0.63929999999999998 0.70220000000000005 -0.44891646000000007;0 0 0 0 0 0 0 5.0209000000000001 4.8766999999999996 -24.48542303;0 0 0 0 0 0 0 5.0508571664482487 4.8485832002131497 -24.489501203917168];
                    F1 = [2 1;3 1;1 4;4 2;2 5;5 3;6 1;3 6];
                    V2 = [4.8626517477548994 5.0358888919949525;4.862766497443328 5.0360121587191697;7.7776999999999994 8.1672999999999973;5.5490479130246575 4.3055955795575649;-162.51531203039792 -173.29768956226653;-162.51519728070951 -173.29756629554231;4.1829323362229154 5.773354525097627;8.4640961652697566 7.4370066875626097;7.0978658387795868 8.9046423663784537];
                    E2 = [1 2 1;2 3 1;1 4 0;1 5 0;2 6 0;2 7 0;3 8 0;3 9 0];
                    Ec2 = zeros(8,6);
                    f2 = [0 0 0 0 0.53197653876051176 0.5 0.46994553666312439 -0.053903925517631707 0.050663818411265056 0.0027309786940796447;0 0 0 0 31722.257097758426 -29773.45514788543 27944.374472199357 -4313.4260064326645 4058.027296471364 269.45620328945313;0 0 0 0 0.54229576008272973 0.5 0.46100305110602618 -0.134163733195448 0.12369978070175279 0.01659602437440539;0 0 0 0 0 0 0 5.0508571664482487 4.8485832002131497 -24.489501203917168;0 0 0 0 0 0 0 5.0209000000000001 4.8766999999999996 -24.48542303;0 0 0 0 0 0 0 8.1672999999999991 7.7777000000000003 -63.522809209999977];
                    F2 = [2 1;3 1;1 4;4 2;2 5;5 3;6 1;3 6];
                case 'orphan4'
                    V1 = [8.8806172589539738 6.0726039565465406;8.8805347317338441 6.0725636717187719;6.3675999999999986 4.8459000000000003;10.016833062172575 5.5144556605804347;227.07236783144538 113.69348060256797;227.07228530422523 113.6934403177402;7.7371159848282529 6.6271959235796771;7.5038158032185995 4.2877517040338944;5.2241812530944074 5.4005322518609056];
                    E1 = [1 2 1;2 3 1;1 4 0;1 5 0;2 6 0;2 7 0;3 8 0;3 9 0];
                    Ec1 = zeros(8,6);
                    f1 = [0 0 0 0 0.24561720334509465 0.50000000000000022 1.0178440133476629 0.85895789597977557 -1.7485703042814369 1.5019482695382986;0 0 0 0 9131.7208390658197 -18513.751399447738 37534.98348462626 31337.438803697161 -63512.676309823823 53695.771578137807;0 0 0 0 0.24253242889444232 0.50000000000000011 1.0307899901864581 0.87860050577174931 -1.8113052134445562 1.5914136766393929;0 0 0 0 0 0 0 6.0764922493755016 8.8727019077977758 -53.9149043737524;0 0 0 0 0 0 0 6.0686999999999998 8.8885000000000005 -53.941639949999995;0 0 0 0 0 0 0 4.8459000000000003 6.3676000000000004 -30.856752840000002];
                    F1 = [1 2;1 3;4 1;2 4;5 2;3 5;1 6;6 3];
                    V2 = [9.9668999999999972 6.6139999999999981;8.8806989311188911 6.0726446614065894;8.8806172589539738 6.0726039565465406;11.095926971599567 6.0593831054690597;8.8540903152847932 7.176699481709198;7.7678892464036871 6.6353441431157894;-207.93055423826206 -100.86731558727523;10.009644230553544 5.5179870620156022;-207.93063591042699 -100.86735629213527];
                    E2 = [1 2 1;2 3 1;1 4 0;1 5 0;2 6 0;2 7 0;3 8 0;3 9 0];
                    Ec2 = zeros(8,6);
                    f2 = [0 0 0 0 0.24561720334509565 0.50000000000000022 1.0178440133476587 0.85895789597976724 -1.748570304281414 1.5019482695382651;0 0 0 0 0.25282826409495579 0.50000000000000011 0.98881349715752731 0.7870859745919857 -1.5565624701998835 1.2251484888705835;0 0 0 0 9131.7208976216461 -18513.751518163463 37534.983725309809 31337.439004598895 -63512.676717130795 53695.771922779095;0 0 0 0 0 0 0 6.6139999999999999 9.9669000000000008 -65.921076599999992;0 0 0 0 0 0 0 6.0686999999999998 8.8885000000000005 -53.941639949999988;0 0 0 0 0 0 0 6.0764922493755016 8.8727019077977758 -53.914904373752407];
                    F2 = [1 2;1 3;4 1;2 4;5 2;3 5;1 6;6 3];
                case 'orphan5'
                    V1 = [2.7772999999999999 0.59089999999999998;3.3778098698941186 2.194688752587926;3.3779128239894884 2.1949637132979238;3.3051442067605667 -0.8511693977199708;2.2373430355894803 2.0006197138001331;2.837852905483599 3.6044084663880591;73.501167000250575 192.32575019180749;3.9057570307500553 0.75289431557795305;73.501269954345943 192.3260251525175];
                    E1 = [1 2 1;2 3 1;1 4 0;1 5 0;2 6 0;2 7 0;3 8 0;3 9 0];
                    Ec1 = zeros(8,6);
                    f1 = [0 0 0 0 1.365999076290044 0.5 0.18301622917563176 -3.4983392345803392 1.2805057101801187 4.479643366027271;0 0 0 0 1.3054000658544593 0.49999999999999983 0.19151217051330577 -3.3300376028975909 1.2754854584436872 4.2474145385665514;0 0 0 0 23745.659154165034 -8757.7764756524739 3230.0071478133214 -60985.551629286667 22496.628159700958 78312.273654939025;0 0 0 0 0 0 0 0.59089999999999998 2.7772999999999999 -1.6411065699999998;0 0 0 0 0 0 0 2.1766999999999999 3.3847 -7.3674764899999987;0 0 0 0 0 0 0 2.2133684194265588 3.3711761041599919 -7.461654725273184];
                    F1 = [2 1;3 1;1 4;4 2;2 5;5 3;6 1;3 6];
                    V2 = [3.3779128239894884 2.1949637132979229;3.3780142397183028 2.195242876822538;5.1589999999999998 7.0976999999999997;3.9118890367690979 0.73614168646231737;-67.560072638877642 -190.14486461726798;-67.559971223148835 -190.14458545374336;2.8480470842084249 3.6651004951870245;5.6929762127796089 5.6388779731643943;4.6290328444901219 8.5675576183644857];
                    E2 = [1 2 1;2 3 1;1 4 0;1 5 0;2 6 0;2 7 0;3 8 0;3 9 0];
                    Ec2 = zeros(8,6);
                    f2 = [0 0 0 0 1.3659990762900454 0.50000000000000022 0.18301622917563173 -3.498339234580341 1.2805057101801196 4.4796433660272772;0 0 0 0 23745.659151606036 -8757.7764747081146 3230.007147464818 -60985.551622715451 22496.62815727594 78312.273646501926;0 0 0 0 1.3867440680831877 0.50000000000000011 0.18027839869945145 -3.6053626472411655 1.2999380095509045 4.6867479433638604;0 0 0 0 0 0 0 2.2133684194265588 3.3711761041599919 -7.4616547252731831;0 0 0 0 0 0 0 2.1766999999999999 3.3847 -7.367476489999996;0 0 0 0 0 0 0 7.0976999999999997 5.1589999999999998 -36.617034299999993];
                    F2 = [2 1;3 1;1 4;4 2;2 5;5 3;6 1;3 6];
                case 'passthrough'
                    V1 = [0.71300000000000008 0.3892000000000001;1.9441537896452004 2.6974104784584174;1.9441684996906563 2.6974380573696353;1.2780501287534205 -0.67627872288494217;0.14469473507363417 1.4485758789621155;1.3758485247188346 3.7567863574205327;273.81695893777595 514.53897385285154;2.5092186284440769 1.631959334484693;273.81697364782138 514.53900143176281];
                    E1 = [1 2 1;2 3 1;1 4 0;1 5 0;2 6 0;2 7 0;3 8 0;3 9 0];
                    Ec1 = zeros(8,6);
                    f1 = [0 0 0 0 0.94281787461542288 0.49999999999999972 0.26516255867759742 -0.47762914460079664 0.25329873216267895 0.12098285677132666;0 0 0 0 0.9320482708350204 0.49999999999999983 0.2682264511644073 -0.46995041710536978 0.25210626520681262 0.11847744448881857;0 0 0 0 217575.17582441517 -115568.52279054454 61386.063044099647 -111261.1674565822 59101.519560471461 28443.882764459358;0 0 0 0 0 0 0 0.38919999999999999 0.71299999999999997 -0.27749960000000012;0 0 0 0 0 0 0 2.6907999999999999 1.9477 -5.2408711600000002;0 0 0 0 0 0 0 2.7040866968566215 1.9406425591744674 -5.2476657276174663];
                    F1 = [2 1;3 1;1 4;4 2;2 5;5 3;6 1;3 6];
                    V2 = [1.9441684996906556 2.6974380573696353;1.9441831487172756 2.6974657513399483;4.6783999999999999 7.8665000000000003;2.5105381263491182 1.6294712420639106;-270.56351158718707 -510.33937314992914;-270.56349693816043 -510.33934545595883;1.3792704933645643 3.7681869694558268;5.2447696266584618 6.7985331846942758;4.1134873446472886 8.9372212181158783];
                    E2 = [1 2 1;2 3 1;1 4 0;1 5 0;2 6 0;2 7 0;3 8 0;3 9 0];
                    Ec2 = zeros(8,6);
                    f2 = [0 0 0 0 0.94281787461542332 0.49999999999999989 0.26516255867759747 -0.47762914460079658 0.25329873216267884 0.12098285677132656;0 0 0 0 217575.17743250669 -115568.52364472141 61386.063497816445 -111261.16827889373 59101.519997262214 28443.882974707427;0 0 0 0 0.94768740615959302 0.5 0.26380006569159714 -0.50041076097703929 0.26401678323705025 0.13211683941036226;0 0 0 0 0 0 0 2.7040866968566215 1.9406425591744674 -5.2476657276174663;0 0 0 0 0 0 0 2.6907999999999999 1.9477 -5.2408711600000002;0 0 0 0 0 0 0 7.8665000000000003 4.6783999999999999 -36.802633600000007];
                    F2 = [2 1;3 1;1 4;4 2;2 5;5 3;6 1;3 6];
                otherwise
                    error('frozenG1G2:unknownKey', 'unknown fixture key %s', key);
            end
            g1 = QuaPar(V1, E1, Ec1, f1, F1);
            g2 = QuaPar(V2, E2, Ec2, f2, F2);
        end

        function [Vt, frow] = extractTriFace(r, k)
            % r: a 2-face RatPol (as produced by convEnvCPLQ's 3-convex-edge split, or its
            % 2-convex-edge split when a genuine split is needed). Returns face k's 3 vertices in
            % CCW order and its coefficient row -- but only for a PLAIN QUADRATIC face (den =
            % [0 0 1]): QuaPoly (unlike RatPol) has no denominator, and conjPieceCPLQ cannot yet
            % conjugate a genuinely rational (quadratic/linear) piece (see convEnvCPLQ.m's
            % 2-convex-edge fix, 2026 session: the "other" sub-triangle of a genuine split is
            % exactly a 1-convex-edge Appendix A.3 RATIONAL envelope, not a plain quadratic like it
            % used to be before that fix) -- silently dropping a nonzero denominator here would
            % build the WRONG QuaPoly rather than fail loudly, so this errors clearly instead.
            if any(abs(r.den(k,:) - [0 0 1]) > 1e-9)
                error('extractTriFace:rationalFaceNotSupported', ...
                    ['face %d is a genuinely rational envelope (den=[%g %g %g]); ' ...
                     'conjPieceCPLQ does not yet support conjugating a rational piece ' ...
                     '(see DESIGN.md/convEnvCPLQ.m).'], k, r.den(k,1), r.den(k,2), r.den(k,3));
            end
            edgeIdx = find(r.F(:,1)==k | r.F(:,2)==k);
            vids = unique(r.E(edgeIdx,1:2));
            Vt = r.V(vids,:);
            area2 = (Vt(2,1)-Vt(1,1))*(Vt(3,2)-Vt(1,2)) - (Vt(2,2)-Vt(1,2))*(Vt(3,1)-Vt(1,1));
            if area2 < 0, Vt = Vt([1 3 2],:); end
            frow = r.f(k,:);
        end

        function idx = findRow(fmat, expected6, testCase)
            % Locate the row of fmat (1x10, cubic coefficients 1:4 assumed zero for these
            % quadratic conjugates) whose quadratic/linear/constant part (columns 5:10) matches
            % expected6 within tolerance; fails the test if no such row exists.
            for i = 1:size(fmat,1)
                if all(abs(fmat(i,5:10) - expected6) < 1e-6)
                    idx = i; return
                end
            end
            testCase.verifyFail(sprintf('no face row matched the expected coefficients [%s]', ...
                num2str(expected6)));
            idx = -1;
        end

        function [delta, Delta] = conicInvariants(diffRow)
            % delta = b^2-4ac (quadratic-part discriminant, decides hyperbolic/parabolic/elliptic
            % TYPE); Delta = full 3x3 conic determinant (decides IRREDUCIBILITY: Delta==0 means
            % the conic degenerates to a point/line/pair of lines, regardless of the sign of
            % delta). See maxQuaPar.m's splitCell for the same construction used at runtime.
            a = diffRow(5)/2; b = diffRow(6); c = diffRow(7)/2;
            d = diffRow(8); e = diffRow(9); f = diffRow(10);
            delta = b^2 - 4*a*c;
            M = [a, b/2, d/2; b/2, c, e/2; d/2, e/2, f];
            Delta = det(M);
        end

        function verifyAgainstGroundTruth(testCase, g)
            % h(s) should equal sup_{(x,y) in T} [s1 x + s2 y - x y], T=conv{(0,0),(3,3),(1,2)},
            % everywhere, including at points inside both former "hyperbola" cells and at (-3,2),
            % which used to be wrong (see maxQuaParResolvesBothHyperbolaCellsWithoutMisclassifying).
            T = [0 0; 3 3; 1 2];
            testPts = [1.90 2.50; 1.70 2.00; 2.3431 1.9; 2.05 1.95; 0 0; 5 5; -3 2];
            for i = 1:size(testPts,1)
                s = testPts(i,:);
                testCase.verifyEqual(g.eval(s), maxQuaParTest.supBilinearOverPoly(s, T), ...
                    'AbsTol', 1e-8, sprintf('s=(%.4f,%.4f)', s(1), s(2)));
            end
        end

        function h = supBilinearOverPoly(s, T)
            % Exact sup_{(x,y) in T} [s1 x + s2 y - x y]: the Hessian of the objective is
            % indefinite (eigenvalues +-1), so no interior point can be a local max, and the sup is
            % attained on T's boundary -- checked in closed form (quadratic-in-t along each edge).
            s1 = s(1); s2 = s(2);
            best = -inf;
            n = size(T,1);
            for i = 1:n
                va = T(i,:); vb = T(mod(i,n)+1,:);
                dx = vb(1)-va(1); dy = vb(2)-va(2);
                A = -dx*dy;
                B = s1*dx + s2*dy - va(1)*dy - va(2)*dx;
                C = s1*va(1) + s2*va(2) - va(1)*va(2);
                cand = [0 1];
                if abs(A) > 1e-14
                    tstar = -B/(2*A);
                    if tstar > 0 && tstar < 1, cand(end+1) = tstar; end %#ok<AGROW>
                end
                for t = cand
                    val = A*t^2 + B*t + C;
                    if val > best, best = val; end
                end
            end
            h = best;
        end
    end
end
