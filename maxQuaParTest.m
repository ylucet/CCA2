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
            T = [0 0; 9.31 7.63; 3.80 7.40];
            [g1, g2] = maxQuaParTest.buildG1G2ForTriangle(T);

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
            T = [0 0; 2.11 1.43; 8.84 4.50];
            [g1, g2] = maxQuaParTest.buildG1G2ForTriangle(T);

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
            T = [6.0365 4.9504; 9.8960 6.3015; 1.4908 3.3753];
            [g1, g2] = maxQuaParTest.buildG1G2ForTriangle(T);

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
            Ts = { [8.5697 2.6142; 5.0151 1.8051; 1.3296 0.9185]
                   [7.4090 3.9129; 4.3476 2.4669; 1.6019 1.1161]
                   [8.1673 7.7777; 5.0209 4.8767; 0.6393 0.7022]
                   [6.6140 9.9669; 6.0687 8.8885; 4.8459 6.3676]
                   [2.1767 3.3847; 0.5909 2.7773; 7.0977 5.1590] };
            for ti = 1:numel(Ts)
                T = Ts{ti};
                [g1, g2] = maxQuaParTest.buildG1G2ForTriangle(T);
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
            T = [7.8665 4.6784; 2.6908 1.9477; 0.3892 0.7130];
            [g1, g2] = maxQuaParTest.buildG1G2ForTriangle(T);
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
