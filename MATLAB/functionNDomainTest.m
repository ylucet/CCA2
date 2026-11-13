classdef functionNDomainTest < matlab.unittest.TestCase
% Regression tests for functionNDomain.conjugateOfPiecePoly's entry adapter.
%
% conjugateOfPiecePoly is edge-indexed end to end -- the scatter that re-indexes constraints by
% edge number, the isQuad chord rewrite, getNormalConeVertex and getSubdiffVertexT1/T2/T2Q all
% address d.ineqs BY EDGE. A constraint that bounds no edge therefore has no slot, and the
% adapter at the top of that routine drops such constraints from its local copy of the domain.
%
% These cases are carved out of testcPLQ/testRectBiconj, which exercises the same code but takes
% ~22 minutes because it runs triangulate + maximum + biconjugateF first. Piece 23 of that test's
% biconjugate is reproduced here directly, so the defect is pinned in ~10 s.

    methods (Test)

        function tangentQuadraticDoesNotEvictAnEdge (testCase)
        % A quadratic constraint active at exactly ONE vertex used to destroy a genuine edge.
        %
        % getEdgeNosInf reports 0 for a constraint with no vertex at all, but for one with
        % exactly one vertex it reports THAT VERTEX'S OWN INDEX -- the same slot the real edge
        % leaving that vertex claims. The scatter is last-write-wins, so the intruder evicted the
        % edge and left the evicted constraint's old slot holding a stale duplicate. Here that
        % turned the triangle below into [quad, s1+2*s2-4, (9*s2)/5-s1+5, s1+2*s2-4]: the edge
        % -s1-7*s2-4 gone, s1+2*s2-4 twice. The isQuad branch then chorded d0.vx(1) to d0.vx(2),
        % two vertices the quadratic does not join, solve() returned a complex conjugate pair,
        % and symbolicFunction.gtd (a bare `if (obj1.f>obj2)`) raised
        % 'Conversion to logical from sym is not possible'.
            [pc, d] = functionNDomainTest.piece23;

            % The collision is a precondition of this test, not an incidental detail: if
            % getEdgeNosInf ever stops producing it, this test stops covering the defect.
            d0 = d.removeInfV;
            if d0.nv == size(d0.ineqs,2)
                d0 = d0.poly2order;
            else
                d0 = d0.poly2orderUnbounded;
            end
            edgeNo = d0.getEdgeNosInf(d0.vars);
            testCase.verifyNotEqual(numel(unique(edgeNo)), numel(edgeNo), ...
                'expected two constraints to claim the same edge slot');

            testCase.verifyEqual(size(pc,2), 5);   % 3 vertex pieces + 2 surviving edge pieces
        end

        function tangentQuadraticConjugateIsExact (testCase)
        % The dropped constraint must be redundant, i.e. the conjugate must be unchanged.
        %
        % 196*s2-148*s1-14*s1*s2-s1^2-49*s2^2+684 is -(s1+7*s2)^2-148*s1+196*s2+684, which in
        % u=s1+7*s2 is -u^2-148*u+1232*s2+684. Its u-derivative -2*u-148 is negative throughout
        % the triangle's u-range [-4, 51/19], so it is maximized on the edge u=-4, where it is
        % 1260+1232*s2 and so maximal at s2=-45/44 -- the vertex (139/44,-45/44), where it is
        % exactly 0. It is <= 0 on the whole triangle: redundant, and safe to drop.
        %
        % Both values below are attained at the vertex (36/5,-8/5), which the evicted edge
        % -s1-7*s2-4 bounds -- so they are exactly what the defect corrupted.
            pc = functionNDomainTest.piece23;
            x = sym('x'); y = sym('y');

            [v, n] = functionNDomainTest.evalAt(pc, [3,-1], x, y);
            testCase.verifyEqual(n, 1, 'the pieces must partition, not overlap, at (3,-1)');
            testCase.verifyTrue(isAlways(v == sym(44)/5));

            [v, n] = functionNDomainTest.evalAt(pc, [2,-2], x, y);
            testCase.verifyEqual(n, 1, 'the pieces must partition, not overlap, at (2,-2)');
            testCase.verifyTrue(isAlways(v == sym(16)/5));
        end

        function twoEdgesOnOneVertexPairAreBothKept (testCase)
        % A LENS -- two genuine edges joining the SAME pair of vertices -- is what neither slot
        % convention in conjugateOfPiecePoly can express, and it is the precondition of the test
        % below. getEdgeNosInf numbers an edge by one of its ENDPOINT vertices, so an arc and its
        % chord get the same number, and the scatter that re-indexes constraints by edge number
        % is last-write-wins: one of the two is destroyed outright.
        %
        % Pinned here so that if the numbering ever stops colliding, the next test stops silently
        % covering nothing.
            [~, d] = functionNDomainTest.halfLens;
            d0 = d.removeInfV;
            if d0.nv == size(d0.ineqs,2)
                d0 = d0.poly2order;
            else
                d0 = d0.poly2orderUnbounded;
            end
            edgeNo = d0.getEdgeNosInf(d0.vars);
            nOn = zeros(size(edgeNo));
            for k = 1:size(d0.ineqs,2)
                nOn(k) = d0.vertexOfEdge(k);
            end
            testCase.verifyTrue( ...
                functionNDomain.edgesStillCollide(edgeNo, true(size(edgeNo)), nOn), ...
                'expected two two-vertex constraints to claim the same edge slot');

            % And the geometry does settle which constraint bounds which edge.
            [eIdx, ok] = functionNDomain.edgeIndexList(d0, true);
            testCase.verifyTrue(ok, 'the edge list must be derivable from the geometry');
            testCase.verifyEqual(numel(eIdx), 2, 'a lens has exactly two edges');
            testCase.verifyEqual(numel(unique(eIdx)), 2, 'and they are distinct constraints');
            % One of them has to be the ARC. Reading the edge list off the slot count instead
            % never selected it, which is the whole defect: the conjugate came out as that of the
            % straight CHORD, finite only on a strip.
            isQ = arrayfun(@(k) d0.ineqs(k).isQuad, eIdx);
            testCase.verifyEqual(nnz(isQ), 1, 'exactly one of the lens edges is the arc');
        end

        function halfLensConjugateIsFiniteEverywhereAndExact (testCase)
        % The half-lens {(s1+s2)^2 <= 4*s1, (s1+s2)^2 <= 4*s2, s2 <= s1} carrying s1 -- piece 1 of
        % f* for f = x*y over the unit square given as TWO faces. It is BOUNDED with a curved
        % edge, so its conjugate is finite on all of R^2 and genuinely curved.
        %
        % It used to come back as two cells carrying max(0, x+y-1), the conjugate of the straight
        % segment conv{(0,0),(1,1)}: +inf off a strip, where the truth is finite. Because f** is a
        % MAX over the pieces, its domain is the INTERSECTION of theirs, so this one piece
        % answering on a strip collapsed the whole biconjugate -- it is the entirety of
        % biconjugateTest/biconjugateOverATwoFaceSubdivisionIsTheEnvelope, at ~10 s instead of the
        % 10-40 minutes that test needs to reach the same piece.
        %
        % The three points checked first are exactly the ones a brute-force sup over the lens
        % showed wrong: +inf where the truth is 0, and 0 where the truth is 1/2.
            pc = functionNDomainTest.halfLens;
            x = sym('x'); y = sym('y');
            testCase.verifyEqual(size(pc,2), 3, ...
                '2 vertex cones + the arc; the chord''s cell is a ray and drops out');

            want = {[0 0], sym(0); [-1 sym(1)/2], sym(0); [2 -1], sym(1)/2; ...
                    [sym(9)/10 sym(3)/5], sym(1)/2; [sym(3)/5 sym(3)/5], sym(1)/5; ...
                    [3 3], sym(5)};
            for i = 1:size(want,1)
                pt = want{i,1};
                [v, n] = functionNDomainTest.evalAt(pc, pt, x, y);
                testCase.verifyEqual(n, 1, sprintf( ...
                    'the cells must partition, not overlap, at (%s,%s)', ...
                    char(string(pt(1))), char(string(pt(2)))));
                testCase.verifyTrue(isAlways(v == want{i,2}), sprintf( ...
                    'half-lens conjugate at (%s,%s)', char(string(pt(1))), ...
                    char(string(pt(2)))));
            end

            % FINITE EVERYWHERE is the property the strip broke; check it away from the region
            % the old answer covered, in all four quadrants.
            for pt = [4 -4; -4 4; -4 -4; 6 6]'
                [~, n] = functionNDomainTest.evalAt(pc, pt', x, y);
                testCase.verifyEqual(n, 1, sprintf( ...
                    'a bounded piece''s conjugate is finite everywhere -- at (%g,%g) it is not', ...
                    pt(1), pt(2)));
            end
        end

        function aBoundedPieceWithATangentVertexConjugatesOntoTheWholePlane (testCase)
        % A BOUNDED piece's conjugate is finite EVERYWHERE. This one was finite on part of the
        % plane only, and that hole is what collapsed the parallelogram's biconjugate to nothing
        % (QuaParCPLQ:conj:emptyResult).
        %
        % Piece 9 of f* for x*y over conv{(0,0),(2,0),(2.5,1),(0.5,1)}, verbatim: 3 vertices, 3
        % genuine edges -- a chord, a line and an ARC -- plus a conic touching at one vertex only.
        % TWO defects met here, and both are of this codebase's recurring kind.
        %   1. 4 constraints for 3 vertices, so the COUNT `size(ineqs,2) == nv` called a bounded
        %      region unbounded, endNv came out nv-1, and it was built one edge cell SHORT.
        %      conjugateOfPiecePoly now derives the edge list from the geometry for any bounded
        %      piece the count mislabels, not only for a lens.
        %   2. At (1/4,7/8) the arc is TANGENT to the chord, so the two half-planes of the vertex
        %      cone coincide and the cone is a HALF-PLANE -- which has no finite vertex, which
        %      region.simplifyUnboundedRegion read as empty. See regionTest/aHalfPlaneIsNotEmpty.
        %
        % Measured against a brute-force sup over the piece: 6 of 10 probe points wrong or
        % uncovered before, 2 after -- and 2026-08-17 the SINGULAR-quadratic overlap behind those
        % two was closed as well (functionNDomain.singularEdgeCut). f is singular here (constant
        % along one whole edge), so the elimination in getInterior returned the gradient map's
        % image LINE, which separates nothing; the KKT multiplier does separate them, is affine
        % in s, and needs no inversion. See the header of getInterior.
        % Do not attack the isQuad chord rewrite for any of this -- it was tried twice and is
        % recorded in DECISIONS.md.
            pc = functionNDomainTest.parallelogramPiece9;
            x = sym('x'); y = sym('y');

            % COVERAGE is the property a bounded piece guarantees, and the one that was broken.
            for pt = [0.6 0.95; 2.2 0.9; 1.5 0.6; 3 3; -2 -2; -1 0.5; 0 0]'
                [~, n] = functionNDomainTest.evalAt(pc, pt', x, y);
                testCase.verifyGreaterThanOrEqual(n, 1, sprintf( ...
                    'a bounded piece''s conjugate is finite at (%g,%g)', pt(1), pt(2)));
            end

            % VALUES, against the closed-form sup (each attained at a vertex of the piece).
            want = {[3 3], sym(23)/8; [sym(11)/5 sym(9)/10], sym(67)/80; ...
                    [-2 -2], sym(1)/4; [-1 sym(1)/2], sym(5)/16; [0 0], sym(0)};
            for i = 1:size(want,1)
                pt = want{i,1};
                v = functionNDomainTest.evalAt(pc, pt, x, y);
                testCase.verifyTrue(isAlways(v == want{i,2}), sprintf( ...
                    'conjugate of the parallelogram''s piece 9 at (%s,%s)', ...
                    char(string(pt(1))), char(string(pt(2)))));
            end

            % THE SINGULAR-QUADRATIC OVERLAP, pinned by the property that fails without the fix:
            % the cells must not OVERLAP, and the value must be the sup. At (1/2,1/4) two cells
            % claimed the point and the first one checked answered 0.00195 where the sup over the
            % piece is 0.03845 (measured on a 1400^2 grid of the piece's own bounding box).
            % f = s1*s2/2 + s1^2/8 + s2^2/2 has a SINGULAR Hessian, so getInterior's elimination
            % returns the gradient map's image line and separates nothing.
            [v, nCells] = functionNDomainTest.evalAt(pc, [sym(1)/2, sym(1)/4], x, y);
            testCase.verifyEqual(nCells, 1, ...
                'exactly one cell may claim (1/2,1/4): a singular f made two of them overlap');
            testCase.verifyEqual(double(v), 0.038452, 'AbsTol', 1e-4, ...
                'conjugate of piece 9 at (1/2,1/4), against a brute-force sup over the piece');
        end

    end

    methods (Static)

        function pc = parallelogramPiece9
        % Piece 9 of f* for x*y over the parallelogram conv{(0,0),(2,0),(2.5,1),(0.5,1)}, exactly
        % as the pipeline delivers it. Its three edges are the chord s1 = 2*s2/3 - 1/3, the line
        % s1 + 2*s2 = 0 and the arc; the fourth constraint, s1 + 2*s2 <= 2, touches one vertex and
        % bounds no edge. Bounded, so its conjugate must be finite everywhere.
            s1 = sym('s_1'); s2 = sym('s_2');
            f = symbolicFunction((s1*s2)/2 + s1^2/8 + s2^2/2);
            d = region([ (2*s2)/3 - s1 - sym(1)/3, ...
                         -s1 - 2*s2, ...
                         s1 + 2*s2 - 2, ...
                         16*s1 - 4*s1*s2 - s1^2 - 4*s2^2 ], [s1, s2]);
            pc = functionNDomain(f, d).conjugateOfPiecePoly;
        end

        function [pc, d] = halfLens
        % Piece 1 of f* for f = x*y over the unit square as two faces, verbatim as the pipeline
        % delivers it: nv = 2, five constraints, of which two are the conics (s1+s2)^2 = 4*s1 and
        % (s1+s2)^2 = 4*s2. With s2 <= s1 the second is the ARC and the first meets the region
        % only at its two vertices; -s1-s2 <= 0 and s1+s2-2 <= 0 are implied by the others but no
        % LP can see that past a conic, so they arrive and stay.
            s1 = sym('s_1'); s2 = sym('s_2');
            f = symbolicFunction(s1);
            d = region([ -s1 - s2, ...
                          s1 + s2 - 2, ...
                          2*s1*s2 - 4*s1 + s1^2 + s2^2, ...
                          2*s1*s2 - 4*s2 + s1^2 + s2^2, ...
                          s2 - s1 ], [s1, s2]);
            pc = functionNDomain(f, d).conjugateOfPiecePoly;
        end

        function [pc, d] = piece23
        % Piece 23 of testcPLQ/testRectBiconj's biconjugate, verbatim: a triangle with vertices
        % (139/44,-45/44), (86/19,-5/19) and (36/5,-8/5) cut by three lines, plus a quadratic
        % that touches only (139/44,-45/44).
            s1 = sym('s_1'); s2 = sym('s_2');
            f = symbolicFunction(s1 - 2*s2 + (s1*s2)/2 + s1^2/8 + s2^2/2 + 2);
            d = region([ (9*s2)/5 - s1 + 5, ...
                         -s1 - 7*s2 - 4, ...
                         196*s2 - 148*s1 - 14*s1*s2 - s1^2 - 49*s2^2 + 684, ...
                         s1 + 2*s2 - 4 ], [s1, s2]);
            pc = functionNDomain(f, d).conjugateOfPiecePoly;
        end

        function [val, nHit] = evalAt (pc, pt, x, y)
        % The conjugate's value at pt, plus how many pieces claim pt (must be 1).
            val = sym(NaN);
            nHit = 0;
            for k = 1:size(pc,2)
                r = pc(k).d;
                if isempty(r) || ~r.ptFeasible(r.vars, pt)
                    continue
                end
                nHit = nHit + 1;
                val = subs(pc(k).f.f, [x y], pt);
            end
        end

    end
end
