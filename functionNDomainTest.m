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

    end

    methods (Static)

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
