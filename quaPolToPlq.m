function p = quaPolToPlq(obj)
% quaPolToPlq  Convert a CCA2 QuaPol (possibly multi-face, genuinely nonconvex) into a cPLQ
%   `plq` object, so the integrated cPLQ pipeline (triangulate/convexEnvelope/conjugate/maximum/
%   biconjugateF) can be run on it. Part of the Phase 1 cPLQ integration (DESIGN.md II.5.1,
%   .claude/SESSION_HANDOFF.md).
%
% [input]  obj : QuaPol, operable (degree<=2), bounded faces (any vertex count -- a
%                non-triangular face is passed through as-is; cPLQ's own plq.triangulate fan-
%                splits it later, same as it already does for a non-triangular plq_1p domain
%                built directly from cPLQ's own test files).
% [output] p   : plq, one plq_1p piece per face (domain = the face's own polygon, function =
%                the face's quadratic re-expressed as a sym expression in x,y).
%
% Each face's quadratic q(x,y) = 1/2 x'Qx + L'x + c is rebuilt as a symbolic expression from
% QuaPol.matrixForm's own [L,Q,C] decomposition (same convention QuaPol/RatPol/conjPieceCPLQ
% already use), then wrapped as symbolicFunction(expr) over domain(V,x,y) (V = the face's own
% vertices, in QuaPol.orderEdges' own clockwise order -- matching domain.m's own expectation,
% see its "Fix order to clockwise" comment).
%
% BOUNDED ONLY, and the check below is deliberately here rather than left to conjCPLQ's own
% isDomBounded gate, because BOTH failure modes on an unbounded face are SILENT.
%
%   (1) In this function. faceVertexIndices takes one vertex per edge, E(j,1), and never consults
%       QuaPol's ray flag E(:,3)==0 -- for which E(j,1) is the base point and E(j,2) is the
%       DIRECTION point, not a second corner. A face bounded by two rays off the same apex
%       (e.g. the second-quadrant cone of the 4-cone fan V=[0 0;-1 0;0 1;1 0;0 -1],
%       E=[1 2 0;1 3 0;1 4 0;1 5 0]) therefore yields V=[(0,0);(0,0)], and domain() -- the
%       BOUNDED constructor, which closes the vertex loop -- turns that into `NaN <= 0` twice.
%       The ray direction is simply thrown away. Building faces from half-planes via
%       domain.domainEdge instead (it takes inequalities, and region then derives the
%       ray/intmax vertices itself) is what this needs, one inequality per edge, segments and
%       rays alike being "the line through E(j,1) and E(j,2)".
%
%   (2) Downstream, and worse. Fixing (1) alone is NOT enough and must not be done alone: cPLQ's
%       plq_1p reads region's infinity markers as ORDINARY FINITE COORDINATES. Measured on the
%       second quadrant built correctly via domainEdge, with f=(x^2+y^2)/2, whose conjugate is
%       min(s1,0)^2/2 + max(s2,0)^2/2:
%         - triangulate then maximum ERRORS (plq_1p.conjugateFunction -> region.getEdgeNos ->
%           symbolicFunction.getLinearCoeffs, 'Index exceeds the number of array elements');
%         - skipping triangulate RUNS and returns garbage -- pieces carrying 2147483647*s_2 and
%           the constant 4611686014132420609, which is exactly intmax^2, i.e. intmax substituted
%           into the closed-form triangle formulas. Max error 1.15e18.
%       So plq_1p needs a genuine unbounded-piece case (a fan emitting WEDGES as well as
%       triangles, plus an envelope/conjugate branch for a quadratic over a wedge) before any
%       unbounded face may be let through. Erroring here keeps that ordering enforced even if
%       conjCPLQ's gate is relaxed first.
    obj.assertOperable();
    if ~obj.isDomBounded
        error('quaPolToPlq:unboundedFace', ...
            ['quaPolToPlq requires a bounded domain: this QuaPol has a ray edge (E(:,3)==0). ' ...
             'Converting it would silently discard the ray direction, and cPLQ''s plq_1p has ' ...
             'no unbounded-piece case behind that. See this file''s header.']);
    end
    x = sym('x'); y = sym('y');
    pieces = plq_1p.empty();
    for k = 1:obj.nf
        iVs = faceVertexIndices(obj, k);
        V = obj.V(iVs, :);
        [L, Q, C] = QuaPol.matrixForm(obj.f(k,:));
        if ~isempty(C)
            error('quaPolToPlq:cubic', 'Face %d is cubic; quaPolToPlq requires degree<=2.', k);
        end
        expr = 0.5*(Q(1,1)*x^2 + 2*Q(1,2)*x*y + Q(2,2)*y^2) + L(1)*x + L(2)*y + obj.f(k,end);
        d = domain(V, x, y);
        f = symbolicFunction(expr);
        pieces(k) = plq_1p(d, f);
    end
    p = plq(pieces);
end

function iVs = faceVertexIndices(obj, k)
% Vertex indices around face k, in the order of its ordered edge list obj.P{k} (same convention
% as convEnvCPLQ.m's/conjCPLQ.m's own file-local helper of the same name, duplicated here since
% it is file-local there too).
    face = obj.P{k}; iVs = zeros(1, numel(face));
    for i = 1:numel(face)
        j = face(i);
        if j > 0, iVs(i) = obj.E(j,1); else, iVs(i) = obj.E(-j,2); end
    end
end
