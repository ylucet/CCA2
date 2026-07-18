function p = quaPolyToPlq(obj)
% quaPolyToPlq  Convert a CCA2 QuaPoly (possibly multi-face, genuinely nonconvex) into a cPLQ
%   `plq` object, so the integrated cPLQ pipeline (triangulate/convexEnvelope/conjugate/maximum/
%   biconjugateF) can be run on it. Part of the Phase 1 cPLQ integration (DESIGN.md II.5.1,
%   .claude/SESSION_HANDOFF.md).
%
% [input]  obj : QuaPoly, operable (degree<=2), any number of TRIANGULAR faces.
% [output] p   : plq, one plq_1p piece per face (domain = the face's own triangle, function =
%                the face's quadratic re-expressed as a sym expression in x,y).
%
% Each face's quadratic q(x,y) = 1/2 x'Qx + L'x + c is rebuilt as a symbolic expression from
% QuaPoly.matrixForm's own [L,Q,C] decomposition (same convention QuaPoly/RatPol/conjPieceCPLQ
% already use), then wrapped as symbolicFunction(expr) over domain(V3,x,y) (V3 = the face's own
% vertices, in QuaPoly.orderEdges' own clockwise order -- matching domain.m's own expectation,
% see its "Fix order to clockwise" comment).
    obj.assertOperable();
    x = sym('x'); y = sym('y');
    pieces = plq_1p.empty();
    for k = 1:obj.nf
        iVs = faceVertexIndices(obj, k);
        if numel(iVs) ~= 3
            error('quaPolyToPlq:notImplemented', ...
                'quaPolyToPlq currently supports only triangular faces (face %d has %d vertices).', ...
                k, numel(iVs));
        end
        V3 = obj.V(iVs, :);
        [L, Q, C] = QuaPoly.matrixForm(obj.f(k,:));
        if ~isempty(C)
            error('quaPolyToPlq:cubic', 'Face %d is cubic; quaPolyToPlq requires degree<=2.', k);
        end
        expr = 0.5*(Q(1,1)*x^2 + 2*Q(1,2)*x*y + Q(2,2)*y^2) + L(1)*x + L(2)*y + obj.f(k,end);
        d = domain(V3, x, y);
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
