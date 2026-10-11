classdef EnvelopeMin
% EnvelopeMin  A convex envelope held as a MINIMUM of upper bounds, each attained somewhere: the
%   atoms' pairwise envelopes (envelopeAssemble face lists) and their triple tangent planes on the
%   contact triangles. `biconjAtoms` builds it; its header says why the minimum IS co f.
%
% This is NOT a node of the AlgCon lattice: it is not meshed into one subdivision. It answers
% `eval`, which is what every consumer of biconjQ's result used; meshing it is a separate item.

    properties
        pairs   % cell array of envelopeAssemble results
        tris    % struct array .T (3x2 contact triangle), .s (dual point), .c (g(s))
        atoms   % the atoms, for inspection
    end

    methods
        function obj = EnvelopeMin(pairs, tris, atoms)
            obj.pairs = pairs;  obj.tris = tris;  obj.atoms = atoms;
        end

        function val = eval(obj, X)
        % co f at the points X (n x 2): the minimum over the pairwise envelopes and the triangle
        % faces that contain each point; +Inf off the domain.
            n = size(X,1);
            val = inf(n,1);
            for p = 1:numel(obj.pairs)
                val = min(val, envelopeEval(obj.pairs{p}, X));
            end
            for t = 1:numel(obj.tris)
                T = obj.tris(t).T;
                in = inTriangle(X, T);
                if any(in)
                    val(in) = min(val(in), X(in,:)*obj.tris(t).s(:) - obj.tris(t).c);
                end
            end
        end
    end
end

function in = inTriangle(X, T)
    c = zeros(size(X,1), 3);
    for i = 1:3
        a = T(i,:);  b = T(mod(i,3)+1,:);
        c(:,i) = (b(1)-a(1))*(X(:,2)-a(2)) - (b(2)-a(2))*(X(:,1)-a(1));
    end
    tol = 1e-9 * max(1, max(abs(T(:))))^2;
    in = all(c >= -tol, 2) | all(c <= tol, 2);
end
