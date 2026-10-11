classdef EnvelopeMin
% EnvelopeMin  A convex envelope held as a MINIMUM of upper bounds, each attained somewhere: the
%   atoms' pairwise envelopes (envelopeAssemble face lists) and their triple tangent planes on the
%   contact triangles. `biconjAtoms` builds it; its header says why the minimum IS co f.
%
% This is NOT a node of the AlgCon lattice: it is not meshed into one subdivision. It answers
% `eval`, which is what every consumer of biconjQ's result used; meshing it is a separate item.

    properties
        pairs   % cell array of envelopeAssemble results
        tris    % struct array .T (3x2 contact triangle), .s (dual point), .c (g(s)), and
                % optionally .G (2 x k): generators of the normal cone of dom f* at s. The face
                % <s,x> - c holds on T + cone(G) -- erratum Case E's T + N_F -- since moving along
                % an outward normal of D at s costs exactly <s, r> (empty G: T alone)
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
                if isfield(obj.tris, 'G') && ~isempty(obj.tris(t).G)
                    in = inTrianglePlusCone(X, T, obj.tris(t).G);   % T + N_D(s), erratum Case E
                else
                    in = inTriangle(X, T);
                end
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

function in = inTrianglePlusCone(X, T, G)
% x in T + cone(G): is there nu >= 0 with x - G nu in T? Over (nu, lambda) the set is a polygon
% in the nonnegative orthant of nu, so when it is nonempty it has a vertex, and a vertex is where
% two of its lines meet: test every such point.
    k = size(G,2);
    in = false(size(X,1),1);
    tol = 1e-9 * max(1, max(abs(T(:))))^2;
    for i = 1:size(X,1)
        x = X(i,:).';
        % unknowns nu (k of them); constraints: nu >= 0 and the three edge rows of T at x - G nu
        A = zeros(0,k);  b = zeros(0,1);              % rows A nu <= b
        for e = 1:3
            a = T(e,:);  c = T(mod(e,3)+1,:);
            % (c-a) x (y - a) >= 0 on one orientation; use the triangle's own orientation
            nrm = [-(c(2)-a(2)), c(1)-a(1)];          % left normal of edge a -> c
            sg = sign(nrm * (T(mod(e+1,3)+1,:) - a).');
            if sg == 0, sg = 1; end
            % sg * nrm . (x - G nu - a) >= 0   <=>   sg*nrm*G nu <= sg*nrm*(x - a)
            A(end+1,:) = sg*nrm*G;  b(end+1,1) = sg*nrm*(x - a.'); %#ok<AGROW>
        end
        A = [A; -eye(k)];  b = [b; zeros(k,1)]; %#ok<AGROW>
        if k == 1
            cand = num2cell(b(A ~= 0) ./ A(A ~= 0)).';
            cand{end+1} = 0;
        else
            cand = {};
            for p = 1:size(A,1)
                for q = p+1:size(A,1)
                    M = A([p q],:);
                    if abs(det(M)) < 1e-14, continue, end
                    cand{end+1} = M \ b([p q]); %#ok<AGROW>
                end
            end
        end
        for c = 1:numel(cand)
            if all(A*cand{c} <= b + tol), in(i) = true; break, end
        end
    end
end
