function vertexCost()
% What does a VERTEX of an AlgCrv cost? A vertex is where two edges meet, and this object's edges are
% conics (degree 2) and algebraic curves (degree 4), so the three kinds of pair need resultants of
% three different sizes. Measured before any table is built: the resultant's DEGREE, the size of its
% coefficients, and the time -- which together say whether the exact route is affordable at all.
    L = load(fullfile(tempdir, 'branchRowPruning_env.mat'));  env = L.env;
    g = assembleAlgCrv(envelopeMeshAlg(env));
    fprintf('AlgCrv: %d conic edges, %d algebraic edges\n', size(g.EcQ,1), numel(g.EcA));

    % the algebraic edges as explicit curves in (x,y)
    A = cell(1, numel(g.EcA));
    for j = 1:numel(g.EcA)
        A{j} = bridgeQuadCut(g.EcA(j).FP, g.EcA(j).num, g.EcA(j).den);
    end
    fprintf('algebraic curve degrees: %s\n', mat2str(cellfun(@(c) max(size(c))-1, A)));

    kinds = {'conic  x conic ', 'conic  x alg   ', 'alg    x alg   '};
    for kind = 1:3
        for trial = 1:2
            [P, Q] = pickPair(g, A, kind, trial);
            if isempty(P), continue, end
            t0 = tic;
            try
                Rb = ratQ.resultantVarBig(P, Q, 2, 2);        % eliminate y
                R  = ratQ.trimN(ratQ.fromBig(Rb));
                fits = true;
            catch ME
                fprintf('%s trial %d: %s after %.1f s\n', kinds{kind}, trial, ME.identifier, toc(t0));
                continue
            end
            el = toc(t0);
            fprintf('%s trial %d: resultant degree %2d | max |coef| %.3e | %s | %.1f s\n', ...
                    kinds{kind}, trial, numel(R)-1, max(abs(R(:))), ...
                    tern(max(abs(R(:))) < 2^53, 'inside 2^53', 'OVER 2^53 (big only)'), el);
        end
    end
end

function [P, Q] = pickPair(g, A, kind, trial)
    P = [];  Q = [];
    cs = {};
    for i = 1:size(g.EcQ,1), cs{end+1} = conicPoly(g.EcQ(i,:)); end %#ok<AGROW>
    switch kind
        case 1
            if numel(cs) >= trial+1, P = cs{trial};  Q = cs{trial+1}; end
        case 2
            if numel(cs) >= trial && numel(A) >= trial, P = cs{trial};  Q = A{trial}; end
        case 3
            if numel(A) >= trial+1, P = A{trial};  Q = A{trial+1}; end
    end
end

function M = conicPoly(c)
% [a b c d e f] as a dense polyN in (x,y): M(i+1,j+1) is the coefficient of x^i y^j.
    M = zeros(3,3);
    M(3,1) = c(1);  M(2,2) = c(2);  M(1,3) = c(3);
    M(2,1) = c(4);  M(1,2) = c(5);  M(1,1) = c(6);
end

function s = tern(c, a, b)
    if c, s = a; else, s = b; end
end
