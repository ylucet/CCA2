function g = assembleAlgCrv(mesh)
% assembleAlgCrv  Turn `envelopeMeshAlg`'s faces into an AlgCrv -- the typed object `biconj` returns
%   for an input whose pieces meet along a CURVED tie.
%
% [input]  mesh : as `envelopeMeshAlg` returns it -- .kind, .num/.den or .FP/.FB, .con, .alg
% [output] g    : AlgCrv
%
% DEDUPLICATION IS BITWISE for the conic edges, which is the whole reason those rows are canonical
% integers -- `assembleQuaConCells`' header carries the measurement behind that. The ALGEBRAIC edges
% are deduplicated on the same principle, by the integer data that defines them: the face's
% polynomial, its branch, and the quadratic's numerator and denominator.
%
% THE VERTEX TABLE IS FILLED BY `algCrvVertices`, and NOT here, because it costs minutes where the
% rest of this costs seconds -- the degree-16 resultant of two quartics, measured at 1.5 to 2.6 s a
% pair. A caller that wants the TOPOLOGY (an edge list, a figure, the paper's tables) asks for it;
% a caller that wants VALUES does not, since `eval` reads FC and the two edge tables and neither
% mentions a vertex. `withVertices` below is the one-line way to ask.
    if nargin ~= 1
        error('assembleAlgCrv:nargin', 'expects (mesh).');
    end
    if ~isstruct(mesh)
        error('assembleAlgCrv:shape', 'mesh must be envelopeMeshAlg''s struct array.');
    end

    EcQ = zeros(0,6);
    EcA = struct('FP', {}, 'FB', {}, 'num', {}, 'den', {});
    FC  = cell(1, numel(mesh));
    fN  = zeros(numel(mesh), 10);
    fD  = ones(numel(mesh), 1);
    FP  = cell(numel(mesh), 1);
    FB  = zeros(numel(mesh), 1);

    for k = 1:numel(mesh)
        f = mesh(k);
        rows = zeros(0,4);
        for t = 1:size(f.con,1)
            c = f.con(t,1:6);
            if all(c(1:5) == 0), continue, end          % a constant is not a curve
            hit = find(all(EcQ == c, 2), 1);
            if isempty(hit), EcQ(end+1,:) = c;  hit = size(EcQ,1); end %#ok<AGROW>
            br = 0;
            if size(f.con,2) >= 8, br = f.con(t,8); end
            rows(end+1,:) = [hit, f.con(t,7), br, 0]; %#ok<AGROW>
        end
        for t = 1:numel(f.alg)
            a = f.alg(t);
            hit = findAlg(EcA, a);
            if hit == 0
                EcA(end+1) = struct('FP', a.FP, 'FB', a.FB, 'num', a.num, 'den', a.den); %#ok<AGROW>
                hit = numel(EcA);
            end
            rows(end+1,:) = [hit, a.side, 0, 1]; %#ok<AGROW>
        end
        FC{k} = rows;

        if strcmp(f.kind, 'bridge')
            FP{k} = f.FP;  FB(k) = f.FB;
        else
            FP{k} = [];    FB(k) = 0;
            fN(k,:) = f.num;  fD(k) = f.den;
        end
    end

    g = AlgCrv(zeros(0,3), zeros(0,3), EcQ, fN, fD, zeros(0,2), FC, FP, FB, EcA);
end

function g = withVertices(g)  %#ok<DEFNU>
% The same object with its vertices named -- separate because it costs minutes, and most callers
% want the function rather than its topology.
    g.Vname = algCrvVertices(g);
    g.nv = size(g.Vname, 1);
end

function hit = findAlg(EcA, a)
% The same edge, by the integer data that defines it.
    hit = 0;
    for j = 1:numel(EcA)
        if EcA(j).FB == a.FB && ratQ.polyEqual(EcA(j).FP, a.FP) && ...
           isequal(EcA(j).num, a.num) && EcA(j).den == a.den
            hit = j;  return
        end
    end
end
