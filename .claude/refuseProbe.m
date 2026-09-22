function refuseProbe()
% WHICH error refuses a pair, and on what shape of input? Stop after a few -- the identifier and the
% degrees are what decide whether the vertex route holds.
    L = load(fullfile(tempdir, 'branchRowPruning_env.mat'));  env = L.env;
    g = assembleAlgCrv(envelopeMeshAlg(env));
    A = cell(1, numel(g.EcA));
    for j = 1:numel(g.EcA), A{j} = bridgeQuadCut(g.EcA(j).FP, g.EcA(j).num, g.EcA(j).den); end
    C = cell(1, size(g.EcQ,1));
    for i = 1:size(g.EcQ,1)
        c = g.EcQ(i,:);  M = zeros(3,3);
        M(3,1)=c(1); M(2,2)=c(2); M(1,3)=c(3); M(2,1)=c(4); M(1,2)=c(5); M(1,1)=c(6);
        C{i} = ratQ.trimN(M);
    end
    shown = 0;  counts = containers.Map('KeyType','char','ValueType','double');
    pairs = 0;
    for i = 1:numel(C)
        for j = 1:numel(A)
            pairs = pairs + 1;
            try
                R = ratQ.trimN(ratQ.fromBig(ratQ.resultantVarBig(C{i}, A{j}, 2, 2)));
            catch ME
                id = ME.identifier;
                if counts.isKey(id), counts(id) = counts(id)+1; else, counts(id) = 1; end
                if shown < 3
                    shown = shown + 1;
                    fprintf('REFUSED %s | conic %d (deg %d) x alg %d (deg %d)\n', id, ...
                            i, max(size(C{i}))-1, j, max(size(A{j}))-1);
                    fprintf('   message: %s\n', ME.message(1:min(160,end)));
                end
                continue
            end
        end
        if i >= 8, break, end
    end
    k = counts.keys;
    for t = 1:numel(k), fprintf('count %-26s %d of %d pairs tried\n', k{t}, counts(k{t}), pairs); end
    if isempty(k), fprintf('no refusals in %d conic x algebraic pairs\n', pairs); end
end
