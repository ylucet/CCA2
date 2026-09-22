function ovProbe()
% WHERE does the overflow come from? One failing pair, with the stack.
    L = load(fullfile(tempdir, 'branchRowPruning_env.mat'));  env = L.env;
    g = assembleAlgCrv(envelopeMeshAlg(env));
    A = cell(1, numel(g.EcA));
    for j = 1:numel(g.EcA), A{j} = bridgeQuadCut(g.EcA(j).FP, g.EcA(j).num, g.EcA(j).den); end
    for i = 1:numel(A)-1
        for j = i+1:numel(A)
            try
                Rb = ratQ.resultantVarBig(A{i}, A{j}, 2, 2);
            catch ME
                fprintf('PAIR %d,%d -> %s: %s\n', i, j, ME.identifier, ME.message(1:min(120,end)));
                for s = 1:min(4, numel(ME.stack))
                    fprintf('   at %s line %d\n', ME.stack(s).name, ME.stack(s).line);
                end
                fprintf('   max |coef| of the inputs: %.3e and %.3e\n', ...
                        max(abs(A{i}(:))), max(abs(A{j}(:))));
                return
            end
        end
    end
    fprintf('no overflow among %d algebraic pairs\n', numel(A)*(numel(A)-1)/2);
end
