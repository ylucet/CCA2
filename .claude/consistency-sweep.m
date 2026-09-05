function consist_sweep()
% How widespread is the overlapping-cells defect? Run checkQuaConConsistent over every fixture the
% exact conjugate accepts, so the quarantined examples(12) is placed in context rather than assumed
% to be unique.
    fx = {};
    for n = {'energy', 'oneNorm', 'oneNormConjugate'}
        try, fx{end+1} = {n{1}, QuaPol.(n{1})()}; catch, end %#ok<AGROW>
    end
    for src = {'examples', 'examples2'}
        try, P = QuaPol.(src{1})(); catch, continue, end
        for k = 1:numel(P)
            if isa(P{k}, 'QuaPol'), fx{end+1} = {sprintf('%s(%d)', src{1}, k), P{k}}; end %#ok<AGROW>
        end
    end
    % plus the shapes the coverage probe exercises
    sq = [0 0; 1 0; 1 1; 0 1];  Esq = [1 2 1; 2 3 1; 3 4 1; 4 1 1];  Fsq = [1 0;1 0;1 0;1 0];
    Hs = {[2 0;0 2], [1 0;0 0], [1 0;0 -1], [-2 0;0 -2], [-1 0;0 0], [0 0;0 0]};
    nm = {'PD','PSDsing','indef','ND','NSDsing','affine'};
    for i = 1:numel(Hs)
        H = Hs{i};
        fx{end+1} = {sprintf('square-%s', nm{i}), ...
            QuaPol(sq, Esq, [0 0 0 0, H(1,1), H(1,2), H(2,2), 1, -1, 2], Fsq)}; %#ok<AGROW>
    end
    fx{end+1} = {'wedge-PD', QuaPol([0 0; 1 0; 0 1], [1 2 0; 1 3 0], ...
        [0 0 0 0 1 0 1 0 0 0], [1 0; 0 1])};

    nbad = 0;  nrun = 0;
    for i = 1:numel(fx)
        try
            g = conjQ(fx{i}{2});
        catch
            continue
        end
        nrun = nrun + 1;
        [ok, rep] = checkQuaConConsistent(g, 3000);
        if ok
            fprintf('%-22s nf=%-4d CONSISTENT\n', fx{i}{1}, g.nf);
        else
            nbad = nbad + 1;
            fprintf('%-22s nf=%-4d OVERLAP  pairs %s  worst gap %.4g  exact-pairs %d\n', ...
                fx{i}{1}, g.nf, mat2str(rep.pairs(1:min(3,size(rep.pairs,1)),:)), ...
                rep.worst, size(rep.exactPairs,1));
        end
    end
    fprintf('CONSISTENCY: %d of %d conjugates have overlapping cells with different values\n', ...
        nbad, nrun);
end
