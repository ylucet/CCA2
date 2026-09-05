function biconj_legacy()
% Differential test of the EXACT envelope against the LEGACY one (biconjCPLQ), over the library's
% own fixture corpus. Same discipline as legacy-diff.m for the conjugate: biconjCPLQ is an
% independent implementation, so where both answer they must agree, and where they differ one is
% wrong and the fixture names it.
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
    nBoth = 0; nDis = 0; onlyE = 0; onlyL = 0; neither = 0; worst = 0;
    for i = 1:numel(fx)
        nm = fx{i}{1};  f = fx{i}{2};
        [he, ee] = tryRun(@() biconjQ(f));
        [hl, el] = tryRun(@() f.biconj());
        if isempty(he) && isempty(hl), neither = neither + 1; continue, end
        if isempty(he), onlyL = onlyL + 1;  fprintf('%-22s LEGACY ONLY (%s)\n', nm, ee);  continue, end
        if isempty(hl), onlyE = onlyE + 1;  fprintf('%-22s EXACT ONLY  (%s)\n', nm, el);  continue, end
        nBoth = nBoth + 1;
        X = samplePts(f, 200);
        ve = he.eval(X);  vl = hl.eval(X);
        both = isfinite(ve) & isfinite(vl);
        d = abs(ve(both) - vl(both)) ./ max(1, abs(vl(both)));
        nv = sum(d > 1e-7);  nd = sum(isfinite(ve) ~= isfinite(vl));
        if nv > 0 || nd > 0
            nDis = nDis + 1;
            fprintf('%-22s DISAGREE %d val, %d dom  (worst %.3g)\n', nm, nv, nd, max([0;d]));
        else
            fprintf('%-22s agree (%d pts)\n', nm, sum(both));
        end
        if ~isempty(d), worst = max(worst, max(d)); end
    end
    fprintf(['BICONJ-LEGACY  both %d (disagree %d, worst %.3g) | exact only %d | legacy only %d ' ...
             '| neither %d\n'], nBoth, nDis, worst, onlyE, onlyL, neither);
end

function [g, err] = tryRun(fn)
    g = [];  err = '';
    try, g = fn(); catch e, id = e.identifier; err = id(max(1,find(id==':',1,'last')+1):end); end
end

function X = samplePts(f, n)
    if isempty(f.V), X = randn(n,2)*2; return, end
    lo = min(f.V,[],1) - 1;  hi = max(f.V,[],1) + 1;
    rng(20260905);
    X = lo + rand(n,2) .* (hi - lo);
    X = [X; f.V];
end
