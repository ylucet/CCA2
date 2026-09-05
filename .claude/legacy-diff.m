function legacy_diff()
% Differential test of the EXACT conjugate against the LEGACY one, over the library's own fixture
% corpus (QuaPol.examples / examples2) plus the named factories.
%
% This is the check the exact work has not had at scale: conjCPLQ is an INDEPENDENT implementation
% -- different algorithm, different arithmetic, written years earlier -- so where both answer, they
% must agree. Where they differ, one of them is wrong and the disagreement names the input.
%
% Only inputs both routes accept are compared; the rest are counted by reason, which doubles as a
% measurement of how far the exact route now reaches relative to the legacy one.
    fx = collectFixtures();
    fprintf('%-28s %-22s %-22s %s\n', 'fixture', 'exact', 'legacy', 'agreement');
    fprintf('%s\n', repmat('-', 1, 96));
    nBoth = 0; nDisagree = 0; worst = 0; onlyExact = 0; onlyLegacy = 0; neither = 0;
    for i = 1:numel(fx)
        name = fx{i}{1};  f = fx{i}{2};
        [ge, ee] = tryRun(@() conjQ(f));
        [gl, el] = tryRun(@() f.conj());
        if isempty(ge) && isempty(gl)
            neither = neither + 1;
            fprintf('%-28s %-22s %-22s -\n', name, shorten(ee), shorten(el));
            continue
        elseif isempty(ge)
            onlyLegacy = onlyLegacy + 1;
            fprintf('%-28s %-22s %-22s LEGACY ONLY\n', name, shorten(ee), 'ok');
            continue
        elseif isempty(gl)
            onlyExact = onlyExact + 1;
            fprintf('%-28s %-22s %-22s EXACT ONLY\n', name, 'ok', shorten(el));
            continue
        end
        nBoth = nBoth + 1;
        rng(20260904);
        S = [randn(200,2)*2; 0 0; 1 1; -1 -1];
        ve = ge.eval(S);  vl = gl.eval(S);
        both = isfinite(ve) & isfinite(vl);
        d = abs(ve(both) - vl(both)) ./ max(1, abs(vl(both)));
        mism = sum(d > 1e-7);
        infDiff = sum(isfinite(ve) ~= isfinite(vl));
        if mism > 0 || infDiff > 0
            nDisagree = nDisagree + 1;
            fprintf('%-28s %-22s %-22s DISAGREE %d val, %d dom\n', name, 'ok', 'ok', mism, infDiff);
        else
            fprintf('%-28s %-22s %-22s agree (%d pts)\n', name, 'ok', 'ok', sum(both));
        end
        if ~isempty(d), worst = max(worst, max(d)); end
    end
    fprintf('%s\n', repmat('-', 1, 96));
    fprintf(['LEGACY-DIFF  both %d (disagree %d, worst rel %.2e) | exact only %d | ' ...
             'legacy only %d | neither %d\n'], nBoth, nDisagree, worst, onlyExact, onlyLegacy, neither);
end

function [g, err] = tryRun(fn)
    g = [];  err = '';
    try
        g = fn();
    catch e
        id = e.identifier;
        err = id(max(1, find(id == ':', 1, 'last')+1):end);
    end
end

function s = shorten(x)
    if isempty(x), s = 'ok'; else, s = x; end
end

function fx = collectFixtures()
    fx = {};
    named = {'energy', 'oneNorm', 'oneNormConjugate'};
    for i = 1:numel(named)
        try
            fx{end+1} = {named{i}, QuaPol.(named{i})()}; %#ok<AGROW>
        catch
        end
    end
    for src = {'examples', 'examples2'}
        try
            P = QuaPol.(src{1})();
        catch
            continue
        end
        % examples/examples2 return a CELL array, not an object array
        for k = 1:numel(P)
            e = P{k};
            if ~isa(e, 'QuaPol'), continue, end
            fx{end+1} = {sprintf('%s(%d)', src{1}, k), e}; %#ok<AGROW>
        end
    end
end
