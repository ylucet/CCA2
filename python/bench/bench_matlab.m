function bench_matlab(fixFile, outFile, reps, slowCap)
% bench_matlab  Time conjQ / biconjQ on the benchmark fixtures (python/bench/fixtures.json).
%   Each fixture: the FIRST call, then `reps` more; a fixture whose first call takes longer than
%   slowCap seconds is not repeated. Times are wall-clock (tic/toc) inside one MATLAB process, so
%   start-up is excluded. Results: one JSON line per fixture in outFile.
    if nargin < 3, reps = 3; end
    if nargin < 4, slowCap = 20; end
    p = strsplit(path, pathsep); m = p(contains(p, [filesep 'toolbox' filesep 'maple']));
    if ~isempty(m), rmpath(m{:}); addpath(fullfile(matlabroot,'toolbox','symbolic','symbolic')); end
    F = jsondecode(fileread(fixFile));
    if isstruct(F), F = num2cell(F); end
    fid = fopen(outFile, 'w');
    % warm the interpreter (JIT, class loading) on one cheap call that is not a fixture
    conjQ(QuaPol([1 0 1 0 0 0]));  biconjQ(QuaPol([1 0 1 0 0 0]));
    for k = 1:numel(F)
        x = F{k};
        f = build(x);
        if strcmp(x.fn, 'conjQ'), call = @() conjQ(f); else, call = @() biconjQ(f); end
        t = zeros(1, 0);  err = '';
        try
            t0 = tic;  call();  t(end+1) = toc(t0); %#ok<AGROW>
            if t(1) <= slowCap
                for r = 1:reps
                    t0 = tic;  call();  t(end+1) = toc(t0); %#ok<AGROW>
                end
            end
        catch e
            err = e.identifier;
        end
        fprintf(fid, '{"id":%d,"first":%.6f,"best":%.6f,"worst":%.6f,"n":%d,"err":"%s"}\n', ...
                x.id, firstOr(t), min([t, inf]), max([t, -inf]), numel(t), err);
        fprintf('%4d %-8s %-11s first %9.4f  best %9.4f  %s\n', x.id, x.fn, x.kind, firstOr(t), min([t inf]), err);
    end
    fclose(fid);
end

function v = firstOr(t)
    if isempty(t), v = -1; else, v = t(1); end
end

function f = build(x)
    ff = mat(x.f);
    if isempty(x.V)
        f = QuaPol(ff);
        return
    end
    E = mat(x.E);  Fm = mat(x.F);
    if isempty(E), E = zeros(0,3); end
    if isempty(Fm), Fm = zeros(0,2); end
    f = QuaPol(mat(x.V), E, ff, Fm);
end

function M = mat(A)
% jsondecode gives a column for a single-row array of scalars written as [[...]]? No: [[a,b]]
% decodes to 1x2 and [[a,b],[c,d]] to 2x2; a cell appears only for ragged input.
    if iscell(A), M = cell2mat(cellfun(@(r) r(:).', A, 'UniformOutput', false)); else, M = A; end
    if isvector(M) && size(M,1) > 1 && size(M,2) == 1, M = M.'; end
end
