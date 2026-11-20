function reviewDiffDriver(fixFile, outFile, shard, nshard)
% reviewDiffDriver  2026-11-19 code review: run conjQ / biconjQ on the review fixtures
% (python/tests/review/plq_gen.py output) and write one line per (fixture, op):
%   {"id":..,"op":"conj"|"biconj","err":"" or identifier,"kind":..,"vals":[..]}
% Not a test; called by python/tests/review/run_matlab.sh. Fixtures are taken k = shard,
% shard + nshard, ... so several MATLABs can share the list.
    if nargin < 3, shard = 1; nshard = 1; end
    p = strsplit(path, pathsep); m = p(contains(p, [filesep 'toolbox' filesep 'maple']));
    if ~isempty(m), rmpath(m{:}); addpath(fullfile(matlabroot,'toolbox','symbolic','symbolic')); end
    addpath(fileparts(fileparts(mfilename('fullpath'))));     % MATLAB/
    F = jsondecode(fileread(fixFile));
    if isstruct(F), F = num2cell(F); end
    fid = fopen(outFile, 'w');
    for k = shard:nshard:numel(F)
        x = F{k};
        f = build(x);
        ops = {'conj', 'biconj'};
        for o = 1:2
            if o == 1, P = mat(x.S); else, P = mat(x.X); end
            err = ''; kind = ''; vals = [];
            try
                if o == 1, g = conjQ(f); else, g = biconjQ(f); end
                kind = class(g);
                vals = g.eval(P);
            catch e
                err = e.identifier;
                if isempty(err), err = ['CRASH ' strrep(e.message, '"', '''')]; end
                err = strrep(err, newline, ' ');
            end
            fprintf(fid, '{"id":%d,"op":"%s","err":"%s","kind":"%s","vals":[%s]}\n', ...
                x.id, ops{o}, err, kind, numlist(vals));
        end
        fprintf('%d %s\n', x.id, x.name);
    end
    fclose(fid);
end

function s = numlist(v)
    v = v(:).';
    c = cell(1, numel(v));
    for i = 1:numel(v)
        if isinf(v(i)) && v(i) > 0, c{i} = '"inf"';
        elseif isinf(v(i)), c{i} = '"-inf"';
        elseif isnan(v(i)), c{i} = '"nan"';
        else, c{i} = sprintf('%.17g', v(i)); end
    end
    s = strjoin(c, ',');
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
    V = mat(x.V);
    f = QuaPol(V, E, ff, Fm);
end

function M = mat(A)
    if iscell(A), M = cell2mat(cellfun(@(r) r(:).', A, 'UniformOutput', false)); else, M = A; end
    if isvector(M) && size(M,1) > 1 && size(M,2) == 1, M = M.'; end
end
