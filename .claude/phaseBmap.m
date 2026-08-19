% Phase B, B1: the direct-formula / symbolic map, MEASURED.
%
% One row per input SHAPE: which route biconj takes, what class comes back, and how long it
% takes. Timings are the MINIMUM of NREP repetitions (CLAUDE.md section 3 -- the machine is
% shared, so a single run settles nothing and the mean is worse than the minimum). The ROUTE is
% read off the returned class plus the dispatch in biconjCPLQ.m, both of which are deterministic
% given the shape; only the cost column is contention-sensitive.
%
% The general-POLYGON row is deliberately absent: it is the A.4/A.5 quadrilateral, measured at
% 43 min in SCIP_READINESS.md, and re-running it here would dominate everything else.
%
%   run:  CCA2DIR=<repo> matlab -batch "run('.claude/phaseBmap.m')"

DIR = getenv('CCA2DIR'); if isempty(DIR), DIR = pwd; end
cd(DIR);
NREP = 3;

E4 = [1 2 1; 2 3 1; 3 4 1; 4 1 1];   F4 = [1 0; 1 0; 1 0; 1 0];
E3 = [1 2 1; 2 3 1; 3 1 1];          F3 = [1 0; 1 0; 1 0];
box = [0 0; 1 0; 1 1; 0 1];

% NOT an anonymous helper: MATLAB captures a closure's variables BY VALUE at creation, so an
% `add = @(nm,p) [cases, {{nm,p}}]` sees the EMPTY cases forever and every row overwrites the
% last. Measured the hard way -- the first run of this script printed one row.
cases = {};

% ---- boxes ---------------------------------------------------------------------------------
cases{end+1} = mkrow('bilinear over a BOX (1 face)',        QuaPol(box, E4, [0 1 0 0 0 0], F4));
cases{end+1} = mkrow('bilinear over a BOX (2 triangles)', ...
            QuaPol(box, [1 2 1; 2 3 1; 1 3 1; 3 4 1; 4 1 1], ...
                   [0 1 0 0 0 0; 0 1 0 0 0 0], [1 0; 1 0; 2 1; 2 0; 2 0]));
cases{end+1} = mkrow('diagonal indefinite over a BOX (x^2-y^2)', QuaPol(box, E4, [1 0 -1 0 0 0], F4));
cases{end+1} = mkrow('convex quadratic over a BOX',              QuaPol(box, E4, [1 0 1 0 0 0], F4));
cases{end+1} = mkrow('bilinear over a DIAMOND',                  QuaPol([1 0; 0 1; -1 0; 0 -1], E4, ...
                                                               [0 1 0 0 0 0], F4));

% ---- triangles, one per Step 1 branch ------------------------------------------------------
tri = { 'affine',                   [0 0; 1 0; 0 1], [0 0 0 -1 2 3]; ...
        'convex',                   [0 0; 1 0; 0 1], [1 0 1 0 0 0]; ...
        'concave',                  [0 0; 1 0; 0 1], [-2 0 -2 0 0 0]; ...
        'indefinite 1-convex-edge', [0 0; 2 0; 1 1], [0 1 0 0 0 0]; ...
        'indefinite 2-convex-edge', [2 1; 0 0; 1 0], [0 1 0 0 0 0]; ...
        'indefinite 3-convex-edge', [0 0; 1 1; 3 2], [0 1 0 0 0 0] };
for t = 1:size(tri,1)
    V = tri{t,2};
    if det([V(2,:)-V(1,:); V(3,:)-V(1,:)]) < 0, V = V([1 3 2],:); end
    cases{end+1} = {sprintf('TRIANGLE, %s', tri{t,1}), QuaPol(V, E3, tri{t,3}, F3)};
end

% ---- unbounded -----------------------------------------------------------------------------
cases{end+1} = mkrow('UNBOUNDED piecewise affine (3 wedges)', ...
            QuaPol([0 0; -1 0; 1 1; 0 -1], [1 2 0; 1 3 0; 1 4 0], ...
                   [0 0 0 0 0 0; 0 0 0 0 1 0; 0 0 0 1 0 0], [1 2; 2 3; 3 1]));

fprintf('\n%-42s %-12s %-10s %s\n', 'shape', 'biconj kind', 'min of 3', 'note');
fprintf('%s\n', repmat('-', 1, 92));
for c = 1:numel(cases)
    nm = cases{c}{1}; p = cases{c}{2};
    best = inf; k = ''; note = '';
    for rep = 1:NREP
        try
            t = tic; h = p.biconj('cplq'); e = toc(t);
            best = min(best, e); k = kind(h);
        catch ME
            note = ME.identifier; best = NaN; break
        end
    end
    fprintf('%-42s %-12s %8.2f s  %s\n', nm, k, best, note);
end

function r = mkrow(nm, p)
    r = {nm, p};
end
