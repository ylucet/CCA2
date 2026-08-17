% Step 3's cost, fold by fold, on the case that measures it: f = x*y over the general convex
% quadrilateral conv{(0,0),(2,0),(2.5,1.5),(0.5,1)}, which the A.4/A.5 split turns into 6 pieces.
%
% WHAT THIS IS FOR. Steps 1 and 2 take ~25 s for all six pieces; functionNDomain.maxOfList then
% takes 73 minutes, with the cell count running 5, 14, 29, 45, 70, 86 -- roughly ten times what
% the answer needs (f* of x*y over a convex quadrilateral has a cone per vertex and a cell per
% edge). This script reports, per fold:
%
%   * the cell count after mtimes (the pairing), after maximumP's split loop, and after each of
%     the two mergeL passes -- so growth can be attributed to cells being MADE or merges REFUSED;
%   * the number of DISTINCT functions among the surviving cells, which is the count the answer
%     actually needs: any excess over it is cells that ought to have merged;
%   * region.mergeTally's refusal reasons for that fold.
%
% Run:  CCA2_STEP3_FOLDS=3 matlab -batch "run('.claude/step3cost.m')"
% CCA2_STEP3_FOLDS bounds the number of folds (default: all of them).

DIR = getenv('CCA2DIR');
if isempty(DIR), DIR = pwd; end
cd(DIR);
fprintf('=== repo: %s ===\n', DIR);

setenv('CCA2_A45_SPLIT', '1');
setenv('CCA2_TRACE_MAXP', '1');

maxFolds = str2double(getenv('CCA2_STEP3_FOLDS'));
if isnan(maxFolds), maxFolds = inf; end

V = [0 0; 2 0; 2.5 1.5; 0.5 1];
E = [1 2 1; 2 3 1; 3 4 1; 4 1 1]; F = [1 0; 1 0; 1 0; 1 0];
q = QuaPol(V, E, [0 1 0 0 0 0], F);

tAll = tic;
p = quaPolToPlq(q);
p = p.triangulate;
fprintf('pieces: %d\n', p.nPieces);
for i = 1:p.nPieces
    tt = tic;
    p.pieces(i) = p.pieces(i).convexEnvelope;
    p.pieces(i) = p.pieces(i).conjugate;
    p.pieces(i) = p.pieces(i).maximumConjugate;
    fprintf('piece %d: %d cells (%.0f s)\n', i, size(p.pieces(i).maxConjugate,2), toc(tt));
end
fprintf('Steps 1+2 total %.0f s\n', toc(tAll));

acc = p.pieces(1).maxConjugate;
fprintf('FOLD 0 (seed): cells=%d distinctF=%d\n', size(acc,2), nDistinct(acc));

for j = 2:min(p.nPieces, 1 + maxFolds)
    region.mergeTally('reset');
    t1 = tic;  acc = acc * p.pieces(j).maxConjugate;              tMul = toc(t1);
    nPaired = size(acc,2);
    t2 = tic;  acc = acc.maximumP(true);                          tMax = toc(t2);
    fprintf('FOLD %d: paired=%3d -> cells=%3d distinctF=%3d   mtimes %.0f s, maximumP %.0f s\n', ...
            j-1, nPaired, size(acc,2), nDistinct(acc), tMul, tMax);
    printTally(region.mergeTally('get'));
end
fprintf('TOTAL %.0f s\n', toc(tAll));

function n = nDistinct(acc)
% How many DIFFERENT functions the surviving cells carry. mergeL groups by exactly this test, so
% it is also the number of cells a perfect merge would leave.
    n = 0; reps = {};
    for i = 1:numel(acc)
        f = acc(i).f(1).f;
        isNew = true;
        for k = 1:numel(reps)
            if isAlways(simplifyFraction(f - reps{k}) == 0, 'Unknown', 'false')
                isNew = false; break
            end
        end
        if isNew, n = n + 1; reps{end+1} = f; end %#ok<AGROW>
    end
end

function printTally(T)
    if isempty(fieldnames(T)), fprintf('  merge: no attempts\n'); return, end
    k = fieldnames(T);
    s = '';
    for i = 1:numel(k)
        s = [s sprintf('%s=%d ', k{i}, T.(k{i}))]; %#ok<AGROW>
    end
    fprintf('  merge: %s\n', strtrim(s));
end
