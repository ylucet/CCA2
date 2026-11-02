function ex2Overlay(outDir)
% ex2Overlay  Erratum Example 2's pairs subdivided alone, then overlaid -- the run TODO item 1 c'
%   measures (DECISIONS.md 2026-11-02). Caches under outDir (default tempdir): ex2h.mat (biconjAtoms,
%   slow), pair<g>.mat (2-4 min each), overlay.mat (~4 h on 2026-11-02, ~50 s and ~40 faces a tuple).
%   Run from outside the repo (MATLAB's current folder shadows addpath), e.g.
%       cd(tempdir); addpath('D:\AI\CCA2'); addpath('D:\AI\CCA2\.claude'); ex2Overlay
    if nargin < 1, outDir = tempdir; end
    f = erratumExamplesTest.example2();
    hf = fullfile(outDir, 'ex2h.mat');
    if isfile(hf), load(hf, 'h'); else, h = biconjAtoms(f); save(hf, 'h'); end
    V = f.VN ./ f.VD;  k = convhull(V(:,1), V(:,2));  dom = zeros(0,3);
    for i = 1:numel(k)-1
        p = V(k(i),:);  q = V(k(i+1),:);
        n = [-(q(2)-p(2)), q(1)-p(1)];  r = [n, -n*p.'];  dom(end+1,:) = r / gcd(gcd(r(1), r(2)), r(3)); %#ok<AGROW>
    end
    subs = cell(1, numel(h.pairs));  infos = subs;
    for g = 1:numel(h.pairs)
        pf = fullfile(outDir, sprintf('pair%d.mat', g));
        if isfile(pf), load(pf, 'sub', 'info');
        else
            t = tic;  [sub, info] = envelopeSubdivisionAlg(envelopePrune(h.pairs{g}), dom, zeros(0,3));
            fprintf('pair %d: %d faces, %.0f s\n', g, info.faces, toc(t));  save(pf, 'sub', 'info');
        end
        subs{g} = sub;  infos{g} = info;
    end
    % the overlay is local to envelopeSubdivisionAlg, so this call subdivides the pairs again (~8 min)
    % before overlaying them; the pair caches above are for looking at the pairs alone
    t = tic;  [sub, info] = envelopeSubdivisionAlg(cellfun(@(x) envelopePrune(x), h.pairs, 'UniformOutput', false), dom, zeros(0,3));
    fprintf('pairs + overlay: %d faces, %.0f s\n', info.faces, toc(t));  disp(info)
    save(fullfile(outDir, 'overlay.mat'), 'sub', 'info');
end
