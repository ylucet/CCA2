function branchRowTies()
% WHAT are the zero-margin cells? The exact test only has to decide those, so their SHAPE decides
% whether it is a combinatorial test or a whole arithmetic layer. Classified, not assumed.
    L = load(fullfile(tempdir, 'branchRowPruning_env.mat'));  env = L.env;
    C = curveTable(env);
    cells = {zeros(0,2)};
    kinds = containers.Map('KeyType','char','ValueType','double');
    shown = 0;
    for j = 1:size(C,1)
        next = {};
        for k = 1:numel(cells)
            for sd = [1 -1]
                cand = [cells{k}; j, sd];
                if ~rationalLive(C, cand), continue, end     % the arrangement's own exact test
                H = halfPlanes(C, cand);
                if size(H,1) < 2, next{end+1} = cand; continue, end %#ok<AGROW>
                [m, tight] = margin(H);
                if abs(m) <= 1e-9
                    kk = classify(C, cand, tight);
                    if kinds.isKey(kk), kinds(kk) = kinds(kk) + 1; else, kinds(kk) = 1; end
                    if shown < 3 && contains(kk, 'other')
                        shown = shown + 1;
                        fprintf('  other: rows %s\n', mat2str(cand));
                    end
                elseif m > 0
                    next{end+1} = cand; %#ok<AGROW>
                end
            end
        end
        cells = next;
        if numel(cells) > 40000, fprintf('stopped at curve %d\n', j); break, end
    end
    ks = kinds.keys;
    for i = 1:numel(ks), fprintf('%-28s %6d\n', ks{i}, kinds(ks{i})); end
    fprintf('live at the end: %d\n', numel(cells));
end

% --------------------------------------------------------------------------------------------------
function kk = classify(C, rows, tight)
% What makes this cell tight: the SAME conic twice (two branches of one pair, or a line against its
% own pair), two DIFFERENT curves, or three.
    idx = tight(tight <= size(rows,1));
    if isempty(idx), kk = 'no tight set'; return, end
    cs = rows(idx,1);
    conics = C(cs, 1:6);
    nDistinctCurve = size(unique(cs), 1);
    nDistinctConic = size(unique(conics, 'rows'), 1);
    anyBranch = any(C(cs,7) ~= 0);
    if nDistinctConic == 1 && nDistinctCurve > 1
        kk = 'same conic, two branches';
    elseif numel(idx) == 2 && nDistinctConic == 2
        kk = sprintf('two curves%s', tern(anyBranch, ' (with branch)', ''));
    elseif numel(idx) >= 3
        kk = sprintf('three curves%s', tern(anyBranch, ' (with branch)', ''));
    else
        kk = 'other';
    end
end

function s = tern(c, a, b)
    if c, s = a; else, s = b; end
end

% --------------------------------------------------------------------------------------------------
function tf = rationalLive(C, rows)
% The rational rows only, by the arrangement's own exact tests, so the ties are classified on the
% cells that actually survive -- which is where the branch curves are.
    lin = zeros(0,3);  curved = zeros(0,6);  sides = [];
    for t = 1:size(rows,1)
        c = C(rows(t,1), 1:6);  br = C(rows(t,1), 7);  sd = rows(t,2);
        if br ~= 0, continue, end
        if all(c(1:3) == 0), lin(end+1,:) = sd * c(4:6); %#ok<AGROW>
        else, curved(end+1,:) = c;  sides(end+1) = sd; %#ok<AGROW>
        end
    end
    tf = true;
    if ~isempty(lin) && ~ratQ.feasible2(lin, true), tf = false; return, end
    if ~isempty(curved)
        tf = cellHasInterior([curved; zeros(size(lin,1),3), lin], ...
                             [(1:size(curved,1)).', sides(:); ...
                              (size(curved,1)+(1:size(lin,1))).', ones(size(lin,1),1)]);
    end
end

% --------------------------------------------------------------------------------------------------
function H = halfPlanes(C, rows)
    H = zeros(0,3);
    for t = 1:size(rows,1)
        c = C(rows(t,1), 1:6);  br = C(rows(t,1), 7);  sd = rows(t,2);
        if br ~= 0
            H(end+1,:) = sd * branchLine(c, br); %#ok<AGROW>
        elseif all(c(1:3) == 0)
            H(end+1,:) = sd * c(4:6); %#ok<AGROW>
        end
    end
end

function L = branchLine(c, b)
    X = [0 0; 1 0; 0 1];
    v = ratQ.pairLineEval(c, b, X);
    L = [v(2) - v(1), v(3) - v(1), v(1)];
    n = hypot(L(1), L(2));
    if n > 0, L = L / n; end
end

function [m, tight] = margin(H)
    R = 50;  nH = size(H,1);
    for i = 1:nH
        n = hypot(H(i,1), H(i,2));
        if n > 0, H(i,:) = H(i,:) / n; end
    end
    M = [H(:,1:2), -ones(nH,1); 1 0 0; -1 0 0; 0 1 0; 0 -1 0];
    b = [-H(:,3); -R; -R; -R; -R];
    m = -inf;  tight = [];  k = size(M,1);
    for i = 1:k-2, for j = i+1:k-1, for l = j+1:k
        S = M([i j l], :);
        if abs(det(S)) < 1e-12, continue, end
        z = S \ b([i j l]);
        if all(M*z >= b - 1e-9) && z(3) > m, m = z(3);  tight = [i j l]; end
    end, end, end
    if isinf(m), m = -R;  tight = []; end
end

function C = curveTable(env)
    C = zeros(0,7);
    for i = 1:numel(env.quad)
        con = env.quad(i).con;
        for t = 1:size(con,1), C = add(C, [con(t,1:6), 0]); end
    end
    for k = 1:numel(env.bridge)
        f = env.bridge{k};
        if ~strcmp(f.kind, 'corner') || isempty(f.con), continue, end
        for t = 1:size(f.con,1), C = add(C, [f.con(t,1:6), f.con(t,8)]); end
    end
end
function C = add(C, r)
    if all(r(1:5) == 0), return, end
    if isempty(find(all(C == r, 2), 1)), C(end+1,:) = r; end
end
