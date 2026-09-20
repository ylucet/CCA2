function cornerWins()
% Does a corner face ever win only PART of a cell?
%
% This is the question that decides whether an algebraic face can be emitted with the rows the type
% already has. A corner's value is affine with quadratic-irrational coefficients, so the boundary
% between it and a quadratic face is a conic whose coefficients are NOT rational -- one member of a
% conjugate pair of conics, which is not an edge any conic subdivision can hold. If a corner face
% always wins a WHOLE cell or none of it, that boundary never has to be drawn and the arrangement's
% existing curves are enough. If it splits a cell, the type needs a new kind of edge.
%
% Sampled, because a sample can only UNDERSTATE splitting: finding two winners in one cell proves a
% split, while finding one is evidence and is reported as such.
    L = load(fullfile(tempdir, 'branchRowPruning_env.mat'));  env = L.env;
    [C, req, isCorner, val] = faceTable(env);
    cells = {zeros(0,2)};
    for j = 1:size(C,1)
        next = {};
        for k = 1:numel(cells)
            for sd = [1 -1]
                cand = [cells{k}; j, sd];
                if cellInterior(C, cand), next{end+1} = cand; end %#ok<AGROW>
            end
        end
        cells = next;
    end

    rng(7);
    X = [rand(20000,1)*6-3, rand(20000,1)*6-3];
    nSeen = 0;  nCornerWins = 0;  nSplit = 0;  nSplitWithCorner = 0;  nNoSample = 0;
    for k = 1:numel(cells)
        inCell = pointsIn(C, cells{k}, X);
        if ~any(inCell), nNoSample = nNoSample + 1; continue, end
        nSeen = nSeen + 1;
        app = applicable(req, cells{k});
        P = X(inCell,:);
        V = zeros(size(P,1), numel(app));
        for a = 1:numel(app), V(:,a) = val{app(a)}(P); end
        [~, who] = max(V, [], 2);
        u = unique(who);
        if numel(u) > 1
            nSplit = nSplit + 1;
            if any(isCorner(app(u))), nSplitWithCorner = nSplitWithCorner + 1; end
        end
        if any(isCorner(app(u))), nCornerWins = nCornerWins + 1; end
    end
    fprintf(['%d cells, %d sampled (%d without a sample point)\n' ...
             'a corner wins somewhere in %d | cells with MORE THAN ONE winner %d ' ...
             '(of those, a corner among the winners in %d)\n'], ...
            numel(cells), nSeen, nNoSample, nCornerWins, nSplit, nSplitWithCorner);
end

% --------------------------------------------------------------------------------------------------
function tf = pointsIn(C, rows, X)
    tf = true(size(X,1), 1);
    for t = 1:size(rows,1)
        c = C(rows(t,1), 1:6);  br = C(rows(t,1), 7);  sd = rows(t,2);
        if br ~= 0
            v = ratQ.pairLineEval(c, br, X);
        else
            v = c(1)*X(:,1).^2 + c(2)*X(:,1).*X(:,2) + c(3)*X(:,2).^2 + ...
                c(4)*X(:,1) + c(5)*X(:,2) + c(6);
        end
        tf = tf & (sd * v > 1e-9);
    end
end

function app = applicable(req, rows)
    app = [];
    for i = 1:numel(req)
        r = req{i};
        ok = true;
        for t = 1:size(r,1)
            s = rows(rows(:,1) == r(t,1), 2);
            if isempty(s) || s(1) ~= r(t,2), ok = false; break, end
        end
        if ok, app(end+1) = i; end %#ok<AGROW>
    end
end

% --------------------------------------------------------------------------------------------------
function [C, req, isCorner, val] = faceTable(env)
    C = zeros(0,7);  req = {};  isCorner = [];  val = {};
    for i = 1:numel(env.quad)
        [C, r] = rowsOf(C, env.quad(i).con);
        req{end+1} = r;  isCorner(end+1) = false; %#ok<AGROW>
        val{end+1} = quadValue(env.quad(i).num, env.quad(i).den); %#ok<AGROW>
    end
    for k = 1:numel(env.bridge)
        f = env.bridge{k};
        if ~strcmp(f.kind, 'corner') || isempty(f.con), continue, end
        [C, r] = rowsOf(C, f.con);
        req{end+1} = r;  isCorner(end+1) = true; %#ok<AGROW>
        val{end+1} = cornerValue(f); %#ok<AGROW>
    end
    isCorner = logical(isCorner);
end

function h = quadValue(num, den)
    h = @(X) ratQ.evalFace(num, den, X);       % the WEIGHTED basis, not a hand-rolled one
end

function h = cornerValue(f)
% <s, x> - g(s) at the FIXED corner s, in doubles.
    n = f.src.n;  d = f.src.d;  s = f.s;
    N = [n(5) n(6); n(6) n(7)];  m = [n(8); n(9)];  q = n(10);
    gs = (s.'*N*s/2 + m.'*s + q) / d;
    h = @(X) X*s - gs;
end

function [C, r] = rowsOf(C, con)
    r = zeros(0,2);
    for t = 1:size(con,1)
        c = con(t,1:6);
        if all(c(1:5) == 0), continue, end
        br = 0;
        if size(con,2) >= 8, br = con(t,8); end
        key = [c, br];
        hit = find(all(C == key, 2), 1);
        if isempty(hit), C(end+1,:) = key;  hit = size(C,1); end %#ok<AGROW>
        r(end+1,:) = [hit, con(t,7)]; %#ok<AGROW>
    end
end
