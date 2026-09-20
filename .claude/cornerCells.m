function cornerCells()
% How much of the mesh can an ALGEBRAIC face be emitted into without a new curve type?
%
% A corner face's value is affine with quadratic-irrational coefficients, so where it has to be
% COMPARED with a quadratic face the boundary is a conic with coefficients in Q(sigma) -- one member
% of a conjugate PAIR of conics, which the arrangement has no row for (it names one line of a pair,
% not one conic of one). Where a corner face is the ONLY candidate on a cell, no such curve is
% needed and the face can be emitted as it stands.
%
% So the question is a COUNT: of the cells of the real arrangement, how many carry a corner face at
% all, and of those, how many carry something else as well.
    L = load(fullfile(tempdir, 'branchRowPruning_env.mat'));  env = L.env;
    [C, req, isCorner] = curveTable(env);
    fprintf('%d curves, %d faces (%d of them corners)\n', size(C,1), numel(req), sum(isCorner));

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

    nNone = 0;  nQuadOnly = 0;  nCornerAlone = 0;  nMixed = 0;  nCornerMulti = 0;
    for k = 1:numel(cells)
        app = applicable(req, cells{k});
        c = isCorner(app);
        if isempty(app),        nNone = nNone + 1;
        elseif ~any(c),         nQuadOnly = nQuadOnly + 1;
        elseif numel(app) == 1, nCornerAlone = nCornerAlone + 1;
        elseif all(c),          nCornerMulti = nCornerMulti + 1;
        else,                   nMixed = nMixed + 1;
        end
    end
    fprintf(['%d cells | no face %d | quadratic only %d | ONE corner alone %d | ' ...
             'corners only, several %d | MIXED %d\n'], ...
            numel(cells), nNone, nQuadOnly, nCornerAlone, nCornerMulti, nMixed);
end

% --------------------------------------------------------------------------------------------------
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
function [C, req, isCorner] = curveTable(env)
    C = zeros(0,7);  req = {};  isCorner = [];
    for i = 1:numel(env.quad)
        [C, r] = rowsOf(C, env.quad(i).con);
        req{end+1} = r;  isCorner(end+1) = false; %#ok<AGROW>
    end
    for k = 1:numel(env.bridge)
        f = env.bridge{k};
        if ~strcmp(f.kind, 'corner') || isempty(f.con), continue, end
        [C, r] = rowsOf(C, f.con);
        req{end+1} = r;  isCorner(end+1) = true; %#ok<AGROW>
    end
    isCorner = logical(isCorner);
end

function [C, r] = rowsOf(C, con)
    r = zeros(0,2);
    for t = 1:size(con,1)
        c = con(t,1:6);
        if all(c(1:5) == 0), continue, end           % a constant is not a curve
        br = 0;
        if size(con,2) >= 8, br = con(t,8); end
        key = [c, br];
        hit = find(all(C == key, 2), 1);
        if isempty(hit), C(end+1,:) = key;  hit = size(C,1); end %#ok<AGROW>
        r(end+1,:) = [hit, con(t,7)]; %#ok<AGROW>
    end
end
