function [fN, fD, FC, info] = compactCells(fN, fD, FC, EcQ)
% compactCells  Drop the rows a cell does not need, so that cells carrying one function can MERGE.
%
% objective: the arrangement emits one cell per sign pattern, and the refinement that picks a winner
%   adds a row per comparison -- so the same function arrives spread over many cells that differ in
%   SEVERAL rows. `mergeAdjacentCells` joins two cells only when they differ in ONE, which is an
%   exact set identity and is all a pairwise fold generates; it cannot touch this. Measured on the
%   convex fixture: 321 faces where a handful is right.
%
% [input]  fN, fD : the face functions;  FC : cell array of [edgeIndex, side, branch] rows
%          EcQ    : the shared curve table, k x 6
% [output] the same, compacted, and info with .dropped and .merged
%
% ==================================================================================================
% THE RULE: A ROW IMPLIED BY THE OTHERS IS NOT A BOUNDARY.
%
% Row r of a cell can go exactly when the rest of the cell already lies on r's side, i.e. when
%
%       { other rows } and { r VIOLATED }   has no interior,
%
% which `cellInterior` decides exactly -- including for a BRANCH row, one named line of a conjugate
% pair, which no earlier test could take. Dropping r then changes the cell by a set with empty
% interior, the same standard by which cells with empty interior are dropped in the first place, and
% it is what lets `mergeAdjacentCells` see two cells as differing in a single row.
%
% The two steps alternate to a fixpoint: dropping a row can expose a merge, and a merge drops a row
% and can expose another redundancy. Nothing here is approximate and no convexity is assumed.
%
% WHAT IT DOES NOT DO. Two cells carrying one function whose union is genuinely not an H-form cell
% -- an L-shape, say -- still come out as two. That is not a defect of the merge but a fact about
% the union, and describing it needs a face to hold a union of cells, which the type does not.

    if nargin ~= 4
        error('compactCells:nargin', 'expects (fN, fD, FC, EcQ).');
    end
    info = struct('dropped', 0, 'merged', 0);
    if isempty(FC), return, end
    before = numel(FC);

    for pass = 1:8
        drops = 0;
        for k = 1:numel(FC)
            [FC{k}, d] = dropImplied(FC{k}, EcQ);
            drops = drops + d;
        end
        info.dropped = info.dropped + drops;
        n0 = numel(FC);
        [fN, fD, FC] = mergeAdjacentCells(fN, fD, FC);
        if drops == 0 && numel(FC) == n0, break, end
    end
    info.merged = before - numel(FC);
end

% ==================================================================================================
function [rows, dropped] = dropImplied(rows, EcQ)
% Every row whose violation is already impossible, removed one at a time -- each decision is made
% against the rows that REMAIN, so no two rows can be dropped on the strength of each other.
    dropped = 0;
    if size(rows,1) < 2, return, end

    % ONLY A CELL WITH INTERIOR MAY LOSE A ROW. A THIN face -- dom f* is a line when the primal is a
    % positive semidefinite singular quadratic, a single point when it is affine -- has no interior
    % at all, so EVERY row's violation is vacuously interior-free and the rule below would strip it
    % to a half-plane. Measured as a red before this guard existed: five fast-bucket tests, among
    % them `conjQ`'s single-point dual domain and `biconjQ`'s collinear chain.
    try
        [C0, local0] = localTable(rows, EcQ);
        if ~cellInterior(C0, local0), return, end
    catch
        return                                       % undecidable: change nothing
    end

    t = 1;
    while t <= size(rows,1)
        if size(rows,1) < 2 || rows(t,2) == 0        % an EQUALITY row is the face, never redundant
            t = t + 1;  continue
        end
        try
            flipped = rows;
            flipped(t,2) = -flipped(t,2);
            [C, local] = localTable(flipped, EcQ);
            if ~cellInterior(C, local)
                rows(t,:) = [];
                dropped = dropped + 1;
                continue                              % t now indexes the next row
            end
        catch ME
            if ~any(strcmp(ME.identifier, {'ratQ:overflow', 'cellInterior:constant'}))
                rethrow(ME)
            end
            % Undecidable here: KEEP the row. A kept row costs a face, a dropped one costs area.
        end
        t = t + 1;
    end
end

% --------------------------------------------------------------------------------------------------
function [C, local] = localTable(rows, EcQ)
% `cellInterior` keys a curve on (conic, branch), while FC carries the branch per ROW -- so the same
% conic used on two branches is two curves here.
    C = zeros(0,7);
    local = zeros(size(rows,1), 2);
    for t = 1:size(rows,1)
        br = 0;
        if size(rows,2) >= 3, br = rows(t,3); end
        key = [EcQ(rows(t,1), :), br];
        hit = find(all(C == key, 2), 1);
        if isempty(hit), C(end+1,:) = key;  hit = size(C,1); end %#ok<AGROW>
        local(t,:) = [hit, rows(t,2)];
    end
end
