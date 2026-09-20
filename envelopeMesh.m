function [g, info] = envelopeMesh(env)
% envelopeMesh  Turn an assembled envelope's FACE LIST into one mesh, by building the ARRANGEMENT of
%   its regions once rather than folding the faces together.
%
% objective: the last step between a verified envelope and an object `biconj` can return. Both fold
%   orders are ruled out by measurement -- left to right reaches 3726 faces by its twelfth cell pair,
%   a balanced tree 570 by its third level -- and for a reason no ordering fixes: a max over a
%   SUBSET of the pairs is a partial max, and a partial max of piecewise quadratics carries pieces
%   that only the full max deletes. An arrangement never forms a partial max at all.
%
% [input]  env : as `envelopeAssemble` returns it, ideally after `envelopePrune`
% [output] g   : QuaCon -- the envelope, where every face of it is quadratic
%          info : .curves, .cells, .refined, .faces
%
% ==================================================================================================
% THE CONSTRUCTION, AND WHY IT IS SMALL
%
% 1. COLLECT THE DISTINCT CURVES the faces' regions are built from. Measured on the standard
%    fixtures after pruning: 13, 25 and 16 of them, against 52, 94 and 38 faces -- the same curve
%    bounds many faces, which is what makes an arrangement the cheap way round.
%
% 2. BUILD THE ARRANGEMENT INCREMENTALLY. Start from the whole plane and add one curve at a time,
%    splitting every cell into its two sides and DROPPING the sides that are empty. Emptiness is
%    decided exactly -- `ratQ.feasible2` for a cell that is all half-planes, `cellHasInterior`'s
%    degree-4 sign kernel when a conic is present -- so the cell list stays the size the geometry
%    actually is.
%
% 3. A CELL'S APPLICABLE FACES ARE READ OFF ITS SIGNS, with no further test. Every face's region is
%    built from these same curves, so a face applies on a cell exactly when the cell's sign for each
%    of that face's rows matches what the row demands. That is a comparison of two sign vectors.
%
% 4. REFINE EACH CELL BY THE WINNER. Where several faces apply, which is largest varies across the
%    cell, and the boundary is their difference -- a CONIC, because both are quadratics over
%    positive denominators. Splitting by it is one more row, and the faces are taken one at a time so
%    the refinement is bounded by how many actually apply (two to four here, not fifty).
%
% CORNER FACES ARE COVERAGE, NOT VALUE -- which is why they can be here at all.
%
% A corner face's value is affine with quadratic-irrational coefficients, so it is not a face any
% conic subdivision can store, and the boundary between it and a quadratic face would be a conic
% with irrational coefficients -- not an edge either. It turns out neither is needed. MEASURED on
% the c = -2 fixture, 4000 points: wherever a corner face owns a point, its value NEVER EXCEEDS the
% best quadratic face owning the same point, and the largest excess is exactly 0.0e0 -- equality,
% not domination, because the corner lies ON the tie and therefore belongs to the neighbouring pair
% too, whose own candidate through that same dual point takes the same value.
%
% So a corner face is carried for its REGION alone. Its rows join the arrangement, and what they
% decide is COVERAGE: by `envelopeEval`'s rule the envelope is the max over faces, and a PAIR with
% no face owning x makes it +infinity there. A cell where every pair is covered emits the largest
% QUADRATIC face; a cell where some pair is not is left out, and a hole in the mesh is exactly
% +infinity. A cell covered ONLY by a corner face would need an algebraic face after all, and is
% refused by name (`envelopeMesh:algebraicFace`) rather than silently mis-valued -- it did not occur
% on any fixture here.
%
% WHAT IT STILL DOES NOT COVER: a BRIDGE face, whose region is a predicate with no rows at all
% (DECISIONS.md 2026-09-17). `envelopeMesh:bridgeFaces` names that, and `envelopeEval` answers for
% such an envelope.

    if nargin ~= 1
        error('envelopeMesh:nargin', 'expects (env).');
    end
    if ~isstruct(env) || ~all(isfield(env, {'quad', 'bridge'}))
        error('envelopeMesh:shape', 'env must be what envelopeAssemble returns.');
    end
    corners = {};
    for k = 1:numel(env.bridge)
        f = env.bridge{k};
        if strcmp(f.kind, 'corner') && ~isempty(f.con)
            corners{end+1} = f; %#ok<AGROW>
        else
            error('envelopeMesh:bridgeFaces', ...
                ['this envelope has a face whose region is a PREDICATE with no rows (kind %s), ' ...
                 'so it has no signs to place in an arrangement. Mesh the quadratic part by ' ...
                 'passing an envelope without it, or keep using envelopeEval, which answers for ' ...
                 'both.'], f.kind);
        end
    end
    if isempty(env.quad)
        g = QuaCon(zeros(0,3), zeros(0,3), zeros(0,6), zeros(0,10), zeros(0,1), zeros(0,2), {});
        info = struct('curves', 0, 'cells', 0, 'refined', 0, 'faces', 0);
        return
    end

    % ---- 1. the distinct curves, from the quadratic faces AND the corners' wedges -----------------
    [C, req, pairOf] = curveTable(env.quad, corners);
    nQuad = numel(env.quad);
    allPairs = unique(pairOf);

    % ---- 2. the arrangement, one curve at a time ------------------------------------------------
    cells = {zeros(0,2)};                       % each cell: rows [curveIndex, side]
    for j = 1:size(C,1)
        next = {};
        for k = 1:numel(cells)
            for sd = [1 -1]
                cand = [cells{k}; j, sd];
                if cellIsLive(C, cand), next{end+1} = cand; end %#ok<AGROW>
            end
        end
        cells = next;
    end
    nArr = numel(cells);

    % ---- 3 and 4. the applicable faces, the COVERAGE test, and the winner ------------------------
    out = struct('num', {}, 'den', {}, 'con', {});
    nRef = 0;  nUncovered = 0;
    for k = 1:numel(cells)
        app = applicable(req, cells{k});
        if isempty(app), nUncovered = nUncovered + 1; continue, end
        % EVERY pair must be covered here, or the envelope is +infinity and the cell is left out.
        % This is `envelopeEval`'s rule, and it is what the corner faces are carried for.
        if ~all(ismember(allPairs, pairOf(app))), nUncovered = nUncovered + 1; continue, end
        quadApp = app(app <= nQuad);
        if isempty(quadApp)
            error('envelopeMesh:algebraicFace', ...
                ['a cell is covered only by CORNER faces, whose value is affine with ' ...
                 'quadratic-irrational coefficients and cannot be stored on a conic subdivision. ' ...
                 'Measured never to occur on the fixtures here -- the corner''s value always ties ' ...
                 'a quadratic face -- so reaching it means the output type needs an algebraic ' ...
                 'face after all.']);
        end
        [sub, C] = winnerCells(C, cells{k}, env.quad, quadApp);
        nRef = nRef + numel(sub);
        for t = 1:numel(sub)
            out(end+1) = sub(t); %#ok<AGROW>
        end
    end

    g = assembleQuaConCells(out);
    info = struct('curves', size(C,1), 'cells', nArr, 'refined', nRef, 'faces', g.nf, ...
                  'uncovered', nUncovered, 'corners', numel(corners));
end

% ==================================================================================================
function [C, req, pairOf] = curveTable(faces, corners)
% The distinct curves, each face's demand on them (req{i}, an m x 2 list [curveIndex, side]), and
% which PAIR each face belongs to -- the quadratic faces first, then the corners, so that a face
% index at most numel(faces) is one whose VALUE can be emitted.
% Bitwise deduplication, because the rows are canonical integers -- the same argument
% assembleQuaConCells' own curve table rests on.
    C = zeros(0,7);                                 % [conic(1:6), branch] -- branch 0 for a curve
    n = numel(faces) + numel(corners);
    req = cell(1, n);
    pairOf = zeros(1, n);
    for i = 1:n
        if i <= numel(faces)
            con = faces(i).con;
            pairOf(i) = faces(i).pair(1);
        else
            con = corners{i - numel(faces)}.con;
            pairOf(i) = corners{i - numel(faces)}.pairIndex;
        end
        r = zeros(size(con,1), 2);
        keep = true(size(con,1), 1);
        for t = 1:size(con,1)
            c = con(t,1:6);
            if all(c(1:5) == 0)
                % A CONSTANT row is not a curve. It is vacuous when it holds and empties the face
                % when it does not, and admitting one makes the exact kernel meet the "curve" 1 = 0.
                % Measured as reachable on the c = -2 fixture (DECISIONS.md 2026-09-18).
                if c(6) < 0, req{i} = [];  keep = false(size(con,1),1);  break, end
                keep(t) = false;  continue
            end
            key = [c, branchOf(con, t)];
            hit = find(all(C == key, 2), 1);
            if isempty(hit), C(end+1,:) = key;  hit = size(C,1); end %#ok<AGROW>
            r(t,:) = [hit, con(t,7)];
        end
        req{i} = r(keep,:);
    end
end

% --------------------------------------------------------------------------------------------------
function tf = cellIsLive(C, rows)
% Does this cell have a two-dimensional interior? `cellInterior` decides it, including the case a
% plain Fourier-Motzkin cannot reach -- a BRANCH row, one named line of a conjugate pair. Measured:
% ignoring such a row costs 31104 live cells against 486 (DECISIONS.md 2026-09-18), because every
% cell kept is split again by every later curve.
    tf = cellInterior(C, rows);
end

% --------------------------------------------------------------------------------------------------
function b = branchOf(con, t)
% The row's branch, where one is carried. Seven-column rows predate the branch column and name an
% ordinary curve.
    if size(con,2) >= 8, b = con(t,8); else, b = 0; end
end

% --------------------------------------------------------------------------------------------------
function app = applicable(req, rows)
% Which faces apply on this cell: those every one of whose demands the cell's signs already meet.
% No feasibility test is needed -- the face's rows ARE curves of the arrangement, so the cell has
% already chosen a side of each.
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
function [out, C] = winnerCells(C, rows, faces, app)
% Split this cell until ONE of the applicable faces is the largest on each piece, taking the faces
% one at a time. The boundary between two of them is their DIFFERENCE, a conic, because both are
% quadratics over positive denominators -- so each split is one more row, and the difference joins
% the curve table like any other curve so that a cell is always just a list of [curve, side].
    out = struct('num', {}, 'den', {}, 'con', {});
    open = {struct('rows', rows, 'win', app(1))};
    for a = 2:numel(app)
        j = app(a);
        nxt = {};
        for k = 1:numel(open)
            cur = open{k};
            [d, sgn] = diffConicOf(faces(j), faces(cur.win));   % d > 0 (as spelled) means j wins
            if isempty(d)
                nxt{end+1} = cur;  continue                     %#ok<AGROW> identical functions
            end
            [C, idx] = addCurve(C, d);
            for sd = [1 -1]
                cand = cur;
                cand.rows = [cur.rows; idx, sd*sgn];
                if sd > 0, cand.win = j; end
                if cellIsLive(C, cand.rows), nxt{end+1} = cand; end %#ok<AGROW>
            end
        end
        open = nxt;
    end
    for k = 1:numel(open)
        f = faces(open{k}.win);
        out(end+1) = struct('num', f.num, 'den', f.den, ...
                            'con', rowsToCon(C, open{k}.rows)); %#ok<AGROW>
    end
end

% --------------------------------------------------------------------------------------------------
function [c, sgn] = diffConicOf(fj, fi)
% The curve { v_j = v_i }, and the SIGN that means "j is the larger" once the conic is spelled
% canonically. Both denominators are positive integers, so the difference's sign is the sign of
% num_j*den_i - num_i*den_j; `ratQ.conic` may flip the whole conic to make its leading entry
% positive, and `sgn` carries that flip.
    num = ratQ.chk(fj.num * fi.den - fi.num * fj.den, 'face difference');
    raw = [num(5), 2*num(6), num(7), 2*num(8), 2*num(9), 2*num(10)];
    if all(raw == 0) || any(num(1:4) ~= 0), c = [];  sgn = 0;  return, end
    nz = find(raw ~= 0, 1);
    c = ratQ.conic(raw);
    sgn = sign(raw(nz));
end

% --------------------------------------------------------------------------------------------------
function [C, idx] = addCurve(C, c)
% Find or add a curve in the table, so every cell row stays a plain [curveIndex, side]. A difference
% conic is an ordinary curve, so its branch column is 0.
    key = [c, 0];
    idx = find(all(C == key, 2), 1);
    if isempty(idx), C(end+1,:) = key;  idx = size(C,1); end
end

% --------------------------------------------------------------------------------------------------
function con = rowsToCon(C, rows)
% The cell's rows as assembleQuaConCells' 8-column form.
    con = zeros(size(rows,1), 8);
    for t = 1:size(rows,1)
        con(t,1:6) = C(rows(t,1), 1:6);
        con(t,7)   = rows(t,2);
        con(t,8)   = C(rows(t,1), 7);
    end
end
