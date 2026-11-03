function [subdivision, info] = envelopeSubdivisionAlg(env, dom, cuts, opts)
% envelopeSubdivisionAlg  The envelope as a subdivision whose EDGES may be algebraic -- the arrangement that
%   `envelopeSubdivision` refuses when a bridge face competes.
%
% objective: a bridge face carries the envelope's value, and where it competes with a quadratic face
%   the boundary is `FP(x, y, q(x,y)) = 0`, of degree four. `envelopeSubdivision` raises
%   `envelopeSubdivision:bridgeCompetes` there because a conic subdivision cannot hold that edge. This
%   builds the same arrangement with such edges allowed, so the GEOMETRY can be verified before the
%   output TYPE is grown to store it.
%
% [input]  env  : as `envelopeAssemble` returns it, ideally pruned
%          dom  : (optional) rows [a b c], a x + b y + c >= 0, the subdivision is clipped to
%          cuts : (optional) integer lines [a b c] added to the arrangement as curves, so that every
%                 face lies on one side of each and records it -- cutTriangles' rational lines
%          opts : (optional) struct; .overlayBudget (default 4000, Inf lifts it) bounds the pieces
%                 the overlay of several pairs may visit -- see minOverlay
%                 .maxCells (default 3000) bounds the arrangement's cells -- arrangementSize
% [output] subdivision : struct array, one entry per face, with
%                   .kind  'quad' | 'bridge'
%                   .num, .den   the value, when kind is 'quad'
%                   .FP, .FB     the value, when kind is 'bridge'
%                   .con   conic rows [a b c d e f sign branch], as everywhere else
%                   .alg   algebraic rows, a struct array with .FP .FB .num .den .side, meaning
%                          side * (q(x) - value_FB(x)) >= 0
%                   .ii    (several pairs only) rows ordering TWO irrational values, a struct array
%                          with .FP1 .FB1 .FP2 .FB2, meaning value_1(x) >= value_2(x)
%          info : .cells, .faces, .algRows, .refused, .skipped; one pair also .cellCon (each cell
%                 with an interior point, as con rows) and .faceCell (the cell of each face)
%
% ==================================================================================================
% WHAT IS NEW HERE, AND WHAT IS DELIBERATELY THE SAME
%
% THE SAME: the arrangement is built from the CONIC curves exactly as before -- the quadratic faces'
% rows, the corner wedges, the bridges' outer bounds -- and every cell's emptiness is decided by
% `cellInterior`. Corner faces are carried for COVERAGE only, because a corner never carries the
% value (measured, DECISIONS.md 2026-09-19), and a pair with no face owning a cell makes the
% envelope +infinity there, which is a hole in the subdivision.
%
% NEW: where two applicable faces cannot be separated by a conic, the cell is split by an ALGEBRAIC
% row instead of being refused. Such a row is `(FP, FB, num, den, side)` -- the curve is the rational
% `bridgeQuadCut` polynomial, and `side` picks one of its two branches, because that polynomial names
% the PAIR and not one of them (DECISIONS.md 2026-09-20). `algCurveSide` decides it exactly at a
% rational point; `subdivisionAlgEval` reads it in doubles, which is what evaluation is allowed to do.
%
% SEVERAL IRRATIONAL faces on one cell (bridges, and irrational CORNERS, which carry value too --
% 2026-10-17) are not cut against each other: that curve is not of this shape. Each is cut against
% the quadratics only, so where two of them both beat every quadratic their pieces OVERLAP; the
% value is still right because every evaluator takes the max over the faces owning a point.
% `info.overlaps` counts such cells (a fixed-grid read). Until 2026-10-17 this was refused as
% `envelopeSubdivisionAlg:twoBridges`; erratum Example 1 reaches it with two corners from Q(sqrt 26) and
% Q(sqrt 6), and 18 of its cells overlap.

    if nargin < 1 || nargin > 4
        error('envelopeSubdivisionAlg:nargin', 'expects (env), (env, dom), (env, dom, cuts) or (env, dom, cuts, opts).');
    end
    if nargin < 2, dom = zeros(0,3); end
    if nargin < 3, cuts = zeros(0,3); end
    if nargin < 4 || isempty(opts), opts = struct(); end
    if ~isfield(opts, 'overlayBudget'), opts.overlayBudget = 4000; end
    if ~isfield(opts, 'maxCells'), opts.maxCells = 3000; end
    % SEVERAL PAIRS: a cell array of pair envelopes is the MINIMUM over them -- an EnvelopeMin with no
    % triangle face. Each pair is subdivided ALONE and the subdivisions are OVERLAID (minOverlay,
    % 2026-10-31). Until then all the pairs' curves went into one arrangement resolved per cell as a
    % min of maxes; on Example 2 that passed 60000 cells, and even arranged locally (2978 cells) a cell
    % met 41-47 faces and took 6-270 s, while each pair alone took 33-69 s and 44-85 cells.
    if iscell(env), envs = env; else, envs = {env}; end
    for g = 1:numel(envs)
        if ~isstruct(envs{g}) || ~all(isfield(envs{g}, {'quad', 'bridge'}))
            error('envelopeSubdivisionAlg:shape', 'env must be what envelopeAssemble returns.');
        end
    end
    if numel(envs) > 1
        subs = cell(1, numel(envs));  infos = cell(1, numel(envs));
        for g = 1:numel(envs)
            [subs{g}, infos{g}] = envelopeSubdivisionAlg(envs{g}, dom, cuts, opts);
        end
        [subdivision, info] = minOverlay(subs, infos, opts.overlayBudget);
        return
    end
    env = combinePairs(envs);
    [corners, bridges] = splitBridge(env.bridge);
    [C, req, pairOf] = conicTable(env.quad, corners);
    nQuad = numel(env.quad);
    allPairs = unique([pairOf{:}, cellfun(@(f) f.pairIndex, bridges)]);
    [C, cEnds] = addCornerCrossings(C, corners);
    [C, bEnds] = addBridgeCurves(C, bridges);   % without them a bridge's membership can change inside a cell
    [C, lEnds] = splitLineCurves(C);
    for i = 1:size(cuts,1)
        C = addKey(C, [ratQ.conic([0 0 0 cuts(i,:)]), 0]);
    end
    owners = curveOwners(C, req, nQuad, bridges, cEnds, bEnds, lEnds);

    % ---- the arrangement over the CONIC curves ---------------------------------------------------
    % THREE THINGS KEEP IT SMALL (2026-10-17, the poly2 kink at c = -2: 85 curves, >90 min before).
    %   * CLIPPED to `dom` when the caller knows it -- dom f** is conv(dom f) for a bounded f, and
    %     every cell outside it is a hole anyway (316 of 332 cells on an indefinite triangle);
    %   * a row is RECORDED only where its curve splits the cell. A curve missing a cell has one sign
    %     on it, read at the cell's interior point when a face demands it (`applicable`), so it adds
    %     nothing to every later emptiness test the cell goes through;
    %   * a CELL budget, not a curve budget: what costs is the cells, and clipped, a curve count says
    %     little about them.
    cells = {zeros(0,2)};
    if ~isempty(dom)
        base = zeros(0,2);
        for i = 1:size(dom,1)
            raw = [0 0 0 dom(i,:)];
            key = [ratQ.conic(raw), 0];
            hit = find(all(C == key, 2), 1);
            if isempty(hit), C(end+1,:) = key;  hit = size(C,1); end %#ok<AGROW>
            base(end+1,:) = [hit, sign(raw(find(raw ~= 0, 1)))]; %#ok<AGROW>
        end
        cells = {base};
    end
    % A LOCAL ARRANGEMENT (2026-10-31). A curve splits a cell only where it can change something
    % there: some face (or bridge block) it bounds must be satisfiable on the cell by its OTHER
    % rows. Where none is, every face it bounds is inapplicable on the whole cell whichever side of
    % the curve a point lies, so the cell is not split, and the curve is remembered as SKIPPED -- it
    % may cross the cell, so it is never emitted as a row of it (fullRows). Why applicability stays
    % constant on every final cell: a face whose row j was skipped on an ancestor A has no point of A
    % meeting its other rows, so it applies nowhere on A; every other face has all its rows recorded
    % or one-signed. Emptiness is liveCell's, which is sound; an undecided test splits.
    % Example 2's three pairs passed 60000 cells globally (DECISIONS 2026-10-30 (seventh)).
    skips = {zeros(1,0)};
    nSkip = 0;
    maxCells = opts.maxCells;
    for j = 1:size(C,1)
        next = {};  nextSk = {};
        for k = 1:numel(cells)
            if any(cells{k}(:,1) == j), next{end+1} = cells{k};  nextSk{end+1} = skips{k}; continue, end %#ok<AGROW>
            live = [liveCell(C, [cells{k}; j, 1]), liveCell(C, [cells{k}; j, -1])];
            if all(live) && ~relevantHere(C, owners{j}, cells{k})
                next{end+1} = cells{k};  nextSk{end+1} = [skips{k}, j]; %#ok<AGROW>
                nSkip = nSkip + 1;
            elseif all(live)
                next{end+1} = [cells{k}; j, 1];  next{end+1} = [cells{k}; j, -1]; %#ok<AGROW>
                nextSk{end+1} = skips{k};  nextSk{end+1} = skips{k}; %#ok<AGROW>
            elseif any(live)
                next{end+1} = cells{k};  nextSk{end+1} = skips{k}; %#ok<AGROW>   % the curve misses this cell
            end
        end
        cells = next;  skips = nextSk;
        if numel(cells) > maxCells
            error('envelopeSubdivisionAlg:arrangementSize', ...
                'the arrangement passed %d cells at curve %d of %d; envelopeEval answers for the face list.', ...
                maxCells, j, size(C,1));
        end
    end

    % ---- one face per cell, with the winner decided and the ties cut algebraically ----------------
    subdivision = struct('kind', {}, 'num', {}, 'den', {}, 'FP', {}, 'FB', {}, 'con', {}, 'alg', {});
    nAlg = 0;  nUncovered = 0;  nOverlap = 0;  nNoPoint = 0;  nBelow = 0;
    cellCon = {};  faceCell = zeros(1, 0);      % the cells, for an overlay of several pairs
    for k = 1:numel(cells)
        rows = cells{k};
        [pt, okPt] = interiorPoint(C, rows);
        if ~okPt, nUncovered = nUncovered + 1;  nNoPoint = nNoPoint + 1; continue, end
        cellCon{end+1} = rowsToCon(C, rows); %#ok<AGROW>

        app = applicable(req, rows, C, pt);
        here = [pairOf{app}];
        bHere = [];
        % A BRIDGE'S MEMBERSHIP BY MAJORITY over five spread interior points (2026-10-17). One point,
        % rounded onto a dyadic grid, can land on a special rational line where a DEGENERATE bridge
        % (two rank-1 ends) owns a set of measure zero; the poly2 kink gave such a bridge a whole
        % cell it owns no part of (0 of 1922 samples), 2.7e-4 too high. The predicate is constant on
        % a cell up to measure zero, so a majority reads the cell's value.
        [Q, Pw] = votePoints(C, rows, pt);
        % THE FULL SIGN VECTOR before anything is emitted: the arrangement records only the curves
        % that split a cell, but a face's rows are what its boundary is traced from, and with the
        % sparse rows algCrvSubdivisionTest found a face vertex ending other than two arcs (2026-10-17). A
        % curve the cell does not record has one sign on it -- read at the cell's points.
        rowsOut = fullRows(C, rows, Q, skips{k});      % emitted; the exact tests below keep the sparse rows
        for b = 1:numel(bridges)
            if sum(bridgeMember(bridges{b}.pair, Q)) * 2 > size(Q,1)
                here(end+1) = bridges{b}.pairIndex; %#ok<AGROW>
                bHere(end+1) = b; %#ok<AGROW>
            end
        end
        if ~all(ismember(allPairs, here)), nUncovered = nUncovered + 1; continue, end

        quadApp = app(app <= nQuad);
        % THE IRRATIONAL CANDIDATES: the bridges owning the cell, and -- since 2026-10-17 -- the
        % CORNERS that apply to it. "A corner never carries the value" was measured on one fixture
        % (DECISIONS.md 2026-09-19) and is false in general: on every one-pair atom envelope a
        % corner is the maximum in slivers near a domain vertex, by up to 8e-6, and ignoring it left
        % the subdivision that much too LOW there.
        irr = cell(1, 0);
        for b = bHere, irr{end+1} = struct('FP', bridges{b}.FP, 'FB', bridges{b}.FB); end %#ok<AGROW>
        for j = app(app > nQuad)
            [fc, okc] = cornerValueFace(corners{j - nQuad}, pt);
            if okc, irr{end+1} = fc; end %#ok<AGROW>
        end
        if isempty(quadApp) && isempty(irr), nUncovered = nUncovered + 1; continue, end
        if numel(irr) > 1
            [irr, nDrop] = dropBelow(irr, C, rows, Pw, env.quad, quadApp);
            nBelow = nBelow + nDrop;
        end
        if numel(irr) > 1 && ~dominantIrrational(irr, C, rows, pt, env.quad, quadApp)
            nOverlap = nOverlap + 1;
        end

        sub = winnerFaces(env.quad, quadApp, irr, rows, C, pt, Pw, rowsOut);
        for t = 1:numel(sub)
            nAlg = nAlg + numel(sub(t).alg);
            subdivision(end+1) = sub(t); %#ok<AGROW>
            faceCell(end+1) = numel(cellCon); %#ok<AGROW>
        end
    end
    info = struct('cells', numel(cells), 'faces', numel(subdivision), 'algRows', nAlg, ...
                  'refused', nUncovered, 'noPoint', nNoPoint, 'overlaps', nOverlap, 'below', nBelow, ...
                  'skipped', nSkip, ...
                  'cellCon', {cellCon}, 'faceCell', faceCell);
end

% ==================================================================================================
function out = winnerFaces(quad, quadApp, irr, rows, C, pt, Q, rowsOut)
% Which applicable face is largest, and the rows that say where. A quad-versus-quad boundary is a
% CONIC and joins the ordinary rows; a quad-versus-bridge boundary is ALGEBRAIC and joins `.alg`.
    out = struct('kind', {}, 'num', {}, 'den', {}, 'FP', {}, 'FB', {}, 'con', {}, 'alg', {});
    base = rowsToCon(C, rowsOut);

    if isempty(irr)
        % No bridge here: the largest quadratic wins, and ties between two of them are conics --
        % one piece PER quadratic (2026-10-17): the winner read at the cell's one point owns only
        % part of the cell, and the rest was left to no face, a hole.
        for j = quadApp(:).'
            extra = quadOver(quad, quadApp, j);
            if ~pieceLive(C, rows, extra, Q), continue, end
            out(end+1) = struct('kind', 'quad', 'num', quad(j).num, 'den', quad(j).den, ...
                                'FP', [], 'FB', 0, 'con', [base; extra], 'alg', emptyAlg()); %#ok<AGROW>
        end
        return
    end

    % IRRATIONAL faces (bridges, corners) AND quadratics: each irrational face owns where it is above
    % every quadratic -- one algebraic row per comparison -- and each quadratic owns where it is the
    % largest quadratic and above every irrational face. Two irrational faces are NOT cut against
    % each other (that curve is not of the irrational-quadratic shape); where both beat every
    % quadratic their pieces overlap, and the value is still right because AlgCrv.eval takes the max
    % over the faces owning a point (counted as `overlaps`, 2026-10-17).
    for f = irr
        algRows = emptyAlg();
        for j = quadApp(:).'
            algRows(end+1) = struct('FP', f{1}.FP, 'FB', f{1}.FB, 'num', quad(j).num, ...
                                    'den', quad(j).den, 'side', -1); %#ok<AGROW>
        end
        out(end+1) = struct('kind', 'bridge', 'num', [], 'den', 0, 'FP', f{1}.FP, 'FB', f{1}.FB, ...
                            'con', base, 'alg', algRows); %#ok<AGROW>
    end

    % ... and, for EACH quadratic, the piece where it is the largest quadratic and beats every
    % irrational face (one piece for the winner at the cell's point left the rest to no face).
    for j = quadApp(:).'
        extra = quadOver(quad, quadApp, j);
        if ~pieceLive(C, rows, extra, Q), continue, end
        algRows = emptyAlg();
        for f = irr
            algRows(end+1) = struct('FP', f{1}.FP, 'FB', f{1}.FB, 'num', quad(j).num, ...
                                    'den', quad(j).den, 'side', 1); %#ok<AGROW>
        end
        if algEmpty(C, rowsOut, extra, algRows, Q), continue, end    % the face's own rows, see algEmpty
        out(end+1) = struct('kind', 'quad', 'num', quad(j).num, 'den', quad(j).den, ...
                            'FP', [], 'FB', 0, 'con', [base; extra], 'alg', algRows); %#ok<AGROW>
    end
end

function rows = fullRows(C, rows, P, skipped)
% Every curve's side on the cell: recorded ones as recorded, the rest by the sign at the cell's points
% (the majority of the nonzero signs; a curve tangent at one point does not decide it). A SKIPPED
% curve may cross the cell (it bounds nothing that applies here) and is left out.
    have = false(size(C,1), 1);  have(rows(:,1)) = true;  have(skipped) = true;
    for j = find(~have).'
        c = C(j, 1:6);
        if C(j,7) ~= 0
            v = sign(ratQ.pairLineEval(c, C(j,7), P));
        else
            v = sign(c(1)*P(:,1).^2 + c(2)*P(:,1).*P(:,2) + c(3)*P(:,2).^2 + c(4)*P(:,1) + c(5)*P(:,2) + c(6));
        end
        v = v(v ~= 0);
        if isempty(v), continue, end                  % on the curve throughout: leave it unrecorded
        rows(end+1,:) = [j, sign(sum(v))]; %#ok<AGROW>
    end
end

function [Q, P] = votePoints(C, rows, pt)
% The cell's interior point and four more spread over the cell (from cellSamplePoint's grid), and
% the whole grid sample P for a caller looking for witnesses.
    Q = pt;
    [~, ok, P] = cellSamplePoint(C, rows, straightRows(C, rows), 1e3);
    if ok && size(P,1) >= 4
        Q = [pt; P(round(linspace(1, size(P,1), 4)), :)];
    end
    P = [pt; P];
end

function tf = liveCell(C, rows)
% cellInterior, except that a cell it KEEPS ONLY BECAUSE IT CANNOT DECIDE (a tie mixing two quadratic
% fields, `unsureTie`) must show an interior point. Kept undecided, such cells were 2938 of 2974 on the
% 3-convex-edge triangle's two-pair arrangement -- slivers pinched at the domain's corners, each split
% again by every later curve -- and the cost of the arrangement was theirs (2026-10-18).
    %
    % AND A CELL MIXING CURVED AND BRANCH ROWS likewise. cellInterior decides the curved rows with the
    % straight ones (exactly) and the branch rows with the straight ones (a margin), but never the
    % curved and branch rows TOGETHER -- so a cell empty only because a conic and a branch line
    % exclude each other passes both. All 2206 refused cells of that arrangement were such (path
    % 'margin', no interior point at a 192-step grid).
    [tf, info] = cellInterior(C, rows);
    if ~tf, return, end
    mixed = any(C(rows(:,1), 7) ~= 0) && any(C(rows(:,1), 7) == 0 & any(C(rows(:,1), 1:3) ~= 0, 2));
    if strcmp(info.path, 'unsureTie') || mixed
        [~, tf] = cellSamplePoint(C, rows, straightRows(C, rows), 1e3);
    end
end

function tf = algEmpty(C, rows, conRows, algRows, Q)
% Is the piece PROVABLY empty because of one of its algebraic rows? pieceLive decides the conic rows
% only, and on the 3-convex-edge triangle 2040 of 2172 faces were empty by an algebraic row.
% A row failing at every sample point of the piece's conic part is a CANDIDATE; algRowFails then
% certifies it false at every point of that part, by Bernstein signs on sub-boxes (2026-10-29).
%
% Until 2026-10-29 the certificate was R = FP(x,y,q) != 0 and D of one sign on the CELL, read at the
% sample points (2026-10-20). That left the ~1600 faces whose cell R = 0 crosses, it relied on a
% sample reaching every connected part of the piece (one sample per part had opened five holes),
% and bridgeQuadCut's R dropped the basis' 1/2 on x^2 and y^2, so on faces with a square term it was
% not the cut at all. The pointwise certificate needs no sample and no connectivity argument.
%
% `rows` are the FULL rows the face is emitted with (fullRows: every curve's side on the cell), not
% the sparse ones the arrangement recorded (2026-10-30). The region is the same set; but the rows
% the cell does not record -- x >= 0, x + y >= 0 through the triangle's corner -- are the factors
% algRowFails divides out, and with the sparse rows 132 faces at that corner stayed.
    tf = false;
    if isempty(algRows), return, end
    M = true(size(Q,1),1);
    for r = 1:size(conRows,1), M = M & conRows(r,7) * QuaPar.evalConic(conRows(r,1:6), Q) > 0; end
    Qk = Q(M,:);
    % A WITNESS must be STRICT (2026-10-30): read with piecePoints' 1e-9 tolerance, a point where a
    % touching row holds with equality kept zero-area faces. And NO sample in the piece is no reason
    % to stop: the certificate needs none, and returning there left ~250 faces failing by a clear
    % margin (each certified at depth 6 even on the whole domain's box) never tried.
    holds = true(size(Qk,1), 1);
    nFail = zeros(numel(algRows), 1);
    for t = 1:numel(algRows)
        a = algRows(t);
        z = algBranch(ratQ.bigToDouble(a.FP), Qk, a.FB);  q = ratQ.evalFace(a.num, a.den, Qk);   % doubles: see branchAt
        ok = isfinite(z) & a.side * (q - z) > 1e-9 * max(1, abs(z));
        holds = holds & ok;
        nFail(t) = sum(~ok);
    end
    if any(holds), return, end                       % a witness: not empty
    for t = 1:numel(algRows)
        if isempty(Qk) || nFail(t) < size(Qk,1), continue, end   % a lone row must fail at every point
        if rowFails(C, rows, conRows, algRows(t)), tf = true; return, end
    end
    % JOINTLY (2026-10-30): no sample has every row holding, but no one row fails at all of them --
    % 274 faces of the 3-convex-edge triangle. A sub-box is decided by whichever row fails on it; the
    % rows failing most often are tried first. With no sample in the piece, every row together.
    if numel(algRows) > 1 || isempty(Qk)
        [~, o] = sort(nFail, 'descend');
        tf = rowFails(C, rows, conRows, algRows(o));
    end
end

function ok = rowFails(C, rows, conRows, a)
% algRowFails (one row, or several jointly) on the piece's conic part (cell rows plus conRows), cached on the rows' CURVES, not
% their indices, so that a cache entry cannot outlive the arrangement it was computed on.
    persistent cache
    if isempty(cache), cache = containers.Map(); end
    ok = false;
    try
        key = sprintf('%s|%s', mat2str(sortrows([C(rows(:,1),:), rows(:,2)])), mat2str(sortrows(conRows)));
        for t = 1:numel(a)
            key = [key, sprintf('|%s|%s|%d|%d|%d', mat2str(a(t).FP(:).'), mat2str(a(t).num), ...
                                a(t).den, a(t).FB, a(t).side)]; %#ok<AGROW>
        end
    catch
        key = '';                                    % a BigInteger FP: not cached
    end
    if ~isempty(key) && isKey(cache, key), ok = cache(key); return, end
    [box, inBox, factors, corners] = regionBox(C, rows, conRows);
    if ~isempty(box)
        try
            ok = algRowFails(a, box, inBox, 6, factors, corners);
        catch
            ok = false;
        end
    end
    if ~isempty(key), cache(key) = ok; end
end

function [box, inBox, rowPolys, V] = regionBox(C, rows, K)
% A box around the cell (its straight rows, within the subdividers' bounding box) and a test rejecting
% sub-boxes that provably miss the cell's or the piece's conic rows or branch rows. rowPolys: the
% conic rows as polynomials >= 0 on the region, which algRowFails divides out of a row's F(q); V the
% straight rows' polygon vertices, the corners algRowFails blows up at (exact ones only certify).
    H = straightRows(C, rows);
    A = [H; 1 0 1e3; -1 0 1e3; 0 1 1e3; 0 -1 1e3];
    V = zeros(0,2);
    for i = 1:size(A,1)
        for j = i+1:size(A,1)
            Mx = A([i j], 1:2);
            if abs(det(Mx)) < 1e-12, continue, end
            x = (Mx \ -A([i j], 3)).';
            if all(A(:,1:2) * x.' + A(:,3) >= -1e-9 * max(1, norm(x))), V(end+1,:) = x; end %#ok<AGROW>
        end
    end
    if size(V,1) < 3, box = []; inBox = []; rowPolys = {}; V = zeros(0,2); return, end
    box = [min(V(:,1)) max(V(:,1)) min(V(:,2)) max(V(:,2))];
    rowPolys = {};
    for t = 1:size(rows,1)
        if C(rows(t,1), 7) ~= 0, continue, end
        rowPolys{end+1} = rows(t,2) * conicPoly(C(rows(t,1), 1:6)); %#ok<AGROW>
    end
    for r = 1:size(K,1)
        if size(K,2) >= 8 && K(r,8) ~= 0, continue, end
        rowPolys{end+1} = K(r,7) * conicPoly(K(r,1:6)); %#ok<AGROW>
    end
    % BRANCH rows too (2026-10-29): without them the box test kept the far side of a line pair's
    % chosen line, where an algebraic row may hold, and all 896 faces left undecided on the
    % 3-convex-edge triangle had one. A branch row is AFFINE, so its maximum on a box is at a corner;
    % the line is a double (pairLineEval, roots of a quadratic), so a box is rejected only when it is
    % clear of the line by a margin -- keeping a box is always safe.
    lines = zeros(0,8);
    for t = 1:size(rows,1)
        if C(rows(t,1), 7) ~= 0, lines(end+1,:) = [C(rows(t,1), 1:6), rows(t,2), C(rows(t,1), 7)]; end %#ok<AGROW>
    end
    if size(K,2) >= 8, lines = [lines; K(K(:,8) ~= 0, 1:8)]; end
    inBox = @(b) ~any(cellfun(@(P) polySignOn(P, b, @(~) true, 0) == -1, rowPolys)) && ...
                 ~lineMisses(lines, b);
end

function tf = lineMisses(lines, b)
% Does box b lie strictly on the wrong side of one of the branch rows?
    tf = false;
    X = [b(1) b(3); b(2) b(3); b(1) b(4); b(2) b(4)];
    for r = 1:size(lines,1)
        v = lines(r,7) * ratQ.pairLineEval(lines(r,1:6), lines(r,8), X);
        if max(v) < -1e-7 * max([1; abs(v)]), tf = true; return, end
    end
end

function P = conicPoly(c)
% [a b c d e f] as a dense polyN: P(i+1,j+1) multiplies x^i y^j.
    P = zeros(3,3);
    P(3,1) = c(1);  P(2,2) = c(2);  P(1,3) = c(3);  P(2,1) = c(4);  P(1,2) = c(5);  P(1,1) = c(6);
end

function tf = pieceLive(C, rows, extra, Q)
% Does the cell, cut to where quadratic j is the largest, keep an interior? Decided exactly by
% cellInterior on the cell's rows plus `extra` (2026-10-17): one piece per applicable quadratic
% made 9153 faces on the poly2 kink, most of them empty. A cell point strictly inside every extra
% row is a WITNESS of interior and spares the exact test (the points are inside the cell already).
    % A CONSTANT row (two quadratics differing by a constant) holds everywhere or nowhere: resolve
    % it here, since cellInterior takes curves only.
    k0 = all(extra(:,1:5) == 0, 2);
    if any(k0)
        if any(extra(k0,7) .* extra(k0,6) < 0), tf = false; return, end
        extra = extra(~k0, :);
    end
    if isempty(extra), tf = true; return, end
    inside = true(size(Q,1), 1);
    for r = 1:size(extra,1)
        inside = inside & extra(r,7) * QuaPar.evalConic(extra(r,1:6), Q) > 1e-9 * max(abs(extra(r,1:6)));
    end
    if any(inside), tf = true; return, end
    % ONE dominating quadratic empties the piece, and one extra row is a far cheaper exact test than
    % all of them: try the rows in order of how often they fail on the samples (2026-10-17).
    fails = zeros(size(extra,1), 1);
    for r = 1:size(extra,1)
        fails(r) = sum(extra(r,7) * QuaPar.evalConic(extra(r,1:6), Q) <= 0);
    end
    [~, o] = sort(fails, 'descend');
    for r = o(1:min(3, end)).'
        if ~liveCell([C; extra(r, [1:6 8])], [rows; size(C,1) + 1, extra(r,7)])
            tf = false;  return
        end
    end
    k = size(extra,1);
    tf = liveCell([C; extra(:, [1:6 8])], [rows; size(C,1) + (1:k).', extra(:,7)]);
end

function extra = quadOver(quad, quadApp, j)
% The conic rows saying quadratic j is at least every other applicable one.
    extra = zeros(0,8);
    for i = quadApp(:).'
        if i == j, continue, end
        [c, sgn] = diffConic(quad(j), quad(i));
        if isempty(c), continue, end
        extra(end+1,:) = [c, sgn, 0]; %#ok<AGROW>
    end
end

% --------------------------------------------------------------------------------------------------
function [c, sgn] = diffConic(fw, fj)
% The curve {v_w = v_j}, with the sign that means "w is the larger", as `envelopeSubdivision` spells it.
    num = ratQ.chkB(fw.num * fj.den - fj.num * fw.den, abs(fw.num) * abs(fj.den) + abs(fj.num) * abs(fw.den), 'face difference');
    raw = [num(5), 2*num(6), num(7), 2*num(8), 2*num(9), 2*num(10)];
    if all(raw == 0) || any(num(1:4) ~= 0), c = [];  sgn = 0;  return, end
    nz = find(raw ~= 0, 1);
    c = ratQ.conic(raw);
    sgn = sign(raw(nz));
end

% --------------------------------------------------------------------------------------------------
function a = emptyAlg()
    a = struct('FP', {}, 'FB', {}, 'num', {}, 'den', {}, 'side', {});
end

function [C, ends] = addCornerCrossings(C, corners)
% An IRRATIONAL corner's value is one root of its FP, the other root being its conjugate corner's
% plane; the two planes cross where E1(x) = 0 in cornerFaceAlg's valueForm (K z = E0 + E1 sigma),
% a RATIONAL line. As a curve of the arrangement it makes the corner's branch index constant on
% every cell, which is what lets a corner carry the value with one FB per cell. `ends` pairs each
% such line's key with its corner, the face it matters for (curveOwners).
    ends = cell(0, 2);
    for j = 1:numel(corners)
        c = corners{j};
        if ~isfield(c, 'sigma') || c.sigma.root == 0, continue, end    % rational corner: one root
        n = c.src.n;  d = c.src.d;
        N = [n(5) n(6); n(6) n(7)];  m = [n(8); n(9)];
        al = c.sigma.quad(1);  be = c.sigma.quad(2);
        s0n = c.sigma.s0n;  w = c.sigma.w;  na = c.sigma.na;
        T1 = ratQ.chkB(2*al*(s0n.'*N*w) - be*(w.'*N*w) + 2*na*al*(m.'*w), ...
                       2*abs(al)*(abs(s0n).'*abs(N)*abs(w)) + abs(be)*(abs(w).'*abs(N)*abs(w)) ...
                       + 2*na*abs(al)*(abs(m).'*abs(w)), 'corner crossing');
        E1 = ratQ.chk([2*al*na*d*w(1), 2*al*na*d*w(2), -T1], 'corner crossing');
        if all(E1(1:2) == 0), continue, end
        key = [ratQ.conic([0 0 0 E1]), 0];
        if isempty(find(all(C == key, 2), 1)), C(end+1,:) = key; end %#ok<AGROW>
        ends(end+1,:) = {key, j}; %#ok<AGROW>
    end
end

function [fc, ok] = cornerValueFace(c, pt)
% The corner as an implicit face on this cell: its FP, and the index of the root that is its own
% plane <s,x> - g(s) at the cell's point (constant on the cell: addCornerCrossings).
    n = c.src.n;  d = c.src.d;  s = c.s(:);
    z = pt(:).'*s - ((s.'*[n(5) n(6); n(6) n(7)]*s)/2 + [n(8) n(9)]*s + n(10)) / d;
    [~, info] = algBranch(c.FP, pt, -1);
    rr = info.allReal{1};
    fc = [];  ok = false;
    if isempty(rr), return, end
    [gap, j] = min(abs(rr - z));
    if gap > 1e-6 * max(1, abs(z)), return, end
    fc = struct('FP', c.FP, 'FB', j - numel(rr) - 1);
    ok = true;
end

function [irr, nDrop] = dropBelow(irr, C, rows, P, quad, quadApp)
% The irrational faces CERTIFIED below another one on the whole cell (algFaceBelow), dropped
% (2026-11-02). Each is cut against the same quadratics, so a face below another owns a subset of
% what that one owns and the max is unchanged. On erratum Example 2 a pair's cell carried 6-11
% irrational faces, one of them ever the maximum (300 samples, every cell), and the overlay of the
% pairs multiplied them. The candidates for "above" are read at the cell's points, most often the
% largest first; reading them is a choice, the certificate is the decision.
    nDrop = 0;
    % THE SAME FACE TWICE first -- equal FP and FB, from two bridges or corners: one is enough.
    % Exact (ratQ.polyEqual); on Example 2 the remaining faces of a cell came in such pairs.
    keep = true(1, numel(irr));
    for i = 2:numel(irr)
        for j = 1:i-1
            if keep(j) && irr{i}.FB == irr{j}.FB && isequal(size3(irr{i}.FP), size3(irr{j}.FP)) && ...
               ratQ.polyEqual(irr{i}.FP, irr{j}.FP)
                keep(i) = false;  nDrop = nDrop + 1;  break
            end
        end
    end
    irr = irr(keep);
    n = numel(irr);
    if n < 2, return, end
    V = -inf(size(P,1), n);
    for i = 1:n
        v = algBranch(ratQ.bigToDouble(irr{i}.FP), P, irr{i}.FB);
        v(~isfinite(v)) = -inf;  V(:,i) = v;
    end
    [~, top] = max(V, [], 2);
    [~, order] = sort(accumarray(top, 1, [n 1]), 'descend');
    above = order(1:min(2, n)).';
    [box, inBox, factors, corners] = regionBox(C, rows, zeros(0,8));
    if isempty(box), return, end
    % EACH FACE'S OWN ROWS (2026-11-02): face i owns x only where z_i beats every applicable
    % quadratic (winnerFaces' rows), so z_i <= max z_j is needed there only -- and the z_j above it
    % then beats them too, so owns x. On erratum Example 2 the faces a cell keeps beside its winner
    % meet it ON the cell's edge, where a box straddling the edge was never decided.
    own = cell(1, n);
    for i = 1:n
        own{i} = struct('FP', {}, 'FB', {}, 'num', {}, 'den', {}, 'side', {});
        if nargin < 6, continue, end
        for j = quadApp(:).'
            own{i}(end+1) = struct('FP', irr{i}.FP, 'FB', irr{i}.FB, 'num', quad(j).num, ...
                                   'den', quad(j).den, 'side', -1);
        end
    end
    keep = true(1, n);
    for i = 1:n
        for j = above
            if j == i || ~keep(j), continue, end
            if all(V(:,i) <= V(:,j) + 1e-9 * max(1, abs(V(:,j)))) && algFaceBelow(irr{i}, irr{j}, box, inBox, 5, factors, corners, own{i})
                keep(i) = false;  nDrop = nDrop + 1;  break
            end
        end
    end
    % BELOW THE MAX OF THE OTHERS, and below no single one (2026-11-01): a face that is the maximum
    % at no sample is tried against every face still kept, box by box. Pointwise z_i <= max_j z_j on
    % the cell, and the z_j attaining it owns x whenever z_i does (same quadratics), so the max is
    % unchanged. Faces winning nowhere are tried fewest-near-the-top first.
    wins = accumarray(top, 1, [n 1]).';
    for i = find(keep & wins == 0)
        others = find(keep & (1:n) ~= i);
        if isempty(others), break, end
        hi = cellfun(@(f) struct('FP', {f.FP}, 'FB', f.FB), irr(others), 'UniformOutput', false);  hi = [hi{:}];
        if algFaceBelow(struct('FP', {irr{i}.FP}, 'FB', irr{i}.FB), hi, box, inBox, 5, factors, corners, own{i})
            keep(i) = false;  nDrop = nDrop + 1;
        end
    end
    irr = irr(keep);
end

function sz = size3(P)
    if ratQ.isBigPoly(P), sz = P.sz; else, sz = size(P); end
    sz = [sz, ones(1, 3 - numel(sz))];
end

function ok = dominantIrrational(irr, C, rows, pt, quad, quadApp)
% Do two irrational faces ever BOTH beat every quadratic on this cell? Only there do their pieces
% overlap. Read on a fixed grid over the cell: a diagnostic for the subdivision's topology, never for its
% values (those are a max, and right either way).
    [~, okS, P] = cellSamplePoint(C, rows, straightRows(C, rows), 1e3);
    if ~okS, P = pt; end
    q = -inf(size(P,1), 1);
    for j = quadApp(:).', q = max(q, ratQ.evalFace(quad(j).num, quad(j).den, P)); end
    above = zeros(size(P,1), 1);
    for f = irr
        v = algBranch(f{1}.FP, P, f{1}.FB);
        above = above + (isfinite(v) & v > q + 1e-12 * max(1, abs(q)));
    end
    ok = ~any(above > 1);
end

function H = straightRows(C, rows)
% The cell's straight rows as normalised half-planes [a b c] >= 0.
    H = zeros(0,3);
    for t = 1:size(rows,1)
        c = C(rows(t,1), 1:6);  br = C(rows(t,1), 7);  sd = rows(t,2);
        if br == 0 && all(c(1:3) == 0)
            h = sd * c(4:6);  nr = hypot(h(1), h(2));
            if nr > 0, h = h / nr; end
            H(end+1,:) = h; %#ok<AGROW>
        end
    end
end

function [env, grpQ, grpB] = combinePairs(envs)
% Several pair envelopes as ONE face list: rational corners made faces per pair, every cell-pair
% index offset by the pair's group (pairGroup reads it back), and each face tagged with its group.
    env = struct('quad', struct('num', {}, 'den', {}, 'con', {}, 'pair', {}), 'bridge', {{}});
    grpQ = zeros(1, 0);  grpB = zeros(1, 0);
    for g = 1:numel(envs)
        e = envelopeRationalCorners(envs{g});
        off = (g - 1) * pairStride();
        for k = 1:numel(e.quad)
            q = e.quad(k);
            env.quad(end+1) = struct('num', q.num, 'den', q.den, 'con', q.con, 'pair', q.pair + off);
            grpQ(end+1) = g; %#ok<AGROW>
        end
        for k = 1:numel(e.bridge)
            b = e.bridge{k};  b.pairIndex = b.pairIndex + off;
            env.bridge{end+1} = b;
            grpB(end+1) = g; %#ok<AGROW>
        end
    end
end

function n = pairStride()
    n = 1e6;
end


function [out, info] = minOverlay(subs, infos, budget)
% THE MINIMUM OVER PAIRS, from each pair's own subdivision (2026-10-31). Every pair's cells partition
% the domain; a cell either carries faces (the pair is finite there, its value the max over them) or
% none (the pair is +infinity there). The overlay visits every TUPLE of cells, one per pair, that
% meets -- decided by liveCell on the union of their rows -- and on it emits, for every choice of one
% face per covered pair and every covered pair w, the piece "the face picked for w, where it is <= the
% face picked for every other covered pair".
%
% WHY THE VALUE IS EXACT UNDER A MAX-EVALUATOR. A face of pair g owning x is <= E_g(x), and a pair whose
% cell is a hole is +infinity at x. So a piece's value is <= min_g E_g = f**, and the piece whose picks
% ARE the pairs' maxima, with w the minimising pair, owns x with value f**. No condition WITHIN a pair
% is needed -- each pair's subdivision already took care of its max.
%
% THE ACROSS-PAIR CONDITIONS. Two quadratics: a conic row. A quadratic and an irrational face: an
% algebraic row, as within a pair. TWO IRRATIONAL faces: an `.ii` row, value_1 >= value_2, read in
% doubles by the evaluator -- the curve where two implicit values tie is not stored (minMaxFaces
% refused this case by name, `minIrrational`; Example 2 reaches it). A piece is dropped only when its
% conic rows are PROVABLY empty; one whose algebraic rows may be empty is kept, which costs a face of
% no area, never a value.
    nP = numel(subs);
    keys = zeros(0,7);
    cellRows = cell(1, nP);  covered = cell(1, nP);  facesOf = cell(1, nP);
    for g = 1:nP
        cc = infos{g}.cellCon;
        cellRows{g} = cell(1, numel(cc));
        for k = 1:numel(cc)
            [keys, cellRows{g}{k}] = keyRows(keys, cc{k});
        end
        facesOf{g} = arrayfun(@(k) find(infos{g}.faceCell == k), 1:numel(cc), 'UniformOutput', false);
        covered{g} = ~cellfun(@isempty, facesOf{g});
    end

    % ---- the tuples of cells that meet ---------------------------------------------------------
    tup = num2cell((1:numel(cellRows{1})).');
    trow = cellRows{1};
    for g = 2:nP
        nt = {};  nr = {};
        for a = 1:numel(trow)
            for k = 1:numel(cellRows{g})
                [r, clash] = mergeRows(trow{a}, cellRows{g}{k});
                if clash || ~liveCell(keys, r), continue, end
                nt{end+1} = [tup{a}, k]; %#ok<AGROW>
                nr{end+1} = r; %#ok<AGROW>
            end
        end
        tup = nt;  trow = nr;
    end

    % ---- A BUDGET on the pieces (2026-11-01) ------------------------------------------------------
    % Each tuple visits up to (faces per covered pair, multiplied) x (covered pairs) pieces, each with
    % exact tests. Erratum Example 2 bounds at 11567 (14563 before that day's pair pruning) and ran
    % for hours (DECISIONS 2026-11-01); a
    % caller with no time for that gets a refusal by name, and biconjQ the EnvelopeMin, which answers
    % eval. Inf lifts it.
    bound = 0;
    for t = 1:numel(tup)
        G = find(arrayfun(@(g) covered{g}(tup{t}(g)), 1:nP));
        bound = bound + prod(arrayfun(@(g) numel(facesOf{g}{tup{t}(g)}), G)) * numel(G);
    end
    if bound > budget
        error('envelopeSubdivisionAlg:overlaySize', ['the overlay of %d pairs would visit up to %d ' ...
              'pieces over %d tuples, past the budget of %d; envelopeEval answers for the minimum.'], ...
              nP, bound, numel(tup), budget);
    end

    % ---- the pieces of each tuple --------------------------------------------------------------
    out = struct('kind', {}, 'num', {}, 'den', {}, 'FP', {}, 'FB', {}, 'con', {}, 'alg', {}, 'ii', {});
    nUnsure = 0;  nII = 0;  nHole = 0;  nCombo = 0;  nTie = 0;
    for t = 1:numel(tup)
        G = find(arrayfun(@(g) covered{g}(tup{t}(g)), 1:nP));
        if isempty(G), nHole = nHole + 1; continue, end
        base = zeros(0,8);                       % every pair's cell, holes included
        for g = 1:nP, base = [base; infos{g}.cellCon{tup{t}(g)}]; end %#ok<AGROW>
        [~, ~, PX] = cellSamplePoint(keys, trow{t}, straightRows(keys, trow{t}), 1e3, true);
        P = struct('X', PX, 'Z', containers.Map());  % branch values at PX, computed once per face
        % The picks, one face per covered pair, grown pair by pair and pruned as soon as they
        % cannot meet: a common sample point keeps a choice, else the conic rows and then algEmpty
        % decide. Without the pruning the picks multiply -- 1184 faces on the 3-convex-edge
        % triangle, 240 with it, against 138 for the single arrangement it replaced.
        combos = {struct('pick', {{}}, 'con', base, 'alg', emptyAlg())};
        for i = 1:numel(G)
            nc = {};
            for c = 1:numel(combos)
                for k = facesOf{G(i)}{tup{t}(G(i))}
                    fk = subs{G(i)}(k);
                    con = unique([combos{c}.con; fk.con(:, 1:8)], 'rows', 'stable');
                    alg = [combos{c}.alg, fk.alg];
                    if ~overlayLive(keys, trow{t}, con, alg, emptyII(), P), continue, end
                    nc{end+1} = struct('pick', {[combos{c}.pick, {fk}]}, 'con', con, 'alg', alg); %#ok<AGROW>
                end
            end
            combos = nc;
        end
        nCombo = nCombo + numel(combos);
        for c = 1:numel(combos)
            pick = combos{c}.pick;
            for w = 1:numel(G)
                f = pick{w};
                con = combos{c}.con;  alg = combos{c}.alg;  ii = emptyII();
                for i = 1:numel(G)
                    if i == w, continue, end
                    [con, alg, ii] = leRows(f, pick{i}, con, alg, ii);
                end
                con = unique(con, 'rows', 'stable');
                if ~overlayLive(keys, trow{t}, con, alg, ii, P)
                    continue
                end
                if ~witnessed(con, alg, ii, P)
                    if tieEmpty(keys, trow{t}, con, ii, P), nTie = nTie + 1; continue, end
                    nUnsure = nUnsure + 1;
                end
                nII = nII + numel(ii);
                out(end+1) = struct('kind', f.kind, 'num', f.num, 'den', f.den, 'FP', f.FP, ...
                                    'FB', f.FB, 'con', con, 'alg', alg, 'ii', ii); %#ok<AGROW>
            end
        end
    end
    nAlg = sum(arrayfun(@(f) numel(f.alg), out));
    info = struct('cells', numel(tup), 'faces', numel(out), 'algRows', nAlg, 'iiRows', nII, 'bound', bound, ...
                  'refused', nHole, 'unsure', nUnsure, 'tieEmpty', nTie, 'combos', nCombo, ...
                  'pairCells', cellfun(@(x) x.cells, infos), 'pairFaces', cellfun(@(x) x.faces, infos), ...
                  'skipped', sum(cellfun(@(x) x.skipped, infos)));
end

function [keys, rows] = keyRows(keys, con)
% con rows [a b c d e f sign branch] as [index side] into a shared curve table.
    rows = zeros(0,2);
    for t = 1:size(con,1)
        key = con(t, [1:6 8]);
        hit = find(all(keys == key, 2), 1);
        if isempty(hit), keys(end+1,:) = key;  hit = size(keys,1); end %#ok<AGROW>
        rows(end+1,:) = [hit, con(t,7)]; %#ok<AGROW>
    end
end

function [r, clash] = mergeRows(a, b)
    clash = false;  r = a;
    for t = 1:size(b,1)
        s = a(a(:,1) == b(t,1), 2);
        if isempty(s), r(end+1,:) = b(t,:); %#ok<AGROW>
        elseif s(1) ~= b(t,2), clash = true; return
        end
    end
end

function [con, alg, ii] = leRows(f, o, con, alg, ii)
% The rows saying face f <= face o.
    if strcmp(f.kind, 'quad') && strcmp(o.kind, 'quad')
        [c, sgn] = diffConic(o, f);                 % o >= f
        if ~isempty(c), con(end+1,:) = [c, sgn, 0]; end
    elseif strcmp(f.kind, 'quad')                   % quad f <= irrational o
        alg(end+1) = struct('FP', o.FP, 'FB', o.FB, 'num', f.num, 'den', f.den, 'side', -1);
    elseif strcmp(o.kind, 'quad')                   % irrational f <= quad o
        alg(end+1) = struct('FP', f.FP, 'FB', f.FB, 'num', o.num, 'den', o.den, 'side', 1);
    elseif ~(isequal(f.FP, o.FP) && f.FB == o.FB)
        ii(end+1) = struct('FP1', o.FP, 'FB1', o.FB, 'FP2', f.FP, 'FB2', f.FB);
    end
end

function ii = emptyII()
    ii = struct('FP1', {}, 'FB1', {}, 'FP2', {}, 'FB2', {});
end

function tf = overlayLive(keys, rows, con, alg, ii, P)
% A witness among the tuple's sample points spares the exact tests. Otherwise the conic rows decide,
% then algEmpty certifies the algebraic ones on the conic part -- on the rows without a branch,
% a larger region, so an emptiness proved there holds on the piece. Tie rows are tieEmpty's.
% The FINAL pieces get algEmpty too since 2026-11-01: skipping it removed 2 faces of 240 on the
% 3-convex-edge triangle, but on erratum Example 2's tuple 7 it removes 25 of 124 -- the final
% pieces add "pick_w <= pick_i" rows, which the picks' own certificate never saw.
    if witnessed(con, alg, ii, P), tf = true; return, end
    k0 = all(con(:,1:5) == 0, 2);
    if any(con(k0,7) .* con(k0,6) < 0), tf = false; return, end
    [K, r] = keyRows(keys, con(~k0,:));
    tf = liveCell(K, unique([rows; r], 'rows'));
    if ~tf || isempty(alg), return, end
    plain = ~k0 & con(:,8) == 0;
    tf = ~algEmpty(keys, rows, con(plain, 1:7), alg, P.X);
end

function tf = tieEmpty(keys, rows, con, ii, PP)
% Is the piece PROVABLY of no area because one of its TIE rows, value_1 >= value_2, fails on it --
% value_1 <= value_2 certified on the piece's conic part, strictly off a finite union of curves
% (algFaceBelow)? overlayLive leaves tie rows uncertified, and on erratum Example 2's first tuple 13
% of 21 pieces had no witness (DECISIONS 2026-11-02 (second)). A piece so emptied is a curve at
% most, where the closures of the pieces around it carry the value. A sample point of the conic
% part with value_1 > value_2 rules the row out before any certificate is tried.
    tf = false;
    if isempty(ii), return, end
    k0 = all(con(:,1:5) == 0, 2);
    K = con(~k0, :);
    X = PP.X;  M = true(size(X,1), 1);
    for r = 1:size(K,1)
        if K(r,8) ~= 0, v = ratQ.pairLineEval(K(r,1:6), K(r,8), X); else, v = QuaPar.evalConic(K(r,1:6), X); end
        M = M & K(r,7) * v >= -1e-9 * max(1, max(abs(K(r,1:6))));
    end
    box = [];
    for t = 1:numel(ii)
        if any(M)
            z1 = branchAt(PP, ii(t).FP1, ii(t).FB1, M);  z2 = branchAt(PP, ii(t).FP2, ii(t).FB2, M);
            if any(isfinite(z1) & isfinite(z2) & z1 > z2), continue, end
        end
        if isempty(box)
            [box, inBox, factors, corners] = regionBox(keys, rows, K);
            if isempty(box), return, end
        end
        lo = struct('FP', ii(t).FP1, 'FB', ii(t).FB1);  hi = struct('FP', ii(t).FP2, 'FB', ii(t).FB2);
        try
            tf = algFaceBelow(lo, hi, box, inBox, 5, factors, corners);
        catch
            tf = false;
        end
        if tf, return, end
    end
end

function tf = witnessed(con, alg, ii, PP)
% A sample point of the tuple strictly inside every row. The branch values are the costly part --
% an exact fibre per point on Example 2's integers -- so each face's are cached on the tuple's
% samples and computed only at the points still in play.
    P = PP.X;
    if isempty(P), tf = false; return, end
    M = true(size(P,1), 1);
    for r = 1:size(con,1)
        if con(r,8) ~= 0
            v = ratQ.pairLineEval(con(r,1:6), con(r,8), P);
        else
            v = QuaPar.evalConic(con(r,1:6), P);
        end
        M = M & con(r,7) * v > 1e-9 * max(1, max(abs(con(r,1:6))));
        if ~any(M), tf = false; return, end
    end
    % At most K points go on to the branch rows: a witness is an economy, not the test, and a
    % piece none of them witnesses is decided exactly by the caller.
    K = 96;  in = find(M);
    if numel(in) > K, M(:) = false;  M(in(round(linspace(1, numel(in), K)))) = true; end
    for t = 1:numel(alg)
        a = alg(t);
        z = branchAt(PP, a.FP, a.FB, M);  q = ratQ.evalFace(a.num, a.den, P(M,:));
        M(M) = isfinite(z) & a.side * (q - z) > 1e-9;
        if ~any(M), tf = false; return, end
    end
    for t = 1:numel(ii)
        z1 = branchAt(PP, ii(t).FP1, ii(t).FB1, M);  z2 = branchAt(PP, ii(t).FP2, ii(t).FB2, M);
        M(M) = isfinite(z1) & isfinite(z2) & z1 - z2 > 1e-9;
        if ~any(M), tf = false; return, end
    end
    tf = any(M);
end

function z = branchAt(PP, FP, FB, M)
% algBranch(FP, X(M,:), FB), each point computed at most once per tuple -- IN DOUBLES, from the
% nearest-double FP (2026-11-02). A witness only KEEPS a piece, so a value off by rounding can at
% worst keep a piece of no area, never change a value; and the exact fibre a big FP took was 10 ms
% a point. (Until then a big FP also missed the cache: double() of ratQ's big form raises, and the
% fallback recomputed every point exactly -- 10344 exact fibres on one Example 2 tuple.)
    FP = ratQ.bigToDouble(FP);
    w = FP(:);
    key = sprintf('%d|%d|%.17g|%.17g', FB, numel(w), sum(w), sum(w .* (1:numel(w)).'));
    if isKey(PP.Z, key)
        e = PP.Z(key);
        if ~isequal(e.FP, FP), z = algBranch(FP, PP.X(M,:), FB); return, end
    else
        e = struct('FP', FP, 'v', nan(size(PP.X,1), 1), 'done', false(size(PP.X,1), 1));
    end
    need = M & ~e.done;
    if any(need)
        e.v(need) = algBranch(FP, PP.X(need,:), FB);  e.done(need) = true;
        PP.Z(key) = e;
    end
    z = e.v(M);
end

function [corners, bridges] = splitBridge(bl)
    corners = {};  bridges = {};
    for k = 1:numel(bl)
        f = bl{k};
        if strcmp(f.kind, 'corner') && ~isempty(f.con)
            corners{end+1} = f; %#ok<AGROW>
        elseif isfield(f, 'FP') && isfield(f, 'pair') && isfield(f, 'cells')   % cells: its outer bound
            bridges{end+1} = f; %#ok<AGROW>
        else
            error('envelopeSubdivisionAlg:faceKind', ...
                  'a face of kind %s has neither rows nor an implicit value.', f.kind);
        end
    end
end

% --------------------------------------------------------------------------------------------------
function [C, req, pairOf] = conicTable(faces, corners)
    C = zeros(0,7);
    n = numel(faces) + numel(corners);
    req = cell(1, n);  pairOf = cell(1, n);
    for i = 1:n
        if i <= numel(faces)
            con = faces(i).con;  pairOf{i} = faces(i).pair(:).';
        else
            con = corners{i - numel(faces)}.con;  pairOf{i} = corners{i - numel(faces)}.pairIndex;
        end
        r = zeros(size(con,1), 2);  keep = true(size(con,1), 1);
        for t = 1:size(con,1)
            c = con(t,1:6);
            if all(c(1:5) == 0)
                if ~constHolds(c(6), con(t,7)), req{i} = [NaN NaN];  break, end   % see envelopeSubdivision
                keep(t) = false;  continue
            end
            br = 0;
            if size(con,2) >= 8, br = con(t,8); end
            key = [c, br];
            hit = find(all(C == key, 2), 1);
            if isempty(hit), C(end+1,:) = key;  hit = size(C,1); end %#ok<AGROW>
            r(t,:) = [hit, con(t,7)];
        end
        if isequaln(req{i}, [NaN NaN]), continue, end
        req{i} = r(keep,:);
    end
end

function tf = constHolds(c6, side)
    tf = (side > 0 && c6 >= 0) || (side < 0 && c6 <= 0) || (side == 0 && c6 == 0);
end

function app = applicable(req, rows, C, pt)
% A face applies where the cell meets every one of its demands. A demanded curve the cell does not
% record does not cross the cell, so its sign is the one at the cell's interior point.
    app = [];
    for i = 1:numel(req)
        r = req{i};  ok = true;
        for t = 1:size(r,1)
            if isnan(r(t,1)), ok = false; break, end
            s = rows(rows(:,1) == r(t,1), 2);
            if isempty(s), s = sideAt(C, r(t,1), pt); end
            if s(1) ~= r(t,2), ok = false; break, end
        end
        if ok, app(end+1) = i; end %#ok<AGROW>
    end
end

function s = sideAt(C, j, p)
    c = C(j, 1:6);  br = C(j, 7);
    if br ~= 0
        v = ratQ.pairLineEval(c, br, p);
    else
        v = c(1)*p(1)^2 + c(2)*p(1)*p(2) + c(3)*p(2)^2 + c(4)*p(1) + c(5)*p(2) + c(6);
    end
    s = sign(v);
end

function [C, ends] = addBridgeCurves(C, bridges)
% A bridge's outer bound (`cells`) joins the curve table as curves with no demand, as in envelopeSubdivision:
% then no cell straddles the face's boundary and one membership test per cell is sound.
%
% AND THE ENDS THE OUTER BOUND OMITS (2026-10-17). bridgeFaceAlg gives a sweep end a row only when
% grad g_k is invertible (its image of the tie is then a conic). For a RANK-1 grad g_k that end lies
% on a rational LINE -- the image of the whole plane -- and it was simply missing: on the poly2 kink
% a bridge with an EMPTY outer bound was handed a whole cell it owns only part of, 2.7e-4 too high.
% The line joins as a curve; a double-line row's underlying line likewise (splitLineCurves).
% `ends` pairs each sweep-end line's key with its bridge: it matters only inside the outer bound.
    ends = cell(0, 2);
    for b = 1:numel(bridges)
        blocks = bridges{b}.cells;
        if ~iscell(blocks), blocks = {blocks}; end
        for t = 1:numel(blocks)
            con = blocks{t};
            for r = 1:size(con,1)
                c = con(r,1:6);
                if all(c(1:5) == 0), continue, end
                br = 0;
                if size(con,2) >= 8, br = con(r,8); end
                C = addKey(C, [c, br]);
            end
        end
        pr = bridges{b}.pair;
        if ~isfield(pr, 'N1'), continue, end                      % a rational face: exact rows
        for k = 1:2
            if k == 1, Nk = pr.N1;  mk = pr.m1; else, Nk = pr.N2;  mk = pr.m2; end
            if all(Nk(:) == 0) || ratQ.detExact(Nk) ~= 0, continue, end
            u = Nk(:, find(any(Nk ~= 0, 1), 1));                  % the range of a rank-1 Nk
            w = [-u(2); u(1)];
            raw = ratQ.chkB([0 0 0, pr.d*w(1), pr.d*w(2), -(w.'*mk)], ...
                            [0 0 0, abs(pr.d*w(1)), abs(pr.d*w(2)), abs(w).'*abs(mk)], 'sweep end line');
            C = addKey(C, [ratQ.conic(raw), 0]);
            ends(end+1,:) = {[ratQ.conic(raw), 0], b}; %#ok<AGROW>
        end
    end
end

function C = addKey(C, key)
    if isempty(find(all(C == key, 2), 1)), C(end+1,:) = key; end
end

function [C, ends] = splitLineCurves(C)
% A DOUBLE LINE L^2 (a perfect-square conic, branch 0) has one sign on both sides of L, so as a curve
% it splits no cell; its line L does. Added as its own curve; `ends` pairs L's key with the double
% line's index, whose regions L matters for.
    ends = cell(0, 2);
    for j = 1:size(C,1)
        c = C(j,1:6);
        if C(j,7) ~= 0 || all(c(1:3) == 0), continue, end
        L = squareRoot(c);
        if ~isempty(L)
            C = addKey(C, [ratQ.conic([0 0 0 L]), 0]);
            ends(end+1,:) = {[ratQ.conic([0 0 0 L]), 0], j}; %#ok<AGROW>
        end
    end
end

function owners = curveOwners(C, req, nQuad, bridges, cEnds, bEnds, lEnds)
% For each curve, the regions where it can matter, as row blocks [curve side]: for a face's row, the
% face's OTHER rows; for a bridge's outer-bound row, the other rows of that block; for a corner's
% crossing line, the corner's rows; for a bridge's sweep-end line, each block of its outer bound; for
% a double line's line, the double line's own regions. A curve with no recorded owner (a cut, a
% domain row) gets the empty block -- it matters everywhere. Duplicates are dropped.
    n = size(C,1);
    owners = cell(1, n);
    for i = 1:numel(req)
        r = req{i};
        if isempty(r) || any(isnan(r(:))), continue, end
        for t = 1:size(r,1)
            owners{r(t,1)}{end+1} = r(r(:,1) ~= r(t,1), :);
        end
    end
    for b = 1:numel(bridges)
        blocks = bridgeBlocks(C, bridges{b});
        for t = 1:numel(blocks)
            r = blocks{t};
            for u = 1:size(r,1)
                owners{r(u,1)}{end+1} = r(r(:,1) ~= r(u,1), :);
            end
        end
    end
    for e = 1:size(cEnds,1)
        j = find(all(C == cEnds{e,1}, 2), 1);
        r = req{nQuad + cEnds{e,2}};
        if any(isnan(r(:))), continue, end
        owners{j}{end+1} = r(r(:,1) ~= j, :);
    end
    for e = 1:size(bEnds,1)
        j = find(all(C == bEnds{e,1}, 2), 1);
        blocks = bridgeBlocks(C, bridges{bEnds{e,2}});
        for t = 1:numel(blocks), owners{j}{end+1} = blocks{t}(blocks{t}(:,1) ~= j, :); end
    end
    for e = 1:size(lEnds,1)                       % after the faces, so the double line's are there
        j = find(all(C == lEnds{e,1}, 2), 1);
        src = owners{lEnds{e,2}};
        if isempty(src), src = {zeros(0,2)}; end
        owners{j} = [owners{j}, src];
    end
    for j = 1:n
        if isempty(owners{j}), owners{j} = {zeros(0,2)}; continue, end
        keys = cellfun(@(r) mat2str(sortrows(r)), owners{j}, 'UniformOutput', false);
        [~, ia] = unique(keys, 'stable');
        owners{j} = owners{j}(ia);
    end
end

function blocks = bridgeBlocks(C, bridge)
% A bridge's outer bound as row blocks over C. A constant row is dropped: an unsatisfiable one only
% makes the block empty, and keeping it is the safe side.
    raw = bridge.cells;
    if ~iscell(raw), raw = {raw}; end
    blocks = cell(1, numel(raw));
    for t = 1:numel(raw)
        con = raw{t};  r = zeros(0,2);
        for u = 1:size(con,1)
            c = con(u,1:6);
            if all(c(1:5) == 0), continue, end
            br = 0;
            if size(con,2) >= 8, br = con(u,8); end
            r(end+1,:) = [find(all(C == [c, br], 2), 1), con(u,7)]; %#ok<AGROW>
        end
        blocks{t} = r;
    end
end

function tf = relevantHere(C, owners, rows)
% Can some region the curve bounds meet the cell? Decided per region by liveCell on the cell's rows
% with the region's; a region demanding the opposite side of a curve the cell records is empty.
    tf = true;
    for o = 1:numel(owners)
        r = owners{o};
        if isempty(r), return, end
        clash = false;  keep = true(size(r,1), 1);
        for t = 1:size(r,1)
            s = rows(rows(:,1) == r(t,1), 2);
            if isempty(s), continue, end
            if s(1) ~= r(t,2), clash = true; break, end
            keep(t) = false;
        end
        if clash, continue, end
        if ~any(keep) || liveCell(C, [rows; r(keep,:)]), return, end
    end
    tf = false;
end

function L = squareRoot(c)
% [d e f] with (d x + e y + f)^2 a positive multiple of the conic c = [a b c d e f], else [].
    L = [];
    a = c(1);  b = c(2);  cc = c(3);
    if b^2 ~= 4*a*cc, return, end                                 % quadratic part not a square
    g = gcd(gcd(abs(a), abs(cc)), abs(b));
    if g == 0, return, end
    A = a/g;  B = b/g;  Cq = cc/g;
    if A < 0 || Cq < 0, A = -A;  B = -B;  Cq = -Cq;  sg = -1; else, sg = 1; end
    p = round(sqrt(A));  q = round(sqrt(Cq));
    if p^2 ~= A || q^2 ~= Cq, return, end                         % not a square over the integers
    if B < 0, q = -q; end                                         % (p x + q y)^2 has cross term 2pq
    if 2*p*q ~= B, return, end
    % the full conic must be  sg*g*(p x + q y + r)^2 : match the linear and constant terms
    k = sg*g;
    if p ~= 0, r = c(4) / (2*k*p); else, r = c(5) / (2*k*q); end
    if r ~= fix(r), r2 = 2*r; else, r2 = []; end
    if ~isempty(r2)                                               % allow a half-integer r: scale by 2
        Lc = [2*p 2*q r2];  kk = k/4;
    else
        Lc = [p q r];  kk = k;
    end
    want = kk * [Lc(1)^2, 2*Lc(1)*Lc(2), Lc(2)^2, 2*Lc(1)*Lc(3), 2*Lc(2)*Lc(3), Lc(3)^2];
    if all(abs(want - c) == 0), L = Lc; end
end


function con = rowsToCon(C, rows)
    con = zeros(size(rows,1), 8);
    for t = 1:size(rows,1)
        con(t,1:6) = C(rows(t,1), 1:6);
        con(t,7)   = rows(t,2);
        con(t,8)   = C(rows(t,1), 7);
    end
end

% --------------------------------------------------------------------------------------------------
function [pt, ok] = interiorPoint(C, rows)
% A point strictly inside the cell, and a SIMPLE rational one, because the exact side test clears its
% denominators through a degree-4 polynomial. The Chebyshev centre is rounded onto a dyadic grid and
% kept only if it still satisfies every row.
    pt = [0 0];  ok = false;
    H = zeros(0,3);
    for t = 1:size(rows,1)
        c = C(rows(t,1), 1:6);  br = C(rows(t,1), 7);  sd = rows(t,2);
        if br == 0 && all(c(1:3) == 0), H(end+1,:) = sd * c(4:6); end %#ok<AGROW>
    end
    R = 1e3;
    if isempty(H)
        cand = [0 0];
    else
        for i = 1:size(H,1)
            nrm = hypot(H(i,1), H(i,2));
            if nrm > 0, H(i,:) = H(i,:) / nrm; end
        end
        M = [H(:,1:2), -ones(size(H,1),1); 1 0 0; -1 0 0; 0 1 0; 0 -1 0];
        b = [-H(:,3); -R; -R; -R; -R];
        best = -inf;  cand = [0 0];  k = size(M,1);
        for i = 1:k-2
            for j = i+1:k-1
                for l = j+1:k
                    S = M([i j l], :);
                    if abs(det(S)) < 1e-12, continue, end
                    z = S \ b([i j l]);
                    if all(M*z >= b - 1e-9) && z(3) > best, best = z(3);  cand = z(1:2).'; end
                end
            end
        end
        if best <= 0, return, end
    end
    for pass = 1:2
        for g = [8 16 64 256 1024]
            p = round(cand * g) / g;
            if holdsAll(C, rows, p), pt = p;  ok = true;  return, end
        end
        if holdsAll(C, rows, cand), pt = cand;  ok = true;  return, end
        if pass == 1                                  % a curved row excludes the centre
            [cand, found] = cellSamplePoint(C, rows, H, R);
            if ~found, return, end
        end
    end
end

function tf = holdsAll(C, rows, p)
    tf = true;
    for t = 1:size(rows,1)
        c = C(rows(t,1), 1:6);  br = C(rows(t,1), 7);  sd = rows(t,2);
        if br ~= 0
            v = ratQ.pairLineEval(c, br, p);
        else
            v = c(1)*p(1)^2 + c(2)*p(1)*p(2) + c(3)*p(2)^2 + c(4)*p(1) + c(5)*p(2) + c(6);
        end
        if sd * v <= 0, tf = false; return, end
    end
end
