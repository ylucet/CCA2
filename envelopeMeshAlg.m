function [mesh, info] = envelopeMeshAlg(env, dom)
% envelopeMeshAlg  The envelope as a mesh whose EDGES may be algebraic -- the arrangement that
%   `envelopeMesh` refuses when a bridge face competes.
%
% objective: a bridge face carries the envelope's value, and where it competes with a quadratic face
%   the boundary is `FP(x, y, q(x,y)) = 0`, of degree four. `envelopeMesh` raises
%   `envelopeMesh:bridgeCompetes` there because a conic subdivision cannot hold that edge. This
%   builds the same arrangement with such edges allowed, so the GEOMETRY can be verified before the
%   output TYPE is grown to store it.
%
% [input]  env  : as `envelopeAssemble` returns it, ideally pruned
% [output] mesh : struct array, one entry per face, with
%                   .kind  'quad' | 'bridge'
%                   .num, .den   the value, when kind is 'quad'
%                   .FP, .FB     the value, when kind is 'bridge'
%                   .con   conic rows [a b c d e f sign branch], as everywhere else
%                   .alg   algebraic rows, a struct array with .FP .FB .num .den .side, meaning
%                          side * (q(x) - value_FB(x)) >= 0
%          info : .cells, .faces, .algRows, .refused
%
% ==================================================================================================
% WHAT IS NEW HERE, AND WHAT IS DELIBERATELY THE SAME
%
% THE SAME: the arrangement is built from the CONIC curves exactly as before -- the quadratic faces'
% rows, the corner wedges, the bridges' outer bounds -- and every cell's emptiness is decided by
% `cellInterior`. Corner faces are carried for COVERAGE only, because a corner never carries the
% value (measured, DECISIONS.md 2026-09-19), and a pair with no face owning a cell makes the
% envelope +infinity there, which is a hole in the mesh.
%
% NEW: where two applicable faces cannot be separated by a conic, the cell is split by an ALGEBRAIC
% row instead of being refused. Such a row is `(FP, FB, num, den, side)` -- the curve is the rational
% `bridgeQuadCut` polynomial, and `side` picks one of its two branches, because that polynomial names
% the PAIR and not one of them (DECISIONS.md 2026-09-20). `algCurveSide` decides it exactly at a
% rational point; `meshAlgEval` reads it in doubles, which is what evaluation is allowed to do.
%
% SEVERAL IRRATIONAL faces on one cell (bridges, and irrational CORNERS, which carry value too --
% 2026-10-17) are not cut against each other: that curve is not of this shape. Each is cut against
% the quadratics only, so where two of them both beat every quadratic their pieces OVERLAP; the
% value is still right because every evaluator takes the max over the faces owning a point.
% `info.overlaps` counts such cells (a fixed-grid read). Until 2026-10-17 this was refused as
% `envelopeMeshAlg:twoBridges`; erratum Example 1 reaches it with two corners from Q(sqrt 26) and
% Q(sqrt 6), and 18 of its cells overlap.

    if nargin < 1 || nargin > 2
        error('envelopeMeshAlg:nargin', 'expects (env) or (env, dom).');
    end
    if nargin < 2, dom = zeros(0,3); end
    % SEVERAL PAIRS (2026-10-18): a cell array of pair envelopes is the MINIMUM over them -- an
    % EnvelopeMin with no triangle face. Their faces go into one arrangement, tagged by pair, and each
    % cell is resolved as a min of maxes (minMaxFaces). One pair is exactly the path below.
    if iscell(env), envs = env; else, envs = {env}; end
    for g = 1:numel(envs)
        if ~isstruct(envs{g}) || ~all(isfield(envs{g}, {'quad', 'bridge'}))
            error('envelopeMeshAlg:shape', 'env must be what envelopeAssemble returns.');
        end
    end
    [env, grpQ, grpB] = combinePairs(envs);
    nGrp = numel(envs);
    [corners, bridges] = splitBridge(env.bridge);
    grpCorner = grpB(cellfun(@(f) strcmp(f.kind, 'corner') && ~isempty(f.con), env.bridge));
    grpBridge = grpB(cellfun(@(f) ~(strcmp(f.kind, 'corner') && ~isempty(f.con)), env.bridge));
    [C, req, pairOf] = conicTable(env.quad, corners);
    nQuad = numel(env.quad);
    allPairs = unique([pairOf{:}, cellfun(@(f) f.pairIndex, bridges)]);
    C = addCornerCrossings(C, corners);
    C = addBridgeCurves(C, bridges);   % without them a bridge's membership can change inside a cell
    C = splitLineCurves(C);

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
    maxCells = 3000;
    for j = 1:size(C,1)
        next = {};
        for k = 1:numel(cells)
            if any(cells{k}(:,1) == j), next{end+1} = cells{k}; continue, end %#ok<AGROW>
            live = [liveCell(C, [cells{k}; j, 1]), liveCell(C, [cells{k}; j, -1])];
            if all(live)
                next{end+1} = [cells{k}; j, 1];  next{end+1} = [cells{k}; j, -1]; %#ok<AGROW>
            elseif any(live)
                next{end+1} = cells{k}; %#ok<AGROW>     % the curve misses this cell
            end
        end
        cells = next;
        if numel(cells) > maxCells
            error('envelopeMeshAlg:arrangementSize', ...
                'the arrangement passed %d cells at curve %d of %d; envelopeEval answers for the face list.', ...
                maxCells, j, size(C,1));
        end
    end

    % ---- one face per cell, with the winner decided and the ties cut algebraically ----------------
    mesh = struct('kind', {}, 'num', {}, 'den', {}, 'FP', {}, 'FB', {}, 'con', {}, 'alg', {});
    nAlg = 0;  nUncovered = 0;  nOverlap = 0;  nNoPoint = 0;
    for k = 1:numel(cells)
        rows = cells{k};
        [pt, okPt] = interiorPoint(C, rows);
        if ~okPt, nUncovered = nUncovered + 1;  nNoPoint = nNoPoint + 1; continue, end

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
        % sparse rows algCrvMeshTest found a face vertex ending other than two arcs (2026-10-17). A
        % curve the cell does not record has one sign on it -- read at the cell's points.
        rowsOut = fullRows(C, rows, Q);      % emitted; the exact tests below keep the sparse rows
        for b = 1:numel(bridges)
            if sum(bridgeMember(bridges{b}.pair, Q)) * 2 > size(Q,1)
                here(end+1) = bridges{b}.pairIndex; %#ok<AGROW>
                bHere(end+1) = b; %#ok<AGROW>
            end
        end
        if nGrp > 1
            % EACH PAIR is finite here only where all of ITS cell pairs are covered; the minimum is
            % over the pairs that are.
            S = struct('quad', {}, 'irr', {});
            for g = 1:nGrp
                pg = allPairs(pairGroup(allPairs) == g);
                if ~all(ismember(pg, here)), continue, end
                qa = app(app <= nQuad);  qa = qa(grpQ(qa) == g);
                ir = cell(1, 0);
                for b = bHere(grpBridge(bHere) == g), ir{end+1} = struct('FP', bridges{b}.FP, 'FB', bridges{b}.FB); end %#ok<AGROW>
                for j = app(app > nQuad)
                    if grpCorner(j - nQuad) ~= g, continue, end
                    [fc, okc] = cornerValueFace(corners{j - nQuad}, pt);
                    if okc, ir{end+1} = fc; end %#ok<AGROW>
                end
                if isempty(qa) && isempty(ir), continue, end
                S(end+1) = struct('quad', qa, 'irr', {ir}); %#ok<AGROW>
            end
            if isempty(S), nUncovered = nUncovered + 1; continue, end
            sub = minMaxFaces(env.quad, S, rows, C, Pw, rowsOut);
            for t = 1:numel(sub)
                nAlg = nAlg + numel(sub(t).alg);
                mesh(end+1) = sub(t); %#ok<AGROW>
            end
            continue
        end
        if ~all(ismember(allPairs, here)), nUncovered = nUncovered + 1; continue, end

        quadApp = app(app <= nQuad);
        % THE IRRATIONAL CANDIDATES: the bridges owning the cell, and -- since 2026-10-17 -- the
        % CORNERS that apply to it. "A corner never carries the value" was measured on one fixture
        % (DECISIONS.md 2026-09-19) and is false in general: on every one-pair atom envelope a
        % corner is the maximum in slivers near a domain vertex, by up to 8e-6, and ignoring it left
        % the mesh that much too LOW there.
        irr = cell(1, 0);
        for b = bHere, irr{end+1} = struct('FP', bridges{b}.FP, 'FB', bridges{b}.FB); end %#ok<AGROW>
        for j = app(app > nQuad)
            [fc, okc] = cornerValueFace(corners{j - nQuad}, pt);
            if okc, irr{end+1} = fc; end %#ok<AGROW>
        end
        if isempty(quadApp) && isempty(irr), nUncovered = nUncovered + 1; continue, end
        if numel(irr) > 1 && ~dominantIrrational(irr, C, rows, pt, env.quad, quadApp)
            nOverlap = nOverlap + 1;
        end

        sub = winnerFaces(env.quad, quadApp, irr, rows, C, pt, Pw, rowsOut);
        for t = 1:numel(sub)
            nAlg = nAlg + numel(sub(t).alg);
            mesh(end+1) = sub(t); %#ok<AGROW>
        end
    end
    info = struct('cells', numel(cells), 'faces', numel(mesh), 'algRows', nAlg, ...
                  'refused', nUncovered, 'noPoint', nNoPoint, 'overlaps', nOverlap);
end

% ==================================================================================================
function out = winnerFaces(quad, quadApp, irr, rows, C, pt, Q, rowsOut) %#ok<INUSD>
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
        if algEmpty(C, rows, extra, algRows, Q), continue, end
        out(end+1) = struct('kind', 'quad', 'num', quad(j).num, 'den', quad(j).den, ...
                            'FP', [], 'FB', 0, 'con', [base; extra], 'alg', algRows); %#ok<AGROW>
    end
end

function rows = fullRows(C, rows, P)
% Every curve's side on the cell: recorded ones as recorded, the rest by the sign at the cell's points
% (the majority of the nonzero signs; a curve tangent at one point does not decide it).
    have = false(size(C,1), 1);  have(rows(:,1)) = true;
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
        z = algBranch(a.FP, Qk, a.FB);  q = ratQ.evalFace(a.num, a.den, Qk);
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
    [box, inBox] = regionBox(C, rows, conRows);
    if ~isempty(box)
        try
            ok = algRowFails(a, box, inBox, 6);
        catch
            ok = false;
        end
    end
    if ~isempty(key), cache(key) = ok; end
end

function [box, inBox] = regionBox(C, rows, K)
% A box around the cell (its straight rows, within the meshers' bounding box) and a test rejecting
% sub-boxes that provably miss the cell's or the piece's conic rows (branch rows are not used).
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
    if size(V,1) < 3, box = []; inBox = []; return, end
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
% The curve {v_w = v_j}, with the sign that means "w is the larger", as `envelopeMesh` spells it.
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

function C = addCornerCrossings(C, corners)
% An IRRATIONAL corner's value is one root of its FP, the other root being its conjugate corner's
% plane; the two planes cross where E1(x) = 0 in cornerFaceAlg's valueForm (K z = E0 + E1 sigma),
% a RATIONAL line. As a curve of the arrangement it makes the corner's branch index constant on
% every cell, which is what lets a corner carry the value with one FB per cell.
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

function ok = dominantIrrational(irr, C, rows, pt, quad, quadApp)
% Do two irrational faces ever BOTH beat every quadratic on this cell? Only there do their pieces
% overlap. Read on a fixed grid over the cell: a diagnostic for the mesh's topology, never for its
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

function g = pairGroup(p)
    g = floor((p - 1) / pairStride()) + 1;
end

function out = minMaxFaces(quad, S, rows, C, Q, rowsOut)
% THE MINIMUM OVER PAIRS OF THE MAXIMUM OVER EACH PAIR'S FACES, on one cell (2026-10-18). A piece is a
% choice of one face per pair and of the pair w whose choice is the smallest; its value is that face.
%
% WHY THE VALUE IS EXACT UNDER A MAX-EVALUATOR. Any face of pair g that applies is <= E_g, the pair's
% value; a piece requires pick_w <= pick_g for every g, so its value is <= min_g E_g = f**. And the
% piece whose picks ARE the pairs' maxima, with w the minimising pair, owns the point with value f**.
% So the max over the pieces owning a point is f** -- the across-pair conditions are what correctness
% needs. The within-pair ones ("pick_g is the largest of g") only prune pieces, so where they would
% compare two irrational faces of one pair they are simply left out. Two irrational faces from
% DIFFERENT pairs must be ordered: if their order is not the same at every sample point of the cell
% that is refused by name, and biconjQ keeps the minimum unmeshed.
    out = struct('kind', {}, 'num', {}, 'den', {}, 'FP', {}, 'FB', {}, 'con', {}, 'alg', {});
    base = rowsToCon(C, rowsOut);
    cand = cell(1, numel(S));  whole = cell(1, numel(S));
    for g = 1:numel(S)
        cg = {};
        for j = S(g).quad(:).', cg{end+1} = struct('kind', 'quad', 'j', j, 'FP', [], 'FB', 0); end %#ok<AGROW>
        for f = S(g).irr, cg{end+1} = struct('kind', 'irr', 'j', 0, 'FP', f{1}.FP, 'FB', f{1}.FB); end %#ok<AGROW>
        % VIABLE PICKS ONLY (2026-10-18): a face can be pair g's pick only where it is the largest
        % of pair g's faces, and that is decided ONCE per face -- exactly, on its conic rows --
        % instead of once per combination of picks (2440 faces on the 3-convex-edge triangle,
        % most of them combinations no point realises).
        keep = true(1, numel(cg));
        for o = 1:numel(cg)
            if ~strcmp(cg{o}.kind, 'quad'), continue, end
            r = zeros(0,8);  a = emptyAlg();
            for o2 = 1:numel(cg)
                if o2 == o || ~strcmp(cg{o2}.kind, 'quad'), continue, end
                [r, a] = geRows(quad, cg{o}, cg{o2}, r, a, Q);
            end
            keep(o) = pieceLive(C, rows, r, Q);
        end
        cand{g} = cg(keep);
        if isempty(cand{g}), cand{g} = cg; end         % something is the max: never leave a hole
        whole{g} = cg;                                 % every face of pair g, for the rows
    end
    sizes = cellfun(@numel, cand);
    for idx = 1:prod(sizes)
        sub = cell(1, numel(S));  [sub{:}] = ind2sub([sizes 1], idx);
        pick = cellfun(@(c, i) c{i}, cand, sub, 'UniformOutput', false);
        for w = 1:numel(S)
            conRows = zeros(0,8);  algRows = emptyAlg();  ok = true;
            for g = 1:numel(S)                         % pick{g} is the largest of pair g
                for o = 1:numel(whole{g})
                    other = whole{g}{o};
                    if isequal(other, pick{g}), continue, end
                    if ~strcmp(pick{g}.kind, 'quad') && ~strcmp(other.kind, 'quad')
                        continue                       % within a pair: pruning only, see above
                    end
                    [conRows, algRows, ok] = geRows(quad, pick{g}, other, conRows, algRows, Q);
                    if ~ok, break, end
                end
                if ~ok, break, end
            end
            if ~ok, continue, end
            deferred = {};
            for g = setdiff(1:numel(S), w)             % ... and pair w's is the smallest
                if ~strcmp(pick{g}.kind, 'quad') && ~strcmp(pick{w}.kind, 'quad')
                    deferred{end+1} = {pick{g}, pick{w}}; %#ok<AGROW>
                    continue
                end
                [conRows, algRows, ok] = geRows(quad, pick{g}, pick{w}, conRows, algRows, Q);
                if ~ok, break, end
            end
            if ~ok, continue, end
            % TWO IRRATIONAL PICKS from different pairs: their order matters only where the piece's
            % other conditions hold, so it is read on the cell's points that satisfy them.
            if ~isempty(deferred)
                M = piecePoints(quad, Q, conRows, algRows);
                if ~any(M), continue, end              % no point of the cell: nothing to own
                for d = 1:numel(deferred)
                    a = deferred{d}{1};  b = deferred{d}{2};    % need a >= b
                    va = algBranch(a.FP, Q(M,:), a.FB);  vb = algBranch(b.FP, Q(M,:), b.FB);
                    tol = 1e-12 * max(1, abs(vb));
                    if all(va >= vb - tol), continue, end
                    if all(va <= vb + tol), ok = false; break, end
                    error('envelopeMeshAlg:minIrrational', ...
                        'two irrational faces from different pairs cross inside one piece of a minimum over pairs.');
                end
                if ~ok, continue, end
            end
            if ~pieceLive(C, rows, conRows, Q), continue, end
            if algEmpty(C, rows, conRows, algRows, Q), continue, end
            f = pick{w};
            if strcmp(f.kind, 'quad')
                out(end+1) = struct('kind', 'quad', 'num', quad(f.j).num, 'den', quad(f.j).den, ...
                                    'FP', [], 'FB', 0, 'con', [base; conRows], 'alg', algRows); %#ok<AGROW>
            else
                out(end+1) = struct('kind', 'bridge', 'num', [], 'den', 0, 'FP', f.FP, 'FB', f.FB, ...
                                    'con', [base; conRows], 'alg', algRows); %#ok<AGROW>
            end
        end
    end
end

function M = piecePoints(quad, Q, conRows, algRows) %#ok<INUSL>
% Which of the cell's points satisfy a piece's conic and algebraic rows.
    M = true(size(Q,1), 1);
    for r = 1:size(conRows,1)
        M = M & conRows(r,7) * QuaPar.evalConic(conRows(r,1:6), Q) >= 0;
    end
    for t = 1:numel(algRows)
        a = algRows(t);
        z = algBranch(a.FP, Q, a.FB);  q = ratQ.evalFace(a.num, a.den, Q);
        M = M & isfinite(z) & a.side * (q - z) >= -1e-9;
    end
end

function [conRows, algRows, ok] = geRows(quad, a, b, conRows, algRows, Q)
% The rows saying face a >= face b. ok = false when they PROVABLY cannot hold on the cell's points
% (two irrational faces ordered the wrong way at every sample) -- the piece is then empty.
    ok = true;
    if strcmp(a.kind, 'quad') && strcmp(b.kind, 'quad')
        if a.j == b.j, return, end
        [c, sgn] = diffConic(quad(a.j), quad(b.j));
        if ~isempty(c), conRows(end+1,:) = [c, sgn, 0]; end
    elseif strcmp(a.kind, 'quad')                       % quad a >= irrational b
        algRows(end+1) = struct('FP', b.FP, 'FB', b.FB, 'num', quad(a.j).num, 'den', quad(a.j).den, 'side', 1);
    elseif strcmp(b.kind, 'quad')                       % irrational a >= quad b
        algRows(end+1) = struct('FP', a.FP, 'FB', a.FB, 'num', quad(b.j).num, 'den', quad(b.j).den, 'side', -1);
    else
        if isequal(a.FP, b.FP) && a.FB == b.FB, return, end
        va = algBranch(a.FP, Q, a.FB);  vb = algBranch(b.FP, Q, b.FB);
        tol = 1e-12 * max(1, abs(vb));
        if all(va >= vb - tol), return, end             % a >= b on the whole cell (sampled)
        if all(va <= vb + tol), ok = false; return, end % never: the piece is empty
        error('envelopeMeshAlg:minIrrational', ...
            'two irrational faces from different pairs cross inside one cell of a minimum over pairs.');
    end
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
            error('envelopeMeshAlg:faceKind', ...
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
                if ~constHolds(c(6), con(t,7)), req{i} = [NaN NaN];  break, end   % see envelopeMesh
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

function C = addBridgeCurves(C, bridges)
% A bridge's outer bound (`cells`) joins the curve table as curves with no demand, as in envelopeMesh:
% then no cell straddles the face's boundary and one membership test per cell is sound.
%
% AND THE ENDS THE OUTER BOUND OMITS (2026-10-17). bridgeFaceAlg gives a sweep end a row only when
% grad g_k is invertible (its image of the tie is then a conic). For a RANK-1 grad g_k that end lies
% on a rational LINE -- the image of the whole plane -- and it was simply missing: on the poly2 kink
% a bridge with an EMPTY outer bound was handed a whole cell it owns only part of, 2.7e-4 too high.
% The line joins as a curve; a double-line row's underlying line likewise (splitLineCurves).
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
        end
    end
end

function C = addKey(C, key)
    if isempty(find(all(C == key, 2), 1)), C(end+1,:) = key; end
end

function C = splitLineCurves(C)
% A DOUBLE LINE L^2 (a perfect-square conic, branch 0) has one sign on both sides of L, so as a curve
% it splits no cell; its line L does. Added as its own curve.
    for j = 1:size(C,1)
        c = C(j,1:6);
        if C(j,7) ~= 0 || all(c(1:3) == 0), continue, end
        L = squareRoot(c);
        if ~isempty(L), C = addKey(C, [ratQ.conic([0 0 0 L]), 0]); end
    end
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
