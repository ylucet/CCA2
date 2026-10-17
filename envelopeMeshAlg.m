function [mesh, info] = envelopeMeshAlg(env)
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

    if nargin ~= 1
        error('envelopeMeshAlg:nargin', 'expects (env).');
    end
    if ~isstruct(env) || ~all(isfield(env, {'quad', 'bridge'}))
        error('envelopeMeshAlg:shape', 'env must be what envelopeAssemble returns.');
    end

    env = envelopeRationalCorners(env);          % a rational corner is a value face (2026-10-17)
    [corners, bridges] = splitBridge(env.bridge);
    [C, req, pairOf] = conicTable(env.quad, corners);
    nQuad = numel(env.quad);
    allPairs = unique([pairOf{:}, cellfun(@(f) f.pairIndex, bridges)]);
    C = addCornerCrossings(C, corners);
    % THE SAME CURVE BUDGET AS envelopeMesh, for the same reason: one global arrangement, whose cost
    % grows far faster than its curve count. Without it the poly2 kink at c = -2 (85 curves), once
    % biconjQ sent irrational-face envelopes here first (2026-10-17), ran >70 min instead of being
    % declined. envelopeEval answers for the face list.
    if size(C,1) > 40
        error('envelopeMeshAlg:arrangementSize', ...
            'the arrangement has %d curves (budget 40), and a global arrangement of that size does not finish.', ...
            size(C,1));
    end

    % ---- the arrangement over the CONIC curves ---------------------------------------------------
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

    % ---- one face per cell, with the winner decided and the ties cut algebraically ----------------
    mesh = struct('kind', {}, 'num', {}, 'den', {}, 'FP', {}, 'FB', {}, 'con', {}, 'alg', {});
    nAlg = 0;  nUncovered = 0;  nOverlap = 0;
    for k = 1:numel(cells)
        rows = cells{k};
        [pt, okPt] = interiorPoint(C, rows);
        if ~okPt, nUncovered = nUncovered + 1; continue, end

        app = applicable(req, rows);
        here = [pairOf{app}];
        bHere = [];
        for b = 1:numel(bridges)
            if bridgeMember(bridges{b}.pair, pt)
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

        sub = winnerFaces(env.quad, quadApp, irr, rows, C, pt);
        for t = 1:numel(sub)
            nAlg = nAlg + numel(sub(t).alg);
            mesh(end+1) = sub(t); %#ok<AGROW>
        end
    end
    info = struct('cells', numel(cells), 'faces', numel(mesh), 'algRows', nAlg, ...
                  'refused', nUncovered, 'overlaps', nOverlap);
end

% ==================================================================================================
function out = winnerFaces(quad, quadApp, irr, rows, C, pt)
% Which applicable face is largest, and the rows that say where. A quad-versus-quad boundary is a
% CONIC and joins the ordinary rows; a quad-versus-bridge boundary is ALGEBRAIC and joins `.alg`.
    out = struct('kind', {}, 'num', {}, 'den', {}, 'FP', {}, 'FB', {}, 'con', {}, 'alg', {});
    base = rowsToCon(C, rows);

    if isempty(irr)
        % No bridge here: the largest quadratic wins, and ties between two of them are conics --
        % one piece PER quadratic (2026-10-17): the winner read at the cell's one point owns only
        % part of the cell, and the rest was left to no face, a hole.
        for j = quadApp(:).'
            out(end+1) = struct('kind', 'quad', 'num', quad(j).num, 'den', quad(j).den, ...
                                'FP', [], 'FB', 0, 'con', [base; quadOver(quad, quadApp, j)], ...
                                'alg', emptyAlg()); %#ok<AGROW>
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
        algRows = emptyAlg();
        for f = irr
            algRows(end+1) = struct('FP', f{1}.FP, 'FB', f{1}.FB, 'num', quad(j).num, ...
                                    'den', quad(j).den, 'side', 1); %#ok<AGROW>
        end
        out(end+1) = struct('kind', 'quad', 'num', quad(j).num, 'den', quad(j).den, ...
                            'FP', [], 'FB', 0, 'con', [base; quadOver(quad, quadApp, j)], ...
                            'alg', algRows); %#ok<AGROW>
    end
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

function app = applicable(req, rows)
    app = [];
    for i = 1:numel(req)
        r = req{i};  ok = true;
        for t = 1:size(r,1)
            s = rows(rows(:,1) == r(t,1), 2);
            if isempty(s) || s(1) ~= r(t,2), ok = false; break, end
        end
        if ok, app(end+1) = i; end %#ok<AGROW>
    end
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
