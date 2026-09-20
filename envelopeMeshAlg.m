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
% STILL REFUSED: two BRIDGE faces competing on one cell. Their difference is not of this shape --
% both sides are irrational -- and it has not occurred on any fixture here, so it is named
% (`envelopeMeshAlg:twoBridges`) rather than guessed at.

    if nargin ~= 1
        error('envelopeMeshAlg:nargin', 'expects (env).');
    end
    if ~isstruct(env) || ~all(isfield(env, {'quad', 'bridge'}))
        error('envelopeMeshAlg:shape', 'env must be what envelopeAssemble returns.');
    end

    [corners, bridges] = splitBridge(env.bridge);
    [C, req, pairOf] = conicTable(env.quad, corners);
    nQuad = numel(env.quad);
    allPairs = unique([pairOf, cellfun(@(f) f.pairIndex, bridges)]);

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
    nAlg = 0;  nUncovered = 0;
    for k = 1:numel(cells)
        rows = cells{k};
        [pt, okPt] = interiorPoint(C, rows);
        if ~okPt, nUncovered = nUncovered + 1; continue, end

        app = applicable(req, rows);
        here = pairOf(app);
        bHere = [];
        for b = 1:numel(bridges)
            if bridgeMember(bridges{b}.pair, pt)
                here(end+1) = bridges{b}.pairIndex; %#ok<AGROW>
                bHere(end+1) = b; %#ok<AGROW>
            end
        end
        if ~all(ismember(allPairs, here)), nUncovered = nUncovered + 1; continue, end

        quadApp = app(app <= nQuad);                 % corners never carry the value
        if numel(bHere) > 1
            error('envelopeMeshAlg:twoBridges', ...
                ['two bridge faces own one cell, so the mesh would have to separate two ' ...
                 'IRRATIONAL values. That curve is not of the bridge-quadratic shape and has not ' ...
                 'occurred on any fixture here.']);
        end
        if isempty(quadApp) && isempty(bHere), nUncovered = nUncovered + 1; continue, end

        sub = winnerFaces(env.quad, quadApp, bridges, bHere, rows, C, pt);
        for t = 1:numel(sub)
            nAlg = nAlg + numel(sub(t).alg);
            mesh(end+1) = sub(t); %#ok<AGROW>
        end
    end
    info = struct('cells', numel(cells), 'faces', numel(mesh), 'algRows', nAlg, ...
                  'refused', nUncovered);
end

% ==================================================================================================
function out = winnerFaces(quad, quadApp, bridges, bHere, rows, C, pt)
% Which applicable face is largest, and the rows that say where. A quad-versus-quad boundary is a
% CONIC and joins the ordinary rows; a quad-versus-bridge boundary is ALGEBRAIC and joins `.alg`.
    out = struct('kind', {}, 'num', {}, 'den', {}, 'FP', {}, 'FB', {}, 'con', {}, 'alg', {});
    base = rowsToCon(C, rows);

    if isempty(bHere)
        % No bridge here: the largest quadratic wins, and ties between two of them are conics.
        [win, extra] = quadWinner(quad, quadApp, pt);
        out(end+1) = struct('kind', 'quad', 'num', quad(win).num, 'den', quad(win).den, ...
                            'FP', [], 'FB', 0, 'con', [base; extra], 'alg', emptyAlg());
        return
    end

    f = bridges{bHere};
    if isempty(quadApp)
        out(end+1) = struct('kind', 'bridge', 'num', [], 'den', 0, ...
                            'FP', f.FP, 'FB', f.FB, 'con', base, 'alg', emptyAlg());
        return
    end

    % A bridge AND quadratics: the bridge is the winner exactly where it is above every one of them,
    % and each comparison is one algebraic row. The complementary pieces are the quadratic winners.
    [win, extra] = quadWinner(quad, quadApp, pt);
    algRows = emptyAlg();
    for j = quadApp(:).'
        algRows(end+1) = struct('FP', f.FP, 'FB', f.FB, 'num', quad(j).num, ...
                                'den', quad(j).den, 'side', -1); %#ok<AGROW>
    end
    out(end+1) = struct('kind', 'bridge', 'num', [], 'den', 0, 'FP', f.FP, 'FB', f.FB, ...
                        'con', base, 'alg', algRows);

    % ... and the piece where that one quadratic beats the bridge
    out(end+1) = struct('kind', 'quad', 'num', quad(win).num, 'den', quad(win).den, ...
                        'FP', [], 'FB', 0, 'con', [base; extra], ...
                        'alg', struct('FP', f.FP, 'FB', f.FB, 'num', quad(win).num, ...
                                      'den', quad(win).den, 'side', 1));
end

% --------------------------------------------------------------------------------------------------
function [win, extra] = quadWinner(quad, quadApp, pt)
% The largest quadratic at the cell's interior point, and the conic rows that pin it against the
% others. Reading the winner at ONE point is sound only because the cell is then cut by every
% difference conic, which is what `extra` is.
    win = quadApp(1);  extra = zeros(0,8);
    best = ratQ.evalFace(quad(win).num, quad(win).den, pt);
    for j = quadApp(2:end)
        v = ratQ.evalFace(quad(j).num, quad(j).den, pt);
        if v > best, best = v;  win = j; end
    end
    for j = quadApp(:).'
        if j == win, continue, end
        [c, sgn] = diffConic(quad(win), quad(j));
        if isempty(c), continue, end
        extra(end+1,:) = [c, sgn, 0]; %#ok<AGROW>
    end
end

% --------------------------------------------------------------------------------------------------
function [c, sgn] = diffConic(fw, fj)
% The curve {v_w = v_j}, with the sign that means "w is the larger", as `envelopeMesh` spells it.
    num = ratQ.chk(fw.num * fj.den - fj.num * fw.den, 'face difference');
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

function [corners, bridges] = splitBridge(bl)
    corners = {};  bridges = {};
    for k = 1:numel(bl)
        f = bl{k};
        if strcmp(f.kind, 'corner') && ~isempty(f.con)
            corners{end+1} = f; %#ok<AGROW>
        elseif isfield(f, 'FP') && isfield(f, 'pair')
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
    req = cell(1, n);  pairOf = zeros(1, n);
    for i = 1:n
        if i <= numel(faces)
            con = faces(i).con;  pairOf(i) = faces(i).pair(1);
        else
            con = corners{i - numel(faces)}.con;  pairOf(i) = corners{i - numel(faces)}.pairIndex;
        end
        r = zeros(size(con,1), 2);  keep = true(size(con,1), 1);
        for t = 1:size(con,1)
            c = con(t,1:6);
            if all(c(1:5) == 0)
                if c(6) < 0, keep(:) = false; break, end
                keep(t) = false;  continue
            end
            br = 0;
            if size(con,2) >= 8, br = con(t,8); end
            key = [c, br];
            hit = find(all(C == key, 2), 1);
            if isempty(hit), C(end+1,:) = key;  hit = size(C,1); end %#ok<AGROW>
            r(t,:) = [hit, con(t,7)];
        end
        req{i} = r(keep,:);
    end
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
    for g = [8 16 64 256 1024]
        p = round(cand * g) / g;
        if holdsAll(C, rows, p), pt = p;  ok = true;  return, end
    end
    if holdsAll(C, rows, cand), pt = cand;  ok = true; end
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
