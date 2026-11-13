function g = conjAffinePLQ(faces)
% conjAffinePLQ  The conjugate of a PIECEWISE-AFFINE function, in one construction, with no
%   triangulation, no symbolic call and no Step 3 fold.
%
% objective: `TODO.md` G2. An affine face over an UNBOUNDED polygon conjugates to a SUPPORT
%   function -- +inf off a cone -- and `maxQuaPar.assertFullDomain` refuses any operand that is not
%   finite everywhere, so `max(0,x,y)` and its family have always gone to the symbolic Case C.
%   `TODO.md` prices the way out: for an input whose pieces are ALL affine, f* is the max of
%   finitely many affine functions restricted to one polyhedron, so it can be built directly and
%   never enter `maxQuaPar` at all. This is that.
%
% [input]  faces : struct array, one entry per face of dom f, with
%                    .W       n x 2 vertices in CCW order (interior on the LEFT of the walk)
%                    .dFirst  1 x 2 recession direction leaving W(1,:),   or [] if bounded
%                    .dLast   1 x 2 recession direction leaving W(end,:), or [] if bounded
%                    .a       1 x 2 gradient, .b scalar: the face carries f(x) = <a,x> + b
%                  -- the same per-face description `conjConvexPolygon` takes, produced by
%                  `conjCPLQ`'s own faceBoundary.
% [output] g     : QuaPol. Every face function is AFFINE and every edge is straight, because a
%                  max of affine functions is affine on each cell of its own subdivision.
%
% THE MATHEMATICS, and it is three lines. On face i, f is <a_i,x> + b_i, so
%
%       f*(s) = max_i  sup_{x in P_i} <s - a_i, x> - b_i
%             = max_i  [ sigma_{P_i}(s - a_i) - b_i ]
%
% and for a polyhedron P_i = conv(W_i) + cone(D_i) the support function is
%
%       sigma_{P_i}(t) = max_j <t, w_j>   when <t,d> <= 0 for every d in D_i,   +inf otherwise.
%
% Two consequences, and together they ARE the algorithm:
%
%   * DOMAIN. f*(s) is finite exactly where EVERY face's support function is, since a max with
%     +inf in it is +inf. So dom f* is the intersection over i and over d in D_i of the half-planes
%     <d, s> <= <d, a_i>. A bounded face contributes nothing. This is why the answer has a
%     BOUNDED domain even though the input does not: `max(0,x,y)`'s three cones contribute
%     s1 >= 0, s2 >= 0 and s1 + s2 <= 1, whose intersection is the unit simplex.
%   * VALUE. On that domain f* is the maximum of the affine functions
%           L_{i,j}(s) = <s, w_j^i> - <a_i, w_j^i> - b_i,
%     one per (face, vertex of that face). No other candidate can attain the sup: the sup of a
%     linear functional over a polyhedron on which it is bounded is attained at a vertex.
%
% THE SUBDIVISION is then the ordinary one for a max of affine functions: cell k is where
% candidate k beats every other, which is an intersection of half-planes,
%       cell_k = { s : <p_k - p_j, s> >= q_j - q_k  for every j }  intersected with the domain,
% computed here by successive half-plane clipping. Duplicate candidates (two faces sharing a
% vertex give the same L, which is the normal case at every interior edge of dom f) are removed
% first, so those constraints are never degenerate.
%
% SCOPE, stated rather than discovered. dom f* must come out BOUNDED; a piecewise-affine f whose
% conjugate has an unbounded domain is refused by name. That is not the interesting case for G2 --
% dom f* is bounded exactly when f grows in every direction, which is what the family this exists
% for does -- and it keeps the clipping below finite. Every face must be affine; a caller with a
% quadratic face has `conjConvexPolygon` or the triangulation route instead.

    if isempty(faces)
        error('PLQ:conjAffinePLQ:notImplemented', 'no faces to conjugate.');
    end

    % ---- the domain: one half-plane per recession direction ------------------------------------
    cons = zeros(0,3);                      % rows [n1 n2 c] meaning n*s' <= c
    for i = 1:numel(faces)
        for d = {faces(i).dFirst, faces(i).dLast}
            if isempty(d{1}), continue, end
            nd = d{1}(:).' / norm(d{1});
            cons(end+1,:) = [nd, nd * faces(i).a(:)]; %#ok<AGROW>
        end
    end

    % ---- the candidate affine functions --------------------------------------------------------
    P = zeros(0,2); Q = zeros(0,1);
    for i = 1:numel(faces)
        W = faces(i).W;
        for j = 1:size(W,1)
            w = W(j,:);
            P(end+1,:) = w;                                        %#ok<AGROW>
            Q(end+1,1) = -(faces(i).a(:).' * w(:)) - faces(i).b;    %#ok<AGROW>
        end
    end
    [P, Q] = dedupCandidates(P, Q);

    % ---- the domain polygon --------------------------------------------------------------------
    dom = boxPolygon(P, cons);
    for r = 1:size(cons,1)
        dom = clipHalfPlane(dom, cons(r,1:2), cons(r,3));
        if size(dom,1) < 3
            error('PLQ:conjAffinePLQ:emptyDomain', ...
                ['dom f* is empty: some face''s support function is +inf everywhere, i.e. f is ' ...
                 'unbounded below along a recession direction of its own domain.']);
        end
    end
    if touchesBox(dom, P, cons)
        error('PLQ:conjAffinePLQ:unboundedDual', ...
            ['dom f* is not bounded, which this construction does not cover. It is bounded ' ...
             'exactly when f grows in every direction; see this file''s SCOPE note.']);
    end

    % ---- one cell per surviving candidate ------------------------------------------------------
    cells = {}; frows = zeros(0,6);
    for k = 1:size(P,1)
        c = dom;
        for j = 1:size(P,1)
            if j == k, continue, end
            % keep  <p_j - p_k, s> <= q_k - q_j
            c = clipHalfPlane(c, P(j,:) - P(k,:), Q(k) - Q(j));
            if size(c,1) < 3, break, end
        end
        if size(c,1) < 3, continue, end
        if polyArea(c) <= 1e-12 * max(1, max(abs(c(:))))^2, continue, end
        cells{end+1} = c; %#ok<AGROW>
        frows(end+1,:) = [0 0 0, P(k,1), P(k,2), Q(k)]; %#ok<AGROW>
    end
    if isempty(cells)
        error('PLQ:conjAffinePLQ:notImplemented', ...
            'the subdivision came out empty, which should not happen on a non-empty domain.');
    end

    g = assembleAffineCells(cells, frows);
end

% ================================================================================================
function [P, Q] = dedupCandidates(P, Q)
% Two faces meeting along an edge share its vertices, and a shared vertex gives the SAME affine
% function from both sides (that is what continuity of f across the edge means). Keeping both
% would make the pairwise constraint <p_j - p_k, s> >= q_k - q_j the trivial 0 >= 0 and split one
% cell into two coincident ones.
    keep = true(1, size(P,1));
    for i = 1:size(P,1)
        if ~keep(i), continue, end
        for j = i+1:size(P,1)
            if ~keep(j), continue, end
            if norm(P(i,:) - P(j,:)) <= 1e-12*(1+norm(P(i,:))) && ...
               abs(Q(i) - Q(j)) <= 1e-12*(1+abs(Q(i)))
                keep(j) = false;
            end
        end
    end
    P = P(keep,:); Q = Q(keep);
end

function W = boxPolygon(P, cons)
% A CCW box large enough that clipping it by `cons` cannot cut a genuine bounded cell, used as the
% starting polygon for the half-plane intersections. Its size is read off the data rather than
% fixed, and `touchesBox` afterwards is what turns "I used a box" into a checked statement.
    R = 10 * max([1; abs(P(:)); abs(cons(:))]);
    W = [-R -R; R -R; R R; -R R];
end

function tf = touchesBox(W, P, cons)
    R = 10 * max([1; abs(P(:)); abs(cons(:))]);
    tf = any(abs(abs(W(:)) - R) <= 1e-9*R);
end

function W = clipHalfPlane(W, n, c)
% Sutherland-Hodgman clip of the convex polygon W by {s : n*s' <= c}. Returns fewer than 3 rows
% when nothing survives.
    if isempty(W), return, end
    nn = norm(n);
    if nn <= 1e-14, return, end                 % a degenerate constraint clips nothing
    n = n / nn; c = c / nn;
    m = size(W,1);
    v = W * n(:) - c;
    tol = 1e-12 * (1 + max(abs(W(:))));
    out = zeros(0,2);
    for i = 1:m
        j = mod(i, m) + 1;
        if v(i) <= tol, out(end+1,:) = W(i,:); end %#ok<AGROW>
        if (v(i) > tol && v(j) < -tol) || (v(i) < -tol && v(j) > tol)
            t = v(i) / (v(i) - v(j));
            out(end+1,:) = W(i,:) + t * (W(j,:) - W(i,:)); %#ok<AGROW>
        end
    end
    W = dedupRows(out);
end

function W = dedupRows(W)
    if size(W,1) < 2, return, end
    keep = true(1, size(W,1));
    for i = 1:size(W,1)
        j = mod(i, size(W,1)) + 1;
        if keep(i) && norm(W(i,:) - W(j,:)) <= 1e-11*(1+max(abs(W(:)))), keep(j) = false; end
    end
    W = W(keep,:);
end

function a = polyArea(W)
    x = W(:,1); y = W(:,2);
    a = abs(0.5 * sum(x .* y([2:end 1]) - x([2:end 1]) .* y));
end

function g = assembleAffineCells(cells, frows)
% Glue the cells into one QuaPol: dedup the vertices by coordinate, emit each cell's edges, and
% pair the ones that coincide. Safe to do by coordinate here, unlike maxQuaPar's assembly: these
% cells all come from clipping ONE box by exact half-planes, so a shared corner is computed by the
% same arithmetic from both sides and the near-degenerate scales that forced maxQuaPar into
% half-edge provenance do not arise.
    V = zeros(0,2); nf = numel(cells); idx = cell(1,nf);
    tolV = 0;
    for k = 1:nf, tolV = max(tolV, max(abs(cells{k}(:)))); end
    tolV = 1e-9 * (1 + tolV);
    for k = 1:nf
        W = ensureCCWlocal(cells{k});
        id = zeros(1, size(W,1));
        for r = 1:size(W,1)
            hit = 0;
            for q = 1:size(V,1)
                if norm(V(q,:) - W(r,:)) <= tolV, hit = q; break, end
            end
            if hit == 0, V(end+1,:) = W(r,:); hit = size(V,1); end %#ok<AGROW>
            id(r) = hit;
        end
        idx{k} = id;
    end
    E = zeros(0,3); F = zeros(0,2);
    key = containers.Map('KeyType','char','ValueType','double');
    for k = 1:nf
        id = idx{k};
        for r = 1:numel(id)
            a = id(r); b = id(mod(r, numel(id)) + 1);
            kk = sprintf('%d-%d', min(a,b), max(a,b));
            if isKey(key, kk)
                e = key(kk);
                if F(e,2) == 0, F(e,2) = k; end
            else
                E(end+1,:) = [a b 1]; %#ok<AGROW>
                F(end+1,:) = [k 0];   %#ok<AGROW>
                key(kk) = size(E,1);
            end
        end
    end
    g = QuaPol(V, E, frows, F);
end

function W = ensureCCWlocal(W)
    x = W(:,1); y = W(:,2);
    if 0.5 * sum(x .* y([2:end 1]) - x([2:end 1]) .* y) < 0, W = flipud(W); end
end
