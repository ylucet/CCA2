function env = biconjAtoms(f)
% biconjAtoms  The convex envelope of a bounded QuaPol as the MINIMUM over its ATOMS' pairwise
%   envelopes and their triple tangent planes -- the route for non-convex pieces that curve UP
%   along an edge, which biconjQ's lower-hull path cannot take.
%
% [input]  f   : QuaPol, exact (fN/fD), bounded (no rays), operable.
% [output] env : EnvelopeMin -- `env.eval(X)` is co f at the points X.
%
% ==================================================================================================
% WHY ATOMS. Karmarkar-Lucet erratum, Theorem 1 (D:\AI\67191168115b2b8bc24afaab\erratum.tex): on a
% bounded piece, sup_x <s,x> - f_i(x) is attained at a VERTEX, at an interior point of an edge
% along which f is strictly convex, or -- only when the piece is convex -- in its interior. So a
% NON-convex piece has exactly the conjugate of
%
%       min( its vertices, each with its value ;  its convex edges, each with f restricted to it ),
%
% and the two functions have the same envelope. Each of those -- a point, a segment carrying a
% convex 1-D quadratic, or a whole CONVEX piece kept as it is -- is a convex function with a
% POLYHEDRAL conjugate (conjQ's KKT active set). That is exactly what `envelopeAssemble` accepts,
% and its refusal of curved cells (`envelopeCellPairs:curvedCell`) is what a non-convex piece's own
% conjugate runs into: the vertex-against-edge tie inside one piece is a parabola.
%
% WHY A MINIMUM. With atoms h_1..h_M, co f = co(min_j h_j). Every co(min(h_a, h_b)) is >= co f, and
% so is every chord: a plane through three points of the atoms' graphs, restricted to their
% triangle. By Caratheodory (dimension 2) co f(x) is attained by at most three points, from at most
% three atoms -- so co f is the MINIMUM over the pairwise envelopes and the triple chords. No
% candidate can undershoot: each one is a convex combination of values of f. That makes dropping a
% candidate unsafe (the answer could only rise), so any refusal inside propagates.
%
%   PAIRS    envelopeAssemble(conjQ(h_a), conjQ(h_b)) -- exact faces, evaluated by envelopeEval.
%   TRIPLES  the common tie s0 of three conjugate FORMS (vertex affine, edge rank-1, convex-piece
%            interior), i.e. a plane tangent to all three, kept when each contact point m(s0) lies
%            in its own atom; the face is <s0,x> - g(s0) on the contact triangle. These are the
%            erratum's Case A polygons (e.g. conv{G1,G2,G3} in its Example 2). s0 is generally
%            irrational (two conics meeting), so it is computed in DOUBLES and refined by Newton;
%            this is a face's DATA, not a mesh decision, and like every face here it is evaluated
%            in doubles.
%
% WHAT IS NOT HERE: rays (the unbounded case needs the erratum's Case E), and the face list is not
% meshed into one AlgCon -- `EnvelopeMin.eval` is the function. Verified against the erratum's
% closed forms in erratumExamplesTest.
% ==================================================================================================
    if ~isa(f, 'QuaPol'), error('biconjAtoms:input', 'expects a QuaPol; got %s.', class(f)); end
    if any(f.E(:,3) == 0)
        error('biconjAtoms:unbounded', ['a piece is unbounded; the atom route covers bounded ' ...
              'domains only (an unbounded one needs the erratum''s Case E).']);
    end
    [atoms, forms] = collectAtoms(f);
    if numel(atoms) < 2
        error('biconjAtoms:oneAtom', ['the input reduces to %d atom(s); a single atom is convex, ' ...
              'so biconjQ''s convex path is the one to take.'], numel(atoms));
    end

    g = cellfun(@(a) conjQ(a.q), atoms, 'UniformOutput', false);
    pairs = {};
    for a = 1:numel(atoms)
        for b = a+1:numel(atoms)
            pairs{end+1} = envelopeAssemble(g{a}, g{b}); %#ok<AGROW>
        end
    end
    tris = tripleFaces(forms);
    env = EnvelopeMin(pairs, tris, atoms);
end

% ==================================================================================================
function [atoms, forms] = collectAtoms(f)
% Atoms (QuaPols) and FORMS (closed-form conjugate pieces for the triples), deduplicated.
    atoms = {};
    forms = struct('kind', {}, 'A', {}, 'b', {}, 'c', {}, 'u', {}, 'd', {}, 'al', {}, 'be', {}, ...
                   'H', {}, 'L', {}, 'poly', {});
    cvx = {};                                    % polygons of the convex-piece atoms
    edgeAt = zeros(0,4);                         % [u w] of each edge atom, for dedup
    vtx = zeros(0,2);  vval = zeros(0,1);        % every vertex that gets a form
    for k = 1:f.nf
        [n, d] = f.faceQ(k);
        H = [n(5) n(6); n(6) n(7)] / d;  L = [n(8); n(9)] / d;  k0 = n(10) / d;
        [poly, eidx] = piecePolygon(f, k);
        isCvx = n(5) >= 0 && n(7) >= 0 && n(5)*n(7) - n(6)^2 >= 0;   % Sylvester on integers
        q = @(X) 0.5*sum((X*H).*X, 2) + X*L + k0;
        for i = 1:size(poly,1)
            [vtx, vval] = addVertex(vtx, vval, poly(i,:), q(poly(i,:)));
        end
        if isCvx
            cvx{end+1} = poly; %#ok<AGROW>
            atoms{end+1} = struct('kind', 'piece', 'q', singlePiece(f, k, eidx)); %#ok<AGROW>
            if n(5)*n(7) - n(6)^2 > 0
                Hi = inv(H);
                forms(end+1) = mkForm('piece', Hi, -Hi*L, 0.5*L.'*Hi*L - k0, [], [], [], [], ...
                                      H, L, poly); %#ok<AGROW>
            end
        end
        for i = 1:size(poly,1)
            u = poly(i,:);  w = poly(mod(i, size(poly,1)) + 1, :);  dv = w - u;
            al = 0.5 * dv*H*dv.';
            if al <= 0, continue, end
            be = (u*H + L.') * dv.';  ga = q(u);
            forms(end+1) = mkForm('edge', dv.'*dv/(2*al), u.' - be*dv.'/(2*al), ...
                                  -ga + be^2/(4*al), u, dv, al, be, [], [], []); %#ok<AGROW>
            if isCvx, continue, end              % the convex piece's atom already holds it
            key = sortrows([u; w]);  key = [key(1,:), key(2,:)];
            if any(all(abs(edgeAt - key) < 1e-12, 2)), continue, end
            edgeAt(end+1,:) = key; %#ok<AGROW>
            atoms{end+1} = struct('kind', 'edge', 'u', u, 'w', w, ...
                'q', QuaPol([u; w], [1 2 1], f.f(k,:), [0 0])); %#ok<AGROW>
        end
    end
    for i = 1:size(vtx,1)
        forms(end+1) = mkForm('vertex', zeros(2), vtx(i,:).', -vval(i), vtx(i,:), [], [], [], ...
                              [], [], []); %#ok<AGROW>
        p = vtx(i,:);
        covered = any(cellfun(@(P) inPoly(p, P), cvx));
        for j = 1:size(edgeAt,1)
            covered = covered || norm(p - edgeAt(j,1:2)) < 1e-12 || norm(p - edgeAt(j,3:4)) < 1e-12;
        end
        if ~covered
            atoms{end+1} = struct('kind', 'vertex', 'u', p, 'q', ...
                QuaPol(p, zeros(0,3), [0 0 0 0 0 0 0 0 0 vval(i)], zeros(0,2))); %#ok<AGROW>
        end
    end
    % an EDGE atom lying in a convex piece is held by that piece already
    keep = true(1, numel(atoms));
    for i = 1:numel(atoms)
        if strcmp(atoms{i}.kind, 'edge')
            keep(i) = ~any(cellfun(@(P) inPoly(atoms{i}.u, P) && inPoly(atoms{i}.w, P), cvx));
        end
    end
    atoms = atoms(keep);
end

function fm = mkForm(kind, A, b, c, u, d, al, be, H, L, poly)
    fm = struct('kind', kind, 'A', A, 'b', b, 'c', c, 'u', u, 'd', d, 'al', al, 'be', be, ...
                'H', H, 'L', L, 'poly', poly);
end

function [vtx, vval] = addVertex(vtx, vval, p, v)
    if ~isempty(vtx) && any(all(abs(vtx - p) < 1e-12, 2)), return, end
    vtx(end+1,:) = p;  vval(end+1,1) = v;
end

function [poly, eidx] = piecePolygon(f, k)
% Face k's vertices in boundary order (counter-clockwise), and its edge indices.
    eidx = find(any(f.F == k, 2));
    E = f.E(eidx, 1:2);
    for j = 1:numel(eidx)                    % orient every edge with face k on its LEFT
        if f.F(eidx(j),1) ~= k, E(j,:) = E(j,[2 1]); end
    end
    order = E(1,1);  nxt = E(1,2);
    while nxt ~= order(1)
        order(end+1) = nxt; %#ok<AGROW>
        nxt = E(E(:,1) == nxt, 2);
        if numel(nxt) ~= 1, error('biconjAtoms:mesh', 'face %d is not a simple polygon.', k); end
    end
    poly = f.V(order, :);
end

function q = singlePiece(f, k, eidx)
    vs = unique(reshape(f.E(eidx,1:2), [], 1));
    idx = zeros(max(vs),1);  idx(vs) = 1:numel(vs);
    E = f.E(eidx,:);  E(:,1:2) = idx(E(:,1:2));
    F = double(f.F(eidx,:) == k);
    q = QuaPol(f.V(vs,:), E, f.f(k,:), F);
end

function tf = inPoly(p, P)
    n = size(P,1);  c = zeros(n,1);
    for i = 1:n
        a = P(i,:);  b = P(mod(i,n)+1,:);
        c(i) = (b(1)-a(1))*(p(2)-a(2)) - (b(2)-a(2))*(p(1)-a(1));
    end
    tf = all(c >= -1e-9 * max(1, max(abs(P(:))))^2);
end

% ==================================================================================================
function tris = tripleFaces(forms)
% For every triple of forms, every common tie s0 whose three contact points lie in their atoms.
    tris = struct('T', {}, 's', {}, 'c', {});
    n = numel(forms);
    for a = 1:n
        for b = a+1:n
            for c = b+1:n
                F3 = forms([a b c]);
                for s0 = commonTies(F3)
                    M = zeros(3,2);  ok = true;
                    for t = 1:3
                        [M(t,:), ok1] = contact(F3(t), s0);
                        ok = ok && ok1;
                    end
                    if ~ok, continue, end
                    ar = (M(2,1)-M(1,1))*(M(3,2)-M(1,2)) - (M(2,2)-M(1,2))*(M(3,1)-M(1,1));
                    if abs(ar) <= 1e-9 * max(1, max(abs(M(:))))^2, continue, end
                    tris(end+1) = struct('T', M, 's', s0.', 'c', gval(F3(1), s0)); %#ok<AGROW>
                end
            end
        end
    end
end

function v = gval(fm, s)
    v = 0.5 * s.'*fm.A*s + fm.b.'*s + fm.c;
end

function [m, ok] = contact(fm, s)
% The maximiser of <s,x> - f over the form's atom, and whether it lies in it.
    switch fm.kind
        case 'vertex'
            m = fm.u;  ok = true;
        case 'edge'
            t = (s.'*fm.d.' - fm.be) / (2*fm.al);
            m = fm.u + t*fm.d;  ok = t >= -1e-9 && t <= 1 + 1e-9;
        case 'piece'
            m = (fm.H \ (s - fm.L)).';  ok = inPoly(m, fm.poly);
    end
end

function S = commonTies(F3)
% Real points where g_1 = g_2 = g_3, as columns.
    c1 = conicOf(F3(1), F3(2));  c2 = conicOf(F3(1), F3(3));
    S = zeros(2,0);
    if all(abs(c1) < 1e-14) || all(abs(c2) < 1e-14), return, end
    for swap = [false true]
        if swap, p1 = c1([3 2 1 5 4 6]); p2 = c2([3 2 1 5 4 6]); else, p1 = c1; p2 = c2; end
        cand = conicMeetDouble(p1, p2);
        if swap, cand = cand([2 1], :); end
        S = [S, cand]; %#ok<AGROW>
        if ~isempty(cand), break, end
    end
    % Newton on (g1-g2, g1-g3), then dedupe
    keep = true(1, size(S,2));
    for i = 1:size(S,2)
        s = S(:,i);
        for it = 1:30
            r = [conicAt(c1, s); conicAt(c2, s)];
            J = [conicGrad(c1, s).'; conicGrad(c2, s).'];
            if rcond(J) < 1e-14, break, end
            ds = -J \ r;  s = s + ds;
            if norm(ds) <= 1e-15 * max(1, norm(s)), break, end
        end
        S(:,i) = s;
        sc = max(1, norm(s))^2 * max(1, max(abs([c1 c2])));
        keep(i) = all(isfinite(s)) && abs(conicAt(c1, s)) < 1e-8*sc && abs(conicAt(c2, s)) < 1e-8*sc;
        for j = 1:i-1
            if keep(j) && norm(S(:,j) - s) < 1e-9 * max(1, norm(s)), keep(i) = false; end
        end
    end
    S = S(:, keep);
end

function c = conicOf(fa, fb)
% g_a - g_b as [a b c d e f]: a s1^2 + b s1 s2 + c s2^2 + d s1 + e s2 + f.
    A = fa.A - fb.A;  b = fa.b - fb.b;
    c = [A(1,1)/2, A(1,2), A(2,2)/2, b(1), b(2), fa.c - fb.c];
end

function v = conicAt(c, s)
    v = c(1)*s(1)^2 + c(2)*s(1)*s(2) + c(3)*s(2)^2 + c(4)*s(1) + c(5)*s(2) + c(6);
end

function g = conicGrad(c, s)
    g = [2*c(1)*s(1) + c(2)*s(2) + c(4); c(2)*s(1) + 2*c(3)*s(2) + c(5)];
end

function S = conicMeetDouble(c1, c2)
% Intersections of two conics through the resultant in s2, in doubles.
    A1 = c1(3);  B1 = [c1(2) c1(5)];  C1 = [c1(1) c1(4) c1(6)];
    A2 = c2(3);  B2 = [c2(2) c2(5)];  C2 = [c2(1) c2(4) c2(6)];
    if A1 == 0 && A2 == 0
        R = padd(conv(B1, C2), -conv(B2, C1));
    else
        X = padd(A1*C2, -A2*C1);  Y = padd(A1*B2, -A2*B1);  Z = padd(conv(B1, C2), -conv(B2, C1));
        R = padd(conv(X, X), -conv(Y, Z));
    end
    R = trimLead(R);
    S = zeros(2,0);
    if isempty(R) || all(R == 0), return, end
    r = roots(R);
    r = real(r(abs(imag(r)) <= 1e-7 * max(1, abs(r))));
    for s1 = r(:).'
        a = A1;  bb = polyval(B1, s1);  cc = polyval(C1, s1);
        if abs(a) > 0, s2c = roots([a bb cc]); elseif bb ~= 0, s2c = -cc/bb; else, s2c = []; end
        a2 = A2;  b2 = polyval(B2, s1);  cc2 = polyval(C2, s1);
        if isempty(s2c) || (abs(a) == 0 && bb == 0)
            if abs(a2) > 0, s2c = roots([a2 b2 cc2]); elseif b2 ~= 0, s2c = -cc2/b2; end
        end
        s2c = real(s2c(abs(imag(s2c)) <= 1e-7 * max(1, abs(s2c))));
        for s2 = s2c(:).'
            S(:,end+1) = [s1; s2]; %#ok<AGROW>
        end
    end
end

function p = padd(a, b)
    n = max(numel(a), numel(b));
    p = [zeros(1, n-numel(a)), a] + [zeros(1, n-numel(b)), b];
end

function p = trimLead(p)
    i = find(abs(p) > 1e-12 * max(1, max(abs(p))), 1);
    if isempty(i), p = []; else, p = p(i:end); end
end
