function env = biconjAtoms(f)
% biconjAtoms  The convex envelope of a bounded QuaPol as the MINIMUM over its ATOMS' pairwise
%   envelopes and their triple tangent planes -- the route for non-convex pieces that curve UP
%   along an edge, which biconjQ's lower-hull path cannot take.
%
% [input]  f   : QuaPol, exact (fN/fD), operable; bounded or not.
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
% ONLY THE ACTIVE CANDIDATES (2026-10-16, erratum Step 1): a pair is assembled only where two of
% its conjugate forms tie for f* along an arc, a triangle kept only where its three forms need three
% atoms -- `activeAtomPairs` says why that is enough. Erratum Example 2: 6 pairs to 3, 60 triangles
% to 2, same values, 780 s to 281 s (one run each, shared machine).
%
% RAYS (erratum Case E) are below, at the restriction to D. WHAT IS NOT HERE: the face list is
% not meshed into one AlgCon -- `EnvelopeMin.eval` is the function. Verified against the erratum's closed forms in erratumExamplesTest.
% ==================================================================================================
    if ~isa(f, 'QuaPol'), error('biconjAtoms:input', 'expects a QuaPol; got %s.', class(f)); end
    % AN UNBOUNDED INPUT (erratum Case E). Its rays carry two things. A ray along which f is
    % strictly convex is an ATOM like a convex edge, with t in [0, inf). A ray along which f grows
    % at most linearly is not an atom at all: its sup sits at the base vertex, and what it does is
    % shrink dom f* to a polyhedron D. So every pair is assembled over D -- sup over s in D of
    % <s,x> - max(g_a, g_b) -- which is how a value like erratum Example 3's x1^2, NOT attained by
    % finitely many points, comes out exact: it is the recession part, read off D. Each such
    % candidate is still >= f** (max(g_a, g_b) + I_D <= f*), so the minimum stays an upper bound
    % that cannot undershoot; that it is attained is verified on Example 3, not proved here.
    D = zeros(0,3);
    if any(f.E(:,3) == 0)
        [~, D] = conjQ(f);
        D = unique(D, 'rows');
        if any(all(D(:,1:2) == 0, 2) & D(:,3) < 0)
            error('biconjAtoms:emptyDomain', ['dom f* is empty: f has no affine minorant, so ' ...
                  'co f is -infinity everywhere.']);
        end
    end
    [atoms, forms] = collectAtoms(f);
    if numel(atoms) == 1 && ~isempty(D)
        % ONE ATOM, BUT THE RECESSION STILL MATTERS: (g + I_D)* is not the atom itself. Paired with
        % itself, the assembly sees identical faces and returns exactly that.
        g1 = conjQ(atoms{1}.q);
        env = EnvelopeMin({envelopeAssemble(g1, g1, D)}, rmfield(tripleFaces(forms, D), 'forms'), atoms);
        return
    end
    if numel(atoms) < 2
        error('biconjAtoms:oneAtom', ['the input reduces to %d atom(s); a single atom is convex, ' ...
              'so biconjQ''s convex path is the one to take.'], numel(atoms));
    end

    % ONLY THE ACTIVE PAIRS (erratum Step 1). See activeAtomPairs: a pair is assembled only where
    % two of its conjugate forms tie for f* along an arc, and a triangle is kept only where its
    % three forms need three different atoms.
    tris = tripleFaces(forms, D);
    [ab, tris] = activeAtomPairs(atoms, forms, tris, D);
    g = cell(1, numel(atoms));
    for a = unique(ab(:)).'
        g{a} = conjQ(atoms{a}.q);
    end
    pairs = cell(1, size(ab,1));
    for k = 1:size(ab,1)
        pairs{k} = envelopeAssemble(g{ab(k,1)}, g{ab(k,2)}, D);
    end
    env = EnvelopeMin(pairs, tris, atoms);
end

% ==================================================================================================
function [ab, tris] = activeAtomPairs(atoms, forms, tris, D)
% The atom pairs, and the triangle faces, that the MINIMUM needs -- erratum Step 1.
%
% WHY THIS IS ENOUGH. f* = max over the forms g (vertex, convex edge, convex-piece interior), and at
% a point x where co f is attained, x is in the subdifferential of f* at a maximiser s -- the
% convex hull of the contact points m_g(s) of the forms ACTIVE at s (plus N_D(s) when D bounds
% dom f*). By Caratheodory three contact points suffice. A contact point of a form lies in an atom
% A, and then it is a subgradient of g_A at s (g_A(s') >= <s',m> - f(m) = g_A(s) + <s'-s, m>), so:
%
%   * if two atoms A, B hold the forms involved, x is in the subdifferential of max(g_A, g_B) at s
%     and that pair's envelope is <= <s,x> - f*(s) = co f(x) there: the PAIR (A,B) attains;
%   * if three atoms are needed, the TRIANGLE face of that tie attains.
%
% And only points x of positive measure matter: every candidate is convex and >= co f, so a minimum
% that equals co f on a dense set equals it everywhere on the interior (each candidate is continuous
% on its polyhedral domain, Rockafellar 10.2). So a form pair counts only if the two forms tie for
% the max along an ARC (an isolated tie point is a vertex of f*'s mesh, which is a triangle's job),
% and only if their contact segments sweep an area -- two vertices, or a vertex on an edge's line,
% or two edges on one line, sweep only that line.
%
% THE TIE TEST IS IN DOUBLES, like the triangle faces themselves, and it errs toward KEEPING: a
% near-tie counts as a tie, and a competing form must beat the tie by a margin to exclude it. Dropping
% a pair that was needed would RAISE the minimum, so every doubt keeps the pair.
    nA = numel(atoms);
    holds = formHolders(atoms, forms);
    req = {};                                % each: struct('pairs', [A B; ...], 'any', atoms)
    nF = numel(forms);
    for i = 1:nF
        for j = i+1:nF
            if sweepsNoArea(forms(i), forms(j)), continue, end
            if ~tieArc(i, j, forms, D), continue, end
            req{end+1} = requirement(holds([i j]), nA); %#ok<AGROW>
        end
    end
    % a convex piece's interior form contacts a 2-D set by ITSELF, and with D so does an edge form
    % (m_g(s) + N_D(s) along an edge of D): those need some pair containing their atom
    for i = 1:nF
        if strcmp(forms(i).kind, 'piece') || (~isempty(D) && strcmp(forms(i).kind, 'edge'))
            req{end+1} = struct('pairs', zeros(0,2), 'any', holds{i}); %#ok<AGROW>
        end
    end
    keepT = true(1, numel(tris));
    for t = 1:numel(tris)
        if ~triangleActive(tris(t), forms, D)
            keepT(t) = false;
            continue
        end
        r = requirement(holds(tris(t).forms), nA);
        if ~isempty(r.pairs) || ~isempty(r.any)
            keepT(t) = false;                % two atoms hold all three: the pair attains it
            req{end+1} = r; %#ok<AGROW>
        end
    end
    tris = rmfield(tris(keepT), 'forms');
    if isempty(req)
        ab = nchoosek(1:nA, 2);              % no arc tie found at all: keep everything (safe)
        return
    end
    ab = coverRequirements(req, nA);
end

function r = requirement(H, nA)
% The atom pairs that hold every form in H (a cell of holder lists): one atom holding some forms
% and another holding the rest, or one atom holding them all (then any pair containing it).
    all1 = H{1};
    for k = 2:numel(H), all1 = intersect(all1, H{k}); end
    P = zeros(0,2);
    for A = 1:nA
        for B = A+1:nA
            if all(cellfun(@(h) any(h == A) || any(h == B), H)), P(end+1,:) = [A B]; end %#ok<AGROW>
        end
    end
    r = struct('pairs', P, 'any', all1(:).');
end

function ab = coverRequirements(req, nA)
% Greedy: keep the pair satisfying the most unsatisfied requirements until none is left.
    cand = nchoosek(1:nA, 2);
    sat = false(numel(req), size(cand,1));
    for r = 1:numel(req)
        for c = 1:size(cand,1)
            sat(r,c) = ismember(cand(c,:), req{r}.pairs, 'rows') || any(ismember(cand(c,:), req{r}.any));
        end
    end
    if any(~any(sat, 2))
        error('biconjAtoms:uncovered', 'a tie no atom pair holds -- the form holders are inconsistent.');
    end
    ab = zeros(0,2);  open = true(numel(req),1);
    while any(open)
        [~, c] = max(sum(sat(open,:), 1));
        ab(end+1,:) = cand(c,:); %#ok<AGROW>
        open = open & ~sat(:,c);
    end
    ab = sortrows(ab);
end

function H = formHolders(atoms, forms)
% H{i}: the atoms whose region contains form i's whole contact set -- for those atoms the contact
% point is a subgradient of the atom's conjugate wherever the form is active.
    H = cell(1, numel(forms));
    for i = 1:numel(forms)
        fm = forms(i);
        h = [];
        for a = 1:numel(atoms)
            if atomHolds(atoms{a}, fm), h(end+1) = a; end %#ok<AGROW>
        end
        H{i} = h;
    end
end

function tf = atomHolds(at, fm)
    tol = 1e-12;
    switch fm.kind
        case 'vertex'
            pts = fm.u;  dirs = zeros(0,2);
        case 'edge'
            if isinf(fm.tmax), pts = fm.u;  dirs = fm.d;
            else, pts = [fm.u; fm.u + fm.d];  dirs = zeros(0,2); end
        case 'piece'
            tf = strcmp(at.kind, 'piece') && isequal(at.region, fm.poly);
            return
    end
    switch at.kind
        case 'vertex'
            tf = isempty(dirs) && all(vecnorm(pts - at.u, 2, 2) < tol);
        case 'edge'
            tf = isempty(dirs) && all(arrayfun(@(k) onSegment(pts(k,:), at.u, at.w), 1:size(pts,1)));
        case 'ray'
            dr = at.q.V(2,:) - at.q.V(1,:);          % a ray runs THROUGH its second vertex
            tf = all(arrayfun(@(k) onRay(pts(k,:), at.u, dr), 1:size(pts,1))) && ...
                 all(arrayfun(@(k) parallelSame(dirs(k,:), dr), 1:size(dirs,1)));
        case 'piece'
            P = at.region;
            tf = all(arrayfun(@(k) inPiece(pts(k,:), P), 1:size(pts,1)));
            if tf && ~isempty(dirs)
                if ~isstruct(P), tf = false; return, end
                R = P.rows;
                tf = all(R(:,1:2) * dirs.' >= -1e-9 * max(1, max(abs(R(:)))) * max(1, norm(dirs)), 'all');
            end
    end
end

function tf = onSegment(p, u, w)
    d = w - u;  L2 = d*d.';
    t = (p - u)*d.' / L2;
    tf = t >= -1e-12 && t <= 1 + 1e-12 && norm(u + t*d - p) <= 1e-12 * max(1, norm(p));
end

function tf = onRay(p, u, d)
    t = (p - u)*d.' / (d*d.');
    tf = t >= -1e-12 && norm(u + t*d - p) <= 1e-12 * max(1, norm(p));
end

function tf = parallelSame(a, b)
    tf = abs(a(1)*b(2) - a(2)*b(1)) <= 1e-12 * norm(a) * norm(b) && a*b.' > 0;
end

function tf = sweepsNoArea(f1, f2)
% Both contact sets on one line: the segments between them sweep only that line.
    k = {f1.kind, f2.kind};
    if any(strcmp(k, 'piece')), tf = false; return, end
    if all(strcmp(k, 'vertex')), tf = true; return, end
    if strcmp(f1.kind, 'vertex'), [f1, f2] = deal(f2, f1); end      % f1 is an edge now
    onLine = @(p) abs(f1.d(1)*(p(2)-f1.u(2)) - f1.d(2)*(p(1)-f1.u(1))) <= ...
                  1e-12 * norm(f1.d) * max(1, norm(p - f1.u));
    if strcmp(f2.kind, 'vertex'), tf = onLine(f2.u);
    else, tf = onLine(f2.u) && onLine(f2.u + f2.d); end
end

% ---- the tie test ----------------------------------------------------------------------------------
function tf = tieArc(i, j, forms, D)
% Do forms i and j tie for the maximum of f* along an arc of positive length? The tie {g_i = g_j}
% is a conic in s; parametrise it (a rational parametrisation through one real point, or lines for a
% degenerate conic), take every breakpoint where some condition changes sign, and test a point
% inside each interval between them.
    Q = conicOf(forms(i), forms(j));
    for par = conicParams(Q)
        bp = [];
        for c = 1:numel(forms)
            R = domRows(forms(c));
            for r = 1:size(R,1), bp = [bp; realRoots(compose(par, [0 0 0 R(r,:)]))]; end %#ok<AGROW>
            if c ~= i && c ~= j
                bp = [bp; realRoots(compose(par, conicOf(forms(c), forms(i))))]; %#ok<AGROW>
            end
        end
        for r = 1:size(D,1), bp = [bp; realRoots(compose(par, [0 0 0 D(r,:)]))]; end %#ok<AGROW>
        if strcmp(par.kind, 'rat'), bp = [bp; realRoots(par.Dn)]; end %#ok<AGROW>
        bp = unique(bp);
        if isempty(bp)
            samp = 0;
        else
            samp = [bp(1) - 1 - abs(bp(1)); (bp(1:end-1) + bp(2:end))/2; bp(end) + 1 + abs(bp(end))];
        end
        for t = samp.'
            [s, ok] = pointAt(par, t);
            if ok && bothActive(s, i, j, forms, D), tf = true; return, end
        end
    end
    tf = false;
end

function tf = bothActive(s, i, j, forms, D)
% s in dom g_i, dom g_j and D (leniently), and no form strictly beats them there.
    tf = inDom(s, domRows(forms(i)), -1) && inDom(s, domRows(forms(j)), -1) && inDom(s, D, -1);
    if ~tf, return, end
    vi = gval(forms(i), s);
    for c = 1:numel(forms)
        if c == i || c == j || ~inDom(s, domRows(forms(c)), +1), continue, end
        if gval(forms(c), s) > vi + 1e-9 * max(1, abs(vi)), tf = false; return, end
    end
end

function tf = triangleActive(tr, forms, D)
% A triangle's tie s0 is a vertex of f*'s mesh only if no other form beats it there.
    s = tr.s(:);
    tf = inDom(s, D, -1);
    v = tr.c;
    for c = setdiff(1:numel(forms), tr.forms)
        if ~tf, return, end
        if inDom(s, domRows(forms(c)), +1) && gval(forms(c), s) > v + 1e-9 * max(1, abs(v))
            tf = false;
        end
    end
end

function tf = inDom(s, R, strictness)
% strictness -1: within tolerance counts as inside; +1: must be inside by the tolerance.
    if isempty(R), tf = true; return, end
    r = R(:,1:2)*s(:) + R(:,3);
    tol = 1e-9 * max(1, max(abs(R(:)))) * max(1, norm(s));
    tf = all(r >= -strictness * tol);
end

function R = domRows(fm)
% dom of a form in s, as rows [a b c] with a s1 + b s2 + c >= 0.
    switch fm.kind
        case 'vertex'
            R = zeros(0,3);
        case 'edge'
            r = [fm.d, -fm.be] / (2*fm.al);          % t(s) >= 0
            R = r;
            if isfinite(fm.tmax), R = [R; -r(1:2), fm.tmax - r(3)]; end
        case 'piece'
            Hi = inv(fm.H);  P = fm.poly;            % m(s) = Hi (s - L) in the piece
            if isstruct(P)
                Rm = P.rows;
            else
                n = size(P,1);  Rm = zeros(n,3);
                for e = 1:n
                    a = P(e,:);  b = P(mod(e,n)+1,:);
                    nr = [-(b(2)-a(2)), b(1)-a(1)];
                    Rm(e,:) = [nr, -nr*a.'];
                end
            end
            R = [Rm(:,1:2)*Hi, Rm(:,3) - Rm(:,1:2)*Hi*fm.L];
    end
end

function P = conicParams(Q)
% Parametrisations of the real conic Q = 0 ([a b c d e f]), as structs: 'rat' (a rational curve
% s = (S1, S2)/Dn in m, through a real point) or 'line' (s = p + t w). Empty or point conics give
% none: they hold no arc.
    P = struct('kind', {}, 'S1', {}, 'S2', {}, 'Dn', {}, 'p', {}, 'w', {});
    sc = max(abs(Q));
    if sc == 0, return, end
    Q = Q / sc;
    A = [Q(1) Q(2)/2; Q(2)/2 Q(3)];  g = Q(4:5).';  f0 = Q(6);
    if norm(A) <= 1e-12
        if norm(g) <= 1e-12, return, end
        P(1) = lineP(-f0*g.'/(g.'*g), [-g(2) g(1)]);
        return
    end
    M = [A g/2; g.'/2 f0];
    if abs(det(M)) > 1e-10 * max(1, norm(M))^3
        p0 = realPoint(A, g, f0);
        if isempty(p0), return, end
        G = 2*A*p0 + g;
        Dn = [A(2,2), 2*A(1,2), A(1,1)];
        N  = [G(2), G(1)];
        P(1) = struct('kind', 'rat', 'S1', padd(p0(1)*Dn, -N), 'S2', padd(p0(2)*Dn, -[N 0]), ...
                      'Dn', Dn, 'p', [], 'w', []);
        return
    end
    dA = det(A);
    if dA > 1e-12, return, end                           % a single point
    if dA < -1e-12                                       % two crossing lines
        c0 = -(A \ g)/2;
        [V, L] = eig(A);  l = diag(L);
        [~, o] = sort(l, 'descend');  V = V(:,o);  l = l(o);
        for sg = [1 -1]
            w = V(:,1)/sqrt(l(1)) + sg*V(:,2)/sqrt(-l(2));
            P(end+1) = lineP(c0.', w.'); %#ok<AGROW>
        end
        return
    end
    [V, L] = eig(A);  [~, k] = max(abs(diag(L)));  l = L(k,k);  v = V(:,k);  vp = [-v(2); v(1)];
    for r = roots([l, g.'*v, f0]).'                      % parallel lines v.s = r
        if abs(imag(r)) > 1e-9 * max(1, abs(r)), continue, end
        P(end+1) = lineP(real(r)*v.', vp.'); %#ok<AGROW>
    end
end

function p = lineP(p0, w)
    p = struct('kind', 'line', 'S1', [], 'S2', [], 'Dn', [], 'p', p0, 'w', w);
end

function p0 = realPoint(A, g, f0)
% One real point of a nondegenerate conic, or [] when it has none.
    p0 = [];
    [V, ~] = eig(A);
    base = {zeros(2,1)};
    if abs(det(A)) > 1e-12, base{end+1} = -(A \ g)/2; end
    dirs = [V, eye(2), null(A)];
    for b = 1:numel(base)
        for k = 1:size(dirs,2)
            w = dirs(:,k);  p = base{b};
            c = [w.'*A*w, (2*A*p + g).'*w, p.'*A*p + g.'*p + f0];
            t = roots(c);
            t = real(t(abs(imag(t)) <= 1e-9 * max(1, abs(t))));
            if ~isempty(t), p0 = p + t(1)*w; return, end
        end
    end
end

function c = compose(par, Q)
% Q = [a b c d e f] along the parametrisation, as a polynomial in the parameter whose real roots
% are where Q changes sign along it (for 'rat', Q times Dn^2, which does not change its sign).
    if strcmp(par.kind, 'line')
        A = [Q(1) Q(2)/2; Q(2)/2 Q(3)];  g = Q(4:5).';  p = par.p(:);  w = par.w(:);
        c = [w.'*A*w, (2*A*p + g).'*w, p.'*A*p + g.'*p + Q(6)];
    else
        S1 = par.S1;  S2 = par.S2;  Dn = par.Dn;
        c = padd(padd(padd(Q(1)*conv(S1,S1), Q(2)*conv(S1,S2)), Q(3)*conv(S2,S2)), ...
                 padd(conv(padd(Q(4)*S1, Q(5)*S2), Dn), Q(6)*conv(Dn,Dn)));
    end
end

function r = realRoots(c)
    c = trimLead(c);
    if numel(c) < 2, r = zeros(0,1); return, end
    r = roots(c);
    r = real(r(abs(imag(r)) <= 1e-6 * max(1, abs(r))));
    r = r(:);
end

function [s, ok] = pointAt(par, t)
    if strcmp(par.kind, 'line')
        s = par.p(:) + t*par.w(:);  ok = true;
    else
        dn = polyval(par.Dn, t);
        ok = abs(dn) > 1e-12 * max(1, max(abs(par.Dn))) * max(1, t^2);
        s = [polyval(par.S1, t); polyval(par.S2, t)] / dn;
    end
end

% ==================================================================================================
function [atoms, forms] = collectAtoms(f)
% Atoms (QuaPols) and FORMS (closed-form conjugate pieces for the triples), deduplicated.
    atoms = {};
    forms = struct('kind', {}, 'A', {}, 'b', {}, 'c', {}, 'u', {}, 'd', {}, 'al', {}, 'be', {}, ...
                   'H', {}, 'L', {}, 'poly', {}, 'tmax', {});
    cvx = {};                                    % the convex-piece atoms' regions (see inPiece)
    edgeAt = zeros(0,4);                         % [u w] of each edge atom, for dedup
    rayAt = zeros(0,4);                          % [u d] of each RAY atom, d primitive
    vtx = zeros(0,2);  vval = zeros(0,1);        % every vertex that gets a form
    for k = 1:f.nf
        [n, d] = f.faceQ(k);
        H = [n(5) n(6); n(6) n(7)] / d;  L = [n(8); n(9)] / d;  k0 = n(10) / d;
        isCvx = n(5) >= 0 && n(7) >= 0 && n(5)*n(7) - n(6)^2 >= 0;   % Sylvester on integers
        q = @(X) 0.5*sum((X*H).*X, 2) + X*L + k0;
        % the piece's edges as (u, dv, isRay): a segment u -> u + dv, or a ray u + t dv, t >= 0
        own = find(any(f.F == k, 2));
        if any(f.E(own,3) == 0)
            ed = zeros(0,5);
            for j = own(:).'
                u = f.V(f.E(j,1),:);  w = f.V(f.E(j,2),:);
                if f.E(j,3) == 0, ed(end+1,:) = [u, w - u, 1]; %#ok<AGROW>  % a ray runs THROUGH V(b)
                else, ed(end+1,:) = [u, w - u, 0]; end %#ok<AGROW>
            end
            % AN UNBOUNDED PIECE HAS NO POLYGON: membership is by its half-planes instead.
            poly = struct('rows', pieceRows(ed));
            eidx = own;
        else
            [poly, eidx] = piecePolygon(f, k);
            m = size(poly,1);
            ed = [poly, poly([2:m 1],:) - poly, zeros(m,1)];
        end
        for i = 1:size(ed,1)
            [vtx, vval] = addVertex(vtx, vval, ed(i,1:2), q(ed(i,1:2)));
            if ~ed(i,5), [vtx, vval] = addVertex(vtx, vval, ed(i,1:2) + ed(i,3:4), ...
                                                 q(ed(i,1:2) + ed(i,3:4))); end
        end
        if isCvx
            cvx{end+1} = poly; %#ok<AGROW>
            atoms{end+1} = struct('kind', 'piece', 'q', singlePiece(f, k, eidx), 'region', poly); %#ok<AGROW>
            if n(5)*n(7) - n(6)^2 > 0
                Hi = inv(H);
                forms(end+1) = mkForm('piece', Hi, -Hi*L, 0.5*L.'*Hi*L - k0, [], [], [], [], ...
                                      H, L, poly); %#ok<AGROW>
            end
        end
        for i = 1:size(ed,1)
            u = ed(i,1:2);  dv = ed(i,3:4);  isRay = ed(i,5) == 1;
            al = 0.5 * dv*H*dv.';
            if al <= 0, continue, end            % a ray of linear growth only shapes dom f* (D)
            be = (u*H + L.') * dv.';  ga = q(u);
            fm = mkForm('edge', dv.'*dv/(2*al), u.' - be*dv.'/(2*al), ...
                        -ga + be^2/(4*al), u, dv, al, be, [], [], []);
            if isRay, fm.tmax = inf; end
            forms(end+1) = fm; %#ok<AGROW>
            if isCvx, continue, end              % the convex piece's atom already holds it
            if isRay
                key = [u, dv / norm(dv)];               % a ray is its base and its direction
                if any(all(abs(rayAt - key) < 1e-12, 2)), continue, end
                rayAt(end+1,:) = key; %#ok<AGROW>
                atoms{end+1} = struct('kind', 'ray', 'u', u, 'w', [], ...
                    'q', QuaPol([u; u + dv], [1 2 0], f.f(k,:), [0 0])); %#ok<AGROW>
                continue
            end
            w = u + dv;
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
        covered = any(cellfun(@(P) inPiece(p, P), cvx));
        for j = 1:size(edgeAt,1)
            covered = covered || norm(p - edgeAt(j,1:2)) < 1e-12 || norm(p - edgeAt(j,3:4)) < 1e-12;
        end
        for j = 1:size(rayAt,1)
            covered = covered || norm(p - rayAt(j,1:2)) < 1e-12;   % a ray atom holds its base
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
            keep(i) = ~any(cellfun(@(P) inPiece(atoms{i}.u, P) && inPiece(atoms{i}.w, P), cvx));
        end
    end
    atoms = atoms(keep);
end

function fm = mkForm(kind, A, b, c, u, d, al, be, H, L, poly)
    fm = struct('kind', kind, 'A', A, 'b', b, 'c', c, 'u', u, 'd', d, 'al', al, 'be', be, ...
                'H', H, 'L', L, 'poly', poly, 'tmax', 1);
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

function tf = inPiece(p, P)
% p in a piece's region: a polygon (bounded) or a struct of half-planes (unbounded).
    if isstruct(P)
        R = P.rows;
        tf = all(R(:,1:2)*p(:) + R(:,3) >= -1e-9 * max(1, max(abs(R(:)))) * max(1, norm(p)));
    else
        tf = inPoly(p, P);
    end
end

function R = pieceRows(ed)
% The half-planes [a b c] (a x + b y + c >= 0) of a CONVEX unbounded piece, from its edges
% (u, dv, isRay), oriented by an interior point: the mean of the base points pushed along the sum
% of the unit ray directions, which lies strictly inside a 2-D convex piece with a vertex.
    rays = ed(ed(:,5) == 1, 3:4);
    ip = mean(ed(:,1:2), 1) + sum(rays ./ vecnorm(rays, 2, 2), 1);
    R = zeros(size(ed,1), 3);
    for i = 1:size(ed,1)
        nrm = [-ed(i,4), ed(i,3)];
        R(i,:) = [nrm, -nrm*ed(i,1:2).'];
        if R(i,1:2)*ip.' + R(i,3) < 0, R(i,:) = -R(i,:); end
    end
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
function tris = tripleFaces(forms, D)
% For every triple of forms, every common tie s0 whose three contact points lie in their atoms.
% With an unbounded input, a tie ON the boundary of D = dom f* also carries the normal cone there:
% the face <s0,x> - g(s0) holds on T + N_D(s0) (erratum Case E), G its generators. Off the
% boundary G is empty. A tie outside D keeps its triangle alone -- still a chord of f, so still
% an upper bound.
    tris = struct('T', {}, 's', {}, 'c', {}, 'G', {}, 'forms', {});
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
                    tris(end+1) = struct('T', M, 's', s0.', 'c', gval(F3(1), s0), ...
                                         'G', normalCone(D, s0), 'forms', [a b c]); %#ok<AGROW>
                end
            end
        end
    end
end

function G = normalCone(D, s0)
% Outward normals of the rows of D active at s0, as columns; empty unless s0 is on D's boundary.
    G = zeros(2,0);
    if isempty(D), return, end
    r = D(:,1:2)*s0 + D(:,3);
    tol = 1e-9 * max(1, max(abs(D(:))));
    if any(r < -tol), return, end                 % s0 is outside D: the triangle alone
    act = abs(r) <= tol * max(1, norm(s0));
    G = -D(act,1:2).';
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
            m = fm.u + t*fm.d;  ok = t >= -1e-9 && t <= fm.tmax + 1e-9;   % tmax = inf on a ray
        case 'piece'
            m = (fm.H \ (s - fm.L)).';  ok = inPiece(m, fm.poly);
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
