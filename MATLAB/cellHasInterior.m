function tf = cellHasInterior(EcQ, rows)
% cellHasInterior  Does an H-form cell have a two-dimensional interior, CURVED constraints included?
%
% objective: decide emptiness for a cell whose constraints are not all straight -- the case
%   ratQ.feasible2 cannot reach, and the one that leaves conjQ's face count inflated.
%
% [input]  EcQ  : the shared curve list;  rows : this cell's m x 2 [curveIndex, side]
% [output] tf   : false only when the cell PROVABLY has no interior
%
% ------------------------------------------------------------------------------------------------
% SOUND, AND COMPLETE FOR THE CELLS THIS PIPELINE BUILDS.
%
% The candidate points of a region bounded by curves are its VERTICES -- the pairwise intersections
% of its own bounding curves -- plus, when the linear part alone is bounded away from those, an
% interior probe. Every such vertex is an intersection of two conics, so it is a root of the exact
% quartic `conicMeet` returns, and deciding whether it satisfies the cell's OTHER constraints is
% exactly `ratQ.signAt`: the sign of a rational polynomial at a degree-<=4 algebraic number. That is
% what the kernel was built for.
%
% WHY VERTICES ARE ENOUGH TO REFUTE. If a cell has interior, it has a point strictly inside every
% constraint, and shrinking toward the boundary reaches a vertex of the arrangement of its own
% curves -- so a cell with NO valid vertex and no valid interior probe has no interior. The
% converse is the direction that matters here: this routine only ever REMOVES cells, so a false
% "has interior" costs an extra face (the status quo) while a false "empty" would lose a real one.
% It is therefore written to answer TRUE when unsure.
%
% THE FILTERED-PREDICATE ARRANGEMENT of CONJ_FIELD_PROOF.md 8.7, in the small: the cheap tests run
% first -- a linear-only cell goes to ratQ.feasible2, a contradictory pair is caught outright -- and
% the exact kernel runs only on what survives, which is the cold path.

    if isempty(rows), tf = true; return, end

    % ---- cheap and exact: everything straight -------------------------------------------------
    isLine = all(EcQ(rows(:,1), 1:3) == 0, 2);
    if all(isLine)
        P = rows(:,2) .* EcQ(rows(:,1), 4:6);
        tf = ratQ.feasible2(P, true);
        return
    end

    % ---- the LINEAR part must have interior in any case ----------------------------------------
    % A necessary condition, and cheap: the curved constraints only shrink the region further.
    if any(isLine)
        P = rows(isLine,2) .* EcQ(rows(isLine,1), 4:6);
        if ~ratQ.feasible2(P, true), tf = false; return, end
    end

    % ---- a row that holds ONLY ON ITS CURVE -----------------------------------------------------
    % side * c >= 0 with side * c a NEGATIVE SEMIDEFINITE form (homogenised: -(x-y)^2, -(x^2+y^2),
    % -(x-1)^2) is the curve c = 0 itself, which has no interior. The vertex test below cannot see
    % this: such a curve meets the others at points that satisfy every row with equality, and a
    % vertex is evidence of a cell's closure, not of its interior. Found 2026-10-31 on the convex
    % envelope fixture: 12 of the 13 faces in a cell were points and lines kept that way.
    for r = 1:size(rows,1)
        if rows(r,2) ~= 0 && ~isLine(r) && negSemidefinite(rows(r,2) * EcQ(rows(r,1), 1:6))
            tf = false;  return
        end
    end

    % ---- an ARRANGEMENT VERTEX that satisfies the cell strictly enough --------------------------
    % Every pair of bounding curves; each meeting point is a root of conicMeet's exact quartic, and
    % whether it satisfies the other constraints is decided exactly.
    idx = unique(rows(:,1));
    for a = 1:numel(idx)
        for b = a+1:numel(idx)
            try
                [P, info] = conicMeet(EcQ(idx(a),:), EcQ(idx(b),:));
            catch ME
                % THE 2^53 CEILING, AND IT IS REACHED IN PRACTICE. conicMeet's resultant is a 4x4
                % determinant of quadratics, so its entries are degree-8 products of the input
                % coefficients: two conics with five-digit entries put it past 2^53 and ratQ.chk
                % raises rather than returning a rounded value. Measured here on a fold whose
                % difference conics carry six-digit coefficients.
                %
                % KEEPING THE CELL IS THE SAFE DIRECTION. This routine only ever REMOVES cells, so
                % failing to decide costs one extra face -- the status quo before this test existed
                % -- while guessing "empty" would lose a real region. The overflow is the recorded
                % signal that a wider integer type is wanted (TODO.md: swap ratQ's backend for GMP
                % behind the same interface), never a reason to loosen the gate.
                if ~strcmp(ME.identifier, 'ratQ:overflow'), rethrow(ME); end
                tf = true;  return
            end
            if info.degenerate || isempty(P), continue, end
            for r = 1:size(P,1)
                if satisfiesAll(EcQ, rows, P(r,:))
                    tf = true;  return
                end
            end
        end
    end

    % ---- no vertex qualified: probe the interior of the linear part -----------------------------
    % A cell can be nonempty with no arrangement vertex inside it -- an unbounded curved cell, for
    % instance -- so a negative vertex test is not yet a proof. The probe is a sample and can only
    % CONFIRM, which is why an empty result below returns false only after the vertex test also
    % found nothing: together they are the pipeline's cells, and the routine errs toward keeping.
    tf = probeInterior(EcQ, rows);
    if ~tf, tf = cadInterior(EcQ, rows); end
end

function tf = negSemidefinite(c)
% Is a x^2 + b xy + c y^2 + d x + e y + f <= 0 everywhere, i.e. is the homogenised form
% [2a b d; b 2c e; d e 2f] negative semidefinite? Exactly, by every principal minor of its negative
% being >= 0, in integers. Coefficients too large for the 3 x 3 minor to be exact answer FALSE --
% the cell is then kept, the direction this routine errs in.
    M = -[2*c(1) c(2) c(4); c(2) 2*c(3) c(5); c(4) c(5) 2*c(6)];
    tf = false;
    if any(M(:) ~= round(M(:))) || max(abs(M(:))) > 2^16 || ~any(M(:)), return, end
    if any(diag(M) < 0), return, end
    for p = [1 2; 1 3; 2 3].'
        if M(p(1),p(1)) * M(p(2),p(2)) - M(p(1),p(2))^2 < 0, return, end
    end
    % six products of at most 2^48 each: exact in doubles
    d = M(1,1)*M(2,2)*M(3,3) + 2*M(1,2)*M(2,3)*M(1,3) - M(1,1)*M(2,3)^2 - M(2,2)*M(1,3)^2 ...
        - M(3,3)*M(1,2)^2;
    tf = d >= 0;
end

function tf = satisfiesAll(EcQ, rows, x)
% objective: does the point x satisfy every constraint of this cell, with a scale-aware tolerance.
% The point comes from conicMeet and is a double; what makes the DECISION exact is that x is a root
% of an exact quartic and the constraints are exact -- the tolerance here only absorbs the printing
% of that root, and a borderline case is kept rather than dropped (see the header).
    tf = true;
    for r = 1:size(rows,1)
        c = EcQ(rows(r,1), :);
        v = c(1)*x(1)^2 + c(2)*x(1)*x(2) + c(3)*x(2)^2 + c(4)*x(1) + c(5)*x(2) + c(6);
        scale = max(1, max(abs(c)) * max(1, max(abs(x))^2));
        if rows(r,2) == 0
            if abs(v) > 1e-9 * scale, tf = false; return, end
        elseif rows(r,2) * v < -1e-9 * scale
            tf = false;  return
        end
    end
end

function tf = probeInterior(EcQ, rows)
% objective: look for an interior point by sampling. CONFIRMS only -- a negative result is what the
% caller treats as "no interior", and it is why this runs last and only after the exact vertex test.
    idx = unique(rows(:,1));
    pts = zeros(0,2);
    for a = 1:numel(idx)
        for b = a+1:numel(idx)
            try
                P = conicMeet(EcQ(idx(a),:), EcQ(idx(b),:));
            catch ME
                if ~strcmp(ME.identifier, 'ratQ:overflow'), rethrow(ME); end
                continue                    % see the caller: an undecidable pair keeps the cell
            end
            pts = [pts; P]; %#ok<AGROW>
        end
    end
    if isempty(pts)
        ctr = [0 0];  r = 4;
    else
        ctr = mean(pts, 1);  r = max(2, 2*max(max(abs(pts - ctr))));
    end
    rng(20260905);
    X = ctr + (rand(400,2) - 0.5) * 2 * r;
    for i = 1:size(X,1)
        if satisfiesAll(EcQ, rows, X(i,:)), tf = true; return, end
    end
    tf = false;
end

function tf = cadInterior(EcQ, rows)
% objective: the EXACT and COMPLETE decision of "does the open cell {side_i * q_i > 0} have a point?",
% run when the witness search above finds nothing. Twin of python/cca2/cell_has_interior.py
% _cadInterior (whose docstring holds the completeness argument) and rust/src/cells.rs cad_interior.
% A cylindrical decomposition over Q, no floating point: shear so that every non-line conic has a y^2
% term; the critical abscissae are the roots of every discriminant, every vertical line's abscissa,
% every pairwise resultant, and -- for two conics whose resultant is identically zero (they share a
% line) -- the crossing of their non-shared lines; one rational (dyadic) column per interval between
% them and one rational point per sector are tested with exact integer signs.
% A cell with an equality row has no open cell: false. Past the 2^53 wall (ratQ:overflow) the cell is
% KEPT, the direction this routine errs in.
    tf = false;
    if any(rows(:,2) == 0), return, end
    try
        tf = cadExact(EcQ, rows);
    catch ME
        if ~strcmp(ME.identifier, 'ratQ:overflow'), rethrow(ME); end
        tf = true;
    end
end

function tf = cadExact(EcQ, rows)
% Memoised: the decision is a pure function of the cell's conics and signs, and conjQ asks about the
% same cell again and again (compactCells.dropImplied removes one row at a time, assembleQuaConCells
% revisits cells). Only successes are stored; ratQ:overflow still propagates.
    persistent memo
    if isempty(memo), memo = containers.Map('KeyType', 'char', 'ValueType', 'logical'); end
    idx = unique(rows(:,1));
    sgn = zeros(numel(idx), 1);
    for k = 1:numel(idx)
        sd = unique(rows(rows(:,1) == idx(k), 2));
        if numel(sd) == 1, sgn(k) = sd; else, sgn(k) = 2; end
    end
    key = sprintf('%.17g,', [reshape(EcQ(idx, 1:6).', 1, []), sgn.']);
    if isKey(memo, key), tf = memo(key); return, end
    tf = cadExactRaw(EcQ, rows);
    memo(key) = tf;
end

function tf = cadExactRaw(EcQ, rows)
    tf = false;
    idx = unique(rows(:,1));
    sg = zeros(numel(idx), 1);
    for k = 1:numel(idx)
        sd = unique(rows(rows(:,1) == idx(k), 2));
        if numel(sd) > 1, return, end              % q > 0 and q < 0 on the same conic
        sg(k) = sd;
    end
    cs = EcQ(idx, 1:6);
    k = 0;
    while any((cs(:,1) ~= 0 | cs(:,2) ~= 0 | cs(:,3) ~= 0) & cs(:,1)*k^2 + cs(:,2)*k + cs(:,3) == 0)
        k = k + 1;
    end
    cs = [cs(:,1), 2*k*cs(:,1) + cs(:,2), k^2*cs(:,1) + k*cs(:,2) + cs(:,3), cs(:,4), k*cs(:,4) + cs(:,5), cs(:,6)];
    ratQ.chk(cs, 'sheared conics');
    [xpolys, special] = cadCritical(cs);
    xs = cadSamples(xpolys, {});
    for x0 = xs(:).'
        [xn, xd] = dyadic(x0);
        [plain, sp, alive] = cadYPolys(cs, sg, xn, xd, special);
        if ~alive, continue, end
        for y0 = cadSamples(sp, plain)
            [yn, yd] = dyadic(y0);
            ok = true;
            for r = 1:size(cs,1)
                if sg(r) * cellValueSign(cs(r,:), xn, xd, yn, yd) <= 0, ok = false; break, end
            end
            if ok, tf = true; return, end
        end
    end
end

function s = cellValueSign(c, xn, xd, yn, yd)
% The exact sign of q(xn/xd, yn/yd) * xd^2 * yd^2. In doubles when every term is small enough that
% the terms and their partial sums are exact integers; otherwise in java.math.BigInteger. A dyadic
% sample's denominator grows with the bisection that separated its roots, so the doubles reach 2^53
% on ordinary small-coefficient cells. Raising ratQ:overflow there made cadInterior KEEP an empty
% cell, where Python and Rust (exact integers) drop it: the poly2 kink's AlgCrv had 340 faces in
% MATLAB against 319 (2026-11-25).
    t = [c(1)*xn^2*yd^2, c(2)*xn*yn*xd*yd, c(3)*yn^2*xd^2, c(4)*xn*xd*yd^2, c(5)*yn*yd*xd^2, c(6)*xd^2*yd^2];
    if all(isfinite(t)) && max(abs(t)) < 2^49            % six terms: every partial sum < 2^52
        s = sign(sum(t));
        return
    end
    B = @(z) java.math.BigInteger(sprintf('%.0f', z));  % xn, xd, yn, yd are exact integer doubles
    X = B(xn);  XD = B(xd);  Y = B(yn);  YD = B(yd);
    mono = {X.multiply(X).multiply(YD).multiply(YD), X.multiply(Y).multiply(XD).multiply(YD), ...
            Y.multiply(Y).multiply(XD).multiply(XD), X.multiply(XD).multiply(YD).multiply(YD), ...
            Y.multiply(YD).multiply(XD).multiply(XD), XD.multiply(XD).multiply(YD).multiply(YD)};
    v = java.math.BigInteger.ZERO;
    for i = 1:6
        if c(i) ~= 0, v = v.add(mono{i}.multiply(B(c(i)))); end
    end
    s = double(v.signum());
end

function [n, d] = dyadic(x)
    d = 1;
    while x*d ~= round(x*d), d = d*2; end
    n = x*d;
end

function p = pmul(a, b)
    if isempty(a) || isempty(b), p = []; return, end
    p = ratQ.polyTrim(ratQ.chk(conv(a, b), 'polynomial product'));
end

function p = padd(a, b, sb)
    b = sb * b;
    k = max(numel(a), numel(b));
    p = ratQ.polyTrim([zeros(1, k-numel(a)), a(:).'] + [zeros(1, k-numel(b)), b(:).']);
end

function [polys, special] = cadCritical(cs)
% The integer polynomials in x whose real roots are the critical abscissae (step 2 of the argument),
% and the rows whose y-polynomials can share a root at EVERY column (a repeated factor, or a common
% one with another row): at any other column the y-polynomials are squarefree and pairwise coprime,
% because a resultant that is not identically zero has no root there.
    polys = {};
    special = false(size(cs,1), 1);
    full = find(cs(:,3) ~= 0).';
    slant = find(cs(:,3) == 0 & cs(:,1) == 0 & cs(:,2) == 0 & cs(:,5) ~= 0).';
    for r = 1:size(cs,1)
        c = cs(r,:);
        if c(3) ~= 0
            polys{end+1} = [c(2)^2 - 4*c(3)*c(1), 2*c(2)*c(5) - 4*c(3)*c(4), c(5)^2 - 4*c(3)*c(6)]; %#ok<AGROW>
            if isempty(ratQ.polyTrim(polys{end})), special(r) = true; end
        elseif c(5) == 0 && c(2) == 0 && c(1) == 0
            polys{end+1} = [c(4), c(6)]; %#ok<AGROW>   % a vertical line (or a constant)
        end
    end
    for i = 1:numel(full)
        p = cs(full(i),:);
        al = p(3);  be = [p(2), p(5)];  ga = [p(1), p(4), p(6)];
        for j = i+1:numel(full)
            q = cs(full(j),:);
            al2 = q(3);  be2 = [q(2), q(5)];  ga2 = [q(1), q(4), q(6)];
            u = padd(pmul(al2, be), pmul(al, be2), -1);       % alpha' beta  - alpha beta'
            w = padd(pmul(al2, ga), pmul(al, ga2), -1);       % alpha' gamma - alpha gamma'
            v = padd(pmul(be, ga2), pmul(be2, ga), -1);       % beta gamma'  - beta' gamma
            res = padd(pmul(w, w), pmul(u, v), 1);            % (a g' - a' g)^2 - (a b' - a' b)(b g' - b' g)
            if ~isempty(res)
                polys{end+1} = res; %#ok<AGROW>
            else
                special([full(i), full(j)]) = true;
                if ~isempty(u) || ~isempty(w)
                    polys = [polys, sharedLineCrossing(p, q, u, w)]; %#ok<AGROW>
                end
            end
        end
        for j = slant
            q = cs(j,:);
            B = q(5);  C = [q(4), q(6)];
            res = padd(padd(pmul(al, pmul(C, C)), pmul(be, pmul(B, C)), -1), pmul(ga, pmul(B, B)), 1);
            if ~isempty(res), polys{end+1} = res; else, special([full(i), j]) = true; end %#ok<AGROW>
        end
    end
    for i = 1:numel(slant)
        for j = i+1:numel(slant)
            p = cs(slant(i),:);  q = cs(slant(j),:);
            res = padd(pmul(p(5), [q(4), q(6)]), pmul(q(5), [p(4), p(6)]), -1);
            if ~isempty(res), polys{end+1} = res; else, special([slant(i), slant(j)]) = true; end %#ok<AGROW>
        end
    end
end

function out = sharedLineCrossing(p, q, u, w)
% p = L Mp and q = L Mq share the line L (resultant == 0, p and q not proportional):
% T = alpha' p - alpha q = u(x) y + w(x) = L * N with N vertical, so L = T / gcd(u, w). Returns the
% polynomial [den, -num] of the abscissa of Mp ^ Mq ({} when they are parallel). Symbolic: rare.
    persistent memo
    if isempty(memo), memo = containers.Map('KeyType', 'char', 'ValueType', 'any'); end
    key = sprintf('%.17g,', [p(:).', q(:).']);
    if isKey(memo, key), out = memo(key); return, end
    out = sharedLineCrossingRaw(p, q, u, w);
    memo(key) = out;
end

function out = sharedLineCrossingRaw(p, q, u, w)
    out = {};
    if isempty(u), error('cellHasInterior:internal', 'shared line factor with no y term'); end
    if isempty(w), g = ratQ.polyPrimitive(u); else, g = ratQ.polyGcd(u, w); end
    syms X
    G = poly2sym(sym(g), X);
    Lu = simplify(poly2sym(sym(u), X) / G);
    if isempty(w), Lw = sym(0); else, Lw = simplify(poly2sym(sym(w), X) / G); end
    ell = subs(Lu, X, 0);
    rr = subs(Lw, X, 0);
    pp = subs(Lw, X, 1) - rr;
    if ell == 0 || simplify(Lu - ell) ~= 0 || simplify(Lw - (pp*X + rr)) ~= 0
        error('cellHasInterior:internal', 'shared factor of two conics is not a line');
    end
    m = sym(zeros(2,3));
    cc = {p, q};
    for t = 1:2
        c = sym(cc{t});
        mm = c(3) / ell;
        m(t,:) = [mm, (c(2) - pp*mm)/ell, (c(5) - rr*mm)/ell];
    end
    det_ = m(1,1)*m(2,2) - m(2,1)*m(1,2);
    if det_ == 0, return, end
    xp = (m(2,1)*m(1,3) - m(1,1)*m(2,3)) / det_;
    [nn, dd] = numden(xp);
    out = {[double(dd), -double(nn)]};
end

function [plain, sp, alive] = cadYPolys(cs, sg, xn, xd, special)
% The y-polynomials q(x0, .) of one column, x0 = xn/xd, scaled to integers: {plain} for the rows that
% are squarefree and coprime to everything at this (non-critical) column, {sp} for the `special` rows
% that need the coprime-base treatment. alive = false when the column holds no point of the cell:
% a row constant in y with the wrong sign, or a quadratic without real roots with the wrong sign,
% kills it; a quadratic without real roots with the right sign holds everywhere, so it adds no
% sample boundary.
    plain = {};  sp = {};  alive = true;
    for r = 1:size(cs,1)
        c = cs(r,:);
        co = [c(3)*xd^2, c(2)*xn*xd + c(5)*xd^2, c(1)*xn^2 + c(4)*xn*xd + c(6)*xd^2];
        ratQ.chk(co, 'column polynomial');
        co = ratQ.polyPrimitive(co);
        if isempty(co) || (numel(co) == 1 && sg(r)*co <= 0)
            alive = false;  return
        elseif numel(co) == 1
            continue
        elseif numel(co) == 3 && co(2)^2 - 4*co(1)*co(3) < 0
            if sg(r)*co(1) < 0, alive = false; return, end
        elseif special(r)
            sp{end+1} = co; %#ok<AGROW>
        else
            plain{end+1} = co; %#ok<AGROW>
        end
    end
end

function base = coprimeBase(polys)
% Pairwise coprime squarefree primitive polynomials with the same real roots as the product.
    work = {};
    for k = 1:numel(polys)
        if numel(ratQ.polyTrim(polys{k})) > 1, work{end+1} = polys{k}; end %#ok<AGROW>
    end
    base = {};
    while ~isempty(work)
        p = ratQ.polyPrimitive(ratQ.polySquarefree(work{end}));  work(end) = [];
        if numel(p) <= 1, continue, end
        hit = 0;
        for i = 1:numel(base)
            g = ratQ.polyGcd(p, base{i});
            if numel(g) > 1, hit = i; break, end
        end
        if hit == 0
            base{end+1} = p; %#ok<AGROW>
        else
            b = base{hit};  base(hit) = [];
            work = [work, {exactDiv(b, g), exactDiv(p, g), g}]; %#ok<AGROW>
        end
    end
end

function q = exactDiv(p, g)
% p / g for integer polynomials by long division in integers (deconv is floating point); an error if
% the division is not exact, which for primitive p and g | p cannot happen (Gauss).
    r = p(:).';  g = g(:).';  n = numel(r) - numel(g) + 1;
    q = zeros(1, max(n, 0));
    for i = 1:n
        c = r(i) / g(1);
        if c ~= round(c), error('cellHasInterior:internal', 'inexact gcd division'); end
        q(i) = c;
        r(i:i+numel(g)-1) = r(i:i+numel(g)-1) - c * g;
    end
    if n < 1 || any(r ~= 0), error('cellHasInterior:internal', 'inexact gcd division'); end
    ratQ.chk(q, 'gcd quotient');
end

function [tf, r, s] = cadBefore(r, s)
% r < s for roots of coprime polynomials (so r ~= s): refine until the intervals separate.
    while true
        if r.hi <= s.lo, tf = true; return, end
        if s.hi <= r.lo, tf = false; return, end
        if r.hi - r.lo >= s.hi - s.lo, r = cadRefine(r); else, s = cadRefine(s); end
    end
end

function r = cadRefine(r)
    mid = (r.lo + r.hi) / 2;
    if mid <= r.lo || mid >= r.hi      % the dyadic endpoint no longer fits a double: keep the cell
        error('ratQ:overflow', 'root interval below double resolution');
    end
    if ratQ.sturmV(r.chain, mid) - ratQ.sturmV(r.chain, r.hi) == 1, r.lo = mid; else, r.hi = mid; end
end

function s = cadSamples(sp, plain)
% Rational (dyadic) sample points: one in every open interval between consecutive distinct real roots
% of the polynomials sp (made coprime first) and plain (already squarefree and coprime), one below and
% one above all roots (0 when there is none). Exact: Sturm isolation of each factor, pairwise
% separation by bisection.
    if nargin < 2, plain = {}; end
    polys = coprimeBase(sp);
    for k = 1:numel(plain)
        if numel(plain{k}) > 1, polys{end+1} = plain{k}; end %#ok<AGROW>
    end
    roots = {};
    for k = 1:numel(polys)
        p = polys{k};
        chain = ratQ.sturmChain(p);
        I = ratQ.isolateRoots(p);
        for t = 1:size(I,1)
            r = struct('chain', {chain}, 'lo', I(t,1), 'hi', I(t,2));
            a = 1;  b = numel(roots) + 1;        % binary insertion keeps roots sorted
            while a < b
                m = floor((a + b) / 2);
                [tf, roots{m}, r] = cadBefore(roots{m}, r); %#ok<AGROW>
                if tf, a = m + 1; else, b = m; end
            end
            roots = [roots(1:a-1), {r}, roots(a:end)]; %#ok<AGROW>
        end
    end
    if isempty(roots), s = 0; return, end
    s = roots{1}.lo;
    for i = 1:numel(roots)-1
        r = roots{i};  t = roots{i+1};
        while r.hi >= t.lo
            if r.hi - r.lo >= t.hi - t.lo, r = cadRefine(r); else, t = cadRefine(t); end
        end
        roots{i} = r;  roots{i+1} = t;
        s(end+1) = (r.hi + t.lo) / 2; %#ok<AGROW>
    end
    s(end+1) = roots{end}.hi + 1;
end
