function tf = algRowFails(a, box, inBox, depth, factors, corners)
% algRowFails  Is an algebraic row side * (q(x) - z(x)) >= 0 FALSE at every point of a region,
%   CERTIFIED by Bernstein bounds on sub-boxes: true only when that was shown, false otherwise.
%   Given SEVERAL rows, is at least one of them false at every point -- their conjunction empty.
%
% [input]  a     : the row(s), as envelopeSubdivisionAlg stores them -- .FP (polyN in x, y, z of degree 1
%                  or 2 in z), .FB (algBranch's index), .num/.den (the quadratic face, weighted
%                  basis), .side (+1 or -1); a struct array for several rows. A point where the
%                  row's branch does not exist counts as failing, as every reader takes it
%          box   : [x0 x1 y0 y1], a box containing the region
%          inBox : function handle, inBox(b) false only when sub-box b PROVABLY misses the region
%          depth : subdivision levels allowed (default 6)
%          factors : optional cell of dense polyN, each >= 0 on the region (the region's own rows,
%                  integer coefficients); see A ROW VANISHING ON THE REGION'S EDGE
%          corners : optional k x 2, candidate points where every row may vanish (the region's
%                  vertices); see A ROW ZERO AT A CORNER
% [output] tf    : true when at every point of the region in box, side * (q - z) < 0 for some row --
%                  except on the zero set of a squared factor (A TOUCHING ROW below), a curve
%
% SEVERAL ROWS (2026-10-30): 274 faces of the 3-convex-edge triangle were empty only JOINTLY -- each
% row holds somewhere on the piece, never all together. A sub-box is decided when ANY row is
% certified false on it, so different rows may decide different boxes.
%
% WHY POINTWISE (2026-10-29). The first test (envelopeSubdivisionAlg, 2026-10-20) proved the
% row's sign CONSTANT on a cell -- R = FP(x, y, q) != 0 there -- and read it at sample points, so a
% cell the curve R = 0 crosses could never be decided: ~1600 faces of the 3-convex-edge triangle.
% But R = 0 is where q meets EITHER root, and a box the OTHER root's curve crosses can still have
% one sign of q - z. With a, b, c the coefficients of z^2, z, 1 and m = -b / (2a) the roots' mean,
%
%       q < z_hi  <=>  q < m  or  F(q)/a < 0        q > z_hi  <=>  q > m  and  F(q)/a > 0
%       q < z_lo  <=>  q < m  and F(q)/a > 0        q > z_lo  <=>  q > m  or   F(q)/a < 0
%
% (F(q)/a < 0 is "q between the roots"), with D = b^2 - 4ac > 0 for both roots to exist. Each of
% sign(a), sign(q - m) = sign(a) sign(2aq + b), sign(F(q)) and sign(D) is a polynomial sign on the
% box, which polySignOn certifies; a box that is not decided is split, and the region is covered
% when every sub-box meeting it is decided. No connectivity, no sample point: the answer holds at
% every point of the region, which is what dropping a face needs.
%
% Degree 1 in z: z = -c/b and sign(q - z) = sign(b) sign(b q + c).
%
% A TOUCHING ROW (2026-10-30). Where q meets z without crossing it -- a tangential tie, on the
% 3-convex-edge triangle along x = y, x = 0 and the like -- F(q) carries a SQUARED factor, G^2, and
% no box containing a point of G = 0 ever gets a strict sign. Divided out, F(q) / G^2 has F(q)'s sign
% everywhere off G = 0, so a certificate then says the row fails everywhere EXCEPT ON A CURVE: the
% piece has no interior, which is what dropping a face needs (its points on that curve lie on the
% boundary of faces with area, whose rows are closed). The same for b q + c on a linear fibre.
%
% A ROW VANISHING ON THE REGION'S EDGE (2026-10-30). 299 faces of the 3-convex-edge triangle were left
% undecided only at boxes straddling the hypotenuse, where a row's q - z is a multiple of
% 1 - x - y: it fails inside, is zero on the edge and holds outside, so no box the edge crosses has
% one sign. A polynomial g >= 0 on the region that divides F(q) exactly, F(q) = g M, leaves
% sign F(q) = sign M wherever g > 0, i.e. on the region off the curve g = 0 -- the same claim as a
% touching row. So the region's rows are divided out like the squared factors. And not only from
% F(q): at the triangle's corner (0, 0), b of a linear fibre and D of a quadratic one are multiples
% of the row x + y >= 0 -- zero at the corner, no box containing it had their sign (228 faces). Every
% polynomial of the test is reduced the same way; each keeps its sign off a curve.
%
% A ROW ZERO AT A CORNER (2026-10-30). Left after that: 20 faces at that corner where a row's
% polynomial is 0 at (0, 0) and positive on the region otherwise, its zero curve tangent to the region
% from outside -- every box containing the corner also holds points of the wrong sign, at any depth.
% BLOWN UP at the corner p, each polynomial is P(p + s (1, t)) = s^k G(s, t) on the chart
% |y - py| <= |x - px|, x > px (and the three other charts: x < px, and y-major), with G's sign that
% of P for s > 0 -- and G(0, t) is P's lowest-order form along direction (1, t), which the tangency
% does not make zero inside the region. A tangency leaves undecided boxes at every scale (~ their
% distance to p, measured out to 6 boxes), so it is the SQUARE of half-width e around p that is
% certified: every row test passes on the four charts [0, e] x [-1, 1] covering it (p itself
% excepted, a point), the region's factors blown up the same way rejecting chart boxes; a box inside
% the square is then decided. The substitution is exact or refused (dyadic p, bounded integers), and only exactly
% zero lowest orders are divided out, so any p is sound -- a wrong p only fails to certify.
    if nargin < 4 || isempty(depth), depth = 6; end
    if nargin < 5, factors = {}; end
    if nargin < 6, corners = zeros(0,2); end
    tf = false;
    rows = cell(1, numel(a));
    for t = 1:numel(a)
        rows{t} = rowPolys(a(t), factors);
    end
    rows = rows(~cellfun(@isempty, rows));
    if isempty(rows), return, end
    tf = boxRec(rows, box, inBox, depth, zeros(0,3));
    if tf || isempty(corners), return, end
    % CORNERS, only once the plain subdivision has failed: around each corner in the box, the largest
    % square (a quarter of the box, halved down to 1/256) whose blown-up charts all pass; a box
    % inside one is decided.
    W = max(box(2) - box(1), box(4) - box(3));
    sq = zeros(0,3);
    for c = 1:size(corners, 1)
        p = corners(c,:);
        if p(1) < box(1) || p(1) > box(2) || p(2) < box(3) || p(2) > box(4), continue, end
        if ~anyVanishes(rows, p), continue, end      % nothing to blow up: G(0, t) = P(p), as plain
        for e = W ./ 2.^(2:8)
            if cornerOK(rows, p, e, factors), sq(end+1,:) = [p, e]; break, end %#ok<AGROW>
        end
    end
    if ~isempty(sq), tf = boxRec(rows, box, inBox, depth, sq); end
end

function r = rowPolys(a, factors)
% The polynomials whose signs decide row a on a box, and the test reading them; [] if the row
% cannot be certified at all (shape, branch index, non-finite coefficients).
    r = [];
    FP = a.FP;
    num = double(a.num(:).');
    if numel(num) ~= 10 || any(num(1:4) ~= 0) || ~(a.den > 0), return, end
    % q = Q / Dn: twice the weighted numerator, in the plain monomial basis (as bridgeQuadCut)
    Q = zeros(3,3);
    Q(3,1) = num(5);    Q(2,2) = 2*num(6);  Q(1,3) = num(7);
    Q(2,1) = 2*num(8);  Q(1,2) = 2*num(9);  Q(1,1) = 2*num(10);
    Dn = 2 * double(a.den);
    want = -a.side;                                  % the sign q - z must have for the row to fail
    if ratQ.isBigPoly(FP), nz = FP.sz(3); else, nz = size(FP, 3); end
    switch nz - 1
        case 2
            hi = any(a.FB == [-1 2]);
            if ~hi && ~any(a.FB == [1 -2]), return, end
            test = @(s) quadOK(s, hi, want);
        case 1
            if ~any(a.FB == [1 -1]), return, end
            test = @(s) s(1) ~= 0 && s(1) * s(2) == want;
        otherwise
            return
    end
    % BIG INTEGERS (2026-11-02). Erratum Example 2's bridges have FP past 2^53, and ratQ.fromBig
    % raised here -- so no row of theirs was ever certified. And a double FP can overflow in the
    % products below, silently. Both go through BigInteger: the same formulas (signPolys) on bigZ,
    % the content divided out, and the result back in doubles -- exactly when it fits, otherwise
    % ROUNDED to the nearest double. That rounding is sound for what reads it: polySignOn counts a
    % Bernstein coefficient only past a margin of 1e-10 of the largest, and a relative error of
    % 2^-53 per coefficient moves each Bernstein coefficient by far less (they are convex
    % combinations). The exact reductions below (squared factors, region factors) need integers
    % within 2^53 and simply refuse a rounded polynomial, which keeps it whole -- sign-correct.
    % REMEMBERED by the row's integers (2026-11-01): the products and the square-free reduction do
    % not depend on the region, and on an erratum Example 2 overlay the same row came back for
    % piece after piece. The key is exact -- every coefficient, the big ones as decimal strings.
    persistent memo
    if isempty(memo), memo = containers.Map('KeyType', 'char', 'ValueType', 'any'); end
    key = rowKey(FP, num, a.den);
    if ~isempty(key) && isKey(memo, key)
        P = memo(key);
    else
        P = squareFreePolys(FP, Q, Dn, nz);
        if isempty(P), return, end
        if ~isempty(key)
            if memo.Count > 5000, memo.remove(memo.keys()); end
            memo(key) = P;
        end
    end
    % EVERY sign is read with its squared factors and the region's factors divided out: each then
    % has the original sign off a finite union of curves, so the conclusion holds off them too
    P = cellfun(@(p) divideOut(p, factors), P, 'UniformOutput', false);
    r = struct('P', {P}, 'test', test);
end

function P = squareFreePolys(FP, Q, Dn, nz)
% The row's sign polynomials with their squared factors divided out; {} when they cannot be formed.
    big = ratQ.isBigPoly(FP);
    if ~big
        bound = signPolys(abs(FP), abs(Q), Dn, nz - 1, 4);  % every partial sum of the products, bounded
        big = any(cellfun(@(p) any(p(:) >= 2^53), bound));
    end
    try
        if big
            P = signPolys(toZ(FP), Q, bigZ(Dn), nz - 1);
            P = cellfun(@fromZ, P, 'UniformOutput', false);
        else
            P = signPolys(FP, Q, Dn, nz - 1);
        end
    catch
        P = {};  return
    end
    if any(cellfun(@(p) ~all(isfinite(p(:))), P)), P = {}; return, end
    P = cellfun(@squareFreeSign, P, 'UniformOutput', false);
end

function key = rowKey(FP, num, den)
% An exact text key for (FP, num, den); '' when one cannot be formed cheaply.
    try
        if ratQ.isBigPoly(FP)
            c = FP.c;  e = cellfun(@isempty, c);
            s = cell(size(c));  s(e) = {'0'};
            s(~e) = cellfun(@(v) char(v.toString()), c(~e), 'UniformOutput', false);
            k1 = [sprintf('%d,', FP.sz), '|', strjoin(s(:).', ',')];
        else
            k1 = [sprintf('%d,', size(FP)), '|', sprintf('%.17g,', FP(:))];
        end
        key = [k1, '|', sprintf('%.17g,', double(num)), '|', sprintf('%.17g', double(den))];
    catch
        key = '';
    end
end

function P = signPolys(FP, Q, Dn, deg, k4)
% The sign polynomials of a row, from FP's z-coefficients: in doubles, or in bigZ (same code).
% k4 = -4 (D = b^2 - 4ac); a bound on |every partial sum| is the same with |inputs| and k4 = +4.
    if nargin < 5, k4 = -4; end
    if deg == 2
        A = FP(:,:,3);  B = FP(:,:,2);  C0 = FP(:,:,1);
        P = {A, ...                                                       % a
             padAdd(2 * conv2g(A, Q), Dn * B), ...                        % Dn (2aq + b)
             padAdd(padAdd(conv2g(A, conv2g(Q, Q)), Dn * conv2g(B, Q)), Dn^2 * C0), ...  % Dn^2 F(q)
             padAdd(conv2g(B, B), k4 * conv2g(A, C0))};                   % D
    else
        A = FP(:,:,2);  C0 = FP(:,:,1);
        P = {A, padAdd(conv2g(A, Q), Dn * C0)};                          % b, Dn (b q + c)
    end
end

function Z = toZ(FP)
% ratQ's big form (or a double array) as a bigZ array of the same shape.
    if ~ratQ.isBigPoly(FP), Z = bigZ(FP); return, end
    c = FP.c;
    c(cellfun(@isempty, c)) = {java.math.BigInteger.ZERO};
    Z = bigZ(reshape(c, FP.sz));
end

function p = fromZ(p)
% Content divided out; exact doubles when they fit, the nearest ones otherwise (see BIG INTEGERS).
    if ~isa(p, 'bigZ'), return, end
    p = contentFree(p);
    try
        p = exactDouble(p, 'row sign polynomial');
    catch
        p = double(p);
    end
end

function C = conv2g(A, B)
% conv2, also for bigZ operands.
    if isnumeric(A) && isnumeric(B), C = conv2(A, B); return, end
    [a1, a2] = size(A);  [b1, b2] = size(B);
    C = bigZ(zeros(a1 + b1 - 1, a2 + b2 - 1));
    for i = 1:a1
        for j = 1:a2
            C(i:i+b1-1, j:j+b2-1) = C(i:i+b1-1, j:j+b2-1) + A(i,j) .* B;
        end
    end
end

function tf = boxRec(rows, b, inBox, depth, sq)
% sq: k x 3, squares [px py e] proved (cornerOK) to fail some row everywhere but at p.
% The row tests come BEFORE inBox (2026-11-02): a box they decide is decided whatever inBox says, and
% inBox -- every region factor's sign on the box -- was 75 of 100 s of algRowFails on an Example 2
% overlay tuple, most of it on chart boxes the rows decide at once.
    for c = 1:size(sq, 1)
        if b(1) >= sq(c,1) - sq(c,3) && b(2) <= sq(c,1) + sq(c,3) && ...
           b(3) >= sq(c,2) - sq(c,3) && b(4) <= sq(c,2) + sq(c,3)
            tf = true;  return
        end
    end
    for t = 1:numel(rows)
        P = rows{t}.P;  s = zeros(1, numel(P));
        for k = 1:numel(P)
            v = polySignOn(P{k}, b, @(~) true, 0);
            if ~isnan(v), s(k) = v; end
        end
        if rows{t}.test(s), tf = true; return, end
    end
    if ~inBox(b), tf = true; return, end             % misses the region: nothing to prove here
    if depth == 0, tf = false; return, end
    xm = (b(1) + b(2)) / 2;  ym = (b(3) + b(4)) / 2;
    tf = boxRec(rows, [b(1) xm b(3) ym], inBox, depth-1, sq) && ...
         boxRec(rows, [xm b(2) b(3) ym], inBox, depth-1, sq) && ...
         boxRec(rows, [b(1) xm ym b(4)], inBox, depth-1, sq) && ...
         boxRec(rows, [xm b(2) ym b(4)], inBox, depth-1, sq);
end

function tf = anyVanishes(rows, p)
% Does some polynomial of some row vanish at p (to rounding)? Only then can a blow-up at p decide
% what the plain subdivision did not: otherwise every chart's lowest order is P(p) itself, a
% nonzero constant, which a small enough plain box reads as well (2026-11-02: the region's
% vertices are many on an overlay tuple, and the charts were tried at every one).
    tf = false;
    for t = 1:numel(rows)
        for k = 1:numel(rows{t}.P)
            if vanishesAt(rows{t}.P{k}, p), tf = true; return, end
        end
    end
end

function tf = vanishesAt(P, p)
% P(p) = 0 to rounding (a dense polyN in x, y).
    T = P .* (p(1) .^ (0:size(P,1)-1).') .* (p(2) .^ (0:size(P,2)-1));
    tf = abs(sum(T(:))) <= 1e-9 * max(sum(abs(T(:))), realmin);
end

function ok = cornerOK(rows, p, e, factors)
% Every point of the square of half-width e around p (p excepted) fails some row: the four blown-up
% charts, each certified by boxRec on [0, e] x [-1, 1]. See A ROW ZERO AT A CORNER.
    ok = false;
    if ~(e > 0), return, end
    for major = 1:2
        for sg = [1 -1]
            rowsC = rows;
            for t = 1:numel(rows)
                for k = 1:numel(rows{t}.P)
                    G = blowUp(rows{t}.P{k}, p, sg, major);
                    if isempty(G), return, end
                    rowsC{t}.P{k} = G;
                end
            end
            facC = cell(1, numel(factors));
            for k = 1:numel(factors)
                % only the region's rows THROUGH p: one positive at p rejects nothing near it, and
                % leaving a factor out only keeps more chart boxes -- safe. On an Example 2 overlay
                % tuple the charts tested ~60 factors per box and that was 49 of their 53 s.
                if ~vanishesAt(factors{k}, p), continue, end
                facC{k} = blowUp(factors{k}, p, sg, major);
            end
            facC = facC(~cellfun(@isempty, facC));
            inC = @(cb) ~any(cellfun(@(G) polySignOn(G, cb, @(~) true, 0) == -1, facC));
            if ~boxRec(rowsC, [0 e -1 1], inC, 6, zeros(0,3)), return, end
        end
    end
    ok = true;
end

function G = blowUp(P, p, sg, major)
% P(p + sg s (1, t)) (major 1) or P(p + sg s (t, 1)) (major 2) divided by the largest power of s that
% divides it EXACTLY, as a polyN in (s, t); [] when the shift could not be done exactly in doubles.
    G = [];
    for it = 1:30                                    % a dyadic P scaled to integers (sign kept)
        if all(P(:) == round(P(:))), break, end
        P = 2 * P;
    end
    sc = 1;                                          % p's denominator, a power of 2 up to 2^8
    while any(p * sc ~= round(p * sc)) && sc < 2^8, sc = 2 * sc; end
    if any(p * sc ~= round(p * sc)) || any(P(:) ~= round(P(:))), return, end
    [n1, n2] = size(P);
    deg = n1 + n2 - 2;
    bound = sum(abs(P(:))) * (abs(p(1)) + 1)^(n1-1) * (abs(p(2)) + 1)^(n2-1) * 2^deg * sc^deg;
    if bound >= 2^52, return, end                    % every partial sum an integer multiple of sc^-deg
    Ps = shiftTo(shiftTo(P, p(1)).', p(2)).';        % Ps(i+1, j+1): coefficient of u^i v^j
    G = zeros(n1 + n2 - 1, max(n1, n2));
    for i = 0:n1-1
        for j = 0:n2-1
            c = Ps(i+1, j+1) * sg^(i+j);
            if major == 1, G(i+j+1, j+1) = G(i+j+1, j+1) + c;
            else,          G(i+j+1, i+1) = G(i+j+1, i+1) + c;
            end
        end
    end
    k = find(any(G ~= 0, 2), 1);
    if isempty(k), G = 0; return, end
    G = G(k:end, :);
end

function S = shiftTo(P, x0)
% Rows of P re-expanded about x = x0: S(k+1, :) = sum_i C(i, k) x0^(i-k) P(i+1, :).
    n1 = size(P, 1);
    S = zeros(size(P));
    for i = 0:n1-1
        for k = 0:i
            S(k+1, :) = S(k+1, :) + nchoosek(i, k) * x0^(i-k) * P(i+1, :);
        end
    end
end

function ok = quadOK(s, hi, want)
% s = signs of [a, Dn(2aq + b), Dn^2 F(q), D] on the box, 0 where not certified.
    if s(1) == 0, ok = false; return, end
    u = s(2) * s(1);                                 % sign(q - m)
    v = s(3) * s(1);                                 % sign(F(q) / a): -1 is "between the roots"
    % WHERE THE BRANCH DOES NOT EXIST (D < 0) THE ROW FAILS -- every reader takes it so (the
    % evaluator's algHolds and algEmpty require a finite value) -- so no conclusion needs a sign of
    % D: each is right wherever the roots exist, and D < 0 fails the row by itself (2026-11-02).
    % Until then u = -1 and u = +1 alone required D > 0, and a box on a face's discriminant curve --
    % erratum Example 2's faces end on theirs at the cells' edges -- never had that sign.
    if hi && want == -1,      ok = u == -1 || v == -1;
    elseif hi,                ok = u == 1 && v == 1;
    elseif want == -1,        ok = u == -1 && v == 1;
    else,                     ok = u == 1 || v == -1;
    end
end

function S = squareFreeSign(R)
% R with the squares of its repeated factors divided out: a polynomial of R's sign at every point
% off the zero set of what was removed. Exact integer arithmetic (ratQ); on any refusal (overflow,
% a gcd whose square does not divide) the reduction so far is kept, which is always sign-correct.
    S = R;
    try
        S = ratQ.contentFree(ratQ.trimN(R));
        for it = 1:4
            G = 1;
            for d = 1:2
                Sd = diffPoly(S, d);
                if any(Sd(:)), G = ratQ.trimN(ratQ.polyGcdBiv(S, Sd)); end
                if numel(G) > 1, break, end
            end
            if numel(G) <= 1, return, end
            S = ratQ.contentFree(ratQ.trimN(ratQ.polyDivExact(S, ratQ.polyMulN(G, G), 2)));
        end
    catch
    end
end

function S = divideOut(S, factors)
% S divided by every factor (>= 0 on the region) as often as it divides EXACTLY: S's sign wherever the
% factors are positive. A dyadic factor is scaled to integers first (a positive scale keeps its sign);
% one that is constant, not dyadic, or does not divide is skipped.
    for k = 1:numel(factors)
        g = factors{k};
        for sc = 1:30
            if all(g(:) == round(g(:))), break, end
            g = 2 * g;
        end
        if numel(g) <= 1 || any(g(:) ~= round(g(:))), continue, end
        try
            g = ratQ.contentFree(ratQ.trimN(g));
            if numel(g) <= 1, continue, end
            if ~mayDivide(S, g), continue, end
            for it = 1:4
                S = ratQ.trimN(ratQ.polyDivExact(S, g, 2));
            end
        catch
        end
    end
end

function tf = mayDivide(S, g)
% false only when g PLAINLY does not divide S: S is clearly nonzero at a point of g = 0 (complex
% points too -- divisibility holds there as well). A filter in doubles before the exact division
% (2026-11-01): on an erratum Example 2 overlay tuple 1.1 million exact divisions, nearly all
% failing, were 60% of the time. Wrong only toward "may", which costs the division, never a sign.
    tf = true;
    if isa(S, 'bigZ') || ~isnumeric(S), return, end
    T = S;  T(1) = 0;
    if ~any(T(:)), tf = (S(1) == 0); return, end     % a constant: only 0 has a non-constant divisor
    swap = size(g,2) < 2;                            % g has no y: solve for x instead
    if swap, g = g.';  S = S.'; end
    if size(g,2) < 2, return, end
    xs = [0.6180339887; -1.3247179572; 2.7182818285];
    for i = 1:numel(xs)
        x = xs(i);
        cy = (x .^ (0:size(g,1)-1)) * g;             % g(x, y) as a polynomial in y, ascending
        cy = cy(1:find(cy ~= 0, 1, 'last'));
        if numel(cy) < 2, continue, end
        ys = roots(fliplr(cy));
        for y = ys(:).'
            T = S .* (x .^ (0:size(S,1)-1)).' .* (y .^ (0:size(S,2)-1));
            if abs(sum(T(:))) > 1e-6 * max(sum(abs(T(:))), realmin), tf = false; return, end
        end
    end
end

function D = diffPoly(S, d)
% dS/dx (d = 1) or dS/dy (d = 2) of a dense polyN, S(i+1,j+1) the coefficient of x^i y^j.
    if d == 2, D = diffPoly(S.', 1).'; return, end
    if size(S,1) < 2, D = zeros(1, size(S,2)); return, end
    D = S(2:end,:) .* (1:size(S,1)-1).';
end

function S = padAdd(A, B)
    n = max(size(A), size(B));
    S = zeros(n);
    if isa(A, 'bigZ') || isa(B, 'bigZ'), S = bigZ(S); end
    S(1:size(A,1), 1:size(A,2)) = A;
    S(1:size(B,1), 1:size(B,2)) = S(1:size(B,1), 1:size(B,2)) + B;
end
