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
%
% FACE AGAINST FACE (2026-11-02). A row may instead carry .W = struct('FP', 'FB'), another implicit
% face, quadratic in z like a's: the row is z_k > z_w (z_k = a's branch), and failing everywhere
% means z_k <= z_w wherever both exist. algFaceBelow uses it where no fitted quadratic separates
% the two -- they TOUCH along a curve no quadratic follows (erratum Example 2's pair cells). With
% t = z_w = (-B_w + e sqrt D_w) / (2 A_w), e = sign A_w for the upper root and -sign A_w for the
% lower, the quadOK disjunctions above need sign(t - m_k) and sign(F_k(t) / A_k), and
%       t - m_k  = (P1 + e A_k sqrt D_w) / (2 A_w A_k),     P1 = A_w B_k - A_k B_w
%       F_k(t)   = (P2 + e P1  sqrt D_w) / (2 A_w^2),       P2 = 2 A_w (A_w C_k - A_k C_w) - P1 B_w
% (A_w t^2 = -B_w t - C_w reduces F_k(t) to degree 1 in t). sign(P + Q sqrt D) is sign P when P, Q
% agree or E = P^2 - Q^2 D > 0, sign Q when E < 0. Where the faces touch, F_k(t) has a double zero
% and E2 = P2^2 - P1^2 D_w carries it SQUARED -- divided out by squareFreeSign, or, past 2^53 (on
% Example 2 E2 is degree ~20 with coefficients ~1e23), by an exact factorisation
% (symSquareFreeSign). A branch missing on a box (D < 0 there) fails the row by itself.
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
        % A vertex computed in doubles lands a few ulps off its dyadic value (erratum Example 2:
        % 47.5 stored as 47.49999999999999, and every blow-up was then refused as not dyadic).
        % Snapped: the certificate is exact about the snapped point, so a wrong snap only fails.
        q = round(p * 2^8) / 2^8;
        if all(abs(q - p) <= 1e-9 * max(1, abs(p))), p = q; end
        if p(1) < box(1) || p(1) > box(2) || p(2) < box(3) || p(2) > box(4), continue, end
        if ~anyVanishes(rows, p), continue, end      % nothing to blow up: G(0, t) = P(p), as plain
        for e = W ./ 2.^(2:8)
            if cornerOK(rows, p, e, factors)
                % around an IRRATIONAL corner the square is certified about the exact point, which
                % the double p misses by ~1e-16: record it a hair smaller
                if ~isDyadic(p), e = e * (1 - 1e-9); end
                sq(end+1,:) = [p, e]; break %#ok<AGROW>
            end
        end
    end
    if ~isempty(sq), tf = boxRec(rows, box, inBox, depth, sq); end
end

function r = rowPolys(a, factors)
% The polynomials whose signs decide row a on a box, and the test reading them; [] if the row
% cannot be certified at all (shape, branch index, non-finite coefficients).
    r = [];
    persistent tmemo
    if isfield(a, 'W') && ~isempty(a.W), r = faceRowPolys(a, factors); return, end
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
    if nz - 1 >= 3
        % THE TOP ROOT OF A POLYNOMIAL OF DEGREE n >= 3 in z (2026-11-04): erratum Example 2's
        % pair 2 has a bridge face quartic in z, and no row of it was ever formed -- a quadratic
        % below it was never shown empty. Sufficient conditions, by Taylor at q (G_n = lead in z):
        % q < z_top when sign(G_n) G(q) < 0 (a root above q); q > z_top when sign(G_n) G^(j)(q) > 0
        % for every j < n (G has no root at or above q). Each Dn^(n-j) G^(j)(Q / Dn), exact.
        if ~isequal(a.FB, -1) || Dn ~= round(Dn), return, end
        if isempty(tmemo), tmemo = containers.Map('KeyType', 'char', 'ValueType', 'any'); end
        key = rowKey(FP, num, a.den);
        if ~isempty(key), key = sprintf('top%d|%s', want, key); end
        if ~isempty(key) && isKey(tmemo, key)
            P = tmemo(key);
        else
            P = topRowPolys(FP, Q, Dn, nz - 1, want);
            if ~isempty(key)
                if tmemo.Count > 5000, tmemo.remove(tmemo.keys()); end
                tmemo(key) = P;
            end
        end
        if isempty(P), return, end
        P = cellfun(@(p) divideOut(p, factors), P, 'UniformOutput', false);
        r = struct('P', {P}, 'test', @(s) s(1) ~= 0 && all(s(1) * s(2:end) == want));
        return
    end
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

function P = topRowPolys(FP, Q, Dn, n, want)
% {G_n, H_0} (want -1) or {G_n, H_0, ..., H_(n-1)} (want +1), H_j = Dn^(n-j) G^(j)(Q / Dn) by Horner,
% exact (zMul/zLin), back in doubles as fromZ does, squared factors divided out; {} on failure.
    P = {};
    try
        if ratQ.isBigPoly(FP) || max(abs(FP(:))) >= 2^53, Z = toZ(FP); else, Z = FP; end
        G = cell(1, n+1);
        for i = 0:n, G{i+1} = zNorm(Z(:,:,i+1)); end
        D = java.math.BigInteger.valueOf(Dn);
        nj = 0;  if want == 1, nj = n - 1; end
        P = cell(1, nj + 2);
        P{1} = fromZ(bigZ(G{n+1}));
        for j = 0:nj
            c = @(i) factorial(i) / factorial(i - j);
            R = zLin(c(n), G{n+1}, 0, 0);
            for i = n-1:-1:j
                R = zLin(1, zMul(R, Q), c(i), zLin(D.pow(n - i), G{i+1}, 0, 0));
            end
            R = fromZ(bigZ(R));
            if ~all(isfinite(R(:))), P = {}; return, end
            P{j+2} = squareFreeSign(R);
        end
    catch
        P = {};
    end
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
            s(k) = polySign(P{k}, b);
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
% P(p) = 0 to rounding (a dense polyN in x, y). A polynomial over Q(sqrt k), U + sqrt(k) V (a
% SURD FACE ROW's), vanishes at a rational p only where U and V both do.
    if isstruct(P), tf = vanishesAt(P.U, p) && vanishesAt(P.V, p); return, end
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
            % A row some polynomial of which cannot be blown up is LEFT OUT of the charts (2026-11-03),
            % not the corner: every point needs only SOME row to fail, so fewer rows is a stronger
            % claim -- sound. On erratum Example 2 the cell's own rows have big coefficients and
            % refused at every dyadic corner, which refused the corner whatever the face row did.
            rowsC = rows;  keep = true(1, numel(rows));
            for t = 1:numel(rows)
                if isfield(rows{t}, 'corner') && ~rows{t}.corner, keep(t) = false; continue, end
                for k = 1:numel(rows{t}.P)
                    G = blowUp(rows{t}.P{k}, p, sg, major);
                    if isempty(G), keep(t) = false; break, end
                    rowsC{t}.P{k} = G;
                end
            end
            rowsC = rowsC(keep);
            if isempty(rowsC), return, end
            facC = cell(1, numel(factors));
            for k = 1:numel(factors)
                % only the region's rows THROUGH p: one positive at p rejects nothing near it, and
                % leaving a factor out only keeps more chart boxes -- safe. On an Example 2 overlay
                % tuple the charts tested ~60 factors per box and that was 49 of their 53 s.
                if ~vanishesAt(factors{k}, p), continue, end
                facC{k} = blowUp(factors{k}, p, sg, major);
            end
            facC = facC(~cellfun(@isempty, facC));
            inC = @(cb) ~any(cellfun(@(G) polySign(G, cb) == -1, facC));
            if ~boxRec(rowsC, [0 e -1 1], inC, 6, zeros(0,3)), return, end
        end
    end
    ok = true;
end

function G = blowUp(P, p, sg, major)
% P(p + sg s (1, t)) (major 1) or P(p + sg s (t, 1)) (major 2) divided by the largest power of s that
% divides it EXACTLY, as a polyN in (s, t); [] when the shift could not be done exactly in doubles.
    G = [];
    if isstruct(P), G = blowUpPair(P, p, sg, major); return, end
    for it = 1:30                                    % a dyadic P scaled to integers (sign kept)
        if all(P(:) == round(P(:))), break, end
        P = 2 * P;
    end
    if any(P(:) ~= round(P(:))), return, end
    if ~isDyadic(p), G = blowUpSurd(P, p, sg, major); return, end
    sc = 1;                                          % p's denominator, a power of 2 up to 2^8
    while any(p * sc ~= round(p * sc)) && sc < 2^8, sc = 2 * sc; end
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

function G = blowUpPair(P, p, sg, major)
% blowUp of a SURD FACE ROW's U + sqrt(k) V at a DYADIC p, exactly (BigInteger): U and V shifted with
% the same positive scale and the power of s dividing BOTH divided out (one of them alone would
% change their ratio); the reduced norm and the reduced U, V (read for their signs only) each with
% their own power of s out. [] at a non-dyadic p or without the exact parts.
    G = [];
    if ~isfield(P, 'EZ') || ~isDyadic(p), return, end
    try
        GU = chartZ(P.UZ, p, sg, major);  GV = chartZ(P.VZ, p, sg, major);
        n = max(size(GU), size(GV));
        GU = padZ(GU, n);  GV = padZ(GV, n);
        k = 0;
        for i = 1:n(1)
            if ~(rowZero(GU, i) && rowZero(GV, i)), k = i; break, end
        end
        if k == 0, G = 0; return, end
        GE = trimS(chartZ(P.EZ, p, sg, major));
        GUs = trimS(chartZ(P.UsZ, p, sg, major));  GVs = trimS(chartZ(P.VsZ, p, sg, major));
        G = makeSurd(zNorm(GU(k:end, :)), zNorm(GV(k:end, :)), P.k, GE, GUs, GVs);
    catch
        G = [];
    end
end

function Z = trimS(Z)
% The largest power of s dividing a chart polynomial (rows = powers of s) divided out, the content
% too (both positive); 0 for the zero polynomial (its sign, 0, decides nothing).
    for i = 1:size(Z, 1)
        if ~rowZero(Z, i), Z = zNorm(contentFree(Z(i:end, :))); return, end
    end
    Z = 0;
end

function tf = rowZero(Z, i)
    tf = true;
    for j = 1:size(Z, 2)
        v = Z(i, j);
        if v.c{1}.signum() ~= 0, tf = false; return, end
    end
end

function Z = padZ(Z, n)
    if all(size(Z) == n), return, end
    W = bigZ(zeros(n));  W(1:size(Z,1), 1:size(Z,2)) = Z;  Z = W;
end

function G = chartZ(Z, p, sg, major)
% sc^deg Z(p + sg s (1, t)) (major 1) or (t, 1) (major 2), exactly, as a bigZ polyN in (s, t) --
% p = a / sc dyadic, the positive scale sc^deg keeping the sign. Not trimmed. Z exact (double or bigZ).
    Z = bigZ(Z);
    sc = 1;
    while any(p * sc ~= round(p * sc)) && sc < 2^8, sc = 2 * sc; end
    a = round(p * sc);
    [n1, n2] = size(Z);  deg = n1 + n2 - 2;
    Ps = bigZ(zeros(n1, n2));                        % coefficient of u^k v^l, u = x - px, v = y - py
    for k = 0:n1-1
        for l = 0:n2-1
            acc = bigZ(0);
            for i = k:n1-1
                for j = l:n2-1
                    zij = Z(i+1, j+1);
                    if zij.c{1}.signum() == 0, continue, end
                    m = bigZ(nchoosek(i, k) * nchoosek(j, l)) .* bigZ(a(1))^(i-k) .* bigZ(a(2))^(j-l) ...
                        .* bigZ(sc)^(k + l + deg - i - j);
                    acc = acc + zij .* m;
                end
            end
            Ps(k+1, l+1) = acc;
        end
    end
    G = bigZ(zeros(n1 + n2 - 1, max(n1, n2)));
    for k = 0:n1-1
        for l = 0:n2-1
            c = Ps(k+1, l+1);
            if sg^(k+l) < 0, c = -c; end
            if major == 1, G(k+l+1, l+1) = G(k+l+1, l+1) + c;
            else,          G(k+l+1, k+1) = G(k+l+1, k+1) + c;
            end
        end
    end
end

function S = surdPoly(U, V, k)
% U + sqrt(k) V for polySign, from EXACT integer polyN U, V (doubles within 2^53, or bigZ), a joint
% positive content divided out. Its norm E = U^2 - k V^2 and U, V themselves are also kept with their
% EVEN-multiplicity factors removed (oddPartZ), for their signs: where the faces TOUCH along a
% curve, U + sqrt(k) V = H^2 M, and E carries norm(H)^2 (U and V carry H^2 when H is rational) --
% only the reduced ones have a sign beside that curve.
    n = max(size(U), size(V));
    U = zPad(U, n);  V = zPad(V, n);
    g = zGcd(zGcd(0, U), V);
    if g > 1, U = zDivScalar(U, g);  V = zDivScalar(V, g); end
    E = zLin(1, zMul(U, U), -k, zMul(V, V));
    S = makeSurd(U, V, k, oddPartZ(E), oddPartZ(U), oddPartZ(V));
end

function S = makeSurd(UZ, VZ, k, EZ, UsZ, VsZ)
% The struct polySign and blowUpPair read: exact parts (*Z), their doubles (rounded past 2^53, see
% BIG INTEGERS), and Lo, Hi = U + r V at the two doubles r bracketing sqrt(k).
    U = double(UZ);  V = double(VZ);  E = fromZ(EZ);
    if ~any(E(:)), E = []; end                       % U = -sqrt(k) V is impossible unless both vanish
    rd = sqrt(k);
    S = struct('U', U, 'V', V, 'E', E, 'Us', fromZ(UsZ), 'Vs', fromZ(VsZ), ...
               'Lo', U + rd * (1 - 4*eps) * V, 'Hi', U + rd * (1 + 4*eps) * V, ...
               'UZ', UZ, 'VZ', VZ, 'EZ', EZ, 'UsZ', UsZ, 'VsZ', VsZ, 'k', k);
end

function O = oddPartZ(Z)
% sign(c) times the product of the ODD-multiplicity factors of Z = c prod f_i^m_i -- Z's sign off the
% removed factors' zeros -- exact: Yun in doubles (oddPart) within 2^53, an exact factorisation
% (Symbolic Math Toolbox) past it. Z itself on any failure (sign-correct, only less decisive).
    O = Z;
    if isnumeric(Z)
        try
            D = ratQ.trimN(Z);
            if numel(D) <= 1, O = sign(D); return, end
            R = ratQ.contentFree(ratQ.trimN(oddPart(D, 1, 0)));
            sg = relSign(D, R);
            if sg ~= 0, O = sg * R; return, end
        catch
        end
        Z = bigZ(Z);                                 % Yun's gcds left 2^53: factor exactly instead
    end
    try
        x = sym('x');  y = sym('y');
        [n1, n2] = size(Z);
        C = sym(zeros(n1, n2));
        for i = 1:n1
            for j = 1:n2
                v = Z(i, j);  v = v.c{1};
                if v.signum() ~= 0, C(i,j) = str2sym(char(v.toString())); end
            end
        end
        p = (x.^(0:n1-1)) * C * (y.^(0:n2-1)).';
        if p == 0, return, end
        F = factor(p);
        keys = arrayfun(@char, F, 'UniformOutput', false);
        [~, first, j] = unique(keys);
        m = accumarray(j(:), 1);
        q = sym(1);  cst = sym(1);
        for i = 1:numel(first)
            f = F(first(i));
            if isempty(symvar(f)), cst = cst * f^m(i); continue, end
            if mod(m(i), 2) == 1, q = q * f; end
        end
        q = expand(sign(cst) * q);
        R = zNorm(symToZ(q, x, y));
        % the identity Z = (unit) O * square is factor's; checked at a point as relSign does
        if relSign(double(Z), double(R)) == 1, O = R; end
    catch
    end
    O = zNorm(O);
end

% ---- exact integer polyN: doubles while every partial sum is within 2^53, bigZ past it ----------

function Z = zNorm(Z)
% bigZ back to doubles when every entry fits exactly.
    if ~isa(Z, 'bigZ'), return, end
    c = Z.c;
    for i = 1:numel(c)
        if c{i}.abs().bitLength() > 52, return, end
    end
    Z = double(Z);
end

function C = zMul(A, B)
    if isnumeric(A) && isnumeric(B) && max(max(conv2(abs(A), abs(B)))) < 2^53, C = conv2(A, B); return, end
    C = zNorm(conv2g(bigZ(A), bigZ(B)));
end

function C = zLin(a, A, b, B)
% a A + b B, scalars a, b exact doubles or BigInteger, padded to the larger size.
    if isa(a, 'java.math.BigInteger'), a = zScalar(a); end
    if isa(b, 'java.math.BigInteger'), b = zScalar(b); end
    n = max(size(A), size(B));
    if isnumeric(A) && isnumeric(B) && isnumeric(a) && isnumeric(b)
        Ap = zPad(A, n);  Bp = zPad(B, n);
        if max(max(abs(a) * abs(Ap) + abs(b) * abs(Bp))) < 2^53, C = a * Ap + b * Bp; return, end
    end
    C = zNorm(zPad(bigZ(a) .* bigZ(A), n) + zPad(bigZ(b) .* bigZ(B), n));
end

function s = zScalar(v)
% a BigInteger as an exact double when it fits, else as a 1x1 bigZ.
    if v.abs().bitLength() <= 52, s = v.doubleValue(); else, s = bigZ({v}); end
end

function Z = zPad(Z, n)
    if isa(Z, 'bigZ'), Z = padZ(Z, n); return, end
    if all(size(Z) == n), return, end
    W = zeros(n);  W(1:size(Z,1), 1:size(Z,2)) = Z;  Z = W;
end

function g = zGcd(g, Z)
% gcd of g (a double) and every entry of Z; -1 when Z is bigZ (then nothing is divided).
    if ~isnumeric(Z) || g < 0, g = -1; return, end
    for v = abs(Z(:)).'
        g = gcd(g, v);
    end
end

function Z = zDivScalar(Z, g)
    Z = Z / g;
end

function Z = symToZ(q, x, y)
% An integer sym polynomial in x, y as a bigZ polyN.
    if isempty(symvar(q))
        Z = bigZ({java.math.BigInteger(char(q))});  return
    end
    [cf, t] = coeffs(q, [x y]);
    Z = bigZ(zeros(1, 1));
    for i = 1:numel(t)
        dx = double(polynomialDegree(t(i), x));  dy = double(polynomialDegree(t(i), y));
        Z(dx+1, dy+1) = bigZ({java.math.BigInteger(char(cf(i)))});
    end
end

function Z = divZ(Z, g)
    c = Z.c;  gg = g;
    for i = 1:numel(c), c{i} = c{i}.divide(gg); end
    Z = bigZ(c);
end

function tf = isDyadic(p)
    tf = all(p * 2^8 == round(p * 2^8));
end

function v = polySign(P, b)
% The certified sign of P on box b (0 when undecided): P a dense polyN, or a blown-up polynomial
% over Q(sqrt d) (blowUpSurd's struct, U + sqrt(d) V).
    if isnumeric(P)
        v = polySignOn(P, b, @(~) true, 0);
        if isnan(v), v = 0; end
        return
    end
    % U + sqrt(d) V lies between U + r V at the two doubles r bracketing sqrt(d), pointwise
    lo = polySignOn(P.Lo, b, @(~) true, 0);
    if lo ~= 0 && ~isnan(lo) && polySignOn(P.Hi, b, @(~) true, 0) == lo, v = lo; return, end
    v = 0;
    if isempty(P.E), return, end
    if isfield(P, 'Us')                              % a surd face row's: reduced, for their signs
        sU = polySignOn(P.Us, b, @(~) true, 0);  sV = polySignOn(P.Vs, b, @(~) true, 0);
    else
        sU = polySignOn(P.U, b, @(~) true, 0);  sV = polySignOn(P.V, b, @(~) true, 0);
    end
    sE = polySignOn(P.E, b, @(~) true, 0);
    v = sqSign(sU, sV, sE);
end

function r = surdOf(p)
% p = ([a1 a2] + [b1 b2] sqrt(d)) / q with small integers, recognised from the doubles: r =
% [a1 b1 a2 b2 q d], or [] when none fits. A WRONG guess is sound: the blow-up is exact about the
% point guessed, so it only fails to certify (the corner left undecided).
    r = [];
    D = [2 3 5 6 7 10 11 13 14 15 17 19 21 22 23 26 29 30 31 33 34 35 37 38 39 41 42 43 46 47];
    bs = -64:64;
    for q = 1:64
        for d = D
            ab = zeros(2,2);  ok = true;
            for c = 1:2
                a = round(q * p(c) - bs * sqrt(d));
                err = abs(q * p(c) - bs * sqrt(d) - a);
                k = find(err < 1e-11 * q * max(1, abs(p(c))));
                if isempty(k), ok = false; break, end
                [~, j] = min(abs(bs(k)));  k = k(j);
                ab(c,:) = [a(k), bs(k)];
            end
            if ok && any(ab(:,2) ~= 0)
                r = [ab(1,:), ab(2,:), q, d];  return
            end
        end
    end
end

function G = blowUpSurd(P, p, sg, major)
% blowUp at a corner in Q(sqrt d)^2 (2026-11-03): the triangle's (1/sqrt5, 1 - 1/sqrt5). With
% p = (A, B) / q, A, B in Z[sqrt d], q^deg P(p + sg s (1, t)) = U + sqrt(d) V exactly, U, V integer
% polyN in (s, t), the powers of s that divide both exactly divided out; q^deg > 0 keeps the sign.
% [] when p is not recognised or a coefficient could pass 2^52.
    G = [];
    r = surdOf(p);
    if isempty(r), return, end
    A = r(1:2);  B = r(3:4);  q = r(5);  d = r(6);
    [n1, n2] = size(P);
    deg = n1 + n2 - 2;
    NA = abs(A(1)) + sqrt(d) * abs(A(2));  NB = abs(B(1)) + sqrt(d) * abs(B(2));
    bound = 0;
    for i = 0:n1-1
        for j = 0:n2-1
            bound = bound + abs(P(i+1,j+1)) * (NA + q)^i * (NB + q)^j * q^(deg-i-j);
        end
    end
    if bound >= 2^52, return, end                    % each component of every partial sum is below
    pA = surdPow(A, n1 - 1, d);  pB = surdPow(B, n2 - 1, d);
    Su = zeros(n1, n2);  Sv = zeros(n1, n2);          % coefficient of u^k v^l, u = x - px, v = y - py
    for k = 0:n1-1
        for l = 0:n2-1
            acc = [0 0];
            for i = k:n1-1
                for j = l:n2-1
                    if P(i+1,j+1) == 0, continue, end
                    m = surdMul(pA(i-k+1,:), pB(j-l+1,:), d);
                    acc = acc + P(i+1,j+1) * nchoosek(i,k) * nchoosek(j,l) * q^(k+l+deg-i-j) * m;
                end
            end
            Su(k+1,l+1) = acc(1);  Sv(k+1,l+1) = acc(2);
        end
    end
    U = zeros(n1 + n2 - 1, max(n1, n2));  V = U;
    for k = 0:n1-1
        for l = 0:n2-1
            f = sg^(k+l);
            if major == 1, U(k+l+1, l+1) = U(k+l+1, l+1) + f * Su(k+1,l+1);  V(k+l+1, l+1) = V(k+l+1, l+1) + f * Sv(k+1,l+1);
            else,          U(k+l+1, k+1) = U(k+l+1, k+1) + f * Su(k+1,l+1);  V(k+l+1, k+1) = V(k+l+1, k+1) + f * Sv(k+1,l+1);
            end
        end
    end
    k = find(any(U ~= 0, 2) | any(V ~= 0, 2), 1);
    if isempty(k), G = 0; return, end
    U = U(k:end, :);  V = V(k:end, :);
    E = [];                                          % U^2 - d V^2, when it is exact in doubles
    if max(max(conv2(abs(U), abs(U)) + d * conv2(abs(V), abs(V)))) < 2^53
        E = conv2(U, U) - d * conv2(V, V);
    end
    rd = sqrt(d);
    G = struct('U', U, 'V', V, 'E', E, 'Lo', U + rd * (1 - 4*eps) * V, 'Hi', U + rd * (1 + 4*eps) * V);
end

function M = surdPow(A, n, d)
% Rows: A^0 .. A^n in Z[sqrt d], as [rational part, sqrt(d) part].
    M = zeros(n + 1, 2);  M(1,:) = [1 0];
    for k = 1:n, M(k+1,:) = surdMul(M(k,:), A, d); end
end

function m = surdMul(x, y, d)
    m = [x(1)*y(1) + d*x(2)*y(2), x(1)*y(2) + x(2)*y(1)];
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

function r = faceRowPolys(a, factors)
% A FACE ROW: z_k > z_w, z_k = a's branch, z_w = a.W's, both quadratic in z. See FACE AGAINST FACE.
    r = [];
    Fk = a.FP;  Fw = a.W.FP;
    nk = nzOf(Fk);  nw = nzOf(Fw);
    % z_w THE TOP ROOT OF A POLYNOMIAL OF DEGREE >= 3 in z (2026-11-04): erratum Example 2's pair 2
    % compares every face it keeps with one bridge face quartic in z, and no row of either kind
    % below was built for it -- every one of those 89 comparisons failed by construction.
    if nw >= 4 && any(nk == [2 3]) && isequal(a.W.FB, -1)
        hiK = false;
        if nk == 3
            hiK = any(a.FB == [-1 2]);  if ~hiK && ~any(a.FB == [1 -2]), return, end
        elseif ~any(a.FB == [1 -1]), return
        end
        r = topFaceRow(Fk, Fw, nk, hiK, factors);
        return
    end
    % ONE FACE RATIONAL (2026-11-03): z = -C/B, degree 1 in z. The 3-convex-edge triangle's last
    % zero-area faces are cut off by a tie between a quadratic-fibre face and such a face, and this
    % returned [] for it -- the row was dropped and the tie never certified.
    if nk == 3 && nw == 2, mode = 'ql'; elseif nk == 2 && nw == 3, mode = 'lq';
    elseif nk == 3 && nw == 3, mode = 'qq'; else, return
    end
    hiK = false;  hiW = false;
    if nk == 3
        hiK = any(a.FB == [-1 2]);  if ~hiK && ~any(a.FB == [1 -2]), return, end
    elseif ~any(a.FB == [1 -1]), return
    end
    if nw == 3
        hiW = any(a.W.FB == [-1 2]);  if ~hiW && ~any(a.W.FB == [1 -2]), return, end
    elseif ~any(a.W.FB == [1 -1]), return
    end
    r = surdFaceRow(Fk, Fw, nk, nw, hiK, hiW, factors);
    if ~isempty(r), return, end
    persistent memo
    if isempty(memo), memo = containers.Map('KeyType', 'char', 'ValueType', 'any'); end
    k1 = rowKey(Fk, 0, 1);  k2 = rowKey(Fw, 0, 1);
    key = '';
    if ~isempty(k1) && ~isempty(k2), key = ['face|', mode, '|', k1, '#', k2]; end
    if ~isempty(key) && isKey(memo, key)
        P = memo(key);
    else
        P = faceSignPolys(Fk, Fw, mode);
        if isempty(P), return, end
        if ~isempty(key)
            if memo.Count > 2000, memo.remove(memo.keys()); end
            memo(key) = P;
        end
    end
    P = cellfun(@(p) divideOut(p, factors), P, 'UniformOutput', false);
    switch mode
        case 'qq', test = @(s) faceOK(s, hiK, hiW);
        case 'ql', test = @(s) s(2) ~= 0 && quadOK([s(1), s(2) * s(3), s(4), s(5)], hiK, 1);    % z_w >= z_k
        case 'lq', test = @(s) s(2) ~= 0 && quadOK([s(1), s(2) * s(3), s(4), s(5)], hiW, -1);   % z_k <= z_w
    end
    r = struct('P', {P}, 'test', test);
end

function r = surdFaceRow(Fk, Fw, nk, nw, hiK, hiW, factors)
% A SURD FACE ROW (2026-11-03): z_k > z_w when each face is a polynomial over ONE quadratic field --
% its discriminant D = c S^2, S in Z[x, y], c > 0 not a square -- or rational (degree 1 in z).
% Erratum Example 2's pair cells are made of such faces (every D = 55 times a square), and the
% general face row could not decide them at the corner where two touch: its E2 = P2^2 - P1^2 D_w is
% the product of z_w - z_k with its CONJUGATE, and the conjugate vanishes along the cell's edge,
% which leaves E2's lowest form at the corner -c L^2, degenerate along the edge L = 0 -- no chart
% decides the direction of the edge. z_w - z_k itself vanishes there only to order two, with a
% DEFINITE lowest form, so it is read directly: with z = (a0 + eta tau sqrt(c) S) / den
% (den = 2A, a0 = -B, eta = +1 for the upper root, tau = sign(A) sign(S) -- constant on a box where
% both are certified, the region's factors divided out of S first), and kappa the field's
% square-free part, sqrt(c) = (mu / kappa) sqrt(kappa) with mu = sqrt(c kappa) an integer:
%       kappa den_k den_w (z_w - z_k) = U + sqrt(kappa) V,
%       U = kappa (den_k a0_w - den_w a0_k),  V = tau_w eta_w mu_w den_k S_w - tau_k eta_k mu_k den_w S_k
% (mu = 0 for a rational face), one (U, V) per (tau_k, tau_w). The row fails where
% sign(U + sqrt(kappa) V) sign(den_k den_w) = +1. [] when the faces are not of this kind.
    r = [];
    persistent memo
    if isempty(memo), memo = containers.Map('KeyType', 'char', 'ValueType', 'any'); end
    key = ['surd|', rowKey(Fk, 0, 1), '#', rowKey(Fw, 0, 1), sprintf('|%d%d', hiK, hiW)];
    if isKey(memo, key)
        P = memo(key);
        if isempty(P), return, end
        P(1:4) = cellfun(@(p) divideOut(p, factors), P(1:4), 'UniformOutput', false);
        P(5:8) = cellfun(@(T) divideOutSurd(T, factors), P(5:8), 'UniformOutput', false);
        r = struct('P', {P}, 'test', @surdOK);
        return
    end
    if memo.Count > 2000, memo.remove(memo.keys()); end
    memo(key) = {};
    try
        fk = surdFace(Fk, nk, hiK);  fw = surdFace(Fw, nw, hiW);
        if isempty(fk) || isempty(fw) || (fk.lin && fw.lin), return, end
        if ~fw.lin, kap = squareFreeInt(fw.c); else, kap = squareFreeInt(fk.c); end
        if kap.bitLength() > 52, return, end
        muK = fieldMu(fk, kap);  muW = fieldMu(fw, kap);
        if isempty(muK) || isempty(muW), return, end
        k = kap.doubleValue();
        U = zLin(k, zLin(1, zMul(fk.den, fw.a0), -1, zMul(fw.den, fk.a0)), 0, 0);
        Wk = zMul(fw.den, fk.S);  Ww = zMul(fk.den, fw.S);         % den_w S_k, den_k S_w
        T = cell(1, 4);
        for tk = [1 -1]
            for tw = [1 -1]
                V = zLin(scaleBy(muW, tw * fw.eta), Ww, scaleBy(muK, -tk * fk.eta), Wk);
                T{(tk < 0) * 2 + (tw < 0) + 1} = surdPoly(U, V, k);
            end
        end
        cf = @(p) fromZ(bigZ(p));                      % content divided out (positive), as doubles
        P = [{cf(fk.den), cf(fw.den), cf(fk.S), cf(fw.S)}, T];
        memo(key) = P;
        P(1:4) = cellfun(@(p) divideOut(p, factors), P(1:4), 'UniformOutput', false);
        P(5:8) = cellfun(@(T) divideOutSurd(T, factors), P(5:8), 'UniformOutput', false);
        r = struct('P', {P}, 'test', @surdOK);
    catch
        r = [];
    end
end

function r = topFaceRow(Fk, Fw, nk, hiK, factors)
% A FACE ROW AGAINST A TOP ROOT (2026-11-04): z_k > z_w with z_w the LARGEST real root of G = Fw, of
% degree n >= 3 in z, and z_k a face over one quadratic field (surdFace) or rational. Sufficient,
% not equivalent: sign(G_n) G(z_k) < 0 puts a real root of G above z_k (G -> sign(G_n) inf), so
% z_k < z_w; the row fails there. With z_k = N / M, N = kappa a0 + eta tau mu S sqrt(kappa) = X + Y
% sqrt(kappa), M = kappa den (see surdFaceRow),
%       M^n G(N / M) = sum_i G_i N^i M^(n-i) = U + sqrt(kappa) V,
% by Horner over Z[sqrt kappa]: no conjugate of z_k enters, so where the two faces touch only
% z_w - z_k vanishes. One (U, V) per tau_k; the row fails where sign(G_n) sign(den)^n sign(T) = -1.
% Where z_k lies below two further roots of G the test reads +1 and only fails to certify.
    r = [];
    persistent memo
    if isempty(memo), memo = containers.Map('KeyType', 'char', 'ValueType', 'any'); end
    key = ['top|', rowKey(Fk, 0, 1), '#', rowKey(Fw, 0, 1), sprintf('|%d', hiK)];
    if isKey(memo, key)
        P = memo(key);
    else
        if memo.Count > 2000, memo.remove(memo.keys()); end
        memo(key) = {};  P = {};
        try
            fk = surdFace(Fk, nk, hiK);
            if ~isempty(fk)
                if fk.lin, kap = java.math.BigInteger.ONE; else, kap = squareFreeInt(fk.c); end
                mu = fieldMu(fk, kap);
                if kap.bitLength() <= 52 && ~isempty(mu)
                    k = kap.doubleValue();
                    if ratQ.isBigPoly(Fw) || max(abs(Fw(:))) >= 2^53, Z = toZ(Fw); else, Z = Fw; end
                    n = nzOf(Fw) - 1;
                    G = cell(1, n+1);
                    for i = 0:n, G{i+1} = zNorm(Z(:,:,i+1)); end
                    X = zLin(k, fk.a0, 0, 0);  M = zLin(k, fk.den, 0, 0);
                    T = cell(1, 2);
                    for tk = [1 -1]
                        Y = zLin(scaleBy(mu, tk * fk.eta), fk.S, 0, 0);
                        U = G{n+1};  V = 0;  Mp = 1;
                        for i = n-1:-1:0
                            Mp = zMul(Mp, M);
                            Un = zLin(1, zMul(U, X), k, zMul(V, Y));
                            V = zLin(1, zMul(U, Y), 1, zMul(V, X));
                            U = zLin(1, Un, 1, zMul(G{i+1}, Mp));
                        end
                        T{(tk < 0) + 1} = surdPoly(U, V, k);
                    end
                    cf = @(p) fromZ(bigZ(p));
                    P = [{cf(fk.den), cf(fk.S), cf(G{n+1})}, T, {mod(n, 2)}];
                    memo(key) = P;
                end
            end
        catch
            P = {};
        end
    end
    if isempty(P), return, end
    odd = P{end};  P = P(1:end-1);
    P(1:3) = cellfun(@(p) divideOut(p, factors), P(1:3), 'UniformOutput', false);
    P(4:5) = cellfun(@(T) divideOutSurd(T, factors), P(4:5), 'UniformOutput', false);
    % NOT blown up at corners (cornerOK leaves it out): its U, V are the quartic's Horner images, of
    % high degree with big coefficients, and the exact BigInteger charts of them ran for hours
    % on the poly2 kink (kinkSubdivisionTest, 2026-11-04). Leaving a row out is sound, as there.
    r = struct('P', {P}, 'test', @(s) topOK(s, odd), 'corner', false);
end

function ok = topOK(s, odd)
% s = signs of {den_k, S_k, G_n, T(+), T(-)}: z_k < z_w (the top root) on the box.
    ok = false;
    if any(s(1:3) == 0), return, end
    t = s(4 + (s(1) * s(2) < 0));
    dn = 1;  if odd, dn = s(1); end
    ok = t * s(3) * dn == -1;
end

function v = scaleBy(mu, sg)
% sg * mu for a BigInteger mu (sg = +-1).
    if sg < 0, v = mu.negate(); else, v = mu; end
end

function ok = surdOK(s)
% s = signs of {den_k, den_w, S_k, S_w, T(++), T(+-), T(-+), T(--)}: z_k <= z_w on the box.
    ok = false;
    if any(s(1:4) == 0), return, end
    tk = s(1) * s(3);  tw = s(2) * s(4);
    t = s(4 + (tk < 0) * 2 + (tw < 0) + 1);
    ok = t ~= 0 && t * s(1) * s(2) == 1;
end

function f = surdFace(F, n, hi)
% A face as z = (a0 + eta tau sqrt(c) S) / den (see surdFaceRow): exact den, a0, S, BigInteger c,
% eta; a rational face has S = 1, c = 0. [] when D is not c S^2 (memoised per face).
    persistent memo
    if isempty(memo), memo = containers.Map('KeyType', 'char', 'ValueType', 'any'); end
    key = [rowKey(F, 0, 1), sprintf('|%d', hi)];
    if isKey(memo, key), f = memo(key); return, end
    f = [];
    if ratQ.isBigPoly(F) || max(abs(F(:))) >= 2^53, Z = toZ(F); else, Z = F; end
    sl = @(i) zNorm(Z(:,:,i));
    if n == 2
        f = struct('lin', true, 'den', sl(2), 'a0', zLin(-1, sl(1), 0, 0), 'S', 1, ...
                   'c', java.math.BigInteger.ZERO, 'eta', 1);
    else
        A = sl(3);  B = sl(2);  C = sl(1);
        D = zLin(1, zMul(B, B), -4, zMul(A, C));
        g = zGcd(0, D);
        if g > 0
            c = java.math.BigInteger.valueOf(g);  Dp = D / g;
        else
            [Dp, c] = contentFree(bigZ(D));  Dp = zNorm(Dp);
        end
        S = zSqrt(Dp);
        if ~isempty(S) && c.signum() > 0
            f = struct('lin', false, 'den', zLin(2, A, 0, 0), 'a0', zLin(-1, B, 0, 0), 'S', S, 'c', c, ...
                       'eta', 2 * hi - 1);
        end
    end
    if memo.Count > 2000, memo.remove(memo.keys()); end
    memo(key) = f;
end

function mu = fieldMu(f, kap)
% sqrt(c kappa), an integer when c is in kappa's square class; 0 for a rational face; [] otherwise.
    if f.lin, mu = java.math.BigInteger.ZERO; return, end
    mu = bigIsqrt(f.c.multiply(kap));
end

function s = bigIsqrt(n)
% The exact integer square root of a BigInteger, [] when n is not a perfect square.
    s = [];
    if n.signum() < 0, return, end
    if n.signum() == 0, s = n; return, end
    x = java.math.BigInteger.ONE.shiftLeft(floor(n.bitLength() / 2) + 1);
    while true
        y = x.add(n.divide(x)).shiftRight(1);
        if y.compareTo(x) >= 0, break, end
        x = y;
    end
    if x.multiply(x).equals(n), s = x; end
end

function k = squareFreeInt(c)
% The square-free part of a positive BigInteger: trial division by small primes, then the rest
% kept whole unless it is a perfect square. A non-minimal k is still a valid field generator (c k
% stays a square for every c of the field only when k IS minimal up to squares -- so a k left
% with a square factor makes fieldMu refuse, which only refuses the row).
    k = java.math.BigInteger.ONE;  r = c;
    for p = primes(5000)
        P = java.math.BigInteger.valueOf(p);
        e = 0;
        while r.mod(P).signum() == 0, r = r.divide(P); e = e + 1; end
        if mod(e, 2) == 1, k = k.multiply(P); end
        if r.equals(java.math.BigInteger.ONE), return, end
    end
    if isempty(bigIsqrt(r)), k = k.multiply(r); end
end

function S = zSqrt(Dp)
% S in Z[x, y], content-free, with S^2 = Dp exactly (Dp a content-free polyN in doubles), or []: Yun's
% square-free decomposition (yunParts, exact in doubles), every part of EVEN multiplicity, and
% S^2 = Dp checked. Refused past 2^53 (bigZ) -- only refuses the surd row.
    S = [];
    if ~isnumeric(Dp), return, end
    try
        P = ratQ.trimN(Dp);
        if numel(P) <= 1
            if P > 0 && sqrt(P) == round(sqrt(P)), S = sqrt(P); end
            return
        end
        L = yunParts(P, 1, 0);
        Q = 1;
        for i = 1:size(L, 1)
            if mod(L{i,2}, 2) == 1, return, end
            for k = 1:L{i,2} / 2, Q = ratQ.trimN(ratQ.polyMulN(Q, L{i,1})); end
        end
        Q = ratQ.contentFree(ratQ.trimN(Q));
        Q2 = ratQ.trimN(ratQ.polyMulN(Q, Q));
        if isequal(size(Q2), size(P)) && all(Q2(:) == P(:)), S = Q; end
    catch
        S = [];
    end
end

function L = yunParts(S, d, lvl)
% Yun's square-free decomposition, S = c prod Z_i^i: rows {Z_i, i} of the non-constant parts, in x
% (d = 1), then what is left free of x in y (d = 2) -- as oddPart, keeping every part.
    L = cell(0, 2);
    S = ratQ.trimN(S);
    if numel(S) <= 1 || lvl > 2, return, end
    Sd = diffPoly(S, d);
    if ~any(Sd(:))
        if d == 1, L = yunParts(S, 2, lvl + 1); end
        return
    end
    C = ratQ.trimN(ratQ.polyGcdBiv(S, Sd));
    W = ratQ.trimN(ratQ.polyDivExact(S, C, 2));
    i = 1;
    while numel(W) > 1
        if i > 60, error('algRowFails:yun', 'square-free decomposition did not end'); end
        Y = ratQ.trimN(ratQ.polyGcdBiv(W, C));
        Z = ratQ.trimN(ratQ.polyDivExact(W, Y, 2));
        if numel(Z) > 1, L(end+1, :) = {Z, i}; end %#ok<AGROW>
        C = ratQ.trimN(ratQ.polyDivExact(C, Y, 2));
        W = Y;  i = i + 1;
    end
    if d == 1 && numel(C) > 1, L = [L; yunParts(C, 2, lvl + 1)]; end
end

function n = nzOf(FP)
    if ratQ.isBigPoly(FP), sz = FP.sz; else, sz = size(FP); end
    sz(end+1:3) = 1;  n = sz(3);
end

function P = faceSignPolys(Fk, Fw, mode)
% {A_w, A_k, D_w, D_k, P1, E1, P2, E2} with their squared factors divided out; {} on failure.
% Exactly in doubles when every partial sum is bounded by 2^53, else in bigZ; the square-free
% reduction of a polynomial past 2^53 is then SYMBOLIC (symSquareFreeSign) -- the tangency's G^2
% is a factor of E2, and on erratum Example 2 E2 has degree ~20 and coefficients ~1e23.
    if nargin < 3, mode = 'qq'; end
    big = ratQ.isBigPoly(Fk) || ratQ.isBigPoly(Fw);
    if ~big
        bound = faceRaw(abs(Fk), abs(Fw), 1, mode);
        big = any(cellfun(@(p) any(p(:) >= 2^53), bound));
    end
    try
        if big
            Z = faceRaw(toZ(Fk), toZ(Fw), -1, mode);
            P = cell(size(Z));
            for i = 1:numel(Z)
                Z{i} = contentFree(Z{i});
                try
                    P{i} = exactDouble(Z{i}, 'face row polynomial');
                    P{i} = squareFreeSign(P{i});
                catch
                    P{i} = symSquareFreeSign(Z{i});
                end
            end
        else
            P = faceRaw(Fk, Fw, -1, mode);
            P = cellfun(@squareFreeSign, P, 'UniformOutput', false);
        end
    catch
        P = {};  return
    end
    if any(cellfun(@(p) isempty(p) || ~all(isfinite(p(:))), P)), P = {}; end
end

function P = faceRaw(Fk, Fw, sg, mode)
% The raw polynomials; sg = -1 for the true ones, +1 (with |inputs|) for a bound on every partial sum.
% mode 'ql' (w rational, z_w = -C_w/B_w) / 'lq' (k rational): {A, B_lin, B B_lin - 2 A C_lin,
% A C_lin^2 - B C_lin B_lin + C B_lin^2, D} of the quadratic one -- q = z_lin in quadOK, scaled by
% B_lin (the second, sign carried by B_lin) and B_lin^2 (the third, positive).
    if nargin < 4, mode = 'qq'; end
    if ~strcmp(mode, 'qq')
        if strcmp(mode, 'ql'), Q = Fk; L = Fw; else, Q = Fw; L = Fk; end
        A = Q(:,:,3);  B = Q(:,:,2);  C = Q(:,:,1);  Bl = L(:,:,2);  Cl = L(:,:,1);
        G = padAdd(conv2g(B, Bl), 2*sg * conv2g(A, Cl));
        H = padAdd(padAdd(conv2g(A, conv2g(Cl, Cl)), sg * conv2g(B, conv2g(Cl, Bl))), conv2g(C, conv2g(Bl, Bl)));
        D = padAdd(conv2g(B, B), 4*sg * conv2g(A, C));
        P = {A, Bl, G, H, D};
        return
    end
    Ak = Fk(:,:,3);  Bk = Fk(:,:,2);  Ck = Fk(:,:,1);
    Aw = Fw(:,:,3);  Bw = Fw(:,:,2);  Cw = Fw(:,:,1);
    Dw = padAdd(conv2g(Bw, Bw), 4*sg * conv2g(Aw, Cw));
    Dk = padAdd(conv2g(Bk, Bk), 4*sg * conv2g(Ak, Ck));
    P1 = padAdd(conv2g(Aw, Bk), sg * conv2g(Ak, Bw));
    E1 = padAdd(conv2g(P1, P1), sg * conv2g(conv2g(Ak, Ak), Dw));
    P2 = padAdd(sg * conv2g(P1, Bw), 2 * conv2g(Aw, padAdd(conv2g(Aw, Ck), sg * conv2g(Ak, Cw))));
    E2 = padAdd(conv2g(P2, P2), sg * conv2g(conv2g(P1, P1), Dw));
    P = {Aw, Ak, Dw, Dk, P1, E1, P2, E2};
end

function S = symSquareFreeSign(Z)
% An integer polynomial past 2^53 (bigZ) with every factor of EVEN multiplicity removed and every
% odd one kept once, by an exact factorisation (Symbolic Math Toolbox `factor`): the original's
% sign off the removed factors' zero sets. Back in doubles, rounded when past 2^53 (sound for
% polySignOn, see BIG INTEGERS). Any failure returns the polynomial unreduced -- sign-correct,
% only less decisive. The sign is checked against the original at a point (relSign).
    S = double(Z);
    try
        x = sym('x');  y = sym('y');
        [n1, n2] = size(Z);
        C = sym(zeros(n1, n2));
        for i = 1:n1
            for j = 1:n2
                v = Z.c{i,j};
                if v.signum() ~= 0, C(i,j) = str2sym(char(v.toString())); end
            end
        end
        p = (x.^(0:n1-1)) * C * (y.^(0:n2-1)).';
        F = factor(p);                               % irreducible factors, repeated by multiplicity
        keys = cell(1, numel(F));
        for i = 1:numel(F), keys{i} = char(F(i)); end
        [u, first, j] = unique(keys);
        m = accumarray(j(:), 1);
        q = sym(1);
        for i = 1:numel(u)
            f = F(first(i));
            if isempty(symvar(f)), continue, end    % constants: the sign is fixed by relSign
            if mod(m(i), 2) == 1, q = q * f; end
        end
        q = expand(q);
        if isempty(symvar(q)), S = 1; return, end
        [cf, t] = coeffs(q, [x y]);
        R = zeros(1, 1);
        for i = 1:numel(t)
            dx = double(polynomialDegree(t(i), x));  dy = double(polynomialDegree(t(i), y));
            R(dx+1, dy+1) = double(cf(i));
        end
        sg = relSign(S, R);
        if sg ~= 0, S = sg * R; end
    catch
    end
end

function ok = faceOK(s, hiK, hiW)
% s = signs of {A_w, A_k, D_w, D_k, P1, E1, P2, E2} on the box (0 where not certified). True when
% z_k <= z_w at every point of the box where both branches exist (off the removed curves).
    ok = false;
    if s(3) == -1 || s(4) == -1, ok = true; return, end          % a branch exists nowhere here
    if s(1) == 0 || s(2) == 0, return, end
    e = s(1);  if ~hiW, e = -e; end                              % z_w = (-B_w + e sqrt D_w) / (2 A_w)
    u = s(1) * s(2) * sqSign(s(5), e * s(2), s(6));              % sign(z_w - m_k)
    v = s(2) * sqSign(s(7), e * s(5), s(8));                     % sign(F_k(z_w) / A_k)
    if hiK, ok = u == 1 && v == 1;                               % z_w >= the upper root
    else,   ok = u == 1 || v == -1;                              % z_w >= the lower root
    end
end

function r = sqSign(sP, sQ, sE)
% sign(P + Q sqrt(D)) from sign P, sign Q and sign E = P^2 - Q^2 D; 0 when they do not decide it.
    r = 0;
    if sP ~= 0 && sP == sQ, r = sP;
    elseif sE == 1 && sP ~= 0, r = sP;
    elseif sE == -1 && sQ ~= 0, r = sQ;
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
%
% THE ODD PART FIRST (2026-11-02): Yun's square-free decomposition, R = c prod Z_i^i, and the
% product of the Z_i of ODD i -- R's sign times sign c off the Z_i's zeros. The plain "divide by
% gcd(R, R')^2" below cannot reduce a factor of multiplicity 3 or more (gcd^2 then has too much of
% it), and two faces TOUCHING along a curve give exactly that: (x - y)^4 on the synthetic pair of
% algFaceBelowTest. The plain loop stays as the fallback when Yun's steps refuse.
    S = R;
    try
        S0 = ratQ.contentFree(ratQ.trimN(R));
        if numel(S0) <= 1, S = S0; return, end
        O = oddPart(S0, 1, 0);
        O = ratQ.contentFree(ratQ.trimN(O));
        sg = relSign(S0, O);
        if sg ~= 0, S = sg * O; return, end
    catch
    end
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

function O = oddPart(S, d, lvl)
% The product of the odd-multiplicity square-free factors of S, up to a constant: Yun's algorithm
% with the derivative in x (d = 1); what is left free of x is then decomposed in y (d = 2).
    O = 1;
    S = ratQ.trimN(S);
    if numel(S) <= 1 || lvl > 2, return, end
    Sd = diffPoly(S, d);
    if ~any(Sd(:))
        if d == 1, O = oddPart(S, 2, lvl + 1); end
        return
    end
    C = ratQ.trimN(ratQ.polyGcdBiv(S, Sd));
    W = ratQ.trimN(ratQ.polyDivExact(S, C, 2));
    i = 1;
    while numel(W) > 1
        if i > 60, error('algRowFails:yun', 'square-free decomposition did not end'); end
        Y = ratQ.trimN(ratQ.polyGcdBiv(W, C));
        Z = ratQ.trimN(ratQ.polyDivExact(W, Y, 2));
        if mod(i, 2) == 1, O = ratQ.trimN(ratQ.polyMulN(O, Z)); end
        C = ratQ.trimN(ratQ.polyDivExact(C, Y, 2));
        W = Y;  i = i + 1;
    end
    if d == 1 && numel(C) > 1, O = ratQ.trimN(ratQ.polyMulN(O, oddPart(C, 2, lvl + 1))); end
end

function sg = relSign(S, O)
% sign(S / O), constant off their zeros (S / O = c times a square); read at a point where both are
% clearly nonzero, 0 when no such point is found among the fixed ones.
    sg = 0;
    pts = [0.6180339887 0.3247179572; -0.7182818285 0.4142135624; 0.2360679775 -0.8660254038; ...
           1.7320508076 1.4142135624; -1.2599210499 -0.5772156649; 3.1415926536 -2.7182818285];
    for k = 1:size(pts, 1)
        [vs, ms] = evalAt(S, pts(k,:));  [vo, mo] = evalAt(O, pts(k,:));
        if abs(vs) > 1e-6 * ms && abs(vo) > 1e-6 * mo
            sg = sign(vs) * sign(vo);  return
        end
    end
end

function [v, m] = evalAt(P, p)
    T = P .* (p(1) .^ (0:size(P,1)-1)).' .* (p(2) .^ (0:size(P,2)-1));
    v = sum(T(:));  m = max(sum(abs(T(:))), realmin);
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

function T = divideOutSurd(T, factors)
% divideOut for U + sqrt(k) V (makeSurd's struct, 2026-11-04): U and V JOINTLY by each region factor
% as often as it divides both exactly (their ratio, which the bracket Lo/Hi reads, is kept), and the
% norm and the reduced U, V -- read for their signs only -- each on its own. Two faces TOUCHING ALONG
% THE REGION'S EDGE (the 3-convex-edge triangle's hypotenuse) leave U + sqrt(k) V = g^m M with g the
% edge's row; undivided, no chart box beside the edge had a sign. Exact doubles only; otherwise kept.
    if ~isnumeric(T.UZ) || ~isnumeric(T.VZ), return, end
    U = T.UZ;  V = T.VZ;  ch = false;
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
            for it = 1:4
                if ~mayDivide(U, g) || (any(V(:)) && ~mayDivide(V, g)), break, end
                U2 = ratQ.trimN(ratQ.polyDivExact(U, g, 2));
                if any(V(:)), V2 = ratQ.trimN(ratQ.polyDivExact(V, g, 2)); else, V2 = 0; end
                U = U2;  V = V2;  ch = true;
            end
        catch
        end
    end
    E = T.EZ;  Us = T.UsZ;  Vs = T.VsZ;
    if isnumeric(E), E2 = divideOut(E, factors); ch = ch || ~isequal(E2, E); E = E2; end
    if isnumeric(Us), U3 = divideOut(Us, factors); ch = ch || ~isequal(U3, Us); Us = U3; end
    if isnumeric(Vs) && any(Vs(:)), V3 = divideOut(Vs, factors); ch = ch || ~isequal(V3, Vs); Vs = V3; end
    if ~ch, return, end
    n = max(size(U), size(V));
    T = makeSurd(zPad(U, n), zPad(V, n), T.k, E, Us, Vs);
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
