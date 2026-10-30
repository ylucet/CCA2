function tf = algRowFails(a, box, inBox, depth)
% algRowFails  Is an algebraic row side * (q(x) - z(x)) >= 0 FALSE at every point of a region,
%   CERTIFIED by Bernstein bounds on sub-boxes: true only when that was shown, false otherwise.
%   Given SEVERAL rows, is at least one of them false at every point -- their conjunction empty.
%
% [input]  a     : the row(s), as envelopeMeshAlg stores them -- .FP (polyN in x, y, z of degree 1
%                  or 2 in z), .FB (algBranch's index), .num/.den (the quadratic face, weighted
%                  basis), .side (+1 or -1); a struct array for several rows
%          box   : [x0 x1 y0 y1], a box containing the region
%          inBox : function handle, inBox(b) false only when sub-box b PROVABLY misses the region
%          depth : subdivision levels allowed (default 6)
% [output] tf    : true when at every point of the region in box, side * (q - z) < 0 for some row --
%                  except on the zero set of a squared factor (A TOUCHING ROW below), a curve
%
% SEVERAL ROWS (2026-10-30): 274 faces of the 3-convex-edge triangle were empty only JOINTLY -- each
% row holds somewhere on the piece, never all together. A sub-box is decided when ANY row is
% certified false on it, so different rows may decide different boxes.
%
% WHY POINTWISE (2026-10-29). The first test (envelopeMeshAlg, 2026-10-20) proved the
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
    if nargin < 4, depth = 6; end
    tf = false;
    rows = cell(1, numel(a));
    for t = 1:numel(a)
        rows{t} = rowPolys(a(t));
    end
    rows = rows(~cellfun(@isempty, rows));
    if isempty(rows), return, end
    tf = boxRec(rows, box, inBox, depth);
end

function r = rowPolys(a)
% The polynomials whose signs decide row a on a box, and the test reading them; [] if the row
% cannot be certified at all (shape, branch index, non-finite coefficients).
    r = [];
    FP = a.FP;
    if ratQ.isBigPoly(FP), FP = ratQ.fromBig(FP); end
    num = double(a.num(:).');
    if numel(num) ~= 10 || any(num(1:4) ~= 0) || ~(a.den > 0), return, end
    % q = Q / Dn: twice the weighted numerator, in the plain monomial basis (as bridgeQuadCut)
    Q = zeros(3,3);
    Q(3,1) = num(5);    Q(2,2) = 2*num(6);  Q(1,3) = num(7);
    Q(2,1) = 2*num(8);  Q(1,2) = 2*num(9);  Q(1,1) = 2*num(10);
    Dn = 2 * double(a.den);
    want = -a.side;                                  % the sign q - z must have for the row to fail
    switch size(FP, 3) - 1
        case 2
            A = FP(:,:,3);  B = FP(:,:,2);  C0 = FP(:,:,1);
            hi = any(a.FB == [-1 2]);
            if ~hi && ~any(a.FB == [1 -2]), return, end
            P = {A, ...                                                   % a
                 padAdd(2 * conv2(A, Q), Dn * B), ...                     % Dn (2aq + b)
                 squareFreeSign(padAdd(padAdd(conv2(A, conv2(Q, Q)), Dn * conv2(B, Q)), Dn^2 * C0)), ...  % Dn^2 F(q)
                 padAdd(conv2(B, B), -4 * conv2(A, C0))};                 % D
            test = @(s) quadOK(s, hi, want);
        case 1
            A = FP(:,:,2);  C0 = FP(:,:,1);
            if ~any(a.FB == [1 -1]), return, end
            P = {A, squareFreeSign(padAdd(conv2(A, Q), Dn * C0))};       % b, Dn (b q + c)
            test = @(s) s(1) ~= 0 && s(1) * s(2) == want;
        otherwise
            return
    end
    if any(cellfun(@(p) ~all(isfinite(p(:))), P)), return, end
    r = struct('P', {P}, 'test', test);
end

function tf = boxRec(rows, b, inBox, depth)
    if ~inBox(b), tf = true; return, end             % misses the region: nothing to prove here
    for t = 1:numel(rows)
        P = rows{t}.P;  s = zeros(1, numel(P));
        for k = 1:numel(P)
            v = polySignOn(P{k}, b, @(~) true, 0);
            if ~isnan(v), s(k) = v; end
        end
        if rows{t}.test(s), tf = true; return, end
    end
    if depth == 0, tf = false; return, end
    xm = (b(1) + b(2)) / 2;  ym = (b(3) + b(4)) / 2;
    tf = boxRec(rows, [b(1) xm b(3) ym], inBox, depth-1) && ...
         boxRec(rows, [xm b(2) b(3) ym], inBox, depth-1) && ...
         boxRec(rows, [b(1) xm ym b(4)], inBox, depth-1) && ...
         boxRec(rows, [xm b(2) ym b(4)], inBox, depth-1);
end

function ok = quadOK(s, hi, want)
% s = signs of [a, Dn(2aq + b), Dn^2 F(q), D] on the box, 0 where not certified.
    if s(1) == 0, ok = false; return, end
    u = s(2) * s(1);                                 % sign(q - m)
    v = s(3) * s(1);                                 % sign(F(q) / a): -1 is "between the roots"
    dPos = s(4) == 1;
    if hi && want == -1,      ok = (u == -1 && dPos) || v == -1;
    elseif hi,                ok = u == 1 && v == 1 && dPos;
    elseif want == -1,        ok = u == -1 && v == 1 && dPos;
    else,                     ok = (u == 1 && dPos) || v == -1;
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

function D = diffPoly(S, d)
% dS/dx (d = 1) or dS/dy (d = 2) of a dense polyN, S(i+1,j+1) the coefficient of x^i y^j.
    if d == 2, D = diffPoly(S.', 1).'; return, end
    if size(S,1) < 2, D = zeros(1, size(S,2)); return, end
    D = S(2:end,:) .* (1:size(S,1)-1).';
end

function S = padAdd(A, B)
    n = max(size(A), size(B));
    S = zeros(n);  S(1:size(A,1), 1:size(A,2)) = A;
    S(1:size(B,1), 1:size(B,2)) = S(1:size(B,1), 1:size(B,2)) + B;
end
