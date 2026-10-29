function tf = algRowFails(a, box, inBox, depth)
% algRowFails  Is an algebraic row side * (q(x) - z(x)) >= 0 FALSE at every point of a region,
%   CERTIFIED by Bernstein bounds on sub-boxes: true only when that was shown, false otherwise.
%
% [input]  a     : the row, as envelopeMeshAlg stores it -- .FP (polyN in x, y, z of degree 1 or 2
%                  in z), .FB (algBranch's index), .num/.den (the quadratic face, weighted basis),
%                  .side (+1 or -1)
%          box   : [x0 x1 y0 y1], a box containing the region
%          inBox : function handle, inBox(b) false only when sub-box b PROVABLY misses the region
%          depth : subdivision levels allowed (default 6)
% [output] tf    : true when side * (q - z) < 0 holds on every point of the region in box
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
    if nargin < 4, depth = 6; end
    tf = false;
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
                 padAdd(padAdd(conv2(A, conv2(Q, Q)), Dn * conv2(B, Q)), Dn^2 * C0), ...  % Dn^2 F(q)
                 padAdd(conv2(B, B), -4 * conv2(A, C0))};                 % D
            test = @(s) quadOK(s, hi, want);
        case 1
            A = FP(:,:,2);  C0 = FP(:,:,1);
            if ~any(a.FB == [1 -1]), return, end
            P = {A, padAdd(conv2(A, Q), Dn * C0)};                       % b, Dn (b q + c)
            test = @(s) s(1) ~= 0 && s(1) * s(2) == want;
        otherwise
            return
    end
    if any(cellfun(@(p) ~all(isfinite(p(:))), P)), return, end
    tf = boxRec(P, test, box, inBox, depth);
end

function tf = boxRec(P, test, b, inBox, depth)
    if ~inBox(b), tf = true; return, end             % misses the region: nothing to prove here
    s = zeros(1, numel(P));
    for k = 1:numel(P)
        v = polySignOn(P{k}, b, @(~) true, 0);
        if ~isnan(v), s(k) = v; end
    end
    if test(s), tf = true; return, end
    if depth == 0, tf = false; return, end
    xm = (b(1) + b(2)) / 2;  ym = (b(3) + b(4)) / 2;
    tf = boxRec(P, test, [b(1) xm b(3) ym], inBox, depth-1) && ...
         boxRec(P, test, [xm b(2) b(3) ym], inBox, depth-1) && ...
         boxRec(P, test, [b(1) xm ym b(4)], inBox, depth-1) && ...
         boxRec(P, test, [xm b(2) ym b(4)], inBox, depth-1);
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

function S = padAdd(A, B)
    n = max(size(A), size(B));
    S = zeros(n);  S(1:size(A,1), 1:size(A,2)) = A;
    S(1:size(B,1), 1:size(B,2)) = S(1:size(B,1), 1:size(B,2)) + B;
end
