function s = polySignOn(P, box, inBox, depth)
% polySignOn  The sign of a bivariate polynomial on a region, CERTIFIED by Bernstein bounds: +1 or -1
%   when P keeps that strict sign on every point of the region, 0 when that could not be shown.
%
% [input]  P     : dense polyN in (x, y) -- P(i+1,j+1) multiplies x^i y^j (bridgeQuadCut's R)
%          box   : [x0 x1 y0 y1], a box containing the region
%          inBox : function handle, inBox(b) false only when sub-box b PROVABLY misses the region
%                  (b is [x0 x1 y0 y1]); sub-boxes it rejects are not examined
%          depth : subdivision levels allowed (default 6)
% [output] s     : +1 / -1 / 0
%
% WHY BERNSTEIN (2026-10-20). A polynomial on a box is a convex combination of its Bernstein
% coefficients, so all of them positive proves it positive there -- one linear map, no root finding.
% Where they mix, the box is split in four and each quarter asked again, and only quarters that may
% meet the region are kept. The coefficients are computed in doubles; a coefficient counts only when
% it clears a margin relative to the largest one, so rounding cannot turn a zero into a sign.
%
% WHAT IT IS FOR: proving a subdivision piece EMPTY when its emptiness is an algebraic row (a quadratic
% against an irrational face) -- see algRowFails and envelopeSubdivisionAlg's algEmpty. Answering 0 only costs a face.
    if nargin < 4, depth = 6; end
    s = signRec(P, box, inBox, depth);
end

function s = signRec(P, b, inBox, depth)
    if ~inBox(b), s = NaN; return, end              % misses the region: no constraint from here
    B = bernstein(P, b);
    m = max(abs(B(:)));
    tol = 1e-10 * max(m, eps);
    if all(B(:) > tol), s = 1; return, end
    if all(B(:) < -tol), s = -1; return, end
    if depth == 0, s = 0; return, end
    xm = (b(1) + b(2)) / 2;  ym = (b(3) + b(4)) / 2;
    q = [signRec(P, [b(1) xm b(3) ym], inBox, depth-1), signRec(P, [xm b(2) b(3) ym], inBox, depth-1), ...
         signRec(P, [b(1) xm ym b(4)], inBox, depth-1), signRec(P, [xm b(2) ym b(4)], inBox, depth-1)];
    q = q(~isnan(q));
    if isempty(q), s = NaN; return, end
    if all(q == 1), s = 1; elseif all(q == -1), s = -1; else, s = 0; end
end

function B = bernstein(P, b)
% Bernstein coefficients of P on box b: substitute x = x0 + hx u, y = y0 + hy v, then convert the
% power basis in (u, v) on [0,1]^2. All as matrix products with cached binomial tables -- the loop
% form called nchoosek ~2e9 times on one subdivision (2026-10-20).
    [n1, n2] = size(P);
    A = shiftMat(n1, b(1), b(2) - b(1)) * P * shiftMat(n2, b(3), b(4) - b(3)).';
    B = convMat(n1) * A * convMat(n2).';
end

function S = shiftMat(n1, x0, h)
% S(k+1, i+1) = C(i,k) x0^(i-k) h^k : coefficients of u^k in (x0 + h u)^i.
    Cb = binom(n1);
    i = 0:n1-1;  k = (0:n1-1).';
    E = i - k;                                    % exponent of x0, negative where k > i
    X = x0 .^ max(E, 0);  X(E < 0) = 0;
    S = Cb.' .* X .* (h .^ k);
end

function M = convMat(n1)
% M(i+1, k+1) = C(i,k) / C(n,k) for k <= i: power basis on [0,1] to Bernstein, degree n = n1 - 1.
    persistent cache
    if isempty(cache), cache = containers.Map('KeyType', 'double', 'ValueType', 'any'); end
    if isKey(cache, n1), M = cache(n1); return, end
    Cb = binom(n1);  n = n1 - 1;
    M = Cb ./ Cb(n+1, :);                         % C(i,k) / C(n,k)
    M = tril(M);
    cache(n1) = M;
end

function Cb = binom(n1)
% Cb(i+1, k+1) = C(i, k), 0 <= i, k < n1 (zero above the diagonal).
    persistent cache
    if isempty(cache), cache = containers.Map('KeyType', 'double', 'ValueType', 'any'); end
    if isKey(cache, n1), Cb = cache(n1); return, end
    Cb = zeros(n1);
    for i = 0:n1-1
        for k = 0:i, Cb(i+1, k+1) = nchoosek(i, k); end
    end
    cache(n1) = Cb;
end
