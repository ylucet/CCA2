function tf = algFaceBelow(lo, hi, box, inBox, depth, factors, corners)
% algFaceBelow  Is one implicit face BELOW another on a region, CERTIFIED: true only when
%   z_lo(x) <= z_hi(x) was shown at every point of the region in box WHERE BOTH BRANCHES EXIST --
%   strictly, except on a finite union of curves -- and false otherwise.
%
% [input]  lo, hi : implicit faces, .FP (polyN in x, y, z of degree 1 or 2 in z; double or ratQ's big
%                   form) and .FB (algBranch's index). hi may be a STRUCT ARRAY: then the claim is
%                   z_lo <= max_j z_hi(j), certified box by box with whichever hi(j) separates there
%                   (2026-11-01: on erratum Example 2 most faces a pair's cell keeps lie below the
%                   max of the others and below no single one)
%          box    : [x0 x1 y0 y1], a box containing the region
%          inBox  : function handle, inBox(b) false only when sub-box b PROVABLY misses the region
%          depth  : subdivision levels allowed (default 6)
%          factors, corners : optional, passed to algRowFails -- the region's rows as polynomials
%                   >= 0 on it, and its vertices. A face often ENDS on the region's edge (its
%                   discriminant vanishes there), and only with the edge divided out does a box
%                   on that edge get a sign (measured on Example 2: y = 5/2 left every such box open)
% [output] tf     : true when certified. "Except on curves" is algRowFails' own reading (a squared
%                   factor of a row divided out); continuity then gives <= on them, which is all
%                   that dropping the lower face needs. "Where both exist" is algRowFails' reading
%                   (no branch, the row fails): a subdivision's face has its branch on its whole
%                   cell -- its evaluator reads it there -- so on a cell both faces own, it is all.
%
% WHAT IT IS FOR (2026-11-02). envelopeSubdivisionAlg cuts each irrational face (bridge, corner)
% against the quadratics only, so where several own a cell they all own it -- and on erratum
% Example 2 a pair's cell carried 6-11 such faces of which ONE was ever the maximum (300 samples
% per cell, every cell read). The overlay of the pairs then chose one face per pair in every
% combination. A face below another on the whole cell owns a subset of what that one owns (both are
% cut against the same quadratics), so under the max-evaluator dropping it changes no value.
%
% HOW. A QUADRATIC q between the two on a box separates them: z_lo < q and q < z_hi are each an
% algebraic row against q, side -1 and side +1, certified false by algRowFails on that box. q is
% fitted to the two branches at a grid of the box in doubles -- a choice, not a decision: a bad q
% only fails to certify. A box with no q is split. (A constant separator was tried first: where the
% gap is small and the faces slope, the boxes it needs are smaller than the gap over the slope.)
% Two faces TOUCHING inside the region are certified only when the fitted q passes through the
% touching curve with both (each row then carries a squared factor); otherwise they are not, which
% only keeps a face.
    if nargin < 5 || isempty(depth), depth = 6; end
    if nargin < 6, factors = {}; end
    if nargin < 7, corners = zeros(0,2); end
    ctx = struct('inBox', inBox, 'factors', {factors}, 'corners', corners);
    loD = struct('FP', ratQ.bigToDouble(lo.FP), 'FB', lo.FB);
    hiD = arrayfun(@(h) struct('FP', ratQ.bigToDouble(h.FP), 'FB', h.FB), hi, 'UniformOutput', false);  hiD = [hiD{:}];
    tf = rec(lo, hi, loD, hiD, box, ctx, depth);
end

function tf = rec(lo, hi, loD, hiD, b, ctx, depth)
    if ~ctx.inBox(b), tf = true; return, end             % misses the region: nothing to prove here
    [u, v] = meshgrid(linspace(b(1), b(2), 5), linspace(b(3), b(4), 5));
    X = [u(:), v(:)];
    zl = algBranch(loD.FP, X, loD.FB);
    % Only the grid points where the two are ordered: a box straddling the region's edge has points
    % outside it, where the faces may cross or end, and those must not veto the separator. The rows
    % are then certified with a few levels of their own, so the sub-boxes outside the region drop out.
    % With several hi, each is tried on the box, the one above lo at most grid points first.
    n = numel(hi);  ZH = zeros(numel(zl), n);  cnt = zeros(1, n);
    for j = 1:n
        ZH(:,j) = algBranch(hiD(j).FP, X, hiD(j).FB);
        cnt(j) = sum(isfinite(zl) & isfinite(ZH(:,j)) & zl < ZH(:,j));
    end
    [~, order] = sort(cnt, 'descend');
    inner = 3;  if depth == 0, inner = 6; end
    for j = order
        zh = ZH(:,j);
        ok = isfinite(zl) & isfinite(zh) & zl < zh;
        if sum(ok) < 6, break, end
        num = separator(X(ok,:), zl(ok), zh(ok));
        if isempty(num), continue, end
        rl = struct('FP', lo.FP, 'FB', lo.FB, 'num', num(1:10), 'den', num(11), 'side', -1);
        rh = struct('FP', hi(j).FP, 'FB', hi(j).FB, 'num', num(1:10), 'den', num(11), 'side', 1);
        % A LEAF BOX gets 6 levels of row certificate, not 3 (2026-11-01). Near a face's
        % discriminant edge "q < z" for the upper root is "q below the roots' mean OR between the
        % roots", a disjunction a box reads only where one side holds throughout, and the curve
        % where q crosses the mean runs through the box; erratum Example 2's faces kept beside the
        % winner of a pair's cell failed only there. Only leaves pay for it.
        if algRowFails(rl, b, ctx.inBox, inner, ctx.factors, ctx.corners) && ...
           algRowFails(rh, b, ctx.inBox, inner, ctx.factors, ctx.corners)
            tf = true;  return
        end
    end
    if depth == 0, tf = false; return, end
    xm = (b(1) + b(2)) / 2;  ym = (b(3) + b(4)) / 2;
    tf = rec(lo, hi, loD, hiD, [b(1) xm b(3) ym], ctx, depth-1) && ...
         rec(lo, hi, loD, hiD, [xm b(2) b(3) ym], ctx, depth-1) && ...
         rec(lo, hi, loD, hiD, [b(1) xm ym b(4)], ctx, depth-1) && ...
         rec(lo, hi, loD, hiD, [xm b(2) ym b(4)], ctx, depth-1);
end

function num = separator(X, zl, zh)
% A quadratic strictly between zl and zh at the points X, in the weighted basis with a power-of-two
% denominator: [num(1:10), den]. The least-squares quadratic through the midpoints, its coefficients
% rounded to k bits for the fewest k that keeps it strictly between; [] if none up to 2^40.
    num = [];
    x = X(:,1);  y = X(:,2);
    V = [x.^2, x.*y, y.^2, x, y, ones(size(x))];
    ws = warning('off', 'MATLAB:rankDeficientMatrix');  ws2 = warning('off', 'MATLAB:singularMatrix');
    c = V \ ((zl + zh) / 2);
    warning(ws);  warning(ws2);                          % q = c1 x^2 + c2 xy + c3 y^2 + c4 x + c5 y + c6
    for k = 0:40
        r = round(c * 2^k);
        if any(abs(r) >= 2^50), return, end
        q = V * r / 2^k;
        if all(q > zl & q < zh)
            % (num5 x^2 + 2 num6 xy + num7 y^2 + 2 num8 x + 2 num9 y + 2 num10) / (2 den)
            num = [0 0 0 0, 2*r(1), r(2), 2*r(3), r(4), r(5), r(6), 2^k];
            return
        end
    end
end
