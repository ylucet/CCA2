function [pt, ok, P] = cellSamplePoint(C, rows, H, R)
% cellSamplePoint  A point strictly inside an arrangement cell whose straight rows' Chebyshev centre
%   fails one of its CURVED rows -- the fallback the meshers' interiorPoint needs.
%
% [input]  C    : the curve table, rows [a b c d e f branch]
%          rows : the cell, [curveIndex side]
%          H    : the cell's straight rows as [a b c] with a x + b y + c >= 0, normalised
%          R    : the half-width of the bounding box the meshers use for unbounded cells
% [output] pt   : 1 x 2, the sampled point with the largest relative slack over EVERY row
%          ok   : false when no sample is strictly inside
%          P    : every sample strictly inside (m x 2), for a caller that checks a property on the cell
%
% WHY IT IS NEEDED (2026-10-16). Both meshers ask a bridge face's membership predicate once per
% cell, at one interior point, and took that point to be the Chebyshev centre of the STRAIGHT rows.
% A cell cut by a conic is often a crescent that misses that centre, and such a cell was dropped as
% undecided -- measured on an indefinite triangle's one-pair envelope: 327 of 332 cells dropped, the
% mesh +infinity on 902 of 2000 interior points where the face list is finite. cellInterior has
% already proved every such cell has an interior, so a point exists; this finds one by a FIXED grid
% (deterministic: the mesh must not depend on a random seed) over a fan of the straight polygon.
    pt = [0 0];  ok = false;  P = zeros(0,2);
    V = polygonVertices(H, R);
    if size(V,1) < 3, return, end
    c0 = mean(V, 1);
    [~, o] = sort(atan2(V(:,2) - c0(2), V(:,1) - c0(1)));
    V = V(o,:);
    n = 12;
    [I, J] = meshgrid(0:n, 0:n);
    keep = I + J < n & I > 0 & J > 0;
    B = [I(keep), J(keep)] / n;                      % strictly inside each fan triangle
    best = 0;
    for e = 1:size(V,1)
        a = V(e,:);  b = V(mod(e, size(V,1)) + 1, :);
        Q = c0 + B(:,1) * (a - c0) + B(:,2) * (b - c0);
        for i = 1:size(Q,1)
            s = slack(C, rows, Q(i,:));
            if s > 0, P(end+1,:) = Q(i,:); end %#ok<AGROW>
            if s > best, best = s;  pt = Q(i,:);  ok = true; end
        end
    end
end

function s = slack(C, rows, p)
% The smallest relative margin by which p satisfies the cell's rows (<= 0: not strictly inside).
    s = inf;
    for t = 1:size(rows,1)
        c = C(rows(t,1), 1:6);  br = C(rows(t,1), 7);  sd = rows(t,2);
        if br ~= 0
            v = ratQ.pairLineEval(c, br, p);
            sc = max(1, norm(p));
        else
            v = c(1)*p(1)^2 + c(2)*p(1)*p(2) + c(3)*p(2)^2 + c(4)*p(1) + c(5)*p(2) + c(6);
            sc = max(abs(c)) * max(1, norm(p))^2;
        end
        s = min(s, sd * v / sc);
        if s <= 0, return, end
    end
end

function V = polygonVertices(H, R)
% Vertices of {H rows >= 0} inside the box [-R, R]^2.
    A = [H; 1 0 R; -1 0 R; 0 1 R; 0 -1 R];
    V = zeros(0,2);
    for i = 1:size(A,1)
        for j = i+1:size(A,1)
            M = A([i j], 1:2);
            if abs(det(M)) < 1e-12, continue, end
            x = (M \ -A([i j], 3)).';
            if all(A(:,1:2) * x.' + A(:,3) >= -1e-9 * max(1, norm(x))), V(end+1,:) = x; end %#ok<AGROW>
        end
    end
    if isempty(V), return, end
    V = uniquetol(V, 1e-12, 'ByRows', true, 'DataScale', 1);
end
