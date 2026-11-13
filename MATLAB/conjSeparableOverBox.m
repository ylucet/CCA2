function cj = conjSeparableOverBox(ax, dx, ay, dy, c, lo, hi, dualVars)
% conjSeparableOverBox  The conjugate of a SEPARABLE quadratic over an axis-aligned BOX.
%
% [input]  ax,dx : f's x-part is ax*x^2 + dx*x        (NO cross term -- that is the hypothesis)
%          ay,dy : f's y-part is ay*y^2 + dy*y
%          c     : the constant
%          lo,hi : 1x2 each, the box [lo(1),hi(1)] x [lo(2),hi(2)]
%          dualVars : 1x2 sym, [s_1 s_2]
% [output] cj    : functionNDomain array tiling R^2 -- 4, 6 or 9 cells.
%
% WHY THIS EXISTS. When f separates AND the domain is a product, the conjugate separates too:
%
%       f*(s) = sup_{x in X, y in Y} (s1*x + s2*y) - f1(x) - f2(y) - c
%             = [ sup_x s1*x - f1(x) ] + [ sup_y s2*y - f2(y) ] - c
%             = f1*(s1) + f2*(s2) - c,
%
% because the two suprema share no variable and the box imposes no coupling constraint. So a
% 2-D problem becomes two 1-D problems and a product -- **no 2-D region arithmetic at all**: no
% triangulation, no normal cones, no `solve`, no `isAlways`, no Step 3 cross-piece maximum.
%
% This matters for the QPLIB family the SCIP spike targets, whose objectives are sums of
% x_i*x_j over a box. The DIAGONAL terms (i = j) are separable by construction and land here.
% The off-diagonal ones are not -- the cross term IS the coupling -- and they have their own
% closed form (McCormick / Al-Khayyal-Falk).
% Measured on f = (x^2+y^2)/2 over the unit square: 456 s through triangulate + Step 3, 4.2 s
% through the 2-D KKT route, and this instead. See SCIP_READINESS.md.
%
% NOTE that separability is a property of the FUNCTION AND THE DOMAIN TOGETHER. x*y separates in
% rotated coordinates (u = x+y, v = x-y), but the box is not a box there, so that does not help
% and is deliberately not attempted.
%
% THE 1-D CONJUGATE, which is the whole of the mathematics here. For phi(t) = a*t^2 + d*t on
% [l,h], g(sigma) = sup_{l<=t<=h} sigma*t - phi(t), and psi(t) = sigma*t - phi(t) has
% psi'' = -2a:
%
%   a > 0  psi is CONCAVE, so an interior maximiser is possible: t* = (sigma-d)/(2a), clamped.
%          THREE pieces, breaking at sigma = 2*a*l + d and sigma = 2*a*h + d, with the middle
%          value (sigma-d)^2/(4a).
%   a = 0  psi is affine in t: the max is at an endpoint, switching at sigma = d. TWO pieces.
%   a < 0  psi is CONVEX, so the max is again at an endpoint -- a concave phi never contributes
%          an interior piece. TWO pieces, switching where the two endpoint values cross, at
%          sigma = a*(l+h) + d. (Which agrees with the a = 0 formula, as it must.)
%
% Every breakpoint and value above is a closed-form expression in the box bounds and the
% coefficients; nothing here is solved for.

    s1 = dualVars(1); s2 = dualVars(2);
    px = onedPieces(sym(ax), sym(dx), sym(lo(1)), sym(hi(1)), s1);
    py = onedPieces(sym(ay), sym(dy), sym(lo(2)), sym(hi(2)), s2);

    cj = functionNDomain.empty();
    for i = 1:numel(px)
        for j = 1:numel(py)
            ineqs = [px(i).ineqs, py(j).ineqs];
            r = region(ineqs, dualVars);
            if isempty(r)
                continue
            end
            f = symbolicFunction(simplify(px(i).val + py(j).val - sym(c)));
            cj = [cj, functionNDomain(f, r)]; %#ok<AGROW>
        end
    end
end

% ------------------------------------------------------------------------------------------------
function P = onedPieces(a, d, l, h, s)
% The 1-D conjugate of a*t^2 + d*t over [l,h], as pieces: each carries its VALUE (a sym in s) and
% the INEQUALITIES in s that delimit it (each written as g <= 0, this class's convention).
    phi = @(t) a*t^2 + d*t;
    atLo = s*l - phi(l);
    atHi = s*h - phi(h);

    if isAlways(a > 0, 'Unknown', 'false')
        bLo = 2*a*l + d;        % below this, the unconstrained maximiser is left of l
        bHi = 2*a*h + d;        % above this, it is right of h
        P(1).val = atLo;                 P(1).ineqs = s - bLo;
        P(2).val = (s - d)^2/(4*a);      P(2).ineqs = [bLo - s, s - bHi];
        P(3).val = atHi;                 P(3).ineqs = bHi - s;
    else
        % a <= 0: psi is affine or convex in t, so the maximiser is an endpoint and the two
        % endpoint values cross at a*(l+h)+d.
        bMid = a*(l + h) + d;
        P(1).val = atLo;                 P(1).ineqs = s - bMid;
        P(2).val = atHi;                 P(2).ineqs = bMid - s;
    end
end
