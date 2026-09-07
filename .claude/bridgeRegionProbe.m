function bridgeRegionProbe()
% IS A BRIDGE FACE'S PRIMAL DOMAIN POLYHEDRAL? The measurement the type decision waits on.
%
% WHY IT DECIDES THE TYPE. CCA2's piecewise types are a 2 x 2 lattice: a FUNCTION axis (Rat ⊃ Qua)
% crossed with a SUBDIVISION axis (Par ⊃ Pol, curved cell boundaries or straight ones). The open
% question is whether the bridge face needs a THIRD value on the function axis -- an implicit
% p(x,y,z) = 0 plus a branch selector. That estimate is only valid if the face's DOMAIN stays
% polyhedral; if the domain is curved too, the change touches both axes and the option being
% weighed is not the real option.
%
% WHAT THE BRIDGE REGION IS. f** (x) = sup_s <s,x> - f*(s) with f* = max(g1, g2). The sup is
% attained at a unique s*(x), and x belongs to the BRIDGE exactly when s*(x) lies on the TIE
% {g1 = g2} rather than strictly inside one piece's cell -- there the envelope is supported by both
% pieces at once, which is what a bridge face is.
%
% THE CLOSED FORM THIS TESTS AGAINST. At a tie point the subdifferential of a max is the segment
% between the two gradients, so
%
%       bridge region  =  union over s on the tie of  conv{ grad g1(s), grad g2(s) } .
%
% Each g_k is quadratic, so each grad g_k is AFFINE, and the tie is a CONIC. An affine image of a
% conic is a conic -- so the region's two boundary curves are conics, and they are straight only in
% the degenerate case where the tie itself is a line (which DECISIONS 2026-09-05 pins as exactly the
% equal-Hessian case). The prediction is therefore that the domain is NOT polyhedral, and for the
% isotropic instance it is sharper still: with g1 = |s|^2/2 and g2 = <s,a> + |s|^2/4 - c the tie is
% the circle of centre 2a and radius R = 2*sqrt(|a|^2 - c), grad g1(s) = s and grad g2(s) = a + s/2,
% so the two boundaries are CONCENTRIC CIRCLES about 2a of radii R and R/2 and the region is an
% ANNULUS.
%
% This measures the region without using any of that: it maximises numerically and classifies by
% where the maximiser landed. Agreement with the annulus is then evidence, not a restatement.
%
% NOTE ON LOCATELLI, because the opposite conclusion is recorded elsewhere in this project.
% DESIGN.md treats an empirically curved transition boundary as proof of an INCOMPLETE SEARCH, on
% the strength of Locatelli's theorem that the convex envelope of a quadratic over a polyhedron
% admits a polyhedral subdivision. That theorem is about ONE quadratic over ONE polyhedron. A bridge
% spans TWO quadratics across a kink, so the theorem does not apply here and a curved boundary is
% the expected answer rather than an artifact. The two findings do not conflict.

    fprintf('\n=== ISOTROPIC: g1 = |s|^2/2,  g2 = <s,a> + |s|^2/4 - c ===\n');
    a = [1; 0.5];  c = -0.25;
    A1 = eye(2);        b1 = [0;0];      k1 = 0;
    A2 = 0.5*eye(2);    b2 = a;          k2 = -c;
    R  = 2*sqrt(a.'*a - c);
    fprintf('predicted annulus: centre (%.4f, %.4f), radii %.6f and %.6f\n', 2*a(1), 2*a(2), R/2, R);
    report(A1,b1,k1, A2,b2,k2, 2*a);

    fprintf('\n=== ANISOTROPIC: the same, with g2 anisotropic (no closed form to lean on) ===\n');
    A2b = [0.5 0.1; 0.1 0.9];
    report(A1,b1,k1, A2b,b2,k2, 2*a);
end

% ------------------------------------------------------------------------------------------------
function report(A1,b1,k1, A2,b2,k2, ctr)
% Trace the bridge region's two boundaries along rays from ctr, then ask whether each traced
% boundary is a LINE or a CONIC. Both fits are least squares on the same points, so the comparison
% is between two models on one data set rather than between two data sets.
    nRay = 240;
    th = linspace(0, 2*pi, nRay+1);  th(end) = [];
    inner = nan(nRay,1);  outer = nan(nRay,1);
    for i = 1:nRay
        u = [cos(th(i)); sin(th(i))];
        [lo, hi] = bracketBridge(A1,b1,k1, A2,b2,k2, ctr, u);
        inner(i) = lo;  outer(i) = hi;
    end
    ok = isfinite(inner) & isfinite(outer);
    fprintf('rays with a bridge interval: %d of %d\n', sum(ok), nRay);
    if sum(ok) < 20, fprintf('  too few to fit; the region is not being found.\n'); return, end

    for which = 1:2
        if which == 1, r = inner; nm = 'inner'; else, r = outer; nm = 'outer'; end
        P = [ctr(1) + r(ok).*cos(th(ok).'), ctr(2) + r(ok).*sin(th(ok).')];
        [eL, eC] = fitBoth(P);
        span = max(range(P(:,1)), range(P(:,2)));
        fprintf(['  %s boundary: %d pts, extent %.4f | best LINE residual %.3e | ' ...
                 'best CONIC residual %.3e\n'], nm, size(P,1), span, eL, eC);
        fprintf('      radius from centre: min %.6f max %.6f  (spread %.2e)\n', ...
                min(r(ok)), max(r(ok)), max(r(ok)) - min(r(ok)));
    end
end

function [lo, hi] = bracketBridge(A1,b1,k1, A2,b2,k2, ctr, u)
% Along the ray ctr + t*u, find the t-interval on which the maximiser lies on the TIE. Bisection on
% the classifier, which is monotone in t here; a ray that never enters the bridge returns NaN.
    T = linspace(0, 6, 1200);
    isB = false(size(T));
    for i = 1:numel(T)
        isB(i) = onTie(A1,b1,k1, A2,b2,k2, ctr + T(i)*u);
    end
    j = find(isB, 1, 'first');  k = find(isB, 1, 'last');
    if isempty(j) || k == numel(T), lo = NaN; hi = NaN; return, end
    lo = refine(A1,b1,k1, A2,b2,k2, ctr, u, T(max(j-1,1)), T(j));
    hi = refine(A1,b1,k1, A2,b2,k2, ctr, u, T(k), T(k+1));
end

function t = refine(A1,b1,k1, A2,b2,k2, ctr, u, tA, tB)
    for it = 1:60
        tm = (tA + tB)/2;
        if onTie(A1,b1,k1, A2,b2,k2, ctr + tm*u) == onTie(A1,b1,k1, A2,b2,k2, ctr + tA*u)
            tA = tm;
        else
            tB = tm;
        end
    end
    t = (tA + tB)/2;
end

function tf = onTie(A1,b1,k1, A2,b2,k2, x)
% Is the maximiser of <s,x> - max(g1(s), g2(s)) a point where the two pieces TIE?
%
% The objective is concave, so the smooth candidate of each branch is available in closed form
% (grad g_k(s) = x, i.e. s = A_k \ (x - b_k)) and is the true maximiser exactly when the OTHER
% piece is not larger there. If neither smooth candidate is feasible, the maximiser is on the ridge
% -- which is the bridge. Deciding it this way rather than by running an optimiser keeps the
% classifier free of a solver's stopping tolerance.
    s1 = A1 \ (x - b1);
    if g(A1,b1,k1,s1) >= g(A2,b2,k2,s1) - 1e-12 * sc(s1), tf = false; return, end
    s2 = A2 \ (x - b2);
    if g(A2,b2,k2,s2) >= g(A1,b1,k1,s2) - 1e-12 * sc(s2), tf = false; return, end
    tf = true;
end

function v = g(A,b,k,s), v = 0.5*(s.'*A*s) + b.'*s + k; end
function v = sc(s),      v = max(1, s.'*s);              end

function [eL, eC] = fitBoth(P)
% Least-squares residuals of the best LINE and the best general CONIC through P, both normalised by
% the extent of P so the two numbers are comparable and scale-free. Each is the smallest singular
% value of the corresponding design matrix, i.e. the algebraic residual of the best unit-norm fit.
    P = P - mean(P,1);
    s = max(vecnorm(P, 2, 2));
    P = P / s;
    ML = [P, ones(size(P,1),1)];
    eL = min(svd(ML)) / sqrt(size(P,1));
    MC = [P(:,1).^2, P(:,1).*P(:,2), P(:,2).^2, P(:,1), P(:,2), ones(size(P,1),1)];
    eC = min(svd(MC)) / sqrt(size(P,1));
end
