function bridgeRegionCross()
% ITEM 2: DOES THE REGION'S BOUNDARY CONSIST OF EXACTLY THE TWO CONSTRUCTED CURVES?
%
% THE QUESTION, AND WHY IT IS NOT RHETORICAL. `.claude/bridgeRegionBuild.m` constructs two conics
% -- the images of the tie under `grad g1` and `grad g2` -- and both are certainly IN the boundary,
% since they are the two ends of every sweeping segment. But the region is the union of a
% ONE-PARAMETER FAMILY of conics,
%
%       x = M_t s + v_t ,   M_t = (1-t) A1 + t A2 ,   v_t = (1-t) b1 + t b2 ,   t in [0,1] ,
%
% and the union of a family of curves is bounded by its extreme members AND, in general, by the
% family's ENVELOPE. If an envelope arc appears, the region is NOT the intersection of two sign
% conditions and an H-form face cannot hold it -- which would change the type decision, not just
% the arithmetic. So this is tested rather than assumed.
%
% THE TEST. Trace the region's boundary numerically, by exactly the independent classifier that
% `.claude/bridgeRegionProbe.m` uses -- maximise `<s,x> - max(g1,g2)` and ask whether the maximiser
% landed on the tie -- and then evaluate the two CONSTRUCTED conics at every traced boundary point.
% Every point lying on one of them means the boundary is exactly those two curves. A traced point
% on neither is an envelope arc, and it would be visible as a residual many orders larger.
%
% The two probes share no arithmetic: one maximises and bisects in doubles, the other substitutes
% an integer matrix into an integer conic.

    fprintf('\n=== isotropic (the same instance bridgeRegionProbe traced) ===\n');
    a = [1; 0.5];  c = -0.25;
    g1 = struct('A', eye(2),     'b', [0;0], 'k', 0);
    g2 = struct('A', 0.5*eye(2), 'b', a,     'k', -c);
    cross(g1, g2, 2*a);

    fprintf('\n=== anisotropic ===\n');
    g2b = struct('A', [0.5 0.1; 0.1 0.9], 'b', a, 'k', -c);
    cross(g1, g2b, 2*a);
end

% ------------------------------------------------------------------------------------------------
function cross(g1, g2, ctr)
    q = diffConic(g1, g2);
    C = {imageConic(q, g1), imageConic(q, g2)};
    for k = 1:2
        if isempty(C{k})
            fprintf('  boundary %d degenerates (singular A); skipping the cross-check.\n', k);
            return
        end
        fprintf('  constructed boundary %d: [%s]  %s\n', k, num2str(C{k}), circleInfo(C{k}));
    end

    [P, nRay] = traceBoundary(g1, g2, ctr);
    if size(P,1) < 20
        fprintf('  only %d traced boundary points; not enough to conclude.\n', size(P,1));  return
    end
    fprintf('  %d traced boundary points over %d rays\n', size(P,1), nRay);

    v1 = abs(evalConic(C{1}, P));   v2 = abs(evalConic(C{2}, P));
    sc = max(1, max(abs(P(:)))^2 * max(max(abs(C{1})), max(abs(C{2}))));
    d  = min(v1, v2) / sc;                    % distance to the NEARER of the two constructed curves
    fprintf('  relative residual to the nearer constructed conic: max %.3e, median %.3e\n', ...
            max(d), median(d));
    if max(d) < 1e-9
        fprintf(['  => EVERY traced boundary point lies on one of the two constructed conics.\n' ...
                 '     No envelope arc: the region IS bounded by exactly those two curves, so an\n' ...
                 '     H-form face over two sign conditions can hold it.\n']);
    else
        fprintf(['  => SOME traced boundary point lies on NEITHER constructed conic (max %.3e).\n' ...
                 '     That is an envelope arc, and the region is not two sign conditions.\n'], max(d));
    end
end

function s = circleInfo(c)
% If the conic is a circle, say so with its centre and radius -- that is what makes the annulus
% recognisable across the two probes.
    if c(2) ~= 0 || c(1) ~= c(3) || c(1) == 0, s = sprintf('kind %s', ratQ.conicKind(ratQ.conic(c))); return, end
    cx = -c(4)/(2*c(1));  cy = -c(5)/(2*c(1));  r2 = cx^2 + cy^2 - c(6)/c(1);
    if r2 < 0, s = 'an empty circle'; return, end
    s = sprintf('circle centre (%.6f, %.6f) radius %.6f', cx, cy, sqrt(r2));
end

function [P, nRay] = traceBoundary(g1, g2, ctr)
% The independent oracle: classify by WHERE THE MAXIMISER LANDS, then bisect for the switch.
    nRay = 360;
    th = linspace(0, 2*pi, nRay+1);  th(end) = [];
    P = zeros(0,2);
    for i = 1:nRay
        u = [cos(th(i)); sin(th(i))];
        T = linspace(0, 8, 1600);
        isB = false(size(T));
        for j = 1:numel(T), isB(j) = onTie(g1, g2, ctr + T(j)*u); end
        sw = find(diff(isB) ~= 0);
        for j = sw
            t = refine(g1, g2, ctr, u, T(j), T(j+1));
            P(end+1,:) = (ctr + t*u).'; %#ok<AGROW>
        end
    end
end

function t = refine(g1, g2, ctr, u, tA, tB)
    for it = 1:60
        tm = (tA + tB)/2;
        if onTie(g1, g2, ctr + tm*u) == onTie(g1, g2, ctr + tA*u), tA = tm; else, tB = tm; end
    end
    t = (tA + tB)/2;
end

function tf = onTie(g1, g2, x)
    s1 = g1.A \ (x - g1.b);
    if gval(g1,s1) >= gval(g2,s1) - 1e-12*max(1, s1.'*s1), tf = false; return, end
    s2 = g2.A \ (x - g2.b);
    if gval(g2,s2) >= gval(g1,s2) - 1e-12*max(1, s2.'*s2), tf = false; return, end
    tf = true;
end

function v = gval(g, s), v = 0.5*(s.'*g.A*s) + g.b.'*s + g.k; end

% ---- the exact half, copied deliberately from bridgeRegionBuild so the two stay comparable -----
function [cimg] = imageConic(q, gk)
    A = gk.A;  b = gk.b;
    if abs(det(A)) < 1e-12, cimg = []; return, end
    [An, ad] = intify(A);  [bn, bd] = intify(b);
    L = lcm(ad, bd);  An = An * (L/ad);  bn = bn * (L/bd);
    adjA = [An(2,2), -An(1,2); -An(2,1), An(1,1)];
    dAn  = ratQ.detExact(An);
    [n2, ~] = ratQ.substAffine(faceOf(q), 1, ratQ.chk(L*adjA), ratQ.chk(-adjA*bn), dAn);
    cimg = ratQ.conic([n2(5), 2*n2(6), n2(7), 2*n2(8), 2*n2(9), 2*n2(10)]);
end

function n = faceOf(q), n = [0 0 0 0, 2*q(1), q(2), 2*q(3), q(4), q(5), q(6)]; end

function q = diffConic(g1, g2)
    A = g1.A - g2.A;  b = g1.b - g2.b;  k = g1.k - g2.k;
    [q, ~] = intify([A(1,1)/2, A(1,2), A(2,2)/2, b(1), b(2), k]);
end

function [N, d] = intify(X)
    d = 1;
    for v = X(:).'
        [~, dv] = rat(v, 1e-12);
        d = lcm(d, abs(dv));
    end
    N = round(X * d);
end

function v = evalConic(c, X)
    v = c(1)*X(:,1).^2 + c(2)*X(:,1).*X(:,2) + c(3)*X(:,2).^2 + ...
        c(4)*X(:,1)      + c(5)*X(:,2)        + c(6);
end
