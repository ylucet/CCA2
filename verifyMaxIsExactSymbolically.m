function [ok, report] = verifyMaxIsExactSymbolically(g1, g2, g)
% verifyMaxIsExactSymbolically  Check g == max(g1,g2) over WHOLE REGIONS by closed-form
%   minimisation, never by sampling points. Companion to maxQuaPar's own piece invariants
%   (FARFIELD_FIX_PLAN.md Phase 4): those constrain the pieces before assembly, this one
%   constrains the assembled result, which is what a caller actually evaluates.
%
% [input]  g1, g2 : QuaPar operands, finite everywhere (conjPieceCPLQ output).
%          g      : QuaPar, the claimed max.
% [output] ok     : true iff every check below passes.
%          report : one string per failure (empty when ok).
%
% WHAT IS PROVED. Overlay the result with each operand. For every face k of g, every face i of g1
% and every face j of g2:
%
%   1. IDENTITY. q_k is exactly one of the operands' face quadratics -- g never invents one.
%   2. DOMINATION. On the region R_k n F_i (and R_k n G_j), min(q_k - q_i) >= 0.
%   3. ATTAINMENT. On R_k n F_i n G_j, q_k equals q_i or q_j (whichever it was identified with in
%      1), so with 2 it IS the max there.
%
% Every point of every face is covered by finitely many minimisations, which no ring of samples can
% do. It also settles overlaps for free: if two faces of g both admit a point, 2 gives q_k >= q_l
% and q_l >= q_k there, so they agree and eval cannot pick wrongly.
%
% A NOTE ON A PREMISE THAT LOOKS TRUE AND IS NOT. It is tempting to skip the intersections and check
% "q_k >= every operand face quadratic, globally", on the grounds that each face's quadratic is the
% sup of the same objective over a subset of the primal domain and so a lower bound for the operand.
% It is not: the stored formula for an EDGE face is the max along that edge's whole LINE, not along
% the segment, so outside its own face it can exceed the operand. Measured on every fixture here --
% face 1's quadratic beats face 3's by +inf along one of face 3's rays. Hence the restriction to
% R_k n F_i, which is where the comparison is meaningful.
%
% HOW THE MINIMA ARE COMPUTED (closed form). Along a straight edge or ray, a quadratic and every
% conic constraint are quadratics in the parameter; along a parabolic edge they are quartics in the
% arc's own frame parameter (parabolaArcFrame.conicCoeffs). So the feasible part of an edge is a
% finite union of intervals whose endpoints are polynomial roots, and the minimum on each is an
% endpoint, a stationary point, or -inf on an unbounded one. The minimum over an intersection of
% two closed regions is attained on the boundary of one of them (or at an interior stationary point
% of a strictly convex difference, checked separately), so enumerating both boundaries is complete.
%
% WHAT IS NOT PROVED. (i) That the faces COVER the plane: every test here quantifies over a face, so
% a hole is invisible to all of them; maxQuaPar's partitionReport samples for that, which is
% evidence, not proof. (ii) The arithmetic is floating point -- the STRUCTURE is exhaustive (all
% extrema, no probing) but the comparisons carry a relative tolerance, because conjPieceCPLQ
% produces numeric coefficients to begin with. Both gaps are recorded in FARFIELD_FIX_PLAN.md.
    report = {};
    QF = [g1.f; g2.f];
    owner = [repmat(1, size(g1.f,1), 1); repmat(2, size(g2.f,1), 1)]; %#ok<REPMAT>
    ownIdx = [(1:size(g1.f,1))'; (1:size(g2.f,1))'];

    for k = 1:size(g.f,1)
        which = 0;
        for r = 1:size(QF,1)
            if norm(g.f(k,:) - QF(r,:), Inf) <= 1e-9*(1 + norm(QF(r,:), Inf)), which = r; break, end
        end
        if which == 0
            report{end+1} = sprintf(['face %d carries a quadratic that is neither operand''s ' ...
                '(%s)'], k, mat2str(g.f(k,:), 4)); %#ok<AGROW>
            continue
        end
        for r = 1:size(QF,1)
            d = g.f(k,:) - QF(r,:);
            if norm(d, Inf) == 0, continue, end
            if owner(r) == 1, h2 = g1; else, h2 = g2; end
            [v, where] = minOverIntersection(g, k, h2, ownIdx(r), d);
            sc = 1 + norm(d, Inf) * max(1, max(abs(g.V(:))))^2;
            if v < -1e-7*sc
                report{end+1} = sprintf(['face %d (carrying g%d face %d) is beaten on its own ' ...
                    'region by g%d face %d, by %.4g %s'], k, owner(which), ownIdx(which), ...
                    owner(r), ownIdx(r), -v, where); %#ok<AGROW>
            end
        end
    end
    ok = isempty(report);
end

function [vmin, where] = minOverIntersection(hA, kA, hB, kB, diffRow)
% Minimum of the quadratic diffRow over (face kA of hA) n (face kB of hB), in closed form.
    vmin = inf; where = '';
    [v, w] = minOverBoundaryInside(hA, kA, hB, kB, diffRow);
    if v < vmin, vmin = v; where = w; end
    [v, w] = minOverBoundaryInside(hB, kB, hA, kA, diffRow);
    if v < vmin, vmin = v; where = w; end
    % interior candidate: only a strictly convex difference can take its minimum off the boundary
    H = [diffRow(5), diffRow(6); diffRow(6), diffRow(7)];
    if all(eig(H) > 1e-12*max(1, norm(H)))
        xs = (-H\[diffRow(8); diffRow(9)])';
        if faceAdmits(hA, kA, xs) && faceAdmits(hB, kB, xs)
            v = QuaPar.evalPoly(diffRow, xs);
            if v < vmin, vmin = v; where = sprintf('at the interior point (%.4g,%.4g)', xs(1), xs(2)); end
        end
    end
end

function [vmin, where] = minOverBoundaryInside(hA, kA, hB, kB, diffRow)
% Minimum of diffRow over the part of face kA's BOUNDARY that lies inside face kB.
    vmin = inf; where = '';
    dc = [0.5*diffRow(5), diffRow(6), 0.5*diffRow(7), diffRow(8), diffRow(9), diffRow(10)];
    ECa = hA.edgeConics(); %#ok<NASGU>
    ECb = hB.edgeConics();
    Pb = hB.P{kB};
    Pa = hA.P{kA};
    for t = 1:numel(Pa)
        j = abs(Pa(t));
        a = hA.E(j,1); b = hA.E(j,2);
        A0 = hA.V(a,:); B0 = hA.V(b,:);
        isArc = ~isempty(hA.Ec) && any(hA.Ec(j,:) ~= 0);
        if isArc
            fr = parabolaArcFrame(hA.Ec(j,:), 'verifyMax');
            lo = fr.uOf(A0); hi = fr.uOf(B0);
            if hi < lo, tmp = lo; lo = hi; hi = tmp; end
            cons = cell(1, numel(Pb));
            for m = 1:numel(Pb)
                cons{m} = sign(Pb(m)) * fr.conicCoeffs(ECb(abs(Pb(m)),:));
            end
            fpoly = fr.conicCoeffs(dc);
            iv = feasibleIntervals(cons, lo, hi);
            for z = 1:size(iv,1)
                [v, uarg] = minPolyOnInterval(fpoly, iv(z,1), iv(z,2));
                if v < vmin, vmin = v; where = sprintf('on arc edge %d (u = %.4g)', j, uarg); end
            end
        else
            d = B0 - A0;
            if hA.E(j,3) == 0, tmax = inf; else, tmax = 1; end
            cons = cell(1, numel(Pb));
            for m = 1:numel(Pb)
                cons{m} = sign(Pb(m)) * conicOnLine(ECb(abs(Pb(m)),:), A0, d);
            end
            [qa, qb, qc] = quadOnLine(diffRow, A0, d);
            iv = feasibleIntervals(cons, 0, tmax);
            for z = 1:size(iv,1)
                t0 = iv(z,1); t1 = iv(z,2);
                A2 = qa; B2 = 2*qa*t0 + qb; C2 = qa*t0^2 + qb*t0 + qc;   % shift to tau = t - t0
                [v, targ] = minQuadOn(A2, B2, C2, t1 - t0, ...
                    max(1, norm(diffRow(5:10), Inf))*max(1, norm(A0))^2);
                if v < vmin, vmin = v; where = sprintf('on edge %d (t = %.4g)', j, t0 + targ); end
            end
        end
    end
end

function iv = feasibleIntervals(cons, lo, hi)
% Sub-intervals of [lo,hi] (hi may be inf) on which every polynomial in `cons` is <= 0. Endpoints
% are the polynomials' own real roots, so this is exact in structure: no scanning, no sampling.
    brk = [lo];
    for m = 1:numel(cons)
        c = cons{m};
        c = c(find(abs(c) > 1e-14, 1):end);           % strip leading zeros
        if isempty(c) || numel(c) == 1, continue, end
        r = roots(c);
        r = real(r(abs(imag(r)) < 1e-9*(1+abs(r))));
        brk = [brk, r(r > lo & r < hi)']; %#ok<AGROW>
    end
    if isfinite(hi), brk(end+1) = hi; end
    brk = unique(sort(brk));
    iv = zeros(0,2);
    for z = 1:numel(brk)-1
        if okOn(cons, 0.5*(brk(z) + brk(z+1))), iv(end+1,:) = [brk(z), brk(z+1)]; end %#ok<AGROW>
    end
    if ~isfinite(hi)
        if okOn(cons, brk(end) + 1), iv(end+1,:) = [brk(end), inf]; end
    end
end

function tf = okOn(cons, t)
    tf = true;
    for m = 1:numel(cons)
        if polyval(cons{m}, t) > 1e-9*max(1, abs(t))^2, tf = false; return, end
    end
end

function tf = faceAdmits(h, k, s)
    tf = true;
    EC = h.edgeConics();
    Pe = h.P{k};
    for t = 1:numel(Pe)
        cv = EC(abs(Pe(t)),:); sc = max(1, max(abs(cv)));
        if QuaPar.evalConic(cv, s)*sign(Pe(t))/sc > 1e-9*max(1, norm(s)^2), tf = false; return, end
    end
end

function p = conicOnLine(cv, apex, dir)
% The conic cv restricted to apex + t*dir, as [A B C] (highest power first).
    a = cv(1); b = cv(2); c = cv(3); d = cv(4); e = cv(5);
    A = a*dir(1)^2 + b*dir(1)*dir(2) + c*dir(2)^2;
    B = 2*a*apex(1)*dir(1) + b*(apex(1)*dir(2) + apex(2)*dir(1)) + 2*c*apex(2)*dir(2) ...
        + d*dir(1) + e*dir(2);
    C = QuaPar.evalConic(cv, apex);
    p = [A B C];
end

function [A,B,C] = quadOnLine(diffRow, apex, dir)
% t -> diffRow(apex + t*dir) = A t^2 + B t + C, in matrixForm's convention (entries 5 and 7 of a
% 10-wide row are TWICE the x^2 and y^2 coefficients).
    Q = [diffRow(5), diffRow(6); diffRow(6), diffRow(7)];
    L = [diffRow(8); diffRow(9)]; K = diffRow(10);
    A = 0.5*(dir*Q*dir');
    B = apex*Q*dir' + L'*dir';
    C = 0.5*(apex*Q*apex') + L'*apex' + K;
end

function [vmin, targ] = minQuadOn(A, B, C, tmax, sc)
    if nargin < 5 || ~isfinite(sc) || sc <= 0, sc = max([abs(A), abs(B), abs(C), 1]); end
    tolA = 1e-9*sc;
    vmin = C; targ = 0;
    if isinf(tmax)
        if A < -tolA || (abs(A) <= tolA && B < -tolA), vmin = -inf; targ = inf; return, end
    else
        v1 = A*tmax^2 + B*tmax + C;
        if v1 < vmin, vmin = v1; targ = tmax; end
    end
    if A > tolA
        ts = -B/(2*A);
        if ts > 0 && ts < tmax
            v = A*ts^2 + B*ts + C;
            if v < vmin, vmin = v; targ = ts; end
        end
    end
end

function [mn, targ] = minPolyOnInterval(coeffs, lo, hi)
    cand = [lo, hi];
    d = polyder(coeffs);
    if ~isempty(d) && any(d ~= 0)
        r = roots(d);
        r = real(r(abs(imag(r)) < 1e-9*(1+abs(r))));
        cand = [cand, r(r > lo & r < hi)'];
    end
    vals = polyval(coeffs, cand);
    [mn, i1] = min(vals);
    targ = cand(i1);
end
