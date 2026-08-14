function [ok, report] = verifyFacesCoverThePlane(g)
% verifyFacesCoverThePlane  Prove that the faces of a full-domain QuaPar leave NO HOLE, by
%   deciding whole curve segments in closed form rather than by probing points.
%
% This is the half of FARFIELD_FIX_PLAN.md Phase 4 that verifyMaxIsExactSymbolically cannot
% reach. Every check in that file quantifies over a FACE ("on R_k n F_i, q_k dominates"), so a
% region belonging to no face at all is invisible to all of them; until now the only evidence for
% coverage was maxQuaPar's partitionReport, which SAMPLES. A ring of samples can miss a hole
% forever -- the hole fixed on 2026-08-13 was one point wide at the sampling density used, and was
% found by accident.
%
% [input]  g  : QuaPar, claimed finite everywhere (V/E/Ec/F/P populated).
% [output] ok : true iff every check below passes.
%          report : one string per failure, empty when ok.
%
% ------------------------------------------------------------------------------------------------
% WHAT IS PROVED, AND WHY THESE THREE CHECKS SUFFICE
%
% Write R_k for face k's CONSTRAINT region -- exactly the set QuaPar.eval tests, i.e. every
% bounding conic sign-oriented <= 0, plus chordCuts -- and A = union of the R_k. A is closed and
% nonempty. Let T_k be the region actually bounded by face k's edge cycle. The claim is A = R^2.
%
%   (A) EVERY EDGE SEPARATES TWO FACES. F(j,:) has no zero entry. A zero is a boundary of the
%       domain, and a full-domain result has none.
%   (B) EVERY EDGE LIES IN BOTH ITS FACES. For each face k and each edge j of k, every constraint
%       of face k is <= 0 along the WHOLE of edge j. Decided by maximising each constraint along
%       the edge in closed form: a conic restricted to a segment or a ray is a quadratic in the
%       parameter, and restricted to a parabolic arc is a quartic in that parabola's own frame
%       (parabolaArcFrame.conicCoeffs). No point is probed.
%   (C) NO FACE HAS BOUNDARY ANYWHERE ELSE. For each face k and each constraint c of face k, the
%       set {c = 0} n R_k is contained in the union of face k's own edges that lie on the curve
%       {c = 0}. Decided by splitting the whole curve at the roots of every other constraint --
%       between consecutive roots each constraint has CONSTANT sign, so one evaluation settles a
%       whole interval exactly -- and requiring each feasible interval to sit inside an edge's own
%       parameter range.
%   (D) NO FACE IS COLLAPSED. At most one constraint of a face may vanish identically along one of
%       its edges. Two would mean the face is locally squeezed onto a curve, and (C)'s "finitely
%       many extra-active points" step would not hold.
%
% THE ARGUMENT. Suppose A ~= R^2. Since R^2 is connected and A is closed and nonempty, the
% boundary dA is nonempty. Take z in dA. Then z is in R_k for some k, and z is not interior to
% R_k, so some constraint c of face k is ACTIVE at z; that is, z is in {c = 0} n R_k. By (C), z
% lies on an edge j of face k. By (A), edge j has a second face m, and by (B) edge j lies in R_m
% too. dA is the boundary of a nonempty open set, so it is not a finite point set: it contains a
% non-degenerate sub-arc of some edge j. Pick z in that sub-arc's relative interior and off the
% finitely many points where a second constraint of k or of m is also active -- finitely many by
% (D), since two distinct conics sharing no component meet in at most four points. Near such a z,
% R_k is exactly {c <= 0} and R_m is exactly {c' <= 0}, where c' is the SAME edge conic oriented
% for the other side (that is what F(j,:) listing two faces means). Two opposite closed half-sides
% of one curve cover a full neighbourhood of z, so z is interior to A -- contradicting z in dA.
% Hence A = R^2. QED
%
% WHAT A FAILURE MEANS. (B) failing says an edge is not on its own face's boundary. (C) failing
% says a face's constraint region has boundary that is not one of its edges -- if the extra piece
% runs to infinity, that is precisely the far-field over-extension of SUPPORT_MATRIX.md 4.1, and
% this check names the curve and the parameter range instead of waiting for a probe to land there.
%
% WHAT IS STILL NOT PROVED. The arithmetic is floating point. The STRUCTURE is exhaustive -- every
% root of every restriction, no scanning and no probing -- but the comparisons carry relative
% tolerances, because conjPieceCPLQ produces numeric coefficients in the first place. That caveat
% is identical to verifyMaxIsExactSymbolically's and is recorded in FARFIELD_FIX_PLAN.md.
    report = {};
    if isempty(g.P) || isempty(g.E)
        ok = isempty(report); return
    end
    EC = g.edgeConics();
    ne = size(g.E,1);
    nf = numel(g.P);
    scV = 1 + max(abs(g.V(:)));

    % ---- (A) every edge separates two faces --------------------------------------------------
    for j = 1:ne
        if g.F(j,1) == 0 || g.F(j,2) == 0
            report{end+1} = sprintf(['edge %d borders only one face (F = [%d %d]); a result that ' ...
                'is finite everywhere has no domain boundary, so the other side is a HOLE'], ...
                j, g.F(j,1), g.F(j,2)); %#ok<AGROW>
        end
    end

    for k = 1:nf
        Pe = g.P{k};
        if isempty(Pe)
            report{end+1} = sprintf('face %d has no bounding edges', k); %#ok<AGROW>
            continue
        end
        cons = faceConstraints(g, k, EC);

        % ---- (B) and (D), edge by edge -------------------------------------------------------
        for t = 1:numel(Pe)
            j = abs(Pe(t));
            [ok1, msg] = edgeLiesInFace(g, k, j, cons, EC, scV);
            if ~ok1, report{end+1} = msg; end %#ok<AGROW>
        end

        % ---- (C) constraint by constraint ----------------------------------------------------
        for t = 1:numel(cons)
            msgs = boundaryStaysOnEdges(g, k, t, cons, EC, scV);
            for z = 1:numel(msgs), report{end+1} = msgs{z}; end %#ok<AGROW>
        end
    end
    ok = isempty(report);
end

% =================================================================================================

function cons = faceConstraints(g, k, EC)
% Face k's constraints as sign-oriented conics: x is in R_k iff evalConic(cons(t).ec, x) <= 0 for
% every t. This is EXACTLY what QuaPar.eval tests -- the sign-oriented edge conics of P{k} plus the
% derived chords -- so a violation reported here is a violation of the real point-location rule and
% not of a paraphrase of it.
    Pe = g.P{k};
    cons = struct('ec', {}, 'edge', {}, 'what', {});
    for t = 1:numel(Pe)
        j = abs(Pe(t));
        ec = sign(Pe(t)) * EC(j,:);
        ec = ec / max(1, max(abs(ec)));
        cons(end+1) = struct('ec', ec, 'edge', j, ...
            'what', sprintf('edge %d', j)); %#ok<AGROW>
    end
    cc = g.chordCuts(k, EC);
    for m = 1:size(cc,1)
        % chordCuts returns [n1 n2 c] with "inside" = {n*x' <= c}; as a conic that is
        % [0 0 0 n1 n2 -c] <= 0, which is the same orientation the loop above uses.
        ec = [0 0 0 cc(m,1:2) -cc(m,3)];
        ec = ec / max(1, max(abs(ec)));
        cons(end+1) = struct('ec', ec, 'edge', 0, ...
            'what', sprintf('chord cut %d', m)); %#ok<AGROW>
    end
end

function [ok, msg] = edgeLiesInFace(g, k, j, cons, EC, scV)
% (B) every constraint of face k is <= 0 along the WHOLE of edge j, and (D) at most one of them
% vanishes identically along it.
    ok = true; msg = '';
    [par, rng] = edgeParametrisation(g, j, EC);
    if isempty(par), return, end                 % a degenerate edge is reported elsewhere
    tol = 1e-7 * scV^2;
    zeroEc = zeros(0,6);
    for t = 1:numel(cons)
        p = restrictConic(par, cons(t).ec);
        if isIdenticallyZero(p)
            % This constraint IS the curve being walked. Several may be: a face that a passthrough
            % vertex left with two COLLINEAR edges has one such constraint per edge, and that is
            % harmless as long as they all keep the same side. Two with OPPOSITE orientation is the
            % collapse (D) exists to catch, so record the oriented conic and compare below rather
            % than counting.
            zeroEc(end+1,:) = cons(t).ec / max(abs(cons(t).ec)); %#ok<AGROW>
            continue
        end
        [mx, arg] = maxPolyOn(p, rng(1), rng(2));
        if mx > tol
            ok = false;
            msg = sprintf(['face %d: its own edge %d leaves the face''s constraint region -- %s ' ...
                'reaches %+.4g (tol %.2g) at parameter %.6g along that edge. The edge is then not ' ...
                'on this face''s boundary, and the region on the far side of it belongs to no ' ...
                'face.'], k, j, cons(t).what, mx, tol, arg);
            return
        end
    end
    for a = 2:size(zeroEc,1)
        if norm(zeroEc(a,:) - zeroEc(1,:), Inf) > 1e-7
            ok = false;
            msg = sprintf(['face %d: two of its constraints lie on edge %d''s own curve with ' ...
                'OPPOSITE orientations, so the face is squeezed onto that curve rather than ' ...
                'bounded by it.'], k, j);
            return
        end
    end
end

function msgs = boundaryStaysOnEdges(g, k, t, cons, EC, scV)
% (C) {cons(t) = 0} n R_k is contained in the face's own edges lying on that same curve.
    msgs = {};
    [par, ~] = curveParametrisation(cons(t).ec);
    if isempty(par), return, end                 % not a curve this frame can hold; skip, loudly below
    tol = 1e-7 * scV^2;

    % restrictions of the OTHER constraints to this curve
    rest = {};
    for m = 1:numel(cons)
        if m == t, continue, end
        p = restrictConic(par, cons(m).ec);
        if isIdenticallyZero(p), continue, end   % a duplicate of this same curve constrains nothing
        rest{end+1} = p; %#ok<AGROW>
    end
    iv = feasibleIntervalsOnLine(rest, tol);

    allowed = allowedRanges(g, k, cons(t).ec, par, EC);
    finiteEnds = [iv(isfinite(iv)); allowed(isfinite(allowed))];
    scP = 1; if ~isempty(finiteEnds), scP = max(1, max(abs(finiteEnds))); end
    for z = 1:size(iv,1)
        a = iv(z,1); b = iv(z,2);
        if b - a <= 1e-9*scP, continue, end                      % a single touch point is fine
        if coveredBy(a, b, allowed, 1e-6*scP)
            continue
        end
        if isempty(allowed)
            where = 'the face has NO edge on that curve at all';
        else
            where = sprintf('its own edges cover only %s', rangeStr(allowed));
        end
        msgs{end+1} = sprintf(['face %d: its constraint region has boundary on %s over the ' ...
            'parameter range %s, but %s. The face therefore %s.'], k, cons(t).what, ...
            rangeStr([a b]), where, ...
            ternary(~isfinite(a) || ~isfinite(b), ...
                'runs to infinity along a curve that is not one of its edges -- the far-field over-extension', ...
                'extends past its own edge')); %#ok<AGROW>
    end
end

% ---- parametrisations ---------------------------------------------------------------------------

function [par, rng] = edgeParametrisation(g, j, EC) %#ok<INUSD>
% A parametrisation of EDGE j together with the parameter range the edge itself occupies.
    par = []; rng = [0 0];
    A0 = g.V(g.E(j,1),:); B0 = g.V(g.E(j,2),:);
    isArc = ~isempty(g.Ec) && any(g.Ec(j,:) ~= 0);
    if isArc
        try
            fr = parabolaArcFrame(g.Ec(j,:), 'verifyFacesCover');
        catch
            return
        end
        par = struct('kind', 'arc', 'fr', fr);
        lo = fr.uOf(A0); hi = fr.uOf(B0);
        rng = [min(lo,hi), max(lo,hi)];
    else
        d = B0 - A0;
        if norm(d) < 1e-12, return, end
        par = struct('kind', 'line', 'base', A0, 'dir', d);
        if g.E(j,3) == 0, rng = [0, inf]; else, rng = [0, 1]; end
    end
end

function [par, ok] = curveParametrisation(ec)
% A parametrisation of the WHOLE curve {ec = 0} -- not of an edge, so the range is all of R. A
% straight conic gives a line; a parabolic one gives its own frame, whose u is a global monotone
% parameter for every point of the conic (see parabolaArcFrame).
    par = []; ok = false;
    if max(abs(ec(1:3))) <= 1e-12
        n = ec(4:5);
        if norm(n) < 1e-12, return, end
        base = (-ec(6)/ (n*n')) * n;              % the point of {n*x' + ec(6) = 0} closest to 0
        par = struct('kind', 'line', 'base', base, 'dir', [-n(2), n(1)]/norm(n));
        ok = true; return
    end
    try
        fr = parabolaArcFrame(ec, 'verifyFacesCover');
    catch
        return                                     % not a parabola: no QuaPar edge is one
    end
    par = struct('kind', 'arc', 'fr', fr);
    ok = true;
end

function p = restrictConic(par, ec)
% The conic ec restricted to the parametrisation par, as a polynomial (highest power first).
    if strcmp(par.kind, 'line')
        a = ec(1); b = ec(2); c = ec(3); d = ec(4); e = ec(5);
        ap = par.base; dr = par.dir;
        A = a*dr(1)^2 + b*dr(1)*dr(2) + c*dr(2)^2;
        B = 2*a*ap(1)*dr(1) + b*(ap(1)*dr(2) + ap(2)*dr(1)) + 2*c*ap(2)*dr(2) + d*dr(1) + e*dr(2);
        C = QuaPar.evalConic(ec, ap);
        p = [A B C];
    else
        p = par.fr.conicCoeffs(ec);
    end
end

function r = allowedRanges(g, k, ec, par, EC)
% The parameter ranges, in par's own parameter, of face k's edges that lie on the curve {ec = 0}.
% Merged, so a chord that coincides with a straight edge and a curve carrying several collinear
% edges both come out as one interval.
    r = zeros(0,2);
    Pe = g.P{k};
    for t = 1:numel(Pe)
        j = abs(Pe(t));
        if ~sameCurve(EC(j,:), ec), continue, end
        A0 = g.V(g.E(j,1),:); B0 = g.V(g.E(j,2),:);
        if strcmp(par.kind, 'line')
            % s is the parameter of base + s*dir, so project and divide by |dir|^2 -- dir is a unit
            % vector as curveParametrisation builds it, but do not rely on that here.
            dd = par.dir * par.dir';
            s0 = ((A0 - par.base) * par.dir') / dd;
            s1 = ((B0 - par.base) * par.dir') / dd;
            if g.E(j,3) == 0
                if s1 >= s0, r(end+1,:) = [s0, inf]; else, r(end+1,:) = [-inf, s0]; end %#ok<AGROW>
            else
                r(end+1,:) = [min(s0,s1), max(s0,s1)]; %#ok<AGROW>
            end
        else
            u0 = par.fr.uOf(A0); u1 = par.fr.uOf(B0);
            if g.E(j,3) == 0
                % A curved RAY is not representable (assertCurvedEdgesAreArcs), so an edge on a
                % parabola is always a bounded arc; treat a stray one as its own u-range anyway.
                r(end+1,:) = [min(u0,u1), max(u0,u1)]; %#ok<AGROW>
            else
                r(end+1,:) = [min(u0,u1), max(u0,u1)]; %#ok<AGROW>
            end
        end
    end
    r = mergeRanges(r);
end

function tf = sameCurve(c1, c2)
% Do the two conics define the same curve? Compared after normalising each to unit max-norm, in
% both signs, because a face's own orientation flips the sign of its edge conic.
    if max(abs(c1)) < 1e-14 || max(abs(c2)) < 1e-14, tf = false; return, end
    a = c1/max(abs(c1)); b = c2/max(abs(c2));
    tf = norm(a - b, Inf) < 1e-7 || norm(a + b, Inf) < 1e-7;
end

% ---- exact interval algebra ---------------------------------------------------------------------

function iv = feasibleIntervalsOnLine(rest, tol)
% Sub-intervals of the WHOLE real line on which every polynomial in `rest` is <= 0. Endpoints are
% the polynomials' own real roots; between two consecutive roots of ALL of them each polynomial has
% constant sign, so ONE evaluation per interval decides the whole interval exactly. That is the same
% argument feasibleIntervals in verifyMaxIsExactSymbolically makes, extended to a two-sided
% infinite range.
    brk = [];
    for m = 1:numel(rest)
        c = rest{m};
        i0 = find(abs(c) > 1e-14, 1);
        if isempty(i0), continue, end
        c = c(i0:end);
        if numel(c) < 2, continue, end
        rr = roots(c);
        rr = real(rr(abs(imag(rr)) < 1e-9*(1 + abs(rr))));
        brk = [brk, rr(:)']; %#ok<AGROW>
    end
    brk = uniqueSorted(brk);
    nodes = [-inf, brk, inf];
    iv = zeros(0,2);
    for z = 1:numel(nodes)-1
        a = nodes(z); b = nodes(z+1);
        s = probeInside(a, b);
        if allNonPositive(rest, s, tol)
            iv(end+1,:) = [a b]; %#ok<AGROW>
        end
    end
    iv = mergeRanges(iv);
end

function s = probeInside(a, b)
% One parameter strictly inside (a,b); any point does, because every constraint has constant sign
% there. Both ends may be infinite.
    if isfinite(a) && isfinite(b), s = 0.5*(a+b);
    elseif isfinite(a),            s = a + 1;
    elseif isfinite(b),            s = b - 1;
    else,                          s = 0;
    end
end

function tf = allNonPositive(rest, s, tol)
    tf = true;
    for m = 1:numel(rest)
        if polyval(rest{m}, s) > tol*max(1, abs(s))^4, tf = false; return, end
    end
end

function v = uniqueSorted(x)
    if isempty(x), v = []; return, end
    v = sort(x(:)');
    keep = [true, diff(v) > 1e-12*max(1, max(abs(v)))];
    v = v(keep);
end

function r = mergeRanges(r0)
    r = zeros(0,2);
    if isempty(r0), return, end
    [~, ord] = sort(r0(:,1));
    r0 = r0(ord,:);
    cur = r0(1,:);
    for z = 2:size(r0,1)
        if r0(z,1) <= cur(2) + 1e-9*max(1, abs(cur(2)))
            cur(2) = max(cur(2), r0(z,2));
        else
            r(end+1,:) = cur; %#ok<AGROW>
            cur = r0(z,:);
        end
    end
    r(end+1,:) = cur;
end

function tf = coveredBy(a, b, allowed, tol)
    tf = false;
    for z = 1:size(allowed,1)
        if allowed(z,1) <= a + tol && b <= allowed(z,2) + tol, tf = true; return, end
    end
end

% ---- closed-form extrema ------------------------------------------------------------------------

function [mx, arg] = maxPolyOn(p, lo, hi)
% Maximum of the polynomial p on [lo,hi], either end possibly infinite. Exact: the candidates are
% the two endpoints and the stationary points, plus the leading term's behaviour at an infinite
% end -- as t -> +inf the sign of p is the sign of its leading coefficient, and as t -> -inf it is
% that sign times (-1)^degree.
    p = stripLeadingZeros(p);
    if isempty(p), mx = 0; arg = 0; return, end
    deg = numel(p) - 1;
    if deg == 0, mx = p(1); arg = 0; return, end
    if ~isfinite(hi) && p(1) > 0, mx = inf; arg = inf; return, end
    if ~isfinite(lo) && p(1)*(-1)^deg > 0, mx = inf; arg = -inf; return, end
    cand = [];
    if isfinite(lo), cand(end+1) = lo; end
    if isfinite(hi), cand(end+1) = hi; end
    d = polyder(p);
    if ~isempty(d) && any(abs(d) > 1e-14)
        rr = roots(d);
        rr = real(rr(abs(imag(rr)) < 1e-9*(1 + abs(rr))));
        cand = [cand, rr(rr > lo & rr < hi)'];
    end
    if isempty(cand)
        % Both ends infinite and no stationary point: the polynomial is monotone and, by the two
        % tests above, tends to -inf at both ends, which only a constant can do -- handled already.
        mx = 0; arg = 0; return
    end
    vals = polyval(p, cand);
    [mx, i1] = max(vals);
    arg = cand(i1);
end

function tf = isIdenticallyZero(p)
% Does the polynomial vanish identically? Every coefficient below the noise floor -- this is the
% test that recognises a constraint which IS the curve being walked (its own edge conic), whose
% restriction is 0 everywhere rather than merely <= 0.
    tf = isempty(stripLeadingZeros(p));
end

function q = stripLeadingZeros(p)
    if isempty(p), q = []; return, end
    sc = max(abs(p));
    if sc <= 1e-14, q = []; return, end
    i0 = find(abs(p) > 1e-12*sc, 1);
    q = p(i0:end);
end

% ---- reporting ----------------------------------------------------------------------------------

function s = rangeStr(r)
    s = '';
    for z = 1:size(r,1)
        if z > 1, s = [s ' U ']; end %#ok<AGROW>
        s = [s sprintf('[%s, %s]', numStr(r(z,1)), numStr(r(z,2)))]; %#ok<AGROW>
    end
    if isempty(s), s = '(empty)'; end
end

function s = numStr(v)
    if v == inf, s = '+inf'; elseif v == -inf, s = '-inf'; else, s = sprintf('%.6g', v); end
end

function s = ternary(c, a, b)
    if c, s = a; else, s = b; end
end
