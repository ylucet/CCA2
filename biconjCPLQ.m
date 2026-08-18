function h = biconjCPLQ(obj)
% biconjCPLQ  Biconjugate f** via the 'cplq' engine.
%
% objective: compute f** = (f*)*, the closed convex envelope of a 2D piecewise linear-quadratic
%   function, dispatching by input shape exactly as conjCPLQ.m does (its Cases A/B/C).
%
% [input]  obj : QuaPol, or a polyhedral QuaPar (all-zero Ec). Operable (degree<=2).
% [output] h   : the biconjugate, a RatPar. Concretely
%                  Case A (full-domain strictly convex quadratic) -> QuaPol   (= f itself)
%                  Case B (single bounded triangle)               -> RatPol   (see below)
%                  Case C (general bounded domain)                -> QuaParCPLQ
%                Call kind() on the result to learn which; see RETURN_TYPE.md.
%
% WHY CASE B IS NOT conj-OF-conj (this function's reason to exist)
% ---------------------------------------------------------------
% For a bounded-domain input, f* is finite everywhere, i.e. it lives on an UNBOUNDED multi-face
% parabolic subdivision -- a shape conjCPLQ still rejects (conjCPLQ.m:103, SUPPORT_MATRIX.md
% section 1.2). So the literal double conjugation `conj(conj(f))` cannot close the loop on the
% numeric path, and `biconj` used to fail for every single bounded triangle -- half of this
% project's stated goal (SUPPORT_MATRIX.md section 8, blocker 1).
%
% It does not have to go that way. For f = q + I_T with T a compact convex set and q continuous
% on T, f is proper and bounded below, so
%
%       f** = cl conv f = conv f            (conv f is already closed: the convex envelope of a
%                                            continuous function on a compact convex polytope is
%                                            continuous on that polytope)
%
% and conv(q + I_T) is EXACTLY what Step 1 of the [COAP]/[JOGO] algorithm computes -- convEnvCPLQ.
% So for Case B the biconjugate is Step 1's own output, in closed form, with no second conjugation
% and no unbounded domain anywhere. This is the "reformulate rather than extend" option
% SUPPORT_MATRIX.md section 3 lists; it leaves the general unbounded multi-face conjugate
% (section 1.2) open, which is a separate, larger task and still the right long-term fix for
% `conj` itself.
%
% This route is also STRICTLY STRONGER than conj-of-conj would have been on Case B, not merely
% equal: it succeeds on triangles whose CONJUGATE conjCPLQ still cannot compute at all (an
% indefinite piece whose Step 1 envelope contains a rational face -- Step 2 rejects those,
% conjPieceCPLQ.m:107, reached via conjCPLQ's own conjMaxOfSubTriangles). Step 1 produces the
% envelope regardless, so f** is available even where f* is not.
%
% VALIDATION (2026-07-28): checked against an INDEPENDENT ground truth that uses none of the
% conjugate pipeline -- f*(s) computed by exact maximization of <s,x> - q(x) over T (vertices,
% per-edge 1D critical points, interior critical point), then f**(x) = sup_s <s,x> - f*(s) by
% grid + Nelder-Mead. Agreement was <= 8.9e-16 at 6 interior points on each of 7 triangles
% covering every Step 1 branch: affine, convex, concave, indefinite 0/1-convex-edge (1-face
% envelope), indefinite 2-convex-edge (2-face split envelope) and indefinite 3-convex-edge
% (4-face split envelope). The last two are exactly the inputs whose conjugate errors. See
% biconjCPLQTest.m, which re-runs the same check.
%
% Cases A and C keep the literal conj-of-conj: A's conjugate is again a full-domain quadratic
% (Case A twice), and C's is a QuaParCPLQ, whose conj drives cPLQ's own symbolic machinery
% (QuaParCPLQ.m) and handles the unbounded parabolic domain. Both already worked.

    obj.assertOperable();   % degree<=2 (cubic rejected; cubic is for isConvex only)

    % ---- f IS ALREADY CONVEX -> f** = f, and there is nothing whatever to compute -----------
    % biconj IS the closed-convex-envelope operator, so a convex f is its own answer. Taken
    % first, before any dispatch, because it is the largest short-circuit here and it costs one
    % eigenvalue test per face.
    %
    % Two ways to reach it, and they differ in what they trust:
    %   * PROVEN, for a single piece -- per-face convexity plus a convex domain IS convexity when
    %     there is only one face, so nothing is assumed;
    %   * ASSERTED, via the caller's fIsConvex flag, for several faces, where the honest test
    %     needs the gradient jump across every shared edge to be consistent. The free NECESSARY
    %     condition is still enforced there (see convexEnough), so a wrong flag is refused loudly
    %     rather than silently answered.
    if convexEnough(obj)
        h = obj;
        return
    end

    % ---- f SEPARABLE over a BOX -> one 1-D envelope per axis, no conjugation at all ---------
    % The convex envelope passes through a separable sum on a product domain:
    %       f = f1(x) + f2(y) + c  on  X x Y   =>   co f = co(f1 + I_X) + co(f2 + I_Y) + c,
    % which follows from conjugating twice -- (f1 (+) f2)* = f1* (+) f2*, applied to f and then
    % to f* -- and needs the domain to be a PRODUCT so the two suprema share no constraint.
    %
    % Each 1-D envelope is immediate: a quadratic with a >= 0 is already convex on the interval
    % and is its own envelope; with a < 0 it is concave, so its envelope is the CHORD through the
    % two endpoint values. No conjugation, no cells, no Step 3.
    %
    % WHY THIS IS SEPARATE FROM THE CONVEXITY SHORT-CIRCUIT ABOVE: it catches the INDEFINITE
    % diagonal, which is not convex and so never reaches that branch. f = x^2 - y^2 over the unit
    % box has co f = x^2 - y, and cost 29 s to reach through conj-of-conj -- the first conjugation
    % took the separable route in about 2 s and the second, whose domain is the whole plane and
    % therefore not a box, did not.
    [okSep, envCoef] = separableEnvelopeCoefs(obj);
    if okSep
        h = obj;
        h.f = envCoef;
        return
    end

    % ---- Case B: a single bounded TRIANGLE -> Step 1's convex envelope ----------------------
    % THE TRIANGLE RESTRICTION IS REAL, and widening it to "any bounded piece" was tried and
    % REVERTED on 2026-08-18. convEnvCPLQ's classification is by Hessian: PSD is convex and
    % returns q on any convex P, but the NEGATIVE SEMIDEFINITE and INDEFINITE branches are stated
    % for a TRIANGLE and raise convEnvCPLQ:notImplemented otherwise. Widening therefore turned
    % all four bilinear box rows of checkBoxEnvelopeForSCIP from OK into that error.
    %
    % Nothing was lost by reverting: the case the widening was meant to reach -- a CONVEX
    % quadratic over a box -- is caught by the convexity short-circuit above, which returns f
    % itself in 0.05 s without calling Step 1 at all.
    %
    % Same predicate as conjCPLQ.m's own Case B, plus "no curved edge": a QuaPar whose triangle
    % has a parabolic side is NOT what convEnvCPLQ takes (it reads V/E/F as straight edges), so
    % it falls through to the conj-of-conj path below and errors there exactly as it did before.
    if obj.nf == 1 && obj.nv == 3 && obj.ne == 3 && obj.isDomBounded && ~hasCurvedEdge(obj)
        h = convEnvCPLQ(asQuaPol(obj));
        return
    end

    % ---- A BOUNDED multi-face domain: first conjugate in SYMBOLIC form, deliberately ---------
    % f** = (f*)*, and f* is the same function whichever representation conjCPLQ returns -- so
    % this is a representation choice, not a mathematical one. It has to be made because the two
    % representations are not equally conjugable: since 2026-08-13 the numeric path completes here
    % and returns a MESHED QuaPar with parabolic edges (which is what the SCIP bridge wants from
    % conj), and handing that to a second conjugation dies at quaPolToPlq:curvedEdge, because
    % quaPolToPlq requires a polyhedral domain. The symbolic QuaParCPLQ conjugates again through
    % cPLQ's own machinery, which handles a curved domain.
    %
    % So conj keeps returning the mesh and biconj asks for the symbolic form. Measured on f = x*y
    % over the unit square given as TWO triangles: this route returns the McCormick envelope
    % max(0, x+y-1); the meshed route errors before producing anything.
    if obj.isDomBounded && ~hasCurvedEdge(obj)
        g = conj(obj, 'cplq');
        % ONLY when the mesh is CURVED. Taking the symbolic route unconditionally here was
        % measured to be a regression: for a CONVEX q on this domain the numeric first conjugate
        % is polyhedral, the second conjugation handles it, and forcing the symbolic form instead
        % puts it back on the functionNDomain.getInterior chain that SUPPORT_MATRIX.md section 7
        % records as broken -- conjCPLQTest/biconjCoverageByInputCase pins exactly that, and it
        % went red. A curved mesh is the only thing the second conjugation cannot take.
        if isMeshed(g) && hasCurvedEdge(g)
            g = conjCPLQ(asQuaPol(obj), [], 'symbolic');
        end
        h = conj(g, 'cplq');
        return
    end

    % ---- Cases A and C (and every unsupported shape): unchanged double conjugation ----------
    h = conj(conj(obj, 'cplq'), 'cplq');
end

% ================================================================================================
function tf = hasCurvedEdge(obj)
% True iff any edge carries a nonzero conic. Ec is all-zero on a QuaPol (RatPar pins it there),
% so this only ever fires for a genuinely curved QuaPar.
    tf = ~isempty(obj.Ec) && any(obj.Ec(:) ~= 0);
end

function p = asQuaPol(obj)
% convEnvCPLQ documents its input as a QuaPol; a polyhedral QuaPar carries the identical mesh, so
% rebuild it as one rather than relying on convEnvCPLQ tolerating the other class.
    if isa(obj, 'QuaPol')
        p = obj;
    else
        p = QuaPol(obj.V, obj.E, obj.f, obj.F);
    end
end

% ------------------------------------------------------------------------------------------------
function tf = convexEnough(obj)
% Is f convex -- either PROVEN cheaply, or ASSERTED by the caller and not contradicted?
% Answers false whenever it cannot say yes; a false answer only costs the caller the short-cut.
    tf = false;
    [allPSD, readable] = everyPieceConvex(obj);
    if ~readable
        return                      % cubic, or coefficients we cannot read: say nothing
    end

    asserted = ~isempty(obj.fIsConvex) && any(obj.fIsConvex);
    if asserted && ~allPSD
        % THE FREE NECESSARY CONDITION, and it fails. Refuse loudly: honouring the flag here
        % would return a non-convex f as its own convex envelope, silently.
        error('PLQ:biconj:notConvexDespiteFlag', ...
              ['fIsConvex is set, but some piece has a Hessian that is not positive ' ...
               'semidefinite, so f cannot be convex. Clear the flag or fix the input.']);
    end
    if ~allPSD
        return
    end

    % Every piece is convex. With ONE piece on a convex domain that is already convexity -- no
    % edges to be inconsistent across -- so this needs no flag and assumes nothing.
    if obj.nf == 1
        try
            tf = logical(obj.dom.isConvex);
        catch
            tf = false;             % domain convexity not recorded: do not guess
        end
        return
    end

    % SEVERAL pieces: per-piece convexity is necessary and NOT sufficient. Only the caller's
    % assertion gets us the rest of the way.
    tf = asserted;
end

% ------------------------------------------------------------------------------------------------
function [tf, readable] = everyPieceConvex(obj)
% Is every piece's quadratic convex? `readable` is false when the coefficients are not a
% quadratic this can classify, in which case tf means nothing.
    tf = false; readable = false;
    try
        if obj.nf < 1, return, end
        for k = 1:obj.nf
            [~, Q, C] = QuaPol.matrixForm(obj.f(k,:));
            if sum(abs(C), 'all') > sqrt(eps)
                return              % genuinely cubic: outside this test
            end
            if ~QuaPol.isPositiveSemidefinite(Q)
                readable = true;    % readable and NOT convex -- a definite no
                return
            end
        end
    catch
        return
    end
    tf = true; readable = true;
end

% ------------------------------------------------------------------------------------------------
function [ok, coefs] = separableEnvelopeCoefs(obj)
% The convex envelope of a SEPARABLE quadratic over an axis-aligned BOX, as a coefficient row in
% this class's own cubic basis. ok is false whenever the hypothesis does not hold, and the caller
% falls through to the general path.
%
% The basis, from RatPar's weight vector [1/6 1/2 1/2 1/6 1/2 1 1/2 1 1 1]: columns 5..10 give
%       f = c5/2 * x^2 + c6 * x*y + c7/2 * y^2 + c8 * x + c9 * y + c10.
    ok = false; coefs = [];
    try
        if obj.nf ~= 1 || obj.nv ~= 4 || obj.ne ~= 4, return, end
        if ~obj.isDomBounded, return, end
        if any(obj.Ec(:) ~= 0), return, end        % a curved edge: not a box
        c = obj.f(1,:);
        if any(abs(c(1:4)) > sqrt(eps)), return, end   % genuinely cubic
        if abs(c(6)) > sqrt(eps), return, end          % a CROSS TERM: not separable
    catch
        return
    end

    % THE DOMAIN MUST BE A BOX: four vertices at the corners of an axis-aligned rectangle.
    V = obj.V;
    xs = unique(round(V(:,1), 12));
    ys = unique(round(V(:,2), 12));
    if numel(xs) ~= 2 || numel(ys) ~= 2, return, end
    if size(V,1) ~= 4, return, end
    for i = 1:4
        if ~any(abs(V(i,1) - xs) < 1e-12) || ~any(abs(V(i,2) - ys) < 1e-12), return, end
    end
    lo = [min(xs), min(ys)]; hi = [max(xs), max(ys)];

    ax = c(5)/2; ay = c(7)/2; dx = c(8); dy = c(9); k0 = c(10);
    [ax2, dx2, kx] = oned(ax, dx, lo(1), hi(1));
    [ay2, dy2, ky] = oned(ay, dy, lo(2), hi(2));

    coefs = zeros(1,10);
    coefs(5)  = 2*ax2;
    coefs(7)  = 2*ay2;
    coefs(8)  = dx2;
    coefs(9)  = dy2;
    coefs(10) = k0 + kx + ky;
    ok = true;
end

% ------------------------------------------------------------------------------------------------
function [a2, d2, k2] = oned(a, d, l, h)
% The convex envelope of a*t^2 + d*t on [l,h], returned as a2*t^2 + d2*t + k2.
    if a >= -sqrt(eps)
        a2 = a; d2 = d; k2 = 0;             % already convex on the interval
    else
        phiL = a*l^2 + d*l;
        phiH = a*h^2 + d*h;
        m = (phiH - phiL)/(h - l);          % the CHORD: a concave function's envelope
        a2 = 0; d2 = m; k2 = phiL - m*l;
    end
end
