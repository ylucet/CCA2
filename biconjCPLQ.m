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

    % ---- Case B: a single bounded triangle -> Step 1's convex envelope ----------------------
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
