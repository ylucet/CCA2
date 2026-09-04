function g = conjQ(obj)
% conjQ  The Legendre-Fenchel conjugate computed EXACTLY over the rationals, returning a QuaCon.
%
% objective: f*(s) = sup_x <s,x> - f(x), with every coefficient an exact rational and every
%   decision made by exact integer arithmetic. This is the sym-free, tolerance-free counterpart of
%   conjCPLQ, and the route CONJ_FIELD_PROOF.md 8.6/8.7 recommends.
%
% [input]  obj : QuaPol, operable (degree <= 2), and EXACT -- it must carry fN/fD, which a QuaPol
%                built from data a caller wrote down always does. One built by the legacy float
%                pipeline out of computed doubles does not, and is refused by name: computing
%                exactly on coefficients that are already one ULP wrong yields exactly the wrong
%                number, which carries no warning (QuaPol.assertExact has the reference).
% [output] g   : QuaCon -- exact faces, primitive integer edge conics, NAMED vertices.
%
% ------------------------------------------------------------------------------------------------
% EXACT OR A NAMED REFUSAL. NEVER A FALLBACK.
%
% Where this routine cannot answer it raises, by name, and the caller may then ask for the slow
% symbolic engine DELIBERATELY (conj(f,'symbolic')). It never falls back on its own. The reason is
% recorded and measured: `isAlways` returns *unknown* rather than merely slow, the symbolic form has
% no canonical shape so equal quantities compare unequal, `solve`/`simplify` cannot be interrupted
% from M-code, and one input took 43 minutes. A silent fallback turns all of that into a mystery
% about why a run is slow; a named refusal turns it into one line in SUPPORT_MATRIX.md.
%
% WHAT IS IMPLEMENTED, AND WHAT IS NOT YET
%
%   * Case A, full-domain strictly convex quadratic -> full-domain quadratic. Exact.
%   * Everything else raises PLQ:conjQ:notImplemented, naming the case.
%
% The order the remaining cases should land in is the plan's Phase 2a: the per-piece closed forms
% (convEnvCPLQ, conjPieceCPLQ, conjConvexOverPiece, conjConvexPolygon, conjAffinePLQ) are ALREADY
% 100% sym-free -- `checkConjSymFree` measures that -- so porting them is replacing double
% arithmetic with ratQ calls, not rewriting an algorithm. Step 3 (maxQuaPar) is the large item and
% comes after.
%
% ONE TOLERANCE IS REMOVED HERE, AND IT IS NOT COSMETIC. conjCPLQ decides strict convexity with
% `all(eig(Q) > sqrt(eps))`. That is a floating-point test with a fixed absolute threshold on a
% quantity that scales with the problem, so it REFUSES genuinely strictly convex quadratics whose
% smaller eigenvalue is small -- Q = [1 16384; 16384 268435457] has determinant exactly 1 and is
% positive definite, and `eig` puts its smaller eigenvalue near 3.7e-9, below the threshold. Here
% the test is Sylvester's criterion on exact integers (a > 0 and det > 0), which is DECIDED. The
% pin is conjQTest/theExactTestAcceptsAStrictlyConvexQuadraticEigWouldRefuse.

    if ~isa(obj, 'QuaPol')
        error('PLQ:conjQ:input', ...
            'conjQ takes a QuaPol (quadratic on polyhedral); got a %s.', class(obj));
    end
    obj.assertOperable();
    obj.assertExact();          % refuses a QuaPol the legacy float pipeline computed -- see there

    fN = obj.fN;  fD = obj.fD;  nf = size(fN, 1);

    % ---- Case A: full domain, one face -------------------------------------------------------
    if obj.nv == 0 && nf == 1
        g = caseAFullDomain(fN, fD);
        return
    end

    error('PLQ:conjQ:notImplemented', ...
        ['the exact conjugate is implemented for a full-domain quadratic only; this input has ' ...
         '%d vertices and %d faces. Use conj(f,''symbolic'') to reach the slow engine ' ...
         'deliberately, or extend conjQ -- see its header for the order the cases should land.'], ...
        obj.nv, nf);
end

% ==================================================================================================

function g = caseAFullDomain(fN, fD)
% objective: the conjugate of f(x) = 1/2 x'H x + L'x + kappa over ALL of R^2, exactly.
%
% For H positive definite the sup is attained at the unique stationary point x = H^-1 (s - L), so
%       f*(s) = 1/2 (s-L)' H^-1 (s-L) - kappa,
% again a full-domain quadratic -- so the answer is ONE face with NO edges and NO vertices, which
% is why the QuaCon below has empty E/EcQ/F and an empty constraint list. An empty constraint list
% is not a degenerate case here: a face with no sign conditions IS the whole plane.
%
% STAYING INTEGRAL. With H = Hn/fD, L = Ln/fD and kappa = fN(10)/fD, writing A = adj(Hn) and
% D = det(Hn) gives H^-1 = fD*A/D, and then
%       quadratic part   1/2 s' (fD A / D) s          -> weighted slots c5,c6,c7 = fD*A/D
%       linear part      -(H^-1 L)  = -(A Ln)/D
%       constant         1/2 L'H^-1 L - kappa = (Ln' A Ln - 2 D fN(10)) / (2 fD D)
% so every slot sits over the common denominator 2*fD*D and `canon` reduces it.

    Hn = [fN(5) fN(6); fN(6) fN(7)];
    Ln = [fN(8); fN(9)];
    D  = ratQ.detExact(Hn);

    % Sylvester's criterion, on exact integers: H > 0 iff Hn(1,1) > 0 and det(Hn) > 0 (fD > 0
    % always, since ratQ.canon normalises the denominator positive, so H and Hn have the same
    % definiteness). DECIDED, not estimated -- see this file's header for what the tolerance test
    % it replaces gets wrong.
    if ~(Hn(1,1) > 0 && D > 0)
        error('PLQ:conjQ:notStrictlyConvex', ...
            ['the full-domain quadratic is not strictly convex (leading minor %d, determinant ' ...
             '%d, both of which must be positive). Its conjugate is not a full-domain quadratic ' ...
             '-- for the affine, rank-deficient and indefinite cases the mathematics is three ' ...
             'lines each and the obstacle is representational; conjCPLQ.m classifies them.'], ...
            Hn(1,1), D);
    end

    A   = [Hn(2,2), -Hn(1,2); -Hn(1,2), Hn(1,1)];    % adjugate of a 2x2, exactly
    ALn = A * Ln;
    den = ratQ.chk(2 * fD * D, 'conjugate denominator');

    c5  = ratQ.chk(2 * fD^2 * A(1,1), 'conjugate c5');
    c6  = ratQ.chk(2 * fD^2 * A(1,2), 'conjugate c6');
    c7  = ratQ.chk(2 * fD^2 * A(2,2), 'conjugate c7');
    c8  = ratQ.chk(-2 * fD * ALn(1),  'conjugate c8');
    c9  = ratQ.chk(-2 * fD * ALn(2),  'conjugate c9');
    c10 = ratQ.chk(Ln' * ALn - 2 * D * fN(10), 'conjugate c10');

    [gN, gD] = ratQ.canon([0 0 0 0, c5, c6, c7, c8, c9, c10], den);

    g = QuaCon(zeros(0,3), zeros(0,3), zeros(0,6), gN, gD, zeros(0,2), {zeros(0,2)});
end
