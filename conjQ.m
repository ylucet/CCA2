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
%   * Case B, one strictly convex quadratic on a BOUNDED convex polygon -> the KKT active-set
%     subdivision: one affine cell per polygon vertex, one quadratic cell per edge, one for the
%     interior. Every cell boundary is a straight line, so the answer is a QuaCon whose edge
%     conics are all lines. Exact.
%   * Case D, EVERY OTHER quadratic on a BOUNDED convex polygon -- indefinite, concave, affine,
%     negative definite or PSD-singular. When H is not positive definite the sup is attained on
%     the BOUNDARY, so the answer is the max of the vertex affines and, for each edge of positive
%     curvature, that edge's clamped one-dimensional maximum. The concave case falls out of this
%     with no qualifying edge, so it is one branch rather than two.
%   * A MULTI-FACE input is the fold of the above over its pieces, via maxQ.
%   * What is still refused by name: an UNBOUNDED piece (a wedge or half-strip -- the sup can be
%     +inf, so dom f* stops being the whole plane) and a domain of dimension < 2.
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

    % ---- one conjugate per PIECE, then the fold ------------------------------------------------
    % f is q_k on P_k, so f*(s) = max_k sup_{x in P_k} <s,x> - q_k(x) = max_k (q_k + I_{P_k})*(s).
    % Conjugation turns a union into a MAX (ALGORITHM.md), which is why the per-piece closed forms
    % compose at all -- and why Step 3 is a fold rather than a special algorithm.
    g = conjPieceQ(obj, 1);
    for k = 2:nf
        g = maxQ(g, conjPieceQ(obj, k));
    end
end

% ==================================================================================================

function g = conjPieceQ(obj, k)
% objective: the exact conjugate of ONE piece -- q_k restricted to its own face P_k -- as a QuaCon.
% [input]  obj : the exact QuaPol; k : face index
%
% Classifying by the EXACT Hessian rather than by eig against a threshold is the whole point: the
% three branches below are decided by two integer comparisons, and conjCPLQ's floating-point
% version of the same decision is what silently claims an empty domain on a badly scaled input
% (DECISIONS.md 2026-09-03).
    [fN, fD] = obj.faceQ(k);
    Hn = [fN(5) fN(6); fN(6) fN(7)];
    D  = ratQ.detExact(Hn);

    if obj.nv == 0
        g = caseAFullDomain(fN, fD);               % the whole plane: no domain constraints at all
        return
    end

    % ---- an UNBOUNDED piece takes a different route, and only the concave/affine half of it ----
    [vs, rays, bounded] = pieceGeometry(obj, k);
    if ~bounded
        if Hn(1,1) <= 0 && Hn(2,2) <= 0 && D >= 0
            [ViU, vdU] = ratQ.combineDen(obj.VN, obj.VD);
            g = caseEUnboundedConcave(fN, fD, ViU, vdU, vs, rays);
            return
        end
        error('PLQ:conjQ:unbounded', ...
            ['piece %d is unbounded and its quadratic is not concave or affine (leading minor ' ...
             '%d, second diagonal %d, determinant %d). On an unbounded piece the finiteness of ' ...
             'the sup is a QUADRATIC condition on the recession CONE rather than a linear one -- ' ...
             'a real case, and a separate one.'], k, Hn(1,1), Hn(2,2), D);
    end

    [Vi, vd, cyc] = polygonExactly(obj, k);
    % A CLEAN DICHOTOMY, and it is the whole case analysis for a bounded piece. If H is positive
    % DEFINITE the objective <s,x> - q(x) is strictly concave and its maximiser can be interior,
    % so the KKT active set is needed. Otherwise the maximiser can always be taken on the
    % BOUNDARY -- caseDBoundaryMax's header proves it -- and the answer is a max over the edges.
    % Sylvester's criterion decides which, on exact integers.
    if Hn(1,1) > 0 && D > 0
        g = caseBConvexOnPolygon(fN, fD, Vi, vd, cyc);      % strictly convex: vertex/edge/interior
    else
        g = caseDBoundaryMax(fN, fD, Vi, vd, cyc);          % everything else: max over the boundary
    end
end

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

function g = caseBConvexOnPolygon(fN, fD, Vi, vd, cyc)
% objective: the conjugate of a STRICTLY CONVEX quadratic over one bounded convex polygon,
%            exactly, as a QuaCon whose faces are the KKT active-set cells.
%
% THE DECOMPOSITION is conjConvexOverPiece's, and its header states it: q* (s) = max{<s,x> - q(x)
% : x in P} is a concave program, so the maximiser lies on exactly one relatively open face of P,
% and which one is decided by where s - grad q(x*) sits. What is new here is only the ARITHMETIC:
% that routine reaches for `sym` ten times to keep its values exact while deciding combinatorics in
% doubles, and every one of those calls becomes integer arithmetic below.
%
%   VERTEX v      x* = v, valid iff s - grad q(v) lies in the normal cone N_P(v), i.e.
%                 <s - grad q(v), e> <= 0 for every edge direction e leaving v.
%                 Value <s,v> - q(v): AFFINE in s.
%
%   EDGE (base a, direction d, outward normal n)
%                 the derivative along d vanishes at t* = <s - grad q(a), d> / (d'Q d), which is
%                 affine in s, and the multiplier condition is <s - grad q(x*), n> >= 0. With
%                 alpha = d'Qd and u = s - grad q(a):
%                       value = <s,a> - q(a) + <u,d>^2 / (2 alpha),   QUADRATIC in s
%                 and the cell is 0 <= t* <= 1 together with the multiplier sign -- all affine.
%
%   INTERIOR      x* = Q^-1 (s - L), valid iff x* lies in P, value 1/2 (s-L)'Q^-1(s-L) - c.
%
% So EVERY cell boundary is a straight line and every face function is a rational quadratic: the
% answer is a QuaCon all of whose edge conics are lines. That is the shape the SCIP bridge wants.
%
% WHAT IS BUILT AND WHAT IS DEFERRED. The faces are complete: exact functions, canonical line
% conics, and the full H-form sign conditions, which is everything `eval` and every predicate
% needs. The INCIDENCE arrays E and F are left empty, and that is stated rather than faked --
% recovering which cell borders which along which segment is an arrangement computation that
% nothing consumes yet, and inventing a plausible-looking incidence would be worse than an honest
% gap. Vertices ARE named, since a corner of a cell is just two of its own bounding lines meeting.

    Hn = [fN(5) fN(6); fN(6) fN(7)];
    Ln = [fN(8); fN(9)];
    kn = fN(10);
    D  = ratQ.detExact(Hn);
    m  = numel(cyc);
    % AN INTERNAL INVARIANT, not a reachable gap: conjQ's dispatch already classified H exactly
    % and only routes a strictly convex one here. Kept because this function's contract is stated
    % in terms of its own input, and a future caller that forgets should be told loudly.
    if ~(Hn(1,1) > 0 && D > 0)
        error('PLQ:conjQ:notStrictlyConvex', ...
            ['Case B needs a STRICTLY convex quadratic (leading minor %d, determinant %d; both ' ...
             'must be positive). A merely semidefinite Q makes the interior cell degenerate -- ' ...
             'the unconstrained sup is finite only on a measure-zero set -- which is a real case ' ...
             'but a different one.'], Hn(1,1), D);
    end

    % Everything below is built over these two denominators, and only at the end is each cell's
    % function reduced by `canon`. q(x) = (1/2 x'Hn x + Ln'x + kn)/fD with x = Vi/vd.

    cells = struct('num', {}, 'den', {}, 'con', {});   % con: k x 7 rows [a b c d e f  sign]

    % ---- one cell per VERTEX -----------------------------------------------------------------
    for i = 1:m
        v  = Vi(cyc(i), :).';                        % over vd
        gv = ratQ.chk(Hn * v + vd * Ln, 'vertex gradient');   % grad q(v) = (Hn v/vd + Ln)/fD
                                                              % numerator over fD*vd
        % q(v) over 2*fD*vd^2
        qv = ratQ.chk(v.' * Hn * v + 2 * vd * (Ln.' * v) + 2 * vd^2 * kn, 'vertex value');

        % value  <s,v> - q(v)  =  (v1 s1 + v2 s2)/vd - qv/(2 fD vd^2)
        %        over the common denominator 2*fD*vd^2
        num = [0 0 0 0, 0, 0, 0, ...
               ratQ.chk(2*fD*vd*v(1), 'c8'), ratQ.chk(2*fD*vd*v(2), 'c9'), -qv];
        den = ratQ.chk(2 * fD * vd^2, 'vertex cell denominator');

        % normal cone: <s - grad q(v), e> <= 0 for each edge direction e leaving v
        con = zeros(0,7);
        dirs = edgeDirsAtVertex(Vi, cyc, i);
        for kk = 1:numel(dirs)
            ed = dirs{kk};                            % integer direction, over vd
            % <s,e> - <grad q(v), e> <= 0, cleared by fD*vd:
            %   fD*vd*(e1 s1 + e2 s2) - <gv, e> <= 0        (gv is over fD*vd)
            row = [0 0 0, ratQ.chk(fD*vd*ed(1),'n1'), ratQ.chk(fD*vd*ed(2),'n2'), ...
                   ratQ.chk(-(gv.' * ed), 'n0')];
            % sgnOf, NOT a hardcoded -1: ratQ.conic normalises the row's overall sign and may
            % have negated it, in which case "<= 0" now names the OTHER half-plane and the cone
            % becomes its own complement -- a wrong answer that still builds a plausible mesh.
            con(end+1, :) = [ratQ.conic(row), sgnOf(row(4:6), -1)]; %#ok<AGROW>
        end
        cells(end+1) = struct('num', num, 'den', den, 'con', con); %#ok<AGROW>
    end

    % ---- one cell per EDGE --------------------------------------------------------------------
    for i = 1:m
        a = Vi(cyc(i), :).';
        b = Vi(cyc(mod(i, m) + 1), :).';
        d = ratQ.chk(b - a, 'edge direction');                 % over vd
        alpha = ratQ.chk(d.' * Hn * d, 'edge curvature');      % d'Hn d, over 1 (d is integral)
        if alpha <= 0
            % q is affine along this edge, so its max sits at an endpoint and the edge carries no
            % 2-D cell. Cannot happen for a strictly convex q with d ~= 0, and is kept as an
            % assertion rather than a branch.
            error('PLQ:conjQ:degenerateEdge', ...
                'edge %d has d''Qd = %d <= 0 under a strictly convex Q, which is impossible.', ...
                i, alpha);
        end
        ga = ratQ.chk(Hn * a + vd * Ln, 'edge base gradient');  % over fD*vd
        qa = ratQ.chk(a.' * Hn * a + 2 * vd * (Ln.' * a) + 2 * vd^2 * kn, 'edge base value');
        n  = outwardNormal(Vi, cyc, i);                          % integer, over vd (scale free)

        % t* = <s - grad q(a), d> / alpha, with <s,d> cleared by fD*vd:
        %      T(s) := fD*vd*(d1 s1 + d2 s2) - <ga,d>   satisfies   t* = T(s) / alpha
        % (d and a are the INTEGER numerators over vd, so the actual direction is d/vd and the
        %  actual curvature is alpha/(fD*vd^2); the two vd's and the fD cancel to exactly this.)
        Td = [ratQ.chk(fD*vd*d(1),'t1'), ratQ.chk(fD*vd*d(2),'t2'), ratQ.chk(-(ga.'*d),'t0')];

        % value = <s,a> - q(a) + <u,d_actual>^2 / (2 alpha_actual), and that last term works out
        % to T(s)^2 / (2 alpha fD vd^2), while <s,a>/vd and q(a) = qa/(2 fD vd^2). So the common
        % denominator is 2*alpha*fD*vd^2 and the three contributions scale by alpha, 1 and 1.
        den = ratQ.chk(2 * alpha * fD * vd^2, 'edge cell denominator');
        lin = [ratQ.chk(2*alpha*fD*vd*a(1), 'c8'), ratQ.chk(2*alpha*fD*vd*a(2), 'c9')];
        num = [0 0 0 0, ...
               ratQ.chk(2*Td(1)^2, 'c5'), ratQ.chk(2*Td(1)*Td(2), 'c6'), ratQ.chk(2*Td(2)^2, 'c7'), ...
               ratQ.chk(lin(1) + 2*Td(1)*Td(3), 'c8'), ratQ.chk(lin(2) + 2*Td(2)*Td(3), 'c9'), ...
               ratQ.chk(Td(3)^2 - alpha*qa, 'c10')];

        con = zeros(0,7);
        % 0 <= t*  ->  T(s) >= 0
        con(end+1,:) = [ratQ.conic([0 0 0, Td(1), Td(2), Td(3)]), sgnOf(Td, +1)]; %#ok<AGROW>
        % t* <= 1   ->  T(s) - alpha <= 0
        Tu = [Td(1), Td(2), ratQ.chk(Td(3) - alpha, 'edge upper')];
        con(end+1,:) = [ratQ.conic([0 0 0, Tu(1), Tu(2), Tu(3)]), sgnOf(Tu, -1)]; %#ok<AGROW>
        % multiplier: <s - grad q(a), n> - t* <Hn d / fD, n> >= 0, cleared to integers.
        %   <s - grad q(a), n> is (fD*vd*<s,n> - <ga,n>)/(fD*vd), and
        %   t* <Q d, n> = T(s) * (d' Hn n) / (fD*vd*alpha*fD) ... clear everything by fD*vd*alpha:
        Hdn = ratQ.chk(d.' * Hn * n, 'edge multiplier coupling');
        Mu  = [ratQ.chk(alpha*fD*vd*n(1) - Td(1)*Hdn, 'm1'), ...
               ratQ.chk(alpha*fD*vd*n(2) - Td(2)*Hdn, 'm2'), ...
               ratQ.chk(-alpha*(ga.'*n) - Td(3)*Hdn, 'm0')];
        con(end+1,:) = [ratQ.conic([0 0 0, Mu(1), Mu(2), Mu(3)]), sgnOf(Mu, +1)]; %#ok<AGROW>

        cells(end+1) = struct('num', num, 'den', den, 'con', con); %#ok<AGROW>
    end

    % ---- the INTERIOR cell ---------------------------------------------------------------------
    % x* = Q^-1(s-L) = fD*A*(s - L/fD)/D with A = adj(Hn); membership is <n_j, x*> <= <n_j, v_j>
    % for every facet j, cleared to integers.
    A = [Hn(2,2), -Hn(1,2); -Hn(1,2), Hn(1,1)];
    den = ratQ.chk(2 * fD * D, 'interior cell denominator');
    ALn = A * Ln;
    numI = [0 0 0 0, ratQ.chk(2*fD^2*A(1,1),'c5'), ratQ.chk(2*fD^2*A(1,2),'c6'), ...
            ratQ.chk(2*fD^2*A(2,2),'c7'), ratQ.chk(-2*fD*ALn(1),'c8'), ...
            ratQ.chk(-2*fD*ALn(2),'c9'), ratQ.chk(Ln.'*ALn - 2*D*kn, 'c10')];
    conI = zeros(0,7);
    for i = 1:m
        n = outwardNormal(Vi, cyc, i);
        a = Vi(cyc(i), :).';
        % <n, fD*A*(s - L/fD)/D>  <=  <n, a/vd>
        %   multiply by D*fD*vd (positive, since D > 0 and fD, vd > 0):
        %   vd*fD*<n, A s> - vd*<n, A Ln>  -  D*fD*<n,a>  <=  0
        An = A.' * n;
        row = [0 0 0, ratQ.chk(vd*fD*An(1), 'i1'), ratQ.chk(vd*fD*An(2), 'i2'), ...
               ratQ.chk(-vd*(An.'*Ln) - D*(n.'*a), 'i0')];
        conI(end+1,:) = [ratQ.conic(row), sgnOf(row(4:6), -1)]; %#ok<AGROW>
    end
    cells(end+1) = struct('num', numI, 'den', den, 'con', conI);

    g = assembleQuaConCells(cells);
end

% --------------------------------------------------------------------------------------------

function s = sgnOf(row, want)
% objective: the H-form sign for a constraint written as `row <= 0` (want = -1) or `row >= 0`
%            (want = +1), AFTER ratQ.conic has normalised the row's overall sign.
%
% ratQ.conic makes the first nonzero entry positive, which may have negated the whole row. The
% face's side must follow that negation, or the cell silently becomes its own complement -- which
% is a wrong answer that still produces a plausible mesh.
%
% Takes either the three LINEAR coefficients [d e f] or a full conic [a b c d e f]; the caller has
% one or the other depending on whether the constraint came from a half-plane or from a difference
% of two quadratics, and getting the padding wrong would look at the wrong leading entry.
    if numel(row) == 3
        r = [0 0 0, row(1), row(2), row(3)];
    else
        r = row;
    end
    nz = find(r ~= 0, 1);
    if isempty(nz)
        error('ratQ:zeroConic', 'a constraint with all-zero coefficients names no half-plane.');
    end
    if r(nz) < 0, s = -want; else, s = want; end
end

function [Vi, vd, cyc] = polygonExactly(obj, k)
% objective: FACE k's vertices as integers over ONE denominator, plus the cycle order.
% [output] Vi : nv x 2 integer (the whole mesh); vd : positive integer;
%          cyc : 1 x m indices into Vi, in boundary order around face k
%
% The cycle is walked from E rather than taken from P, so this does not depend on P's clockwise/
% counter-clockwise convention -- a convention it would be easy to read backwards, and reading it
% backwards flips every outward normal, which turns every cell into its complement.
    [Vi, vd] = ratQ.combineDen(obj.VN, obj.VD);
    own = find(any(obj.F == k, 2));
    E = obj.E(own, :);
    if isempty(E)
        error('PLQ:conjQ:noFace', 'face %d has no edges.', k);
    end
    % An INTERNAL INVARIANT since caseE landed: conjPieceQ tests boundedness first and routes an
    % unbounded piece elsewhere, so this cannot fire from there. Kept because the routine's own
    % contract is about its input -- walking a cycle that does not close would otherwise build a
    % silently wrong polygon rather than raise.
    if any(E(:,3) == 0)
        error('PLQ:conjQ:unbounded', ...
            ['the exact conjugate needs a BOUNDED piece; face %d has a ray (edge %d). The ' ...
             'unbounded polyhedral pieces -- wedges and half-strips -- are a real case and a ' ...
             'separate one: the sup can be +inf, so dom f* stops being the whole plane.'], ...
            k, own(find(E(:,3) == 0, 1)));
    end
    m = size(E,1);
    adj = zeros(size(Vi,1), 2);
    for j = 1:m
        for c = 1:2
            v = E(j,c);
            if adj(v,1) == 0, adj(v,1) = E(j,3-c); else, adj(v,2) = E(j,3-c); end
        end
    end
    cyc = zeros(1, m);
    cyc(1) = E(1,1);  cyc(2) = E(1,2);
    for k = 3:m
        nxt = adj(cyc(k-1), :);
        cyc(k) = nxt(nxt ~= cyc(k-2));
    end
    if numel(unique(cyc)) ~= m
        error('PLQ:conjQ:notSimple', ...
            'the boundary of face %d is not a single simple cycle of %d edges.', k, m);
    end
end

function dirs = edgeDirsAtVertex(Vi, cyc, i)
% objective: the two edge directions leaving cycle position i, as integer column vectors.
    m = numel(cyc);
    v  = Vi(cyc(i), :).';
    p  = Vi(cyc(mod(i-2, m) + 1), :).';
    n  = Vi(cyc(mod(i,   m) + 1), :).';
    dirs = {ratQ.chk(n - v, 'edge direction'), ratQ.chk(p - v, 'edge direction')};
end

function n = outwardNormal(Vi, cyc, i)
% objective: the OUTWARD normal of the edge from cycle position i to i+1, as an integer vector.
%
% The direction is fixed by testing against the polygon's own centroid, exactly: the outward side
% is the one the centroid is NOT on. Orientation conventions are not consulted, so this cannot be
% read backwards -- and reading it backwards would turn every cell into its complement.
    m = numel(cyc);
    a = Vi(cyc(i), :).';
    b = Vi(cyc(mod(i, m) + 1), :).';
    d = b - a;
    n = [d(2); -d(1)];
    ctr = sum(Vi(cyc, :), 1).';                 % m * centroid, same denominator
    if ratQ.chk(n.' * (ctr - m*a), 'centroid side') > 0
        n = -n;
    end
end

function g = maxOverVerticesQuaCon(fN, fD, Vi, vd, cyc)
% objective: the max of the AFFINE functions <s,v> - q(v) over the polygon's vertices, as a QuaCon.
%
% RENAMED 2026-09-04 from caseCConcaveOnPolygon, because the construction was never about
% concavity: it is the max of finitely many affine functions and is meaningful for any q. What
% concavity supplied was the THEOREM that this max IS the conjugate -- and caseDBoundaryMax now
% owns that reasoning, in the more general form that also covers the indefinite and PSD-singular
% cases. For a concave q that routine finds no qualifying edge and returns exactly this object, so
% the concave case is unchanged, not merely similar.
%
% WHY THIS CASE IS EASY, AND IT IS THE MATHEMATICS THAT MAKES IT SO. If q is concave then -q is
% convex, so the objective <s,x> - q(x) is CONVEX in x, and a convex function attains its maximum
% over a polytope at an EXTREME POINT. Hence
%       q*(s) = max_i [ <s, v_i> - q(v_i) ],
% a maximum of finitely many AFFINE functions of s -- no normal cones, no multipliers, no interior
% cell, and no curvature anywhere in the answer. The affine case is the same statement with a zero
% Hessian, which is why one branch serves both.
%
% This is the conjugate half of what ALGORITHM.md calls the concave envelope: co q on the polygon
% is the affine interpolant over the LOWER HULL of the lifted vertices, and the cells below are
% exactly that hull's normal fan. The vertices NOT on the lower hull are the ones whose cell comes
% out empty -- so assembleQuaCon's feasibility filter is doing real work here rather than tidying,
% and that is why it had to be exact.
    Hn = [fN(5) fN(6); fN(6) fN(7)];
    Ln = [fN(8); fN(9)];
    kn = fN(10);
    m  = numel(cyc);

    % q(v_i), over the common denominator 2*fD*vd^2 -- the same clearing as Case B's vertex cells
    qv = zeros(m,1);
    for i = 1:m
        v = Vi(cyc(i), :).';
        qv(i) = ratQ.chk(v.'*Hn*v + 2*vd*(Ln.'*v) + 2*vd^2*kn, 'vertex value');
    end

    den = ratQ.chk(2 * fD * vd^2, 'cell denominator');
    cells = struct('num', {}, 'den', {}, 'con', {});
    for i = 1:m
        vi = Vi(cyc(i), :).';
        num = [0 0 0 0, 0, 0, 0, ...
               ratQ.chk(2*fD*vd*vi(1), 'c8'), ratQ.chk(2*fD*vd*vi(2), 'c9'), -qv(i)];

        % cell i:  <s,v_i> - q(v_i) >= <s,v_j> - q(v_j)  for every other vertex j, cleared by the
        % same positive 2*fD*vd^2:
        %       2*fD*vd*<s, Vi_i - Vi_j>  -  (qv_i - qv_j)  >= 0
        con = zeros(0,7);
        for j = 1:m
            if j == i, continue, end
            vj = Vi(cyc(j), :).';
            row = [0 0 0, ratQ.chk(2*fD*vd*(vi(1)-vj(1)), 'd1'), ...
                          ratQ.chk(2*fD*vd*(vi(2)-vj(2)), 'd2'), ...
                          ratQ.chk(-(qv(i) - qv(j)), 'd0')];
            if all(row(4:6) == 0)
                % Two DISTINCT vertices cannot give a zero direction, so this can only be reached
                % if the polygon carried a repeated vertex. Refuse rather than skip: a repeated
                % vertex means the caller's mesh is malformed, and silently ignoring it would put
                % a duplicated affine piece into the answer.
                error('PLQ:conjQ:repeatedVertex', ...
                    'polygon vertices %d and %d coincide, so the domain is not a simple polygon.', ...
                    cyc(i), cyc(j));
            end
            con(end+1,:) = [ratQ.conic(row), sgnOf(row(4:6), +1)]; %#ok<AGROW>
        end
        cells(end+1) = struct('num', num, 'den', den, 'con', con); %#ok<AGROW>
    end

    g = assembleQuaConCells(cells);
end

function g = caseDBoundaryMax(fN, fD, Vi, vd, cyc)
% objective: the conjugate of a quadratic that is NOT positive definite, over a bounded convex
%            polygon -- indefinite, negative definite, concave, affine, or PSD-singular.
%
% THE ONE FACT THAT MAKES THIS CASE TRACTABLE, and it replaces the whole A.2-A.5 apparatus for the
% conjugate: when H is not positive DEFINITE, the sup of <s,x> - q(x) over P is attained on the
% BOUNDARY of P. Suppose it were attained at an interior point x0. Then x0 is a local maximum of
% the objective, so the objective's Hessian -H must be negative semidefinite there, i.e. H is PSD.
% If H is PSD but SINGULAR the argument still gives the boundary: there is a direction d with
% d'Hd = 0, along which the objective is affine, so from any interior maximiser one can travel to
% the boundary without decreasing. Only a positive DEFINITE H makes the interior stationary point a
% strict maximiser, and that is Case B.
%
% So q*(s) is the largest of finitely many CANDIDATES:
%
%   every vertex v          val_v(s) = <s,v> - q(v),  AFFINE in s, always available;
%   every edge of positive  the clamped stationary point of the edge's 1-D restriction, QUADRATIC
%   curvature d'Hd > 0      in s, and available only where 0 <= t* <= 1 -- which is a pair of
%                           AFFINE conditions, since t* = T(s)/alpha is affine.
%
% An edge of non-positive curvature has a convex (or affine) restriction, so its maximum sits at an
% endpoint and it contributes no candidate at all -- both endpoints are already vertices. When H is
% negative semidefinite NO edge qualifies, the candidate set is just the vertices, and this routine
% degenerates to exactly maxOverVerticesQuaCon. So the concave case is one branch of this, not a
% separate algorithm.
%
% ------------------------------------------------------------------------------------------------
% CLASSIFY FIRST, SPLIT LAST -- and here that is worth a measurement, not just a principle.
%
% The obvious implementation makes each edge's clamped maximum into a total three-piece function
% and folds them pairwise with maxQ. That is CORRECT and was the first version, but each fold
% re-splits cells the previous fold had already separated, and the splits are CONIC, so no exact
% test available today can tell that the resulting cell is empty. Measured on a PSD-singular
% pentagon: 501 faces, of which 75 were ever occupied and which carried 10 distinct functions.
% Three sound filters (linear infeasibility, contradictory sides, constant-sign conics) took it to
% 274, and exact redundancy elimination shrank the lists from 14 constraints to 6.3 without
% removing a single cell -- because what separates those cells is curved.
%
% This version removes the cause instead of filtering the symptom. EVERY availability condition is
% AFFINE, so the arrangement they cut the plane into is a LINE arrangement, whose cells
% ratQ.feasible2 decides exactly and completely. Build that arrangement first, pruning as it grows;
% only then, inside each cell, ask which of the few surviving candidates is largest -- and that is
% the only place a conic ever enters. The exponential 3^m enumeration never materialises because an
% infeasible partial assignment is dropped as soon as it is infeasible.

    Hn = [fN(5) fN(6); fN(6) fN(7)];
    Ln = [fN(8); fN(9)];
    kn = fN(10);
    m  = numel(cyc);

    % ---- the vertex candidates, and the linear cells on which each is the vertex maximum -------
    denV = ratQ.chk(2 * fD * vd^2, 'vertex denominator');
    qv = zeros(m,1);  numV = zeros(m,10);
    for i = 1:m
        v = Vi(cyc(i), :).';
        qv(i) = ratQ.chk(v.'*Hn*v + 2*vd*(Ln.'*v) + 2*vd^2*kn, 'vertex value');
        numV(i,:) = [0 0 0 0, 0, 0, 0, ...
                     ratQ.chk(2*fD*vd*v(1),'c8'), ratQ.chk(2*fD*vd*v(2),'c9'), -qv(i)];
    end

    % ---- the edge candidates ------------------------------------------------------------------
    Ed = struct('Td', {}, 'alpha', {}, 'num', {}, 'den', {});
    for j = 1:m
        a = Vi(cyc(j), :).';
        b = Vi(cyc(mod(j, m) + 1), :).';
        d = ratQ.chk(b - a, 'edge direction');
        alpha = ratQ.chk(d.' * Hn * d, 'edge curvature');
        if alpha <= 0, continue, end
        ga = ratQ.chk(Hn * a + vd * Ln, 'edge base gradient');
        qa = ratQ.chk(a.'*Hn*a + 2*vd*(Ln.'*a) + 2*vd^2*kn, 'edge base value');
        Td = [ratQ.chk(fD*vd*d(1),'t1'), ratQ.chk(fD*vd*d(2),'t2'), ratQ.chk(-(ga.'*d),'t0')];
        lin = [ratQ.chk(2*alpha*fD*vd*a(1),'c8'), ratQ.chk(2*alpha*fD*vd*a(2),'c9')];
        Ed(end+1) = struct('Td', Td, 'alpha', alpha, ...
            'num', [0 0 0 0, ...
                    ratQ.chk(2*Td(1)^2,'c5'), ratQ.chk(2*Td(1)*Td(2),'c6'), ...
                    ratQ.chk(2*Td(2)^2,'c7'), ...
                    ratQ.chk(lin(1) + 2*Td(1)*Td(3),'c8'), ...
                    ratQ.chk(lin(2) + 2*Td(2)*Td(3),'c9'), ...
                    ratQ.chk(Td(3)^2 - alpha*qa,'c10')], ...
            'den', ratQ.chk(2*alpha*fD*vd^2, 'edge denominator')); %#ok<AGROW>
    end

    % ---- the LINE arrangement, built incrementally and pruned exactly at every step -------------
    % A cell records the half-planes defining it, which vertex is the vertex maximum there, and
    % which edge candidates are available. Every condition here is affine, so feasible2 decides
    % emptiness exactly -- this is the pruning that keeps 3^(#edges) from ever being enumerated.
    cells0 = struct('lin', {}, 'vi', {}, 'ei', {});
    for i = 1:m
        rows = zeros(0,3);
        for j = 1:m
            if j == i, continue, end
            r = [ratQ.chk(2*fD*vd*(Vi(cyc(i),1)-Vi(cyc(j),1)), 'd1'), ...
                 ratQ.chk(2*fD*vd*(Vi(cyc(i),2)-Vi(cyc(j),2)), 'd2'), ...
                 ratQ.chk(-(qv(i) - qv(j)), 'd0')];
            if all(r == 0)
                error('PLQ:conjQ:repeatedVertex', ...
                    'polygon vertices %d and %d coincide.', cyc(i), cyc(j));
            end
            rows(end+1,:) = r; %#ok<AGROW>
        end
        if ratQ.feasible2(rows, true)
            cells0(end+1) = struct('lin', rows, 'vi', i, 'ei', []); %#ok<AGROW>
        end
    end

    for e = 1:numel(Ed)
        Td = Ed(e).Td;  alpha = Ed(e).alpha;
        regimes = { [-Td(1) -Td(2) -Td(3)], ...                                   % t* <= 0
                    [Td; -Td(1), -Td(2), alpha - Td(3)], ...                      % 0 <= t* <= 1
                    [Td(1) Td(2) Td(3)-alpha] };                                  % t* >= 1
        nxt = struct('lin', {}, 'vi', {}, 'ei', {});
        for c = 1:numel(cells0)
            for r = 1:3
                lin = [cells0(c).lin; regimes{r}];
                if ~ratQ.feasible2(lin, true), continue, end
                ei = cells0(c).ei;
                if r == 2, ei = [ei, e]; end
                nxt(end+1) = struct('lin', lin, 'vi', cells0(c).vi, 'ei', ei); %#ok<AGROW>
            end
        end
        cells0 = nxt;
    end

    % ---- inside each linear cell, the winner among its few candidates --------------------------
    % This is the ONLY place a conic appears: the difference of a vertex affine and an edge
    % quadratic, or of two edge quadratics.
    out = struct('num', {}, 'den', {}, 'con', {});
    for c = 1:numel(cells0)
        base = [zeros(size(cells0(c).lin,1), 3), cells0(c).lin, ones(size(cells0(c).lin,1),1)];
        for k = 1:size(base,1)
            row = ratQ.conic(base(k,1:6));
            base(k,:) = [row, sgnOf(base(k,4:6), +1)];
        end

        cand = struct('num', {}, 'den', {});
        cand(1) = struct('num', numV(cells0(c).vi,:), 'den', denV);
        for e = cells0(c).ei
            cand(end+1) = struct('num', Ed(e).num, 'den', Ed(e).den); %#ok<AGROW>
        end

        for i = 1:numel(cand)
            con = base;
            for j = 1:numel(cand)
                if j == i, continue, end
                [dn, ~] = ratQ.sub(cand(i).num, cand(i).den, cand(j).num, cand(j).den);
                if all(dn == 0), continue, end        % the same function: no boundary between them
                R = ratQ.chk([dn(5), 2*dn(6), dn(7), 2*dn(8), 2*dn(9), 2*dn(10)], 'winner');
                con(end+1,:) = [ratQ.conic(R), sgnOf(R, +1)]; %#ok<AGROW>
            end
            out(end+1) = struct('num', cand(i).num, 'den', cand(i).den, 'con', con); %#ok<AGROW>
        end
    end

    g = assembleQuaConCells(out);
end

function g = caseEUnboundedConcave(fN, fD, Vi, vd, vs, rays)
% objective: the conjugate of a CONCAVE or AFFINE quadratic over an UNBOUNDED polyhedral piece.
% [input]  fN, fD : the piece's function, exactly;  Vi, vd : the mesh vertices over one denominator
%          vs     : the indices of this piece's own vertices;  rays : k x 2 integer directions
% [output] g      : QuaCon -- and unlike every case so far, its cells DO NOT COVER THE PLANE
%
% ------------------------------------------------------------------------------------------------
% THIS IS THE FIRST CASE WHERE dom f* IS A PROPER SUBSET, and that is the mathematics, not a gap.
%
% Along a recession direction d of the piece, with a any point of it,
%       <s, a+td> - q(a+td)  =  [ <s,a> - q(a) ]  +  t <s - grad q(a), d>  -  (t^2/2) d'Hd,
% so the sup is +infinity as soon as the leading behaviour is upward. With H negative semidefinite
% the quadratic term is -(t^2/2) d'Hd >= 0, so it ALREADY diverges unless d'Hd = 0 -- and for an
% NSD H that is equivalent to Hd = 0. Two consequences, both exact tests on integers:
%
%   * if Hd is nonzero for some recession direction, f* is +infinity EVERYWHERE and dom f* is
%     empty. Refused by name: a function with no domain has no QuaCon to be stored in, which is a
%     representational gap this type shares with conjCPLQ (`conjugateHasEmptyDomain`).
%
%   * otherwise Hd = 0 makes <grad q(a), d> = <Ha + L, d> = a'Hd + <L,d> = <L,d>, independent of
%     which point of the ray is used, so the finiteness condition is
%           <s - L, d>  <=  0
%     -- LINEAR in s, and linear in d as well, which is why testing it on the EXTREME rays settles
%     it for the whole recession cone rather than only its generators.
%
% dom f* is therefore the polar-type cone { s : <s - L, d> <= 0 for every ray direction d }, and on
% it the objective is CONVEX in x (its Hessian is -H, positive semidefinite), so the maximum over
% the piece is attained at an EXTREME POINT -- and the extreme points of a polyhedron are exactly
% its vertices, unbounded or not. So the value is the same vertex max as the bounded concave case;
% all that changes is that every cell also carries the domain constraints.
%
% WHAT THIS COVERS, and it is the reason to build it before the general unbounded case: q AFFINE is
% the H = 0 instance, and that is the whole family of elementary unbounded conjugates -- an
% indicator function, a support function, a norm, max(0,x,y). Those are the inputs a user writes
% first and the ones the legacy conjAffinePLQ exists for.

    Hn = [fN(5) fN(6); fN(6) fN(7)];
    Ln = [fN(8); fN(9)];
    kn = fN(10);
    m  = numel(vs);
    if m == 0
        error('PLQ:conjQ:noVertex', ...
            ['an unbounded piece with no vertex at all is a half-plane or a slab, whose ' ...
             'conjugate is supported on a line or a point -- a domain of dimension < 2, which ' ...
             'this type cannot carry.']);
    end

    % ---- the finiteness conditions, and the refusal when there are none ----------------------
    dom = zeros(0,3);
    for r = 1:size(rays,1)
        d = rays(r,:).';
        if any(ratQ.chk(Hn * d, 'recession curvature') ~= 0)
            error('PLQ:conjQ:emptyDomain', ...
                ['the piece recedes along (%d,%d), where the objective grows without bound for ' ...
                 'every s, so f* is +infinity everywhere and dom f* is EMPTY. That is the right ' ...
                 'answer and there is nowhere to put it -- a QuaCon carries at least one face. ' ...
                 'Same representational gap conjCPLQ records as conjugateHasEmptyDomain.'], ...
                d(1), d(2));
        end
        % <s - L, d> <= 0, cleared by fD > 0 and written as a >= 0 half-plane
        row = [ratQ.chk(-fD*d(1), 'r1'), ratQ.chk(-fD*d(2), 'r2'), ratQ.chk(Ln.'*d, 'r0')];
        if all(row == 0), continue, end            % a zero direction constrains nothing
        dom(end+1,:) = row; %#ok<AGROW>
    end

    % ---- the vertex values, exactly ----------------------------------------------------------
    denV = ratQ.chk(2 * fD * vd^2, 'vertex denominator');
    qv = zeros(m,1);  numV = zeros(m,10);
    for i = 1:m
        v = Vi(vs(i), :).';
        qv(i) = ratQ.chk(v.'*Hn*v + 2*vd*(Ln.'*v) + 2*vd^2*kn, 'vertex value');
        numV(i,:) = [0 0 0 0, 0, 0, 0, ...
                     ratQ.chk(2*fD*vd*v(1),'c8'), ratQ.chk(2*fD*vd*v(2),'c9'), -qv(i)];
    end

    % ---- one cell per vertex that wins somewhere INSIDE dom f* --------------------------------
    cells = struct('num', {}, 'den', {}, 'con', {});
    for i = 1:m
        rows = dom;
        for j = 1:m
            if j == i, continue, end
            r = [ratQ.chk(2*fD*vd*(Vi(vs(i),1)-Vi(vs(j),1)), 'd1'), ...
                 ratQ.chk(2*fD*vd*(Vi(vs(i),2)-Vi(vs(j),2)), 'd2'), ...
                 ratQ.chk(-(qv(i) - qv(j)), 'd0')];
            if all(r == 0)
                error('PLQ:conjQ:repeatedVertex', ...
                    'piece vertices %d and %d coincide.', vs(i), vs(j));
            end
            rows(end+1,:) = r; %#ok<AGROW>
        end
        if ~ratQ.feasible2(rows, true), continue, end     % this vertex never attains the max
        con = zeros(size(rows,1), 7);
        for t = 1:size(rows,1)
            con(t,:) = [ratQ.conic([0 0 0, rows(t,:)]), sgnOf(rows(t,:), +1)];
        end
        cells(end+1) = struct('num', numV(i,:), 'den', denV, 'con', con); %#ok<AGROW>
    end

    g = assembleQuaConCells(cells);
end

function [vs, rays, bounded] = pieceGeometry(obj, k)
% objective: face k's own vertices and its recession directions, exactly.
% [output] vs      : the indices into obj.VN of the vertices on this face
%          rays    : r x 2 INTEGER directions (over the shared vertex denominator, but a direction
%                    is scale-free so the denominator never has to be divided out)
%          bounded : true when the face has no ray edge
%
% A ray edge is stored as [v1 v2 0] meaning the ray [V(v1), V(v2)) -- base then a second point
% giving the direction, both of them real vertices of the mesh (QuaPol.belongToEdge is where that
% convention is written down). So the direction is simply V(v2) - V(v1), and it is exact.
    own = find(any(obj.F == k, 2));
    if isempty(own)
        error('PLQ:conjQ:noFace', 'face %d has no edges.', k);
    end
    E = obj.E(own, :);
    [Vi, ~] = ratQ.combineDen(obj.VN, obj.VD);

    % A RAY EDGE STORES A DIRECTION MARKER, NOT A SECOND VERTEX. `[v1 v2 0]` means the ray
    % [V(v1), V(v2)) -- v1 is the base and V(v2) only fixes the direction (QuaPol.belongToEdge is
    % where that convention is written down). So V(v2) is a point OF the piece but not an extreme
    % point of it, and putting it in the vertex list would ask the max to consider a point that can
    % only ever tie. Harmless for the VALUE -- every point of the domain is bounded by the sup --
    % but it manufactures cells that are empty, which is the thing this pass exists to avoid.
    vs = [];
    rays = zeros(0,2);
    for j = 1:size(E,1)
        if E(j,3) ~= 0
            vs = [vs; E(j,1); E(j,2)]; %#ok<AGROW>
        else
            vs = [vs; E(j,1)]; %#ok<AGROW>
            rays(end+1,:) = ratQ.chk(Vi(E(j,2),:) - Vi(E(j,1),:), 'ray direction'); %#ok<AGROW>
        end
    end
    vs = unique(vs);
    bounded = isempty(rays);
end
