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
%   * Case C, one CONCAVE or AFFINE quadratic on a BOUNDED convex polygon -> the max of the
%     affine functions <s,v> - q(v) over the polygon's vertices, since <s,x> - q(x) is then
%     CONVEX in x and a convex function attains its maximum over a polytope at an extreme point.
%     One cell per vertex that actually wins somewhere; the rest are empty and are dropped.
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

    % ---- one quadratic on a bounded convex polygon: classify H EXACTLY and route ---------------
    if nf == 1 && obj.nv >= 3 && all(obj.E(:,3) == 1)
        Hn = [fN(5) fN(6); fN(6) fN(7)];
        D  = ratQ.detExact(Hn);
        if Hn(1,1) > 0 && D > 0
            g = caseBConvexOnPolygon(obj);          % strictly convex: vertex/edge/interior cells
            return
        elseif Hn(1,1) <= 0 && Hn(2,2) <= 0 && D >= 0
            g = caseCConcaveOnPolygon(obj);         % concave or affine: the max is at a vertex
            return
        end
        error('PLQ:conjQ:notImplemented', ...
            ['on a bounded polygon the exact conjugate covers a strictly convex Q (leading minor ' ...
             'and determinant both positive) and a concave or affine Q (both diagonal entries ' ...
             'nonpositive and determinant nonnegative). This Q has leading minor %d, second ' ...
             'diagonal %d and determinant %d, so it is semidefinite-singular or INDEFINITE -- ' ...
             'the indefinite case needs the x*y frame change and [COAP] A.2-A.5, which is the ' ...
             'next work item.'], Hn(1,1), Hn(2,2), D);
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

function g = caseBConvexOnPolygon(obj)
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

    [fN, fD] = obj.faceQ(1);
    Hn = [fN(5) fN(6); fN(6) fN(7)];
    Ln = [fN(8); fN(9)];
    kn = fN(10);
    D  = ratQ.detExact(Hn);
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

    % ---- the polygon, exactly, as integer vertices over one denominator --------------------
    [Vi, vd, cyc] = polygonExactly(obj);
    m = numel(cyc);

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

    g = assembleQuaCon(cells);
end

% --------------------------------------------------------------------------------------------

function s = sgnOf(row, want)
% objective: the H-form sign for a constraint written as `row <= 0` (want = -1) or `row >= 0`
%            (want = +1), AFTER ratQ.conic has normalised the row's overall sign.
%
% ratQ.conic makes the first nonzero entry positive, which may have negated the whole row. The
% face's side must follow that negation, or the cell silently becomes its own complement -- which
% is a wrong answer that still produces a plausible mesh.
    r  = [0 0 0, row(1), row(2), row(3)];
    nz = find(r ~= 0, 1);
    if isempty(nz)
        error('ratQ:zeroConic', 'a constraint with all-zero coefficients names no half-plane.');
    end
    if r(nz) < 0, s = -want; else, s = want; end
end

function [Vi, vd, cyc] = polygonExactly(obj)
% objective: the domain's vertices as integers over ONE denominator, plus the cycle order.
% [output] Vi : nv x 2 integer; vd : positive integer; cyc : 1 x m vertex indices, in order
%
% The cycle is walked from E rather than taken from P, so this does not depend on P's clockwise/
% counter-clockwise convention -- a convention it would be easy to read backwards, and reading it
% backwards flips every outward normal.
    [Vi, vd] = ratQ.combineDen(obj.VN, obj.VD);
    E = obj.E;
    if any(E(:,3) == 0)
        error('PLQ:conjQ:unbounded', ...
            'Case B needs a BOUNDED polygon; edge %d is a ray.', find(E(:,3) == 0, 1));
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
            'the domain boundary is not a single simple cycle of %d edges.', m);
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

function g = caseCConcaveOnPolygon(obj)
% objective: the conjugate of a CONCAVE or AFFINE quadratic over one bounded convex polygon.
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
    [fN, fD] = obj.faceQ(1);
    Hn = [fN(5) fN(6); fN(6) fN(7)];
    Ln = [fN(8); fN(9)];
    kn = fN(10);

    [Vi, vd, cyc] = polygonExactly(obj);
    m = numel(cyc);

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

    g = assembleQuaCon(cells);
end

function g = assembleQuaCon(cells)
% objective: turn the per-cell constraint lists into a QuaCon -- deduplicating the bounding lines
%            into one canonical edge list, naming the corners, and reducing each face function.
%
% DEDUPLICATION IS BITWISE, which is the whole point of ratQ.conic: two cells that share a facet
% arrived at it by different routes, and in doubles those two spellings differ by an ULP and the
% shared facet becomes invisible (DECISIONS.md 2026-08-17, measured: 57 cells carrying 10 distinct
% functions, 4 merges out of 612 attempts). Here they are the same integer row, and `find` on an
% integer matrix is the whole comparison.
    % ---- drop the cells that describe the empty set, or no 2-D face ------------------------
    % A polygon vertex that never attains the maximum contributes a cell that is EMPTY, not small,
    % and a cell degenerating to a point or a segment carries no face either. Both are decided
    % exactly by ratQ.feasible2 rather than inferred from a sample. Leaving them in would not make
    % `eval` wrong -- no point satisfies them -- but it would put faces into the mesh that bound
    % nothing, and `nf` would stop meaning what it says.
    live = true(1, numel(cells));
    for k = 1:numel(cells)
        rows = cells(k).con;
        live(k) = ratQ.feasible2(rows(:,7) .* rows(:,4:6), true);
    end
    cells = cells(live);
    if isempty(cells)
        error('PLQ:conjQ:noCells', ...
            'every candidate cell is empty, which cannot happen for a bounded domain.');
    end

    EcQ = zeros(0,6);
    FC  = cell(numel(cells), 1);
    for k = 1:numel(cells)
        rows = cells(k).con;
        idx  = zeros(size(rows,1), 1);
        for r = 1:size(rows,1)
            c = rows(r, 1:6);
            hit = find(all(EcQ == c, 2), 1);
            if isempty(hit)
                EcQ(end+1, :) = c; %#ok<AGROW>
                hit = size(EcQ, 1);
            end
            idx(r) = hit;
        end
        FC{k} = [idx, rows(:,7)];
    end

    fN = zeros(numel(cells), 10);  fD = ones(numel(cells), 1);
    for k = 1:numel(cells)
        [fN(k,:), fD(k)] = ratQ.canon(cells(k).num, cells(k).den);
    end

    % ---- name the corners, and ONLY the corners ---------------------------------------------
    % A pair of bounding lines names a vertex when their intersection actually satisfies the rest
    % of that cell's constraints. Two distinct lines meet in at most one point, so the root index
    % is always 1. Every step here is exact: the intersection comes from ratQ.solve2 (Cramer over
    % an exact cofactor determinant) and the membership test is the sign of an integer.
    %
    % Naming every pair of non-parallel lines instead would be much simpler and WRONG in a way
    % that matters: it would list points that exist but bound nothing, so `nv` would stop meaning
    % "corners of this mesh" and any later incidence computation would start from a padded list.
    seen = false(size(EcQ,1));
    Vname = zeros(0,3);
    for k = 1:numel(cells)
        rows = FC{k};
        for p1 = 1:size(rows,1)
            for p2 = p1+1:size(rows,1)
                i = rows(p1,1);  j = rows(p2,1);
                if i == j, continue, end
                A = [EcQ(i,4) EcQ(i,5); EcQ(j,4) EcQ(j,5)];
                b = [-EcQ(i,6); -EcQ(j,6)];
                if ratQ.detExact(A) == 0, continue, end       % parallel: they meet nowhere
                [xn, xd] = ratQ.solve2(A, b);
                if ~cellHoldsAt(EcQ, rows, xn, xd), continue, end
                lo = min(i,j);  hi = max(i,j);
                if ~seen(lo,hi)
                    seen(lo,hi) = true;
                    Vname(end+1, :) = [lo hi 1]; %#ok<AGROW>
                end
            end
        end
    end

    % E and F are left EMPTY on purpose -- see caseBConvexOnPolygon's header. A zero E of the right
    % shape keeps the constructor's arity checks meaningful without claiming an incidence that has
    % not been computed.
    g = QuaCon(Vname, zeros(size(EcQ,1),3), EcQ, fN, fD, zeros(size(EcQ,1),2), FC);
end

function tf = cellHoldsAt(EcQ, rows, xn, xd)
% objective: does the point xn/xd satisfy EVERY constraint of this cell. Exact.
% xd > 0 (ratQ.canon normalises it), so multiplying the constraint through by xd preserves signs.
    tf = true;
    for r = 1:size(rows,1)
        c = EcQ(rows(r,1), :);
        val = ratQ.chk(c(4)*xn(1) + c(5)*xn(2) + c(6)*xd, 'constraint at corner');
        if rows(r,2) * val < 0, tf = false; return, end
    end
end
