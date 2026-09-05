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
%
% THE WHOLE CASE ANALYSIS, and it is two questions asked in order. FIRST, can the sup be finite at
% all: on an unbounded piece that is a condition on the recession CONE, and recessionConditions
% decides it exactly and returns whatever linear conditions on s it imposes. SECOND, is the
% objective strictly concave -- i.e. is H positive DEFINITE -- because that is what decides whether
% the maximiser can be interior (the KKT active set, caseB) or must lie on the boundary (the max
% over edges, caseD). Both questions are integer comparisons; conjCPLQ asks the second with
% eig against a fixed threshold, which is what silently claims an empty domain on a badly scaled
% input (DECISIONS.md 2026-09-03).
    [fN, fD] = obj.faceQ(k);
    Hn = [fN(5) fN(6); fN(6) fN(7)];
    D  = ratQ.detExact(Hn);

    if obj.nv == 0
        g = caseAFullDomain(fN, fD);               % the whole plane: no domain constraints at all
        return
    end

    [Vi, vd] = ratQ.combineDen(obj.VN, obj.VD);

    % ---- a domain of dimension < 2: a needle, a segment, a ray, or a chain of edges ------------
    % THE CONJUGATE OF A LOW-DIMENSIONAL INPUT IS FULL-DIMENSIONAL, which is why this needs no new
    % storage: a needle at p with value c conjugates to <s,p> - c, affine and finite on all of R^2,
    % and a segment to a max over its two endpoints and, where the restriction is concave, its
    % interior stationary point. That is exactly caseD's candidate set, so the only new code is
    % building the shape -- there is no face to read it from, since such a mesh has nf = 0 and an
    % empty F, so every edge belongs to the single piece.
    %
    % caseD unconditionally, never caseB, and the reason is that the problem is ONE-dimensional:
    % what decides the shape is the curvature ALONG the edge, d'Hd, not the definiteness of H in
    % the plane. An interior maximiser of the 1-D restriction is already the edge candidate.
    if obj.dom.dim < 2
        if size(obj.fN,1) > 1
            error('PLQ:conjQ:ambiguousChain', ...
                ['a domain of dimension < 2 carrying %d different functions has no F to say ' ...
                 'which edge takes which, so the mapping is ambiguous. One function only.'], ...
                size(obj.fN,1));
        end
        sh = degenerateShape(obj, Vi);
        [okFinite, dom, ~] = recessionConditions(sh, Hn, [fN(8); fN(9)], fD, Vi, vd);
        if ~okFinite, g = emptyConjugate(); return, end
        g = caseDBoundaryMax(fN, fD, Vi, vd, sh, dom);
        return
    end

    sh = pieceShape(obj, k, Vi, vd);
    if numel(sh.vs) == 0
        error('PLQ:conjQ:noVertex', ...
            ['piece %d has no vertex, so it is a half-plane or a slab and its conjugate is ' ...
             'supported on a line or a point -- a domain of dimension < 2, which this type ' ...
             'cannot carry.'], k);
    end

    % ---- a NON-CONVEX piece is SPLIT, not refused ----------------------------------------------
    % conj turns a union into a MAX: f restricted to P = union of P_i gives
    % f* = max_i (q + I_{P_i})*, so triangulating the face and folding the pieces with maxQ is
    % exact, not an approximation. The triangulation itself is exact too -- ear clipping, whose
    % every decision is the sign of an integer cross product.
    if ~sh.convex
        if ~sh.bounded
            error('PLQ:conjQ:nonConvexPiece', ...
                ['face %d is non-convex AND unbounded. Ear clipping needs a closed boundary, so ' ...
                 'splitting it means cutting along its recession directions first -- a real ' ...
                 'case, and a separate one.'], k);
        end
        tris = earClip(Vi, boundaryCycle(sh, Vi));
        g = [];
        for t = 1:size(tris,1)
            sht = triangleShape(Vi, vd, tris(t,:));
            if Hn(1,1) > 0 && D > 0
                gt = caseBConvexOnPiece(fN, fD, Vi, vd, sht, zeros(0,3));
            else
                gt = caseDBoundaryMax(fN, fD, Vi, vd, sht, zeros(0,3));
            end
            if isempty(g), g = gt; else, g = maxQ(g, gt); end
        end
        return
    end

    [okFinite, dom, why] = recessionConditions(sh, Hn, [fN(8); fN(9)], fD, Vi, vd);
    if ~okFinite
        % dom f* is EMPTY -- `why` says which recession direction does it. That is the right
        % answer and it IS representable: a QuaCon with zero faces evaluates to +infinity
        % everywhere. This used to raise PLQ:conjQ:emptyDomain on the belief that the type could
        % not hold such a function, which was wrong.
        assert(~isempty(why));                 %#ok<NASGU> -- the reason, kept for debugging
        g = emptyConjugate();
        return
    end

    if Hn(1,1) > 0 && D > 0
        g = caseBConvexOnPiece(fN, fD, Vi, vd, sh, dom);     % strictly convex: KKT active set
    else
        g = caseDBoundaryMax(fN, fD, Vi, vd, sh, dom);       % everything else: max over the boundary
    end
end

function g = emptyConjugate()
% objective: the function that is +infinity everywhere, as a QuaCon.
%
% A QuaCon with NO faces: `eval` matches nothing and so answers +inf at every point, which is
% exactly that function. Verified directly rather than assumed. This is the honest answer whenever
% the sup diverges for every s, and it needs no extension to the type -- what was missing was only
% the willingness to return it.
    g = QuaCon(zeros(0,3), zeros(0,3), zeros(0,6), zeros(0,10), zeros(0,1), zeros(0,2), {});
end

function sh = degenerateShape(obj, Vi)
% objective: the shape record for a domain of dimension < 2, where every edge belongs to the one
%            piece and there are no bounding half-planes.
%
% `hp` is empty on purpose: a segment is not the intersection of half-planes, it IS the candidate
% set. The membership constraint that a 2-D piece needs is carried here by the edge's own clamping
% conditions instead.
    sh.vs = [];
    sh.ed = struct('a', {}, 'b', {}, 'd', {}, 'isRay', {});
    sh.rays = zeros(0,2);
    sh.hp = zeros(0,3);
    sh.bounded = true;
    for j = 1:size(obj.E,1)
        a = obj.E(j,1);  b = obj.E(j,2);
        d = ratQ.chk(Vi(b,:) - Vi(a,:), 'edge direction');
        if obj.E(j,3) ~= 0
            sh.vs = [sh.vs; a; b];
            sh.ed(end+1) = struct('a', a, 'b', b, 'd', d, 'isRay', false); %#ok<AGROW>
        else
            sh.vs = [sh.vs; a];
            sh.ed(end+1) = struct('a', a, 'b', 0, 'd', d, 'isRay', true); %#ok<AGROW>
            sh.rays(end+1,:) = d; %#ok<AGROW>
            sh.bounded = false;
        end
    end
    if isempty(sh.ed)
        sh.vs = (1:size(obj.VN,1)).';        % a NEEDLE: the single vertex, and no edge at all
    end
    sh.vs = unique(sh.vs);
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
        % NOT ONE CASE BUT TWO, and only the second is still a gap.
        %
        % If H has a direction of NEGATIVE curvature -- H not positive semidefinite, which covers
        % negative definite, negative semidefinite-singular and indefinite -- then <s,x> - q(x)
        % rises without bound along it for EVERY s, so f* is +infinity everywhere. That is an
        % answer, and a QuaCon with zero faces is exactly it.
        if ~ratQ.isPSD2(Hn)
            g = emptyConjugate();
            return
        end
        % Otherwise H is positive SEMIdefinite and singular, and dom f* is thin rather than empty:
        % a LINE when H is nonzero (the sup is finite only for s in the range of H), and a single
        % POINT when H is zero (q affine, so f* is finite only at s = L). Both are correct answers
        % and neither fits a QuaCon, whose faces are two-dimensional. Representing them needs the
        % H-form's side column to carry 0 for "on the curve" -- one value in an existing field,
        % not a new type. See TODO.md 2026-09-04.
        % A THIN dual domain, and both kinds are now BUILT rather than refused, using the
        % H-form's equality side (a side of 0 on a curve means ON it).
        if all(Hn(:) == 0)
            % q is AFFINE on the plane: f*(s) = -kappa at the single point s = L, +inf elsewhere.
            % Two equalities, fD*s_i - Ln_i = 0, and a constant face.
            con = [ratQ.conic([0 0 0, fD, 0, -Ln(1)]), sgnOf([fD, 0, -Ln(1)], 0); ...
                   ratQ.conic([0 0 0, 0, fD, -Ln(2)]), sgnOf([0, fD, -Ln(2)], 0)];
            g = assembleQuaConCells(struct('num', [0 0 0 0 0 0 0 0 0 -fN(10)], ...
                                           'den', fD, 'con', con));
            return
        end
        % q is positive semidefinite and SINGULAR. With u = fD*s - Ln, the unconstrained sup of
        % u'x - x'Hn x/2 is finite exactly when u lies in range(Hn) -- equivalently u is orthogonal
        % to the null direction n -- and then equals u'pinv(Hn)u/2. For a rank-one PSD Hn the
        % pseudo-inverse is m*m'/(m'Hn m) for any m spanning the range, so
        %       f*(s) = ( (u.m)^2 / (2 m'Hn m) - kappa ) / fD      on the line  u.n = 0,
        % which is rational throughout. Checked against q = x^2/2: n = (0,1) gives s2 = 0, m = (1,0)
        % gives f* = s1^2/2.
        if Hn(1,1) ~= 0
            nl = [-Hn(1,2); Hn(1,1)];               % Hn*nl = 0 because det(Hn) = 0
        else
            nl = [1; 0];                            % Hn(1,1) = 0 forces Hn(1,2) = 0 here
        end
        m = [Hn(1,1); Hn(1,2)];
        if all(m == 0), m = [Hn(1,2); Hn(2,2)]; end
        mHm = ratQ.chk(m.' * Hn * m, 'range scaling');
        P = [ratQ.chk(fD*m(1),'p1'), ratQ.chk(fD*m(2),'p2'), ratQ.chk(-(Ln.'*m),'p3')];
        den = ratQ.chk(2 * mHm * fD, 'line denominator');
        num = [0 0 0 0, ratQ.chk(2*P(1)^2,'c5'), ratQ.chk(2*P(1)*P(2),'c6'), ...
               ratQ.chk(2*P(2)^2,'c7'), ratQ.chk(2*P(1)*P(3),'c8'), ratQ.chk(2*P(2)*P(3),'c9'), ...
               ratQ.chk(P(3)^2 - 2*mHm*fN(10), 'c10')];
        row = [ratQ.chk(fD*nl(1),'e1'), ratQ.chk(fD*nl(2),'e2'), ratQ.chk(-(Ln.'*nl),'e0')];
        g = assembleQuaConCells(struct('num', num, 'den', den, ...
                                       'con', [ratQ.conic([0 0 0, row]), sgnOf(row, 0)]));
        return
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

function g = caseBConvexOnPiece(fN, fD, Vi, vd, sh, dom)
% objective: the conjugate of a STRICTLY CONVEX quadratic over one polyhedral piece, bounded or
%            not, exactly, as the KKT active-set subdivision.
%
% GENERALISED 2026-09-04 from the bounded-polygon version. The decomposition is unchanged and is
% conjConvexOverPiece's: the maximiser of the concave program lies on exactly one relatively open
% face of P, so there is one AFFINE cell per vertex (its normal cone, shifted by grad q(v)), one
% QUADRATIC cell per edge of positive curvature, and one interior cell.
%
% WHAT AN UNBOUNDED PIECE CHANGES, and conjConvexOverPiece's header already said the first part:
% "a RAY leaving v contributes its direction exactly as a bounded edge does, so unbounded pieces
% need no special case" -- the normal cone is built from the directions leaving the vertex, and a
% ray gives one. The second part is the edge cell: a segment is clamped at BOTH ends (0 <= t* <= 1)
% and a ray only at its base (t* >= 0), which is one affine condition instead of two. The third is
% the interior cell, whose membership test is the piece's half-planes evaluated at the maximiser --
% which for a bounded polygon was written from its vertex cycle and is now read off sh.hp, the same
% list for either shape.
    Hn = [fN(5) fN(6); fN(6) fN(7)];
    Ln = [fN(8); fN(9)];
    kn = fN(10);
    D  = ratQ.detExact(Hn);
    if ~(Hn(1,1) > 0 && D > 0)
        error('PLQ:conjQ:notStrictlyConvex', ...
            'caseB needs a strictly convex Q; leading minor %d, determinant %d.', Hn(1,1), D);
    end

    denV = ratQ.chk(2 * fD * vd^2, 'vertex denominator');
    cells = struct('num', {}, 'den', {}, 'con', {});

    % ---- one cell per VERTEX: its normal cone, shifted by grad q(v) ---------------------------
    for i = 1:numel(sh.vs)
        vi = sh.vs(i);
        v  = Vi(vi, :).';
        gv = ratQ.chk(Hn * v + vd * Ln, 'vertex gradient');       % over fD*vd
        qv = ratQ.chk(v.'*Hn*v + 2*vd*(Ln.'*v) + 2*vd^2*kn, 'vertex value');
        num = [0 0 0 0, 0, 0, 0, ...
               ratQ.chk(2*fD*vd*v(1),'c8'), ratQ.chk(2*fD*vd*v(2),'c9'), -qv];

        rows = dom;
        dirs = dirsAtVertex(sh, vi);
        for kk = 1:numel(dirs)
            ed = dirs{kk};
            % <s - grad q(v), e> <= 0, cleared by fD*vd, written as a >= 0 half-plane
            rows(end+1,:) = [ratQ.chk(-fD*vd*ed(1),'n1'), ratQ.chk(-fD*vd*ed(2),'n2'), ...
                             ratQ.chk(gv.' * ed, 'n0')]; %#ok<AGROW>
        end
        [con, ok] = toCon(rows);
        if ok
            cells(end+1) = struct('num', num, 'den', denV, 'con', con); %#ok<AGROW>
        end
    end

    % ---- one cell per EDGE of positive curvature ----------------------------------------------
    for j = 1:numel(sh.ed)
        e = sh.ed(j);
        a = Vi(e.a, :).';
        d = e.d(:);
        alpha = ratQ.chk(d.' * Hn * d, 'edge curvature');
        if alpha <= 0
            % Impossible for a strictly convex Q and a nonzero direction; kept as an assertion.
            error('PLQ:conjQ:degenerateEdge', ...
                'edge %d has d''Qd = %d <= 0 under a strictly convex Q.', j, alpha);
        end
        ga = ratQ.chk(Hn * a + vd * Ln, 'edge base gradient');
        qa = ratQ.chk(a.'*Hn*a + 2*vd*(Ln.'*a) + 2*vd^2*kn, 'edge base value');
        Td = [ratQ.chk(fD*vd*d(1),'t1'), ratQ.chk(fD*vd*d(2),'t2'), ratQ.chk(-(ga.'*d),'t0')];

        den = ratQ.chk(2 * alpha * fD * vd^2, 'edge denominator');
        lin = [ratQ.chk(2*alpha*fD*vd*a(1),'c8'), ratQ.chk(2*alpha*fD*vd*a(2),'c9')];
        num = [0 0 0 0, ...
               ratQ.chk(2*Td(1)^2,'c5'), ratQ.chk(2*Td(1)*Td(2),'c6'), ratQ.chk(2*Td(2)^2,'c7'), ...
               ratQ.chk(lin(1) + 2*Td(1)*Td(3),'c8'), ratQ.chk(lin(2) + 2*Td(2)*Td(3),'c9'), ...
               ratQ.chk(Td(3)^2 - alpha*qa,'c10')];

        rows = [dom; Td];                                  % 0 <= t*
        if ~e.isRay
            rows(end+1,:) = [-Td(1), -Td(2), ratQ.chk(alpha - Td(3), 'edge upper')]; %#ok<AGROW>
        end
        % the multiplier condition, cleared by alpha > 0
        n = outwardNormalOf(sh, j, Vi);
        Hdn = ratQ.chk(d.' * Hn * n, 'edge multiplier coupling');
        rows(end+1,:) = [ratQ.chk(alpha*fD*vd*n(1) - Td(1)*Hdn, 'm1'), ...
                         ratQ.chk(alpha*fD*vd*n(2) - Td(2)*Hdn, 'm2'), ...
                         ratQ.chk(-alpha*(ga.'*n) - Td(3)*Hdn, 'm0')]; %#ok<AGROW>
        [con, ok] = toCon(rows);
        if ok
            cells(end+1) = struct('num', num, 'den', den, 'con', con); %#ok<AGROW>
        end
    end

    % ---- the INTERIOR cell ----------------------------------------------------------------------
    % x* = H^-1(s - L) = fD*A*(s - L/fD)/D with A = adj(Hn); it is in the piece exactly when every
    % half-plane holds there, and each of those is affine in s.
    A = [Hn(2,2), -Hn(1,2); -Hn(1,2), Hn(1,1)];
    ALn = A * Ln;
    numI = [0 0 0 0, ratQ.chk(2*fD^2*A(1,1),'c5'), ratQ.chk(2*fD^2*A(1,2),'c6'), ...
            ratQ.chk(2*fD^2*A(2,2),'c7'), ratQ.chk(-2*fD*ALn(1),'c8'), ...
            ratQ.chk(-2*fD*ALn(2),'c9'), ratQ.chk(Ln.'*ALn - 2*D*kn, 'c10')];
    rows = dom;
    for r = 1:size(sh.hp,1)
        % hp is [p q c] meaning p*x + q*y + c >= 0 at actual points; substitute x = fD*A*s/D - A*L/D
        p = sh.hp(r,1:2);
        Ap = A.' * p.';
        rows(end+1,:) = [ratQ.chk(fD*Ap(1), 'i1'), ratQ.chk(fD*Ap(2), 'i2'), ...
                         ratQ.chk(-(Ap.'*Ln) + D*sh.hp(r,3), 'i0')]; %#ok<AGROW>
    end
    [con, ok] = toCon(rows);
    if ok
        cells(end+1) = struct('num', numI, 'den', ratQ.chk(2*fD*D,'interior denominator'), ...
                              'con', con);
    end

    g = assembleQuaConCells(cells);
end

function n = outwardNormalOf(sh, j, Vi)
% objective: the OUTWARD normal of edge j, as an integer column.
%
% Built from the edge's own direction and then ORIENTED, rather than matched against the
% half-plane list. Matching by perpendicularity is AMBIGUOUS whenever two edges are parallel: on
% the unit square the top and bottom edges have the same normal direction, so the scan returned the
% first one and the top edge's cell came out as s2 <= 1 where it must be s2 >= 1. Measured: 74 of
% 307 dual points wrong and 53 in no cell at all.
%
% Outward means every point and every recession direction of the piece lies on the NON-positive
% side: <n, v - a> <= 0 for each vertex v, and <n, r> <= 0 for each recession direction r. That
% fixes the sign with no reference to any other edge, so parallel edges cannot be confused.
    d = sh.ed(j).d(:);
    n = [d(2); -d(1)];
    a = Vi(sh.ed(j).a, :).';
    for i = 1:numel(sh.vs)
        t = ratQ.chk(n.' * (Vi(sh.vs(i), :).' - a), 'orientation');
        if t > 0, n = -n; return, end
        if t < 0, return, end
    end
    for r = 1:size(sh.rays,1)
        t = ratQ.chk(n.' * sh.rays(r,:).', 'orientation');
        if t > 0, n = -n; return, end
        if t < 0, return, end
    end
end

function [con, ok] = toCon(rows)
% objective: turn [a b c] half-planes into the H-form [conic6, side] rows, dropping the ones that
%            say nothing and reporting the ones that make the cell EMPTY.
%
% A row whose linear part vanishes is not a constraint on s at all: it reads c >= 0, which is
% either vacuous (c >= 0) or unsatisfiable (c < 0). Dropping BOTH would be wrong -- the second
% means the cell is empty -- so the two are separated here rather than silently merged.
    ok = true;  con = zeros(0,7);
    if isempty(rows), return, end
    for r = 1:size(rows,1)
        if all(rows(r,1:2) == 0)
            if rows(r,3) < 0, ok = false; return, end
            continue                                   % vacuous
        end
        con(end+1,:) = [ratQ.conic([0 0 0, rows(r,:)]), sgnOf(rows(r,:), +1)]; %#ok<AGROW>
    end
end

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
    if want == 0
        % An EQUALITY is symmetric: {c = 0} and {-c = 0} are the same set, so ratQ.conic's sign
        % normalisation cannot disturb it and there is nothing to follow.
        s = 0;  return
    end
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

function sh = pieceShape(obj, k, Vi, vd)
% objective: everything the conjugate needs to know about ONE piece's geometry, bounded or not.
% [output] sh : struct with
%     vs    : the piece's true vertex indices (a ray's direction marker is NOT one)
%     ed    : struct array of edges, each with .a (base vertex index), .d (integer direction),
%             .b (the other vertex index, 0 for a ray) and .isRay
%     rays  : r x 2 integer extreme recession directions
%     hp    : h x 3 half-planes [a b c] meaning a*x + b*y + c >= 0 on the piece, in ACTUAL
%             coordinates
%     bounded : logical
%
% ONE ROUTINE FOR BOTH SHAPES, replacing the cycle walk. A bounded face was described by walking its
% boundary in order, which an unbounded face has no way to give -- its boundary is not a cycle. The
% edge list plus the half-planes is what both shapes DO have, and nothing downstream actually needed
% the ordering: the vertex normal cones are built from the edges incident to each vertex and the
% cells from pairwise comparisons.
%
% THE SIDE OF EACH EDGE is fixed by requiring every vertex AND every recession direction of the
% piece to satisfy it, so no orientation convention is consulted and none can be read backwards. A
% recession direction imposes a condition on the LINEAR PART only: the constraint must survive
% travelling to infinity along it, which is <n, d> >= 0.
%
% TWO COORDINATE SYSTEMS, NOT TO BE MIXED. Vi holds NUMERATORS over vd, so one line reads
% n.X + c0 = 0 at numerators and vd*n.x + c0 = 0 at actual points. The side is decided at
% numerators, where the vertices live; the row is emitted for actual points, where eval reads it.
% Using one where the other belongs rescales the offset into a PARALLEL line -- the defect measured
% in biconjQ on a triangle with half-integer vertices, which excluded two of its own three vertices.
    own = find(any(obj.F == k, 2));
    if isempty(own)
        error('PLQ:conjQ:noFace', 'face %d has no edges.', k);
    end
    E = obj.E(own, :);

    sh.vs = [];
    sh.ed = struct('a', {}, 'b', {}, 'd', {}, 'isRay', {});
    sh.rays = zeros(0,2);
    for j = 1:size(E,1)
        a = E(j,1);  b = E(j,2);
        d = ratQ.chk(Vi(b,:) - Vi(a,:), 'edge direction');
        if E(j,3) ~= 0
            sh.vs = [sh.vs; a; b];
            sh.ed(end+1) = struct('a', a, 'b', b, 'd', d, 'isRay', false); %#ok<AGROW>
        else
            sh.vs = [sh.vs; a];
            sh.ed(end+1) = struct('a', a, 'b', 0, 'd', d, 'isRay', true); %#ok<AGROW>
            sh.rays(end+1,:) = d; %#ok<AGROW>
        end
    end
    sh.vs = unique(sh.vs);
    sh.bounded = isempty(sh.rays);
    sh.convex = true;                 % until an edge below finds a vertex on its far side

    % ---- the half-planes -----------------------------------------------------------------------
    sh.hp = zeros(0,3);
    for j = 1:numel(sh.ed)
        p = Vi(sh.ed(j).a, :);
        n = [sh.ed(j).d(2), -sh.ed(j).d(1)];          % a normal to this edge
        c0 = ratQ.chk(-(n * p.'), 'offset');
        sgn = 0;
        for i = 1:numel(sh.vs)
            v = Vi(sh.vs(i), :);
            t = ratQ.chk(n(1)*v(1) + n(2)*v(2) + c0, 'vertex side');
            if t ~= 0, sgn = sign(t); break, end
        end
        if sgn == 0
            for r = 1:size(sh.rays,1)
                t = ratQ.chk(n * sh.rays(r,:).', 'recession side');
                if t ~= 0, sgn = sign(t); break, end
            end
        end
        if sgn == 0, continue, end                     % everything on the line: degenerate

        % IS THE PIECE CONVEX? A face described by half-planes IS an intersection of half-planes,
        % so a REFLEX corner cannot be expressed: every other vertex would have to lie on the
        % chosen side of this edge, and at a reflex corner some do not. Recorded rather than
        % raised, because a non-convex face is now SPLIT rather than refused -- but it must never
        % be conjugated through this description, which would silently take the sup over the
        % smaller set the half-planes cut out (measured on an L: 54 of 303 dual points too low).
        for i = 1:numel(sh.vs)
            v = Vi(sh.vs(i), :);
            if ratQ.chk(sgn * (n(1)*v(1) + n(2)*v(2) + c0), 'convexity') < 0
                sh.convex = false;
            end
        end
        for r = 1:size(sh.rays,1)
            if ratQ.chk(sgn * (n * sh.rays(r,:).'), 'convexity') < 0
                sh.convex = false;
            end
        end

        sh.hp(end+1,:) = sgn * [ratQ.chk(vd*n(1), 'edge normal'), ...
                                ratQ.chk(vd*n(2), 'edge normal'), c0]; %#ok<AGROW>
    end
end

function dirs = dirsAtVertex(sh, i)
% objective: the edge directions leaving vertex index i, as integer column vectors.
% A ray contributes its direction at its BASE only -- its far end is a marker, not a vertex.
    dirs = {};
    for j = 1:numel(sh.ed)
        e = sh.ed(j);
        if e.a == i
            dirs{end+1} = e.d(:); %#ok<AGROW>
        elseif ~e.isRay && e.b == i
            dirs{end+1} = -e.d(:); %#ok<AGROW>
        end
    end
end

function [ok, dom, why] = recessionConditions(sh, Hn, Ln, fD, Vi, vd)
% objective: decide whether the sup can be finite on this piece, and return the LINEAR conditions
%            on s that make it so.
% [output] ok : false when the sup is +infinity for every s, i.e. dom f* is empty
%          dom : d x 3 half-planes on s;  why : the reason when ok is false
%
% ALONG a recession direction d the objective carries  t*slope(x)  -  (t^2/2) d'Hd, so
%   d'Hd > 0   the objective falls away: no condition at all;
%   d'Hd < 0   it rises without bound for every s: dom f* is EMPTY;
%   d'Hd = 0   it is linear in t and finite exactly when the slope is <= 0 -- for EVERY x of the
%              piece, since the ray may be based anywhere on it.
%
% THE NULL CASE IS STILL A HALF-PLANE, and working that out removed a refusal rather than adding
% one. The slope is
%       slope(x) = <s - L, d> - <Hd, x>,
% so the condition is  <s - L, d>  <=  inf over the piece of <Hd, x>  -- and that infimum is a
% LINEAR PROGRAM over a polyhedron, hence: minus infinity when some recession direction r has
% <Hd, r> < 0 (and then no s works at all, so dom f* is empty), and otherwise the minimum over the
% VERTICES. Both are exact integer computations. When Hd = 0 the infimum is 0 and this reduces to
% <s - L, d> <= 0, which is the concave case's own condition.
    ok = true;  dom = zeros(0,3);  why = '';
    R = sh.rays;
    if isempty(R), return, end

    for r = 1:size(R,1)
        if ratQ.chk(R(r,:) * Hn * R(r,:).', 'ray curvature') < 0
            ok = false;
            why = sprintf('the piece recedes along (%d,%d), where the objective rises without bound for every s', ...
                          R(r,1), R(r,2));
            return
        end
    end
    % THE CONE, NOT ONLY ITS GENERATORS. On cone(d1,d2) the form is a*l1^2 + 2b*l1*l2 + c*l2^2, and
    % it is nonnegative for every l >= 0 exactly when a >= 0, c >= 0 and (b >= 0 or b^2 <= a*c).
    % Testing the generators alone would miss a negative direction strictly inside the cone.
    if size(R,1) >= 2
        a = ratQ.chk(R(1,:)*Hn*R(1,:).', 'cone a');
        b = ratQ.chk(R(1,:)*Hn*R(2,:).', 'cone b');
        c = ratQ.chk(R(2,:)*Hn*R(2,:).', 'cone c');
        if b < 0 && ratQ.chk(b^2 - a*c, 'cone discriminant') > 0
            ok = false;
            why = 'some direction strictly inside the recession cone has negative curvature';
            return
        end
    end

    for r = 1:size(R,1)
        d = R(r,:).';
        if ratQ.chk(d.' * Hn * d, 'ray curvature') ~= 0, continue, end
        Hd = ratQ.chk(Hn * d, 'null direction image');

        % inf over the piece of <Hd, x>: -infinity if it decreases along some recession direction
        for q = 1:size(R,1)
            if ratQ.chk(R(q,:) * Hd, 'slope along recession') < 0
                ok = false;
                why = sprintf(['the objective''s slope along the null recession direction ' ...
                               '(%d,%d) is unbounded above on the piece'], d(1), d(2));
                return
            end
        end
        mu = inf;
        for i = 1:numel(sh.vs)
            t = ratQ.chk(Vi(sh.vs(i), :) * Hd, 'slope at vertex');
            if t < mu, mu = t; end          % over vd, like Vi itself
        end

        % <s,d> - (Ln.d)/fD <= mu/vd, cleared by fD*vd > 0 and written as a >= 0 half-plane
        row = [ratQ.chk(-fD*vd*d(1), 'r1'), ratQ.chk(-fD*vd*d(2), 'r2'), ...
               ratQ.chk(fD*mu + vd*(Ln.'*d), 'r0')];
        if all(row == 0), continue, end
        dom(end+1,:) = row; %#ok<AGROW>
    end
end

function g = caseDBoundaryMax(fN, fD, Vi, vd, sh, dom)
% objective: the conjugate of a quadratic that is NOT positive definite, over one polyhedral piece,
%            bounded or not.
%
% THE ONE FACT THAT MAKES THIS CASE TRACTABLE, and it replaces the whole A.2-A.5 apparatus: when H
% is not positive DEFINITE the sup of <s,x> - q(x) over P is attained on the BOUNDARY. An interior
% maximiser would force the objective's Hessian -H to be negative semidefinite, i.e. H PSD; and if H
% is PSD but SINGULAR there is a direction of zero curvature along which one can walk to the
% boundary without decreasing. Only a positive definite H makes the interior stationary point a
% strict maximiser, and that is caseB.
%
% So q*(s) is the largest of finitely many CANDIDATES: every vertex value <s,v> - q(v), which is
% affine; and for every edge of positive curvature d'Hd > 0, that edge's clamped one-dimensional
% maximum, which is quadratic and available only where the clamp is inactive -- an AFFINE condition,
% since t* = T(s)/alpha is affine. An edge of non-positive curvature has a convex restriction whose
% maximum sits at an endpoint, already carried by the vertices, so it contributes nothing.
%
% GENERALISED 2026-09-04 to an unbounded piece. A RAY is clamped only at its base, so its
% availability is the single condition t* >= 0 rather than 0 <= t* <= 1; and the finiteness
% conditions on s that an unbounded piece imposes arrive as `dom` and are affine too, so they simply
% join the arrangement below. When H is negative semidefinite no edge qualifies, the candidate set
% is the vertices alone, and this returns exactly the vertex max -- so the concave case is one
% branch of this rather than a separate algorithm.
%
% CLASSIFY FIRST, SPLIT LAST. Every availability condition is AFFINE, so they cut the plane into a
% LINE arrangement whose cells ratQ.feasible2 decides exactly AND completely. Build that first,
% pruning each partial assignment the moment it is infeasible -- so the 3^(#edges) enumeration never
% materialises -- and only then ask which of a cell's few surviving candidates is largest, which is
% the only place a conic ever enters. Measured when this replaced a pairwise maxQ fold: 274 faces
% to 121, and 10 s to 4 s.
    Hn = [fN(5) fN(6); fN(6) fN(7)];
    Ln = [fN(8); fN(9)];
    kn = fN(10);
    m  = numel(sh.vs);

    denV = ratQ.chk(2 * fD * vd^2, 'vertex denominator');
    qv = zeros(m,1);  numV = zeros(m,10);
    for i = 1:m
        v = Vi(sh.vs(i), :).';
        qv(i) = ratQ.chk(v.'*Hn*v + 2*vd*(Ln.'*v) + 2*vd^2*kn, 'vertex value');
        numV(i,:) = [0 0 0 0, 0, 0, 0, ...
                     ratQ.chk(2*fD*vd*v(1),'c8'), ratQ.chk(2*fD*vd*v(2),'c9'), -qv(i)];
    end

    Ed = struct('Td', {}, 'alpha', {}, 'num', {}, 'den', {}, 'isRay', {});
    for j = 1:numel(sh.ed)
        e = sh.ed(j);
        d = e.d(:);
        alpha = ratQ.chk(d.' * Hn * d, 'edge curvature');
        if alpha <= 0, continue, end
        a  = Vi(e.a, :).';
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
            'den', ratQ.chk(2*alpha*fD*vd^2, 'edge denominator'), 'isRay', e.isRay); %#ok<AGROW>
    end

    % ---- the LINE arrangement: which vertex leads, and which edges are available ---------------
    cells0 = struct('lin', {}, 'vi', {}, 'ei', {});
    for i = 1:m
        rows = dom;
        for j = 1:m
            if j == i, continue, end
            r = [ratQ.chk(2*fD*vd*(Vi(sh.vs(i),1)-Vi(sh.vs(j),1)), 'd1'), ...
                 ratQ.chk(2*fD*vd*(Vi(sh.vs(i),2)-Vi(sh.vs(j),2)), 'd2'), ...
                 ratQ.chk(-(qv(i) - qv(j)), 'd0')];
            if all(r == 0)
                error('PLQ:conjQ:repeatedVertex', ...
                    'piece vertices %d and %d coincide.', sh.vs(i), sh.vs(j));
            end
            rows(end+1,:) = r; %#ok<AGROW>
        end
        if ratQ.feasible2(rows, true)
            cells0(end+1) = struct('lin', rows, 'vi', i, 'ei', []); %#ok<AGROW>
        end
    end

    for e = 1:numel(Ed)
        Td = Ed(e).Td;  alpha = Ed(e).alpha;
        if Ed(e).isRay
            % clamped at the base only: available on t* >= 0, unavailable on t* <= 0
            regimes = { [-Td(1) -Td(2) -Td(3)], Td };
            avail   = [false, true];
        else
            regimes = { [-Td(1) -Td(2) -Td(3)], ...
                        [Td; -Td(1), -Td(2), alpha - Td(3)], ...
                        [Td(1) Td(2) Td(3)-alpha] };
            avail   = [false, true, false];
        end
        nxt = struct('lin', {}, 'vi', {}, 'ei', {});
        for c = 1:numel(cells0)
            for r = 1:numel(regimes)
                lin = [cells0(c).lin; regimes{r}];
                if ~ratQ.feasible2(lin, true), continue, end
                ei = cells0(c).ei;
                if avail(r), ei = [ei, e]; end
                nxt(end+1) = struct('lin', lin, 'vi', cells0(c).vi, 'ei', ei); %#ok<AGROW>
            end
        end
        cells0 = nxt;
    end

    % ---- inside each linear cell, the winner among its few candidates --------------------------
    out = struct('num', {}, 'den', {}, 'con', {});
    for c = 1:numel(cells0)
        [base, ok] = toCon(cells0(c).lin);
        if ~ok, continue, end

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
                if all(dn == 0), continue, end
                R = ratQ.chk([dn(5), 2*dn(6), dn(7), 2*dn(8), 2*dn(9), 2*dn(10)], 'winner');
                con(end+1,:) = [ratQ.conic(R), sgnOf(R, +1)]; %#ok<AGROW>
            end
            out(end+1) = struct('num', cand(i).num, 'den', cand(i).den, 'con', con); %#ok<AGROW>
        end
    end

    g = assembleQuaConCells(out);
end

function cyc = boundaryCycle(sh, Vi)
% objective: the face's boundary as an ORDERED cycle of vertex indices.
% [output] cyc : 1 x m, consecutive entries joined by an edge, walked from the edge list
%
% Walked from the edges rather than taken from P, so it does not depend on P's clockwise /
% counter-clockwise convention -- a convention easy to read backwards, and reading it backwards
% reverses the orientation that the ear test below depends on.
    m = numel(sh.ed);
    adj = containers.Map('KeyType', 'double', 'ValueType', 'any');
    for j = 1:m
        a = sh.ed(j).a;  b = sh.ed(j).b;
        for pair = [a b; b a].'
            u = pair(1);  v = pair(2);
            if isKey(adj, u), adj(u) = [adj(u), v]; else, adj(u) = v; end
        end
    end
    cyc = zeros(1, m);
    cyc(1) = sh.ed(1).a;  cyc(2) = sh.ed(1).b;
    for t = 3:m
        nxt = adj(cyc(t-1));
        nxt = nxt(nxt ~= cyc(t-2));
        if isempty(nxt)
            error('PLQ:conjQ:notSimple', 'the face boundary is not a single closed cycle.');
        end
        cyc(t) = nxt(1);
    end
    if numel(unique(cyc)) ~= m
        error('PLQ:conjQ:notSimple', ...
            'the face boundary revisits a vertex, so it is not a simple polygon.');
    end
    % Orient counter-clockwise, exactly: twice the signed area is an integer sum of cross products.
    if twiceArea(Vi, cyc) < 0, cyc = fliplr(cyc); end
end

function A2 = twiceArea(Vi, cyc)
% twice the SIGNED area of the polygon, as an exact integer (over vd^2, but only its sign is used)
    A2 = 0;
    m = numel(cyc);
    for i = 1:m
        p = Vi(cyc(i), :);  q = Vi(cyc(mod(i, m) + 1), :);
        A2 = A2 + ratQ.chk(p(1)*q(2) - p(2)*q(1), 'signed area');
    end
    A2 = ratQ.chk(A2, 'signed area');
end

function tris = earClip(Vi, cyc)
% objective: triangulate a simple polygon EXACTLY, by ear clipping.
% [input]  cyc : the boundary in counter-clockwise order
% [output] tris : t x 3 vertex indices, one row per triangle
%
% WHY EAR CLIPPING AND NOT A LIBRARY. Every decision is a sign test on an integer: whether a corner
% turns left (a cross product) and whether another vertex lies inside a candidate triangle (three
% cross products). So the triangulation is DECIDED rather than estimated, which is the whole
% premise here -- a floating-point triangulator could return a different mesh for the same polygon
% depending on rounding, and each triangle's conjugate would then be exact about the wrong thing.
%
% Two theorems make the loop terminate: every simple polygon with more than three vertices has at
% least two ears (Meisters), and clipping one leaves a simple polygon. So the list shrinks by one
% each pass and the recursion is m - 2 triangles.
%
% COLLINEAR VERTICES ARE DROPPED FIRST. A vertex whose corner turns neither left nor right lies on
% the segment between its neighbours and changes nothing about the region; keeping it would make a
% degenerate ear of zero area, which is not a face and would divide by zero downstream.
    cyc = dropCollinear(Vi, cyc);
    m = numel(cyc);
    if m < 3
        error('PLQ:conjQ:degenerate', 'the face collapses to fewer than three corners.');
    end
    tris = zeros(0,3);
    guard = 0;
    while numel(cyc) > 3
        guard = guard + 1;
        if guard > m + 2
            error('PLQ:conjQ:noEar', ...
                ['no ear was found on a polygon of %d corners, which cannot happen for a simple ' ...
                 'polygon -- so the boundary is self-intersecting.'], numel(cyc));
        end
        n = numel(cyc);
        clipped = false;
        for i = 1:n
            a = cyc(mod(i-2, n) + 1);  b = cyc(i);  c = cyc(mod(i, n) + 1);
            if ~isEar(Vi, cyc, a, b, c), continue, end
            tris(end+1, :) = [a b c]; %#ok<AGROW>
            cyc(i) = [];
            clipped = true;
            break
        end
        if ~clipped
            error('PLQ:conjQ:noEar', ...
                'no ear was found, so the boundary is self-intersecting or not counter-clockwise.');
        end
    end
    tris(end+1, :) = cyc;
end

function tf = isEar(Vi, cyc, a, b, c)
% objective: is the corner at b an EAR of the polygon -- a left turn whose triangle is empty. Exact.
    if ratQ.chk(cross2(Vi, a, b, c), 'corner turn') <= 0
        tf = false;  return                  % reflex, or collinear: not an ear
    end
    for i = 1:numel(cyc)
        v = cyc(i);
        if v == a || v == b || v == c, continue, end
        if inTriangle(Vi, a, b, c, v)
            tf = false;  return              % another corner is inside it
        end
    end
    tf = true;
end

function z = cross2(Vi, a, b, c)
% the z component of (b-a) x (c-b), as an exact integer: positive means a LEFT turn at b
    p = Vi(a,:);  q = Vi(b,:);  r = Vi(c,:);
    z = (q(1)-p(1))*(r(2)-q(2)) - (q(2)-p(2))*(r(1)-q(1));
end

function tf = inTriangle(Vi, a, b, c, v)
% objective: does vertex v lie inside the counter-clockwise triangle (a,b,c), boundary included.
% Boundary INCLUDED on purpose: a vertex sitting exactly on an ear's edge makes the clip produce
% two pieces that share more than an edge, so such an ear is refused and another is taken instead.
    tf = ratQ.chk(cross2v(Vi, a, b, v), 'in-triangle') >= 0 && ...
         ratQ.chk(cross2v(Vi, b, c, v), 'in-triangle') >= 0 && ...
         ratQ.chk(cross2v(Vi, c, a, v), 'in-triangle') >= 0;
end

function z = cross2v(Vi, a, b, v)
% the z component of (b-a) x (v-a): positive means v is to the LEFT of the directed edge a->b
    p = Vi(a,:);  q = Vi(b,:);  r = Vi(v,:);
    z = (q(1)-p(1))*(r(2)-p(2)) - (q(2)-p(2))*(r(1)-p(1));
end

function cyc = dropCollinear(Vi, cyc)
    changed = true;
    while changed && numel(cyc) > 3
        changed = false;
        n = numel(cyc);
        for i = 1:n
            a = cyc(mod(i-2, n) + 1);  b = cyc(i);  c = cyc(mod(i, n) + 1);
            if ratQ.chk(cross2(Vi, a, b, c), 'collinear test') == 0
                cyc(i) = [];  changed = true;  break
            end
        end
    end
end

function sh = triangleShape(Vi, vd, tri)
% objective: the shape record for one triangle of the split, in the same form pieceShape produces.
% Convex by construction, so the side of each edge is read straight off the third vertex.
%
% The half-plane rows are emitted in ACTUAL coordinates -- linear part scaled by vd, constant not --
% because that is what eval reads them at, while the SIDE is decided at numerators where the
% vertices live. pieceShape's header explains why mixing the two produces a parallel line.
    sh.vs = sort(tri(:));
    sh.ed = struct('a', {}, 'b', {}, 'd', {}, 'isRay', {});
    sh.rays = zeros(0,2);
    sh.hp = zeros(0,3);
    sh.bounded = true;
    sh.convex = true;
    for i = 1:3
        a = tri(i);  b = tri(mod(i,3) + 1);  o = tri(mod(i+1,3) + 1);
        d = ratQ.chk(Vi(b,:) - Vi(a,:), 'edge direction');
        sh.ed(end+1) = struct('a', a, 'b', b, 'd', d, 'isRay', false); %#ok<AGROW>
        n = [d(2), -d(1)];
        c0 = ratQ.chk(-(n * Vi(a,:).'), 'offset');
        t  = ratQ.chk(n(1)*Vi(o,1) + n(2)*Vi(o,2) + c0, 'third vertex side');
        if t == 0
            error('PLQ:conjQ:degenerate', 'a split triangle has zero area.');
        end
        sh.hp(end+1,:) = sign(t) * [ratQ.chk(vd*n(1), 'edge normal'), ...
                                    ratQ.chk(vd*n(2), 'edge normal'), c0]; %#ok<AGROW>
    end
end
