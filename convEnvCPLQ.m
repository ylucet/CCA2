function r = convEnvCPLQ(obj)
% convEnvCPLQ  Step 1 of the 'cplq' conjugate pipeline: convex envelope of ONE quadratic piece.
%
% objective: compute conv(q + I_P), the convex envelope of a single quadratic q over its
%   (convex) polyhedral domain P, returned as a RatPol. Step 1 of the algorithm in
%     [COAP] Karmarkar & Lucet, Comput. Optim. Appl. 94 (2026) 747-780 (Sect. 3, Appendix A);
%     [JOGO] Karmarkar & Lucet, J. Glob. Optim. 94 (2026) 3-34.
%
% [input]  obj : QuaPoly with a SINGLE piece (nf==1: full domain nv==0, or one face).
%                Operable (quadratic numerator; cubic rejected). A bounded-triangle face must be
%                a valid QuaPoly: with F=[k 0] (face on the left of each directed edge) the
%                vertices are listed counter-clockwise. The single-piece result reuses obj's
%                domain, so an inconsistent winding would make eval miss the interior.
% [output] r   : RatPol = conv(q + I_P).
%
% Classification by the (constant) Hessian Q of q = 1/2 x'Q x + L'x + kappa:
%   * Q positive semidefinite  -> q is convex on the convex set P, so conv = q itself.
%   * Q negative semidefinite over a TRIANGLE -> affine interpolation of the 3 vertex values.
%   * Q indefinite over a TRIANGLE -> reduce to the bilinear form u1*u2 via u = M x (eigen-
%       rotation; for an already-bilinear positive q the map is the identity), compute the
%       envelope with [COAP] Appendix A (0/1/2 convex edges -> affine / rational eq.16 /
%       quadratic harmonic-mean; a triangle with 3 convex edges is split by the smooth-fit line
%       through the middle vertex -- see splitThreeConvex's HISTORY -- into two 2-convex-edge
%       sub-triangles), then substitute u = M x
%       back. Yields linear / quadratic / rational (quadratic over linear) pieces on a polyhedral
%       subdivision (1 or 2 faces). Verified vs Appendix A examples, e.g. conv(xy) over
%       conv{(1,1),(0,0),(2,0)} = 2y^2/(y-x+2), and conv(x^2-y^2) over (1,0),(0,0),(1,1) =
%       (x-y)^2/(1-y).
%
% Multi-piece input (nf>1): each bounded face is fan-triangulated (CCW), the single-triangle
% envelope above is computed per triangle, and the triangle-pieces are merged into one connected
% RatPol (convEnvMultiFace). This is Step 1 of the pipeline: for a single NON-convex quadratic
% over a polygon it returns the per-triangle envelopes (a triangulated intermediate), not the
% polygon's true convex envelope (the latter is recovered later by the conjugate/max steps); for
% convex pieces the per-triangle envelope equals the piece itself, so the result is exact.
%
% Not implemented yet: non-convex over a single non-triangular face (nf==1; pass nf>1 or a
% triangle), non-convex full domain, and unbounded faces. See DESIGN.md II.5.1 and codeOld/cPLQ.

    obj.assertOperable();
    if obj.nf > 1
        r = convEnvMultiFace(obj);   % triangulate each face, envelope per triangle, assemble
        return
    end

    [L,Q,~] = QuaPoly.matrixForm(obj.f(1,:));
    ev    = eig(Q);
    isPSD = all(ev >= -sqrt(eps));
    isNSD = all(ev <=  sqrt(eps));
    f6    = obj.f(1,5:10);                 % quadratic numerator [x^2 xy y^2 x y const] (stored)

    % ---- Convex: envelope = q itself (full domain or any single convex face) ----------------
    if isPSD
        if obj.nv == 0, r = RatPol(f6);
        else,           r = RatPol(obj.V, obj.E, f6, obj.F);
        end
        return
    end
    if obj.nv == 0
        error('convEnvCPLQ:notImplemented', ...
            'Convex envelope of a non-convex full-domain quadratic is not finite; not implemented.');
    end
    isTriangle = (obj.nv == 3) && (obj.ne == 3) && obj.isDomBounded;
    if ~isTriangle
        error('convEnvCPLQ:notImplemented', ...
            ['Convex envelope of a non-convex quadratic is implemented only over a bounded ' ...
             'triangle (got nv=%d, ne=%d, bounded=%d).'], obj.nv, obj.ne, obj.isDomBounded);
    end

    % ---- Concave over a triangle: affine interpolation of the 3 vertex values ---------------
    if isNSD
        z   = QuaPoly.evalPoly(obj.f(1,:), obj.V);
        abc = [obj.V, ones(3,1)] \ z;
        r = RatPol(obj.V, obj.E, [0 0 0 abc(1) abc(2) abc(3)], obj.F);
        return
    end

    % ---- Indefinite over a triangle: reduce to bilinear u1*u2 via u = M x --------------------
    kappa = obj.f(1,end);
    if abs(f6(1)) < sqrt(eps) && abs(f6(3)) < sqrt(eps) && f6(2) > sqrt(eps)
        % already bilinear with positive coefficient: identity frame
        M = eye(2); Minv = eye(2); Vbf = obj.V; beta = f6(2); lin = [f6(4) f6(5) f6(6)];
    else
        % general indefinite (incl. negative bilinear): rotate so 1/2 x'Q x = u1*u2
        M = bilinearFrame(Q); Minv = inv(M); %#ok<MINV>
        Vbf = (M * obj.V')';
        Ltil = (M') \ L;                   % linear part in u-coords: M^{-T} L
        beta = 1; lin = [Ltil(1) Ltil(2) kappa];
    end

    ce  = classifyConvexEdges(Vbf);
    nCE = size(ce,1);
    if nCE == 2
        % [COAP] Appendix A.4 fix (see splitTwoConvexEdges): the single-quadratic formula is not
        % always the tightest envelope over the WHOLE triangle -- split into 2 sub-triangles,
        % UNLESS the special (measure-zero, e.g. mirror-symmetric) case where it already is one.
        [needsSplit, numPlain, den3, triU, faces, edgeList] = splitTwoConvexEdges(Vbf, beta, lin, ce);
        if ~needsSplit
            numX = substituteQuadPlain(numPlain, M);
            denX = substituteLin(den3, M);
            r = RatPol(obj.V, obj.E, plainToStored(numX), obj.F, denX);
            return
        end
        % Each sub-triangle is re-classified and solved via envelopeFromClassified (same pattern
        % as the nCE==3 case below): the sub-triangle keeping P is still exactly 2-convex-edge
        % (reproducing q1 unchanged), and the other sub-triangle is exactly 1-convex-edge (a
        % genuine Appendix A.3 rational piece) -- see splitTwoConvexEdges's header for why.
        numX = zeros(2,6); denX = zeros(2,3);
        for k = 1:2
            sub = triU(faces{k}, :);
            [npk, dk] = envelopeFromClassified(sub, beta, lin, classifyConvexEdges(sub));
            numX(k,:) = plainToStored(substituteQuadPlain(npk, M));
            denX(k,:) = substituteLin(dk, M);
        end
        Vx = (Minv * triU')';
        r  = assembleTwoTriangles(Vx, faces, edgeList, numX, denX);
        return
    end
    if nCE <= 1
        [numPlain, den3] = envelopeFromClassified(Vbf, beta, lin, ce);
        numX = substituteQuadPlain(numPlain, M);
        denX = substituteLin(den3, M);
        r = RatPol(obj.V, obj.E, plainToStored(numX), obj.F, denX);
        return
    end

    % nCE == 3: split the triangle (in the bilinear frame) into two 2-convex-edge sub-triangles
    [triU, faces, edgeList] = splitThreeConvex(Vbf);
    num6 = zeros(2,6); den3 = zeros(2,3);
    for k = 1:2
        sub = triU(faces{k}, :);
        [npk, dk] = envelopeFromClassified(sub, beta, lin, classifyConvexEdges(sub));
        num6(k,:) = plainToStored(substituteQuadPlain(npk, M));
        den3(k,:) = substituteLin(dk, M);
    end
    Vx = (Minv * triU')';                  % sub-triangle vertices back in x-coords
    r  = assembleTwoTriangles(Vx, faces, edgeList, num6, den3);
end

% ===================== local functions ======================================================
function M = bilinearFrame(Q)
% Build M (u = M x) so that 1/2 x'Q x = u1*u2 for an indefinite symmetric Q.
% With Q = lam1 r1 r1' + lam2 r2 r2' (lam1>0>lam2) and a = sqrt(lam1/2) r1 + sqrt(-lam2/2) r2,
% b = sqrt(lam1/2) r1 - sqrt(-lam2/2) r2, one has a b' + b a' = Q, hence (a'x)(b'x) = 1/2 x'Q x.
    [R, Lam] = eig(Q); lam = diag(Lam);
    [lp, ip] = max(lam); [ln, in] = min(lam);
    r1 = R(:,ip); r2 = R(:,in);
    a = sqrt(lp/2)*r1 + sqrt(-ln/2)*r2;
    b = sqrt(lp/2)*r1 - sqrt(-ln/2)*r2;
    M = [a'; b'];
end

function ce = classifyConvexEdges(V)
% Convex edges of a triangle for the bilinear form u1*u2: an edge is convex iff its slope > 0
% (u1*u2 restricted to y=m x+q is m x^2+q x, convex iff m>0). Vertical edges are not convex.
% Returns ce: (#convex) x 3 rows [slope, intercept, opposite-vertex-index].
    tol = sqrt(eps); ed = [1 2; 2 3; 3 1]; ce = zeros(0,3);
    for t = 1:3
        i = ed(t,1); j = ed(t,2); opp = 6 - i - j;
        vi = V(i,:); vj = V(j,:); dx = vj(1)-vi(1); dy = vj(2)-vi(2);
        if abs(dx) < tol, continue; end
        m = dy/dx;
        if m > tol
            ce(end+1,:) = [m, vi(2)-m*vi(1), opp]; %#ok<AGROW>
        end
    end
end

function [numPlain, den3] = envelopeFromClassified(V, beta, lin, ce)
% conv(beta*u1*u2 + lin) over the triangle V, given the convex-edge list ce (0/1/2 rows).
% Returns PLAIN numerator coefficients [a b c d e f] (a u1^2+b u1u2+c u2^2+d u1+e u2+f) and a
% linear denominator den3 = [g h k] (den3 = [0 0 1] when the envelope is a polynomial).
    d = lin(1); e = lin(2); f0 = lin(3);
    switch size(ce,1)
        case 0   % Appendix A.2: affine interpolation of the vertex values
            qv  = beta*(V(:,1).*V(:,2)) + d*V(:,1) + e*V(:,2) + f0;
            abc = [V, ones(3,1)] \ qv;
            numPlain = [0 0 0 abc(1) abc(2) abc(3)]; den3 = [0 0 1];
        case 1   % Appendix A.3 eq.16: rational (quadratic over linear)
            m = ce(1,1); q = ce(1,2); x1 = V(ce(1,3),1); y1 = V(ce(1,3),2);
            a = -m*y1; b = q; c = x1;
            dd = -q*y1 + m*x1*y1; ee = -q*x1 - x1*y1; ff = q*x1*y1;
            g = -m; h = 1; k = -y1 + m*x1;
            numB = [a b c dd ee ff]; denB = [g h k];
            linDen = [d*g, d*h+e*g, e*h, d*k+f0*g, e*k+f0*h, f0*k];
            numPlain = beta*numB + linDen; den3 = denB;
        case 2   % Appendix A.4: quadratic (harmonic-mean form)
            q2 = twoEdgeQuadPlain(ce(1,1), ce(1,2), ce(2,1), ce(2,2), V);
            numPlain = beta*q2 + [0 0 0 d e f0]; den3 = [0 0 1];
        otherwise
            error('convEnvCPLQ:internal','envelopeFromClassified expects 0,1,or 2 convex edges.');
    end
end

function [needsSplit, num6, den3, tri, faces, edgeList] = splitTwoConvexEdges(V, beta, lin, ce)
% [COAP] Appendix A.4 FIX (2026 session, corrected in a later session -- see below): the single
% quadratic q1 = twoEdgeQuadPlain(...), which touches u1*u2 along BOTH classified convex edges, is
% a valid convex minorant but NOT always the tightest one over the WHOLE triangle. Along the
% triangle's third ("weak", non-convex) edge, u1*u2 is CONCAVE (slope <= 0), so the true envelope
% there is the AFFINE CHORD between its two endpoint u1*u2 values -- strictly ABOVE q1 in the open
% interior of that edge, since q1 is a single convex quadratic forced through only the 2 shared
% endpoints, not the whole edge. Reproduced and confirmed on the paper's own Appendix A.4.3
% example -- see DESIGN.md / session handoff for the diagnosis this fix resolves.
%
% Fix: split the triangle by a cevian from ONE weak-edge endpoint into the OPPOSITE convex edge.
% The sub-triangle containing the two convex edges' common vertex P keeps q1 UNCHANGED (it still
% touches u1*u2 exactly along both of ITS two convex-edge sub-segments, since each is a
% sub-segment of an original convex edge's own infinite line). The OTHER sub-triangle's own
% remaining edges are: one sub-segment of the OTHER original convex edge (still convex, same
% line), the fully-contained original weak edge, and the new internal seam -- exactly the
% Appendix A.3 one-convex-edge case, NOT a bespoke construction. (An earlier version of this fix
% used a bespoke quadratic here, touching u1*u2 along the remaining convex edge and the affine
% CHORD along the weak edge; that is a valid boundary condition on the triangle's BOUNDARY but not
% enough to pin down the unique tight envelope in the INTERIOR -- confirmed wrong by an
% independent ground-truth check (3D convex hull of a dense sample of the graph over a
% counterexample triangle, gap up to ~2.5): the true envelope over that sub-triangle is not a
% single quadratic at all, so no single boundary-matching quadratic can equal it everywhere. The
% caller (convEnvCPLQ's nCE==2 branch) now re-classifies and re-solves BOTH sub-triangles via
% envelopeFromClassified, exactly like the nCE==3 case below -- this function only returns the
% split's GEOMETRY (needsSplit, tri, faces, edgeList), not a final formula, except in the no-split
% case. Re-verified to match ground truth to grid resolution on the original counterexample and on
% a fresh random sample of 2-convex-edge triangles -- see session handoff.)
%
% The cevian's direction is forced, not free -- exactly analogous to splitThreeConvex: it is the
% unique direction making the two pieces agree identically along the shared internal seam. q1 and
% ANY quadratic touching u1*u2 along the SAME original convex edge (the one being split) agree
% identically along that edge's whole line (both equal u1*u2 there) and therefore agree exactly at
% the triangle's third vertex too (both equal u1*u2 there as well); q1 minus such a quadratic is
% therefore a quadratic vanishing along that whole line, hence factors as (that line) * (a second
% line) that must pass through the shared opposite vertex -- the second line is the seam, found by
% seamPoint. `buildEdgeAffinePiece` supplies one convenient such placeholder quadratic (matching
% u1*u2 along the convex edge and the affine chord along the weak edge) purely to locate this
% seam -- its own value is NOT used as a final envelope piece (see above). Of the two candidate
% cevians (from the weak edge's first endpoint into the second convex edge, or from its second
% endpoint into the first), exactly one lands strictly inside its target edge and the other's
% intersection falls outside it -- verified across ~6500 random 2-convex-edge triangles in the
% original session (never both valid, never neither); used here as the selection rule, mirroring
% how twoEdgeQuadPlain already picks its own +/- branch by validity.
%
% Special (measure-zero) case: if q1 itself is already exactly affine (in the edge's own
% parameter) along the weak edge, it already touches the chord everywhere -- e.g. the mirror-
% symmetric triangle (0,0),(2,1),(1,2), where mh*mw=1 makes q1 constant along the weak edge -- no
% split is needed (needsSplit=false); both candidate cevians degenerate (their seam runs parallel
% to the target edge, so seamPoint's intersection is at infinity) in exactly this case, but the
% weak-edge curvature check below is the direct, cheap way to detect it up front.
    Bidx = ce(1,3); Aidx = ce(2,3); Pidx = setdiff(1:3, [Bidx Aidx]);
    P = V(Pidx,:); A = V(Aidx,:); B = V(Bidx,:);
    mh = ce(1,1); qh = ce(1,2);      % edge P-A
    mw = ce(2,1); qw = ce(2,2);      % edge P-B

    q1 = twoEdgeQuadPlain(mh, qh, mw, qw, V);       % unchanged: touches u1*u2 along both P-A, P-B

    dxAB = B(1)-A(1); dyAB = B(2)-A(2);
    curv = q1(1)*dxAB^2 + q1(2)*dxAB*dyAB + q1(3)*dyAB^2;   % q1's curvature along the weak edge
    scale = (abs(q1(1))+abs(q1(2))+abs(q1(3))+1) * (dxAB^2+dyAB^2);
    if abs(curv) < 1e-9*scale
        needsSplit = false;
        num6 = beta*q1 + [0 0 0 lin(1) lin(2) lin(3)]; den3 = [0 0 1];
        tri = []; faces = {}; edgeList = [];
        return
    end
    needsSplit = true;
    num6 = []; den3 = [];   % recomputed by the caller, one envelopeFromClassified call per face

    qFromA = buildEdgeAffinePiece(mw, qw, A, B);     % seam-finding aid only, see header above
    qFromB = buildEdgeAffinePiece(mh, qh, A, B);     % seam-finding aid only, see header above

    [Ra, ta] = seamPoint(q1 - qFromA, mw, qw, P, B); % candidate: cevian from A into edge P-B
    if ta > 1e-9 && ta < 1 - 1e-9
        tri = [P; A; B; Ra];
        faces = {[1 2 4], [2 4 3]};                  % T1={P,A,R} (2CE); T2={A,R,B} (1CE)
        edgeList = [1 2; 1 4; 4 3; 2 3; 2 4];         % P-A, P-R, R-B, A-B, A-R(internal seam)
    else
        [Rb, tb] = seamPoint(q1 - qFromB, mh, qh, P, A); % candidate: cevian from B into edge P-A
        if ~(tb > 1e-9 && tb < 1 - 1e-9)
            error('convEnvCPLQ:internal', ...
                'splitTwoConvexEdges: neither candidate cevian lands inside its target edge.');
        end
        tri = [P; A; B; Rb];
        faces = {[1 4 3], [4 2 3]};                  % T1={P,R,B} (2CE); T2={R,A,B} (1CE)
        edgeList = [1 3; 1 4; 4 2; 2 3; 3 4];         % P-B, P-R, R-A, A-B, B-R(internal seam)
    end
end

function coef = buildEdgeAffinePiece(m1, q1c, Ap, Bp)
% Seam-finding aid for splitTwoConvexEdges (NOT itself a final envelope piece -- see that
% function's header): conv(u1*u2) that touches u1*u2 exactly along the convex edge y=m1*x+q1c AND
% touches, along the (possibly vertical) line through Ap,Bp, the AFFINE CHORD interpolating
% u1*u2's values at Ap and Bp. Used only so that q1 minus this placeholder vanishes identically
% along the shared convex edge (both this placeholder and q1 equal u1*u2 there), which is what
% makes seamPoint's line-factoring argument locate the correct cevian.
% Both "match along a whole line" conditions are linear in the quadratic's 6 plain coefficients,
% but together are RANK-DEFICIENT BY EXACTLY 1 (same pencil structure as twoEdgeQuadPlain's own
% derivation): any solution plus a multiple of ellLine1*ellLineAB (the product of the two matched
% lines) also satisfies both. Unlike twoEdgeQuadPlain (2 distinct rank-1/PSD-boundary solutions,
% needing a +/- branch choice), this pencil's rank-1 (Hessian-singular) condition is a DOUBLE root
% -- tangent, not crossing -- always picking a UNIQUE convex quadratic (verified symbolically and
% numerically this session), so no branch selection is needed here.
    Ax = Ap(1); Ay = Ap(2); Bx = Bp(1); By = Bp(2);
    dx = Bx - Ax; dy = By - Ay; fA = Ax*Ay; fB = Bx*By;
    Amat = [ 1, m1, m1^2, 0, 0, 0;
             0, q1c, 2*m1*q1c, 1, m1, 0;
             0, 0, q1c^2, 0, q1c, 1;
             dx^2, dx*dy, dy^2, 0, 0, 0;
             2*Ax*dx, Ax*dy+Ay*dx, 2*Ay*dy, dx, dy, 0;
             Ax^2, Ax*Ay, Ay^2, Ax, Ay, 1 ];
    rhs = [m1; q1c; 0; 0; fB-fA; fA];
    sol0 = pinv(Amat) * rhs;                          % minimum-norm particular solution (rank 5)

    K = dx*Ay - dy*Ax;                                 % null-space direction: ellLine1*ellLineAB
    Lnull = [-m1*dy, m1*dx+dy, -dx, -m1*K-q1c*dy, K+q1c*dx, -q1c*K];

    a0=sol0(1); b0=sol0(2); c0=sol0(3);
    aL=Lnull(1); bL=Lnull(2); cL=Lnull(3);
    A2 = 4*aL*cL - bL^2;
    A1 = 4*(a0*cL + aL*c0) - 2*b0*bL;
    lam = -A1/(2*A2);
    coef = sol0' + lam*Lnull;
end

function [R, t] = seamPoint(diffCoef, m0, q0, Pp, otherEnd)
% diffCoef (plain [a b c d e f]) vanishes identically along y=m0*x+q0 (m0 is a classified convex
% edge's slope, always > 0, so never vertical); it therefore factors as
% (y-m0*x-q0)*(p*x+q*y+r) -- solve for [p q r] directly from diffCoef's own coefficients (exact,
% no risk of dividing by a vanishing quantity since m0~=0) -- then intersect that second line with
% y=m0*x+q0 to get the seam's far endpoint R. t is R's fractional position from Pp to otherEnd
% along that same edge (0<t<1 required for R to lie strictly inside it).
    a = diffCoef(1); c = diffCoef(3); dd = diffCoef(4);
    p = -a/m0; q = c; r = -(dd + q0*p)/m0;
    denom = p + q*m0;
    xR = -(q*q0 + r)/denom;
    yR = m0*xR + q0;
    R = [xR, yR];
    if abs(otherEnd(1)-Pp(1)) > abs(otherEnd(2)-Pp(2))
        t = (xR - Pp(1)) / (otherEnd(1) - Pp(1));
    else
        t = (yR - Pp(2)) / (otherEnd(2) - Pp(2));
    end
end

function q = twoEdgeQuadPlain(mh, qh, mw, qw, V)
% conv(u1*u2) over a triangle with two convex edges (Appendix A.4), plain coeffs [a b c d e f].
% Of the two (+/-) solutions, pick the one that touches f=u1*u2 along a convex edge.
    for s = [1 -1]
        cand = buildTwoEdge(mh, qh, mw, qw, s);
        xt = mean(V(:,1)); yt = mh*xt + qh;
        if abs(evalPlain(cand, xt, yt) - xt*yt) < 1e-7*(1+abs(xt*yt))
            q = cand; return;
        end
    end
    q = buildTwoEdge(mh, qh, mw, qw, 1);
end

function c = buildTwoEdge(mh, qh, mw, qw, s)
    rr = sqrt(mh*mw); denom = mh + mw + s*2*rr;
    c = [ mh*mw/denom, s*2*rr/denom, 1/denom, ...
          (mh*qw + mw*qh)/denom, -(qh+qw)/denom, qh*qw/denom ];
end

function v = evalPlain(c, x, y)
    v = c(1)*x.^2 + c(2)*x.*y + c(3)*y.^2 + c(4)*x + c(5)*y + c(6);
end

function s6 = plainToStored(p6)
% plain [A B C D E F] (A x^2+...) -> stored weighted basis [2A B 2C D E F] used by evalPoly.
    s6 = [2*p6(1), p6(2), 2*p6(3), p6(4), p6(5), p6(6)];
end

function p = substituteQuadPlain(N, M)
% Substitute u = M x into a PLAIN quadratic N=[A B C D E F] (in u1,u2); return PLAIN coeffs in x,y.
% u1 = M(1,1) x + M(1,2) y,  u2 = M(2,1) x + M(2,2) y.
    A=N(1); B=N(2); C=N(3); D=N(4); E=N(5); F=N(6);
    m11=M(1,1); m12=M(1,2); m21=M(2,1); m22=M(2,2);
    p = [ A*m11^2 + B*m11*m21 + C*m21^2, ...
          A*2*m11*m12 + B*(m11*m22+m12*m21) + C*2*m21*m22, ...
          A*m12^2 + B*m12*m22 + C*m22^2, ...
          D*m11 + E*m21, ...
          D*m12 + E*m22, ...
          F ];
end

function d3 = substituteLin(den, M)
% Substitute u = M x into a linear denominator den=[g h k] (g u1+h u2+k); return [.. .. k] in x,y.
    g=den(1); h=den(2); k=den(3);
    d3 = [ g*M(1,1) + h*M(2,1), g*M(1,2) + h*M(2,2), k ];
end

function [tri, faces, edgeList] = splitThreeConvex(V)
% Split a triangle (3x2, in the bilinear frame) by the SMOOTH-FIT line through the middle (by
% 2nd-coordinate) vertex, hitting the opposite (low-high) edge at Pnew. Returns 4 vertices
% tri = [vlow; vmid; vhigh; Pnew], the two sub-triangles faces (vertex indices into tri), and a
% directed edge list. Each sub-triangle then has exactly two convex edges.
%
% The split line is NOT simply horizontal through vmid (an earlier version of this function used
% that, which is wrong in general -- see HISTORY below). Both sub-envelopes q1=conv(u1u2+I_lowMidPnew)
% and q2=conv(u1u2+I_midHighPnew) touch u1*u2 exactly along the ENTIRE low-high edge (the "w" edge
% of Appendix A.4's two-edge formula is shared by both sub-triangles), so q1-q2 is a quadratic that
% vanishes identically on that whole line, hence factors as (low-high line) * (a second line). That
% second line necessarily passes through vmid too (q1(vmid)=q2(vmid)=u1u2(vmid), both touching
% there via their OTHER convex edge), and is the unique choice of split direction through vmid that
% makes q1 and q2 agree identically along the shared internal seam (smooth C^1 pasting, required
% for the glued 2-piece function to be the TRUE convex envelope, not just two locally-correct but
% badly-glued pieces). Only for the special case sqrt(mh1)==sqrt(mh2) (mirror-symmetric split) does
% this second line reduce to horizontal.
%
% HISTORY: the original "horizontal line through vmid" (matching [COAP] Appendix A.5's likely
% intent for a specific/symmetric configuration, but not the general condition) produced q1,q2 that
% only agreed with each other at the two seam endpoints (vmid, Pnew), not along the interior of the
% seam -- a real but small (~0.02-0.04) mismatch for the well-tested T=(0,0),(3,3),(1,2) example
% (small enough that it never affected that example's tested conjugate values), but a LARGE
% (~1-13) mismatch for other triangles (e.g. T=(0,0),(7.02,0.67),(8.43,7.63)), making the glued
% 2-piece "envelope" genuinely non-convex (a real jump discontinuity across the interior seam) and
% its conjugate (via maxQuaPar) numerically too high. Found by comparing the two sub-envelopes'
% closed forms along their shared seam for a stress-test triangle where the final conjugate was
% wrong, discovering they only match at the seam's 2 endpoints, then deriving the correct split
% line as the second factor of q1-q2 above (verified to give exact, to machine precision, agreement
% along the WHOLE seam for both the known-good example and every previously-failing triangle).
    [~, ord] = sort(V(:,2));
    vlow = V(ord(1),:); vmid = V(ord(2),:); vhigh = V(ord(3),:);

    mh1 = (vmid(2)-vlow(2))/(vmid(1)-vlow(1));  qh1 = vlow(2) - mh1*vlow(1);
    mw  = (vhigh(2)-vlow(2))/(vhigh(1)-vlow(1)); qw = vlow(2) - mw*vlow(1);
    mh2 = (vhigh(2)-vmid(2))/(vhigh(1)-vmid(1)); qh2 = vmid(2) - mh2*vmid(1);

    c1 = buildTwoEdge(mh1, qh1, mw, qw, 1);
    c2 = buildTwoEdge(mh2, qh2, mw, qw, 1);
    d = c1 - c2;                                % plain [A B C D E F], vanishes on the low-high line
    q = -d(3); p = mw*q - d(2); r = (d(4) - qw*p)/mw;   % second factor p*x+q*y+r=0, through vmid

    xPnew = -(q*qw + r) / (p + q*mw);           % intersect p x+q y+r=0 with y = mw x + qw
    yPnew = mw*xPnew + qw;
    Pnew  = [xPnew, yPnew];

    tri = [vlow; vmid; vhigh; Pnew];          % indices 1=low, 2=mid, 3=high, 4=Pnew
    faces = {[1 2 4], [2 3 4]};               % T1={low,mid,Pnew}, T2={mid,high,Pnew}
    edgeList = [1 2; 2 3; 1 4; 4 3; 2 4];     % low-mid, mid-high, low-Pnew, Pnew-high, mid-Pnew(internal)
end

function r = assembleTwoTriangles(V, faces, edgeList, num6, den3)
% Build a 2-face RatPol from two triangles sharing an edge. Edge orientation (which face is on
% the left/right) is determined by a centroid side test, so F is consistent for the constructor.
    ne = size(edgeList,1); nf = numel(faces);
    E = [edgeList, ones(ne,1)];
    cent = zeros(nf,2);
    for k = 1:nf, cent(k,:) = mean(V(faces{k},:),1); end
    F = zeros(ne,2);
    for j = 1:ne
        a = edgeList(j,1); b = edgeList(j,2); va = V(a,:); dir = V(b,:) - va;
        for k = 1:nf
            if all(ismember([a b], faces{k}))
                cr = dir(1)*(cent(k,2)-va(2)) - dir(2)*(cent(k,1)-va(1)); % (b-a) x (cent-a)
                if cr > 0, F(j,1) = k; else, F(j,2) = k; end
            end
        end
    end
    r = RatPol(V, E, num6, F, den3);
end

% ----- multi-face (Step 1 over a triangulated domain) ---------------------------------------
function r = convEnvMultiFace(obj)
% Fan-triangulate every bounded face of obj, compute the single-triangle convex envelope of the
% face's quadratic on each triangle, and merge the resulting triangle-pieces into one RatPol.
    pieces = struct('V', {}, 'num6', {}, 'den3', {});
    for i = 1:obj.nf
        tris = extractFaceTrianglesCCW(obj, i);
        q6   = obj.f(i, 5:10);                 % the (quadratic) numerator on face i
        for t = 1:numel(tris)
            qT = QuaPoly(tris{t}, [1 2 1; 2 3 1; 3 1 1], q6, [1 0; 1 0; 1 0]);
            rT = convEnvCPLQ(qT);              % single-triangle path (1 or 2 output faces)
            for g = 1:rT.nf
                Vg = ensureCCW(rT.V(faceVertexIndices(rT, g), :));
                pieces(end+1) = struct('V', Vg, 'num6', rT.f(g,5:10), 'den3', rT.den(g,:)); %#ok<AGROW>
            end
        end
    end
    r = assembleTriangles(pieces);
end

function tris = extractFaceTrianglesCCW(obj, i)
% Vertices of bounded face i in boundary order -> CCW polygon -> fan triangulation (convex face).
    ej = find(any(obj.F == i, 2));             % edges incident to face i
    if any(obj.E(ej,3) == 0)
        error('convEnvCPLQ:notImplemented', ...
            'Multi-piece convex envelope requires bounded faces (face %d is unbounded).', i);
    end
    W = obj.V(faceVertexIndices(obj, i), :);
    if signedArea(W) < 0, W = flipud(W); end   % make the polygon CCW
    n = size(W,1);
    if n < 3, error('convEnvCPLQ:degenerateFace', 'Face %d has fewer than 3 vertices.', i); end
    tris = cell(1, n-2);
    for t = 1:n-2, tris{t} = [W(1,:); W(t+1,:); W(t+2,:)]; end   % fan from W(1)
end

function iVs = faceVertexIndices(obj, k)
% Vertex indices around face k, in the order of its ordered edge list obj.P{k}.
    face = obj.P{k}; iVs = zeros(1, numel(face));
    for i = 1:numel(face)
        j = face(i);
        if j > 0, iVs(i) = obj.E(j,1); else, iVs(i) = obj.E(-j,2); end
    end
end

function W = ensureCCW(W)
    if signedArea(W) < 0, W = W([1 3 2], :); end   % 3 vertices: swap to CCW
end

function a = signedArea(W)
    x = W(:,1); y = W(:,2); n = size(W,1); a = 0;
    for i = 1:n, j = mod(i,n)+1; a = a + (x(i)*y(j) - x(j)*y(i)); end
    a = a/2;
end

function r = assembleTriangles(pieces)
% Merge CCW triangle pieces {V(3x2), num6, den3} into one connected RatPol. Vertices are
% deduplicated by coordinate; an undirected edge shared by two CCW triangles is internal (the
% face that traverses it as a->b is on its left, the other on its right); boundary edges keep 0.
    tol = sqrt(eps); nf = numel(pieces);
    allV = zeros(3*nf, 2);
    for k = 1:nf, allV(3*k-2:3*k, :) = pieces(k).V; end
    V = zeros(0,2); idx = zeros(3*nf,1);
    for p = 1:3*nf
        hit = 0;
        for q = 1:size(V,1)
            if norm(V(q,:) - allV(p,:)) < tol, idx(p) = q; hit = 1; break; end
        end
        if ~hit, V(end+1,:) = allV(p,:); idx(p) = size(V,1); end %#ok<AGROW>
    end
    faceV = reshape(idx, 3, nf)';              % nf x 3 CCW vertex indices
    HE = zeros(3*nf, 3); c = 0;                % half-edges [face a b]
    for k = 1:nf
        vs = faceV(k,:);
        for t = 1:3
            c = c+1; HE(c,:) = [k, vs(t), vs(mod(t,3)+1)];
        end
    end
    used = false(3*nf,1); E = zeros(0,3); F = zeros(0,2);
    for c = 1:3*nf
        if used(c), continue; end
        used(c) = true; k = HE(c,1); a = HE(c,2); b = HE(c,3);
        opp = find(~used & HE(:,2)==b & HE(:,3)==a, 1);
        if isempty(opp)
            E(end+1,:) = [a b 1]; F(end+1,:) = [k 0]; %#ok<AGROW>
        else
            used(opp) = true;
            E(end+1,:) = [a b 1]; F(end+1,:) = [k HE(opp,1)]; %#ok<AGROW>
        end
    end
    num6 = zeros(nf,6); den3 = zeros(nf,3);
    for k = 1:nf, num6(k,:) = pieces(k).num6; den3(k,:) = pieces(k).den3; end
    r = RatPol(V, E, num6, F, den3);
end
