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
        [needsSplit, numPlain, den3, triU, faces, ~] = splitTwoConvexEdges(Vbf, beta, lin, ce);
        if ~needsSplit
            numX = substituteQuadPlain(numPlain, M);
            denX = substituteLin(den3, M);
            r = RatPol(obj.V, obj.E, plainToStored(numX), obj.F, denX);
            return
        end
        % Each sub-triangle is itself solved via solveTriangleBF, NOT envelopeFromClassified
        % directly: the sub-triangle keeping P is again exactly 2-convex-edge (same mh,mw), and
        % REMAINS subject to the exact same tightness issue splitTwoConvexEdges itself fixes --
        % solveTriangleBF recurses into it (and the other sub-triangle) exactly as needed. See
        % solveTriangleBF's header and DESIGN.md's Part 2c for why a single further split of the
        % ORIGINAL triangle is not always sufficient in general.
        piecesBF = [solveTriangleBF(triU(faces{1},:), beta, lin), ...
                    solveTriangleBF(triU(faces{2},:), beta, lin)];
        r = assemblePiecesBF(piecesBF, M, Minv);
        return
    end
    if nCE <= 1
        [numPlain, den3] = envelopeFromClassified(Vbf, beta, lin, ce);
        numX = substituteQuadPlain(numPlain, M);
        denX = substituteLin(den3, M);
        r = RatPol(obj.V, obj.E, plainToStored(numX), obj.F, denX);
        return
    end

    % nCE == 3: split the triangle (in the bilinear frame) into two 2-convex-edge sub-triangles,
    % each solved via solveTriangleBF (see above -- either sub-triangle can itself need further
    % splitting, exactly like the nCE==2 case, so this must NOT call envelopeFromClassified
    % directly on them; see DESIGN.md's Part 2c for a repro where skipping this recursion silently
    % gave a wrong answer).
    [triU, faces, ~] = splitThreeConvex(Vbf);
    piecesBF = [solveTriangleBF(triU(faces{1},:), beta, lin), ...
                solveTriangleBF(triU(faces{2},:), beta, lin)];
    r = assemblePiecesBF(piecesBF, M, Minv);
end

function pieces = solveTriangleBF(V, beta, lin)
% Solve conv(beta*u1*u2+lin) over a triangle V (bilinear frame), returning a cell array of pieces
% {V, num (plain), den (plain)} in the SAME bilinear frame -- a single piece if V has 0 or 1
% convex edges, or 2 convex edges but no further split is needed; RECURSIVELY more pieces if a
% 2-convex-edge triangle's own single-quadratic formula (q1) is not tight throughout it (see
% splitTwoConvexEdges). This is what lets EITHER of splitTwoConvexEdges' own two sub-triangles, or
% either of splitThreeConvex's two sub-triangles, need a further split themselves -- both produce
% exactly-2-convex-edge sub-triangles, subject to the identical tightness issue as the original
% triangle, so both must go through this same check rather than assuming a single split suffices.
    ce = classifyConvexEdges(V);
    if size(ce,1) == 2
        [needsSplit, numPlain, den3, triU, faces, ~] = splitTwoConvexEdges(V, beta, lin, ce);
        if needsSplit
            pieces = [solveTriangleBF(triU(faces{1},:), beta, lin), ...
                      solveTriangleBF(triU(faces{2},:), beta, lin)];
            return
        end
    else
        [numPlain, den3] = envelopeFromClassified(V, beta, lin, ce);
    end
    pieces = {struct('V', V, 'num', numPlain, 'den', den3)};
end

function r = assemblePiecesBF(piecesBF, M, Minv)
% Transform a cell array of bilinear-frame {V,num,den} pieces (from solveTriangleBF) back to x,y
% coordinates and assemble them into one connected RatPol via assembleTriangles (which derives
% face adjacency from shared vertices/edges, so it handles any number of pieces, not just 2).
    n = numel(piecesBF);
    pieces = struct('V', {}, 'num6', {}, 'den3', {});
    for i = 1:n
        p = piecesBF{i};
        Vx = ensureCCW((Minv * p.V')');
        numX = plainToStored(substituteQuadPlain(p.num, M));
        denX = substituteLin(p.den, M);
        pieces(end+1) = struct('V', Vx, 'num6', numX, 'den3', denX); %#ok<AGROW>
    end
    r = assembleTriangles(pieces);
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
% [COAP] Appendix A.4 FIX, Parts 1-2c (2026 sessions -- see DESIGN.md for the full derivation
% trail): the single quadratic q1 = twoEdgeQuadPlain(...), which touches u1*u2 along BOTH
% classified convex edges, is a valid convex minorant but NOT always the tightest one over the
% WHOLE triangle -- it is only tight over a genuine SUB-region containing P, not the whole
% triangle. Reproduced and confirmed on the paper's own Appendix A.4.3 example, whose own claim
% ("the domain is the entire triangle") is false there too: q1(0.474343,0) = -0.042780 instead of
% the true value 0 (f=xy is identically 0, hence >=0 is a valid global minorant, along that whole
% weak edge).
%
% Fix: split the triangle by a cevian from ONE convex edge's far vertex into the OPPOSITE convex
% edge. The sub-triangle containing the two convex edges' common vertex P keeps q1 UNCHANGED. The
% OTHER sub-triangle's own remaining edges are: one sub-segment of the OTHER original convex edge
% (still convex, same line), the ORIGINAL weak edge's anchor vertex reclassified as the new
% triangle's far vertex, and the new internal seam (itself always non-convex, checked on many
% random triangles) -- exactly the Appendix A.3 one-convex-edge case, not a bespoke construction.
%
% The cevian's direction is FORCED (not free), and is exactly the unique direction making q1 and
% the OTHER sub-triangle's Appendix A.3 formula agree not just in VALUE but in GRADIENT along their
% shared seam (a genuine C1, tangent contact, matching the smooth-fit pattern splitThreeConvex
% already uses for the 3-convex-edge case). Derived by writing q1 - R_anchor (R_anchor = the
% Appendix A.3 formula anchored at the far vertex of ONE convex edge, using the OTHER edge) as a
% single rational expression and clearing denominators: using the fact that the anchor vertex lies
% on its OWN edge's line, the resulting polynomial factors as (the used edge's line) times a PERFECT
% SQUARE, (sqrt(mh*mw)*(x-x0) + (y-y0))^2, where (x0,y0) is the anchor vertex -- a double root, i.e.
% tangency, not a transversal crossing. This means the correct cevian is simply the LINE THROUGH THE
% FAR VERTEX WITH SLOPE -sqrt(mh*mw), intersected with the other convex edge -- remarkably, this
% slope depends ONLY on mh and mw (beta and the affine shift `lin` cancel out of the derivation
% completely, an exact algebraic fact, not an approximation). Verified against ground truth
% (numerically maximized biconjugate) to solver precision, replacing an EARLIER, closely related but
% still not fully correct criterion from this same session (matching only ONE classified edge's own
% dual-tangency condition, `s_j(x,y) = far-endpoint x-coordinate`): that criterion is a valid
% necessary condition (confirmed: q1 and the Appendix A.3 formula DO agree in value, though not
% gradient, along its line) but not sufficient -- it under-corrects, still leaving a thin residual
% region (found via a full barycentric grid scan) where q1 undershoots truth by as much as ~0.11 in
% one repro triangle, immediately adjacent to that first-pass cevian. The gradient-tangency
% criterion here supersedes it entirely and closes that residual gap (verified exactly, both
% sub-regions matching ground truth on multiple repro triangles).
%
% Of the two candidate cevians (from one convex edge's far vertex into the other edge, or vice
% versa), exactly one lands strictly inside its target edge and the other's intersection falls
% outside it -- same selection rule as before, mirroring how twoEdgeQuadPlain already picks its own
% +/- branch by validity.
%
% Special (measure-zero) case: if q1 itself is already exactly affine (in the edge's own
% parameter) along the weak edge, it already touches the chord everywhere -- e.g. the mirror-
% symmetric triangle (0,0),(2,1),(1,2), where mh*mw=1 makes q1 constant along the weak edge -- no
% split is needed (needsSplit=false); both candidate cevians degenerate (their direction runs
% parallel to the target edge, so the intersection is at infinity) in exactly this case, but the
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

    slope = -sqrt(mh*mw);   % forced tangency direction; independent of beta, lin (see header)

    % candidate: cevian from A (edge h's far vertex), slope -sqrt(mh*mw), into edge P-B
    [Ra, ta] = tangentCevian(A, slope, mw, qw, P, B);
    if ta > 1e-9 && ta < 1 - 1e-9
        tri = [P; A; B; Ra];
        faces = {[1 2 4], [2 4 3]};                  % T1={P,A,R} (2CE); T2={A,R,B} (1CE)
        edgeList = [1 2; 1 4; 4 3; 2 3; 2 4];         % P-A, P-R, R-B, A-B, A-R(internal seam)
    else
        % candidate: cevian from B (edge w's far vertex), slope -sqrt(mh*mw), into edge P-A
        [Rb, tb] = tangentCevian(B, slope, mh, qh, P, A);
        if ~(tb > 1e-9 && tb < 1 - 1e-9)
            error('convEnvCPLQ:internal', ...
                'splitTwoConvexEdges: neither candidate cevian lands inside its target edge.');
        end
        tri = [P; A; B; Rb];
        faces = {[1 4 3], [4 2 3]};                  % T1={P,R,B} (2CE); T2={R,A,B} (1CE)
        edgeList = [1 3; 1 4; 4 2; 2 3; 3 4];         % P-B, P-R, R-A, A-B, B-R(internal seam)
    end
end

function [R, t] = tangentCevian(anchor, slope, mTarget, qTarget, Pp, otherEnd)
% Line through `anchor` (a triangle vertex) with the given `slope`, intersected with the OTHER
% convex edge y = mTarget*x + qTarget -- see splitTwoConvexEdges' header for the derivation of why
% this direction (not `anchor`'s own position alone) locates the correct cevian. t is R's
% fractional position from Pp to otherEnd along the target edge (0<t<1 required for R to lie
% strictly inside it).
    xA = anchor(1); yA = anchor(2);
    xR = (yA - slope*xA - qTarget) / (mTarget - slope);
    yR = mTarget*xR + qTarget;
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
