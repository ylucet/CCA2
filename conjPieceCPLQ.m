function g = conjPieceCPLQ(p)
% conjPieceCPLQ  Step 2 of the 'cplq' pipeline: Fenchel conjugate of a single envelope piece.
%
% objective: conjugate (p + I_P)*(s) = sup_{x in P} <s,x> - p(x) of one convex-envelope piece,
%   returned as a QuaPar (quadratic on a parabolic subdivision; here polyhedral). Step 2 of
%     [COAP] Karmarkar & Lucet, Comput. Optim. Appl. 94 (2026) 747-780 (Appendix B);
%     [JOGO] Karmarkar & Lucet, J. Glob. Optim. 94 (2026) 3-34.
%
% [input]  p : a single piece over a bounded triangle. Currently implemented:
%              - AFFINE piece  -> f* = max_i(<s,v_i>-ell(v_i)), three cones;
%              - CONVEX (positive-definite) QUADRATIC -> seven polyhedral regions;
%              - CONVEX, RANK-1 (PSD, one zero eigenvalue) QUADRATIC -> six-face QuaPar (three
%                PARABOLIC edge strips + three linear vertex cones). This is the case actually
%                produced by convEnvCPLQ's Appendix A.4 two-convex-edge envelope: its harmonic-
%                mean formula is PROVABLY always rank-1 PSD (b^2-4ac=0 identically, any slopes),
%                never strictly indefinite -- see conjPieceCPLQTest for the discriminant proof
%                and NOTES below on why the "raw indefinite bilinear, 2 convex edges" case is a
%                dead end for QuaPar (needs a genuine hyperbola, not a parabola).
%              - pure BILINEAR f = x*y with ZERO convex edges -> three-cone piecewise-linear
%                QuaPar (max of vertex linears, same construction as the affine case);
%              - pure BILINEAR f = x*y with exactly ONE convex edge -> six-face PARABOLIC QuaPar
%                (COAP B.2: edge region E, three L_W regions, two vertex cones; parabola arc).
%              den = 1 (no rational denominator). p may be a QuaPoly, QuaPar, or RatPol.
%              TODO: beta*x*y + linear (scale/shift), rational pieces (parabolic edges).
%
% NOTES (2026-07-07 investigation): the raw indefinite-bilinear "conjugate of x*y over a
% triangle with TWO convex edges" (as opposed to the rank-1 PSD quadratic above) requires, in
% general, comparing the two convex edges' own quadratic conjugate formulas against each other;
% that comparison boundary is PROVABLY a genuine hyperbola whenever the two edges' slopes differ
% (discriminant = (m1-m2)^2/(4 m1 m2) > 0), which QuaPar cannot represent (parabolic/linear edges
% only). This case does not actually arise in the wired pipeline: Step 1 (convEnvCPLQ) always
% convexifies a 2-convex-edge piece into the rank-1 PSD quadratic handled above before Step 2
% ever sees it, so conjPieceCPLQ is never asked to conjugate the raw indefinite piece directly in
% that case. See conjPieceCPLQTest/bilinearTwoConvexEdgesDocumentedLimitation.
% [output] g : QuaPar = the conjugate, a quadratic on a polyhedral subdivision with seven faces:
%              the central gradient image triangle T' = {A v_i + b}, three edge strips, and
%              three vertex cones; the conjugate is finite everywhere (domain = R^2).
%
% Since the gradient of a quadratic is affine, T' is a triangle and the subdivision is
% polyhedral (linear edges) -- no curved-edge handling is needed here. The conjugate of a
% RATIONAL envelope piece (whose gradient is nonlinear) yields genuinely PARABOLIC edges and is
% NOT implemented yet (it needs curved-edge orderEdges; see DESIGN.md II.5.1 and cPLQ Appendix B).
%
% Structure (vertices v1,v2,v3 CCW, edges e1=(v1,v2), e2=(v2,v3), e3=(v3,v1), outward normals
% n1,n2,n3):  s_i = A v_i + b;  on T'  f*(s) = 1/2 (s-b)' inv(A) (s-b) - c;  on the cone at v_i
% f*(s) = <s,v_i> - q(v_i);  on the strip of edge (v_i,v_j) with e = v_j - v_i,
% f*(s) = <s,v_i> - q(v_i) + (<s,e> - <grad q(v_i),e>)^2 / (2 e' A e).

    if ~(isa(p,'QuaPoly') || isa(p,'QuaPar') || isa(p,'RatPol'))
        error('conjPieceCPLQ:type','Input must be a QuaPoly/QuaPar/RatPol piece.');
    end
    if p.nf ~= 1 || p.nv ~= 3 || p.ne ~= 3 || ~p.isDomBounded
        error('conjPieceCPLQ:notImplemented', ...
            'Step 2 currently conjugates a single bounded-triangle piece (got nf=%d, nv=%d).', p.nf, p.nv);
    end
    if isa(p,'RatPol') && any(abs(p.den(1,1:2)) > sqrt(eps))
        error('conjPieceCPLQ:notImplemented', ...
            'Conjugate of a rational piece (-> parabolic) is not implemented yet (needs curved orderEdges).');
    end
    [L,Q,C] = QuaPoly.matrixForm(p.f(1,:));
    if ~isempty(C)
        error('conjPieceCPLQ:notImplemented','Quadratic numerator required (got a cubic).');
    end
    A = Q; b = L; c = p.f(1,end);
    V3 = p.V;
    if triSignedArea(V3) < 0, V3 = V3([1 3 2],:); end     % ensure CCW
    if norm(A) < sqrt(eps)
        g = conjLinearTriangle(V3, b, c);                 % affine piece (e.g. a concave envelope)
    elseif all(eig(A) > sqrt(eps))
        g = conjConvexQuadTriangle(V3, A, b, c);          % strictly convex quadratic
    elseif min(eig(A)) > -sqrt(eps) && max(eig(A)) > sqrt(eps)
        g = conjPSDRank1QuadTriangle(V3, A, b, c);        % convex, rank 1 (COAP A.4 envelope)
    elseif abs(A(1,1)) < sqrt(eps) && abs(A(2,2)) < sqrt(eps) && abs(A(1,2)-1) < sqrt(eps) ...
            && norm(b) < sqrt(eps) && abs(c) < sqrt(eps)
        % pure bilinear f = x*y (already in bilinear frame; general indefinite reduces to this
        % via bilinearFrame + an affine substitution -- to be added). Requires one convex edge.
        ce = convexEdgesXY(V3);
        if size(ce,1) == 0
            g = conjBilinearXYzeroCE(V3);                 % -> 3-cone linear QuaPar
        elseif size(ce,1) == 1
            g = conjBilinearXYoneCE(V3, ce);              % -> parabolic QuaPar (6 faces)
        else
            error('conjPieceCPLQ:notImplemented', ...
                'Bilinear conjugate is implemented for zero or one convex edge (got %d).', size(ce,1));
        end
    else
        error('conjPieceCPLQ:notImplemented', ...
            ['Only affine, positive-definite quadratic, or pure-bilinear (x*y) pieces over a ' ...
             'triangle are supported so far; general indefinite (rotation) and rational are next.']);
    end
end

% ============================================================================================
function g = conjLinearTriangle(V, gvec, h)
% Conjugate of an affine ell(x)=gvec'x+h over the CCW triangle V:
%   f*(s) = max_i [ <s,v_i> - ell(v_i) ],
% a piecewise-linear QuaPar with three unbounded cones meeting at the point s0 where the three
% vertex-linears L_i(s) = <s,v_i> - ell(v_i) are equal ([COAP] Appendix B.1).
    v1 = V(1,:)'; v2 = V(2,:)'; v3 = V(3,:)';
    lv = @(v) gvec'*v + h;
    g = conjVertexMax(V, [lv(v1); lv(v2); lv(v3)]);
end

function g = conjVertexMax(V, ellv)
% Conjugate f*(s) = max_i [ <s,v_i> - ellv(i) ] of arbitrary pointwise vertex values ellv (need
% not come from an affine function) over the CCW triangle V: the same 3-cone piecewise-linear
% QuaPar geometry as conjLinearTriangle, driven only by the three vertex values.
    v1 = V(1,:)'; v2 = V(2,:)'; v3 = V(3,:)';
    lv = @(i) ellv(i);
    s0 = [(v1-v2)'; (v1-v3)'] \ [lv(1)-lv(2); lv(1)-lv(3)];   % L1=L2=L3 here

    d12 = boundaryDir(v1-v2, v1-v3);   % along boundary L1=L2, oriented so L1>=L3
    d23 = boundaryDir(v2-v3, v2-v1);   % along boundary L2=L3, oriented so L2>=L1
    d13 = boundaryDir(v1-v3, v1-v2);   % along boundary L1=L3, oriented so L1>=L2

    Vd = [s0'; (s0+d12)'; (s0+d23)'; (s0+d13)'];
    E  = [1 2 0; 1 3 0; 1 4 0];        % rays from s0: r12, r23, r13
    Ec = zeros(3,6);
    adj = [1 2; 2 3; 1 3];             % r12:{cone1,cone2}, r23:{cone2,cone3}, r13:{cone1,cone3}
    d  = 0.5 * max([norm(d12), norm(d23), norm(d13)]);
    rep = [ (s0 + d*udir(d12) + d*udir(d13))';   % cone1 (between r12,r13)
            (s0 + d*udir(d12) + d*udir(d23))';   % cone2 (between r12,r23)
            (s0 + d*udir(d23) + d*udir(d13))' ];  % cone3 (between r23,r13)
    F = orientFaces(Vd, E, adj, rep);
    f = [0 0 0 v1(1) v1(2) -lv(1);
         0 0 0 v2(1) v2(2) -lv(2);
         0 0 0 v3(1) v3(2) -lv(3)];
    g = QuaPar(Vd, E, Ec, f, F);
end

function d = boundaryDir(w, u)
% Direction perpendicular to w (along the line w'*s = const), oriented so that <d,u> > 0.
    d = [-w(2); w(1)];
    if d'*u < 0, d = -d; end
end

function u = udir(d)
    u = d / norm(d);
end

% ============================================================================================
function g = conjConvexQuadTriangle(V, A, b, c)
% Build the 7-face QuaPar conjugate of q(x)=1/2 x'A x + b'x + c (A PD) over the CCW triangle V.
    v1 = V(1,:)'; v2 = V(2,:)'; v3 = V(3,:)';
    s1 = A*v1 + b; s2 = A*v2 + b; s3 = A*v3 + b;          % dual vertices = grad q(v_i)
    % edge vectors and unit OUTWARD normals (CCW triangle: outward = rotate edge by -90 deg)
    e1 = v2 - v1; e2 = v3 - v2; e3 = v1 - v3;
    n1 = unitNormal(e1); n2 = unitNormal(e2); n3 = unitNormal(e3);

    % vertices: 3 finite duals + 6 ray endpoints (s_i + unit normal direction)
    Vd = [s1'; s2'; s3'; (s1+n1)'; (s1+n3)'; (s2+n1)'; (s2+n2)'; (s3+n2)'; (s3+n3)'];

    % edges: 3 central segments, then 6 rays (isSegment = 0). Columns: [vtxA vtxB isSegment]
    E = [1 2 1;            % c1  s1-s2   (F0 | S1)
         2 3 1;            % c2  s2-s3   (F0 | S2)
         3 1 1;            % c3  s3-s1   (F0 | S3)
         1 4 0;            % r1  s1->s1+n1 (S1 | C1)
         1 5 0;            % r2  s1->s1+n3 (S3 | C1)
         2 6 0;            % r3  s2->s2+n1 (S1 | C2)
         2 7 0;            % r4  s2->s2+n2 (S2 | C2)
         3 8 0;            % r5  s3->s3+n2 (S2 | C3)
         3 9 0];           % r6  s3->s3+n3 (S3 | C3)
    Ec = zeros(9,6);        % all edges linear (polyhedral)

    % faces: 1=F0, 2=S1, 3=S2, 4=S3, 5=C1, 6=C2, 7=C3 ; adjacency (two faces per edge)
    adj = [1 2; 1 3; 1 4; 2 5; 4 5; 2 6; 3 6; 3 7; 4 7];

    % representative interior point of each face (for edge-orientation test)
    d = 0.5 * max([norm(s1-s2), norm(s2-s3), norm(s3-s1)]);   % offset scaled to dual size
    rep = [ ((s1+s2+s3)/3)';                 % F0 centroid
            ((s1+s2)/2 + d*n1)';             % S1
            ((s2+s3)/2 + d*n2)';             % S2
            ((s3+s1)/2 + d*n3)';             % S3
            (s1 + d*(n1+n3))';               % C1
            (s2 + d*(n1+n2))';               % C2
            (s3 + d*(n2+n3))' ];             % C3
    F = orientFaces(Vd, E, adj, rep);

    % functions (stored weighted basis [x^2 xy y^2 x y const]) per face
    M = inv(A);                                          %#ok<MINV>
    grad = -M*b;
    f0 = [M(1,1), M(1,2), M(2,2), grad(1), grad(2), 0.5*(b'*M*b) - c];
    f = zeros(7,6);
    f(1,:) = f0;
    f(2,:) = stripFun(v1, e1, A, b, c);                  % S1: edge e1 from v1
    f(3,:) = stripFun(v2, e2, A, b, c);                  % S2: edge e2 from v2
    f(4,:) = stripFun(v3, e3, A, b, c);                  % S3: edge e3 from v3
    f(5,:) = vertexFun(v1, A, b, c);                     % C1
    f(6,:) = vertexFun(v2, A, b, c);                     % C2
    f(7,:) = vertexFun(v3, A, b, c);                     % C3

    g = QuaPar(Vd, E, Ec, f, F);
end

function s6 = stripFun(vi, e, A, b, c)
% conjugate on the strip of an edge from vi with edge vector e:
%   <s,vi> - q(vi) + (<s,e> - g0)^2 / D,  g0 = <grad q(vi), e>,  D = 2 e'A e.
    D  = 2*(e'*A*e);
    g0 = (A*vi + b)'*e;
    qv = 0.5*(vi'*A*vi) + b'*vi + c;
    e1 = e(1); e2 = e(2);
    % plain coefficients a x^2 + b xy + c y^2 + d x + e y + f  (x=s1, y=s2)
    px2 = e1^2/D;  pxy = 2*e1*e2/D;  py2 = e2^2/D;
    pxc = vi(1) - 2*g0*e1/D;  pyc = vi(2) - 2*g0*e2/D;  pc = -qv + g0^2/D;
    s6 = [2*px2, pxy, 2*py2, pxc, pyc, pc];              % -> stored weighted basis
end

function s6 = vertexFun(vi, A, b, c)
% conjugate on the cone at vertex vi:  <s,vi> - q(vi)  (linear).
    qv = 0.5*(vi'*A*vi) + b'*vi + c;
    s6 = [0, 0, 0, vi(1), vi(2), -qv];
end

function n = unitNormal(e)
% unit outward normal of a CCW-triangle edge vector e: rotate e by -90 degrees, normalize.
    n = [e(2); -e(1)]; n = n / norm(n);
end

function F = orientFaces(V, E, adj, rep)
% Assign F=[leftFace rightFace] per edge using a representative interior point of each face.
    ne = size(E,1); F = zeros(ne,2);
    for j = 1:ne
        a = E(j,1); va = V(a,:); dir = V(E(j,2),:) - va;
        for w = 1:2
            k = adj(j,w);
            cr = dir(1)*(rep(k,2)-va(2)) - dir(2)*(rep(k,1)-va(1));  % (b-a) x (rep-a)
            if cr > 0, F(j,1) = k; else, F(j,2) = k; end
        end
    end
end

function a = triSignedArea(W)
    a = 0.5*((W(2,1)-W(1,1))*(W(3,2)-W(1,2)) - (W(3,1)-W(1,1))*(W(2,2)-W(1,2)));
end

% ============================================================================================
% ----- convex, RANK-1 (PSD, singular) quadratic conjugate over a triangle -------------------
function g = conjPSDRank1QuadTriangle(V, A, b, c)
% Conjugate of q(x)=1/2 x'Ax+b'x+c (A PSD, rank 1) over the CCW triangle V. This is the case
% produced by convEnvCPLQ's Appendix A.4 two-convex-edge envelope (always rank-1 PSD, see file
% header). Six-face QuaPar: three PARABOLIC edge strips (one per triangle edge; q is jointly
% convex so every edge -- not just "convex-for-xy" edges -- contributes an interior-critical-
% point region) and three LINEAR vertex cones, meeting at three COLLINEAR dual points
% s_i = A*v_i+b (collinear because A has rank 1: A*v = lambda*(u'v)*u varies only along u, the
% eigenvector of the nonzero eigenvalue lambda).
%
% Derivation: rotate to t=u'x (the curved direction) and r=uperp'x (the flat direction, where A
% has a zero eigenvalue). q(t,r) = (lambda/2)t^2 + p*t + w*r + c is jointly convex in (t,r), so
% for s off the measure-zero line {s'*uperp = w} the sup over the triangle is attained on its
% boundary: whichever of the triangle's two (t,r) boundary chains ("upper"/"lower") is active,
% selected by the sign of gamma = s'*uperp - w. Each chain segment (endpoints (ti,ri), slope m in
% (t,r)) contributes, where its own interior critical point falls within the segment's t-range, a
% quadratic value that is PROVABLY a parabola in s (its (alpha,gamma) quadratic part is a perfect
% square, same mechanism as the b^2-4ac=0 identity for the COAP B.2 bilinear edge formula).
% Elsewhere, the vertex conjugate f*(s)=<s,v_i>-q(v_i) applies. Validated by numeric sup-sampling
% over 200+ random (triangle, rank-1 PSD quadratic) pairs (see conjPieceCPLQTest).
    [Vec, Dg] = eig(A);
    [lam, idx] = max(diag(Dg));
    u = Vec(:,idx); uperp = [-u(2); u(1)];
    T = V*u; R = V*uperp;
    [Ts, ord] = sort(T); Rs = R(ord); Vs = V(ord,:);
    t1=Ts(1); t2=Ts(2); t3=Ts(3); r1=Rs(1); r2=Rs(2); r3=Rs(3);
    v1=Vs(1,:)'; v2=Vs(2,:)'; v3=Vs(3,:)';
    if min(t2-t1, t3-t2) < sqrt(eps)*(1+abs(t3-t1))
        error('conjPieceCPLQ:notImplemented', ...
            ['conjPSDRank1QuadTriangle: degenerate case (two vertices tie in the rotated t-' ...
             'coordinate, i.e. a triangle edge is exactly perpendicular to the quadratic''s ' ...
             'nonzero eigenvector) is not implemented.']);
    end
    s1 = A*v1+b; s2 = A*v2+b; s3 = A*v3+b;             % collinear dual "vertices"
    qv = @(v) 0.5*(v'*A*v) + b'*v + c;
    q1=qv(v1); q2=qv(v2); q3=qv(v3);

    m13 = (r3-r1)/(t3-t1); m12 = (r2-r1)/(t2-t1); m23 = (r3-r2)/(t3-t2);
    p = b'*u; w = b'*uperp;

    % which (t,r) chain the middle vertex sits on ("kink" chain = 2 segments; the other chain is
    % the direct edge v1-v3, one segment). Ray direction along {alpha+gamma*m = const} always has
    % gamma-rate +1 (d.uperp==1 identically); the sign below picks the ray that increases gamma
    % on the chain that actually needs gamma>0 (the "kink" chain iff it's the upper one).
    kinkHi = r2 > r1 + (r3-r1)*(t2-t1)/(t3-t1);
    if kinkHi, sgn13 = -1; sgn12 =  1; else, sgn13 =  1; sgn12 = -1; end
    sgn23 = sgn12;

    dray = @(m) [-(m*u(1)+u(2)); u(1)-m*u(2)];
    d13 = sgn13*dray(m13); d12 = sgn12*dray(m12); d23 = sgn23*dray(m23);

    E13 = rank1EdgeQuad(t1,r1,m13,u,lam,p,w,c);
    E12 = rank1EdgeQuad(t1,r1,m12,u,lam,p,w,c);
    E23 = rank1EdgeQuad(t2,r2,m23,u,lam,p,w,c);

    Vd = [s1'; s2'; s3'; (s1+d13)'; (s1+d12)'; (s2+d12)'; (s2+d23)'; (s3+d13)'; (s3+d23)'];
    E = [1 2 1; 2 3 1; 1 4 0; 1 5 0; 2 6 0; 2 7 0; 3 8 0; 3 9 0];
    Ec = zeros(8,6);                                   % every edge here is straight
    % faces: 1=E13(strip v1-v3), 2=E12(strip v1-v2), 3=E23(strip v2-v3), 4=L1, 5=L2, 6=L3
    adj = [1 2; 1 3; 1 4; 2 4; 2 5; 3 5; 1 6; 3 6];

    % Representative interior points, picked via the oracle (max of the 6 gated candidates)
    % rather than a hand-derived bisector sign: with two rays d_i,d_j from a shared apex, the
    % wedge between them is bisected by EITHER +(unit(d_i)+unit(d_j)) or its negation depending
    % on which of the two complementary angular sectors is the actual face, so try both and keep
    % whichever the oracle confirms (avoids a fragile by-hand orientation argument).
    winner = @(s) rank1Winner(s, u, uperp, lam, p, w, t1,t2,t3, r1,r2, m13,m12,m23, ...
        v1,v2,v3, q1,q2,q3, c);
    un = @(x) x/norm(x);
    scale = 0.4*norm(s3-s1) + 1;
    dirNonKink = uperp * (1 - 2*kinkHi);               % kinkHi -> -uperp (into gamma<0), else +uperp
    dirKink    = -dirNonKink;
    rep = zeros(6,2);
    rep(1,:) = pickRep((s1+s3)'/2, dirNonKink', 1, winner, scale);   % inside E13
    rep(2,:) = pickRep((s1+s2)'/2, dirKink', 2, winner, scale);      % inside E12
    rep(3,:) = pickRep((s2+s3)'/2, dirKink', 3, winner, scale);      % inside E23
    rep(4,:) = pickRep(s1', un(un(d13)+un(d12))', 4, winner, scale);   % inside L1 wedge
    rep(5,:) = pickRep(s2', un(un(d12)+un(d23))', 5, winner, scale);   % inside L2 wedge
    rep(6,:) = pickRep(s3', un(un(d13)+un(d23))', 6, winner, scale);   % inside L3 wedge
    F = orientFaces(Vd, E, adj, rep);

    f = [E13; E12; E23; ...
         0 0 0 v1(1) v1(2) -q1; ...
         0 0 0 v2(1) v2(2) -q2; ...
         0 0 0 v3(1) v3(2) -q3];
    g = QuaPar(Vd, E, Ec, f, F);
end

function s6 = rank1EdgeQuad(ti, ri, m, u, lam, p, w, c)
% Conjugate value contributed by the (t,r)-chain segment through (ti,ri) with slope m: the
% interior-critical-point value (alpha+gamma*m)^2/(2*lam) + gamma*(ri-m*ti) - c, expressed as a
% plain-then-stored quadratic in s (alpha=s.u-p, gamma=s.uperp-w are affine-orthogonal in s, so
% substituting keeps the quadratic part a perfect square: b^2-4ac=0, i.e. a genuine parabola).
    u1=u(1); u2=u(2);
    k1 = u1-m*u2; k2 = u2+m*u1; k0 = -(p+m*w);
    Ap = k1^2/(2*lam); Bp = k1*k2/lam; Cp = k2^2/(2*lam);
    Dp = k1*k0/lam - u2*(ri-m*ti);
    Ep = k2*k0/lam + u1*(ri-m*ti);
    Fp = k0^2/(2*lam) - w*(ri-m*ti) - c;
    s6 = [2*Ap, Bp, 2*Cp, Dp, Ep, Fp];                 % plain -> stored weighted basis
end

function p = pickRep(base, dir, wantIdx, winner, scale)
% Try base +/- dir*mag for a geometric sequence of magnitudes (dir need not be unit length): the
% correct sign (which of the two complementary sides a shared boundary/apex is on) and the
% correct magnitude (small enough not to overshoot a bounded strip's valid range, large enough to
% clear a wedge's apex) are both easier to search for than to derive by hand.
    for mag = scale * [1, 0.3, 0.1, 0.03, 0.01, 0.003, 0.001, 3, 10]
        for sgn = [1, -1]
            pt = base + sgn*mag*dir;
            if winner(pt) == wantIdx, p = pt; return; end
        end
    end
    error('conjPieceCPLQ:internal', ...
        'conjPSDRank1QuadTriangle: could not locate a representative point for face %d.', wantIdx);
end

function idx = rank1Winner(s, u, uperp, lam, p, w, t1, t2, t3, r1, r2, m13, m12, m23, ...
        v1, v2, v3, q1, q2, q3, c)
% Oracle: index (1=E13,2=E12,3=E23,4=L1,5=L2,6=L3) of the true argmax of the six rank-1-quadratic
% conjugate candidates at dual point s (row vector), gating each edge candidate by its own
% interior-critical-point t*-range (see conjPSDRank1QuadTriangle).
    s = s(:);
    alpha = s'*u - p; gamma = s'*uperp - w;
    tol = 1e-7*(1+abs(t3-t1));
    evalStored = @(row) 0.5*row(1)*s(1)^2 + row(2)*s(1)*s(2) + 0.5*row(3)*s(2)^2 ...
        + row(4)*s(1) + row(5)*s(2) + row(6);
    vals = -inf(6,1);
    vals(4) = s'*v1 - q1; vals(5) = s'*v2 - q2; vals(6) = s'*v3 - q3;
    xs13 = (alpha+gamma*m13)/lam;
    if xs13 >= t1-tol && xs13 <= t3+tol
        vals(1) = evalStored(rank1EdgeQuad(t1,r1,m13,u,lam,p,w,c));
    end
    xs12 = (alpha+gamma*m12)/lam;
    if xs12 >= t1-tol && xs12 <= t2+tol
        vals(2) = evalStored(rank1EdgeQuad(t1,r1,m12,u,lam,p,w,c));
    end
    xs23 = (alpha+gamma*m23)/lam;
    if xs23 >= t2-tol && xs23 <= t3+tol
        vals(3) = evalStored(rank1EdgeQuad(t2,r2,m23,u,lam,p,w,c));
    end
    [~, idx] = max(vals);
end

% ----- bilinear (f=x*y) conjugate over a triangle with ONE convex edge (COAP B.2) -----------
function ce = convexEdgesXY(V)
% Convex edges of the triangle for f = x*y: positive-slope edges. Returns rows [i j m q].
    tol = sqrt(eps); ed = [1 2; 2 3; 3 1]; ce = zeros(0,4);
    for t = 1:3
        i = ed(t,1); j = ed(t,2); vi = V(i,:); vj = V(j,:); dx = vj(1)-vi(1);
        if abs(dx) < tol, continue; end
        m = (vj(2)-vi(2))/dx;
        if m > tol, ce(end+1,:) = [i, j, m, vi(2)-m*vi(1)]; end %#ok<AGROW>
    end
end

function g = conjBilinearXYzeroCE(V)
% Conjugate of f = x*y over a triangle with NO convex edge: every edge is concave for x*y, so the
% restriction of <s,x>-xy to each edge is concave and its max is at an endpoint; hence the overall
% sup over the triangle is attained at a vertex and f*(s) = max_i(<s,v_i> - xy(v_i)), the same
% 3-cone piecewise-linear QuaPar as conjLinearTriangle (COAP B.1 applied to the vertex values).
    xy = @(v) v(1)*v(2);
    ellv = [xy(V(1,:)); xy(V(2,:)); xy(V(3,:))];
    g = conjVertexMax(V, ellv);
end

function g = conjBilinearXYoneCE(V, ce)
% Conjugate of f = x*y over a triangle with one convex edge (endpoints A,B; special vertex W).
% Six convex faces: SE (parabolic quad E), CW_mid/CW_below/CW_above (linear L_W), C_A (L_A),
% C_B (L_B); seven edges: the parabola arc P0-P1 plus six rays. Verified vs numeric sup.
    m = ce(3); q = ce(4);
    A = V(ce(1),:)'; B = V(ce(2),:)';
    W = V(setdiff(1:3, [ce(1) ce(2)]), :)';
    if A(1) > B(1), tmp = A; A = B; B = tmp; end   % canonical order x_A < x_B (m>0 => distinct)
    xy = @(v) v(1)*v(2);
    % parabola (crs) and the region functions (plain coefficients)
    crs = [-1, -2*m, -m^2, 2*q+4*m*W(1), -(2*m*q-4*m*W(2)), -(q^2+4*m*W(1)*W(2))];
    E  = [1/(4*m), 1/2, m/4, -q/(2*m), q/2, q^2/(4*m)];
    LA = [0 0 0 A(1) A(2) -xy(A)]; LB = [0 0 0 B(1) B(2) -xy(B)]; LW = [0 0 0 W(1) W(2) -xy(W)];
    % strip lines x* = x_A, x_B : s1 + m s2 - q - 2 m x = 0
    P0 = lineParabolaInt([1 m (-q-2*m*A(1))], crs);
    P1 = lineParabolaInt([1 m (-q-2*m*B(1))], crs);
    gc = @(s) [2*crs(1)*s(1)+crs(2)*s(2)+crs(4); crs(2)*s(1)+2*crs(3)*s(2)+crs(5)];
    gw = [1; m]; rp = @(w) [-w(2); w(1)]; un = @(d) d/norm(d);
    dSEA = [m;-1]; if gc(P0)'*dSEA > 0, dSEA = -dSEA; end   % crs decreasing along strip line
    dSEB = [m;-1]; if gc(P1)'*dSEB > 0, dSEB = -dSEB; end
    dCAW = rp(A-W);  if gw'*dCAW > 0, dCAW = -dCAW; end     % x* decreasing (C_A side)
    dCBW = rp(B-W);  if gw'*dCBW < 0, dCBW = -dCBW; end     % x* increasing (C_B side)
    dAL = -dSEA; dBL = -dSEB;                               % lower strip rays
    Vd = [P0'; P1'; (P0+un(dSEA))'; (P0+un(dAL))'; (P0+un(dCAW))'; ...
          (P1+un(dSEB))'; (P1+un(dBL))'; (P1+un(dCBW))'];
    E7 = [1 2 1; 1 3 0; 1 4 0; 1 5 0; 2 6 0; 2 7 0; 2 8 0];
    Ec = zeros(7,6); Ec(1,:) = crs;
    % faces: 1=SE, 2=CW_mid, 3=C_A, 4=CW_below, 5=C_B, 6=CW_above
    adj = [1 2; 1 3; 2 4; 3 4; 1 5; 2 6; 5 6];
    sc = 0.6*norm(P1-P0) + 1;
    rep = [ ((P0+P1)/2 + sc*un(dSEA))';                 % SE
            ((P0+P1)/2 + sc*un(-dSEA))';                % CW_mid (parabola interior)
            (P0 + sc*(un(dSEA)+un(dCAW)))';             % C_A
            (P0 + sc*(un(dAL)+un(dCAW)))';              % CW_below
            (P1 + sc*(un(dSEB)+un(dCBW)))';             % C_B
            (P1 + sc*(un(dBL)+un(dCBW)))' ];            % CW_above
    F = orientFaces(Vd, E7, adj, rep);
    if evalConicPt(Ec(1,:), rep(F(1,1),:)) < 0, Ec(1,:) = -Ec(1,:); end   % >0 on left face
    toS = @(p) [2*p(1), p(2), 2*p(3), p(4), p(5), p(6)];
    f = [toS(E); toS(LW); toS(LA); toS(LW); toS(LB); toS(LW)];
    g = QuaPar(Vd, E7, Ec, f, F);
end

function P = lineParabolaInt(l, c)
% Intersection of line l=[1 b c0] (s1 = -c0 - b s2) with conic c. For our strip lines (parallel
% to the parabola axis) this is a single point; solve the induced polynomial in s2.
    b = l(2); c0 = l(3); r = -b; p = -c0;               % s1 = p + r s2  (l(1)=1)
    A2 = c(1)*r^2 + c(2)*r + c(3);
    A1 = c(1)*2*p*r + c(2)*p + c(4)*r + c(5);
    A0 = c(1)*p^2 + c(4)*p + c(6);
    if abs(A2) < 1e-10
        s2 = -A0/A1;
    else
        rr = roots([A2 A1 A0]); rr = rr(abs(imag(rr)) < 1e-9); s2 = real(rr(1));
    end
    P = [p + r*s2; s2];
end

function z = evalConicPt(c, s)
    z = c(1)*s(1)^2 + c(2)*s(1)*s(2) + c(3)*s(2)^2 + c(4)*s(1) + c(5)*s(2) + c(6);
end
