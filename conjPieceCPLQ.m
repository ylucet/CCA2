function g = conjPieceCPLQ(p)
% conjPieceCPLQ  Step 2 of the 'cplq' pipeline: Fenchel conjugate of a single envelope piece.
%
% objective: conjugate (p + I_P)*(s) = sup_{x in P} <s,x> - p(x) of one convex-envelope piece,
%   returned as a QuaPar (quadratic on a parabolic subdivision; here polyhedral). Step 2 of
%     [COAP] Karmarkar & Lucet, Comput. Optim. Appl. 94 (2026) 747-780 (Appendix B);
%     [JOGO] Karmarkar & Lucet, J. Glob. Optim. 94 (2026) 3-34.
%
% [input]  p : a single convex piece over a bounded triangle. Currently implemented:
%              an AFFINE piece (-> f* = max_i(<s,v_i>-ell(v_i)), three cones) and a CONVEX
%              (positive-definite) QUADRATIC (-> seven regions). den = 1 (no rational
%              denominator). p may be a QuaPoly, QuaPar, or RatPol.
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
    else
        error('conjPieceCPLQ:notImplemented', ...
            ['Only an affine or a positive-definite quadratic piece over a triangle is supported; ' ...
             'got a non-PD, non-zero Hessian (indefinite or rank-deficient, or a rational piece).']);
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
    s0 = [(v1-v2)'; (v1-v3)'] \ [lv(v1)-lv(v2); lv(v1)-lv(v3)];   % L1=L2=L3 here

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
    f = [0 0 0 v1(1) v1(2) -lv(v1);
         0 0 0 v2(1) v2(2) -lv(v2);
         0 0 0 v3(1) v3(2) -lv(v3)];
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
