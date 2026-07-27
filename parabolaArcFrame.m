function fr = parabolaArcFrame(Ec, errId)
% parabolaArcFrame  Intrinsic (u,v) frame of a PARABOLA: the shared substrate for every
%   parabola-arc computation in the 'cplq' pipeline (clipArcByHalfPlane.m, and maxQuaPar.m's
%   curved-edge clipping). Factored out of clipArcByHalfPlane.m -- which was its sole original
%   home -- when maxQuaPar.m needed the same frame for arc midpoints/tangents and for its own
%   "does the arc bulge across this clip line" test, so the two cannot drift apart.
%
% [input]  Ec    : 1x6 conic [a b c d e f] for a*x^2+b*xy+c*y^2+d*x+e*y+f=0, a genuine PARABOLA
%                  (b^2-4ac==0, checked).
%          errId : optional identifier COMPONENT for the two errors raised here, so each caller
%                  keeps its OWN public error identifiers (e.g. clipArcByHalfPlane:notParabola)
%                  instead of leaking this helper's name into them. Defaults to the file's name.
% [output] fr    : struct with
%            uDir, nullDir : orthonormal 1x2 rows. nullDir is the parabola's own AXIS direction
%                  (the zero-eigenvalue eigenvector of the quadratic part Q=[a b/2;b/2 c], which
%                  has rank exactly 1 for a genuine parabola); uDir is its perpendicular, hence
%                  guaranteed TRANSVERSE to the axis. That is what makes u a GLOBAL, monotonic
%                  parameter for every point of the conic (not just for one arc), so ordering two
%                  points by u correctly orders them along the finite arc between them.
%            lam, du, dv, f0 : this frame's own coefficients; v = q(u) = -(lam*u^2+du*u+f0)/dv.
%            point(u)   : 1x2 point of the parabola at parameter u.
%            tangent(u) : 1x2 derivative d/du of point(u) (never zero: its uDir component is 1).
%            uOf(X)     : parameter u of a point X (1x2 row) lying on the parabola.
%            lineCoeffs(nrm,c) : [A2 A1 A0] with nrm*point(u)'-c == A2*u^2+A1*u+A0 -- the
%                  restriction of a half-plane's boundary function to the conic, an explicit
%                  quadratic in u. A line/arc crossing is therefore solved directly, with no
%                  line-conic-intersection-then-reproject detour, and with no ambiguity about
%                  which of the (up to 2) roots lies on a GIVEN arc: that is decided by u-range
%                  membership instead.
    if nargin < 2, errId = 'parabolaArcFrame'; end
    Q = [Ec(1), Ec(2)/2; Ec(2)/2, Ec(3)];
    delta = Ec(2)^2 - 4*Ec(1)*Ec(3);
    if abs(delta) > 1e-6*(norm(Q,'fro')+1)^2
        error([errId ':notParabola'], 'Ec is not a parabola (b^2-4ac=%g).', delta);
    end
    [Vq, D] = eig(Q);
    [~, i0] = min(abs(diag(D)));
    iu = 3 - i0;
    nullDir = Vq(:,i0)'/norm(Vq(:,i0));
    uDir = Vq(:,iu)'/norm(Vq(:,iu));
    lam = uDir*Q*uDir'; dv = [Ec(4) Ec(5)]*nullDir'; du = [Ec(4) Ec(5)]*uDir'; f0 = Ec(6);
    if abs(dv) < 1e-10*(abs(lam)+abs(du)+1)
        error([errId ':degenerateAxis'], ...
            'Conic does not depend on v (dv~0): not a genuine parabola in this frame.');
    end
    fr.uDir = uDir; fr.nullDir = nullDir;
    fr.lam = lam; fr.du = du; fr.dv = dv; fr.f0 = f0;
    fr.point      = @(u) u*uDir + (-(lam*u.^2+du*u+f0)/dv)*nullDir;
    fr.tangent    = @(u) uDir + (-(2*lam*u+du)/dv)*nullDir;
    fr.uOf        = @(X) X*uDir';
    fr.lineCoeffs = @(nrm,c) [ -(nrm*nullDir')*lam/dv, ...
                                nrm*uDir' - (nrm*nullDir')*du/dv, ...
                               -(nrm*nullDir')*f0/dv - c ];
end
