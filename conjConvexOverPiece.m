function cj = conjConvexOverPiece(r, Q, L, c, dualVars)
% conjConvexOverPiece  The conjugate of a CONVEX quadratic over one polyhedral piece.
%
% [input]  r        : region with affine facets -- one of the three shapes plq_1p produces,
%                     i.e. a bounded TRIANGLE, a WEDGE (one vertex, two rays) or a HALF-STRIP
%                     (two vertices, two parallel rays). Deliberately not "any polyhedron":
%                     these are the only shapes fanUnboundedFace emits, and QuaPol's own faces
%                     reduce to them.
%          Q, L, c  : q(x) = 1/2 x'Qx + L'x + c, with Q positive SEMIdefinite.
%          dualVars : 1x2 sym, [s_1 s_2].
% [output] cj       : functionNDomain array partitioning dom q* into vertex, edge and interior
%                     cells.
%
% WHY THIS IS NEEDED AT ALL. When q is convex, co(q|P) = q -- there is no Step 1 work to do --
% but Step 2 then has to conjugate a CURVED function, which cPLQ's nCE==0 branch (a support
% function, valid only for an affine envelope) cannot. Measured: conjugateOfPiecePoly on a
% BOUNDED triangle with strictly convex q returns only the 3 vertex pieces, leaving every
% interior dual point uncovered. That is by design for cPLQ, whose Step 1 always hands Step 2 an
% affine or rank-1-PSD envelope, so the edge and interior cells never arise there. They do arise
% here, because for a convex face the envelope IS q.
%
% THE DECOMPOSITION is the KKT active set of the concave program
%       q*(s) = max { <s,x> - q(x) : x in P },
% whose objective is concave, so the maximizer is unique-ish and lies on exactly one relatively
% open face of P. Which face is decided by where s - grad q(x*) sits:
%
%   VERTEX v      x* = v, valid iff s - grad q(v) lies in the normal cone N_P(v), which for a
%                 polyhedron is {u : <u,e> <= 0 for every edge direction e leaving v}. A RAY
%                 leaving v contributes its direction exactly as a bounded edge does, so
%                 unbounded pieces need no special case. Value: <s,v> - q(v), affine in s.
%
%   EDGE (facet with outward normal n, direction d, base x0)
%                 x* = x0 + t* d with t* from d/dt[<s,x> - q(x)] = 0, i.e.
%                     t* = <s - grad q(x0), d> / (d'Qd),
%                 which needs d'Qd > 0; when d'Qd = 0 the objective is affine along the edge and
%                 its max sits at an endpoint, so the edge contributes no cell. Valid iff t* lies
%                 within the edge's extent AND s - grad q(x*) = mu*n with mu >= 0. Both are
%                 affine in s (t* is), so the cell is polyhedral; the value is quadratic in s.
%
%   INTERIOR      x* = Q^{-1}(s - L), needing Q nonsingular (Q PSD and nonsingular = PD). Valid
%                 iff x* lies in P. Value 1/2 (s-L)'Q^{-1}(s-L) - c, the familiar
%                 [GARDINER-13] Fact 3 formula with the inverse in place of the pseudo-inverse.
%                 A singular PSD Q makes the unconstrained sup finite only on the measure-zero
%                 set s - L in range(Q), which carries no 2-D cell, so it is correctly skipped.
%
% Cells are emitted for whichever of the above are nonempty; region() itself returns empty for
% an infeasible constraint set, so an inactive cell drops out rather than needing a test here.

    s1 = dualVars(1); s2 = dualVars(2);
    sv = [s1; s2];
    % NOTE: every transpose of a SYMBOLIC vector below is .' (transpose), never '
    % (ctranspose). MATLAB's ' conjugates, and s_1/s_2 are not declared real, so
    % sv'*x produced literal conj(s_1) terms and simplify then folded s_1*conj(s_1)
    % into abs(s_1)^2. Numerically identical on real input -- which is why the value
    % checks still passed -- but the EXPRESSIONS go on to build region inequalities and
    % be tested with isAlways, where a stray conj/abs is not harmless.
    cj = functionNDomain.empty();

    Q = double(Q); L = double(L(:)); c = double(c);
    tolQ = 1e-9 * max(1, max(abs(Q(:))));
    if any(eig((Q+Q')/2) < -tolQ)
        error('conjConvexOverPiece:notConvex', ...
            'q must be convex here (Q positive semidefinite); this Q is indefinite or negative.');
    end
    Q = (Q + Q')/2;

    [A, b, lin] = r.linearForm;
    if ~all(lin)
        error('conjConvexOverPiece:nonAffineFacet', ...
            'every facet must be affine to read the active set.');
    end
    [~, px, py] = r.finiteVertices;
    V = uniqueRows(double([px(:), py(:)]));

    scale = max(1, max(abs([A(:); b(:)])));
    tolA  = 1e-7 * scale;
    inP   = @(z) all(A*z(:) <= b + tolA);
    gradq = @(z) Q*z(:) + L;
    qAt   = @(z) 0.5*(z(:).'*Q*z(:)) + L.'*z(:) + c;

    % ---- vertex cells ---------------------------------------------------------------------
    for k = 1:size(V,1)
        v = V(k,:)';
        dirs = edgeDirsAt(A, b, v, tolA);
        if isempty(dirs), continue, end
        g = sym.empty(1,0);
        gv = gradq(v);
        for t = 1:size(dirs,1)
            g(t) = dirs(t,1)*(s1 - gv(1)) + dirs(t,2)*(s2 - gv(2));
        end
        rk = region(g, dualVars);
        if isempty(rk), continue, end
        fk = symbolicFunction(simplify(s1*v(1) + s2*v(2) - qAt(v)));
        cj = [cj, functionNDomain(fk, rk)]; %#ok<AGROW>
    end

    % ---- edge cells -----------------------------------------------------------------------
    for i = 1:size(A,1)
        n = A(i,:)';
        if norm(n) <= tolA, continue, end
        d = [-n(2); n(1)];  d = d / norm(d);
        dQd = d'*Q*d;
        if dQd <= tolQ
            continue        % q is affine along this edge: its max is at an endpoint, no cell
        end
        onIdx = [];
        for k = 1:size(V,1)
            if abs(A(i,:)*V(k,:)' - b(i)) <= tolA, onIdx(end+1) = k; end %#ok<AGROW>
        end
        if isempty(onIdx), continue, end
        x0 = V(onIdx(1),:)';
        [tmin, tmax] = edgeExtent(A, b, V, onIdx, x0, d, tolA);

        g0 = gradq(x0);
        tstar = ( d(1)*(s1 - g0(1)) + d(2)*(s2 - g0(2)) ) / dQd;   % affine in s
        xstar = x0 + d*tstar;                                       % affine in s
        gx    = Q*xstar + L;
        mu    = n.'*(sv - gx) / (n.'*n);                            % affine in s

        g = sym.empty(1,0); ng = 0;
        ng = ng+1; g(ng) = -mu;                                     % mu >= 0
        if isfinite(tmin), ng = ng+1; g(ng) = tmin - tstar; end
        if isfinite(tmax), ng = ng+1; g(ng) = tstar - tmax; end
        rk = region(expand(g), dualVars);
        if isempty(rk), continue, end
        val = sv.'*xstar - (0.5*(xstar.'*Q*xstar) + L.'*xstar + c);
        cj = [cj, functionNDomain(symbolicFunction(simplify(expand(val))), rk)]; %#ok<AGROW>
    end

    % ---- interior cell --------------------------------------------------------------------
    if rcond(Q) > 1e-12
        Qi = inv(Q);
        xstar = Qi*(sv - L);
        g = sym.empty(1,0);
        for i = 1:size(A,1)
            g(i) = expand(A(i,:)*xstar - b(i));
        end
        rk = region(g, dualVars);
        if ~isempty(rk)
            val = 0.5*((sv - L).'*Qi*(sv - L)) - c;
            cj = [cj, functionNDomain(symbolicFunction(simplify(expand(val))), rk)];
        end
    end

    if isempty(cj)
        error('conjConvexOverPiece:noCell', ...
            'no dual cell was produced; the piece or the quadratic is degenerate beyond what this handles.');
    end
end

% ------------------------------------------------------------------------------------------
function dirs = edgeDirsAt(A, b, v, tolA)
% Directions of the edges of {A x <= b} leaving vertex v: for each ACTIVE facet, whichever of
% the two directions along it stays feasible. A ray is found the same way a bounded edge is.
    dirs = zeros(0,2);
    for i = 1:size(A,1)
        if abs(A(i,:)*v - b(i)) > tolA, continue, end
        n = A(i,:);
        if norm(n) <= 0, continue, end
        e = [-n(2), n(1)] / norm(n);
        for sg = [1 -1]
            if all(A*(v + 1e-6*sg*e') <= b + tolA)
                dirs(end+1,:) = sg*e; %#ok<AGROW>
            end
        end
    end
    dirs = uniqueRows(dirs);
end

function [tmin, tmax] = edgeExtent(A, b, V, onIdx, x0, d, tolA)
% The t-range of the facet's edge in the parametrization x0 + t*d, with -inf/+inf where the edge
% runs off to infinity (which is exactly when +-d is a recession direction of the region).
    ts = zeros(1, numel(onIdx));
    for k = 1:numel(onIdx)
        ts(k) = (V(onIdx(k),:)' - x0)' * d;
    end
    tmin = min(ts); tmax = max(ts);
    big = 1e8;
    if all(A*(x0 + big*d) <= b + tolA*big), tmax = inf; end
    if all(A*(x0 - big*d) <= b + tolA*big), tmin = -inf; end
end

function V = uniqueRows(V)
    keep = true(size(V,1),1);
    for i = 2:size(V,1)
        for j = 1:i-1
            if keep(j) && norm(V(i,:) - V(j,:)) < 1e-9, keep(i) = false; break, end
        end
    end
    V = V(keep,:);
end
