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

    % EXACT VALUES, NUMERIC DECISIONS. This routine used to open with
    %     Q = double(Q); L = double(L(:)); c = double(c);   and   V = double([px, py])
    % and build every cell out of those. That is a defect, and a subtle one: the numbers are
    % right to 16 digits, every value check passes, and yet two cells that SHARE A FACET come out
    % carrying two different doubles of the same exact number, because each rounds it by its own
    % route. MEASURED 2026-08-17 on x*y over conv{(0,0),(3,3),(1,2)}: two adjacent cells carried
    %     s_2 - 659536895553805/562949953421312   and   s_2 - 5276295164430439/4503599627370496
    % -- both are `4 - 2*sqrt(2)`, ONE ULP APART. region.merge finds a shared facet by asking
    % whether one constraint is the negation of another, and no comparison, structural or
    % symbolic, can identify those two. So the facet is invisible, the cells never merge, and
    % Step 3's cross-piece maximum grows without bound: 57 cells carrying 10 distinct functions,
    % 4 merges out of 612 attempts. DECISIONS.md 2026-08-17 has the chain.
    %
    % The split below is the fix and the rule: anything that becomes part of an EXPRESSION is
    % built from the exact `*S` quantities; the doubles decide only COMBINATORICS -- which facets
    % are active at a vertex, which side of a probe a point falls, whether an edge runs to
    % infinity, whether Q is invertible. Those are decisions, not values, and they are where this
    % file's tolerances belong.
    QS = sym(Q); QS = (QS + QS.')/2;
    LS = sym(L(:));
    cS = sym(c);
    Q = double(QS); L = double(LS); c = double(cS);
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
    % The same rows, exactly: constraint j is ineqs(j) = c0 + c1*x + c2*y <= 0, so
    % A(j,:) = [c1 c2] and b(j) = -c0 -- read by EVALUATION, the convention region.linearForm
    % states and for the reason it gives (it does not care how the expression is written).
    AS = sym(zeros(size(A))); bS = sym(zeros(size(b)));
    for j = 1:size(A,1)
        gj = r.ineqs(j).f;
        e0 = subs(gj, r.vars, [0 0]);
        AS(j,:) = [subs(gj, r.vars, [1 0]) - e0, subs(gj, r.vars, [0 1]) - e0];
        bS(j)   = -e0;
    end

    [~, px, py] = r.finiteVertices;
    VS = [px(:), py(:)];
    [V, keepV] = uniqueRows(double(VS));
    VS = VS(keepV,:);

    scale = max(1, max(abs([A(:); b(:)])));
    tolA  = 1e-7 * scale;
    gradqS = @(z) QS*z(:) + LS;
    qAtS   = @(z) 0.5*(z(:).'*QS*z(:)) + LS.'*z(:) + cS;

    % ---- vertex cells ---------------------------------------------------------------------
    for k = 1:size(V,1)
        v = V(k,:)';
        vS = VS(k,:).';
        [dirIdx, dirSgn] = edgeDirsAt(A, b, v, tolA);
        if isempty(dirIdx), continue, end
        g = sym.empty(1,0);
        gv = gradqS(vS);
        for t = 1:numel(dirIdx)
            % The exact edge direction, deliberately NOT normalized: normalizing would introduce
            % a square root for nothing, and a cone constraint is unchanged by a POSITIVE scaling
            % of its direction.
            e = dirSgn(t) * [-AS(dirIdx(t),2); AS(dirIdx(t),1)];
            g(t) = e(1)*(s1 - gv(1)) + e(2)*(s2 - gv(2));
        end
        rk = region(g, dualVars);
        if isempty(rk), continue, end
        fk = symbolicFunction(simplify(s1*vS(1) + s2*vS(2) - qAtS(vS)));
        cj = [cj, functionNDomain(fk, rk)]; %#ok<AGROW>
    end

    % ---- edge cells -----------------------------------------------------------------------
    for i = 1:size(A,1)
        n = A(i,:)';
        if norm(n) <= tolA, continue, end
        nS = AS(i,:).';
        dS = [-nS(2); nS(1)];               % exact, unnormalized (see below)
        d  = [-n(2); n(1)];                 % the same vector, numerically
        du = d / norm(d);                   % unit, for the run-to-infinity probe only
        dQd = dS.'*QS*dS;
        % `tstar` scales like 1/alpha under d -> alpha*d and `xstar = x0 + dS*tstar` is therefore
        % INVARIANT, so leaving dS unnormalized changes no cell -- it only keeps the square root
        % out of every expression below. The extents tS are put in the same parametrization.
        if double(dQd) <= tolQ * max(1, norm(d)^2)
            continue        % q is affine along this edge: its max is at an endpoint, no cell
        end
        onIdx = [];
        for k = 1:size(V,1)
            if abs(A(i,:)*V(k,:)' - b(i)) <= tolA, onIdx(end+1) = k; end %#ok<AGROW>
        end
        if isempty(onIdx), continue, end
        x0 = V(onIdx(1),:)';
        x0S = VS(onIdx(1),:).';

        tS = sym(zeros(1, numel(onIdx)));
        for k = 1:numel(onIdx)
            tS(k) = ((VS(onIdx(k),:).' - x0S).' * dS) / (dS.'*dS);
        end
        tn = double(tS);
        [~, kmin] = min(tn); [~, kmax] = max(tn);
        tminS = tS(kmin); tmaxS = tS(kmax);
        big = 1e8;
        runsUp   = all(A*(x0 + big*du) <= b + tolA*big);
        runsDown = all(A*(x0 - big*du) <= b + tolA*big);

        g0 = gradqS(x0S);
        tstar = ( dS(1)*(s1 - g0(1)) + dS(2)*(s2 - g0(2)) ) / dQd;   % affine in s
        xstar = x0S + dS*tstar;                                      % affine in s
        gx    = QS*xstar + LS;
        mu    = nS.'*(sv - gx) / (nS.'*nS);                          % affine in s

        g = sym.empty(1,0); ng = 0;
        ng = ng+1; g(ng) = -mu;                                      % mu >= 0
        if ~runsDown, ng = ng+1; g(ng) = tminS - tstar; end
        if ~runsUp,   ng = ng+1; g(ng) = tstar - tmaxS; end
        rk = region(expand(g), dualVars);
        if isempty(rk), continue, end
        val = sv.'*xstar - (0.5*(xstar.'*QS*xstar) + LS.'*xstar + cS);
        cj = [cj, functionNDomain(symbolicFunction(simplify(expand(val))), rk)]; %#ok<AGROW>
    end

    % ---- interior cell --------------------------------------------------------------------
    if rcond(Q) > 1e-12
        QiS = inv(QS);
        xstar = QiS*(sv - LS);
        g = sym.empty(1,0);
        for i = 1:size(A,1)
            g(i) = expand(AS(i,:)*xstar - bS(i));
        end
        rk = region(g, dualVars);
        if ~isempty(rk)
            val = 0.5*((sv - LS).'*QiS*(sv - LS)) - cS;
            cj = [cj, functionNDomain(symbolicFunction(simplify(expand(val))), rk)];
        end
    end

    if isempty(cj)
        error('conjConvexOverPiece:noCell', ...
            'no dual cell was produced; the piece or the quadratic is degenerate beyond what this handles.');
    end
end

% ------------------------------------------------------------------------------------------
function [idx, sgn] = edgeDirsAt(A, b, v, tolA)
% Directions of the edges of {A x <= b} leaving vertex v: for each ACTIVE facet, whichever of
% the two directions along it stays feasible. A ray is found the same way a bounded edge is.
%
% Returns the facet INDEX and the SIGN rather than the direction itself, so the caller can
% rebuild the direction EXACTLY from its symbolic constraint row. The feasibility probe and the
% deduplication stay numeric -- they are decisions, and this file's tolerances belong to them.
%
% HISTORY (vertex-cone coverage gap, 2026-08-30): the probe step used to be the fixed constant
% 1e-6, independent of tolA. At a vertex where the OTHER active facet's cross-sensitivity to a
% move along e is small, the resulting constraint change can come out SMALLER than tolA, so a
% genuinely infeasible direction (moving backward along an edge) spuriously passes `<= b+tolA` --
% both signs of the same facet then look feasible, producing a duplicate direction that reduces
% the vertex's cell from a 2D cone to a single point. Measured on T1=(0,0),(60,10),(15,10),
% q=1/2 x'diag stuff -- Q=[3 -1;-1 5], at v=(0,0): step 1e-6 against tolA=1e-6 left the
% cross-facet violation (4.9e-7) below tolA, so BOTH signs along the (60,10)-facing edge passed,
% and the resulting cell excluded s=(-1,-1) even though vertex (0,0) is provably optimal there
% (independent brute-force oracle: f*=0). The step must be large relative to tolA, not an
% independent constant, since both are decisions about the SAME comparison.
    idx = []; sgn = [];
    seen = zeros(0,2);
    step = max(1e-6, 1e4*tolA);
    for i = 1:size(A,1)
        if abs(A(i,:)*v - b(i)) > tolA, continue, end
        n = A(i,:);
        if norm(n) <= 0, continue, end
        e = [-n(2), n(1)] / norm(n);
        for sg = [1 -1]
            if all(A*(v + step*sg*e') <= b + tolA)
                cand = sg*e;
                dup = false;
                for t = 1:size(seen,1)
                    if norm(seen(t,:) - cand) < 1e-9, dup = true; break, end
                end
                if dup, continue, end
                seen(end+1,:) = cand;   %#ok<AGROW>
                idx(end+1) = i;         %#ok<AGROW>
                sgn(end+1) = sg;        %#ok<AGROW>
            end
        end
    end
end

% edgeExtent is gone: its two jobs are now split by the rule this file runs on. The t-range
% itself is a VALUE and is computed exactly, inline, from the symbolic vertices; whether the edge
% runs off to infinity is a DECISION and stays the same numeric probe it always was.

function [V, keep] = uniqueRows(V)
    keep = true(size(V,1),1);
    for i = 2:size(V,1)
        for j = 1:i-1
            if keep(j) && norm(V(i,:) - V(j,:)) < 1e-9, keep(i) = false; break, end
        end
    end
    V = V(keep,:);
end
