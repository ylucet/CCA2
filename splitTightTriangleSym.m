function ts = splitTightTriangleSym(V)
% splitTightTriangleSym  Sub-triangles of V on each of which cPLQ's own closed form for that
%   sub-triangle's convex-edge count IS the convex envelope of x*y -- [COAP] Appendix A.4 (two
%   convex edges) and A.5 (three), applied recursively. Returns {V} when no split is needed,
%   which is every input cPLQ itself ever had.
%
% [input]  V  : 3x2, the triangle's vertices. SYMBOLIC and exact -- see WHY SYMBOLIC below.
%               The function is x*y; that is the only case this file is for (plq_1p.isCanonicalXY
%               guarantees it on the path that calls this).
% [output] ts : 1xk cell of 3x2 sym vertex matrices whose union is V.
%
% WHY THIS EXISTS. cPLQ's Step 1 applies its two-convex-edge closed form to the WHOLE triangle.
% That form touches x*y along both convex edges and is a valid convex MINORANT, but A.4 shows it
% is tight only over a sub-region -- and plq_1p.convexEnvelope1 never tests. Measured on
% conv{(2.5,1.5),(2,0),(0,0)}: the form returned reaches -0.2835 at (1,0), where both coordinates
% are >= 0, so x*y >= 0, the affine minorant 0 is admissible and the true envelope is >= 0. A
% too-small envelope gives a too-large conjugate. A three-convex-edge triangle has no such form at
% all, which is where that routine falls off the end and leaves an EMPTY envelope.
%
% WHY THE SPLIT IS A DOMAIN SPLIT -- and one claim that USED to be here is now false. It said an
% envelope-level split "cannot work: those faces are RATIONAL, and plq_1p.conjugateFunction reads
% its envelope with coeffs() and dispatches on the PIECE's nCE, so it raises
% symbolic:coeffs:NotAPolynomial". The obstacle was the DISPATCH, not the faces: conjugateFunction
% took nCE from the piece, so a rational A.3 face was read by the nCE == 2 branch. Since 2026-08-20
% it dispatches each face on that FACE's own geometry and the cross-face max can split a RATIONAL
% pair, so plq_1p.convexEnvelope1 calls THIS routine for a piece that did not come through
% triangulate and installs the sub-triangles' envelopes as several faces.
%
% The domain split remains the route for a piece that DOES come through triangulate, and it is the
% cheaper one: splitting the DOMAIN leaves every sub-piece on a path Step 2 already has -- A.4 yields one two-convex-edge
% sub-triangle, on which its own form IS tight, plus one one-convex-edge sub-triangle whose A.3
% rational envelope cPLQ's nCE == 1 branch derives analytically; A.5 yields two two-convex-edge
% ones, which recurse into A.4. And it is sound for Step 3 because a sup over a union is the max of
% the sups. See DECISIONS.md.
%
% WHY SYMBOLIC, and this is the whole reason this file exists beside convEnvCPLQ's own split.
% convEnvCPLQ computes the same geometry in DOUBLE precision. Taking the sub-triangles from it was
% tried and hangs: A.4's cevian has slope -sqrt(mh*mw), so its foot is IRRATIONAL, and sym of a
% double is EXACT -- a denominator near 2^53. Snapping the new vertices to the simplest rational
% within 1e-10 bounds the VERTEX denominators but not the downstream ones, because the conjugate is
% a rational function of those coordinates: a few squarings carry 1e5 to 1e25 and MuPAD's isAlways
% then decides nothing. The first conjugate ran 45+ minutes with no output. Carried symbolically the
% foot is an exact algebraic number, sqrt is something the symbolic layer keeps small, and the
% expressions stay the size the rest of the pipeline expects.
%
% REFUSALS ARE NOT ERRORS. Anything undecidable -- a comparison isAlways cannot settle, a cevian
% that lands outside its target edge -- returns the triangle UNSPLIT, which is exactly today's
% behaviour. So this can only ever subdivide a triangle whose envelope cPLQ was getting wrong.
    ts = {V};
    if size(V,1) ~= 3
        return
    end
    ts = splitRec(sym(V), 0);
end

% ================================================================================================
function ts = splitRec(V, depth)
% A.4 splits into a two-convex-edge half that keeps q1 unchanged (so its own curvature test comes
% out zero and it stops) plus a one-convex-edge half (which stops immediately), and A.5 into two
% two-convex-edge halves that each go through A.4 once. So the recursion is two deep by
% construction; the cap is a backstop against an undecidable comparison cycling, not a real limit.
    ts = {V};
    if depth > 4
        return
    end
    ce = classifyConvexEdgesSym(V);
    switch size(ce,1)
        case 2
            halves = splitTwoConvexSym(V, ce);
        case 3
            halves = splitThreeConvexSym(V);
        otherwise
            return                      % 0 or 1 convex edge: cPLQ's form is the envelope there
    end
    if isempty(halves)
        return                          % no split needed, or none this can certify
    end
    ts = {};
    for k = 1:numel(halves)
        ts = [ts, splitRec(halves{k}, depth+1)];   %#ok<AGROW>
    end
end

% ------------------------------------------------------------------------------------------------
function ce = classifyConvexEdgesSym(V)
% Convex edges of a triangle for the bilinear form x*y: an edge is convex iff its slope is > 0,
% since x*y restricted to y = m*x + q is m*x^2 + q*x. A VERTICAL edge is not convex -- it has no
% slope, and testing the geometry rather than the syntax is this file's own recurring lesson.
% Returns ce: (#convex) x 3 rows [slope, intercept, opposite-vertex-index].
    ed = [1 2; 2 3; 3 1];
    ce = sym(zeros(0,3));
    for t = 1:3
        i = ed(t,1); j = ed(t,2); opp = 6 - i - j;
        dx = V(j,1) - V(i,1);
        dy = V(j,2) - V(i,2);
        if provably(dx == 0)
            continue
        end
        m = dy/dx;
        if provably(m > 0)
            ce(end+1,:) = [m, V(i,2) - m*V(i,1), opp];   %#ok<AGROW>
        end
    end
end

% ------------------------------------------------------------------------------------------------
function halves = splitTwoConvexSym(V, ce)
% [COAP] Appendix A.4. The single quadratic q1, which touches x*y along BOTH convex edges, is a
% valid convex minorant but is tight only over a sub-region. Split by the cevian from one convex
% edge's far vertex into the OTHER convex edge, in the direction that makes q1 and the other half's
% A.3 formula agree in GRADIENT as well as value along the shared seam -- which is the line of
% slope -sqrt(mh*mw) through that far vertex, a direction depending only on mh and mw.
% convEnvCPLQ's splitTwoConvexEdges carries the full derivation; this is the same construction in
% exact arithmetic.
%
% THE "NO SPLIT NEEDED" TEST, IN CLOSED FORM. q1's quadratic part is
% [mh*mw, 2*sqrt(mh*mw), 1]/denom (see buildTwoEdgeSym), so its curvature along the weak edge AB is
%
%     ( sqrt(mh*mw)*dx + dy )^2 / ( sqrt(mh) + sqrt(mw) )^2 ,
%
% a PERFECT SQUARE -- which vanishes exactly when AB is PARALLEL to the cevian direction
% -sqrt(mh*mw). That is the measure-zero mirror-symmetric case where q1 already touches the chord
% along the whole weak edge (e.g. (0,0),(2,1),(1,2), where mh*mw = 1), and it is the same
% degeneracy convEnvCPLQ's header describes as "both candidate cevians degenerate". Testing it
% this way needs neither q1 nor the +/- branch selection.
%
% Empty means "no split": that degeneracy, or neither candidate cevian certifiably landing inside
% its target edge.
    halves = {};
    Bidx = double(ce(1,3)); Aidx = double(ce(2,3));
    Pidx = setdiff(1:3, [Bidx Aidx]);
    P = V(Pidx,:); A = V(Aidx,:); B = V(Bidx,:);
    mh = ce(1,1); qh = ce(1,2);      % edge P-A
    mw = ce(2,1); qw = ce(2,2);      % edge P-B

    g = sqrt(mh*mw);
    if provably(g*(B(1)-A(1)) + (B(2)-A(2)) == 0)
        return                       % q1 already touches the chord along the whole weak edge
    end

    slope = -g;                      % forced tangency direction; independent of the affine part

    [Ra, ta] = tangentCevianSym(A, slope, mw, qw, P, B);
    if provably(ta > 0) && provably(ta < 1)
        halves = {[P; A; Ra], [A; Ra; B]};        % {P,A,R} has 2 convex edges; {A,R,B} has 1
        return
    end
    [Rb, tb] = tangentCevianSym(B, slope, mh, qh, P, A);
    if provably(tb > 0) && provably(tb < 1)
        halves = {[P; Rb; B], [Rb; A; B]};
    end
end

% ------------------------------------------------------------------------------------------------
function halves = splitThreeConvexSym(V)
% [COAP] Appendix A.5. Split by the SMOOTH-FIT line through the middle (by second coordinate)
% vertex, hitting the opposite low-high edge at Pnew; each half then has exactly two convex edges.
%
% The line is NOT horizontal through the middle vertex in general. Both sub-envelopes touch x*y
% along the ENTIRE low-high edge, so their difference is a quadratic vanishing identically on that
% line and therefore factors as (low-high line) x (a second line); that second factor passes
% through the middle vertex, and is the unique direction making the two halves agree along the whole
% seam rather than only at its endpoints. convEnvCPLQ's splitThreeConvex records what the horizontal
% version got wrong and by how much.
    halves = {};
    ord = sortByYSym(V);
    if isempty(ord)
        return
    end
    vlow = V(ord(1),:); vmid = V(ord(2),:); vhigh = V(ord(3),:);

    if provably(vmid(1) - vlow(1) == 0) || provably(vhigh(1) - vlow(1) == 0) || ...
       provably(vhigh(1) - vmid(1) == 0)
        return
    end
    mh1 = (vmid(2)-vlow(2))/(vmid(1)-vlow(1));   qh1 = vlow(2) - mh1*vlow(1);
    mw  = (vhigh(2)-vlow(2))/(vhigh(1)-vlow(1)); qw  = vlow(2) - mw*vlow(1);
    mh2 = (vhigh(2)-vmid(2))/(vhigh(1)-vmid(1)); qh2 = vmid(2) - mh2*vmid(1);

    % Only three components of c1 - c2 are used, so write them out rather than build and simplify
    % two six-vectors. denom_i = mh_i + mw + 2*sqrt(mh_i*mw) = (sqrt(mh_i) + sqrt(mw))^2 > 0.
    den1 = mh1 + mw + 2*sqrt(mh1*mw);
    den2 = mh2 + mw + 2*sqrt(mh2*mw);
    if provably(den1 == 0) || provably(den2 == 0)
        return
    end
    d2 = 2*sqrt(mh1*mw)/den1 - 2*sqrt(mh2*mw)/den2;                  % x*y coefficient
    d3 = 1/den1 - 1/den2;                                            % y^2 coefficient
    d4 = (mh1*qw + mw*qh1)/den1 - (mh2*qw + mw*qh2)/den2;            % x coefficient
    d  = [sym(0), d2, d3, d4];                   % c1 - c2 vanishes on the low-high line
    q  = -d(3);
    pp = mw*q - d(2);
    if provably(pp + q*mw == 0)
        return
    end
    r  = (d(4) - qw*pp)/mw;
    xPnew = -(q*qw + r) / (pp + q*mw);           % p*x + q*y + r = 0 meets y = mw*x + qw
    Pnew  = simplify([xPnew, mw*xPnew + qw]);

    halves = {[vlow; vmid; Pnew], [vmid; vhigh; Pnew]};
end

% ------------------------------------------------------------------------------------------------
% THE +/- BRANCH OF [COAP] A.4's FORM IS VESTIGIAL, and this file therefore does not carry it.
% convEnvCPLQ's twoEdgeQuadPlain builds both signs and keeps the one that TOUCHES x*y along a convex
% edge. Substituting y = mh*x + qh into the s-form, with denom = mh + mw + s*2*sqrt(mh*mw):
%   x^2 coefficient : ( mh*mw + 2*s*sqrt(mh*mw)*mh + mh^2 ) / denom = mh*denom/denom = mh
%   x   coefficient : ( 2*s*sqrt(mh*mw)*qh + mh*qh + mw*qh ) / denom = qh*denom/denom = qh
%   constant        : ( qh^2 - (qh+qw)*qh + qh*qw ) / denom = 0
% -- for BOTH signs. So the form touches x*y along that whole edge either way, the touching test
% cannot separate them, and s = -1 is merely UNDEFINED when mh = mw. s = +1 always wins, which is
% what cPLQ hard-codes. Only the quantities the split actually needs are written out above.

% ------------------------------------------------------------------------------------------------
function [R, t] = tangentCevianSym(anchor, slope, mTarget, qTarget, Pp, otherEnd)
% The line through `anchor` with the given slope, met with the other convex edge
% y = mTarget*x + qTarget. t is R's fractional position from Pp to otherEnd along that edge, so
% 0 < t < 1 says R lies strictly inside it.
    R = sym([NaN NaN]); t = sym(NaN);
    den = mTarget - slope;
    if provably(den == 0)
        return
    end
    xR = (anchor(2) - slope*anchor(1) - qTarget) / den;
    yR = mTarget*xR + qTarget;
    R = simplify([xR, yR]);
    if provably(otherEnd(1) - Pp(1) ~= 0)
        t = simplify((xR - Pp(1)) / (otherEnd(1) - Pp(1)));
    elseif provably(otherEnd(2) - Pp(2) ~= 0)
        t = simplify((yR - Pp(2)) / (otherEnd(2) - Pp(2)));
    end
end

% ------------------------------------------------------------------------------------------------
function ord = sortByYSym(V)
% The vertex order by second coordinate, decided by PROVED comparisons only. Empty when two of them
% cannot be ordered -- in which case the caller leaves the triangle alone rather than guessing.
    ord = [];
    for i = 1:3
        for j = i+1:3
            if ~provably(V(i,2) < V(j,2)) && ~provably(V(j,2) < V(i,2)) && ...
               ~provably(V(i,2) == V(j,2))
                return
            end
        end
    end
    key = zeros(1,3);
    for i = 1:3
        for j = 1:3
            if i ~= j && provably(V(j,2) < V(i,2))
                key(i) = key(i) + 1;
            end
        end
    end
    if numel(unique(key)) ~= 3
        return                       % a tie in the second coordinate: no low/mid/high to speak of
    end
    [~, ord] = sort(key);
end

% ------------------------------------------------------------------------------------------------
function tf = provably(cond)
% isAlways with "unknown" folded into FALSE and no warning. Every caller is written so that a false
% answer means "leave it alone", so an undecidable comparison costs a split that was not made --
% never a wrong one.
    ws = warning('off', 'symbolic:sym:isAlways:TruthUnknown');
    try
        tf = logical(isAlways(cond, 'Unknown', 'false'));
    catch
        tf = false;
    end
    warning(ws);
end
