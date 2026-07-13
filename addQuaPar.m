function h = addQuaPar(f, g)
% addQuaPar  Pointwise sum h = f+g of two QuaPar functions (the geometry behind QuaPar.add).
%
% objective: h(x) = f(x) + g(x) for all x. h's domain is the INTERSECTION of f's and g's domains;
%   outside that intersection h is implicitly +infinity, same convention as QuaPoly/addQuaPoly.m.
%
% [input]  f, g : QuaPar, both operable (quadratic numerator, degree<=2 -- see assertOperable)
% [output] h    : QuaPar with h.f(k,:) = the sum of whichever f-piece and g-piece overlap there
%
% METHOD: generalizes addQuaPoly.m's pairwise convex-polygon clipping (Sutherland-Hodgman, every
%   face of f against every face of g) to allow ONE curved (parabolic-arc, "Ec") edge per face.
%   Every boundary edge -- straight or curved -- and every clipping constraint -- a face's linear
%   half-plane or its one conic (Ec) side -- is represented uniformly as a 6-coefficient row
%   [a b c d e f] for {a x^2+b xy+c y^2+d x+e y+f <= 0} (a=b=c=0 is the ordinary linear case, same
%   convention as QuaPar.edgeConics/evalConic). A straight edge is parametrized affinely in a
%   scalar t (segment: t in [0,1]; ray: t in [0,Inf)); the one curved edge is parametrized via an
%   axis-rotation of its parabola, x(u),y(u) each quadratic in a scalar u (see parabolaParam),
%   with u ranging directly between its two finite endpoints' own u-values. Either parametrization
%   substituted into any constraint row gives a plain polynomial in t or u (degree<=2 for a
%   straight edge, degree<=4 for the curved edge) whose real roots (via the ordinary quadratic
%   formula, or MATLAB's roots() for the quartic) are exactly where that edge crosses that
%   constraint -- found directly, not via an endpoint-sign-disagreement shortcut, since (unlike
%   the pure-linear case in addQuaPoly.m/maxQuaPar.m) a conic constraint can cross even a STRAIGHT
%   edge twice while agreeing in sign at both its endpoints (a dip-and-return), so both endpoints'
%   signs are never sufficient on their own.
%
% STATUS (incremental implementation, see DESIGN.md's Implementation status / II.4-II.6):
%   * IMPLEMENTED: f, g each have, per face, AT MOST ONE curved edge, and that curved edge is
%     always a bounded SEGMENT between two finite vertices (never a curved ray) -- the only shape
%     conjPieceCPLQ.m's conjBilinearXYoneCE actually produces today. A genuinely degenerate conic
%     (parabola collapsing to a double line) is rejected with a clear error, not silently wrong.
%   * NOT IMPLEMENTED (errors clearly rather than silently mishandling): a curved RAY edge; a
%     result cell that would need TWO curved edges (poly's own surviving arc plus a newly
%     introduced cut arc from the other operand's curved constraint) -- both raise
%     'addQuaPar:notImplemented' rather than dropping/misrepresenting geometry.

    if f.nv == 0 && g.nv == 0        % both finite everywhere: h is the single summed quadratic
        h = QuaPar(f.f + g.f);
        return
    end
    if f.nv == 0                     % f finite everywhere: add its row to every g face, keep g's
        h = g; h.f = g.f + f.f;      % own domain/mesh untouched
        return
    end
    if g.nv == 0
        h = f; h.f = f.f + g.f;
        return
    end

    pieces = struct('V', {}, 'dirIn', {}, 'dirOut', {}, 'curveAfter', {}, 'curveEc', {}, 'f', {});
    for k = 1:f.nf
        polyK = facePoly(f, k);
        for l = 1:g.nf
            polyL = facePoly(g, l);
            cell = clipByFace(polyK, polyL);
            if isempty(cell), continue; end
            cell.f = f.f(k,:) + g.f(l,:);
            pieces(end+1) = cell; %#ok<AGROW>
        end
    end
    if isempty(pieces)
        error('addQuaPar:noOverlap', 'add: f''s and g''s domains do not overlap anywhere.');
    end
    h = assemblePiecesAdd(pieces);
end

% ============================================================================================
% ----- extracting a face's boundary as {V, dirIn, dirOut, curveAfter, curveEc} ----------------
function poly = facePoly(obj, k)
% Boundary of face k, exactly like QuaPoly/maxQuaPar's facePoly (CCW finite vertices poly.V, plus
% ray directions dirIn/dirOut for an unbounded face), PLUS: poly.curveAfter (0, or the index i
% such that the edge from poly.V(i) to poly.V(mod(i,nv)+1) is the one curved arc) and poly.curveEc
% (that edge's conic, oriented so evalConic(curveEc,x)<=0 means "inside face k" -- the SAME
% "<=0 means inside" convention QuaPar.eval already uses via EC(...)*sign(P{i}(t)), so it composes
% directly with the linear half-plane constraints below).
    if obj.nv == 0
        poly.V = zeros(0,2); poly.dirIn = []; poly.dirOut = [];
        poly.curveAfter = 0; poly.curveEc = zeros(1,6);
        return
    end
    Pk = obj.P{k};
    n = numel(Pk);
    V = zeros(0,2); dirIn = []; dirOut = [];
    curveEndpoints = []; curveEcRaw = zeros(1,6);
    for t = 1:n
        j = abs(Pk(t)); s = sign(Pk(t));
        a = obj.E(j,1); b = obj.E(j,2);
        if obj.E(j,3) == 0   % ray: column 1 is always the finite apex
            if any(obj.Ec(j,:) ~= 0)
                error('addQuaPar:notImplemented', ...
                    'addQuaPar: curved (parabolic) ray edges are not implemented (edge %d).', j);
            end
            apex = obj.V(a,:); d = obj.V(b,:) - apex; d = d/norm(d);
            if t == 1
                V(end+1,:) = apex; dirIn = d; %#ok<AGROW>
            else
                dirOut = d;
            end
        else
            if s > 0, V(end+1,:) = obj.V(b,:); else, V(end+1,:) = obj.V(a,:); end %#ok<AGROW>
            if any(obj.Ec(j,:) ~= 0)
                if ~isempty(curveEndpoints)
                    error('addQuaPar:notImplemented', ...
                        'addQuaPar: face %d has more than one curved edge; not implemented.', k);
                end
                curveEndpoints = [obj.V(a,:); obj.V(b,:)]; %#ok<AGROW>
                curveEcRaw = obj.Ec(j,:) * s;
            end
        end
    end
    poly.V = flipud(V);
    poly.dirIn = dirOut;
    poly.dirOut = dirIn;
    poly.curveAfter = 0; poly.curveEc = zeros(1,6);
    if ~isempty(curveEndpoints)
        m = size(poly.V,1);
        tol = 1e-7;
        % An unbounded poly's real-vertex chain does NOT wrap (boundaryEdges never treats index m
        % as connecting back to 1 for it -- there is no edge there, just the two rays), so only
        % i=1..m-1 are valid positions there; a bounded poly's chain is a true cycle, i=1..m.
        lastPos = m - 1; if isempty(poly.dirIn), lastPos = m; end   % bounded (no rays): true cycle
        for i = 1:lastPos
            jn = mod(i,m) + 1;
            if (norm(poly.V(i,:)-curveEndpoints(1,:))<tol && norm(poly.V(jn,:)-curveEndpoints(2,:))<tol) || ...
               (norm(poly.V(i,:)-curveEndpoints(2,:))<tol && norm(poly.V(jn,:)-curveEndpoints(1,:))<tol)
                poly.curveAfter = i; poly.curveEc = curveEcRaw;
                break
            end
        end
        if poly.curveAfter == 0
            error('addQuaPar:internal', 'facePoly: could not locate curved edge in assembled boundary.');
        end
    end
end

% ============================================================================================
% ----- per-face boundary constraints, as uniform 6-coefficient rows --------------------------
function cons = polyConstraints(poly)
% One row per boundary edge of poly: [0 0 0 n1 n2 -c] for a straight edge with boundary value
% n.x<=c -- evalConic computes d*x+e*y+f as an ADDITIVE constant f, so the constant must be
% STORED AS -c to make evalConic<=0 equivalent to n.x<=c (the earlier, purely-linear
% clipPolyHalfPlane instead subtracted c explicitly when evaluating, so never needed this flip;
% routing everything through evalConic/evalConic-style substitution here does) -- or poly.curveEc
% itself for the one curved edge (already oriented "<=0 means inside poly", see facePoly).
    nv = size(poly.V,1);
    cons = zeros(0,6);
    last = nv - 1;
    if isempty(poly.dirIn), last = nv; end   % bounded: also include the closing edge (nv,1)
    for i = 1:last
        jn = mod(i,nv) + 1;
        if poly.curveAfter == i
            cons(end+1,:) = poly.curveEc; %#ok<AGROW>
        else
            d = poly.V(jn,:) - poly.V(i,:);
            n = [d(2), -d(1)];
            cons(end+1,:) = [0 0 0, n, -n*poly.V(i,:)']; %#ok<AGROW>
        end
    end
    if ~isempty(poly.dirIn)
        n = [-poly.dirIn(2), poly.dirIn(1)];
        cons(end+1,:) = [0 0 0, n, -n*poly.V(1,:)']; %#ok<AGROW>
        n = [poly.dirOut(2), -poly.dirOut(1)];
        cons(end+1,:) = [0 0 0, n, -n*poly.V(end,:)']; %#ok<AGROW>
    end
end

% ============================================================================================
% ----- parametrizing the one curved (parabolic) edge as x(u),y(u), u a scalar ----------------
function [Xc, Yc, p, q] = parabolaParam(Ec)
% For a parabola Ec=[a b c d e f] (b^2-4ac=0, not all zero): a x^2+bxy+cy^2 factors EXACTLY as
% s*(px+qy)^2 for s=+-1 (this is exactly what b^2=4ac buys you). Writing u=px+qy (the "squared"
% axis coordinate) and v as the complementary rotated coordinate, the conic becomes
% s*u^2 + A*u + B*v + C = 0, i.e. v = -(s*u^2+A*u+C)/B whenever B~=0 (the generic single-branch
% parabola case -- B==0 would mean the conic degenerates to a pair of parallel lines in the
% rotated frame, rejected below since QuaPar's own isParabola allows that degeneracy but this
% pipeline never needs it). Substituting v(u) back gives x(u),y(u) as quadratics in u.
% [output] Xc,Yc: 1x3 coefficient rows (highest power first) s.t. x(u)=Xc*[u^2;u;1], y(u)=Yc*[u^2;u;1]
%          p,q  : the axis-rotation coefficients (u = p*x+q*y), needed by the caller to place a
%                 known point's u-value.
    a=Ec(1); b=Ec(2); c=Ec(3); d=Ec(4); e=Ec(5); f=Ec(6);
    if abs(a) >= abs(c)
        s = sign(a); p = sqrt(s*a); q = b/(2*s*p);
    else
        s = sign(c); q = sqrt(s*c); p = b/(2*s*q);
    end
    nrm = p^2 + q^2;
    if nrm < 1e-12
        error('addQuaPar:degenerateConic', 'parabolaParam: Ec is degenerate (p=q=0).');
    end
    Alin = (d*p+e*q)/nrm; Blin = (-d*q+e*p)/nrm; Clin = f;
    if abs(Blin) < 1e-9*max(1,abs(Alin))
        error('addQuaPar:notImplemented', ...
            ['parabolaParam: conic degenerates to (at most) a pair of parallel lines in the ' ...
             'rotated frame, not a genuine single-branch parabola; not implemented.']);
    end
    v2 = -s/Blin; v1 = -Alin/Blin; v0 = -Clin/Blin;
    Xc = [-q*v2, p - q*v1, -q*v0] / nrm;
    Yc = [ p*v2, q + p*v1,  p*v0] / nrm;
end

function [x,y] = evalParam(Xc, Yc, u)
    x = Xc(1)*u.^2 + Xc(2)*u + Xc(3);
    y = Yc(1)*u.^2 + Yc(2)*u + Yc(3);
end

function coeffs = constraintAlongCurve(ecRow, Xc, Yc)
% Coefficients (highest power first, length 5) of u -> evalConic(ecRow, x(u), y(u)); degree<=4.
    a=ecRow(1); b=ecRow(2); c=ecRow(3); d=ecRow(4); e=ecRow(5); f=ecRow(6);
    X2 = conv(Xc,Xc); XY = conv(Xc,Yc); Y2 = conv(Yc,Yc);   % each length 5
    pad = @(Cc) [0 0 Cc];                                    % length 3 -> length 5
    coeffs = a*X2 + b*XY + c*Y2 + d*pad(Xc) + e*pad(Yc) + [0 0 0 0 f];
end

% ============================================================================================
% ----- restricting a constraint to a straight-edge/ray parametrization apex+t*dir ------------
function [A,B,C] = conicAlongRay(ecRow, apex, dir)
% Coefficients of t -> evalConic(ecRow, apex+t*dir) = A t^2 + B t + C.
    a=ecRow(1); b=ecRow(2); c=ecRow(3); d=ecRow(4); e=ecRow(5);
    A = a*dir(1)^2 + b*dir(1)*dir(2) + c*dir(2)^2;
    B = 2*a*apex(1)*dir(1) + b*(apex(1)*dir(2)+apex(2)*dir(1)) + 2*c*apex(2)*dir(2) + d*dir(1) + e*dir(2);
    C = QuaPar.evalConic(ecRow, apex);
end

function r = solveQuad(A,B,C)
    tol = 1e-9*(1+abs(A)+abs(B)+abs(C));
    if abs(A) <= tol
        if abs(B) <= tol, r = []; return; end
        r = -C/B; return
    end
    disc = B^2 - 4*A*C;
    if disc < -tol, r = []; return; end
    disc = max(disc,0);
    r = [(-B+sqrt(disc))/(2*A), (-B-sqrt(disc))/(2*A)];
end

function s = sign2scalar(v, tol)
    if v > tol, s = 1; elseif v < -tol, s = -1; else, s = 0; end
end

function v = asymptoticVal(A, B, C, tol)
% Sign-representative value of A t^2+B t+C as t->+inf: the leading nonzero coefficient.
    if abs(A) > tol, v = A; return; end
    if abs(B) > tol, v = B; return; end
    v = C;
end

% ============================================================================================
% ----- ordered boundary-edge descriptors for a poly (0 or 1 curved edges) --------------------
function segs = boundaryEdges(poly)
% One struct per boundary edge, in CCW WALK order: for an unbounded poly, the FIRST entry is
% always the dirIn ray and the LAST is always the dirOut ray (kind='ray', parametrized apex+t*dir,
% t in [0,Inf) -- note this is the ray's OWN natural direction, apex outward; dirIn's true WALK
% direction is the reverse of that, handled by the caller, see clipPolyByConstraint).
    nv = size(poly.V,1);
    segs = struct('kind',{}, 'apex',{}, 'dir',{}, 'tMax',{}, 'Xc',{}, 'Yc',{}, 'uA',{}, 'uB',{});
    unbounded = ~isempty(poly.dirIn);
    if unbounded
        segs(end+1) = struct('kind','ray','apex',poly.V(1,:),'dir',poly.dirIn,'tMax',Inf, ...
            'Xc',[],'Yc',[],'uA',[],'uB',[]);
        lastReal = nv - 1;
    else
        lastReal = nv;
    end
    for i = 1:lastReal
        j = mod(i,nv) + 1;
        if poly.curveAfter == i
            [Xc,Yc,p,q] = parabolaParam(poly.curveEc);
            uA = p*poly.V(i,1)+q*poly.V(i,2);
            uB = p*poly.V(j,1)+q*poly.V(j,2);
            segs(end+1) = struct('kind','curve','apex',[],'dir',[],'tMax',[], ...
                'Xc',Xc,'Yc',Yc,'uA',uA,'uB',uB); %#ok<AGROW>
        else
            segs(end+1) = struct('kind','seg','apex',poly.V(i,:),'dir',poly.V(j,:)-poly.V(i,:), ...
                'tMax',1,'Xc',[],'Yc',[],'uA',[],'uB',[]); %#ok<AGROW>
        end
    end
    if unbounded
        segs(end+1) = struct('kind','ray','apex',poly.V(end,:),'dir',poly.dirOut,'tMax',Inf, ...
            'Xc',[],'Yc',[],'uA',[],'uB',[]);
    end
end

function atoms = straightAtoms(seg, ecRow, tol)
% Sub-intervals of t in [0,1] (seg) or [0,Inf) (ray) with constant sign, split at every real root
% of evalConic(ecRow, apex+t*dir)=0. atoms(k) = {sign, ptStart, ptEnd} in NATURAL (t increasing)
% order; ptEnd=[] marks "extends to infinity" for the last atom of a ray.
    isRay = isinf(seg.tMax);
    [A,B,C] = conicAlongRay(ecRow, seg.apex, seg.dir);
    rr = solveQuad(A,B,C);
    if isRay
        rr = sort(rr(rr>tol));
        tsAll = [0, rr];
    else
        rr = sort(rr(rr>tol & rr<1-tol));
        tsAll = [0, rr, 1];
    end
    atoms = struct('sign',{},'ptStart',{},'ptEnd',{});
    n = numel(tsAll);
    last = n; if isRay, last = n; end
    for k = 1:last
        tS = tsAll(k);
        if k < n
            tE = tsAll(k+1);
            tm = (tS+tE)/2;
            val = A*tm^2+B*tm+C;
            ptE = seg.apex + tE*seg.dir;
        elseif isRay
            tE = Inf;
            val = asymptoticVal(A,B,C, tol);
            ptE = [];
        else
            break   % bounded: nothing after the last node
        end
        s = sign2scalar(val, tol);
        ptS = seg.apex + tS*seg.dir;
        atoms(end+1) = struct('sign',s,'ptStart',ptS,'ptEnd',ptE); %#ok<AGROW>
    end
end

function atoms = curveAtoms(seg, ecRow, tol)
% Sub-intervals of u between seg.uA and seg.uB (walk direction, possibly decreasing) with
% constant sign, split at every real root (within range) of evalConic(ecRow,x(u),y(u))=0.
    quart = constraintAlongCurve(ecRow, seg.Xc, seg.Yc);
    rr = roots(quart);
    rr = rr(abs(imag(rr)) < 1e-6*max(1,abs(rr)));
    rr = real(rr);
    uLo = min(seg.uA,seg.uB); uHi = max(seg.uA,seg.uB);
    rr = rr(rr>uLo+tol & rr<uHi-tol);
    if seg.uA <= seg.uB
        uWalk = [seg.uA, sort(rr(:)'), seg.uB];
    else
        uWalk = [seg.uA, sort(rr(:)','descend'), seg.uB];
    end
    atoms = struct('sign',{},'ptStart',{},'ptEnd',{});
    for k = 1:numel(uWalk)-1
        um = (uWalk(k)+uWalk(k+1))/2;
        [xm,ym] = evalParam(seg.Xc, seg.Yc, um);
        val = QuaPar.evalConic(ecRow, [xm ym]);
        s = sign2scalar(val, tol);
        [xS,yS] = evalParam(seg.Xc, seg.Yc, uWalk(k));
        [xE,yE] = evalParam(seg.Xc, seg.Yc, uWalk(k+1));
        atoms(end+1) = struct('sign',s,'ptStart',[xS yS],'ptEnd',[xE yE]); %#ok<AGROW>
    end
end

% ============================================================================================
% ----- locating the single surviving run of "kept" atoms -------------------------------------
% Both poly and (the region behind) ecRow are convex (poly is one face of a QuaPar; ecRow is one
% bounding half-plane/half-conic of another convex face), so their intersection is convex, hence
% connected: exactly one surviving run is expected. More than one signals either a genuine bug or
% an input violating that convexity assumption -- errors clearly either way, never silently picks
% one arbitrarily.
function runs = findRunsLinear(keep)
    runs = zeros(0,2);
    n = numel(keep); i = 1;
    while i <= n
        if keep(i)
            j = i;
            while j < n && keep(j+1), j = j+1; end
            runs(end+1,:) = [i j]; %#ok<AGROW>
            i = j+1;
        else
            i = i+1;
        end
    end
end

function [rStart, rEnd] = findSingleRunLinear(keep)
    runs = findRunsLinear(keep);
    if size(runs,1) ~= 1
        error('addQuaPar:notImplemented', ...
            ['clipPolyByConstraint: clipping produced %d disjoint surviving pieces; only a single ' ...
             'convex surviving piece is implemented.'], size(runs,1));
    end
    rStart = runs(1,1); rEnd = runs(1,2);
end

function [rStart, rEnd] = findSingleRunCyclic(keep)
    n = numel(keep);
    rStart = [];
    for a = 1:n
        prev = mod(a-2,n) + 1;
        if keep(a) && ~keep(prev), rStart = a; break; end
    end
    if isempty(rStart)
        error('addQuaPar:internal', 'findSingleRunCyclic: no discard->keep transition found.');
    end
    len = 0;
    for k = 0:n-1
        idx = mod(rStart-1+k, n) + 1;
        if ~keep(idx), break; end
        len = len + 1;
    end
    if len ~= nnz(keep)
        error('addQuaPar:notImplemented', ...
            ['clipPolyByConstraint: clipping produced multiple disjoint surviving pieces; only a ' ...
             'single convex surviving piece is implemented.']);
    end
    rEnd = mod(rStart-1+len-1, n) + 1;
end

% ============================================================================================
% ----- clipping poly by one constraint row [a b c d e f] (evalConic<=0 means inside) ---------
function poly2 = clipPolyByConstraint(poly, ecRow)
% Generalizes clipPolyHalfPlane (addQuaPoly.m/maxQuaPar.m) to a possibly-conic constraint against
% a poly that may itself have one curved edge. Every boundary edge is split into constant-sign
% "atoms" (straight edges via the ordinary quadratic formula on t; the curved edge via roots() on
% a quartic in u, see curveAtoms) and the single surviving run of kept atoms is rebuilt into the
% new poly. Introducing a NEW ray (case: an unbounded poly's entire original ray is discarded, and
% the surviving piece must extend to infinity along ecRow's own line instead) is only implemented
% for a LINEAR ecRow (reusing clipPolyHalfPlane's own formula) -- a conic ecRow needing an
% unbounded curved replacement ray is out of scope (curved edges are bounded segments only here)
% and errors clearly rather than being silently dropped.
    tol = 1e-9*(1+norm(ecRow));
    isConicConstraint = any(ecRow(1:3) ~= 0);
    unbounded = ~isempty(poly.dirIn);
    segs = boundaryEdges(poly);
    nSeg = numel(segs);

    allAtoms = struct('kind',{},'sign',{},'ptStart',{},'ptEnd',{});
    for si = 1:nSeg
        sg = segs(si);
        if strcmp(sg.kind,'curve')
            atoms = curveAtoms(sg, ecRow, tol);
        else
            atoms = straightAtoms(sg, ecRow, tol);
        end
        if unbounded && si == 1   % dirIn: true walk direction is the reverse of "apex outward"
            atoms = atoms(end:-1:1);
            for a = 1:numel(atoms)
                tmp = atoms(a).ptStart; atoms(a).ptStart = atoms(a).ptEnd; atoms(a).ptEnd = tmp;
            end
        end
        for a = 1:numel(atoms)
            allAtoms(end+1) = struct('kind',sg.kind,'sign',atoms(a).sign, ...
                'ptStart',atoms(a).ptStart,'ptEnd',atoms(a).ptEnd); %#ok<AGROW>
        end
    end

    nA = numel(allAtoms);
    keep = [allAtoms.sign] <= 0;
    if all(keep), poly2 = poly; return; end
    if ~any(keep), poly2 = []; return; end

    % "Asymptotically kept" = the ATOM THAT TOUCHES INFINITY (always allAtoms(1) for dirIn,
    % allAtoms(nA) for dirOut, by construction) is itself kept -- the authoritative signal for
    % "does the result still extend to infinity via this ray", exactly as clipPolyHalfPlane's own
    % st(1)/asymptotic-status branches decide it (NOT whether some root happens to fall on that
    % specific ray -- a ray can be kept near its apex and discarded toward infinity, which must
    % become a CLOSED end, not "ray survives, unchanged").
    startOpen = unbounded && keep(1);
    endOpen   = unbounded && keep(nA);

    if unbounded
        runs = findRunsLinear(keep);
        if size(runs,1) == 1
            idxList = runs(1,1):runs(1,2);
        elseif size(runs,1) == 2 && startOpen && endOpen && runs(1,1) == 1 && runs(2,2) == nA
            % Both original rays survive asymptotically, with a middle bulge discarded --
            % mirrors clipPolyHalfPlane's "both ray ends inside" 2-crossing case (a still-single,
            % still-connected boundary, bridged by one new chord/arc across the discarded gap).
            % Bridging while ALSO keeping both original rays is not implemented here.
            error('addQuaPar:notImplemented', ...
                ['clipPolyByConstraint: both ray ends survive with a discarded middle bulge; ' ...
                 'bridging that gap is not implemented.']);
        else
            error('addQuaPar:notImplemented', ...
                ['clipPolyByConstraint: clipping produced %d disjoint surviving pieces; only a ' ...
                 'single convex surviving piece is implemented.'], size(runs,1));
        end
    else
        [rStart, rEnd] = findSingleRunCyclic(keep);
        len = mod(rEnd-rStart, nA) + 1;
        idxList = mod((rStart:rStart+len-1)-1, nA) + 1;
    end

    V2 = zeros(0,2);
    curveAfter2 = 0; curveEc2 = zeros(1,6);
    if ~startOpen
        V2(end+1,:) = allAtoms(idxList(1)).ptStart; %#ok<AGROW>
    end
    for k = 1:numel(idxList)
        at = allAtoms(idxList(k));
        if k == numel(idxList) && endOpen, continue; end   % ptEnd is the open (infinite) end
        V2(end+1,:) = at.ptEnd; %#ok<AGROW>
        if strcmp(at.kind,'curve')
            if curveAfter2 ~= 0
                error('addQuaPar:notImplemented', ...
                    'clipPolyByConstraint: result would need more than one curved edge; not implemented.');
            end
            curveAfter2 = size(V2,1) - 1;
            curveEc2 = poly.curveEc;
        end
    end

    poly2.dirIn = []; poly2.dirOut = [];
    if startOpen
        % dirIn's own infinite end survives (possibly truncated closer to the apex, if there was
        % an interior root): direction UNCHANGED.
        poly2.dirIn = poly.dirIn;
    elseif unbounded && endOpen
        % This end is closed, but the OTHER end (dirOut) still reaches infinity, so the overall
        % result must remain unbounded: a brand-new ray extends from V2(1) along ecRow's own line.
        if isConicConstraint
            error('addQuaPar:notImplemented', ...
                ['clipPolyByConstraint: the result needs a new UNBOUNDED ray along a curved ' ...
                 'constraint; not implemented (curved edges are bounded segments only here).']);
        end
        nrm = ecRow(4:5);
        poly2.dirIn = [nrm(2), -nrm(1)];
    end
    if endOpen
        poly2.dirOut = poly.dirOut;
    elseif unbounded && startOpen
        if isConicConstraint
            error('addQuaPar:notImplemented', ...
                ['clipPolyByConstraint: the result needs a new UNBOUNDED ray along a curved ' ...
                 'constraint; not implemented (curved edges are bounded segments only here).']);
        end
        nrm = ecRow(4:5);
        poly2.dirOut = [-nrm(2), nrm(1)];
    end

    % Having reached this point, at least one atom was discarded (the "all kept" case already
    % returned early above), so a CLOSED (no-ray) result always has a genuine gap between V2(end)
    % and V2(1) needing its own new closing edge -- never just an accidental reuse of an existing one.
    isClosedShape = isempty(poly2.dirIn) && isempty(poly2.dirOut);
    if isClosedShape
        if isConicConstraint
            if curveAfter2 ~= 0
                error('addQuaPar:notImplemented', ...
                    'clipPolyByConstraint: result would need more than one curved edge; not implemented.');
            end
            curveAfter2 = size(V2,1);
            curveEc2 = ecRow;
        end
        % ecRow linear: the implicit straight closing edge (V2(end) back to V2(1)) is already
        % exactly correct -- no extra bookkeeping needed, same as the original clipPolyHalfPlane.
    end

    poly2.V = V2;
    poly2.curveAfter = curveAfter2;
    poly2.curveEc = curveEc2;
    if size(poly2.V,1) < 1 || (isempty(poly2.dirIn) && size(poly2.V,1) < 3)
        poly2 = []; return
    end
end

% ============================================================================================
% ----- clipping one face's boundary by every constraint of another ---------------------------
function cell = clipByFace(polyK, polyL)
    cell = polyK;
    cons = polyConstraints(polyL);
    for i = 1:size(cons,1)
        cell = clipPolyByConstraint(cell, cons(i,:));
        if isempty(cell), return; end
    end
    cell = insertPassthroughVertices(cell, [polyK.V; polyL.V]);
    if size(cell.V,1) < 1 || (isempty(cell.dirIn) && size(cell.V,1) < 3)
        cell = []; return
    end
end

function poly = insertPassthroughVertices(poly, pts)
% Same fix as addQuaPoly.m/maxQuaPar.m (re-insert a polyK/polyL vertex lying in the open interior
% of one of cell's own STRAIGHT edges, needed when two collinear edges make one half-plane clip a
% geometric no-op -- see maxQuaPar.m header HISTORY); the one curved edge (poly.curveAfter) is
% skipped by this straight-line test and left alone (its own two endpoints are always genuine
% clip-produced vertices already, never mid-arc pass-throughs in this pipeline).
    if isempty(pts), return; end
    tol = 1e-7;
    for pi = 1:size(pts,1)
        p = pts(pi,:);
        again = true;
        while again
            again = false;
            nv = size(poly.V,1);
            if nv == 0 || any(all(abs(poly.V - p) < tol, 2)), break; end
            if isempty(poly.dirIn)
                for i = 1:nv
                    if i == poly.curveAfter, continue; end
                    j = mod(i,nv) + 1;
                    if onOpenSegment(poly.V(i,:), poly.V(j,:), p, tol)
                        poly.V = [poly.V(1:i,:); p; poly.V(i+1:end,:)];
                        if poly.curveAfter > i, poly.curveAfter = poly.curveAfter + 1; end
                        again = true; break
                    end
                end
            else
                for i = 1:nv-1
                    if i == poly.curveAfter, continue; end
                    if onOpenSegment(poly.V(i,:), poly.V(i+1,:), p, tol)
                        poly.V = [poly.V(1:i,:); p; poly.V(i+1:end,:)];
                        if poly.curveAfter > i, poly.curveAfter = poly.curveAfter + 1; end
                        again = true; break
                    end
                end
                if ~again && onOpenRay(poly.V(1,:), poly.dirIn, p, tol)
                    poly.V = [p; poly.V];
                    if poly.curveAfter > 0, poly.curveAfter = poly.curveAfter + 1; end
                    again = true;
                elseif ~again && onOpenRay(poly.V(end,:), poly.dirOut, p, tol)
                    poly.V = [poly.V; p]; again = true;
                end
            end
        end
    end
end

function tf = onOpenSegment(a, b, p, tol)
    d = b - a; L = norm(d);
    if L < tol, tf = false; return; end
    t = dot(p-a, d) / L^2;
    if t <= tol/L || t >= 1 - tol/L, tf = false; return; end
    tf = norm(a + t*d - p) < tol;
end

function tf = onOpenRay(apex, dir, p, tol)
    dn = dir/norm(dir);
    t = dot(p-apex, dn);
    if t <= tol, tf = false; return; end
    tf = norm(apex + t*dn - p) < tol;
end

% ============================================================================================
% ----- reassembling clipped+summed pieces into one QuaPar -------------------------------------
function h = assemblePiecesAdd(pieces)
% Adapted from addQuaPoly.m's assemblePiecesAdd (unmatched edge -> genuine domain-boundary edge of
% h, not an error, since f/g may have a bounded domain -- see that file's header) generalized with
% maxQuaPar.m's Ec bookkeeping (piece.curveAfter/curveEc -> an Ec column in the final edge list) to
% build a QuaPar.
    tol = 1e-6;
    n = numel(pieces);
    allV = cell(1,n); allE = cell(1,n); allEc = cell(1,n);
    for p = 1:n
        piece = pieces(p);
        nv = size(piece.V,1);
        if isempty(piece.dirIn)
            Vp = piece.V;
            Ep = zeros(nv,3); Ecp = zeros(nv,6);
            for i = 1:nv
                Ep(i,:) = [i, mod(i,nv)+1, 1];
            end
        else
            Vp = [piece.V; piece.V(1,:)+piece.dirIn; piece.V(end,:)+piece.dirOut];
            ne = nv+1;
            Ep = zeros(ne,3); Ecp = zeros(ne,6);
            Ep(1,:) = [1, nv+1, 0];
            for i = 1:nv-1
                Ep(i+1,:) = [i, i+1, 1];
            end
            Ep(ne,:) = [nv, nv+2, 0];
        end
        if piece.curveAfter > 0
            % piece.curveEc uses this file's own "evalConic<=0 means inside" convention
            % throughout (needed so cell.curveEc stays valid across the MULTIPLE successive
            % clipPolyByConstraint calls within one clipByFace); the curved edge is always built
            % here as an isSeg=1 edge with `piece` on its LEFT (same Ep(i,:)=[i,i+1,1] convention
            % as every other edge), which per QuaPar's OWN storage convention requires
            % evalConic(Ec,x)>0 on the left -- i.e. >0 for THIS piece's interior -- the OPPOSITE
            % sign. Negate exactly once, here at final assembly, not earlier.
            Ecp(piece.curveAfter,:) = -piece.curveEc;
        end
        allV{p} = Vp; allE{p} = Ep; allEc{p} = Ecp;
    end

    V = zeros(0,2);
    globalIdx = cell(1,n);
    for p = 1:n
        Vp = allV{p}; gi = zeros(size(Vp,1),1);
        for r = 1:size(Vp,1)
            hit = 0;
            for q = 1:size(V,1)
                if norm(V(q,:)-Vp(r,:)) < tol, gi(r) = q; hit = 1; break; end
            end
            if ~hit, V(end+1,:) = Vp(r,:); gi(r) = size(V,1); end %#ok<AGROW>
        end
        globalIdx{p} = gi;
    end

    HE = struct('piece', {}, 'a', {}, 'b', {}, 'isSeg', {}, 'ec', {}, 'rayOut', {});
    for p = 1:n
        Ep = allE{p}; Ecp = allEc{p}; gi = globalIdx{p};
        for e = 1:size(Ep,1)
            rayOut = ~Ep(e,3) && e == size(Ep,1);
            HE(end+1) = struct('piece', p, 'a', gi(Ep(e,1)), 'b', gi(Ep(e,2)), ...
                'isSeg', Ep(e,3), 'ec', Ecp(e,:), 'rayOut', rayOut); %#ok<AGROW>
        end
    end
    haVec = [HE.a]; hbVec = [HE.b]; hsVec = [HE.isSeg];
    used = false(1,numel(HE));
    E = zeros(0,3); Ec = zeros(0,6); F = zeros(0,2);
    for hIdx = 1:numel(HE)
        if used(hIdx), continue; end
        used(hIdx) = true;
        he = HE(hIdx);
        if he.isSeg
            opp = find(~used & hbVec==he.a & haVec==he.b & hsVec, 1);
        else
            opp = find(~used & haVec==he.a & hbVec==he.b & ~hsVec, 1);
        end
        E(end+1,:) = [he.a, he.b, he.isSeg]; %#ok<AGROW>
        if isempty(opp)
            % No neighbour: a genuine domain-boundary edge of h (f's or g's own domain boundary
            % passes through here), same left/right convention as addQuaPoly.m's assemblePiecesAdd.
            Ec(end+1,:) = he.ec; %#ok<AGROW>
            if he.isSeg || he.rayOut
                F(end+1,:) = [he.piece, 0]; %#ok<AGROW>
            else
                F(end+1,:) = [0, he.piece]; %#ok<AGROW>
            end
            continue
        end
        used(opp) = true;
        ecRow = he.ec; if all(ecRow==0), ecRow = HE(opp).ec; end
        Ec(end+1,:) = ecRow; %#ok<AGROW>
        if he.isSeg
            F(end+1,:) = [he.piece, HE(opp).piece]; %#ok<AGROW>
        else
            if he.rayOut
                F(end+1,:) = [he.piece, HE(opp).piece]; %#ok<AGROW>
            else
                F(end+1,:) = [HE(opp).piece, he.piece]; %#ok<AGROW>
            end
        end
    end

    fMat = zeros(n,10);
    for p = 1:n, fMat(p,:) = pieces(p).f; end
    h = QuaPar(V, E, Ec, fMat, F);
end
