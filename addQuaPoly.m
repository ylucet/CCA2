function h = addQuaPoly(f, g)
% addQuaPoly  Pointwise sum h = f+g of two QuaPoly functions (the geometry behind QuaPoly.add).
%
% objective: h(x) = f(x) + g(x) for all x. h's domain is the INTERSECTION of f's and g's domains
%   (finite exactly where both f and g are finite); outside that intersection h is implicitly
%   +infinity, the same convention QuaPoly itself already uses for a bounded domain (F(j,:) has a
%   0 entry on the side where the function is +infinity).
%
% [input]  f, g : QuaPoly, both operable (quadratic numerator, degree<=2 -- see assertOperable)
% [output] h    : QuaPoly with h.f(k,:) = the sum of whichever f-piece and g-piece overlap there
%
% METHOD: overlay f's and g's polyhedral subdivisions by pairwise convex-polygon clipping --
%   every face of f against every face of g (Sutherland-Hodgman half-plane clipping, adapted from
%   maxQuaPar.m's facePoly/clipByFace/clipPolyHalfPlane) -- then SUM the two quadratics on each
%   nonempty overlap cell. Unlike maxQuaPar's pointwise MAX, a sum is always exactly representable
%   as ONE quadratic per cell, so there is no decideWinner/splitCell case analysis at all: every
%   nonempty clip result is immediately a finished piece.
%
%   Reassembly (adapted from maxQuaPar.m's assemblePieces) matches each piece's edges into shared
%   half-edges to build the final V/E/F. The one substantive difference from maxQuaPar: its inputs
%   are always full-domain (finite everywhere), so an edge with no matching neighbour is an
%   internal-consistency ERROR there. Here f or g (or both) may have a genuinely bounded domain
%   (the ordinary QuaPoly case), so an edge with no neighbour is not an error -- it is exactly
%   where f's or g's own domain boundary passes through, and becomes a real boundary edge of h's
%   domain (F=0 on the open side), the same way any bounded QuaPoly already represents its domain.
%
% STATUS: implemented for two general QuaPoly objects (bounded or unbounded domain, any number of
%   faces each). Not yet extended to QuaPar (would need Ec/curved-edge clipping) or RatPol (would
%   need a common-denominator sum) -- see DESIGN.md's Implementation status section.

    % --- degenerate/trivial full-domain cases: no geometry needed at all ---
    if f.nv == 0 && g.nv == 0        % both finite everywhere: h is the single summed quadratic
        h = QuaPoly(f.f + g.f);
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

    % --- general case: pairwise face-vs-face clipping ---
    pieces = struct('V', {}, 'dirIn', {}, 'dirOut', {}, 'f', {});
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
        error('addQuaPoly:noOverlap', 'add: f''s and g''s domains do not overlap anywhere.');
    end
    h = assemblePiecesAdd(pieces);
end

% ============================================================================================
% ----- extracting a face's boundary as {V (CCW finite vertices), dirIn, dirOut} --------------
% (identical to maxQuaPar.m's facePoly, just applied to a QuaPoly instead of a QuaPar -- both
% classes share the same V/E/F/P field layout, so the same code works verbatim)
function poly = facePoly(obj, k)
    if obj.nv == 0
        poly.V = zeros(0,2); poly.dirIn = []; poly.dirOut = [];
        return
    end
    Pk = obj.P{k};
    n = numel(Pk);
    V = zeros(0,2); dirIn = []; dirOut = [];
    for t = 1:n
        j = abs(Pk(t)); s = sign(Pk(t));
        a = obj.E(j,1); b = obj.E(j,2);
        if obj.E(j,3) == 0   % ray: column 1 is always the finite apex
            apex = obj.V(a,:); d = obj.V(b,:) - apex; d = d/norm(d);
            if t == 1
                V(end+1,:) = apex; dirIn = d; %#ok<AGROW>
            else
                dirOut = d;
            end
        else
            if s > 0, V(end+1,:) = obj.V(b,:); else, V(end+1,:) = obj.V(a,:); end %#ok<AGROW>
        end
    end
    poly.V = flipud(V);
    poly.dirIn = dirOut;
    poly.dirOut = dirIn;
end

% ============================================================================================
% ----- clipping one convex poly by every boundary constraint of another (face intersection) --
% (identical to maxQuaPar.m's clipByFace/clipPolyHalfPlane/polyConstraints/
% insertPassthroughVertices/onOpenSegment/onOpenRay/sign2/crossingPoint/dedupConsecutive -- see
% that file's header HISTORY for the derivation and the bugs this exact code already survived)
function cell = clipByFace(polyK, polyL)
    cell = polyK;
    cons = polyConstraints(polyL);
    for i = 1:size(cons,1)
        cell = clipPolyHalfPlane(cell, cons(i,1:2), cons(i,3));
        if isempty(cell), return; end
    end
    cell = insertPassthroughVertices(cell, [polyK.V; polyL.V]);
    if size(cell.V,1) < 1 || (isempty(cell.dirIn) && size(cell.V,1) < 3)
        cell = []; return
    end
end

function cell = insertPassthroughVertices(cell, pts)
    if isempty(pts), return; end
    tol = 1e-7;
    for pi = 1:size(pts,1)
        p = pts(pi,:);
        again = true;
        while again
            again = false;
            nv = size(cell.V,1);
            if nv == 0 || any(all(abs(cell.V - p) < tol, 2)), break; end
            if isempty(cell.dirIn)
                for i = 1:nv
                    j = mod(i,nv)+1;
                    if onOpenSegment(cell.V(i,:), cell.V(j,:), p, tol)
                        cell.V = [cell.V(1:i,:); p; cell.V(i+1:end,:)];
                        again = true; break
                    end
                end
            else
                for i = 1:nv-1
                    if onOpenSegment(cell.V(i,:), cell.V(i+1,:), p, tol)
                        cell.V = [cell.V(1:i,:); p; cell.V(i+1:end,:)];
                        again = true; break
                    end
                end
                if ~again && onOpenRay(cell.V(1,:), cell.dirIn, p, tol)
                    cell.V = [p; cell.V]; again = true;
                elseif ~again && onOpenRay(cell.V(end,:), cell.dirOut, p, tol)
                    cell.V = [cell.V; p]; again = true;
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

function cons = polyConstraints(poly)
% BUGFIX (found while implementing addQuaPar.m): a BOUNDED poly (dirIn empty) is a closed
% cycle of nv edges (1,2),...,(nv-1,nv),(nv,1) -- the old "for i=1:nv-1" loop always dropped
% the closing edge (nv,1), so clipByFace never enforced that one constraint of polyL. Silent
% under-constraining, not caught by the existing test suite (see addQuaPolyTest.m's new
% missingClosingEdgeConstraintIsEnforced regression test). An UNBOUNDED poly's nv-1 real-vertex
% edges (1,2),...,(nv-1,nv) do NOT wrap (the two ends connect to rays, not to each other), so
% that case is unaffected: `last` below is nv-1 for it, exactly as before.
    nv = size(poly.V,1);
    cons = zeros(0,3);
    last = nv - 1;
    if isempty(poly.dirIn), last = nv; end   % bounded: also include the closing edge (nv,1)
    for i = 1:last
        jn = mod(i,nv) + 1;
        d = poly.V(jn,:) - poly.V(i,:);
        n = [d(2), -d(1)];
        cons(end+1,:) = [n, n*poly.V(i,:)']; %#ok<AGROW>
    end
    if ~isempty(poly.dirIn)
        n = [-poly.dirIn(2), poly.dirIn(1)];
        cons(end+1,:) = [n, n*poly.V(1,:)']; %#ok<AGROW>
        n = [poly.dirOut(2), -poly.dirOut(1)];
        cons(end+1,:) = [n, n*poly.V(end,:)']; %#ok<AGROW>
    end
end

function s = sign2(v, tol)
    s = zeros(size(v));
    s(v > tol) = 1; s(v < -tol) = -1;
end

function poly2 = clipPolyHalfPlane(poly, nrm, c)
    nv = size(poly.V,1);
    tol = 1e-9*(1+abs(c)+norm(nrm));
    val = poly.V*nrm' - c;
    unbounded = ~isempty(poly.dirIn);

    if unbounded
        inSign = sign2(poly.dirIn*nrm', tol);
        if inSign == 0, inSign = sign2(val(1), tol); end
        outSign = sign2(poly.dirOut*nrm', tol);
        if outSign == 0, outSign = sign2(val(end), tol); end
        st = [inSign; sign2(val,tol); outSign];
    else
        st = sign2(val,tol);
    end
    m = numel(st);

    if all(st <= 0)
        poly2 = poly; return
    end
    if all(st >= 0) && any(st > 0)
        poly2 = []; return
    end

    pairs = [(1:m-1)', (2:m)'];
    if ~unbounded, pairs = [pairs; m 1]; end
    cross = find((st(pairs(:,1)) > 0) ~= (st(pairs(:,2)) > 0));

    xpt = @(p) crossingPoint(poly, val, nrm, c, p, m, unbounded);

    if ~unbounded
        if numel(cross) ~= 2
            error('addQuaPoly:internal', ...
                'clipPolyHalfPlane: expected 2 crossings on a bounded poly, found %d.', numel(cross));
        end
        p1 = cross(1); p2 = cross(2); X1 = xpt(p1); X2 = xpt(p2);
        midIdx = mod((p1):(p2-1), nv) + 1;
        if isempty(midIdx) || all(st(midIdx) <= 0)
            Vnew = dedupConsecutive([X1; poly.V(midIdx,:); X2]);
        else
            keepIdx = mod((p2):(p1-1+nv), nv) + 1;
            Vnew = dedupConsecutive([X2; poly.V(keepIdx,:); X1]);
        end
        if size(Vnew,1) < 3, poly2 = []; return; end
        poly2.V = Vnew; poly2.dirIn = []; poly2.dirOut = [];
        return
    end

    if numel(cross) == 2
        p1 = cross(1); p2 = cross(2); X1 = xpt(p1); X2 = xpt(p2);
        vIdxMid = (p1):(p2-1);
        if st(1) <= 0
            keepBefore = poly.V(1:p1-1,:);
            keepAfter  = poly.V(p2:nv,:);
            poly2.V = dedupConsecutive([keepBefore; X1; X2; keepAfter]);
            poly2.dirIn = poly.dirIn; poly2.dirOut = poly.dirOut;
        else
            Vnew = dedupConsecutive([X1; poly.V(vIdxMid,:); X2]);
            if size(Vnew,1) < 3, poly2 = []; return; end
            poly2.V = Vnew; poly2.dirIn = []; poly2.dirOut = [];
        end
        return
    end
    if numel(cross) ~= 1
        error('addQuaPoly:internal', ...
            'clipPolyHalfPlane: expected 1 or 2 crossings on an unbounded poly, found %d.', numel(cross));
    end
    p = cross(1); X = xpt(p);
    if st(1) <= 0
        keepV = poly.V(1:min(p-1,nv),:);
        poly2.V = dedupConsecutive([keepV; X]);
        poly2.dirIn = poly.dirIn;
        poly2.dirOut = [-nrm(2), nrm(1)];
    else
        keepV = poly.V(max(p,1):nv,:);
        poly2.V = dedupConsecutive([X; keepV]);
        poly2.dirIn = [nrm(2), -nrm(1)];
        poly2.dirOut = poly.dirOut;
    end
    if isempty(poly2.V), poly2 = []; end
end

function pt = crossingPoint(poly, val, nrm, c, pairIdx, m, unbounded)
    if unbounded && pairIdx == 1
        v0 = poly.V(1,:)*nrm' - c; denom = poly.dirIn*nrm';
        pt = poly.V(1,:) - (v0/denom)*poly.dirIn;
        return
    end
    if unbounded && pairIdx == m-1
        vn = poly.V(end,:)*nrm' - c; denom = poly.dirOut*nrm';
        pt = poly.V(end,:) - (vn/denom)*poly.dirOut;
        return
    end
    if unbounded
        i = pairIdx - 1; j = pairIdx;
        vA = val(i); vB = val(j);
        pt = poly.V(i,:) + (vA/(vA-vB))*(poly.V(j,:)-poly.V(i,:));
    else
        i = pairIdx; j = mod(pairIdx, m) + 1;
        vA = val(i); vB = val(j);
        pt = poly.V(i,:) + (vA/(vA-vB))*(poly.V(j,:)-poly.V(i,:));
    end
end

function V = dedupConsecutive(V)
    tol = sqrt(eps);
    keep = true(size(V,1),1);
    for i = 2:size(V,1)
        if norm(V(i,:)-V(i-1,:)) < tol, keep(i) = false; end
    end
    if size(V,1) > 1 && norm(V(1,:)-V(end,:)) < tol, keep(end) = false; end
    V = V(keep,:);
end

% ============================================================================================
% ----- reassembling clipped+summed pieces into one QuaPoly -----------------------------------
function h = assemblePiecesAdd(pieces)
% Adapted from maxQuaPar.m's assemblePieces: same coordinate-dedup + half-edge pairing, but (a)
% no Ec/curved edges at all (QuaPoly is purely polyhedral) and (b) an edge with no matching
% neighbour is NOT an error -- f or g may have a bounded domain, so an unpaired edge is exactly
% where that domain boundary passes through, and becomes a genuine boundary edge of h (F=0 on the
% open side), using the same left/right convention as a paired edge (see maxQuaPar.m header
% HISTORY for the ray-orientation derivation this reuses).
    tol = 1e-6;
    n = numel(pieces);
    allV = cell(1,n); allE = cell(1,n);
    for p = 1:n
        piece = pieces(p);
        nv = size(piece.V,1);
        if isempty(piece.dirIn)
            Vp = piece.V;
            Ep = zeros(nv,3);
            for i = 1:nv
                Ep(i,:) = [i, mod(i,nv)+1, 1];
            end
        else
            Vp = [piece.V; piece.V(1,:)+piece.dirIn; piece.V(end,:)+piece.dirOut];
            ne = nv+1;
            Ep = zeros(ne,3);
            Ep(1,:) = [1, nv+1, 0];
            for i = 1:nv-1
                Ep(i+1,:) = [i, i+1, 1];
            end
            Ep(ne,:) = [nv, nv+2, 0];
        end
        allV{p} = Vp; allE{p} = Ep;
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

    HE = struct('piece', {}, 'a', {}, 'b', {}, 'isSeg', {}, 'rayOut', {});
    for p = 1:n
        Ep = allE{p}; gi = globalIdx{p};
        for e = 1:size(Ep,1)
            rayOut = ~Ep(e,3) && e == size(Ep,1);   % Ep's last row is always the OUTGOING ray
            HE(end+1) = struct('piece', p, 'a', gi(Ep(e,1)), 'b', gi(Ep(e,2)), ...
                'isSeg', Ep(e,3), 'rayOut', rayOut); %#ok<AGROW>
        end
    end
    haVec = [HE.a]; hbVec = [HE.b]; hsVec = [HE.isSeg];
    used = false(1,numel(HE));
    E = zeros(0,3); F = zeros(0,2);
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
            % passes through here). Same left/right rule as a paired edge below: the piece is on
            % the LEFT of stored (a,b) for a segment or an outgoing ray, on the RIGHT for an
            % incoming ray; the missing side is 0 (infinity), QuaPoly's own convention.
            if he.isSeg || he.rayOut
                F(end+1,:) = [he.piece, 0]; %#ok<AGROW>
            else
                F(end+1,:) = [0, he.piece]; %#ok<AGROW>
            end
            continue
        end
        used(opp) = true;
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
    h = QuaPoly(V, E, fMat, F);
end
