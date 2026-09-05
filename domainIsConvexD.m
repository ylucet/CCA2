function tf = domainIsConvexD(V, E, F, tol)
% domainIsConvexD  Is the domain of a mesh CONVEX. Doubles, for createDom.
%
% objective: decide convexity of `dom f` from the mesh alone, without walking its boundary.
%
% [input]  V : nv x 2 vertices; E : ne x 3 [v1 v2 isSegment]; F : ne x 2 face indices, 0 = +infinity
%          tol : (optional) relative tolerance, default sqrt(eps)
% [output] tf : logical
%
% ------------------------------------------------------------------------------------------------
% WHY THIS EXISTS. `createDom` used to take its `isConvex` from `orderEdges(0)`, which walks the
% boundary, and that is WRONG for an unbounded domain: QuaPol.examples{16} -- the half-star, five
% rays from the origin -- was reported convex although (-1,-1) and (1,0) are in it and their
% midpoint (0,-0.5) is not. Two boundary rays do not say which of the two cones they bound is the
% domain, and the walk assumed the acute one. DECISIONS.md 2026-09-05.
%
% THE ALGORITHM, which is an H-representation containment test rather than a walk. A polyhedral set
% is convex iff it is the intersection of the closed half-planes of its own boundary edges, so:
%
%   1. collect every recession direction of the mesh -- one per RAY edge;
%   2. for each BOUNDARY edge (exactly one finite side), take the line through its endpoints,
%      oriented so the finite side is `l >= 0`;
%   3. require every vertex to satisfy `l >= 0` and every recession direction to satisfy the
%      homogeneous part `<a,u> >= 0`.
%
% Necessity is immediate. Sufficiency: let Omega be that intersection. Each face is
% conv(its vertices) + cone(its ray directions), so step 3 gives dom inside Omega; and for x in
% Omega, the segment from an interior point of dom to x cannot leave dom, because leaving means
% crossing some boundary edge's line, and a straight segment that crosses a line stays on the far
% side -- which would put l_j(x) < 0. So dom = Omega, hence convex.
%
% O(ne * nv), not the O(ne) a boundary walk would be. Deliberate: nv is small here and each
% predicate is a 2x2 determinant. `isConvexQ` runs the same test over the INTEGERS, with no
% tolerance at all, and is what the exact pipeline uses; this is the doubles twin for the legacy
% mesh property.
%
% TOLERANCE. Each test is a signed DISTANCE -- the line's normal is divided out -- so `tol` is
% comparable across edges of different lengths, and it is scaled by the magnitude of the point so
% that a vertex far from the origin is not judged more harshly than a near one.

    if nargin < 4, tol = sqrt(eps); end
    tf = true;
    if isempty(F) || isempty(E) || isempty(V), return, end
    % A THIN mesh -- every edge with +infinity on both sides -- has no boundary edge, so the loop
    % below would find nothing to violate and return true vacuously. That is not an answer: a
    % non-collinear chain of edges has a 1-dimensional NON-convex domain. Refuse it by name;
    % createDom decides those before it gets here.
    if all(F(:) == 0)
        error('PLQ:domainIsConvexD:thinMesh', ...
            ['every edge has +infinity on both sides, so the domain has dimension < 2 and this ' ...
             'test says nothing about it.']);
    end

    % Every recession direction of the mesh, normalised.
    dirs = zeros(0,2);
    for k = 1:size(E,1)
        if E(k,3) == 0
            u = V(E(k,2),:) - V(E(k,1),:);
            n = norm(u);
            if n > 0, dirs(end+1,:) = u / n; end %#ok<AGROW>
        end
    end

    for j = 1:size(F,1)
        iL = F(j,1);  iR = F(j,2);
        if (iL ~= 0 && iR ~= 0) || (iL == 0 && iR == 0)
            continue     % an interior edge is not a domain constraint; a dangling one is not a mesh
        end
        v1 = V(E(j,1),:);  v2 = V(E(j,2),:);
        a = [-(v2(2)-v1(2)), v2(1)-v1(1)];      % the LEFT normal of the edge direction
        na = norm(a);
        if na == 0, continue, end
        c = v1(1)*v2(2) - v2(1)*v1(2);
        s = 1;  if iL == 0, s = -1; end          % the finite side is the RIGHT one
        a = s*a/na;  c = s*c/na;                 % now l is a signed distance

        for i = 1:size(V,1)
            if a*V(i,:)' + c < -tol*max(1, norm(V(i,:)))
                tf = false;  return
            end
        end
        for r = 1:size(dirs,1)
            if a*dirs(r,:)' < -tol
                tf = false;  return
            end
        end
    end
end
