function [g, nMerged] = mergeSameQuadFaces(obj)
% mergeSameQuadFaces  Step 0 of `conj` and `biconj`: NORMALISE the subdivision.
%
% objective: delete every interior edge whose two sides carry the SAME quadratic AND whose two
%   faces have a convex union, fusing the faces it separated. The function is unchanged -- such
%   an edge carries no information -- but the MESH is, and the mesh is what every dispatch
%   downstream reads.
%
% [input]  obj      : any RatPar. Only a QuaPol (quadratic on a polyhedral subdivision) is
%                     normalised; anything else is returned untouched, since a curved edge or a
%                     rational face makes "the same quadratic on both sides" the wrong question.
% [output] g        : the normalised function, same class
%          nMerged  : how many faces disappeared (0 when nothing applied)
%
% WHY THIS IS STEP 0 AND NOT AN OPTIMISATION. `biconj`'s short-circuits are stated over SHAPES:
% McCormick needs one bilinear face on a four-vertex box, the diamond case needs one face on a
% rotated box, separability needs one face on a product domain (biconjCPLQ.m, and ALGORITHM.md's
% ordering). The unit square handed in as two triangles sharing its diagonal is the SAME function
% as the unit square handed in as one face, and the second answers in 0 s while the first used to
% go the whole way round -- triangulate, per-piece conjugate, Step 3, second conjugation -- and
% come out wrong (biconjugateTest.biconjugateOverATwoFaceSubdivisionIsTheEnvelope). The answer
% must not depend on how the caller drew the mesh, so the mesh is normalised before the dispatch
% reads it.
%
% THE CONVEXITY GUARD IS NOT OPTIONAL, and it was measured rather than assumed. Two faces sharing
% an edge can have a NON-CONVEX union -- three wedges round the origin, the outer two carrying the
% same affine piece, fuse into a 225-degree reflex wedge -- and this representation does not
% support one: the merged mesh builds without complaint and then `eval` returns +inf at a point
% inside the reflex face, because a face is read as the intersection of the half-planes its edges
% define. So a candidate merge is BUILT and its merged face checked with `orderEdges`, whose
% isConvex flag is exactly this question; a merge that fails is skipped and its edge kept.
% Collinear boundary vertices pass that test, which is what lets strips merge.
%
% WHAT "THE SAME QUADRATIC" MEANS HERE: the two faces' coefficient rows are EXACTLY equal. No
% tolerance, deliberately. A tolerance would make the normalisation depend on a threshold nobody
% downstream knows about, and two faces that agree to within 1e-12 are two faces -- merging them
% changes the function, however slightly, which is not what a normaliser may do. Callers that
% want fuzzy merging should round their own coefficients first.
%
% WHAT IT DOES NOT DO: it leaves the vertices of the fused faces alone. Deleting an interior edge
% can leave a vertex of degree 2 whose two remaining edges are COLLINEAR -- a boundary point that
% is no longer a corner -- and the shape tests above count vertices (`obj.nv ~= 4`). Two triangles
% sharing a diagonal do not produce one (their four corners all stay corners), which is the case
% this exists for; a fan of strips would, and would then merge without becoming a box. See
% ALGORITHM.md.
%
% An edge with the SAME face on both sides is left alone as well: it is a slit, not a shared
% boundary, and removing it would change the domain rather than the mesh.

    g = obj;
    nMerged = 0;

    % OPT-OUT, and it exists for one regression test rather than for tuning.
    % biconjugateTest.twoFaceConjugateIsExactIncludingTheLensCell pins a silent wrong answer that
    % only appears on the two-face route (the lens where the two triangles' conjugates overlap and
    % s1, s2 tie at every vertex). Normalising that input away would leave the fix unpinned, so the
    % test sets CCA2_NO_STEP0 and keeps exercising the multi-face path deliberately.
    if ~isempty(getenv('CCA2_NO_STEP0')), return, end

    % QuaPol only: needs a mesh, no conics (Pol), one polynomial per face (Qua).
    if ~(isa(obj, 'Qua') && isa(obj, 'Pol')), return, end
    if ~ismethod(obj, 'isDomBounded'), return, end          % QuaParCPLQ carries no mesh
    if isempty(obj.nf) || obj.nf < 2 || obj.ne < 1, return, end

    % ONE MERGE AT A TIME, restarting after each: dropping an edge renumbers the faces, the edges
    % and the vertices, so a batch pass computed against the original numbering would be applying
    % a stale plan. The meshes here are small (a handful of faces) and each step is one
    % construction, so the simple loop is the right trade.
    changed = true;
    while changed
        changed = false;
        for j = 1:g.ne
            a = g.F(j,1); b = g.F(j,2);
            if a > 0 && b > 0 && a ~= b && isequal(g.f(a,:), g.f(b,:))
                trial = dropInteriorEdge(g, j, a, b);
                if ~isempty(trial)
                    g = trial;
                    nMerged = nMerged + 1;
                    changed = true;
                    break
                end
            end
        end
    end
end

% ------------------------------------------------------------------------------------------------
function h = dropInteriorEdge(g, j, a, b)
% Delete edge j and fuse faces a and b (b's rows go away, a keeps the quadratic). Empty when the
% result is not a well-formed mesh with a CONVEX merged face -- the caller then keeps the edge.
    h = [];

    keepF = true(g.nf, 1); keepF(b) = false;
    fNew = g.f(keepF, :);
    fMap = zeros(g.nf, 1);
    fMap(keepF) = 1:nnz(keepF);
    fMap(b) = fMap(a);                       % everything that was face b is now face a

    keepE = true(g.ne, 1); keepE(j) = false;
    E = g.E(keepE, :);
    F = g.F(keepE, :);
    Fnew = zeros(size(F));
    nz = F > 0;
    Fnew(nz) = fMap(F(nz));

    % AN EDGE THAT NOW SEPARATES A FACE FROM ITSELF goes with it. Two half-planes carrying the
    % same quadratic are drawn as two opposite RAYS: dropping one of them leaves the other with
    % the merged face on both sides, which is not a boundary of anything and which orderEdges
    % refuses (a face must carry zero or two rays). Only edges that separated two DIFFERENT faces
    % before this merge are removed this way -- a slit the caller drew stays a slit.
    fused = (F(:,1) > 0) & (F(:,2) > 0) & (F(:,1) ~= F(:,2)) & (Fnew(:,1) == Fnew(:,2));
    E = E(~fused, :); Fnew = Fnew(~fused, :);

    if isempty(E)
        % Every edge was interior: the domain is the whole plane carrying one quadratic, which is
        % the one-argument QuaPol. Convexity is not in question there.
        h = QuaPol(fNew(1,:));
        return
    end

    used = false(g.nv, 1);
    used(E(:,1:2)) = true;
    vMap = zeros(g.nv, 1);
    vMap(used) = 1:nnz(used);
    E(:,1:2) = vMap(E(:,1:2));

    try
        cand = QuaPol(g.V(used,:), E, fNew, Fnew);
        [~, isConv] = cand.orderEdges(fMap(a));
    catch
        return                                % not a mesh this representation can express
    end
    if ~isConv, return, end                   % a reflex face: eval cannot read it (see header)
    h = cand;
end
