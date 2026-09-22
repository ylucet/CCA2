function [E, F, info] = algCrvIncidence(g, Vname)
% algCrvIncidence  The EDGE and FACE incidence of an AlgCrv: which vertices bound each edge, and
%   which faces each edge separates.
%
% objective: `algCrvVertices` names the vertices, which is what an object needs to have topology at
%   all; this is what a figure or an edge list needs on top of it. Without it the mesh is a set of
%   faces that happen to share curves, with nothing saying they are adjacent.
%
% [input]  g     : AlgCrv
%          Vname : as `algCrvVertices` returns it -- nv x 5 [edgeA tableA edgeB tableB rootIdx]
% [output] E     : ne x 5, each row [edgeIndex, table, vA, vB, flag]. vA and vB index ROWS OF Vname
%                  and are 0 where the edge has no such endpoint; flag is 1 for a segment, 0 for a
%                  ray or line (no endpoints found), 2 for a CLOSED curve -- AlgCon's own convention
%                  for the third column, and the reason that value exists (a whole ellipse has no
%                  endpoints to name).
%          F     : ne x 2, the faces on each side: [k1 k2], with k2 = 0 where only one face uses it.
%          info  : .edges, .segments, .rays, .oversubdivided, .shared
%
% ==================================================================================================
% WHAT IS EXACT HERE AND WHAT IS COMBINATORIAL
%
% Nothing is computed geometrically. Every edge is an entry already in `EcQ` or `EcA`, every vertex
% is already NAMED by the pair of edges that meet there, and a face already lists the edges it is
% bounded by. So the incidence is read off those three tables: an edge's endpoints are the vertices
% whose name mentions it, and the faces it separates are the faces whose row list references it.
%
% AN EDGE WITH MORE THAN TWO VERTICES is the one case the classic mesh cannot hold, because `E` has
% room for two endpoints. A conic can meet a face's other edges at three points or more, and then
% the single "edge" really is several arcs. `assembleQuaConCells` has exactly this limitation one
% degree down and resolves it the same way -- keep the first two and say so -- and the count is
% reported in `info.oversubdivided` rather than left for a reader to discover. Splitting such an
% edge into its arcs means ordering the vertices ALONG the curve, which is a comparison of algebraic
% numbers and the next piece of work rather than this one.

    if nargin ~= 2
        error('algCrvIncidence:nargin', 'expects (g, Vname).');
    end
    if ~isa(g, 'AlgCrv')
        error('algCrvIncidence:type', 'expects an AlgCrv; got a %s.', class(g));
    end
    if ~isempty(Vname) && size(Vname,2) ~= 5
        error('algCrvIncidence:shape', ...
              'Vname must be nv x 5 [edgeA tableA edgeB tableB rootIdx]; got %d columns.', ...
              size(Vname,2));
    end

    % ---- every edge either table actually uses, in a stable order ---------------------------------
    used = zeros(0,2);
    for k = 1:g.nf
        rows = g.FC{k};
        for t = 1:size(rows,1)
            tbl = 0;
            if size(rows,2) >= 4, tbl = rows(t,4); end
            key = [rows(t,1), tbl];
            if isempty(find(all(used == key, 2), 1)), used(end+1,:) = key; end %#ok<AGROW>
        end
    end

    ne = size(used,1);
    E = zeros(ne, 5);
    F = zeros(ne, 2);
    info = struct('edges', ne, 'segments', 0, 'rays', 0, 'oversubdivided', 0, 'shared', 0);

    for e = 1:ne
        idx = used(e,1);  tbl = used(e,2);

        % ---- the vertices whose NAME mentions this edge -------------------------------------------
        onIt = [];
        for v = 1:size(Vname,1)
            if (Vname(v,1) == idx && Vname(v,2) == tbl) || (Vname(v,3) == idx && Vname(v,4) == tbl)
                onIt(end+1) = v; %#ok<AGROW>
            end
        end

        flag = 0;
        if numel(onIt) >= 2
            flag = 1;
            info.segments = info.segments + 1;
            if numel(onIt) > 2, info.oversubdivided = info.oversubdivided + 1; end
        else
            info.rays = info.rays + 1;
        end
        vA = 0;  vB = 0;
        if numel(onIt) >= 1, vA = onIt(1); end
        if numel(onIt) >= 2, vB = onIt(2); end
        E(e,:) = [idx, tbl, vA, vB, flag];

        % ---- the faces that use it ---------------------------------------------------------------
        fs = [];
        for k = 1:g.nf
            rows = g.FC{k};
            for t = 1:size(rows,1)
                rt = 0;
                if size(rows,2) >= 4, rt = rows(t,4); end
                if rows(t,1) == idx && rt == tbl, fs(end+1) = k; break, end %#ok<AGROW>
            end
        end
        fs = unique(fs);
        if numel(fs) >= 2, info.shared = info.shared + 1; end
        F(e,1) = fs(1);
        if numel(fs) >= 2, F(e,2) = fs(2); end
    end
end
