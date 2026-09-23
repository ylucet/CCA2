classdef algCrvIncidenceTest < matlab.unittest.TestCase
% [fast] The EDGE and FACE incidence of an AlgCrv -- which vertices bound each edge and which faces
% each edge separates. Everything here is read off the tables, so every assertion is a consistency
% statement about them rather than a stored number.
%
% The fixture is the one AlgCrvTest uses: the cone's upper branch (|x|) against the constant 1, whose
% algebraic edge is the unit circle, cut by the line y = 0 at (-1,0) and (1,0).

    methods (Test)
        function anEDGECarriesTheVerticesThatNAMEIt(t)
            g = fixture();
            V = algCrvVertices(g);
            [E, F, info] = algCrvIncidence(g, V);
            t.verifyEqual(info.edges, 2, 'the circle and the line');
            t.verifyEqual(size(E,1), 2);
            t.verifyEqual(size(F,1), 2);
            t.verifyEqual(info.segments, 2, 'each edge is cut twice, so each is a segment');
            t.verifyEqual(info.rays, 0);
            t.verifyEqual(info.oversubdivided, 0, 'and neither carries a third vertex');

            % every endpoint must be a vertex whose own NAME mentions that edge
            for e = 1:size(E,1)
                for c = [3 4]
                    v = E(e,c);
                    t.verifyGreaterThan(v, 0, 'a segment has both endpoints');
                    named = (V(v,1) == E(e,1) && V(v,2) == E(e,2)) || ...
                            (V(v,3) == E(e,1) && V(v,4) == E(e,2));
                    t.verifyTrue(named, 'an endpoint must name the edge it bounds');
                end
            end
        end

        function aFACEListedForAnEdgeMustUSEThatEdge(t)
            g = fixture();
            V = algCrvVertices(g);
            [E, F] = algCrvIncidence(g, V);
            for e = 1:size(E,1)
                t.verifyEqual(F(e,1), 1, 'one face here, so it owns every edge');
                t.verifyEqual(F(e,2), 0, 'and nothing is shared');
                rows = g.FC{F(e,1)};
                hit = false;
                for r = 1:size(rows,1)
                    tbl = 0;
                    if size(rows,2) >= 4, tbl = rows(r,4); end
                    if rows(r,1) == E(e,1) && tbl == E(e,2), hit = true; break, end
                end
                t.verifyTrue(hit, 'the face must reference the edge it is listed for');
            end
        end

        function anEdgeWithNOVertexIsARAYRatherThanASegment(t)
        % The other branch of the flag, and the reason AlgCon's third column has more than one
        % value: an edge nothing cuts has no endpoints to name.
            g = fixture([0 0 0 0 1 -3]);            % the line y = 3 misses the unit circle
            V = algCrvVertices(g);
            [E, ~, info] = algCrvIncidence(g, V);
            t.verifyEmpty(V, 'nothing meets');
            t.verifyEqual(info.rays, 2, 'so both edges are unbounded as far as this mesh knows');
            t.verifyEqual(unique(E(:,5)).', 0, 'and the flag says so');
        end

        function anEdgeWithTHREEVerticesIsSplitIntoItsARCS(t)
        % The circle bounds TWO quarter-disks and carries three vertices, (1,0), (0,1), (-1,0): one
        % edge, two arcs that bound a face each, and a lower half that bounds neither. Before the
        % split it was one row holding two of the three and `oversubdivided` said so.
            g = quarters();
            V = algCrvVertices(g);
            [E, F, info] = algCrvIncidence(g, V);
            t.verifyEqual(info.oversubdivided, 2, 'the circle and the line y = 0 carry three each');
            t.verifyEqual(info.split, 2, 'and both are split');
            t.verifyEqual(info.untraced, 0);
            t.verifyGreaterThanOrEqual(info.unowned, 1, 'the lower half bounds no face');

            circ = find(E(:,1) == 1 & E(:,2) == 0);
            t.verifyEqual(numel(circ), 2, 'the circle is two arcs of the mesh');
            [~, iv] = algCrvVertices(g);  X = iv.X;
            for r = circ.'
                ends = sortrows(round(X(E(r,3:4),:), 9));
                k = F(r,1);
                t.verifyEqual(F(r,2), 0, 'each quarter arc bounds one face');
                if k == 1, want = [0 1; 1 0]; else, want = [-1 0; 0 1]; end
                t.verifyEqual(ends, want, 'AbsTol', 1e-9, ...
                              'an arc runs between the two vertices ITS face has on the circle');
            end
        end

        function everyVERTEXOfAFaceEndsExactlyTWOOfItsArcs(t)
        % The definition of a face's boundary as a closed chain: at each of its vertex POINTS
        % exactly two of the face's arcs end. One row per edge breaks this on the circle, which
        % would give face 1 the whole circle's row with a third vertex hanging off it.
            g = quarters();
            [V, iv] = algCrvVertices(g);
            [E, F] = algCrvIncidence(g, V, iv.X);
            for k = 1:g.nf
                mine = find(F(:,1) == k | F(:,2) == k);
                P = zeros(0,2);
                for r = mine.'
                    for c = [3 4]
                        if E(r,c) > 0, P(end+1,:) = round(iv.X(E(r,c),:), 9); end %#ok<AGROW>
                    end
                end
                [~, ~, j] = unique(P, 'rows');
                t.verifyEqual(unique(accumarray(j, 1)).', 2, ...
                              sprintf('face %d: every vertex ends exactly two arcs', k));
            end
        end

        function theARGUMENTSAreChecked(t)
            g = fixture();
            t.verifyError(@() algCrvIncidence(g), 'algCrvIncidence:nargin');
            t.verifyError(@() algCrvIncidence(5, zeros(0,5)), 'algCrvIncidence:type');
            t.verifyError(@() algCrvIncidence(g, zeros(1,3)), 'algCrvIncidence:shape');
        end
    end
end

% ==================================================================================================
function g = fixture(line)
    if nargin < 1, line = [0 0 0 0 1 0]; end        % y = 0
    FP = zeros(3,3,3);
    FP(1,1,3) =  1;  FP(3,1,1) = -1;  FP(1,3,1) = -1;
    num = [0 0 0 0, 0, 0, 0, 0, 0, 1];
    EcA = struct('FP', FP, 'FB', 2, 'num', num, 'den', 1);
    FC  = {[1 1 0 1; 1 1 0 0]};
    g = AlgCrv(zeros(0,3), zeros(0,3), ratQ.conic(line), num, 1, zeros(0,2), FC, {[]}, 0, EcA);
end

function g = quarters()
% The unit disk's upper half, split by x = 0: face 1 is the quarter with x >= 0, face 2 the other.
% Edges: 1 the circle (inside is <= 0), 2 the line y = 0, 3 the line x = 0. Constant 1 on both.
    EcQ = [ratQ.conic([1 0 1 0 0 -1]); ratQ.conic([0 0 0 0 1 0]); ratQ.conic([0 0 0 1 0 0])];
    num = [0 0 0 0, 0, 0, 0, 0, 0, 1];
    FC  = {[1 -1 0 0; 2 1 0 0; 3  1 0 0]; ...
           [1 -1 0 0; 2 1 0 0; 3 -1 0 0]};
    g = AlgCrv(zeros(0,3), zeros(0,3), EcQ, [num; num], [1; 1], zeros(0,2), FC, {[], []}, [0; 0]);
end
