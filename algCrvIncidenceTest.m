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
