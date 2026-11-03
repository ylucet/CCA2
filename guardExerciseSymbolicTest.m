classdef guardExerciseSymbolicTest < matlab.unittest.TestCase
% The SYMBOLIC half of guardExerciseTest: the named refusals whose input is a symbolic `region` /
% `plq_1p`, each fired by one constructed input and required to report its own identifier
% (verifyError fails on any other error). guardExerciseTest's pin counts a guard asserted here as
% covered (it reads every other test file for the identifier).
%
% BUCKET: slow -- building a region is the symbolic pipeline (2026-11-02).
%
% Found on the way: conjConvexOverPiece:noCell never fired for an EMPTY piece -- region() returns a
% 0x0 array for an infeasible system and linearForm raised MATLAB:NonScalarInput first.

    methods (Test)
        function convEnvUnboundedNotTriangle(t)                              % reachable
            syms x y
            r = region([-x, -y, x-1, y-1], [x y]);                         % a bounded SQUARE
            t.verifyError(@() convEnvUnbounded(r, x*y, [x y]), 'convEnvUnbounded:notTriangle');
        end

        function convEnvUnboundedConvexAlongATriangleEdge(t)                 % reachable
            syms x y
            r = region([-x, -y, x+y-1], [x y]);                            % x^2 - y^2 bends along y = 0
            t.verifyError(@() convEnvUnbounded(r, x^2 - y^2, [x y]), 'convEnvUnbounded:convexAlongRay');
        end

        function convEnvUnboundedConvexAlongAHalfStripBase(t)                % reachable
            syms x y
            r = region([-x, x-1, -y], [x y]);                              % 0 <= x <= 1, y >= 0
            t.verifyError(@() convEnvUnbounded(r, x*y + x^2, [x y]), 'convEnvUnbounded:convexAlongRay');
        end

        function convEnvUnboundedNotPointed(t)                               % reachable
            % The strip 0 <= y <= 1: x^2 - y^2 is bounded below on it (its only recession
            % direction is x) and not convex, so it reaches the shape switch -- with no apex.
            syms x y
            r = region([-y, y-1], [x y]);
            t.verifyError(@() convEnvUnbounded(r, x^2 - y^2, [x y]), 'convEnvUnbounded:notPointed');
        end

        function conjConvexOverPieceNotConvex(t)                             % reachable
            syms x y s1 s2
            r = region([-x, -y, x+y-2], [x y]);
            t.verifyError(@() conjConvexOverPiece(r, [1 0; 0 -1], [0; 0], 0, [s1 s2]), ...
                          'conjConvexOverPiece:notConvex');
        end

        function conjConvexOverPieceNonAffineFacet(t)                        % reachable
            syms x y s1 s2
            r = region([-x, x^2 - y], [x y]);                              % above y = x^2, x >= 0
            t.verifyError(@() conjConvexOverPiece(r, eye(2), [0; 0], 0, [s1 s2]), ...
                          'conjConvexOverPiece:nonAffineFacet');
        end

        function conjConvexOverPieceEmptyPiece(t)                            % reachable
            syms x y s1 s2
            r = region([-x, -y, x+y+1], [x y]);                            % x, y >= 0, x + y <= -1
            t.verifyError(@() conjConvexOverPiece(r, eye(2), [0; 0], 0, [s1 s2]), ...
                          'conjConvexOverPiece:noCell');
        end

        function plqCheckRegionBuilders(t)                                   % reachable
            syms x y
            t.verifyError(@() plqCheck.triRegion(zeros(2,2), [x y]), 'plqCheck:triRegion');
            t.verifyError(@() plqCheck.wedgeRegion([0 0], [1 0], [0 -1], [x y]), 'plqCheck:wedgeRegion');
            t.verifyError(@() plqCheck.halfStripRegion([0 0], [1 0], [1 0], [x y]), ...
                          'plqCheck:halfStripRegion');
        end

        function plq1pPieceVarsWithNeither(t)                                % reachable
            % A constant piece with no domain: neither carries the two variables.
            t.verifyError(@() pieceVars(plq_1p([], symbolicFunction(sym(3)))), 'plq_1p:pieceVars');
        end
    end
end
