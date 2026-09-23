classdef algCrvVerticesTest < matlab.unittest.TestCase
% [fast] Naming the VERTICES of an AlgCrv, where an edge may be a conic or a curve of higher degree.
% The assertions are geometric: a named vertex must be a point where BOTH its edges vanish and every
% row of its face holds, and the count must be the one the configuration actually has.
%
% The fixture is the cone z^2 = x^2 + y^2 against the constant face q = 1, whose algebraic edge is
% the unit circle, cut by coordinate lines whose vertices can be written down.

    methods (Test)
        function aLINEMeetsTheCURVEWhereItShould(t)
        % y = 0 crosses the unit circle at (-1,0) and (1,0), and nowhere else.
            g = fixture([0 0 0 0 1 0]);            % the line y = 0 as a conic row
            [V, info] = algCrvVertices(g);
            t.verifyEqual(info.refused, 0, 'a line and a circle need no refusal');
            t.verifyEqual(size(V,1), 2, 'exactly two crossings');
            % One edge from each TABLE, whichever order the face happens to list them in -- the
            % name records the tables, not a convention about which comes first.
            for i = 1:size(V,1)
                t.verifyEqual(sort([V(i,2) V(i,4)]), [0 1], ...
                              'each vertex joins a CONIC edge to an ALGEBRAIC one');
            end
            t.verifyEqual(sort(V(:,5)).', [1 2], 'and they are the two roots, in order');
        end

        function aVERTICALLineIsNamedToo(t)
        % x = 0 has no y at all, so eliminating y refuses it -- the resultant needs both arguments
        % to have degree in the variable it removes. Measured before this was handled: 66 of 176
        % conic-by-curve pairs of the real fixture refused for exactly this reason. Eliminating x
        % instead names the same two vertices.
            g = fixture([0 0 0 1 0 0]);            % the line x = 0
            [V, info] = algCrvVertices(g);
            t.verifyEqual(info.refused, 0, 'a vertical line must not refuse');
            t.verifyEqual(size(V,1), 2, '(0,-1) and (0,1)');
        end

        function aLINEThatMISSESTheCurveNamesNothing(t)
        % The other direction, without which the tests above would pass on code that named points
        % it had not checked: y = 3 does not meet the unit circle.
            g = fixture([0 0 0 0 1 -3]);
            [V, ~] = algCrvVertices(g);
            t.verifyEmpty(V, 'no real crossing, no vertex');
        end

        function aCORNEROnlyTheSECONDFaceHasIsNamed(t)
        % Two quarter-disks share the circle and the line y = 0, and (-1,0) is a corner of the
        % SECOND one only. Deduplicating by edge pair skipped the pair on its second face, so the
        % corner went unnamed and the circle looked like it carried two vertices instead of three.
            [V, info] = algCrvVertices(quarters());
            P = sortrows(round(info.X, 9));
            t.verifyTrue(ismember([-1 0], P, 'rows'), 'the second face''s own corner');
            t.verifyEqual(size(unique(V, 'rows'), 1), size(V, 1), 'and no name appears twice');
        end

        function twoAXISParallelLINESMeet(t)
        % y = 0 has no x and x = 0 has no y, so neither variable can be eliminated between them --
        % the resultant needs both to have positive degree in it. They still meet, at the origin,
        % which is a corner of both quarter-disks.
            [~, info] = algCrvVertices(quarters());
            t.verifyTrue(ismember([0 0], round(info.X, 9), 'rows'), 'x = 0 meets y = 0');
        end

        function everyCORNEROfTheQuarterDisksIsNamedAndNothingElse(t)
        % The two faces' corners are exactly (1,0), (0,1), (-1,0) and (0,0); (0,-1), where x = 0
        % meets the circle, is on neither face.
            [~, info] = algCrvVertices(quarters());
            P = unique(round(info.X, 9), 'rows');
            t.verifyEqual(P, [-1 0; 0 0; 0 1; 1 0], 'AbsTol', 1e-9);
        end

        function aTANGENTCornerIsOnePointAndExact(t)
        % x + y = 2 touches the circle x^2 + y^2 = 2 at (1,1): a DOUBLE root of the resultant, which
        % a numeric root finder returns to about half the digits, as a pair of nearby roots. It is one
        % point, and a rational one, so it must come back as exactly (1,1) once.
            EcQ = [ratQ.conic([1 0 1 0 0 -2]); ratQ.conic([0 0 0 1 1 -2])];
            num = [0 0 0 0, 0, 0, 0, 0, 0, 1];
            g = AlgCrv(zeros(0,3), zeros(0,3), EcQ, num, 1, zeros(0,2), ...
                       {[1 -1 0 0; 2 -1 0 0]}, {[]}, 0);
            [~, info] = algCrvVertices(g);
            t.verifyEqual(unique(info.X, 'rows'), [1 1], 'AbsTol', 1e-12);
        end

        function theARGUMENTSAreChecked(t)
            t.verifyError(@() algCrvVertices(), 'algCrvVertices:nargin');
            t.verifyError(@() algCrvVertices(5), 'algCrvVertices:type');
        end
    end
end

% ==================================================================================================
function g = fixture(line)
% One face: inside the unit circle (the algebraic edge) and on the positive side of `line`.
    FP = zeros(3,3,3);
    FP(1,1,3) =  1;  FP(3,1,1) = -1;  FP(1,3,1) = -1;      % z^2 = x^2 + y^2
    num = [0 0 0 0, 0, 0, 0, 0, 0, 1];                     % the constant 1
    EcA = struct('FP', FP, 'FB', 2, 'num', num, 'den', 1);
    EcQ = ratQ.conic(line);
    FC  = {[1 1 0 1; 1 1 0 0]};                            % the curve, then the line
    g = AlgCrv(zeros(0,3), zeros(0,3), EcQ, num, 1, zeros(0,2), FC, {[]}, 0, EcA);
end

function g = quarters()
% The unit disk's upper half, split by x = 0: face 1 is the quarter with x >= 0, face 2 the other.
% Edges: 1 the circle (inside is <= 0), 2 the line y = 0, 3 the line x = 0. algCrvIncidenceTest
% uses the same object.
    EcQ = [ratQ.conic([1 0 1 0 0 -1]); ratQ.conic([0 0 0 0 1 0]); ratQ.conic([0 0 0 1 0 0])];
    num = [0 0 0 0, 0, 0, 0, 0, 0, 1];
    FC  = {[1 -1 0 0; 2 1 0 0; 3  1 0 0]; ...
           [1 -1 0 0; 2 1 0 0; 3 -1 0 0]};
    g = AlgCrv(zeros(0,3), zeros(0,3), EcQ, [num; num], [1; 1], zeros(0,2), FC, {[], []}, [0; 0]);
end
