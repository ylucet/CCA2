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
