classdef compactCellsTest < matlab.unittest.TestCase
% [fast] Does dropping IMPLIED rows preserve each face's region and expose the merges that
% `mergeAdjacentCells` alone cannot see? Both halves are asserted against the definition: a row may
% go only when the rest of the cell already lies on its side, and a merge may happen only when the
% union really is the cell that comes out.

    methods (Test)
        function anIMPLIEDRowIsDroppedAndTheREGIONIsUnchanged(t)
        % x >= 0 and x >= -1: the second is implied by the first, so it must go -- and no point may
        % change its membership, which is checked on a grid rather than argued.
            EcQ = [0 0 0 1 0 0;          % x
                   0 0 0 1 0 1;          % x + 1
                   0 0 0 0 1 0];         % y
            FC = {[1 1 0; 2 1 0; 3 1 0]};
            fN = [0 0 0 0 1 0 0 0 0 0];  fD = 1;
            [fN2, fD2, FC2, info] = compactCells(fN, fD, FC, EcQ);
            t.verifyEqual(info.dropped, 1, 'exactly the implied row goes');
            t.verifyEqual(sortrows(FC2{1}), sortrows([1 1 0; 3 1 0]));
            t.verifyEqual(fN2, fN);  t.verifyEqual(fD2, fD);
            [XX, YY] = meshgrid(linspace(-2, 2, 21));
            X = [XX(:), YY(:)];
            t.verifyEqual(holds(FC2{1}, EcQ, X), holds(FC{1}, EcQ, X), ...
                          'the region must be the same set, point by point');
        end

        function DroppingAnImpliedRowIsWhatLetsTwoCellsMERGE(t)
        % The two cells below carry one function and differ in the side of y AND in an extra row
        % that one of them does not need. `mergeAdjacentCells` sees two differences and refuses;
        % after the implied row goes it sees one, and the union {x >= 0} is exact.
            EcQ = [0 0 0 1 0 0;          % x
                   0 0 0 0 1 0;          % y
                   0 0 0 1 0 1];         % x + 1, implied by x >= 0
            FC = {[1 1 0; 2 1 0; 3 1 0], [1 1 0; 2 -1 0]};
            fN = repmat([0 0 0 0 1 0 0 0 0 0], 2, 1);  fD = [1; 1];
            [~, ~, FC2, info] = compactCells(fN, fD, FC, EcQ);
            t.verifyEqual(numel(FC2), 1, 'the two cells must become one');
            t.verifyEqual(sortrows(FC2{1}), [1 1 0], 'and the union is exactly x >= 0');
            t.verifyGreaterThan(info.merged, 0);
        end

        function aCellWhoseROWSAreAllNeededIsUntouched(t)
        % The other direction, without which the test above would pass on code that drops anything:
        % a triangle needs all three of its edges, and none may be removed.
            EcQ = [0 0 0 1 0 0; 0 0 0 0 1 0; 0 0 0 -1 -1 1];
            FC = {[1 1 0; 2 1 0; 3 1 0]};
            fN = [0 0 0 0 1 0 0 0 0 0];  fD = 1;
            [~, ~, FC2, info] = compactCells(fN, fD, FC, EcQ);
            t.verifyEqual(info.dropped, 0);
            t.verifyEqual(sortrows(FC2{1}), sortrows(FC{1}));
        end

        function aCURVEDRowCanBeDroppedToo(t)
        % What the older linear-only dropper could not do: inside the unit disc, the constraint
        % x^2 + y^2 <= 4 is implied, and deciding that needs the exact curved test.
            EcQ = [1 0 1 0 0 -1;         % x^2 + y^2 - 1
                   1 0 1 0 0 -4];        % x^2 + y^2 - 4
            FC = {[1 -1 0; 2 -1 0]};
            fN = [0 0 0 0 1 0 1 0 0 0];  fD = 1;
            [~, ~, FC2, info] = compactCells(fN, fD, FC, EcQ);
            t.verifyEqual(info.dropped, 1, 'the larger disc is implied by the smaller');
            t.verifyEqual(FC2{1}, [1 -1 0]);
        end

        function aTHINFaceKeepsEveryRow(t)
        % The red this guard exists for. `dom f*` is a LINE when the primal is a positive
        % semidefinite singular quadratic and a single POINT when it is affine, and QuaCon stores
        % both -- a cell with no interior at all. Every row's violation is then vacuously
        % interior-free, so the redundancy rule would strip such a face to a half-plane. Measured
        % before the guard: five fast-bucket reds, among them conjQ's single-point dual domain.
            EcQ = [0 0 0 1 0 0; 0 0 0 -1 0 0; 0 0 0 0 1 0; 0 0 0 0 -1 0];
            FC = {[1 1 0; 2 1 0; 3 1 0; 4 1 0]};      % x >= 0, x <= 0, y >= 0, y <= 0: the ORIGIN
            fN = [0 0 0 0 1 0 0 0 0 0];  fD = 1;
            [~, ~, FC2, info] = compactCells(fN, fD, FC, EcQ);
            t.verifyEqual(info.dropped, 0, 'a face with no interior may not lose a row');
            t.verifyEqual(sortrows(FC2{1}), sortrows(FC{1}));

            % ... and a LINE, the other thin shape, likewise
            FC = {[1 1 0; 2 1 0; 3 1 0]};             % x = 0 and y >= 0: a ray
            [~, ~, FC3, info] = compactCells(fN, fD, FC, EcQ);
            t.verifyEqual(info.dropped, 0);
            t.verifyEqual(sortrows(FC3{1}), sortrows(FC{1}));
        end

        function theARGUMENTSAreChecked(t)
            t.verifyError(@() compactCells(1, 2, 3), 'compactCells:nargin');
        end
    end
end

% ==================================================================================================
function tf = holds(rows, EcQ, X)
    tf = true(size(X,1), 1);
    for t = 1:size(rows,1)
        c = EcQ(rows(t,1), :);
        v = c(1)*X(:,1).^2 + c(2)*X(:,1).*X(:,2) + c(3)*X(:,2).^2 + ...
            c(4)*X(:,1) + c(5)*X(:,2) + c(6);
        tf = tf & (rows(t,2)*v >= 0);
    end
end
