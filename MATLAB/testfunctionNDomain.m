classdef testfunctionNDomain < matlab.unittest.TestCase
% Definition checks for the two list operations `functionNDomain` provides.
%
% WHAT CHANGED AND WHY. Both tests here used to build a five-cell list, call the operation, and
% `printL` the result -- no assertion of any kind, in a suite that sits in the NORMAL bucket, i.e.
% in the pre-commit gate. A suite that can only go red by throwing reads as coverage and is not.
%
% Both operations have the same contract and it is the one asserted below: a `functionNDomain`
% array IS a piecewise function, so `mergeL` and `unique` may change the CELLS as much as they
% like and may not change the FUNCTION. Sample the domain, evaluate before and after, compare.
% Nothing here pins a cell count or a constraint list, so a better merge strategy -- which is an
% open item (TODO G4, the N-ary fan merge) -- makes these tests pass harder, not fail.

    properties
        r
    end

    methods (Static)
        function fL = fixtureList(fs)
        % The five-region fixture both tests were built on, with one function per region. `fs` is
        % a cell array of symbolic expressions, one per region, so a caller can give two cells the
        % same function (which is what makes a merge possible) or different ones (which must
        % prevent it).
            x = sym('x'); y = sym('y');
            ls = { [-x - 7*y - 4, x + 7*y-10, 196*y-148*x-(x+7*y)^2 + 684, 4-2*y-x ], ...
                   [(y/3)-x+ 14/3, 10-7*y-x, y-2], ...
                   [x+7*y+4, 4-2*y-x], ...
                   [-y/3+x-14/3, 10-7*y-x, y-2, 5*y/7-x+25/7], ...
                   [x + 7*y-10, -196*y+148*x+(x+7*y)^2 - 684, 4-2*y-x, x - y/3-14/3, 5*y/7-x+25/7] };
            fL = functionNDomain.empty();
            for k = 1:numel(ls)
                fL = [fL, functionNDomain(symbolicFunction(fs{min(k,numel(fs))}), ...
                                          region(ls{k},[x,y]))];   %#ok<AGROW>
            end
        end

        function verifySamePiecewiseFunction(tc, before, after, name)
        % `before` and `after` describe the same piecewise function: at every point that is
        % STRICTLY inside one of `before`'s cells, both lists must be defined and must agree.
        %
        % Strictly inside, because a point on a shared facet belongs to two cells and which one
        % `evalFunctionNDomain` reaches there is not part of any contract -- the values agree
        % anyway when the cells carry the same function, but the coverage flag need not.
            box = plqCheck.regionBox(arrayfun(@(z) z.d, before));
            P = zeros(0,2);
            for k = 1:numel(before)
                P = [P; plqCheck.regionSample(before(k).d, box, 25, 1e-6)];    %#ok<AGROW>
            end
            tc.verifyGreaterThan(size(P,1), 0, sprintf( ...
                '%s: no cell of the input has a sampled interior point -- the check is vacuous', name));
            nBad = 0; firstBad = [];
            for i = 1:size(P,1)
                [v1, ok1] = plqCheck.safeEval(before, P(i,:));
                [v2, ok2] = plqCheck.safeEval(after,  P(i,:));
                agree = (ok1 == ok2) && (~ok1 || abs(v1 - v2) <= 1e-9 * max(1, abs(v1)));
                if ~agree
                    nBad = nBad + 1;
                    if isempty(firstBad), firstBad = [P(i,:) v1 v2 ok1 ok2]; end
                end
            end
            if nBad == 0
                tc.verifyEqual(nBad, 0);            % records the pass with the count
            else
                tc.verifyEqual(nBad, 0, sprintf( ...
                    ['%s: the operation changed the piecewise function at %d of %d sampled ' ...
                     'points; first at (%g,%g) where before = %g (covered %d) and after = %g ' ...
                     '(covered %d)'], name, nBad, size(P,1), firstBad(1), firstBad(2), ...
                     firstBad(3), firstBad(5), firstBad(4), firstBad(6)));
            end
        end
    end

    methods(Test)

        function testMerge(testCase)
        % mergeL folds cells that carry the SAME function and whose union is exactly convex. All
        % five cells here carry `x`, so merging is permitted everywhere it is geometrically sound.
        %
        % Two assertions, and the second is the one that matters for TODO G4. mergeL may not
        % change the function (first), and it may not GROW the list (second) -- a fold that
        % returns more cells than it was given is not a fold. How many it removes is deliberately
        % not pinned: that number is the open research question, and a test that fixed it would
        % have to be edited by whoever improves the strategy, which is exactly the golden-value
        % trap.
          x = sym('x');
          fL = testfunctionNDomain.fixtureList({x, x, x, x, x});
          fL2 = fL.mergeL;

          testCase.verifyNotEmpty(fL2, 'mergeL returned an empty list');
          testCase.verifyLessThanOrEqual(numel(fL2), numel(fL), sprintf( ...
              'mergeL grew the list from %d cells to %d', numel(fL), numel(fL2)));
          testfunctionNDomain.verifySamePiecewiseFunction(testCase, fL, fL2, 'mergeL');

          % IDEMPOTENT ENOUGH: a second pass may still find something (maximumP itself runs two),
          % but it must not undo the first -- the function is still the same and the list still
          % does not grow. This is the property the "extra mergeL passes find nothing" experiment
          % of 2026-08-30 was measuring by hand.
          fL3 = fL2.mergeL;
          testCase.verifyLessThanOrEqual(numel(fL3), numel(fL2), ...
              'a second mergeL pass grew the list');
          testfunctionNDomain.verifySamePiecewiseFunction(testCase, fL, fL3, 'mergeL twice');

          % A CONTROL, so the test cannot pass by mergeL doing nothing at all: give the same five
          % regions five DIFFERENT functions and no two cells may be folded, because folding cells
          % that disagree would change the function.
          y = sym('y');
          fD = testfunctionNDomain.fixtureList({x, y, x+y, x-y, 2*x});
          fD2 = fD.mergeL;
          testCase.verifyEqual(numel(fD2), numel(fD), ...
              'mergeL folded cells carrying DIFFERENT functions, which cannot be sound');
        end

        function testUnique(testCase)
        % `unique` drops entries that repeat an existing (function, domain) pair exactly. The
        % contract is de-duplication, not merging: the point set and the values are untouched, the
        % duplicates go, and an entry that differs in EITHER coordinate stays.
            x = sym('x');
            fL = testfunctionNDomain.fixtureList({x, x, x, x, x});
            r1 = fL(1).d;

            % five distinct cells, then the first one twice more
            dup = [fL, functionNDomain(symbolicFunction(x), r1), ...
                       functionNDomain(symbolicFunction(x), r1)];
            testCase.verifyEqual(numel(dup), 7, 'fixture construction');

            u = dup.unique;
            testCase.verifyEqual(numel(u), numel(fL), sprintf( ...
                'unique left %d cells; the 7 entries hold %d distinct (f,domain) pairs', ...
                numel(u), numel(fL)));
            testfunctionNDomain.verifySamePiecewiseFunction(testCase, fL, u, 'unique');

            % SAME DOMAIN, DIFFERENT FUNCTION: not a duplicate. This is the half that a
            % de-duplicator keyed on the domain alone would get wrong, and it is silent -- it
            % deletes a piece of the function.
            u2 = [u, functionNDomain(symbolicFunction(x^2), r1)];
            u3 = u2.unique;
            testCase.verifyEqual(numel(u3), numel(u) + 1, ...
                'unique removed a cell that shares a domain but carries a different function');
        end
    end
end
