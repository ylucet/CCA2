classdef symbolicFunctionUnitTest < matlab.unittest.TestCase
% Method-level unit tests for `symbolicFunction`, the rational-expression wrapper the whole
% symbolic pipeline is built on.
%
% WHY. 57 methods, and 35% of its lines under fast+normal on 2026-08-31. `testSymbolicFunction`
% covers the arithmetic and the tangent; what is left dark is the SUPPORTING vocabulary -- the
% predicates the region and conjugate code branch on (isLinear, isQuad, quadterm, isNegativeSqr),
% the accessors they read (getLinearCoeffs), the comparison operators, and the two
% routines that exist because MATLAB's own tools give the wrong answer here
% (removeDenominator2, limitDirectional).
%
% A wrong predicate here is the worst kind of defect this codebase can have: it does not throw, it
% silently selects a different branch several files away. That is exactly the shape a pipeline
% test cannot localise and a unit test can.
%
% BUCKET: fast. Small expressions; no `solve`, and the one `limit` call is on a rational in one
% direction.

    properties (Constant)
        X = sym('x')
        Y = sym('y')
    end

    methods (Static)
        function v = vars()
            v = [symbolicFunctionUnitTest.X, symbolicFunctionUnitTest.Y];
        end
    end

    methods (Test)

        % =========================================================================================
        % isLinear / isQuad -- the branch predicates, asserted as a partition of a family
        % =========================================================================================
        function degreePredicatesClassifyTheWholeFamily(testCase)
        % Every expression the pipeline can produce, with the answer each predicate owes it.
        %
        % TWO CLAUSES A CALLER FORGETS, both pinned here because both decide branches elsewhere:
        %   * a CONSTANT is neither linear nor quadratic. Both predicates test the degree for
        %     EXACT equality (degreeNum == 1, == 2), so a facet that has collapsed to a constant
        %     takes neither branch. Measured, not assumed -- the first version of this test read
        %     "linear" as "affine" and went red on 0 and 7.
        %   * a RATIONAL is neither, however simple its numerator, since both open by refusing a
        %     nonzero denominator degree.
            x = symbolicFunctionUnitTest.X; y = symbolicFunctionUnitTest.Y;
            cases = { sym(0),          false, false, 'zero: degree 0, so neither'
                      sym(7),          false, false, 'a nonzero constant: degree 0, so neither'
                      x,               true,  false, 'x'
                      2*x - 3*y + 4,   true,  false, 'a general affine form'
                      x^2,             false, true,  'x^2'
                      x*y,             false, true,  'x*y'
                      x^2 + y^2 - 1,   false, true,  'a circle'
                      x^3,             false, false, 'a cubic'
                      x/y,             false, false, 'a rational, numerator linear'
                      (x^2+1)/y,       false, false, 'a rational, numerator quadratic' };
            for k = 1:size(cases,1)
                f = symbolicFunction(cases{k,1});
                testCase.verifyEqual(logical(f.isLinear), cases{k,2}, sprintf( ...
                    'isLinear(%s) -- %s', char(cases{k,1}), cases{k,4}));
                testCase.verifyEqual(logical(f.isQuad), cases{k,3}, sprintf( ...
                    'isQuad(%s) -- %s', char(cases{k,1}), cases{k,4}));
                testCase.verifyFalse(logical(f.isLinear) && logical(f.isQuad), sprintf( ...
                    'isLinear and isQuad must be exclusive, and both are true for %s', ...
                    char(cases{k,1})));
            end
        end

        % =========================================================================================
        % getLinearCoeffs -- the coefficients callers rebuild the constraint row from
        % =========================================================================================
        function getLinearCoeffsReconstructsTheAffineForm(testCase)
        % region.slopeAtVertex and region.linearForm both read this and rebuild the form from it,
        % so the check is the reconstruction: c(1)x + c(2)y + c(3) must equal the original at
        % arbitrary points. Asserted that way rather than against a pinned triple, so the ORDER of
        % the coefficients is what is being tested and not the caller's memory of it.
            v = symbolicFunctionUnitTest.vars();
            exprs = { 2*v(1) - 3*v(2) + 4, v(1), v(2), sym(5), sym(0), -v(1) - v(2) };
            P = [0 0; 1 0; 0 1; -2 3; 0.5 -0.25];
            for k = 1:numel(exprs)
                f = symbolicFunction(exprs{k});
                c = double(f.getLinearCoeffs(v));
                testCase.verifyEqual(numel(c), 3, sprintf( ...
                    'getLinearCoeffs(%s) must return [x y const]', char(exprs{k})));
                for i = 1:size(P,1)
                    got  = c(1)*P(i,1) + c(2)*P(i,2) + c(3);
                    want = double(subs(exprs{k}, v, P(i,:)));
                    testCase.verifyEqual(got, want, 'AbsTol', 1e-12, sprintf( ...
                        'getLinearCoeffs(%s) at (%g,%g)', char(exprs{k}), P(i,1), P(i,2)));
                end
            end
        end

        % =========================================================================================
        % quadterm -- "does x appear squared here, and with what coefficient"
        % =========================================================================================
        function quadtermFindsTheSquareTermAndItsCoefficient(testCase)
            x = symbolicFunctionUnitTest.X; y = symbolicFunctionUnitTest.Y;
            cases = { 3*x^2 + x + 1,   x, true,  3,  'x^2 present'
                      x*y + y,         x, false, [], 'x appears, but not squared'
                      y^2 + 1,         x, false, [], 'the square is in the OTHER variable'
                      -2*x^2,          x, true,  -2, 'a negative coefficient'
                      y^2 + 1,         y, true,  1,  'the same expression, asked about y' };
            for k = 1:size(cases,1)
                f = symbolicFunction(cases{k,1});
                [l, c] = f.quadterm(cases{k,2});
                testCase.verifyEqual(logical(l), cases{k,3}, sprintf( ...
                    'quadterm(%s, %s) -- %s', char(cases{k,1}), char(cases{k,2}), cases{k,5}));
                if cases{k,3}
                    testCase.verifyEqual(double(c), cases{k,4}, 'AbsTol', 1e-12, sprintf( ...
                        'quadterm(%s, %s) coefficient', char(cases{k,1}), char(cases{k,2})));
                else
                    % The coefficient must come back EMPTY, not unassigned. Before 2026-08-31 it
                    % was assigned only on the true branch, so this two-output call RAISED
                    % "Output argument c ... not assigned a value" whenever the square was
                    % absent -- the routine could only be called safely on inputs already known
                    % to contain it. This assertion is why that was found.
                    testCase.verifyEmpty(c, sprintf( ...
                        'quadterm(%s, %s): with no square term the coefficient must be empty', ...
                        char(cases{k,1}), char(cases{k,2})));
                end
            end
        end

        % =========================================================================================
        % lt / gt against a number -- comparison operators are only useful if they separate
        % =========================================================================================
        function comparisonWithANumberSeparatesBothWays(testCase)
            a = symbolicFunction(sym(2));
            testCase.verifyTrue(a < 5,  '2 < 5');
            testCase.verifyFalse(a < 1, '2 < 1 is false');
            testCase.verifyTrue(a > 1,  '2 > 1');
            testCase.verifyFalse(a > 5, '2 > 5 is false');
            testCase.verifyFalse(a < 2, 'the comparison must be strict at equality');
            testCase.verifyFalse(a > 2, 'the comparison must be strict at equality');
        end

        % =========================================================================================
        % removeDenominator2 -- the rewrite must preserve the SIGN SET, which is all callers use
        % =========================================================================================
        function removeDenominatorPreservesTheSignWhereTheDenominatorIsPositive(testCase)
        % region calls this on a constraint, and a constraint is only ever read through its sign:
        % {g <= 0}. So the contract that matters is that the rewritten form has the same sign as
        % the original wherever the denominator does not vanish and is positive -- which is the
        % region of the plane the caller has already restricted itself to.
            v = symbolicFunctionUnitTest.vars();
            exprs = { (v(1) + 1)/(v(2)^2 + 1), (v(1)*v(2) - 2)/(v(1)^2 + 3), (v(1) - v(2))/2 };
            [gx, gy] = meshgrid(linspace(-3, 3, 41), linspace(-3, 3, 41));
            P = [gx(:), gy(:)];
            for k = 1:numel(exprs)
                f = symbolicFunction(exprs{k});
                % It answers with a bare SYM, not a symbolicFunction -- worth pinning, since a
                % caller that expected the wrapper would get a silent class change.
                g = f.removeDenominator2;
                testCase.verifyClass(g, 'sym', 'removeDenominator2 answers with a sym');
                hf = matlabFunction(f.f, 'Vars', {v(1), v(2)});
                hg = matlabFunction(g,   'Vars', {v(1), v(2)});
                bad = 0;
                for i = 1:size(P,1)
                    a = double(hf(P(i,1), P(i,2)));
                    b = double(hg(P(i,1), P(i,2)));
                    if ~isfinite(a) || ~isfinite(b), continue, end
                    if abs(a) < 1e-9 || abs(b) < 1e-9, continue, end   % on the zero set: no claim
                    if sign(a) ~= sign(b), bad = bad + 1; end
                end
                testCase.verifyEqual(bad, 0, sprintf( ...
                    ['removeDenominator2(%s) changed the SIGN at %d of %d sampled points -- ' ...
                     'every caller reads this only through {g <= 0}'], char(exprs{k}), bad, size(P,1)));
            end
        end


        % =========================================================================================
        % limitDirectional -- the bivariate 0/0 that MATLAB's own `limit` cannot resolve
        % =========================================================================================
        function limitDirectionalResolvesARemovableBivariateSingularity(testCase)
        % Its header: `limit` takes ITERATED univariate limits, which is the wrong tool for a
        % genuine bivariate 0/0 -- the inner limit is itself indeterminate and the iteration
        % returns NaN. This restricts to straight lines into the point and requires the values to
        % AGREE, which is what makes agreement a proof that the limit exists.
        %
        % The fixture is removable by construction: (x^2 - y^2)/(x - y) = x + y away from the
        % diagonal, so the limit at (1,1) is 2 along every direction. The answer is therefore
        % known independently of the routine.
            v = symbolicFunctionUnitTest.vars();
            f = symbolicFunction((v(1)^2 - v(2)^2) / (v(1) - v(2)));
            got = f.limitDirectional(v, [1 1]);
            testCase.verifyEqual(double(got.f), 2, 'AbsTol', 1e-9, ...
                '(x^2-y^2)/(x-y) simplifies to x+y, so its limit at (1,1) is 2');

            g = symbolicFunction((v(1)^2 - v(2)^2) / (v(1) - v(2)));
            got0 = g.limitDirectional(v, [0 0]);
            testCase.verifyEqual(double(got0.f), 0, 'AbsTol', 1e-9, ...
                'the same expression at the origin has limit 0');
        end

        function limitDirectionalRefusesADirectionDependentSingularity(testCase)
        % The other half of the contract, and the reason the routine compares directions at all:
        % x*y/(x^2 + y^2) at the origin has a DIFFERENT limit along every line (0 along the axes,
        % 1/2 along y = x), so no limit exists and the documented answer is NaN. A routine that
        % returned one arbitrary branch's value would look perfectly reasonable at any single
        % probe -- which is exactly why this is asserted.
            v = symbolicFunctionUnitTest.vars();
            f = symbolicFunction(v(1)*v(2) / (v(1)^2 + v(2)^2));
            got = f.limitDirectional(v, [0 0]);
            val = double(got.f);
            testCase.verifyTrue(isnan(val), sprintf( ...
                ['x*y/(x^2+y^2) has no limit at the origin -- 0 along the axes, 1/2 along ' ...
                 'y = x -- so the answer must be NaN, not %g'], val));
        end

        % =========================================================================================
        % dfdx and normalize1 -- small, and both are read by slopeAtVertex
        % =========================================================================================
        function dfdxIsThePartialDerivative(testCase)
            v = symbolicFunctionUnitTest.vars();
            exprs = { v(1)^2*v(2), v(1) + 3*v(2), sym(4), v(1)*v(2)^3 };
            P = [1 2; -1 0.5; 0 0; 2 -3];
            for k = 1:numel(exprs)
                f = symbolicFunction(exprs{k});
                for j = 1:2
                    d = f.dfdx(v(j));
                    want = diff(exprs{k}, v(j));
                    for i = 1:size(P,1)
                        testCase.verifyEqual(double(subs(d.f, v, P(i,:))), ...
                            double(subs(want, v, P(i,:))), 'AbsTol', 1e-12, sprintf( ...
                            'dfdx(%s, %s) at (%g,%g)', char(exprs{k}), char(v(j)), P(i,1), P(i,2)));
                    end
                end
            end
        end

        function normalize1KeepsTheZeroSetAndTheSide(testCase)
        % normalize1 rescales a constraint. A positive rescaling leaves {g <= 0} alone; a negative
        % one would invert it and silently flip a half-plane. So: same zero set, same side.
            v = symbolicFunctionUnitTest.vars();
            exprs = { 2*v(1) - 4*v(2) + 6, -3*v(1) + 9, v(1) + v(2) };
            [gx, gy] = meshgrid(linspace(-4, 4, 31), linspace(-4, 4, 31));
            P = [gx(:), gy(:)];
            for k = 1:numel(exprs)
                f = symbolicFunction(exprs{k});
                g = f.normalize1;
                hf = matlabFunction(f.f, 'Vars', {v(1), v(2)});
                hg = matlabFunction(sym(g.f), 'Vars', {v(1), v(2)});
                bad = 0;
                for i = 1:size(P,1)
                    a = double(hf(P(i,1), P(i,2)));
                    b = double(hg(P(i,1), P(i,2)));
                    if abs(a) < 1e-9, continue, end
                    if sign(a) ~= sign(b), bad = bad + 1; end
                end
                testCase.verifyEqual(bad, 0, sprintf( ...
                    ['normalize1(%s) flipped the sign at %d of %d points -- a negative rescaling ' ...
                     'turns the half-plane {g <= 0} into its complement'], ...
                    char(exprs{k}), bad, size(P,1)));
            end
        end
    end
end
