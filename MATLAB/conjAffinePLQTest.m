classdef conjAffinePLQTest < matlab.unittest.TestCase
% Tests for conjAffinePLQ, the direct conjugate of a piecewise-AFFINE function (TODO.md G2).
%
% Both fixtures are checked against a CLOSED FORM derived independently, not against today's
% output: an all-affine input's conjugate is a max of support functions, and for these two the max
% collapses to something nameable.

    methods (Test)

        function maxOfZeroXandYIsTheIndicatorOfTheSimplex(testCase)
        % f = max(0,x,y), G2's canonical case, and the reason this file exists. Its three faces are
        % CONES, so every one of them conjugates to a support function that is +inf off a cone --
        % which is exactly the shape `maxQuaPar.assertFullDomain` refuses, and why this input has
        % always gone to the symbolic Case C.
        %
        % f*(s) = sup_x <s,x> - max(0,x,y). Taking x along either axis shows it is +inf unless
        % s >= 0, and pushing along (1,1) shows it is +inf unless s1 + s2 <= 1; on the simplex the
        % sup is attained at the origin and is 0. So f* is the INDICATOR of the unit simplex.
            f = conjAffinePLQTest.maxZeroXY();
            g = conjAffinePLQ(f);
            testCase.verifyEqual(g.kind(), 'QuaPol', 'the conjugate of a PL function is PL');
            inside  = [0 0; 1 0; 0 1; 0.3 0.3; 0.5 0.5; 0.25 0.5];
            outside = [0.6 0.6; -0.1 0.5; 0.5 -0.1; 1 1; -1 -1; 2 0];
            for i = 1:size(inside,1)
                s = inside(i,:);
                testCase.verifyEqual(g.eval(s), 0, 'AbsTol', 1e-9, ...
                    sprintf('f*(%g,%g) must be 0 on the simplex', s(1), s(2)));
            end
            for i = 1:size(outside,1)
                s = outside(i,:);
                v = g.eval(s);
                testCase.verifyTrue(isinf(v) || isnan(v), sprintf( ...
                    'f*(%g,%g) must be +inf off the simplex, got %g', s(1), s(2), v));
            end
        end

        function theClippedInfinityNormGivesTheL1BallAndA4CellSUBDIVISION(testCase)
        % f = max(0, |x|-1, |y|-1). Five faces: the square [-1,1]^2 carrying 0, and four unbounded
        % faces each carrying an affine function and each having TWO vertices -- so unlike the
        % fixture above this one genuinely exercises the subdivision rather than collapsing to a
        % single cell.
        %
        % CLOSED FORM. <s,x> <= ||s||_1 ||x||_inf, and on ||x||_inf = t >= 1 the objective is
        % t*||s||_1 - (t-1) = t(||s||_1 - 1) + 1, which is bounded in t exactly when ||s||_1 <= 1
        % and is then largest at t = 1. So f*(s) = |s1| + |s2| on the L1 ball, +inf outside.
            f = conjAffinePLQTest.clippedInfNorm();
            g = conjAffinePLQ(f);
            testCase.verifyEqual(g.kind(), 'QuaPol');
            testCase.verifyEqual(g.nf, 4, ...
                'the L1 ball splits into one cell per sign pattern of (s1,s2)');
            % STRICTLY INSIDE their cells, and that is deliberate. A point lying exactly ON a
            % cell boundary asks `QuaPol.eval` to break a tie: it locates a point by requiring
            % every bounding conic to be <= tol, and on the boundary that quantity is zero up to
            % rounding. Measured here at s = (-0.5,0.5), which is exactly on the L1 ball's
            % second-quadrant edge: the binding constraint evaluates to +2.220e-16, so the point
            % is excluded, while its mirror image in the first quadrant comes out <= 0 and is
            % included. That is floating-point luck in the point location, not a statement about
            % f*, so asserting it would pin the tie-break rather than the mathematics.
            %
            % The geometry is still covered: the ring below sits 1e-6 INSIDE the boundary points
            % that would otherwise be used, which fails just as loudly if a cell is missing, wrong
            % or mis-oriented, and the +inf assertions below sit strictly outside.
            inside = [0 0; 0.25 0.25; 0.75 0.2; -0.3 -0.6; 0.2 -0.5; -0.6 0.1];
            edge   = [1 0; 0 1; -1 0; 0 -1; 0.5 0.5; -0.5 0.5; 0.5 -0.5; -0.5 -0.5];
            inside = [inside; edge * (1 - 1e-6)];
            for i = 1:size(inside,1)
                s = inside(i,:);
                testCase.verifyEqual(g.eval(s), abs(s(1)) + abs(s(2)), 'AbsTol', 1e-9, ...
                    sprintf('f*(%g,%g) must be |s1|+|s2|', s(1), s(2)));
            end
            for s = {[0.7 0.7], [1.2 0], [0 -1.4], [-0.9 -0.9]}
                v = g.eval(s{1});
                testCase.verifyTrue(isinf(v) || isnan(v), sprintf( ...
                    'f*(%g,%g) must be +inf off the L1 ball, got %g', s{1}(1), s{1}(2), v));
            end
        end

        function aBoundedDOMAINIsRefusedByNameRatherThanApproximated(testCase)
        % dom f bounded means f* is finite everywhere, so dom f* is unbounded -- outside this
        % construction's scope, and it says so instead of returning the clipping box it works in.
            f.W = [0 0; 1 0; 0 1]; f.dFirst = []; f.dLast = []; f.a = [1 2]; f.b = 0;
            testCase.verifyError(@() conjAffinePLQ(f), 'PLQ:conjAffinePLQ:unboundedDual');
        end
    end

    methods (Static)
        function f = maxZeroXY()
            f(1).W = [0 0]; f(1).dFirst = [-1 0]; f(1).dLast = [0 -1]; f(1).a = [0 0]; f(1).b = 0;
            f(2).W = [0 0]; f(2).dFirst = [1 1];  f(2).dLast = [1 0];  f(2).a = [1 0]; f(2).b = 0;
            f(3).W = [0 0]; f(3).dFirst = [0 1];  f(3).dLast = [1 1];  f(3).a = [0 1]; f(3).b = 0;
        end

        function f = clippedInfNorm()
        % max(0,|x|-1,|y|-1): the square, then the four strips where one coordinate dominates.
            f(1).W = [-1 -1; 1 -1; 1 1; -1 1]; f(1).dFirst = []; f(1).dLast = [];
            f(1).a = [0 0]; f(1).b = 0;
            f(2).W = [1 -1; 1 1];  f(2).dFirst = [1 -1]; f(2).dLast = [1 1];   % x >= 1, x >= |y|
            f(2).a = [1 0];  f(2).b = -1;
            f(3).W = [1 1; -1 1];  f(3).dFirst = [1 1];  f(3).dLast = [-1 1];  % y >= 1, y >= |x|
            f(3).a = [0 1];  f(3).b = -1;
            f(4).W = [-1 1; -1 -1]; f(4).dFirst = [-1 1]; f(4).dLast = [-1 -1];
            f(4).a = [-1 0]; f(4).b = -1;
            f(5).W = [-1 -1; 1 -1]; f(5).dFirst = [-1 -1]; f(5).dLast = [1 -1];
            f(5).a = [0 -1]; f(5).b = -1;
        end
    end
end
