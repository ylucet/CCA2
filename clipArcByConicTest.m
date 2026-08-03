classdef clipArcByConicTest < matlab.unittest.TestCase
    % Tests for clipArcByConic -- the arc-vs-arc primitive. One test per GEOMETRIC CONFIGURATION,
    % which is the whole point: the configurations are enumerable, so the tests should enumerate
    % them rather than sample a few arbitrary inputs.
    %
    % All fixtures use two parabolas written as explicit conics, so every expected crossing is
    % available in closed form and no test compares the routine against another routine:
    %   arc   : y = x^2                 ->  Ec   = [1 0 0  0 -1  0]
    %   cut A : y = (x-1)^2 + k         ->  Ecut = [1 0 0 -2 -1  1+k]   (same leading coeff)
    %   cut B : y = 2x^2 - 3            ->  Ecut = [2 0 0  0 -1 -3]     (different leading coeff)
    % The region KEPT is {evalConic(Ecut,.) >= 0}, i.e. BELOW the cutting parabola.
    %
    % Note why both cut families are needed. Two parabolas with the SAME leading coefficient differ
    % by an affine function, so they meet at most once -- those fixtures can produce inside /
    % outside / one crossing but never two. Two crossings require different leading coefficients.
    % That is the geometric reason the 'twice' configuration exists at all, and it is why a test
    % suite built from one cut family would silently miss it.

    methods (Static)
        function Ec = arcParabola()
            Ec = [1 0 0 0 -1 0];                 % y = x^2
        end
        function Ecut = cutSameLeading(k)
            Ecut = [1 0 0 -2 -1 1+k];            % y = (x-1)^2 + k
        end
        function Ecut = cutDifferentLeading()
            Ecut = [2 0 0 0 -1 -3];              % y = 2x^2 - 3
        end
        function X = onArc(x)
            X = [x, x.^2];
        end
    end

    methods (Test)

        function wholeArcInside(testCase)
            % CONFIGURATION 1: no crossing, the arc lies entirely in the kept region.
            % Cut y = (x-1)^2 + 5 meets y = x^2 where x^2 = x^2 - 2x + 6, i.e. x = 3 -- outside
            % the arc's own x-range [-1,2], so there is no crossing on this arc and the arc is
            % below the cutting parabola throughout.
            Ec = clipArcByConicTest.arcParabola();
            X0 = clipArcByConicTest.onArc(-1); X1 = clipArcByConicTest.onArc(2);
            [st, Xn] = clipArcByConic(Ec, X0, X1, clipArcByConicTest.cutSameLeading(5));
            testCase.verifyEqual(st, 'inside');
            testCase.verifyEqual(Xn, [X0; X1], 'AbsTol', 1e-12);
        end

        function wholeArcOutside(testCase)
            % CONFIGURATION 2: no crossing, the arc lies entirely OUTSIDE. Cut y = (x-1)^2 - 10
            % meets y = x^2 at x = -4.5, again off the arc, and the arc is above the cutting
            % parabola throughout, so nothing survives.
            Ec = clipArcByConicTest.arcParabola();
            X0 = clipArcByConicTest.onArc(-1); X1 = clipArcByConicTest.onArc(2);
            [st, Xn] = clipArcByConic(Ec, X0, X1, clipArcByConicTest.cutSameLeading(-10));
            testCase.verifyEqual(st, 'outside');
            testCase.verifyEmpty(Xn);
        end

        function singleCrossingReplacesOneEndpoint(testCase)
            % CONFIGURATION 3: exactly one crossing strictly inside the arc, so one endpoint is
            % replaced and the other is kept. Cut y = (x-1)^2 meets y = x^2 at 2x = 1, x = 0.5,
            % which lies strictly inside [-1,2]. Below the cut (kept) is x < 0.5, so the x = -1
            % end survives and the x = 2 end is replaced by (0.5, 0.25).
            Ec = clipArcByConicTest.arcParabola();
            X0 = clipArcByConicTest.onArc(-1); X1 = clipArcByConicTest.onArc(2);
            [st, Xn] = clipArcByConic(Ec, X0, X1, clipArcByConicTest.cutSameLeading(0));
            testCase.verifyEqual(st, 'cut');
            testCase.verifyEqual(Xn(1,:), X0, 'AbsTol', 1e-9);
            testCase.verifyEqual(Xn(2,:), [0.5, 0.25], 'AbsTol', 1e-9);
        end

        function singleCrossingRespectsEndpointOrder(testCase)
            % CONFIGURATION 3, endpoints supplied the OTHER way round. The frame parameter u is
            % intrinsic to the parabola and does not know which endpoint the caller called X0, so
            % the result must come back in the caller's own orientation. Getting this wrong is
            % silent: the arc is geometrically right and its two ends are swapped, which shows up
            % much later as a self-intersecting cell.
            Ec = clipArcByConicTest.arcParabola();
            X0 = clipArcByConicTest.onArc(2); X1 = clipArcByConicTest.onArc(-1);
            [st, Xn] = clipArcByConic(Ec, X0, X1, clipArcByConicTest.cutSameLeading(0));
            testCase.verifyEqual(st, 'cut');
            testCase.verifyEqual(Xn(1,:), [0.5, 0.25], 'AbsTol', 1e-9);
            testCase.verifyEqual(Xn(2,:), X1, 'AbsTol', 1e-9);
        end

        function twoCrossingsAreReportedNotGuessed(testCase)
            % CONFIGURATION 4: BOTH endpoints on the kept side with the arc leaving and re-entering
            % -- only possible when the two parabolas have different leading coefficients. Cut
            % y = 2x^2 - 3 against y = x^2 gives x^2 - 3 = 0, i.e. x = +/- sqrt(3), both strictly
            % inside the arc's range [-2,2], while the endpoints x = +/-2 give x^2-3 = 1 > 0.
            % The survivor is then two components, which one QuaPar face cannot hold, so this must
            % be reported rather than collapsed to a single cut.
            Ec = clipArcByConicTest.arcParabola();
            X0 = clipArcByConicTest.onArc(-2); X1 = clipArcByConicTest.onArc(2);
            [st, Xn] = clipArcByConic(Ec, X0, X1, clipArcByConicTest.cutDifferentLeading());
            testCase.verifyEqual(st, 'twice');
            testCase.verifyEqual(sort(Xn(:,1))', [-sqrt(3), sqrt(3)], 'AbsTol', 1e-9);
        end

        function tangencyIsNotACrossing(testCase)
            % A DOUBLE root is a tangency: the arc touches the cutting curve and returns to the
            % same side, so nothing is cut. y = 2x^2 meets y = x^2 only at x = 0, to second order.
            % Treating that as a crossing would split a cell along a curve it never crosses.
            Ec = clipArcByConicTest.arcParabola();
            Ecut = [2 0 0 0 -1 0];               % y = 2x^2, kept region y <= 2x^2
            X0 = clipArcByConicTest.onArc(-1); X1 = clipArcByConicTest.onArc(1);
            [st, ~] = clipArcByConic(Ec, X0, X1, Ecut);
            testCase.verifyEqual(st, 'inside');
        end

        function crossingPointsLieOnBothConics(testCase)
            % The invariant that must hold whatever the internals do: a reported crossing is a
            % point of BOTH curves. Checked on the two-crossing fixture, where the crossings are
            % genuinely interior to the arc.
            Ec   = clipArcByConicTest.arcParabola();
            Ecut = clipArcByConicTest.cutDifferentLeading();
            X0 = clipArcByConicTest.onArc(-2); X1 = clipArcByConicTest.onArc(2);
            [~, Xn] = clipArcByConic(Ec, X0, X1, Ecut);
            for i = 1:size(Xn,1)
                testCase.verifyEqual(QuaPar.evalConic(Ec,   Xn(i,:)), 0, 'AbsTol', 1e-9);
                testCase.verifyEqual(QuaPar.evalConic(Ecut, Xn(i,:)), 0, 'AbsTol', 1e-9);
            end
        end
    end
end
