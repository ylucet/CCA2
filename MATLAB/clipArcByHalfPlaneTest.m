classdef clipArcByHalfPlaneTest < matlab.unittest.TestCase
    % Standalone validation of clipArcByHalfPlane (the parabola-arc/half-plane clip primitive
    % maxQuaPar.m's curved-edge TODO needs -- see that file's header and DESIGN.md II.5.1/Phase 2).
    % Not yet exercised through maxQuaPar.m itself; these tests check the primitive in isolation,
    % against hand-derived and independently-constructed (rotated/shifted) parabola cases.

    methods (Test)
        function verticalClipCutsFarEndpoint(testCase)
            % y=x^2 from (-1,1) to (2,4); clip x<=1 cuts off the far endpoint at (1,1).
            [Ec, X0, X1] = clipArcByHalfPlaneTest.yEqualsXSquared();
            [status, Xnew] = clipArcByHalfPlane(Ec, X0, X1, [1 0], 1);
            testCase.verifyEqual(status, 'cut');
            testCase.verifyEqual(Xnew(1,:), X0, 'AbsTol', 1e-10);
            testCase.verifyEqual(Xnew(2,:), [1,1], 'AbsTol', 1e-8);
        end

        function horizontalClipCutsFarEndpoint(testCase)
            [Ec, X0, X1] = clipArcByHalfPlaneTest.yEqualsXSquared();
            [status, Xnew] = clipArcByHalfPlane(Ec, X0, X1, [0 1], 2);
            testCase.verifyEqual(status, 'cut');
            testCase.verifyEqual(Xnew(1,:), X0, 'AbsTol', 1e-10);
            testCase.verifyEqual(Xnew(2,:), [sqrt(2), 2], 'AbsTol', 1e-8);
        end

        function clipCutsNearEndpoint(testCase)
            % {x>=0} (i.e. -x<=0) cuts off the NEAR endpoint (-1,1), replacing it with (0,0).
            [Ec, X0, X1] = clipArcByHalfPlaneTest.yEqualsXSquared();
            [status, Xnew] = clipArcByHalfPlane(Ec, X0, X1, [-1 0], 0);
            testCase.verifyEqual(status, 'cut');
            testCase.verifyEqual(Xnew(1,:), [0,0], 'AbsTol', 1e-8);
            testCase.verifyEqual(Xnew(2,:), X1, 'AbsTol', 1e-10);
        end

        function fullyOutsideReturnsEmpty(testCase)
            [Ec, X0, X1] = clipArcByHalfPlaneTest.yEqualsXSquared();
            [status, Xnew] = clipArcByHalfPlane(Ec, X0, X1, [1 0], -5);
            testCase.verifyEqual(status, 'outside');
            testCase.verifyTrue(isempty(Xnew));
        end

        function fullyInsideReturnsBothEndpointsUnchanged(testCase)
            [Ec, X0, X1] = clipArcByHalfPlaneTest.yEqualsXSquared();
            [status, Xnew] = clipArcByHalfPlane(Ec, X0, X1, [1 0], 10);
            testCase.verifyEqual(status, 'inside');
            testCase.verifyEqual(Xnew, [X0; X1], 'AbsTol', 1e-10);
        end

        function rotatedShiftedParabolaMatchesRigidTransformOfAxisAlignedCase(testCase)
            % Rotate/shift y=x^2 by an arbitrary rigid transform (theta=pi/6, shift=(2,-1)); the
            % clip result must be the SAME rigid transform applied to the axis-aligned answer --
            % an independent check that the (u,v)-frame construction is not accidentally tied to
            % any axis-aligned special case.
            [Ec0, X0, X1] = clipArcByHalfPlaneTest.yEqualsXSquared();
            theta = pi/6; R = [cos(theta) -sin(theta); sin(theta) cos(theta)]; shift = [2,-1];
            Ec = clipArcByHalfPlaneTest.transformConic(Ec0, R, shift);
            X0r = (R*X0(:) + shift(:))'; X1r = (R*X1(:) + shift(:))';

            % Original clip: x<=1 (normal [1 0], through point (1,0)); transform rigidly too.
            nrmNew = (R*[1;0])';
            cNew = nrmNew*(R*[1;0] + shift(:));
            [status, Xnew] = clipArcByHalfPlane(Ec, X0r, X1r, nrmNew, cNew);

            testCase.verifyEqual(status, 'cut');
            expectedCut = (R*[1;1] + shift(:))';
            testCase.verifyEqual(Xnew(1,:), X0r, 'AbsTol', 1e-8);
            testCase.verifyEqual(Xnew(2,:), expectedCut, 'AbsTol', 1e-6);
        end

        function rejectsNonParabolicConic(testCase)
            % x^2+y^2-1=0 is a circle (delta=b^2-4ac=-4, not 0): must error, not silently proceed.
            Ec = [1 0 1 0 0 -1];
            testCase.verifyError(@() clipArcByHalfPlane(Ec, [1,0], [0,1], [1 0], 0.5), ...
                'clipArcByHalfPlane:notParabola');
        end
    end

    methods (Static)
        function [Ec, X0, X1] = yEqualsXSquared()
            % x^2 - y = 0  <=>  y = x^2, arc from (-1,1) to (2,4).
            Ec = [1 0 0 0 -1 0];
            X0 = [-1, 1]; X1 = [2, 4];
        end

        function EcNew = transformConic(Ec, R, shift)
            % Conic coefficients after the rigid map Y = R*X + shift (X = R'*(Y-shift)), by direct
            % symbolic substitution -- an independent derivation from clipArcByHalfPlane's own
            % (u,v)-frame construction, used here only as an external cross-check.
            syms Yx Yy real
            Rinv = R';
            XY = Rinv*([Yx;Yy] - shift(:));
            expr = expand(Ec(1)*XY(1)^2 + Ec(2)*XY(1)*XY(2) + Ec(3)*XY(2)^2 ...
                + Ec(4)*XY(1) + Ec(5)*XY(2) + Ec(6));
            a_ = double(subs(diff(expr,Yx,2),{Yx,Yy},{0,0}))/2;
            b_ = double(subs(diff(diff(expr,Yx),Yy),{Yx,Yy},{0,0}));
            cc_ = double(subs(diff(expr,Yy,2),{Yx,Yy},{0,0}))/2;
            d_ = double(subs(diff(expr,Yx),{Yx,Yy},{0,0}));
            e_ = double(subs(diff(expr,Yy),{Yx,Yy},{0,0}));
            f_ = double(subs(expr,{Yx,Yy},{0,0}));
            EcNew = [a_ b_ cc_ d_ e_ f_];
        end
    end
end
