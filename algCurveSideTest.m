classdef algCurveSideTest < matlab.unittest.TestCase
% [fast] The edge a CONIC subdivision cannot hold: the boundary between a bridge face's value and a
% quadratic face's. Two routines, and the second exists because of what the first cannot do.
%
% `bridgeQuadCut` builds the rational curve R = FP(x, y, q(x,y)), which CONTAINS the boundary. But
% FP is degree 2 in z, so R vanishes wherever q meets EITHER of its roots: R names the PAIR of
% branches, exactly as a degenerate conic names a pair of lines. `algCurveSide` names ONE, by
% counting how many roots of the fibre lie at or below q -- a Sturm count over the integers, with no
% root ever computed and no irrational value ever formed.
%
% The fixture is the CONE z^2 = x^2 + y^2 against the constant face q = 1. Its upper branch is
% |x|, so the boundary is the unit circle and the two branches disagree everywhere outside it --
% which is what separates a correct branch test from R's own sign.

    methods (Test)
        function theCUTIsTheCurveWhereTheVALUESMeet(t)
            [FP, num, den] = fixture();
            [R, info] = bridgeQuadCut(FP, num, den);
            t.verifyEqual(info.zdeg, 2);
            % R = 1 - x^2 - y^2 up to content and sign: the unit circle
            [XX, YY] = meshgrid(linspace(-2, 2, 33));
            X = [XX(:), YY(:)];
            r = evalPoly2(R, X);
            onCircle = abs(hypot(X(:,1), X(:,2)) - 1) < 1e-12;
            t.verifyTrue(all(abs(r(onCircle)) < 1e-9), 'R must vanish where |x| = 1');
            inside = hypot(X(:,1), X(:,2)) < 0.99;
            outside = hypot(X(:,1), X(:,2)) > 1.01;
            t.verifyTrue(all(sign(r(inside)) == sign(r(find(inside,1)))), 'one sign inside');
            t.verifyTrue(all(sign(r(outside)) == sign(r(find(outside,1)))), 'one sign outside');
            t.verifyNotEqual(sign(r(find(inside,1))), sign(r(find(outside,1))), ...
                             'and the two sides differ');
        end

        function theSIDEIsPerBRANCHWhereTheCUTIsNot(t)
        % The distinction the whole design rests on. Against the UPPER branch (|x|) the constant 1
        % is larger inside the circle and smaller outside; against the LOWER branch (-|x|) it is
        % larger EVERYWHERE. R cannot tell those apart -- its sign is the product -- and the branch
        % test must.
            [FP, num, den] = fixture();
            X = [0 0; 0.5 0; 0 0.5; 0.25 0.25; 2 0; 0 2; 1.5 1.5; -3 0; -0.5 -0.5];
            inUnit = hypot(X(:,1), X(:,2)) < 1;

            upper = algCurveSide(FP, 2, num, den, X);       % value = +|x|
            t.verifyEqual(upper(inUnit), ones(sum(inUnit),1), ...
                          'inside the circle the constant 1 is above +|x|');
            t.verifyEqual(upper(~inUnit), -ones(sum(~inUnit),1), ...
                          'outside it is below');

            lower = algCurveSide(FP, 1, num, den, X);       % value = -|x|
            t.verifyEqual(lower, ones(size(X,1),1), ...
                          'the constant 1 is above -|x| everywhere, inside and out');
        end

        function aPointONTheCurveIsReportedAsZero(t)
            [FP, num, den] = fixture();
            % The DYADIC points of the unit circle are its axis points -- (3/5, 4/5) is on the
            % circle but is not a simple rational, and the exact path refuses such a point by name
            % rather than overflowing, which is asserted here too.
            X = [1 0; 0 1; -1 0; 0 -1];
            s = algCurveSide(FP, 2, num, den, X);
            t.verifyEqual(s, zeros(size(X,1),1), 'q IS the branch value there');
            t.verifyError(@() algCurveSide(FP, 2, num, den, [0.6 0.8]), ...
                          'algCurveSide:pointTooFine');
        end

        function theNEGATIVEBranchIndexCountsFromTheTop(t)
        % -1 is the largest real root, which is `algBranch`'s convention and `AlgCon.FB`'s.
            [FP, num, den] = fixture();
            X = [0.5 0; 2 0];
            t.verifyEqual(algCurveSide(FP, -1, num, den, X), ...
                          algCurveSide(FP,  2, num, den, X), ...
                          '-1 and 2 are the same branch on a two-root fibre');
        end

        function theARGUMENTSAreChecked(t)
            [FP, num, den] = fixture();
            t.verifyError(@() algCurveSide(FP, 2, num, den), 'algCurveSide:nargin');
            t.verifyError(@() algCurveSide(zeros(2,2,2), 2, num, den, [0 0]), ...
                          'algCurveSide:zDegree');
            t.verifyError(@() algCurveSide(FP, 0, num, den, [0 0]), 'algCurveSide:branch');
            t.verifyError(@() algCurveSide(FP, 2, num, den, [0 0 0]), 'algCurveSide:points');
            t.verifyError(@() bridgeQuadCut(FP, num), 'bridgeQuadCut:nargin');
            t.verifyError(@() bridgeQuadCut(FP, num, -1), 'bridgeQuadCut:denominator');
            t.verifyError(@() bridgeQuadCut(FP, [1 0 0 0 0 0 0 0 0 1], 1), 'bridgeQuadCut:cubic');
            t.verifyError(@() bridgeQuadCut(zeros(2,2,2), num, 1), 'bridgeQuadCut:zDegree');
        end
    end
end

% ==================================================================================================
function [FP, num, den] = fixture()
% The cone z^2 = x^2 + y^2, whose two branches are +|x| and -|x|, against the constant face q = 1.
    FP = zeros(3,3,3);
    FP(1,1,3) =  1;                 % z^2
    FP(3,1,1) = -1;                 % -x^2
    FP(1,3,1) = -1;                 % -y^2
    num = [0 0 0 0, 0, 0, 0, 0, 0, 1];      % the weighted basis: entry 10 is the constant
    den = 1;
end

function v = evalPoly2(R, X)
    v = zeros(size(X,1), 1);
    for i = 1:size(R,1)
        for j = 1:size(R,2)
            if R(i,j) ~= 0
                v = v + R(i,j) * X(:,1).^(i-1) .* X(:,2).^(j-1);
            end
        end
    end
end
