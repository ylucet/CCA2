classdef algRowFailsTest < matlab.unittest.TestCase
% [fast] algRowFails: an algebraic row side * (q - z) >= 0 certified FALSE on a whole region, by
% Bernstein signs on sub-boxes -- what lets envelopeMeshAlg drop a face whose cell the curve
% R = FP(x, y, q) = 0 crosses (2026-10-29).
%
% BUCKET: fast (Bernstein bounds on small boxes; no symbolic call).
%
% The fixture is z^2 = x^2 + y^2 + 1 -- branches +-sqrt(x^2 + y^2 + 1), real everywhere -- against
% q = 2x. R = 3x^2 - y^2 - 1 vanishes on x = -sqrt((1 + y^2)/3), where q meets the LOWER branch,
% so R crosses every box below, and the old constant-sign certificate cannot decide any of them.

    methods (Test)
        function aRowFailingOnABoxTheCutCROSSESIsCertified(t)
            a = row(-1, 1, [0 0 0 0 0 0 0 2 0 0]);        % q = 2x >= +sqrt(.): false for x < 1/sqrt3
            box = [-3 0.3 -1 1];
            t.verifyEqual(polySignOn(bridgeQuadCut(a.FP, a.num, a.den), box, @(~) true, 6), 0, ...
                          'the fixture must have R = 0 crossing the box');
            t.verifyTrue(algRowFails(a, box, @(~) true));
        end

        function aRowHOLDINGSomewhereIsNotCertified(t)
            a = row(-1, 1, [0 0 0 0 0 0 0 2 0 0]);        % at (1, 0): 2 >= sqrt 2, the row holds
            t.verifyFalse(algRowFails(a, [-3 1 -1 1], @(~) true));
        end

        function theLOWERBranchAndTheOtherSide(t)
            a = row(1, -1, [0 0 0 0 0 0 0 2 0 0]);        % q = 2x <= -sqrt(.): false for x > -1/sqrt3
            t.verifyTrue(algRowFails(a, [-0.3 3 -1 1], @(~) true));
            t.verifyFalse(algRowFails(a, [-1 3 -1 1], @(~) true), 'at (-1, 0): -2 <= -sqrt 2');
        end

        function aSQUARETermIsReadInTheWeightedBasis(t)
        % num(5) = 2 is q = x^2 (the weighted basis halves it). Read unweighted, q = 2 x^2 > 1 at
        % x = 0.9 and the row would hold there.
            FP = zeros(1,3,3);  FP(1,1,3) = 1;  FP(1,3,1) = -1;  FP(1,1,1) = -1;   % z^2 = y^2 + 1
            a = struct('FP', FP, 'FB', -1, 'num', [0 0 0 0 2 0 0 0 0 0], 'den', 1, 'side', 1);
            t.verifyTrue(algRowFails(a, [-0.9 0.9 -1 1], @(~) true));
            t.verifyFalse(algRowFails(a, [-1.2 1.2 -1 1], @(~) true), 'at (1.2, 0): 1.44 >= 1');
        end

        function aLINEARFibre(t)
            FP = zeros(2,1,2);  FP(1,1,2) = 2;  FP(1,1,1) = -1;  FP(2,1,1) = -1;   % 2z = x + 1
            a = struct('FP', FP, 'FB', 1, 'num', [0 0 0 0 0 0 0 1 0 0], 'den', 1, 'side', 1);
            t.verifyTrue(algRowFails(a, [-2 0.9 0 1], @(~) true), 'x >= (x+1)/2 false for x < 1');
            t.verifyFalse(algRowFails(a, [-2 1.5 0 1], @(~) true));
        end

        function aBoxOutsideTheRegionIsNotExamined(t)
            a = row(-1, 1, [0 0 0 0 0 0 0 2 0 0]);
            inBox = @(b) b(1) < 0;                         % the region lies in x < 0
            t.verifyTrue(algRowFails(a, [-3 1 -1 1], inBox));
        end

        function aCERTIFIEDRowFailsAtEveryPointRandomInputs(t)
        % The definition, against the evaluator that shares no code with the certificate: wherever
        % algRowFails says true, algBranch and ratQ.evalFace find the row false at every sample.
            rng(29);  nCert = 0;  nRef = 0;
            for trial = 1:300
                FP = zeros(3,3,3);
                FP(1,1,3) = randi([1 3]) * sign(randn);
                FP(:,:,2) = randi([-2 2], 3) .* fliplr(triu(ones(3)));
                FP(:,:,1) = randi([-3 3], 3) .* fliplr(triu(ones(3)));
                num = [0 0 0 0, randi([-3 3], 1, 6)];
                a = struct('FP', FP, 'FB', subsref([-2 -1 1 2], substruct('()', {randi(4)})), 'num', num, ...
                           'den', randi(3), 'side', 2*randi(2) - 3);
                c = randn(1,2);  h = 0.05 + rand(1,2);
                box = [c(1) - h(1), c(1) + h(1), c(2) - h(2), c(2) + h(2)];
                if ~algRowFails(a, box, @(~) true, 5), continue, end
                nCert = nCert + 1;
                X = [box(1) + (box(2) - box(1)) * rand(400,1), box(3) + (box(4) - box(3)) * rand(400,1)];
                X = [X; box([1 3]); box([2 3]); box([1 4]); box([2 4])]; %#ok<AGROW>
                z = algBranch(a.FP, X, a.FB);  q = ratQ.evalFace(a.num, a.den, X);
                t.verifyTrue(all(isfinite(z)), sprintf('trial %d: the branch must exist on the box', trial));
                t.verifyTrue(all(a.side * (q - z) < 0), sprintf('trial %d: the row holds at a sample', trial));
                if polySignOn(bridgeQuadCut(a.FP, a.num, a.den), box, @(~) true, 5) == 0
                    nRef = nRef + 1;                       % one the old certificate could not give
                end
            end
            t.verifyGreaterThan(nCert, 20, 'too few certified trials to mean anything');
            t.verifyGreaterThan(nRef, 0, 'no certified trial had R = 0 crossing its box');
        end
    end
end

function a = row(FB, side, num)
    FP = zeros(3,3,3);
    FP(1,1,3) = 1;  FP(3,1,1) = -1;  FP(1,3,1) = -1;  FP(1,1,1) = -1;    % z^2 = x^2 + y^2 + 1
    a = struct('FP', FP, 'FB', FB, 'num', num, 'den', 1, 'side', side);
end
