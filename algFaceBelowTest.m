classdef algFaceBelowTest < matlab.unittest.TestCase
% [fast] algFaceBelow: one implicit face certified below another on a whole region -- what lets
% envelopeSubdivisionAlg drop an irrational face another one dominates on a cell (2026-11-02).
%
% BUCKET: fast (Bernstein bounds on small boxes; no symbolic call).

    methods (Test)
        function theLOWERBranchIsBelowTheUPPER(t)
            F = fp([1 1 1]);                                   % z^2 = x^2 + y^2 + 1
            t.verifyTrue(algFaceBelow(face(F, 1), face(F, -1), [-1 1 -1 1], @(~) true));
            t.verifyFalse(algFaceBelow(face(F, -1), face(F, 1), [-1 1 -1 1], @(~) true), 'reversed');
        end

        function aFaceABOVESomewhereIsNotCertified(t)
            lo = face(fp([1 1 4]), -1);  hi = face(fp([1 1 1]), -1);     % sqrt(r + 4) > sqrt(r + 1)
            t.verifyTrue(algFaceBelow(hi, lo, [-2 2 -2 2], @(~) true));
            t.verifyFalse(algFaceBelow(lo, hi, [-2 2 -2 2], @(~) true));
            % z = x against sqrt(x^2 + y^2 + 1): below everywhere, and the gap shrinks as x grows
            L = zeros(2,1,2);  L(1,1,2) = 1;  L(2,1,1) = -1;
            t.verifyTrue(algFaceBelow(face(L, 1), hi, [-3 3 -1 1], @(~) true));
            t.verifyFalse(algFaceBelow(hi, face(L, 1), [-3 3 -1 1], @(~) true));
        end

        function TOUCHINGFacesAreBelowOrEqual(t)
        % z = 2x - 1 <= z = x^2, equal at x = 1 only: the claim is <=, strict off a curve, and the
        % fitted separator (x^2 + 2x - 1)/2 meets both on x = 1 -- each row's squared factor. A face
        % CROSSING the other is never certified, near the crossing or not.
            lo = zeros(2,1,2);  lo(1,1,2) = 1;  lo(2,1,1) = -2;  lo(1,1,1) = 1;   % z - 2x + 1
            hi = zeros(3,1,2);  hi(1,1,2) = 1;  hi(3,1,1) = -1;                    % z - x^2
            t.verifyTrue(algFaceBelow(face(lo, 1), face(hi, 1), [0 2 -1 1], @(~) true, 4));
            t.verifyTrue(algFaceBelow(face(lo, 1), face(hi, 1), [-2 0.75 -1 1], @(~) true));
            x2 = zeros(2,1,2);  x2(1,1,2) = 4;  x2(2,1,1) = -8;  x2(1,1,1) = 3;   % z = 2x - 3/4
            t.verifyFalse(algFaceBelow(face(x2, 1), face(hi, 1), [0 2 -1 1], @(~) true, 4), ...
                          'z = 2x - 3/4 is above x^2 on (1/2, 3/2)');
            inBox = @(b) b(1) < 0.25;                           % boxes left of x = 1/4 in a larger one
            t.verifyTrue(algFaceBelow(face(x2, 1), face(hi, 1), [-2 2 -1 1], inBox), 'clear of the crossing');
        end

        function aBIGIntegerFaceAsItsDoubleTwin(t)
            K = 2^30 + 1;
            F = fp([1 1 1]);
            B = ratQ.bigMul(ratQ.bigMul(ratQ.toBig(F), ratQ.toBig(K)), ratQ.toBig(K));
            t.verifyTrue(algFaceBelow(face(B, 1), face(F, -1), [-1 1 -1 1], @(~) true));
            t.verifyTrue(algFaceBelow(face(F, 1), face(B, -1), [-1 1 -1 1], @(~) true));
            t.verifyFalse(algFaceBelow(face(B, -1), face(F, 1), [-1 1 -1 1], @(~) true));
        end

        function aCERTIFIEDPairIsOrderedAtEverySampleRandomInputs(t)
        % The definition: wherever algFaceBelow says true, z_lo <= z_hi at every sample of the box where
        % both branches exist.
            rng(43);  nCert = 0;  nNot = 0;  nBoth = 0;
            for trial = 1:150
                lo = face(randFP(), 2*randi(2) - 3);  hi = face(randFP(), 2*randi(2) - 3);
                c = randi([-4 4], 1, 2) / 2;  h = 0.25 + rand;
                box = [c(1) - h, c(1) + h, c(2) - h, c(2) + h];
                if ~algFaceBelow(lo, hi, box, @(~) true, 4), nNot = nNot + 1; continue, end
                nCert = nCert + 1;
                X = [box(1) + 2*h * rand(500,1), box(3) + 2*h * rand(500,1)];
                X = [X; box([1 3]); box([2 3]); box([1 4]); box([2 4])]; %#ok<AGROW>
                zl = algBranch(lo.FP, X, lo.FB);  zh = algBranch(hi.FP, X, hi.FB);
                both = isfinite(zl) & isfinite(zh);
                t.verifyTrue(all(zl(both) <= zh(both)), sprintf('trial %d: z_lo > z_hi at a sample', trial));
                nBoth = nBoth + sum(both);
            end
            t.verifyGreaterThan(nCert, 15, 'too few certified trials to mean anything');
            t.verifyGreaterThan(nNot, 15, 'and some must be refused, or the inputs are too easy');
            t.verifyGreaterThan(nBoth, 250 * nCert, 'and most of the 504 samples must have both branches');
        end
    end
end

% ==================================================================================================
function F = fp(c)
% z^2 = c1 x^2 + c2 y^2 + c3
    F = zeros(3,3,3);  F(1,1,3) = 1;  F(3,1,1) = -c(1);  F(1,3,1) = -c(2);  F(1,1,1) = -c(3);
end

function f = face(FP, FB)
    f = struct('FP', FP, 'FB', FB);
end

function F = randFP()
% z^2 + b z + c with b affine and c a quadratic, small integers -- branches over most of a box
    F = zeros(3,3,3);  F(1,1,3) = 1;
    F(1,1,2) = randi([-2 2]);  F(2,1,2) = randi([-1 1]);  F(1,2,2) = randi([-1 1]);
    F(1,1,1) = -randi([1 4]);  F(3,1,1) = randi([-1 0]);  F(1,3,1) = randi([-1 0]);  F(2,2,1) = randi([-1 1]);
end
