classdef algRowFailsTest < matlab.unittest.TestCase
% [fast] algRowFails: an algebraic row side * (q - z) >= 0 certified FALSE on a whole region, by
% Bernstein signs on sub-boxes -- what lets envelopeSubdivisionAlg drop a face whose cell the curve
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

        function aTOUCHINGRowHoldsOnlyOnACurve(t)
        % q = 0 against z = (x - y)^2 holds on the diagonal alone: F(q) = -(x - y)^2 has no strict
        % sign on any box the diagonal crosses, so this needs the squared factor divided out.
            L = zeros(3,3);  L(3,1) = 1;  L(2,2) = -2;  L(1,3) = 1;   % (x - y)^2
            FP = zeros(3,3,2);  FP(1,1,2) = 1;  FP(:,:,1) = -L;        % z = (x - y)^2
            a = struct('FP', FP, 'FB', 1, 'num', zeros(1,10), 'den', 1, 'side', 1);
            t.verifyTrue(algRowFails(a, [-1 1 -1 1], @(~) true), 'linear fibre');
            FP = zeros(3,3,3);  FP(1,1,3) = 1;  FP(:,:,2) = -L;  FP(1,1,2) = 1;  FP(:,:,1) = -L;
            a.FP = FP;  a.FB = -1;                                     % (z - (x-y)^2)(z + 1): upper
            t.verifyTrue(algRowFails(a, [-1 1 -1 1], @(~) true), 'quadratic fibre');
            a.num(10) = 1;  a.den = 100;                               % q = 1/100: a band, not a curve
            t.verifyFalse(algRowFails(a, [-1 1 -1 1], @(~) true), 'a row holding on a band');
        end

        function aRowVanishingOnTheREGIONSEdge(t)
        % q = 0 against z = 1 - x - y on the triangle x, y >= 0, x + y <= 1: q - z = -(1 - x - y)
        % fails inside and is ZERO on the hypotenuse, so every box the hypotenuse crosses has both
        % signs -- the 3-convex-edge triangle's corner faces, 2026-10-30. 1 - x - y >= 0 holds on the
        % region and divides F(q); divided out, the row fails everywhere but on that edge.
            ell = [1 -1; -1 0];                                        % 1 - x - y
            inBox = @(b) b(1) + b(3) <= 1;                             % meets x + y <= 1
            FP = zeros(2,2,2);  FP(1,1,2) = 1;  FP(:,:,1) = -ell;      % z = 1 - x - y
            a = struct('FP', FP, 'FB', 1, 'num', zeros(1,10), 'den', 1, 'side', 1);
            t.verifyFalse(algRowFails(a, [0 1 0 1], inBox, 6), 'without the factor: undecided');
            t.verifyTrue(algRowFails(a, [0 1 0 1], inBox, 6, {ell}), 'linear fibre');
            FP = zeros(2,2,3);  FP(1,1,3) = 1;  FP(:,:,2) = -ell;  FP(1,1,2) = FP(1,1,2) + 1;
            FP(:,:,1) = -ell;                                          % (z - (1-x-y))(z + 1): upper
            a.FP = FP;  a.FB = -1;
            t.verifyTrue(algRowFails(a, [0 1 0 1], inBox, 6, {ell}), 'quadratic fibre');
            a.side = -1;                                               % q <= z_hi holds inside
            t.verifyFalse(algRowFails(a, [0 1 0 1], inBox, 6, {ell}), 'a row holding inside');
            a.side = 1;  a.num(10) = -1;  a.den = 100;                 % q = -1/100 crosses z_hi inside
            t.verifyFalse(algRowFails(a, [0 1 0 1], inBox, 6, {ell}), 'a factor not dividing');
        end

        function aLEADINGCoefficientOrDISCRIMINANTVanishingAtACorner(t)
        % At the 3-convex-edge triangle's corner (0, 0), b of a linear fibre and D of a quadratic one
        % are multiples of the region's row x + y >= 0: zero at the corner, so no box containing it
        % gets their sign, and 228 faces stayed (2026-10-30). Divided out like F(q)'s factors.
            inBox = @(b) b(1) + b(3) <= 1;
            xy = [0 1; 1 0];                                           % x + y
            FP = zeros(2,2,2);  FP(:,:,2) = xy;  FP(1,1,1) = -1;  FP(2,1,1) = -1;   % (x+y) z = 1 + x
            a = struct('FP', FP, 'FB', 1, 'num', zeros(1,10), 'den', 1, 'side', 1);  % 0 >= z: false
            t.verifyFalse(algRowFails(a, [0 1 0 1], inBox, 6, {}), 'without the factor: undecided');
            t.verifyTrue(algRowFails(a, [0 1 0 1], inBox, 6, {xy}), 'b divided by x + y');
            FP = zeros(2,2,3);  FP(1,1,3) = 1;  FP(:,:,1) = -xy;       % z^2 = x + y, D = 4 (x + y)
            a = struct('FP', FP, 'FB', -1, 'num', [0 0 0 0 0 0 0 0 0 -1], 'den', 1, 'side', 1);  % -1 >= +sqrt
            % D's sign is no longer read for this conclusion (q = -1 below the roots' mean 0 puts it
            % below the upper root wherever the roots exist; 2026-11-02), so the factor is not needed
            t.verifyTrue(algRowFails(a, [0 1 0 1], inBox, 6, {}), 'D is not needed');
            t.verifyTrue(algRowFails(a, [0 1 0 1], inBox, 6, {xy}), 'and dividing it out does no harm');
            a.side = -1;                                               % -1 <= +sqrt holds
            t.verifyFalse(algRowFails(a, [0 1 0 1], inBox, 6, {xy}), 'a row holding');
        end

        function aRowZEROAtACorner(t)
        % c = (x - y) + y^2 > 0 on the wedge 0 <= y <= x/2 except at the origin, but every box
        % containing the origin also holds (0, h) where c < 0: no depth decides it (20 faces at the
        % 3-convex-edge triangle's corner, 2026-10-30). Blown up, y = s t x-major: c / s = 1 - t + s t^2 > 0
        % for t <= 1/2.
            FP = zeros(2,3,2);  FP(1,1,2) = -1;                         % -z + c = 0: z = c
            FP(2,1,1) = 1;  FP(1,2,1) = -1;  FP(1,3,1) = 1;
            a = struct('FP', FP, 'FB', 1, 'num', zeros(1,10), 'den', 1, 'side', 1);   % 0 >= c
            fac = {[0 1], [0 -2; 1 0]};                                % y, x - 2y
            inBox = @(b) ~any(cellfun(@(P) polySignOn(P, b, @(~) true, 0) == -1, fac));
            t.verifyFalse(algRowFails(a, [0 1 0 1], inBox, 6, fac), 'without the corner: undecided');
            t.verifyTrue(algRowFails(a, [0 1 0 1], inBox, 6, fac, [0 0]), 'blown up at the corner');
            t.verifyFalse(algRowFails(a, [0 1 0 1], inBox, 6, fac, [0.5 0]), 'a corner elsewhere');
            fac2 = {[0 1], [0 -1; 1 0]};                               % 0 <= y <= x: c's zero curve
            inBox2 = @(b) ~any(cellfun(@(P) polySignOn(P, b, @(~) true, 0) == -1, fac2));
            a.FP(1,3,1) = -1;                                          % c = x - y - y^2 < 0 near y = x
            t.verifyFalse(algRowFails(a, [0 1 0 1], inBox2, 6, fac2, [0 0]), 'a row holding in the wedge');
        end

        function aCornerChartOnTheNEGATIVESide(t)
        % c = x + 32 (x^2 + y^2) <= 0 on the disk of radius 1/64 left of the origin, inside the double
        % wedge x^2 >= 4 y^2: the row 0 >= c HOLDS there, but only within the leaf next to the
        % corner, so only the x < 0 chart sees it. A chart reading P(p + s (1, t)) for both signs
        % (no sg) certified it.
            FP = zeros(3,3,2);  FP(1,1,2) = -1;  FP(2,1,1) = 1;  FP(3,1,1) = 32;  FP(1,3,1) = 32;
            a = struct('FP', FP, 'FB', 1, 'num', zeros(1,10), 'den', 1, 'side', 1);
            W = zeros(3,3);  W(3,1) = 1;  W(1,3) = -4;                 % x^2 - 4 y^2
            inBox = @(b) polySignOn(W, b, @(~) true, 0) ~= -1;
            t.verifyFalse(algRowFails(a, [-1 1 -1 1], inBox, 4, {W}, [0 0]), 'the row holds at (-1/64, 0)');
            a.FP(2,1,1) = -1;                                          % mirrored: holds right of 0
            t.verifyFalse(algRowFails(a, [-1 1 -1 1], inBox, 4, {W}, [0 0]), 'mirrored');
        end

        function aRowZEROAtAnIRRATIONALCorner(t)
        % The 3-convex-edge triangle's last 4 faces (2026-10-30 (third)) meet at (1/sqrt5, 1 - 1/sqrt5),
        % where x + y = 1 meets 5 x^2 = 1. Here c = 5 x^2 + x + y - 2 is < 0 on the region
        % {x + y <= 1, 5 x^2 <= 1} near that corner except AT it, and every box containing it holds
        % points outside where c > 0: the row c >= 0 fails on the region only via a blow-up, which
        % at an irrational corner must be done over Q(sqrt5).
            FP = zeros(3,2,2);  FP(1,1,2) = -1;                         % -z + c = 0: z = c
            FP(3,1,1) = 5;  FP(2,1,1) = 1;  FP(1,2,1) = 1;  FP(1,1,1) = -2;
            a = struct('FP', FP, 'FB', 1, 'num', zeros(1,10), 'den', 1, 'side', -1);  % c >= 0
            fac = {[1 -1; -1 0], [1; 0; -5]};                          % 1 - x - y, 1 - 5 x^2
            inBox = @(b) ~any(cellfun(@(P) polySignOn(P, b, @(~) true, 0) == -1, fac));
            p = [1/sqrt(5), 1 - 1/sqrt(5)];
            box = [p(1) - 0.1, p(1) + 0.1, p(2) - 0.1, p(2) + 0.1];
            t.verifyFalse(algRowFails(a, box, inBox, 6, fac), 'without the corner: undecided');
            t.verifyTrue(algRowFails(a, box, inBox, 6, fac, p), 'blown up over Q(sqrt5)');
            t.verifyFalse(algRowFails(a, box, inBox, 6, fac, p + [1e-6 0]), 'a corner beside it');
            % the row holding instead: -c >= 0 is TRUE on the region, never certified as failing
            a.side = 1;
            t.verifyFalse(algRowFails(a, box, inBox, 6, fac, p), 'a row that holds');
            % the definition, sampled: c < 0 on the region in the box, the corner excepted
            rng(7);  X = [box(1) + 0.2*rand(20000,1), box(3) + 0.2*rand(20000,1); p + 1e-5*randn(2000,2)];
            X = X(X(:,1) + X(:,2) <= 1 & 5*X(:,1).^2 <= 1, :);
            t.verifyLessThan(max(5*X(:,1).^2 + X(:,1) + X(:,2) - 2), 0);
        end

        function aCERTIFIEDRowAtAnIRRATIONALCornerFailsRandomInputs(t)
        % The oracle for the Q(sqrt d) charts: rows vanishing at p = (a + b sqrt d) / q, the region
        % the wedge of two lines through p with coefficients in Q (both a multiple of one surd line
        % b2 = -b1, so x + y = (a1 + a2) / q is rational) and the quadratic (x - p1)(x - p1').
        % Wherever algRowFails with the corner says true, no sample of the region holds every row.
            rng(33);  nCert = 0;  nNeed = 0;
            for trial = 1:400
                d = randsample([2 3 5 7], 1);  qd = randi(4);
                ab = randi([-3 3], 2, 2);  ab(:,2) = ab(1,2) * [1; -1];  if ab(1,2) == 0, ab(:,2) = [1; -1]; end
                p = ((ab(:,1) + ab(:,2) * sqrt(d)) / qd).';
                % (x - p1)(x - p1') qd^2 = qd^2 x^2 - 2 a1 qd x + a1^2 - d b1^2
                s1 = 2*randi(2) - 3;  s2 = 2*randi(2) - 3;
                L1 = s1 * [-(ab(1,1) + ab(2,1)) qd; qd 0];
                L2 = zeros(3,1);  L2(3) = qd^2;  L2(2) = -2*ab(1,1)*qd;  L2(1) = ab(1,1)^2 - d*ab(1,2)^2;  L2 = s2 * L2;
                fac = {L1, L2};
                rows = struct('FP', {}, 'FB', {}, 'num', {}, 'den', {}, 'side', {});
                for r = 1:2
                    c = randi([-3 3]) * [L1; 0 0] + randi([-3 3]) * [L2, zeros(3,1)];  % zero at p
                    FP = zeros(3,2,2);  FP(1,1,2) = 2*randi(2) - 3;  FP(:,:,1) = -FP(1,1,2) * c;
                    rows(r) = struct('FP', FP, 'FB', 1, 'num', zeros(1,10), 'den', 1, 'side', 2*randi(2) - 3);
                end
                h = 0.25 + rand / 2;
                box = [p(1) - h, p(1) + h, p(2) - h, p(2) + h];
                inBox = @(b) ~any(cellfun(@(P) polySignOn(P, b, @(~) true, 0) == -1, fac));
                if ~algRowFails(rows, box, inBox, 4, fac, p), continue, end
                nCert = nCert + 1;
                X = [box(1) + 2*h * rand(3000,1), box(3) + 2*h * rand(3000,1); p + 1e-4 * randn(500, 2)];
                ev = @(P, X) polyval2(P, X);
                X = X(ev(L1, X) >= 0 & ev(L2, X) >= 0 & any(abs(X - p) > 1e-12, 2), :);
                holds = true(size(X,1), 1);
                for r = 1:2
                    z = algBranch(rows(r).FP, X, rows(r).FB);
                    holds = holds & rows(r).side * (0 - z) >= 0;
                end
                t.verifyFalse(any(holds), sprintf('trial %d: every row holds at a sample of the region', trial));
                if ~algRowFails(rows, box, inBox, 4, fac), nNeed = nNeed + 1; end
            end
            t.verifyGreaterThan(nCert, 20, 'too few certified trials to mean anything');
            t.verifyGreaterThan(nNeed, 5, 'too few certified trials needed the corner');
        end

        function aCERTIFIEDRowAtACornerFailsOnTheWedgeRandomInputs(t)
        % The definition at a corner: rows vanishing at a random dyadic corner p, the region a random
        % wedge at p (two lines through p as factors); wherever algRowFails with the corner says
        % true, every sample of the wedge in the box fails some row (algBranch / ratQ.evalFace).
            rng(32);  nCert = 0;  nNeed = 0;
            for trial = 1:1000
                p = randi([-8 8], 1, 2) / 4;
                L = cell(1,2);
                for k = 1:2
                    n = randi([-3 3], 1, 2);  if ~any(n), n = [1 0]; end
                    L{k} = [-n * p.', n(2); n(1), 0];                  % n . (x - p) >= 0
                end
                rows = struct('FP', {}, 'FB', {}, 'num', {}, 'den', {}, 'side', {});
                for r = 1:2
                    C = zeros(3,3);                                    % c(p) = 0
                    C(2,1) = randi([-3 3]);  C(1,2) = randi([-3 3]);
                    C(3,1) = randi([-2 2]);  C(2,2) = randi([-2 2]);  C(1,3) = randi([-2 2]);
                    C(1,1) = -(C(2,1)*p(1) + C(1,2)*p(2) + C(3,1)*p(1)^2 + C(2,2)*p(1)*p(2) + C(1,3)*p(2)^2);
                    FP = zeros(3,3,2);  FP(1,1,2) = 2*randi(2) - 3;  FP(:,:,1) = -FP(1,1,2) * C;  % z = c
                    rows(r) = struct('FP', FP, 'FB', 1, 'num', zeros(1,10), 'den', 1, 'side', 2*randi(2) - 3);
                end
                h = 0.25 + rand / 2;
                box = [p(1) - h, p(1) + h, p(2) - h, p(2) + h];
                inBox = @(b) ~any(cellfun(@(P) polySignOn(P, b, @(~) true, 0) == -1, L));
                if ~algRowFails(rows, box, inBox, 4, L, p), continue, end
                nCert = nCert + 1;
                X = [box(1) + 2*h * rand(3000,1), box(3) + 2*h * rand(3000,1)];
                X = [X; p + 1e-4 * randn(500, 2)]; %#ok<AGROW>
                ev = @(P, X) P(1,1) + P(2,1) * X(:,1) + P(1,2) * X(:,2);
                X = X(ev(L{1}, X) >= 0 & ev(L{2}, X) >= 0 & any(X ~= p, 2), :);
                holds = true(size(X,1), 1);
                for r = 1:2
                    z = algBranch(rows(r).FP, X, rows(r).FB);  q = ratQ.evalFace(rows(r).num, rows(r).den, X);
                    holds = holds & rows(r).side * (q - z) >= 0;
                end
                t.verifyFalse(any(holds), sprintf('trial %d: every row holds at a sample of the wedge', trial));
                if ~algRowFails(rows, box, inBox, 4, L), nNeed = nNeed + 1; end
            end
            t.verifyGreaterThan(nCert, 20, 'too few certified trials to mean anything');
            t.verifyGreaterThan(nNeed, 5, 'too few certified trials needed the corner');
        end

        function aCERTIFIEDRowWithAFactorFailsOffTheEdgeRandomInputs(t)
        % The definition with a factor: on the region ell >= 0, wherever algRowFails says true, the
        % row is false at every sample with ell > 0 (algBranch / ratQ.evalFace). F(q) is built to
        % carry ell, so the factor actually divides.
            rng(31);  nCert = 0;  nNeed = 0;
            for trial = 1:300
                ell = zeros(2,2);  ell(1,1) = randi([-2 2]);  ell(2,1) = randi([-3 3]);  ell(1,2) = randi([-3 3]);
                if ~any(ell(2:end)), continue, end
                m = zeros(2,2);  m(1:3) = randi([-3 3], 1, 3);
                Cq = conv2(ell, m);                                    % F(0) = ell * m
                if randi(2) == 1
                    FP = zeros(3,3,2);  FP(1,1,2) = randi([1 3]) * sign(randn);
                    FP(1:2,1:2,1) = Cq(1:2,1:2);  FP(3,1,1) = Cq(3,1);  FP(1,3,1) = Cq(1,3);
                    FB = 2*randi(2) - 3;
                else
                    FP = zeros(3,3,3);  FP(1,1,3) = randi([1 3]) * sign(randn);
                    FP(:,:,2) = randi([-2 2], 3) .* fliplr(triu(ones(3)));
                    FP(1:2,1:2,1) = Cq(1:2,1:2);  FP(3,1,1) = Cq(3,1);  FP(1,3,1) = Cq(1,3);
                    br = [-2 -1 1 2];  FB = br(randi(4));
                end
                a = struct('FP', FP, 'FB', FB, 'num', zeros(1,10), 'den', 1, 'side', 2*randi(2) - 3);
                c = randn(1,2);  h = 0.05 + rand(1,2);
                box = [c(1) - h(1), c(1) + h(1), c(2) - h(2), c(2) + h(2)];
                ev = @(X) ell(1,1) + ell(2,1) * X(:,1) + ell(1,2) * X(:,2);
                inBox = @(b) max(ev([b(1) b(3); b(2) b(3); b(1) b(4); b(2) b(4)])) >= 0;
                if ~inBox(box) || ~algRowFails(a, box, inBox, 5, {ell}), continue, end
                nCert = nCert + 1;
                X = [box(1) + (box(2) - box(1)) * rand(2000,1), box(3) + (box(4) - box(3)) * rand(2000,1)];
                X = X(ev(X) > 0, :);
                z = algBranch(a.FP, X, a.FB);  q = ratQ.evalFace(a.num, a.den, X);
                ok = isfinite(z);
                t.verifyTrue(all(a.side * (q(ok) - z(ok)) < 0), ...
                             sprintf('trial %d: the row holds at a sample inside the region', trial));
                if ~algRowFails(a, box, inBox, 5), nNeed = nNeed + 1; end
            end
            t.verifyGreaterThan(nCert, 20, 'too few certified trials to mean anything');
            t.verifyGreaterThan(nNeed, 0, 'no certified trial needed the factor');
        end

        function twoRowsEmptyOnlyJOINTLY(t)
        % q = 2x >= +sqrt(.) holds for x > 1/sqrt3, q = 2x <= -sqrt(.) for x < -1/sqrt3: each holds
        % on part of the box, never both at one point.
            a = [row(-1, 1, [0 0 0 0 0 0 0 2 0 0]), row(1, -1, [0 0 0 0 0 0 0 2 0 0])];
            box = [-3 3 -1 1];
            t.verifyFalse(algRowFails(a(1), box, @(~) true), 'the first row alone holds somewhere');
            t.verifyFalse(algRowFails(a(2), box, @(~) true), 'the second too');
            t.verifyTrue(algRowFails(a, box, @(~) true), 'together they hold nowhere');
        end

        function aCERTIFIEDPairOfRowsNeverHoldsTogetherRandomInputs(t)
        % The definition for several rows: wherever algRowFails says true, no sample has every row
        % holding (algBranch and ratQ.evalFace, sharing no code with the certificate).
            rng(30);  nCert = 0;  nJoint = 0;
            for trial = 1:300
                a = [randRow(), randRow()];
                c = randn(1,2);  h = 0.05 + rand(1,2);
                box = [c(1) - h(1), c(1) + h(1), c(2) - h(2), c(2) + h(2)];
                if ~algRowFails(a, box, @(~) true, 5), continue, end
                nCert = nCert + 1;
                X = [box(1) + (box(2) - box(1)) * rand(400,1), box(3) + (box(4) - box(3)) * rand(400,1)];
                holds = true(size(X,1), 1);
                for k = 1:2
                    z = algBranch(a(k).FP, X, a(k).FB);  q = ratQ.evalFace(a(k).num, a(k).den, X);
                    holds = holds & isfinite(z) & a(k).side * (q - z) >= 0;
                end
                t.verifyFalse(any(holds), sprintf('trial %d: both rows hold at a sample', trial));
                if ~algRowFails(a(1), box, @(~) true, 5) && ~algRowFails(a(2), box, @(~) true, 5)
                    nJoint = nJoint + 1;                   % neither row alone was enough
                end
            end
            t.verifyGreaterThan(nCert, 20, 'too few certified trials to mean anything');
            t.verifyGreaterThan(nJoint, 0, 'no trial needed both rows');
        end

        function aCERTIFIEDRowFailsAtEveryPointRandomInputs(t)
        % The definition, against the evaluator that shares no code with the certificate: wherever
        % algRowFails says true, algBranch and ratQ.evalFace find the row false at every sample.
            rng(29);  nCert = 0;  nRef = 0;  nFin = 0;
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
                % the row as the evaluator reads it (subdivisionAlgEval's algHolds): a point with no
                % branch fails it. (Until 2026-11-02 a certificate also implied the branch EXISTS on
                % the box; no caller needs that, and requiring it left every box on a discriminant
                % curve undecided.)
                holds = isfinite(z) & a.side * (q - z) >= 0;
                t.verifyFalse(any(holds), sprintf('trial %d: the row holds at a sample', trial));
                nFin = nFin + sum(isfinite(z));
                if polySignOn(bridgeQuadCut(a.FP, a.num, a.den), box, @(~) true, 5) == 0
                    nRef = nRef + 1;                       % one the old certificate could not give
                end
            end
            t.verifyGreaterThan(nCert, 20, 'too few certified trials to mean anything');
            t.verifyGreaterThan(nRef, 0, 'no certified trial had R = 0 crossing its box');
            t.verifyGreaterThan(nFin, 200 * nCert, 'most samples must have the branch, or this tests nothing');
        end

        function aBIGIntegerFibreIsCertifiedAsItsDoubleTwin(t)
        % Erratum Example 2's bridges carry FP past 2^53 (ratQ's big form). rowPolys called
        % ratQ.fromBig, which raises there, and envelopeSubdivisionAlg's caller turned the raise into
        % "not certified" -- so no algebraic row of those pairs was ever decided (2026-11-02). The
        % sign polynomials are formed in BigInteger instead. Two inputs: z^2 = x^2 + y^2 + 1 scaled by
        % K^2 (content K^2, reduced back exactly), and the same with the constant moved by one unit
        % (no content: the reduced polynomial is past 2^53 too) -- its branch differs by ~1e-37.
            K = 2^30 + 1;
            a = row(-1, 1, [0 0 0 0 0 0 0 2 0 0]);
            box = [-3 0.3 -1 1];
            t.verifyTrue(algRowFails(a, box, @(~) true), 'the double twin must be certified');
            B = ratQ.bigMul(ratQ.bigMul(ratQ.toBig(a.FP), ratQ.toBig(K)), ratQ.toBig(K));
            t.verifyTrue(ratQ.isBigPoly(B) && ratQ.bigBits(B.c{1}) > 53, 'the fixture must be past 2^53');
            a.FP = B;
            t.verifyTrue(algRowFails(a, box, @(~) true), 'scaled by K^2');
            B.c{1} = B.c{1}.subtract(java.math.BigInteger.ONE);
            a.FP = B;
            t.verifyTrue(algRowFails(a, box, @(~) true), 'scaled by K^2, constant moved by one');
            a = row(-1, 1, [0 0 0 0 0 0 0 2 0 0]);  a.FP = B;
            t.verifyFalse(algRowFails(a, [-3 1 -1 1], @(~) true), 'and a row holding at (1, 0) is not');
        end
    end
end

function a = row(FB, side, num)
    FP = zeros(3,3,3);
    FP(1,1,3) = 1;  FP(3,1,1) = -1;  FP(1,3,1) = -1;  FP(1,1,1) = -1;    % z^2 = x^2 + y^2 + 1
    a = struct('FP', FP, 'FB', FB, 'num', num, 'den', 1, 'side', side);
end

function a = randRow()
% A random row on a random quadratic fibre, as aCERTIFIEDRowFailsAtEveryPointRandomInputs draws them.
    FP = zeros(3,3,3);
    FP(1,1,3) = randi([1 3]) * sign(randn);
    FP(:,:,2) = randi([-2 2], 3) .* fliplr(triu(ones(3)));
    FP(:,:,1) = randi([-3 3], 3) .* fliplr(triu(ones(3)));
    br = [-2 -1 1 2];
    a = struct('FP', FP, 'FB', br(randi(4)), 'num', [0 0 0 0, randi([-3 3], 1, 6)], ...
               'den', randi(3), 'side', 2*randi(2) - 3);
end

function v = polyval2(P, X)
% A dense polyN at the rows of X.
    v = zeros(size(X,1), 1);
    for i = 1:size(P,1)
        for j = 1:size(P,2)
            v = v + P(i,j) * X(:,1).^(i-1) .* X(:,2).^(j-1);
        end
    end
end
