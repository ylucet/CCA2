classdef algRowFailsSurdTest < matlab.unittest.TestCase
% [slow] algRowFails' SURD FACE ROW (2026-11-03): z_k > z_w certified false when both faces are
% polynomials over one quadratic field Q(sqrt k) -- discriminant k times a square -- by reading
% z_w - z_k = U + sqrt(k) V directly, with the even-multiplicity factors of its norm U^2 - k V^2
% removed. Erratum Example 2's pair cells are made of such faces.
%
% BUCKET: slow. The norm's square-free part is an exact factorisation (Symbolic Math Toolbox): Yun's
% gcds in doubles overflow on 192 of 200 random degree-8 norms (measured 2026-11-03), so the
% reduction that decides a touching curve is a symbolic call. The rest of algRowFails is in
% algRowFailsTest (fast).

    methods (Test)
        function aSURDFaceRowTOUCHINGAtACornerOnTheEdge(t)
        % Both faces polynomials over Q(sqrt 2) on x >= 0 -- z_w = sqrt2 x, z_k = sqrt2 x - x^2 - y^2
        % (discriminants 8 x^2, each 2 times a square) -- touching at the corner (0, 0) of the
        % region, where the discriminants vanish along its edge x = 0 -- the SHAPE of erratum Example
        % 2's pair cells (2026-11-03), small enough that the general face row decides it too (the
        % next test is the one only the surd row decides). Plus faces touching along a curve.
            [Fw, Fk] = deal(zeros(3,3,3));
            Fw(1,1,3) = 1;  Fw(3,1,1) = -2;                                    % z^2 - 2 x^2
            Fk(1,1,3) = 1;  Fk(3,1,2) = 2;  Fk(1,3,2) = 2;                     % (z + x^2 + y^2)^2 - 2 x^2
            Fk = cat(1, Fk, zeros(2,3,3));  Fk = cat(2, Fk, zeros(5,2,3));
            Fk(5,1,1) = 1;  Fk(3,3,1) = 2;  Fk(1,5,1) = 1;  Fk(3,1,1) = -2;
            fr = @(Fa, ba, Fb, bb) struct('FP', Fa, 'FB', ba, 'num', [], 'den', [], 'side', [], ...
                                          'W', struct('FP', Fb, 'FB', bb));
            X = [0; 1];                                                         % x >= 0
            inBox = @(b) b(2) >= 0;  box = [0 1 -1 1];
            t.verifyTrue(algRowFails(fr(Fk, -1, Fw, -1), box, inBox, 6, {X}, [0 0]), 'z_k > z_w fails');
            t.verifyFalse(algRowFails(fr(Fw, -1, Fk, -1), box, inBox, 6, {X}, [0 0]), 'z_w > z_k holds off (0, 0)');
            % touching along the CURVE y = 0: z_k = sqrt2 x - y^2, the norm carries the square
            G = zeros(3,3,3);  G(1,1,3) = 1;  G(1,3,2) = 2;  G = cat(2, G, zeros(3,2,3));
            G(1,5,1) = 1;  G(3,1,1) = -2;                                       % (z + y^2)^2 - 2 x^2
            t.verifyTrue(algRowFails(fr(G, -1, Fw, -1), box, inBox, 6, {X}, [0 0]), 'touching along y = 0');
        end

        function aSURDFaceRowOnERRATUMExampleTwosCorner(t)
        % Erratum Example 2, pair 1: two faces of one cell (discriminants 55 times a square) that
        % touch at the cell's corner (47.5, 2.5), the cell the wedge 2y >= 5, 5y - 3x + 130 >= 0.
        % Before the surd row nothing decided it: the blow-up was refused (the corner a few ulps off
        % 47.5 in doubles, and partial sums past 2^52), and done exactly the general face row's E2
        % has the lowest form -c (3u - 5v)^2 there, degenerate along the edge. (A 60-digit check:
        % z_k - z_w ~ -c(theta) r^2 on every ray; sampled in doubles it reads up to +2e-9, which is
        % algBranch's rounding at coefficients ~1e9.)
            Fk = zeros(3,3,3);
            Fk(1,1,1) = 1700725625;  Fk(2,1,1) = -77532000;  Fk(3,1,1) = 883600;  Fk(1,2,1) = -31149400;
            Fk(2,2,1) = 718240;  Fk(1,3,1) = -549216;  Fk(1,1,2) = 728400;  Fk(2,1,2) = -16640;
            Fk(1,2,2) = -448;  Fk(1,1,3) = 64;
            Fw = zeros(5,5,3);
            Fw(1,1,1) = -28730000;  Fw(2,1,1) = 5720000;  Fw(3,1,1) = -167400;  Fw(5,1,1) = 27;
            Fw(1,2,1) = 3468400;  Fw(2,2,1) = -586000;  Fw(3,2,1) = 19476;  Fw(4,2,1) = -180;
            Fw(1,3,1) = -52628;  Fw(2,3,1) = -7720;  Fw(3,3,1) = 120;  Fw(1,4,1) = -2880;
            Fw(2,4,1) = 600;  Fw(1,5,1) = 300;  Fw(1,1,2) = 67600;  Fw(3,1,2) = -36;  Fw(1,2,2) = -4184;
            Fw(2,2,2) = 120;  Fw(1,3,2) = -320;  Fw(1,1,3) = 12;
            fr = @(Fa, ba, Fb, bb) struct('FP', Fa, 'FB', ba, 'num', [], 'den', [], 'side', [], ...
                                          'W', struct('FP', Fb, 'FB', bb));
            fac = {[-5 2], [130 5; -3 0]};                                     % 2y - 5, 5y - 3x + 130
            inBox = @(b) all(cellfun(@(g) polySignOn(g, b, @(~) true, 0) ~= -1, fac));
            box = [47.4 47.6 2.5 2.7];  p = [47.5 - 7e-15, 2.5];                % as the vertex is stored
            t.verifyTrue(algRowFails(fr(Fk, -1, Fw, -1), box, inBox, 6, fac, p), 'z_k > z_w fails on the cell');
            t.verifyFalse(algRowFails(fr(Fw, -1, Fk, -1), box, inBox, 6, fac, p), 'z_w > z_k holds off the corner');
        end

        function aCERTIFIEDSurdFaceRowNeverHoldsRandomInputs(t)
        % The oracle: random faces z = P0 + e c sqrt(k) S, (z - P0)^2 - c^2 k S^2 = 0 (P0 quadratic, S
        % linear, k in {2, 3}), either root, sometimes a rational face instead; wherever the face row
        % is certified false, no sample of the box has both branches with z_k > z_w.
            rng(43);  nCert = 0;
            for trial = 1:1000
                kap = 2 + (rand < 0.5);
                F = cell(1, 2);  B = zeros(1, 2);
                for f = 1:2
                    P0 = zeros(3,3);  P0(1,1) = randi([-3 3]);  P0(2,1) = randi([-2 2]);  P0(1,2) = randi([-2 2]);
                    P0(3,1) = randi([-1 1]);  P0(2,2) = randi([-1 1]);  P0(1,3) = randi([-1 1]);
                    if f == 2 && rand < 0.25                                    % rational: z = P0
                        F{f} = cat(3, -P0, ones(3,3) .* [1 0 0; 0 0 0; 0 0 0]);  B(f) = 1;  continue
                    end
                    S = zeros(2,2);  S(1,1) = randi([-2 2]);  S(2,1) = randi([-2 2]);  S(1,2) = randi([-2 2]);
                    c = randi([1 2]);
                    C = conv2(P0, P0) - c^2 * kap * padTo(conv2(S, S), 5);
                    Fq = zeros(5,5,3);  Fq(:,:,1) = C;  Fq(1:3,1:3,2) = -2 * P0;  Fq(1,1,3) = 1;
                    F{f} = Fq;  B(f) = 2 * randi(2) - 3;
                end
                a = struct('FP', F{1}, 'FB', B(1), 'num', [], 'den', [], 'side', [], ...
                           'W', struct('FP', F{2}, 'FB', B(2)));
                box = [-1 1 -1 1] * (0.25 + rand);
                if ~algRowFails(a, box, @(~) true, 5), continue, end
                nCert = nCert + 1;
                X = [box(1) + (box(2) - box(1)) * rand(3000,1), box(3) + (box(4) - box(3)) * rand(3000,1)];
                zk = algBranch(F{1}, X, B(1));  zw = algBranch(F{2}, X, B(2));
                both = isfinite(zk) & isfinite(zw);
                t.verifyFalse(any(zk(both) > zw(both) + 1e-9 * max(1, abs(zw(both)))), ...
                              sprintf('trial %d: z_k > z_w at a sample', trial));
            end
            t.verifyGreaterThan(nCert, 40, 'too few certified trials to mean anything');
        end

    end
end

function Q = padTo(P, n)
    Q = zeros(n, n);  Q(1:size(P,1), 1:size(P,2)) = P;
end
