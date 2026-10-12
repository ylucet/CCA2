classdef erratumExamplesTest < matlab.unittest.TestCase
% The three examples of the Karmarkar-Lucet ERRATUM (COAP 94 (2026) 747-780 and JOGO 94 (2026)
% 3-34; source D:\AI\67191168115b2b8bc24afaab\erratum.tex, Lean companion in its lean/ folder),
% run through the EXACT pipeline: conjQ for f*, biconjQ for f**.
%
% BUCKET: fast (conjQ is integer arithmetic; no symbolic call). The BICONJUGATE tests of the same
% examples are in erratumBiconjTest (normal bucket: biconjQ takes ~45 s on Example 1); the
% fixtures and oracles below are shared with it.
%
% WHY THESE THREE. Each is the smallest witness of one way the original papers' structure theorems
% fail, so together they cover what the corrected theorems add:
%   Example 1  a strictly CONVEX piece next to an indefinite one -- assumption (H) fails. f* has a
%              CIRCULAR cell boundary; f** has a circular cell boundary and a norm-shaped cell
%              (affine + const*|x - v|), so no polyhedral/rational description exists.
%   Example 2  three INDEFINITE pieces with different quadratics -- (H) holds. f* has a HYPERBOLIC
%              boundary (two convex edges in different pieces, delta ~= 0); f** is polyhedral but
%              on one cell is the minimum of a 1-D function whose stationarity is a QUARTIC, and is
%              not a rational function there (erratum Theorem 2 case (b)).
%   Example 3  an UNBOUNDED domain: dom f* is a half-plane and f** = x1^2 is not attained by any
%              finite convex combination.
%
% THE ORACLES, and why they are independent of CCA2:
%   f*   by exact candidate enumeration over each piece -- vertices, the stationary point of each
%        edge of positive curvature, and the stationary point of a strictly convex piece. That is
%        the complete Frank-Wolfe candidate set for a quadratic on a triangle, so the max over it
%        IS the sup; it shares no code with conjQ. (The erratum's own tables agree with it to
%        1e-11 on 3000 random s per example -- scratch check, 2026-10-10.)
%   f**  by the erratum's closed forms per cell. Each was BRACKETED independently before being
%        written here (2026-10-10, scratch check): between a dual lower bound sup_s <s,x> - f*(s)
%        and a primal upper bound min over chords through x, at one interior point of every
%        cell. On top of the values, two DEFINITIONS that need no formula: f** <= f on the domain,
%        and f** >= <s,x> - f*(s) for every s (Fenchel-Young), with f* from the oracle above.

    properties (Constant)
        % Example 2's cell 3 (case (b)): x = W + p h1 + r h2 with these constants, erratum Sec. 5.
        EX2_DELTA = 60375/8
    end

    methods (Test)

        % ============================== Example 1 ===============================================

        function ex1ConjugateMatchesTheDefinition(testCase)
            [f, pieces] = erratumExamplesTest.example1();
            g = conjQ(f);
            S = erratumExamplesTest.ex1DualSample();
            want = erratumExamplesTest.fstarOracle(S, pieces);
            testCase.verifyEqual(g.eval(S), want, 'RelTol', 1e-10, 'AbsTol', 1e-10, ...
                'conjQ disagrees with the exact sup on Example 1');
        end

        function ex1ConjugateIsTheMaxOfTheErratumFunctionsAcrossTheCircle(testCase)
        % The boundary between the convex piece's cell, (s1^2+s2^2)/4, and the vertex v's cell,
        % s1 - 3 s2 + 11, is the circle (s1-2)^2 + (s2+6)^2 = 84. Points straddling the arc, just
        % inside and just outside, are where a polygonal subdivision would have to be wrong.
            [f, ~] = erratumExamplesTest.example1();
            g = conjQ(f);
            th = linspace(1.692, 1.966, 25)';        % the arc between the cell's two corners
            onArc = [2 + sqrt(84)*cos(th), -6 + sqrt(84)*sin(th)];
            for rho = [-1e-3, 0, 1e-3]
                S = [2 + (sqrt(84)+rho)*cos(th), -6 + (sqrt(84)+rho)*sin(th)];
                gI = (S(:,1).^2 + S(:,2).^2)/4;
                gv = S(:,1) - 3*S(:,2) + 11;
                testCase.verifyEqual(g.eval(S), max(gI, gv), 'RelTol', 1e-10, 'AbsTol', 1e-10, ...
                    sprintf('wrong value at radius sqrt(84) %+g', rho));
            end
            testCase.verifyEqual((onArc(:,1).^2 + onArc(:,2).^2)/4, ...
                onArc(:,1) - 3*onArc(:,2) + 11, 'AbsTol', 1e-9, 'fixture self-check: the arc ties');
        end

        % ============================== Example 2 ===============================================

        function ex2ConjugateMatchesTheDefinition(testCase)
            [f, pieces] = erratumExamplesTest.example2();
            g = conjQ(f);
            S = erratumExamplesTest.ex2DualSample();
            want = erratumExamplesTest.fstarOracle(S, pieces);
            testCase.verifyEqual(g.eval(S), want, 'RelTol', 1e-10, 'AbsTol', 1e-8, ...
                'conjQ disagrees with the exact sup on Example 2');
        end

        function ex2ConjugateIsRightAcrossTheHyperbola(testCase)
        % BC | DE: g_BC = g_DE on (s1-335/2)^2/6 - (2497-5s1-2s2)^2/428 = -60375/8, a hyperbola,
        % between s ~ (-5, 104.2) and s ~ (64.9, 88.6). Solve for s2 on that branch and probe
        % just above and below it.
            [f, pieces] = erratumExamplesTest.example2();
            g = conjQ(f);
            s1 = linspace(-4, 64, 30)';
            % (2497 - 5 s1 - 2 s2)^2 = 428 ((s1-335/2)^2/6 + 60375/8); the cell sits on the root
            % with 2497 - 5 s1 - 2 s2 > 0.
            rhs = 428*((s1 - 335/2).^2/6 + erratumExamplesTest.EX2_DELTA);
            s2 = (2497 - 5*s1 - sqrt(rhs))/2;
            gBC = @(S) (S(:,1)+50).^2/6 + 10*S(:,2) - 260;
            gDE = @(S) 35*S(:,1) - 5*S(:,2) + 1285/2 + (892 - 5*S(:,1) - 2*S(:,2)).^2/428;
            testCase.verifyEqual(gBC([s1 s2]), gDE([s1 s2]), 'RelTol', 1e-10, ...
                'fixture self-check: the two edge functions tie on the hyperbola');
            for ds = [-0.05, 0, 0.05]
                S = [s1, s2 + ds];
                testCase.verifyEqual(g.eval(S), erratumExamplesTest.fstarOracle(S, pieces), ...
                    'RelTol', 1e-10, 'AbsTol', 1e-8, sprintf('wrong value at s2 offset %+g', ds));
            end
        end

        function ex2ConjugateValueDoesNotDependOnTheOtherPointsInTheCall(testCase)
        % Found 2026-10-10 by ex2ConjugateMatchesTheDefinition: QuaCon.eval scaled its membership
        % tolerance by the largest |x| in the WHOLE batch and took the first face accepted within
        % it, so a point near a cell boundary got a neighbour's (smaller) value as soon as a far
        % point shared the call. A point's value is a function of that point only.
            [f, pieces] = erratumExamplesTest.example2();
            g = conjQ(f);
            P = [45.25 99; 60.75 91.25; 66.5 85.5; 74.75 96.5; 64.94 88.60];
            want = erratumExamplesTest.fstarOracle(P, pieces);
            alone = arrayfun(@(i) g.eval(P(i,:)), (1:size(P,1))');
            batch = g.eval([P; 1000 1000]);
            testCase.verifyEqual(alone, want, 'RelTol', 1e-10, 'AbsTol', 1e-8, 'one point per call');
            testCase.verifyEqual(batch(1:end-1), want, 'RelTol', 1e-10, 'AbsTol', 1e-8, ...
                'the same points with a far point in the same call');
        end

        % ============================== Example 3 ===============================================

        function ex3ConjugateIsAHalfPlaneParabola(testCase)
        % f = x1^2 + |x1| x2 on {x2 >= 0}: f*(s) = s1^2/4 on {s2 <= 0}, +inf on {s2 > 0}.
            f = erratumExamplesTest.example3();
            g = conjQ(f);
            rng(20261010);
            S = [randn(200,1)*5, -abs(randn(200,1))*5; 3 0; -2 0; 0 0];
            testCase.verifyEqual(g.eval(S), S(:,1).^2/4, 'RelTol', 1e-12, 'AbsTol', 1e-12);
            T = [randn(50,1)*5, abs(randn(50,1))*5 + 1e-3];
            testCase.verifyTrue(all(isinf(g.eval(T)) & g.eval(T) > 0), ...
                'f* must be +inf wherever s2 > 0');
        end

    end

    % ============================== fixtures and oracles =========================================
    methods (Static)

        function [f, pieces] = example1()
        % P1 = conv{A,B,C} with x^2+y^2, P2 = conv{A,B,v} with x^2+xy-y^2.
            V = [-2 0; 2 0; 0 2; 1 -3];               % A B C v
            E = [1 4 1; 4 2 1; 2 3 1; 3 1 1; 1 2 1];
            F = [2 0; 2 0; 1 0; 1 0; 1 2];
            H1 = [2 0; 0 2];  H2 = [2 1; 1 -2];
            f = QuaPol(V, E, [0 0 0 0 2 0 2 0 0 0; 0 0 0 0 2 1 -2 0 0 0], F);
            pieces = {struct('T', V([1 2 3],:), 'H', H1, 'L', [0;0], 'k', 0), ...
                      struct('T', V([1 2 4],:), 'H', H2, 'L', [0;0], 'k', 0)};
        end

        function [f, pieces] = example2()
        % Three indefinite pieces on conv{A,B,C}, conv{A,C,D}, conv{A,D,E}.
            V = [0 0; 15 10; 60 10; 35 -5; 10 -15];  % A B C D E
            E = [1 5 1; 5 4 1; 4 3 1; 3 2 1; 2 1 1; 1 3 1; 1 4 1];
            F = [3 0; 3 0; 2 0; 1 0; 1 0; 1 2; 2 3];
            H = {[3 -5; -5 5], [0 4; 4 5], [2 10; 10 -9]};
            rows = zeros(3,10);
            for k = 1:3, rows(k,:) = [0 0 0 0, H{k}(1,1), H{k}(1,2), H{k}(2,2), 0, 1, 0]; end
            f = QuaPol(V, E, rows, F);
            tri = {[1 2 3], [1 3 4], [1 4 5]};
            pieces = cell(1,3);
            for k = 1:3
                pieces{k} = struct('T', V(tri{k},:), 'H', H{k}, 'L', [0;1], 'k', 0);
            end
        end

        function f = example3()
        % x1^2 + x1 x2 on the first quadrant, x1^2 - x1 x2 on the second; three rays from 0.
            V = [0 0; 1 0; 0 1; -1 0];
            E = [1 2 0; 1 3 0; 1 4 0];
            F = [1 0; 2 1; 0 2];
            f = QuaPol(V, E, [0 0 0 0 2 1 0 0 0 0; 0 0 0 0 2 -1 0 0 0 0], F);
        end

        function v = fstarOracle(S, pieces)
        % max over pieces of the exact sup of <s,x> - q(x) on a triangle, by the complete
        % candidate set (vertices, stationary points of edges of positive curvature, interior
        % stationary point of a strictly convex q).
            v = -inf(size(S,1),1);
            for k = 1:numel(pieces)
                P = pieces{k};  T = P.T;  H = P.H;  L = P.L;  k0 = P.k;
                q = @(X) 0.5*sum((X*H).*X, 2) + X*L + k0;
                for i = 1:3
                    v = max(v, S*T(i,:)' - q(T(i,:)));
                    u = T(i,:);  d = T(mod(i,3)+1,:) - u;
                    al = 0.5*d*H*d';
                    if al > 0
                        beta = (u*H + L')*d';
                        t = (S*d' - beta)/(2*al);
                        ok = t > 0 & t < 1;
                        Y = u + t*d;
                        v(ok) = max(v(ok), sum(S(ok,:).*Y(ok,:), 2) - q(Y(ok,:)));
                    end
                end
                if H(1,1) > 0 && det(H) > 0
                    Y = (H \ (S' - L))';
                    ok = erratumExamplesTest.inTriangle(Y, T);
                    v(ok) = max(v(ok), sum(S(ok,:).*Y(ok,:), 2) - q(Y(ok,:)));
                end
            end
        end

        function tf = inTriangle(X, T)
            c = zeros(size(X,1), 3);
            for i = 1:3
                a = T(i,:);  b = T(mod(i,3)+1,:);
                c(:,i) = (b(1)-a(1))*(X(:,2)-a(2)) - (b(2)-a(2))*(X(:,1)-a(1));
            end
            tf = all(c >= -1e-12, 2) | all(c <= 1e-12, 2);
        end

        function y = fOnPieces(X, pieces)
            y = inf(size(X,1),1);
            for k = numel(pieces):-1:1
                P = pieces{k};
                in = erratumExamplesTest.inTriangle(X, P.T);
                y(in) = 0.5*sum((X(in,:)*P.H).*X(in,:), 2) + X(in,:)*P.L + P.k;
            end
        end

        function S = ex1DualSample()
            rng(20261012);
            S = [(rand(400,2) - 0.5)*30; 0 0; 2 -6; 0 3; -1.55 2.45; 0.9 3.1];
        end

        function S = ex2DualSample()
            rng(20261013);
            S = [(rand(400,1)*900 - 300), (rand(400,1)*1200 - 800); 64.94 88.60; -5 104.21];
        end

        function [X, want] = ex1BiconjCells()
        % One or two interior points per cell of erratum Table 3 / Figure 1(b).
            v = [1 -3];
            X = [0 1.8; -0.1 1.85;               % cell 1: P1 outside the disc |x-v|^2 < 21
                 0 1.0;  0.2 0.5;                % cell 2: conv({v} u Gamma)
                 -1.2 0;  -0.8 -0.6;             % cell 3: conv{v, A, P-}
                 1.4 0;   1.2 -0.8];             % cell 4: conv{v, P+, B}
            x = X(:,1);  y = X(:,2);
            c1 = x.^2 + y.^2;
            c2 = -11 + 2*((x - v(1))*v(1) + (y - v(2))*v(2)) + 2*sqrt(21)*hypot(x - v(1), y - v(2));
            c3 = -(15*x.^2 + 2*x.*y - 2*x + 7*y.^2 + 18*y - 16)./(2*(x - y - 4));
            c4 = (45*x.^2 + 26*x.*y - 56*x + 13*y.^2 + 8*y - 4)./(4*(x + y + 2));
            want = [c1(1:2); c2(3:4); c3(5:6); c4(7:8)];
        end

        function [X, want, names] = ex2BiconjCells()
            % cell 6's point was (32,-5.5) until 2026-10-11, which is NOT in cell 6: it lies in cell 4,
            % conv{G1,G2,G3}, barycentric (0.004, 0.311, 0.685). (35.5,-4) is inside cell 6 =
            % conv{G2,D,G3}, barycentric (0.385, 0.231, 0.384).
            X = [8 0.5; 13 -8; 25 0; 28 3; 50 7; 35.5 -4];
            names = {'cell 1 (affine)', 'cell 2 (form (a), vertex B + edge DE)', ...
                     'cell 3 (case (b), quartic)', 'cell 3 (case (b), quartic)', ...
                     'cell 5 (delta = 0 at W = C)', 'cell 6 (delta = 0 at W = D)'};
            want = zeros(6,1);
            want(1) = (-2125*X(1,1) + 2791*X(1,2))/26;
            mu = (2*X(2,1) - 5*X(2,2) + 20)/115;  nu = (-3*X(2,1) - 4*X(2,2) + 85)/115;
            want(2) = 2675*nu^2/mu - 4460*nu - 1285/2*mu - 305/2*(1 - mu);
            want(3) = erratumExamplesTest.ex2Cell3(X(3,:));
            want(4) = erratumExamplesTest.ex2Cell3(X(4,:));
            x = X(5,:);
            p = x(1)/45 - x(2)/27 - 26/27;  r = 2/3 - x(2)/15;
            want(5) = (sqrt(6075/2)*abs(p) + sqrt(4125/2)*abs(r))^2 + 5850*p - 5365*r + 2660;
            x = X(6,:);
            p = (2*x(1) - 5*x(2) - 95)/25;  r = (-3*x(1) + 5*x(2) + 130)/25;
            want(6) = (sqrt(4125/2)*abs(p) + sqrt(2675)*abs(r))^2 - 1240*p - 4460*r - 1285/2;
        end

        function v = ex2Cell3(x)
        % min over lambda in (0,1) of F_x, erratum (eq:min); its stationary point is the unique
        % root in (0,1) of the quartic Q_x, which is what the minimum is taken at.
            p = x(1)/45 - x(2)/18 - 19/18;  r = 1 - x(2)/10;  dl = erratumExamplesTest.EX2_DELTA;
            Fx = @(l) 6075*p^2./(2*l) + 2675*r^2./(1 - l) - dl*l + 15075/2*p - 12485*r + 48265/4;
            Qx = [-dl, 2*dl, 2675*r^2 - 6075/2*p^2 - dl, 6075*p^2, -6075/2*p^2];
            lam = roots(Qx);
            lam = real(lam(abs(imag(lam)) < 1e-12 & real(lam) > 0 & real(lam) < 1));
            assert(isscalar(lam), 'fixture: Q_x must have exactly one root in (0,1)');
            v = Fx(lam);
        end

        function checkEnvelopeDefinition(testCase, h, f, pieces, S)
        % f** <= f on the domain, and f** >= <s,x> - f*(s) for every sampled s.
            rng(20261014);
            lo = min(f.V, [], 1);  hi = max(f.V, [], 1);
            X = lo + rand(3000,2).*(hi - lo);
            fx = erratumExamplesTest.fOnPieces(X, pieces);
            X = X(isfinite(fx),:);  fx = fx(isfinite(fx));
            hx = h.eval(X);
            testCase.verifyLessThanOrEqual(max(hx - fx), 1e-9 * max(1, max(abs(fx))), ...
                'f** exceeds f somewhere on the domain');
            fs = erratumExamplesTest.fstarOracle(S, pieces);
            lowerBound = max(X*S' - fs', [], 2);
            testCase.verifyGreaterThanOrEqual(min(hx - lowerBound), -1e-9 * max(1, max(abs(fx))), ...
                'f** lies below an affine minorant <s,x> - f*(s): it is not the LARGEST convex minorant');
        end
    end
end
