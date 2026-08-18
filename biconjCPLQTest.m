classdef biconjCPLQTest < matlab.unittest.TestCase
    % Tests for biconjCPLQ (the 'cplq' biconjugate f** = (f*)*).
    %
    % The case that matters here is a SINGLE BOUNDED TRIANGLE: f* is then finite everywhere, i.e.
    % it lives on an unbounded multi-face parabolic domain that conjCPLQ still rejects, so the
    % literal double conjugation cannot close the loop. biconjCPLQ instead returns Step 1's convex
    % envelope, which for a compact domain IS f** -- see biconjCPLQ.m for the derivation.
    %
    % Ground truth below uses NONE of the conjugate pipeline: f*(s) is obtained by maximizing
    % <s,x> - q(x) over the triangle exactly (the max of a quadratic over a polytope is attained at
    % a vertex, at a 1D critical point of an edge restriction, or at an interior critical point),
    % and f**(x) = sup_s <s,x> - f*(s) by a coarse grid plus Nelder-Mead. So a passing test says
    % the envelope really is the biconjugate, not merely that two CCA2 routines agree.

    properties (Constant)
        E3 = [1 2 1; 2 3 1; 3 1 1];
        F3 = [1 0; 1 0; 1 0];
    end

    methods (Test)
        function convexOverABoxIsItsOwnBiconjugate(testCase)
        % biconj IS the closed convex envelope, so a CONVEX f is its own answer and there is
        % nothing to compute. Pinned because the short-circuit that says so is easy to lose and
        % its absence is not visible in any VALUE -- only in the clock: this input cost 436 s
        % through conj-of-conj before, and 0.05 s after, for the same function.
            V = [0 0; 1 0; 1 1; 0 1];
            E = [1 2 1; 2 3 1; 3 4 1; 4 1 1];
            F = [1 0; 1 0; 1 0; 1 0];
            p = QuaPol(V, E, [1 0 1 0 0 0], F);          % (x^2 + y^2)/2, convex
            h = p.biconj('cplq');
            testCase.verifyEqual(h.f(1,5:10), p.f(1,5:10), 'AbsTol', 1e-12, ...
                'co f must be f itself for a convex f');
        end

        function separableOverABoxTakesTheOneDimensionalRoute(testCase)
        % f = x^2 - y^2 over the unit box is SEPARABLE, so its envelope is one 1-D envelope per
        % axis: x^2 is convex in x and survives; -y^2 is concave in y, so it is replaced by the
        % CHORD through its endpoint values, which on [0,1] is -y. Hence co f = x^2 - y.
        %
        % This is the case that used to raise MATLAB:badsubscript, and then cost 29 s once the
        % first conjugation learned the separable route while the second did not. It must not
        % go back to either.
            V = [0 0; 1 0; 1 1; 0 1];
            E = [1 2 1; 2 3 1; 3 4 1; 4 1 1];
            F = [1 0; 1 0; 1 0; 1 0];
            p = QuaPol(V, E, [2 0 -2 0 0 0], F);          % x^2 - y^2, indefinite and DIAGONAL
            h = p.biconj('cplq');
            testCase.verifyEqual(h.f(1,5:10), [2 0 0 0 -1 0], 'AbsTol', 1e-12, ...
                'co(x^2 - y^2) over the unit box is x^2 - y');

            % and it really is the envelope: <= f on the box, and equal at the y-endpoints where
            % the chord touches.
            [uu, vv] = meshgrid(linspace(0,1,40), linspace(0,1,40));
            env = uu.^2 - vv;
            testCase.verifyLessThanOrEqual(max(env(:) - (uu(:).^2 - vv(:).^2)), 1e-12, ...
                'the envelope must be a MINORANT of f on the box');
        end

        function aWrongConvexityFlagIsRefusedNotTrusted(testCase)
        % The flag is trusted, but not blindly: the free NECESSARY condition -- every piece's
        % Hessian positive semidefinite -- is still enforced, because honouring a wrong flag
        % would return a NON-convex f as its own convex envelope, silently.
            V = [0 0; 1 0; 1 1; 0 1];
            E = [1 2 1; 2 3 1; 3 4 1; 4 1 1];
            F = [1 0; 1 0; 1 0; 1 0];
            p = QuaPol(V, E, [0 1 0 0 0 0], F);           % x*y: indefinite, NOT convex
            p.fIsConvex = true;                           % ... and the caller says otherwise
            testCase.verifyError(@() p.biconj('cplq'), 'PLQ:biconj:notConvexDespiteFlag');
        end

        function bilinearOverABoxReturnsMcCormickDirectly(testCase)
        % co(b*x*y + d*x + e*y + k) over a rectangle is Al-Khayyal-Falk's McCormick envelope, the
        % max of two affine functions. It is returned in closed form rather than derived through
        % triangulate -> conjugate -> Step 3 -> conjugate, which cost 40-63 s for an answer that
        % has been in the literature since 1983.
        %
        % Checked on BOTH signs of b, because they use different affine pairs and meet on
        % different diagonals of the box -- the anti-diagonal for b > 0, the main diagonal for
        % b < 0.
            E = [1 2 1; 2 3 1; 3 4 1; 4 1 1];
            F = [1 0; 1 0; 1 0; 1 0];
            cases = { [0 1 0 0 0 0],   [0 1 0 1]        % x*y on the unit box
                      [0 1 0 0 0 0],   [-2 3 -1 4]      % x*y, general bounds
                      [0 3 0 7 -2 5],  [0 1 0 1]        % with a linear part
                      [0 -2 0 0 0 0],  [0 1 0 1] };     % b < 0: the other diagonal
            for c = 1:size(cases,1)
                coef = cases{c,1}; bb = cases{c,2};
                xl = bb(1); xu = bb(2); yl = bb(3); yu = bb(4);
                p = QuaPol([xl yl; xu yl; xu yu; xl yu], E, coef, F);
                h = p.biconj('cplq');
                testCase.verifyEqual(h.nf, 2, 'the McCormick envelope has exactly two faces');

                b = coef(2); d = coef(4); e = coef(5); k0 = coef(6);
                if b > 0
                    mc = @(x,y) max(b*(xl*y + yl*x - xl*yl), b*(xu*y + yu*x - xu*yu)) ...
                                + d*x + e*y + k0;
                else
                    mc = @(x,y) max(b*(xu*y + yl*x - xu*yl), b*(xl*y + yu*x - xl*yu)) ...
                                + d*x + e*y + k0;
                end
                [uu, vv] = meshgrid(linspace(xl,xu,12), linspace(yl,yu,12));
                for i = 1:numel(uu)
                    pt = [uu(i), vv(i)];
                    testCase.verifyEqual(h.eval(pt), mc(pt(1), pt(2)), 'AbsTol', 1e-9, ...
                        sprintf('McCormick at (%g,%g), case %d', pt(1), pt(2), c));
                    % and a MINORANT of f, which is what makes it an envelope and not merely a fit
                    testCase.verifyLessThanOrEqual(mc(pt(1),pt(2)) - ...
                        (b*pt(1)*pt(2) + d*pt(1) + e*pt(2) + k0), 1e-9);
                end
            end
        end

        function bilinearOverADiamondRotatesToTheSeparableCase(testCase)
        % x*y separates in the 45-degree frame -- with u = (x+y)/sqrt(2), v = (x-y)/sqrt(2),
        % x*y = (u^2 - v^2)/2 -- but that alone buys nothing, because the ENVELOPE separates only
        % when the DOMAIN is a product in the same coordinates, and rotating the function rotates
        % the domain too. It pays in exactly one shape: a DIAMOND, which rotates to a box.
        %
        % On the unit diamond |x| + |y| <= 1 the answer is a SINGLE smooth quadratic, because the
        % concave direction contributes its chord, which is affine:
        %       co(x*y) = (x+y)^2/4 - 1/4.
        % Check it at the origin: (1/2,-1/2) and (-1/2,1/2) both give x*y = -1/4 and average to
        % the origin, so co f (0,0) = -1/4, which is what the formula gives.
            E = [1 2 1; 2 3 1; 3 4 1; 4 1 1];
            F = [1 0; 1 0; 1 0; 1 0];
            V = [1 0; 0 1; -1 0; 0 -1];              % the unit diamond, counter-clockwise
            p = QuaPol(V, E, [0 1 0 0 0 0], F);      % f = x*y
            h = p.biconj('cplq');
            testCase.verifyEqual(h.nf, 1, 'the envelope over a diamond is one smooth quadratic');

            X = [0 0; 0.5 -0.5; -0.5 0.5; 0.25 0.25; 0.5 0.25; -0.3 -0.2; 0.9 0; 0 -0.75];
            for i = 1:size(X,1)
                want = (X(i,1) + X(i,2))^2/4 - 0.25;
                testCase.verifyEqual(h.eval(X(i,:)), want, 'AbsTol', 1e-12, ...
                    sprintf('co(x*y) over the unit diamond at (%g,%g)', X(i,1), X(i,2)));
                % and a MINORANT of f on the diamond
                testCase.verifyLessThanOrEqual(want - X(i,1)*X(i,2), 1e-12, ...
                    sprintf('envelope must underestimate x*y at (%g,%g)', X(i,1), X(i,2)));
            end
        end

        function singleBoundedTriangleNoLongerErrors(testCase)
            % The exact input SUPPORT_MATRIX.md section 8 listed as blocker 1: f = xy over the
            % unit triangle. conj works and gives a QuaPar; biconj used to raise
            % PLQ:conjCPLQ:notImplemented on the second conjugation.
            p = QuaPol([0 0; 1 0; 0 1], biconjCPLQTest.E3, [0 1 0 0 0 0], biconjCPLQTest.F3);
            testCase.verifyEqual(p.conj().kind(), 'QuaPar');
            b = p.biconj();
            testCase.verifyEqual(b.kind(), 'RatPol');
            % f = xy >= 0 on this triangle and vanishes at all three vertices, so the largest
            % convex minorant is the zero function.
            S = [0.2 0.2; 0.3 0.5; 0.1 0.7; 1/3 1/3];
            testCase.verifyEqual(b.eval(S), zeros(4,1), 'AbsTol', 1e-12);
        end

        function matchesIndependentBiconjugate(testCase)
            % Every Step 1 branch, checked against the pipeline-free ground truth.
            cases = { ...
                {'affine',                  [0 0; 1 0; 0 1], [0 0 0 -1 2 3]}, ...
                {'convex',                  [0 0; 1 0; 0 1], [1 0 1 0 0 0]}, ...
                {'concave',                 [0 0; 1 0; 0 1], [-2 0 -2 0 0 0]}, ...
                {'indefinite 1-convex-edge',[0 0; 2 0; 1 1], [0 1 0 0 0 0]}, ...
                {'indefinite + linear',     [0 0; 3 0; 1 2], [1 3 -1 2 -1 0.5]}, ...
                {'indefinite 2-convex-edge',[2 1; 0 0; 1 0], [0 1 0 0 0 0]}, ...
                {'indefinite 3-convex-edge',[0 0; 1 1; 3 2], [0 1 0 0 0 0]}};
            for c = 1:numel(cases)
                [name, V, f6] = deal(cases{c}{:});
                p = biconjCPLQTest.triangle(V, f6);
                b = p.biconj();
                % TYPE: a RatPar, and WHICH one is deliberately not pinned here.
                % biconj's own contract is "a RatPar -- QuaPol (Case A), RatPol (Case B) or
                % QuaParCPLQ (Case C); call kind() to learn which". Since 2026-08-18 the convex
                % and affine rows come back as the QuaPol they went in as, because co f = f for
                % a convex f and returning the input unchanged is the faithful answer -- Step 1
                % would only re-derive it and widen the type. The VALUE checks below are what
                % this test is really for, and they are unchanged.
                testCase.verifyTrue(isa(b, 'RatPar'), ...
                    sprintf('%s: expected a RatPar biconjugate, got %s', name, class(b)));
                X = [0.6 0.2 0.2; 0.2 0.2 0.6; 1/3 1/3 1/3] * p.V;
                for i = 1:size(X,1)
                    truth = biconjCPLQTest.biconjTrue(f6, p.V, X(i,:));
                    testCase.verifyEqual(b.eval(X(i,:)), truth, 'AbsTol', 1e-9, ...
                        sprintf('%s at (%g,%g)', name, X(i,1), X(i,2)));
                end
            end
        end

        function succeedsWhereTheConjugateItselfFails(testCase)
            % A triangle whose Step 1 envelope splits into 4 faces: Step 2 completes on them,
            % but cPLQ's Step 3 (the cross-piece maximum) does not yet, so f* is not computable
            % (PLQ:conjCPLQ:cplqFailed -- see
            % conjCPLQTest.indefiniteTriangleThreeConvexEdgesUsesStep3 for the full diagnosis).
            % f** is computable anyway, because it IS Step 1's output and never has to cross into
            % the symbolic layer. This is why biconjCPLQ is not a thin wrapper around conj-of-conj.
            %
            % UPDATE (2026-07-28): the 2-convex-edge triangle conv{(2,1),(0,0),(1,0)} used to be
            % in this list too. Its conjugate now works, via Step 2's fallback to cPLQ's symbolic
            % Step 2/3, so it no longer belongs here -- see
            % conjCPLQTest.indefiniteTriangleTwoConvexEdgesSplitViaCPLQStep2.
            %
            % UPDATE (2026-08-13): the conjugate of THIS triangle now works too. Its four-face
            % envelope's per-piece conjugates are curved, and arc-vs-arc assembly took the
            % numeric route all the way, so there is no longer a cplqFailed to expect. The point
            % the test makes still stands and is checked below -- biconj does not go through
            % conj-of-conj -- so what is pinned here now is that conj SUCCEEDS and agrees with
            % the closed-form sup, which is a stronger statement than the old refusal.
            p = biconjCPLQTest.triangle([0 0; 1 1; 3 2], [0 1 0 0 0 0]);   % f = xy
            gc = p.conj();
            testCase.verifyTrue(isa(gc, 'RatPar'));
            for sPt = {[1 1], [-3 -3], [2 -1], [0 0], [5 5]}
                sv = sPt{1};
                testCase.verifyEqual(conjCPLQTest.evalConjResult(gc, sv), ...
                    convEnvCPLQTest.supBilinearOverPoly(sv, [0 0; 1 1; 3 2]), ...
                    'AbsTol', 1e-6, sprintf('conj at (%g,%g)', sv(1), sv(2)));
            end
            b = p.biconj();
            testCase.verifyEqual(b.kind(), 'RatPol');
            testCase.verifyGreaterThan(b.nf, 1);                % genuinely split envelope
        end

        function isAConvexUnderestimatorTouchingAtTheVertices(testCase)
            % The three defining properties of f**, independent of any formula: it underestimates
            % f on dom f, it is convex, and it agrees with f at the triangle's vertices (each
            % vertex is an extreme point, so no strictly-better convex minorant exists there).
            p = biconjCPLQTest.triangle([0 0; 3 0; 1 2], [1 3 -1 2 -1 0.5]);
            b = p.biconj();
            B = [0.6 0.2 0.2; 0.2 0.6 0.2; 0.2 0.2 0.6; 1/3 1/3 1/3; 0.5 0.3 0.2; 0.1 0.45 0.45];
            X = B * p.V;
            testCase.verifyLessThanOrEqual(b.eval(X) - p.eval(X), 1e-9);
            testCase.verifyEqual(b.eval(p.V), p.eval(p.V), 'AbsTol', 1e-10);
            for i = 1:size(X,1)
                for j = i+1:size(X,1)
                    mid = 0.5*(X(i,:) + X(j,:));
                    testCase.verifyLessThanOrEqual( ...
                        b.eval(mid) - 0.5*(b.eval(X(i,:)) + b.eval(X(j,:))), 1e-9);
                end
            end
        end

        function convexTriangleIsItsOwnBiconjugate(testCase)
            % f** = f for a closed convex f -- on a bounded domain too, not just Case A.
            p = biconjCPLQTest.triangle([0 0; 1 0; 0 1], [1 0 1 0 0 0]);   % 1/2(x^2+y^2)
            b = p.biconj();
            X = [0.2 0.2; 0.6 0.2; 0.2 0.6; 1/3 1/3];
            testCase.verifyEqual(b.eval(X), p.eval(X), 'AbsTol', 1e-12);
        end

        function domainOfTheBiconjugateIsTheTriangle(testCase)
            % dom f** = conv(dom f): finite inside, +Inf outside.
            p = biconjCPLQTest.triangle([0 0; 1 0; 0 1], [0 1 0 0 0 0]);
            b = p.biconj();
            testCase.verifyTrue(all(isfinite(b.eval([0.2 0.2; 0.1 0.5]))));
            testCase.verifyEqual(b.eval([2 2; -1 0.5]), [Inf; Inf]);
        end

        function polyhedralQuaParTriangleTakesTheSameRoute(testCase)
            % QuaPar.biconj routes through biconjCPLQ too: an all-zero-Ec triangle is the same
            % mesh as the QuaPol above and must give the same answer.
            V = [0 0; 2 0; 1 1];
            q = QuaPar(V, biconjCPLQTest.E3, zeros(3,6), [0 1 0 0 0 0], biconjCPLQTest.F3);
            b = q.biconj();
            testCase.verifyEqual(b.kind(), 'RatPol');
            X = [0.6 0.2; 1.0 0.6; 1.0 1/3];
            expected = biconjCPLQTest.triangle(V, [0 1 0 0 0 0]).biconj().eval(X);
            testCase.verifyEqual(b.eval(X), expected, 'AbsTol', 1e-12);
        end

        function fullDomainQuadraticIsUnchangedCaseA(testCase)
            % Case A keeps the literal double conjugation, and its QuaPol return type.
            p = QuaPol.energy();
            b = p.biconj('cplq');
            testCase.verifyEqual(b.kind(), 'QuaPol');
            testCase.verifyEqual(b.f, p.f, 'AbsTol', 1e-12);
        end

        function conjugateOfATriangleBiconjugatesBackToItself(testCase)
            % HISTORY: this was `unsupportedShapesStillErrorAsBefore`, and it asserted that
            % g.biconj() ERRORS for g = (x*y on a triangle)* -- first with
            % PLQ:conjCPLQ:notImplemented (the isDomBounded gate), then with
            % QuaParCPLQ:conj:emptyResult. Neither was a property of the input: the biconjugate
            % simply did not work for ANY input, because the last step of the algorithm --
            % "max of all those conjugates" -- was missing from the second pass, and
            % conjugateOfPiecePoly mis-assigned regions for pieces whose domain carried a
            % redundant constraint. Both are fixed, so the input is now supported and the test
            % asserts the VALUE instead of an error.
            %
            % g = (x*y on conv{(0,0),(1,0),(0,1)})* = max(0,s1,s2), which is convex, so
            % g** = g. That also exercises the shape the second pass found hardest: an
            % unbounded, multi-face, piecewise-affine conjugate.
            p = biconjCPLQTest.triangle([0 0; 1 0; 0 1], [0 1 0 0 0 0]);
            g = p.conj();
            testCase.verifyEqual(g.kind(), 'QuaPar');
            h = g.biconj();
            X = [0 0; 2 1; -1 3; -2 -3; 0.5 0.25; 1 1; -0.5 -0.5];
            for i = 1:size(X,1)
                v = evalFunctionNDomain(h.fnd, X(i,:));
                testCase.verifyEqual(v, max([0, X(i,1), X(i,2)]), 'AbsTol', 1e-9, ...
                    sprintf('g** must equal max(0,x,y) at (%g,%g)', X(i,1), X(i,2)));
            end
        end

        function nonCplqEnginesStillError(testCase)
            % Engine dispatch is untouched: only 'cplq' takes the new route.
            p = biconjCPLQTest.triangle([0 0; 1 0; 0 1], [0 1 0 0 0 0]);
            testCase.verifyError(@() p.biconj('pqp'),   'PLQ:conj:engine');
            testCase.verifyError(@() p.biconj('graph'), 'PLQ:conj:engine');
            testCase.verifyError(@() p.biconj('nope'),  'PLQ:conj:engine');
        end
    end

    methods (Static, Access = private)
        function p = triangle(V, f6)
            % QuaPol on one bounded triangle, vertices forced counter-clockwise so that
            % F=[k 0] (face on the left of each directed edge) is consistent.
            if det([V(2,:)-V(1,:); V(3,:)-V(1,:)]) < 0, V = V([1 3 2],:); end
            p = QuaPol(V, biconjCPLQTest.E3, f6, biconjCPLQTest.F3);
        end

        function v = fstarExact(f6, V, s)
            % f*(s) = max_{x in triangle V} <s,x> - q(x), q = 1/2 x'Qx + L'x + kappa in
            % QuaPol.matrixForm's stored convention (Q = [c5 c6; c6 c7]). The maximum of a
            % quadratic over a polytope is attained at a vertex, at an interior critical point of
            % an edge restriction, or at an interior critical point of the whole quadratic, so
            % enumerating those candidates is exact.
            Q = [f6(1) f6(2); f6(2) f6(3)];  L = [f6(4); f6(5)];  kap = f6(6);
            obj = @(X) X*s(:) - (0.5*sum((X*Q).*X,2) + X*L + kap);
            cand = V;
            for e = 1:3
                A = V(e,:); Bv = V(mod(e,3)+1,:); d = Bv - A;
                a2 = -0.5*(d*Q*d');                     % t^2 coefficient of the edge restriction
                a1 = d*s(:) - (A*Q*d' + d*L);           % t   coefficient
                if a2 < -1e-14
                    t = -a1/(2*a2);
                    if t > 0 && t < 1, cand(end+1,:) = A + t*d; end %#ok<AGROW>
                end
            end
            if rcond(Q) > 1e-12
                xi  = (Q \ (s(:) - L))';
                lam = [V ones(3,1)]' \ [xi 1]';         % barycentric coordinates
                if all(lam > -1e-12), cand(end+1,:) = xi; end
            end
            v = max(obj(cand));
        end

        function v = biconjTrue(f6, V, x)
            % f**(x) = sup_s <s,x> - f*(s), by coarse grid then Nelder-Mead (f* is convex in s).
            obj = @(s) -(s(:)'*x(:) - biconjCPLQTest.fstarExact(f6, V, s));
            [S1,S2] = meshgrid(linspace(-15,15,31), linspace(-15,15,31));
            best = inf; bs = [0 0];
            for k = 1:numel(S1)
                val = obj([S1(k) S2(k)]);
                if val < best, best = val; bs = [S1(k) S2(k)]; end
            end
            opts = optimset('TolX',1e-13,'TolFun',1e-15,'MaxFunEvals',20000,'MaxIter',20000);
            [bs, fv]  = fminsearch(obj, bs, opts);
            [~,  fv2] = fminsearch(obj, bs, opts);
            v = -min(fv, fv2);
        end
    end
end
