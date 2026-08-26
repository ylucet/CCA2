classdef conjSymFreeTest < matlab.unittest.TestCase
% Pins WHICH ROUTE `conj` takes, not only what it returns.
%
% BUCKET: fast. Every case here is answered by the closed-form numeric path in well under a second,
% and the cases that are NOT are asserted to refuse immediately rather than to run the symbolic
% pipeline -- which is what keeps this suite fast even though it covers the fallbacks.
%
% WHY A ROUTE TEST EXISTS AT ALL. "conj is sym-free except as a fallback" is a claim about a RATE,
% and under the default route a numeric refusal is SWALLOWED: the symbolic Case C answers, the
% values are still right, and nothing in the test suite notices that an input moved onto the slow
% path. Measured, that path costs 86-112 s where the numeric one costs 0.01-1 s, and it needs the
% Symbolic Toolbox and the licence server. So a regression that pushes a shape onto it is exactly
% the kind of thing that goes unnoticed for a month.
%
% `route='numeric'` refuses to fall back, so each case below states its route as an assertion. The
% VALUES are checked in conjCPLQTest and conjConvexPolygonTest; what is checked here is the ROUTE,
% plus one definition-oracle spot check per shape family so a route test cannot pass on a wrong
% answer.
%
% HOW TO READ A FAILURE.
%   * a `mustBeNumeric` case gone red = something pushed a supported shape back onto the symbolic
%     path. Find the reason with `checkConjSymFree`, which reports the identifier.
%   * a `stillFallsBack` case gone GREEN is good news, not a failure: move it up and delete its
%     entry from the gap list in TODO.md.

    methods (Test)

        % ---- the shapes the numeric route covers ------------------------------------------------

        function fullDomainQuadraticsAreNumeric(testCase)
            for q = { QuaPol([1 0 1 0 0 0]), QuaPol([2 1 3 -1 4 5]) }
                g = conjCPLQ(q{1}, [], 'numeric');
                testCase.verifyEqual(g.kind(), 'QuaPol');
            end
        end

        function everySingleTriangleIsNumeric(testCase)
        % One per Hessian class: affine, convex PD, rank-1 PSD, concave, indefinite (0, 1 and 3
        % convex edges). The A.5 three-convex-edge triangle is included on purpose -- it is the
        % most expensive numeric case and the one most likely to regress onto the symbolic path.
            T  = [0 0; 1 0; 0 1];
            cases = { {T, [0 0 0 1 2 3]},  {T, [2 0 2 0 0 0]},  {T, [1 0 0 0 0 0]}, ...
                      {T, [-2 0 -2 0 0 0]}, {T, [0 1 0 0 0 0]}, {T, [1 0 -1 0 0 0]}, ...
                      {[1 1; 0 0; 2 0], [0 1 0 0 0 0]}, ...
                      {[5/2 3/2; 0 0; 1/2 1], [0 1 0 0 0 0]} };
            for k = 1:numel(cases)
                q = conjSymFreeTest.tri(cases{k}{1}, cases{k}{2});
                g = conjCPLQ(q, [], 'numeric');
                testCase.verifyTrue(g.isMeshed(), sprintf('triangle case %d returned no mesh', k));
                testCase.verifyFalse(isa(g, 'QuaParCPLQ'), sprintf('triangle case %d went symbolic', k));
            end
        end

        function boundedPolygonsWithOneQuadraticAreNumeric(testCase)
            S = [0 0; 1 0; 1 1; 0 1];
            for c = { {S, [2 0 2 0 0 0]}, {S, [-2 0 -2 0 0 0]}, {S, [0 0 0 1 2 3]} }
                q = conjSymFreeTest.poly4(c{1}{1}, c{1}{2});
                g = conjCPLQ(q, [], 'numeric');
                testCase.verifyFalse(isa(g, 'QuaParCPLQ'));
            end
        end

        function aConvexPolygonTakesTheWHOLEFacePathAndComesBackPolyhedral(testCase)
        % The specific gain of 2026-08-24: a convex face is no longer fan-triangulated, so its
        % conjugate carries no parabolic edge at all. A QuaPol return type IS that statement, and
        % it is what keeps the arc-vs-arc machinery out of Step 3 for convex input.
            q = conjSymFreeTest.poly4([0 0; 1 0; 1 1; 0 1], [2 0 2 0 0 0]);
            g = conjCPLQ(q, [], 'numeric');
            testCase.verifyEqual(g.kind(), 'QuaPol');
            testCase.verifyTrue(all(g.Ec(:) == 0));
        end

        function anUNBOUNDEDConvexMultiFaceDomainIsNumericAndCorrect(testCase)
        % Before 2026-08-24 this shape had NO numeric route -- conjCPLQ's Case B2 was gated on
        % isDomBounded and everything unbounded went to the symbolic Case C. The values are checked
        % against a closed form, because each face is separable on its own quadrant: the sup of
        % s*t - a*t^2/2 over a half-line is s^2/(2a) when s points into it and 0 otherwise.
            V = [0 0; -1 0; 0 1; 1 0; 0 -1];
            E = [1 2 0; 1 3 0; 1 4 0; 1 5 0];
            f = [1 0 1 0 0 0; 1 0 2 0 0 0; 2 0 2 0 0 0; 2 0 1 0 0 0];
            F = [1 2; 2 3; 3 4; 4 1];
            q = QuaPol(V, E, f, F);
            testCase.verifyFalse(q.isDomBounded);
            g = conjCPLQ(q, [], 'numeric');
            testCase.verifyFalse(isa(g, 'QuaParCPLQ'), 'this must not need the symbolic path');
            testCase.verifyTrue(g.isMeshed());

            % Which quadrant each face occupies is DERIVED from the mesh rather than hardcoded:
            % getting the face-to-cone assignment wrong by hand silently compares against the
            % wrong closed form. Same construction as conjCPLQTest's own version of this fixture.
            S = [0 0; 1 1; -1 -1; 2 0.5; -0.5 2; 3 -1; -2 3; 0.25 -0.75];
            for t = 1:size(S,1)
                s = S(t,:);
                best = -inf;
                for k = 1:4
                    ax = [f(k,1), f(k,3)];
                    sg = conjSymFreeTest.quadrantSignOf(V, E, F, k);
                    v = 0;
                    for d = 1:2
                        if sign(s(d)) == sg(d), v = v + s(d)^2/(2*ax(d)); end
                    end
                    best = max(best, v);
                end
                testCase.verifyEqual(g.eval(s), best, 'AbsTol', 1e-9, ...
                    sprintf('at s = (%g,%g)', s(1), s(2)));
            end
        end

        % ---- the boundary of what is covered, pinned so it is not mistaken for a bug ------------

        function anIndefiniteQuadraticOverAPolygonStillFallsBack(testCase)
        % TWO indefinite triangles conjugate to two PARABOLIC QuaPars, and folding them puts one
        % operand's face vertex in the open interior of the other's arc. Since 2026-08-24 that is
        % SPLIT rather than refused (maxQuaPar's bulgeSplit / the passthrough split), but the
        % pieces then fail to pair up in assemblePieces -- the subdivision is consistent per face
        % PAIR and not yet globally. DECISIONS.md 2026-08-24 has the design for closing it.
        %
        % Asserted as an identifier PREFIX, not an exact message: what is being pinned is that the
        % numeric route declines with a maxQuaPar reason, not today's wording of it.
            q = conjSymFreeTest.twoTri([0 0; 2 0; 5/2 1; 1/2 1], [0 1 0 0 0 0]);
            try
                conjCPLQ(q, [], 'numeric');
                testCase.verifyFail(['the indefinite-over-a-polygon case now completes ' ...
                    'numerically. That is GOOD -- promote this test and update TODO.md.']);
            catch ME
                testCase.verifyTrue(strncmp(ME.identifier, 'maxQuaPar:', 10), ...
                    sprintf('declined with %s, which is not a maxQuaPar reason', ME.identifier));
            end
            % The DEFAULT route still answers this correctly through the symbolic fallback -- that
            % is the point of having one -- and conjCPLQTest checks those values. It is NOT checked
            % here on purpose: the symbolic path costs 86 s against 0.01-1 s for the numeric one,
            % which would put this whole suite in the wrong bucket to pin a ROUTE.
        end

        function anAFFINEUnboundedFaceIsNowNUMERICAndCorrect(testCase)
        % max(0,x,y): convex, its own biconjugate, every face affine over an unbounded wedge.
        %
        % TEST CHANGED 2026-08-25 (overnight), and this is the sentence umbrella CLAUDE.md 8 asks
        % for. It was `anAFFINEUnboundedFaceStillFallsBack`, and it asserted that this input must
        % DECLINE on the numeric route -- because the conjugate of an affine function over an
        % unbounded polygon is a SUPPORT function, +inf off a cone, and `maxQuaPar` refuses an
        % operand that is not finite everywhere. That was a true statement about a GAP, TODO.md's
        % G2, and the gap is now closed: `conjAffinePLQ` builds the answer straight from the
        % definition and never enters `maxQuaPar` at all. A test that pins a gap has to change when
        % the gap does, or it forbids the fix.
        %
        % It now asserts the thing worth having, which the old one could not: the numeric route
        % SUCCEEDS and the value is right. f*(s) = sup_x <s,x> - max(0,x,y) is the INDICATOR of the
        % unit simplex -- +inf unless s >= 0 (push along either axis) and unless s1 + s2 <= 1 (push
        % along (1,1)), and 0 on it, attained at the origin.
            V = [0 0; 0 -1; -1 0; 1 1];
            E = [1 2 0; 1 3 0; 1 4 0];
            F = [2 1; 1 3; 3 2];
            f = [0 0 0 0 0 0; 0 0 0 1 0 0; 0 0 0 0 1 0];
            q = QuaPol(V, E, f, F);
            g = conjCPLQ(q, [], 'numeric');
            testCase.verifyFalse(isa(g, 'QuaParCPLQ'), 'the affine route must not go symbolic');
            for s = {[0 0], [0.3 0.3], [0.25 0.5], [0.9 0.05]}
                testCase.verifyEqual(g.eval(s{1}), 0, 'AbsTol', 1e-9, sprintf( ...
                    'f*(%g,%g) must be 0 on the simplex', s{1}(1), s{1}(2)));
            end
            for s = {[0.6 0.6], [-0.1 0.5], [0.5 -0.1], [2 0]}
                v = g.eval(s{1});
                testCase.verifyTrue(isinf(v) || isnan(v), sprintf( ...
                    'f*(%g,%g) must be +inf off the simplex, got %g', s{1}(1), s{1}(2), v));
            end
        end

        function routeSYMBOLICHasSomewhereToGoForASingleTriangle(testCase)
        % TODO.md G12. `route='symbolic'` is documented as skipping the numeric path, and Case B
        % used to ignore it entirely: a single bounded triangle came back as the same numeric mesh
        % whatever the caller asked for. That is not only a testability point -- `biconjCPLQ` asks
        % for the symbolic form exactly when the numeric first conjugate is a CURVED mesh, because
        % the second conjugation cannot take one, and for a triangle it was handed the same curved
        % mesh back. The escape did nothing.
        %
        % Gating Case B on ~forceSymbolic so it fell through to Case C is REFUTED (DECISIONS.md
        % 2026-08-25 G12): Case C does not cover a single triangle and raises cplqFailed after
        % 102 s. The destination is cPLQ's own PER-TRIANGLE form, `conjEnvelopeViaCPLQ`, which is
        % what `conjSingleTriangle` already falls back to.
        %
        % The value is asserted, not just the type: `x*y` over {(0,0),(1,0),(0,1)} conjugates to
        % max(0,s1,s2), whose closed form is exact here.
            q = conjSymFreeTest.tri([0 0; 1 0; 0 1], [0 1 0 0 0 0]);
            gs = conjCPLQ(q, [], 'symbolic');
            testCase.verifyEqual(class(gs), 'QuaParCPLQ', ...
                'route=symbolic on a triangle must return cPLQ''s symbolic form');
            gn = conjCPLQ(q, [], 'numeric');
            testCase.verifyFalse(isa(gn, 'QuaParCPLQ'), 'route=numeric must stay a mesh');
            S = [0 0; 1 1; -1 -1; 2 0.5; -0.5 2; 1 -1];
            for i = 1:size(S,1)
                s = S(i,:);
                testCase.verifyEqual(evalFunctionNDomain(gs.fnd, s), ...
                    max([0, s(1), s(2)]), 'AbsTol', 1e-9, sprintf( ...
                    'the symbolic form must equal max(0,s1,s2) at (%g,%g)', s(1), s(2)));
            end
        end

        function theRouteArgumentItselfIsValidated(testCase)
            q = QuaPol([1 0 1 0 0 0]);
            testCase.verifyError(@() conjCPLQ(q, [], 'nonsense'), 'PLQ:conjCPLQ:route');
        end
    end

    methods (Static)
        function q = tri(V, f6)
            q = QuaPol(V, [1 2 1; 2 3 1; 3 1 1], f6, [1 0; 1 0; 1 0]);
        end

        function q = poly4(V, f6)
            q = QuaPol(V, [1 2 1; 2 3 1; 3 4 1; 4 1 1], f6, [1 0; 1 0; 1 0; 1 0]);
        end

        function q = twoTri(V, f6)
        % A quadrilateral as two triangles carrying the SAME quadratic. Step 0 merges them into one
        % face, so what this exercises is one indefinite quadratic over a four-sided polygon.
            q = QuaPol(V, [1 2 1; 2 3 1; 1 3 1; 3 4 1; 4 1 1], [f6; f6], ...
                       [1 0; 1 0; 2 1; 2 0; 2 0]);
        end

        function sg = quadrantSignOf(V, E, F, k)
        % objective: which quadrant face k of the four-wedge fixture occupies, as a sign pair,
        %   derived from the mesh's own ray directions rather than hardcoded.
            ds = zeros(0,2);
            for j = 1:size(E,1)
                if any(F(j,:) == k)
                    ds(end+1,:) = V(E(j,2),:) - V(E(j,1),:); %#ok<AGROW>
                end
            end
            m = mean(ds, 1);
            sg = sign(m);
            sg(sg == 0) = 1;
        end
    end
end
