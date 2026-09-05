classdef isConvexQTest < matlab.unittest.TestCase
% Unit tests for isConvexQ, the EXACT convexity test for a QuaPol (Singh & Lucet, SIOPT 2021,
% Theorem 3.7, resting on Bauschke, Lucet & Phan, ESAIM COCV 2016, Theorem 5.1).
%
% BUCKET: fast (integer sign tests only; no symbolic call anywhere).
%
% THE ORACLE IS THE DEFINITION, not the legacy routine. `f` is convex iff
% `f(lambda x + (1-lambda) y) <= lambda f(x) + (1-lambda) f(y)` for every pair, and that is what
% `sampledViolation` below recomputes -- in doubles, from `eval` alone, with no reference to
% isConvexQ's algebra. Sampling is ONE-SIDED and used only in the direction where it proves
% something: a violation is a definite refutation, while finding none is not a proof. It is
% therefore applied to every f the test calls CONVEX, which is the dangerous answer -- a wrong
% `true` makes biconjQ hand back a non-convex f as its own convex envelope.
%
% The legacy doubles `QuaPol.isConvex` appears once, as a DIFFERENTIAL ORACLE (umbrella CLAUDE.md
% section 6) rather than as a reference: it has two latent smells recorded in TODO.md 2026-09-05,
% and on `examples{16}` it is simply wrong -- so that fixture is pinned here by an explicit
% witness rather than by agreement.

    methods (Test)

        function theSignConventionIsPinnedByTheOneNormAndItsNegative(testCase)
        % The one orientation that cannot be checked by inspection: which of F(j,1), F(j,2) plays
        % the theorem's P_i, the piece on the {a'x >= beta} side. Swapping them turns the test into
        % a CONCAVITY test that would pass every fixture below just as happily, so it gets a pair
        % whose only difference is the sign.
        %
        % |x| + |y| and its negative have the same mesh and identically ZERO Hessians, so the
        % Hessian test cannot separate them and the domain test cannot either (both are the plane).
        % Only the sign of delta can, which is why the reason string is asserted too.
            p = QuaPol.oneNorm();
            [tf, why] = isConvexQ(p);
            testCase.verifyTrue(tf, sprintf('|x|+|y| must be convex; got "%s"', why));

            q = QuaPol(p.V, p.E, -p.f, p.F);       % -(|x|+|y|), concave, same mesh
            [tf2, why2] = isConvexQ(q);
            testCase.verifyFalse(tf2, '-(|x|+|y|) must not be convex');
            testCase.verifySubstring(why2, 'subdifferential', ...
                'the concave twin must be rejected by the delta sign, not by some other condition');
        end

        function everyConvexAnswerSatisfiesTheDefinitionOnSampledSegments(testCase)
        % The definition, recomputed from eval. See the class header for why this runs on the
        % TRUE answers: a wrong true is the silent failure, a wrong false is only a missed
        % short-circuit.
            rng(20260905);
            L = isConvexQTest.corpus();
            nChecked = 0;
            for k = 1:numel(L)
                f = L{k}{2};
                tf = false;
                try, tf = isConvexQ(f); catch, continue, end   % cubic / inexact / thin: not ours
                if ~tf, continue, end
                nChecked = nChecked + 1;
                [bad, msg] = isConvexQTest.sampledViolation(f, 400);
                testCase.verifyFalse(bad, sprintf('%s is called convex but %s', L{k}{1}, msg));
            end
            testCase.verifyGreaterThan(nChecked, 8, ...
                'the corpus must contain enough convex fixtures for this to mean anything');
        end

        function theEightFixturesBiconjQRefusesAreConvexExceptTheOneThatIsNot(testCase)
        % TODO.md 2026-09-05 measured that all 8 fixtures biconjQ refuses become computable once
        % fIsConvex is asserted. Seven of them are genuinely convex and this test unblocks them.
        % The eighth, examples{16}, is NOT -- see the next test -- so asserting the flag there was
        % producing a wrong answer, and this routine must decline to bless it.
            P  = QuaPol.examples();
            P2 = QuaPol.examples2();
            convexOnes = {'oneNorm', QuaPol.oneNorm(); ...
                          'examples{1}',  P{1};  'examples{5}',  P{5};  'examples{8}',  P{8}; ...
                          'examples{18}', P{18}; 'examples{19}', P{19}; 'examples2{9}', P2{9}};
            for i = 1:size(convexOnes,1)
                [tf, why] = isConvexQ(convexOnes{i,2});
                testCase.verifyTrue(tf, sprintf('%s must be convex; got "%s"', convexOnes{i,1}, why));
            end
        end

        function theHalfStarDomainIsNotConvexAndTheWitnessProvesIt(testCase)
        % examples{16} is where the exact test and the legacy doubles one disagree, so it is
        % settled by the DEFINITION rather than by either of them: two points the mesh itself puts
        % in the domain whose midpoint the mesh puts outside it.
            f = QuaPol.examples();  f = f{16};
            a = [-1 -1];  b = [1 0];  mid = (a + b)/2;
            v = f.eval([a; b; mid]);
            testCase.verifyTrue(isfinite(v(1)) && isfinite(v(2)), ...
                'both witness points must be in the domain for the witness to mean anything');
            testCase.verifyEqual(isfinite(v(3)), false, ...
                'their midpoint must be outside the domain: that is the non-convexity');

            [tf, why] = isConvexQ(f);
            testCase.verifyFalse(tf, 'examples{16} must be rejected');
            testCase.verifySubstring(why, 'domain', ...
                'and rejected for the domain, which is what the witness exhibits');
        end

        function everyNonconvexFixtureIsRejected(testCase)
            L = QuaPol.examplesNonconvex();
            for k = 1:numel(L)
                [tf, why] = isConvexQ(L{k});
                testCase.verifyFalse(tf, ...
                    sprintf('examplesNonconvex{%d} must be rejected; got true', k));
                testCase.verifyNotEmpty(why, 'a false answer must say which condition failed');
            end
        end

        function everyDiscontinuousFixtureIsRejectedForDISCONTINUITY(testCase)
        % Not merely rejected: rejected by condition (b). A discontinuous f that happened to be
        % refused by the Hessian or the domain test would leave continuity untested.
            L = QuaPol.examplesDiscontinuous();
            nChecked = 0;
            for k = 1:numel(L)
                try, [tf, why] = isConvexQ(L{k}); catch, continue, end   % the cubic member
                nChecked = nChecked + 1;
                testCase.verifyFalse(tf, sprintf('examplesDiscontinuous{%d} must be rejected', k));
                testCase.verifySubstring(why, 'discontinuous', ...
                    sprintf('examplesDiscontinuous{%d} must fail the CONTINUITY condition', k));
            end
            testCase.verifyGreaterThan(nChecked, 2, 'too few discontinuous fixtures exercised');
        end

        function deltaMustNotDECREASEAlongARay(testCase)
        % Condition (d) on its own. Everything else about this function passes: both Hessians are
        % positive definite, the domain is the closed upper half-plane, f is continuous across the
        % splitting ray, and delta is nonnegative at BOTH stored endpoints of that ray. It is only
        % that delta decreases along the ray, so it goes negative further out.
        %
        %   f1 = x^2 + y^2 + xy - x  on the second quadrant,  f2 = x^2 + y^2  on the first.
        %
        % Then d = f1 - f2 = x(y-1), which vanishes on the splitting ray x = 0 (continuity), and
        % with a = (-1,0) gives delta proportional to 1 - y: nonnegative at (0,0) and at (0,1),
        % decreasing, hence negative for y > 1. Without condition (d) this f would be called
        % convex, and the witness below shows it is not.
            f = isConvexQTest.upperHalfPlaneSplit([2 1 2 -1 0 0], [2 0 2 0 0 0]);

            [tf, why] = isConvexQ(f);
            testCase.verifyFalse(tf, 'delta decreasing along the ray must be caught');
            testCase.verifySubstring(why, 'unbounded', ...
                'and caught by condition (d), not by the endpoint signs');

            % The definition, independently: at y = 3 the chord across the ray dips BELOW f.
            v = f.eval([-0.5 3; 0.5 3; 0 3]);
            testCase.verifyGreaterThan(v(3), (v(1)+v(2))/2 + 1e-12, ...
                'the witness must actually violate the midpoint inequality');
        end

        function theSameFunctionOnBOTHSidesIsConvexWhenTheQuadraticIs(testCase)
        % The companion of the previous test, and the reason condition (d) is not vacuous: with
        % delta identically zero the same mesh is convex, so the rejection above is about the jump
        % and not about the mesh being unbounded.
            f = isConvexQTest.upperHalfPlaneSplit([2 0 2 0 0 0], [2 0 2 0 0 0]);
            [tf, why] = isConvexQ(f);
            testCase.verifyTrue(tf, sprintf('x^2+y^2 on the half-plane is convex; got "%s"', why));
            [bad, msg] = isConvexQTest.sampledViolation(f, 300);
            testCase.verifyFalse(bad, msg);
        end

        function itAgreesWithTheLegacyDoublesTestWhereverTheLegacyIsRight(testCase)
        % The differential check (umbrella CLAUDE.md section 6). examples{16} is excluded BY NAME
        % and with its reason, not by a silent filter -- the test above proves the legacy answer
        % there is the wrong one.
            L = isConvexQTest.corpus();
            nCompared = 0;
            for k = 1:numel(L)
                name = L{k}{1};  f = L{k}{2};
                if strcmp(name, 'examples{16}'), continue, end   % legacy says convex; it is not
                try, new = isConvexQ(f); catch, continue, end
                try, old = logical(f.isConvex()); catch, continue, end
                nCompared = nCompared + 1;
                testCase.verifyEqual(new, old, ...
                    sprintf('%s: exact says %d, legacy doubles says %d', name, new, old));
            end
            testCase.verifyGreaterThan(nCompared, 15, ...
                'too few fixtures compared for the differential check to mean anything');
        end

        function whatIsOutOfScopeIsREFUSEDByNameRatherThanAnswered(testCase)
        % A cubic (the theorem used here is 3.7, the PLQ one -- 3.13 is a different algorithm), a
        % thin domain, an inexact object and a non-QuaPol. Each must raise its own identifier, so
        % that a caller can tell "not supported" from "not convex".
            P = QuaPol.examples();
            testCase.verifyError(@() isConvexQ(P{4}),  'PLQ:op:unsupportedType');   % cubic1
            testCase.verifyError(@() isConvexQ(P{17}), 'PLQ:QuaPol:notExact');          % unit circle
            P2 = QuaPol.examples2();
            testCase.verifyError(@() isConvexQ(P2{1}), 'PLQ:isConvexQ:thinDomain');
            testCase.verifyError(@() isConvexQ(42),    'PLQ:isConvexQ:input');
        end

    end

    methods (Static)

        function L = corpus()
        % Every fixture in the repository that a convexity test could be asked about, each with
        % the name it is known by, so a failure message names it.
            L = {};
            P = QuaPol.examples();
            for k = 1:numel(P), L{end+1} = {sprintf('examples{%d}',k), P{k}}; end %#ok<AGROW>
            P2 = QuaPol.examples2();
            for k = 1:numel(P2), L{end+1} = {sprintf('examples2{%d}',k), P2{k}}; end %#ok<AGROW>
            PN = QuaPol.examplesNonconvex();
            for k = 1:numel(PN), L{end+1} = {sprintf('examplesNonconvex{%d}',k), PN{k}}; end %#ok<AGROW>
            PD = QuaPol.examplesDiscontinuous();
            for k = 1:numel(PD), L{end+1} = {sprintf('examplesDiscontinuous{%d}',k), PD{k}}; end %#ok<AGROW>
        end

        function [bad, msg] = sampledViolation(f, n)
        % objective: search for a pair of points refuting convexity of f, using ONLY eval.
        % [output] bad : true if a violation was found; msg : which pair, for the failure message.
        %
        % A point outside the domain evaluates to +Inf, so a midpoint that leaves the domain shows
        % up as a violation of the same inequality -- the domain and the function are tested at
        % once, which is right, since f is convex exactly when epi f is.
            bad = false;  msg = 'no violation found on sampled segments';
            X = 6*rand(n,2) - 3;
            v = f.eval(X);
            in = isfinite(v);
            X = X(in,:);  v = v(in);
            m = size(X,1);
            if m < 2, return, end
            for t = 1:min(2000, m*(m-1)/2)
                i = randi(m);  j = randi(m);
                if i == j, continue, end
                lam = rand();
                z = lam*X(i,:) + (1-lam)*X(j,:);
                fz = f.eval(z);
                if isnan(fz), continue, end                    % eval reports a discontinuity
                rhs = lam*v(i) + (1-lam)*v(j);
                if fz > rhs + 1e-7*max(1,abs(rhs))
                    bad = true;
                    msg = sprintf(['convexity fails at lambda=%.4f between (%.4f,%.4f) and ' ...
                                   '(%.4f,%.4f): f=%.6g > %.6g'], lam, X(i,1), X(i,2), ...
                                   X(j,1), X(j,2), fz, rhs);
                    return
                end
            end
        end

        function f = upperHalfPlaneSplit(c1, c2)
        % The closed upper half-plane split by the ray x = 0, y >= 0 into the second quadrant
        % (face 1, coefficients c1) and the first (face 2, c2). Coefficients are in the 6-wide
        % weighted basis [x^2 xy y^2 x y 1] with weights [1/2 1 1/2 1 1 1], so [2 0 2 0 0 0] is
        % x^2 + y^2.
        %
        % Edge 1 runs (0,0) -> (0,1) and is the shared one: its direction is +y, so its LEFT is
        % the -x side and F(1,:) = [1 2]. Edges 2 and 3 are the two halves of the x-axis, with
        % +infinity on the outside.
            V = [0 0; 0 1; -1 0; 1 0];
            E = [1 2 0; 1 3 0; 1 4 0];
            F = [1 2; 0 1; 2 0];
            f = QuaPol(V, E, [c1; c2], F);
        end

    end
end
