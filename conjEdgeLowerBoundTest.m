classdef conjEdgeLowerBoundTest < matlab.unittest.TestCase
% Unit tests for conjEdgeLowerBound and for the refusal conjCPLQ builds on it.
%
% BUCKET: fast (closed form throughout; the only loop is over edges).
%
% WHAT THE BOUND IS. For ANY subset A of dom f, f*(s) = sup_{x in dom f} <s,x> - f(x) is at least
% the sup over A. Taking A to be one EDGE makes that a one-dimensional quadratic maximisation with
% a closed form, and the best over all edges is the bound. It is therefore SOUND by construction --
% the tests below check that it is also TIGHT where it should be, and that the refusal built on it
% fires where it must and nowhere else.

    methods (Test)

        function theBoundIsTIGHTWhereTheSupIsOnTheBOUNDARY(testCase)
        % For a CONCAVE or INDEFINITE piece the maximiser of <s,x> - q(x) is on the boundary, so
        % the edge bound is not merely valid there but EQUAL to f*. That is what makes it a useful
        % check rather than a slack one.
            W = [0 0; 1 0; 1 1; 0 1];
            for f6 = { [-2 0 -2 0 0 0], [0 1 0 0 0 0], [0 0 0 1 -2 0.5] }
                q = conjEdgeLowerBoundTest.poly(W, f6{1});
                g = q.conj('cplq');
                S = conjEdgeLowerBoundTest.probes();
                lb = conjEdgeLowerBound(q, S);
                for i = 1:size(S,1)
                    testCase.verifyEqual(g.eval(S(i,:)), lb(i), 'AbsTol', 1e-9, ...
                        sprintf('at s=(%g,%g)', S(i,1), S(i,2)));
                end
            end
        end

        function theBoundIsVALIDButSLACKWhereTheSupIsINTERIOR(testCase)
        % A CONVEX piece can attain its sup strictly inside the domain, and there the boundary says
        % strictly less than f*. Asserted as an INEQUALITY plus a witness that the gap is real, so
        % the test states the limitation rather than pretending the bound is tight everywhere.
            W = [0 0; 1 0; 1 1; 0 1];
            q = conjEdgeLowerBoundTest.poly(W, [2 0 2 0 0 0]);
            g = q.conj('cplq');
            S = conjEdgeLowerBoundTest.probes();
            lb = conjEdgeLowerBound(q, S);
            gap = 0;
            for i = 1:size(S,1)
                testCase.verifyGreaterThanOrEqual(g.eval(S(i,:)), lb(i) - 1e-9, ...
                    sprintf('the bound must never exceed f*, at s=(%g,%g)', S(i,1), S(i,2)));
                gap = max(gap, g.eval(S(i,:)) - lb(i));
            end
            testCase.verifyGreaterThan(gap, 1e-6, ...
                'on a convex piece the bound should be genuinely slack somewhere');
        end

        function aRAYIsParametrisedToInfinityAndItsMARKERIsNotADomainPoint(testCase)
        % The second stored point of a ray is a DIRECTION MARKER, not a point of the domain
        % (RatCon.m's `E`). Reading it as one would make the bound too large on the marker's side,
        % which would turn a correct conjugate into a reported defect -- the worst possible failure
        % for a check. Pinned on the first quadrant, whose conjugate is known in closed form.
            % F, spelled out because getting it backwards is easy and silent: edge 1 runs along
            % +x with the quadrant ABOVE it, so face 1 is on its LEFT -> [1 0]; edge 2 runs along
            % +y with the quadrant to its RIGHT -> [0 1].
            V = [0 0; 1 0; 0 1];
            q = QuaPol(V, [1 2 0; 1 3 0], [2 0 2 0 0 0], [1 0; 0 1]);
            S = [1 1; -1 -1; 3 0.5; -2 4];
            lb = conjEdgeLowerBound(q, S);
            for i = 1:size(S,1)
                s = S(i,:);
                % sup over the two boundary RAYS of {x>=0,y>=0} of <s,x> - (x^2+y^2)
                b = max(rayMax(s, [1 0]), rayMax(s, [0 1]));
                testCase.verifyEqual(lb(i), b, 'AbsTol', 1e-9, ...
                    sprintf('at s=(%g,%g)', s(1), s(2)));
            end
            function m = rayMax(s, d)
                % max over t >= 0 of  a1*t + a2*t^2  with a1 = <s,d> and a2 = -(d1^2+d2^2) < 0.
                % t = 0 gives 0; the stationary point t* = -a1/(2*a2) is positive exactly when
                % a1 > 0, and there the value is a1^2/(-4*a2).
                a2 = -(d(1)^2 + d(2)^2);
                a1 = s*d.';
                m = 0;
                if a1 > 0, m = a1^2 / (-4*a2); end
            end
        end

        function theKnownBadFixtureIsREFUSEDByTheNumericRoute(testCase)
        % G4's fixture: `conj`'s numeric answer on it is wrong, and the numeric route refuses it
        % by name rather than returning it. That is the property worth pinning, and it is now
        % pinned on the NUMERIC route specifically -- `conjCPLQ(q, [], 'numeric')` -- for a
        % measured reason.
        %
        % HISTORY, 2026-08-25, two changes in one session and the second undid the first's test.
        % This started as `theRefusalFIRESOnTheKnownMinorantAndNowhereElse`, asserting
        % `PLQ:conjCPLQ:belowEdgeBound`. Then the fold cross-check landed and got there first
        % (`PLQ:conjCPLQ:foldDroppedACell`), so it was widened to accept either. Then the actual
        % maxQuaPar defects were fixed -- a recession cone decided by 1e-17 noise, and an
        % unbounded cell whose winner was read from its centroid -- and now maxQuaPar itself
        % refuses this fixture EARLIER than either definition check, with
        % `maxQuaPar:notImplemented`: an unbounded piece genuinely straddles {f1 = f2} and has no
        % single winner to carry.
        %
        % WHY THE `auto` ROUTE IS NOT USED HERE. `conjSingleTriangle` swallows `maxQuaPar:` and
        % falls through to the symbolic Case C, which on this fixture runs for 420 s and then
        % fails in MuPAD (`PLQ:conjCPLQ:cplqFailed`). That is a true statement about the symbolic
        % route and it belongs in the slow bucket, not in a fast test that would spend seven
        % minutes reaching it. The numeric route refuses in about 3 s.
        %
        % NOTE FOR WHOEVER FIXES THE STRADDLE: when maxQuaPar can assemble this cell, the fold
        % cross-check and the edge lower bound become reachable on it again, and the expected
        % identifier here moves back. Widen the list, do not delete the test.
            W = [0.6057047151 0.9300751811; -0.3353947472 0.5251524293; -1.082499617 0.08448609744];
            f6 = [0 1 0 -0.7177913413 -0.6075645347 -0.6781835233];
            q = QuaPol(W, [1 2 1; 2 3 1; 3 1 1], f6, [1 0; 1 0; 1 0]);
            try
                g = conjCPLQ(q, [], 'numeric');
                testCase.verifyFail(sprintf( ...
                    ['the known-bad fixture returned a %s: the numeric route neither refused it ' ...
                     'nor was caught by a definition check'], class(g)));
            catch ME
                testCase.verifyTrue(any(strcmp(ME.identifier, ...
                    {'maxQuaPar:notImplemented', 'PLQ:conjCPLQ:foldDroppedACell', ...
                     'PLQ:conjCPLQ:belowEdgeBound', 'PLQ:conjCPLQ:numericRouteDeclined'})), ...
                    sprintf(['refused with %s, which is not one of the refusals this fixture may ' ...
                             'legitimately produce'], ME.identifier));
            end
        end

        function theDefinitionChecksCanBeTurnedOffWithTheGlobal(testCase)
        % `CCA2_CONJ_VERIFY = 0` turns off the DEFINITION checks -- the edge lower bound and the
        % fold cross-check -- and with them off the old, unchecked answer comes back. That is what
        % the switch is for: a caller who wants the previous behaviour, or who is debugging the
        % wrong answer itself, needs to be able to obtain it.
        %
        % Renamed from `theRefusalCanBeTurnedOffWithTheGlobal` and moved onto the NUMERIC route.
        % Both are for speed and precision, not for a change of meaning: under `auto` this fixture
        % now spends 420 s in the symbolic fallback before failing there, which is a true statement
        % about the symbolic route and does not belong in a fast test. `conjCPLQ(q, [], 'numeric')`
        % asks the question this test is actually about -- what the numeric path returns when the
        % definition checks are off -- in about 3 s.
            prev = getGlobal();
            c = onCleanup(@() setGlobal(prev)); %#ok<NASGU>
            W = [0.6057047151 0.9300751811; -0.3353947472 0.5251524293; -1.082499617 0.08448609744];
            f6 = [0 1 0 -0.7177913413 -0.6075645347 -0.6781835233];
            q = QuaPol(W, [1 2 1; 2 3 1; 3 1 1], f6, [1 0; 1 0; 1 0]);
            setGlobal(0);
            g = conjCPLQ(q, [], 'numeric');
            testCase.verifyTrue(isa(g, 'RatPar'), 'with the checks off the old answer comes back');
        end
    end

    methods (Static)
        function q = poly(W, f6)
            n = size(W,1);
            q = QuaPol(W, [(1:n)', mod((1:n),n)'+1, ones(n,1)], f6, [ones(n,1), zeros(n,1)]);
        end

        function S = probes()
            S = [0 0; 2 2; -1 3; 0.3 0.4; 5 -2; -4 -4; 1.5 0.5; -3 1];
        end
    end
end

function v = getGlobal()
    global CCA2_CONJ_VERIFY %#ok<GVMIS>
    v = CCA2_CONJ_VERIFY;
end

function setGlobal(v)
    global CCA2_CONJ_VERIFY %#ok<GVMIS>
    CCA2_CONJ_VERIFY = v;
end
