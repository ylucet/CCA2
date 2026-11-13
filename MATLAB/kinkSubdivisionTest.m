classdef kinkSubdivisionTest < matlab.unittest.TestCase
% biconjQ on the poly2 kink at c = -2: the two-piece envelope whose single global arrangement had 85
% curves and did not finish in 90 min (TODO item 2). Since 2026-10-17 it is SUBDIVIDED: the arrangement is
% clipped to dom f** = conv(dom f), a curve's row is recorded only where it splits a cell, and a
% bridge's membership is read by majority -- envelopeSubdivisionAlg's header.
%
% BUCKET: slow -- one biconjQ, ~11 min on the shared machine (one run), computed once per session.

    methods (Test)
        function theKinkIsSUBDIVIDEDAndAgreesWithItsFaceList(testCase)
        % Subdivision against the unsubdivided face list, which envelopeAssembleTest verifies against the
        % definition: the subdivision must add no hole and change no value.
            [h, env, V] = kinkSubdivisionTest.cached();
            testCase.verifyClass(h, 'AlgCrv');
            rng(4);  W = [rand(1500,4); rand(500,4).^4];  X = (W ./ sum(W,2)) * V;
            want = envelopeEval(env, X);
            got = h.eval(X);
            testCase.verifyEqual(isinf(got), isinf(want), 'a hole in the subdivision where the face list is finite');
            fin = isfinite(want);
            testCase.verifyEqual(got(fin), want(fin), 'RelTol', 1e-9, 'AbsTol', 1e-9);
        end

        function theKinkSatisfiesTheDEFINITION(testCase)
        % Oracle: the numeric sup over the two pieces' conjugates (moved from envelopeAssembleTest,
        % where this asserted the unsubdivided answer).
            [h, ~, V, gA, gB] = kinkSubdivisionTest.cached();
            rng(29);  W = rand(6,3);  W = W ./ sum(W,2);
            X = [W*V([1 2 3],:); W*V([1 3 4],:)];
            v = h.eval(X);
            for i = 1:size(X,1)
                want = kinkSubdivisionTest.supOracle(gA, gB, X(i,:));
                testCase.verifyLessThan(abs(v(i) - want)/max(1,abs(want)), 1e-7, ...
                    sprintf('biconjQ must equal the definition at (%g,%g)', X(i,1), X(i,2)));
            end
        end
    end

    methods (Test)
        function anIRRATIONALStraightTieGetsItsCorners(testCase)
        % Found 2026-10-11 on the poly2 kink at c = -2 (biconjQ's own fixture): one cell pair ties
        % on a PARALLEL PAIR of IRRATIONAL lines, (s1 - 2 s2) = -6 +- sqrt(80). A straight tie is
        % treated as an ordinary facet, but an irrational line cannot be a row, so nothing produced
        % candidates ON it -- and both faces are linear along it, so the pair's sup sits at its
        % ends, the corners. The pair covered no point at all and the envelope was +inf on 384 of
        % 600 points of its own domain. Oracle: the numeric sup, as above.
            [~, ~, V, gA, gB] = kinkSubdivisionTest.cached();
            env = envelopeAssemble(gA, gB);
            rng(23);  W = rand(60,3);  W = W ./ sum(W,2);
            X = [W*V([1 2 3],:); W*V([1 3 4],:)];
            v = envelopeEval(env, X);
            testCase.verifyTrue(all(isfinite(v)), sprintf( ...
                'the envelope must be finite on its domain; +inf at %d of %d', nnz(~isfinite(v)), numel(v)));
            for i = 1:8:size(X,1)
                want = kinkSubdivisionTest.supOracle(gA, gB, X(i,:));
                testCase.verifyGreaterThanOrEqual(v(i), want - 1e-9*max(1,abs(want)), 'no undershoot');
                testCase.verifyLessThan(abs(v(i) - want)/max(1,abs(want)), 1e-7, ...
                    sprintf('the envelope must equal the definition at (%g,%g)', X(i,1), X(i,2)));
            end
        end
    end

    methods (Static)
        function [h, env, V, gA, gB] = cached()
            persistent C
            if isempty(C)
                V = [-2 0; 2 0; 1 2; -1 2];
                q0 = [0 0 0 0, 2, 0, 2, 0, 0, 0];  q1 = q0;  q1(9) = -2;
                f = QuaPol(V, [1 2 1; 2 3 1; 3 4 1; 4 1 1; 1 3 1], [q0; q1], [1 0; 1 0; 2 0; 2 0; 2 1]);
                gA = conjQ(QuaPol(V([1 2 3],:), [1 2 1; 2 3 1; 3 1 1], q0, [1 0;1 0;1 0]));
                gB = conjQ(QuaPol(V([1 3 4],:), [1 2 1; 2 3 1; 3 1 1], q1, [1 0;1 0;1 0]));
                C = struct('h', biconjQ(f), 'env', envelopePrune(envelopeAssemble(gA, gB)), ...
                           'V', V, 'gA', gA, 'gB', gB);
            end
            h = C.h;  env = C.env;  V = C.V;  gA = C.gA;  gB = C.gB;
        end

        function v = supOracle(gA, gB, x)
        % sup_s [ <s,x> - max(gA,gB)(s) ]: a grid for the basin, then Nelder-Mead from two starts.
            obj = @(s) -( s(:).'*x(:) - max(ev(gA,s), ev(gB,s)) );
            best = inf;  bs = [0;0];
            for s1 = -8:0.1:8
                for s2 = -8:0.1:8
                    o = obj([s1;s2]);
                    if o < best, best = o;  bs = [s1;s2]; end
                end
            end
            op = optimset('TolX',1e-14,'TolFun',1e-14,'MaxFunEvals',20000,'Display','off');
            s1 = fminsearch(obj, bs, op);
            s2 = fminsearch(obj, bs + [1e-3; -1e-3], op);
            v = -min([obj(s1), obj(s2), best]);
        end
    end
end

function v = ev(g, s)
    v = g.eval(s(:).');
    if isempty(v) || ~isfinite(v), v = -inf; end
end
