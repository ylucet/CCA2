classdef kinkSubdivisionTest < matlab.unittest.TestCase
% The poly2 kink at c = -2 (DECISIONS 2026-09-05): q0 = x^2 + y^2 below the diagonal, q1 = q0 - 2x
% above it. REVIEW 2026-11-19 #10: the two pieces do NOT agree on the shared diagonal -- the
% perturbation -2x vanishes on the line x = 0, not on the diagonal 2x - 3y + 4 = 0 (a jump of 4 at
% (-2,0), 2 at (1,2)) -- so f is DISCONTINUOUS, and the exact engine now refuses it by name
% (PLQ:notContinuous) at the entry of BOTH conj and biconj. The suite used to pin the old biconjQ
% answer (the 317-face AlgCrv); that answer is exactly the class finding #10 flags as wrong (12 of
% 24 answered D-family biconjugates above a certified upper bound).
%
% What survives: the two pieces' OWN conjugates are single-face and stay in scope, so the
% envelope-assembly test below still runs on them (it exercises envelopeAssemble/envelopeEval,
% which the fix does not touch).
%
% BUCKET: slow -- the two single-face conjQs are computed once per session.

    methods (Test)
        function biconjQRefusesTheDiscontinuousKink(testCase)
        % Review 2026-11-19 #10: a jump across the shared diagonal is out of the exact engine's
        % scope (the paper's standing assumption is continuity), refused by name at the entry of
        % both operators -- no longer answered (the old answer was the 317-face AlgCrv).
            f = kinkSubdivisionTest.kink();
            testCase.verifyError(@() conjQ(f), 'PLQ:notContinuous', ...
                'a discontinuous two-piece input is out of scope at the entry of conj');
            testCase.verifyError(@() biconjQ(f), 'PLQ:notContinuous', ...
                'a discontinuous two-piece input is out of scope at the entry of biconj');
        end

        function anIRRATIONALStraightTieGetsItsCorners(testCase)
        % Found 2026-10-11 on the poly2 kink at c = -2 (biconjQ's own fixture): one cell pair ties
        % on a PARALLEL PAIR of IRRATIONAL lines, (s1 - 2 s2) = -6 +- sqrt(80). A straight tie is
        % treated as an ordinary facet, but an irrational line cannot be a row, so nothing produced
        % candidates ON it -- and both faces are linear along it, so the pair's sup sits at its
        % ends, the corners. The pair covered no point at all and the envelope was +inf on 384 of
        % 600 points of its own domain. Oracle: the numeric sup, as above.
        % (Unaffected by #10: it tests the envelope of the two SINGLE-FACE conjugates, which stay
        % in scope.)
            [V, gA, gB] = kinkSubdivisionTest.cached();
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
        function f = kink()
        % The discontinuous two-piece kink (see the class header).
            V = [-2 0; 2 0; 1 2; -1 2];
            q0 = [0 0 0 0, 2, 0, 2, 0, 0, 0];  q1 = q0;  q1(9) = -2;
            f = QuaPol(V, [1 2 1; 2 3 1; 3 4 1; 4 1 1; 1 3 1], [q0; q1], [1 0; 1 0; 2 0; 2 0; 2 1]);
        end

        function [V, gA, gB] = cached()
        % The two pieces' OWN conjugates: each face alone is a single-face QuaPol, still in scope
        % (biconjQ/conjQ on the DISCONTINUOUS pair is refused -- biconjQRefusesTheDiscontinuousKink).
            persistent C
            if isempty(C)
                V = [-2 0; 2 0; 1 2; -1 2];
                q0 = [0 0 0 0, 2, 0, 2, 0, 0, 0];  q1 = q0;  q1(9) = -2;
                gA = conjQ(QuaPol(V([1 2 3],:), [1 2 1; 2 3 1; 3 1 1], q0, [1 0;1 0;1 0]));
                gB = conjQ(QuaPol(V([1 3 4],:), [1 2 1; 2 3 1; 3 1 1], q1, [1 0;1 0;1 0]));
                C = struct('V', V, 'gA', gA, 'gB', gB);
            end
            V = C.V;  gA = C.gA;  gB = C.gB;
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
