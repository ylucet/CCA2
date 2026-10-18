classdef minMeshTest < matlab.unittest.TestCase
% biconjQ on an indefinite triangle with THREE convex edges: three edge atoms, two ACTIVE pairs and
% no triangle face, so co f is the minimum of two pair envelopes -- meshed since 2026-10-18 as a min
% of maxes in one arrangement (envelopeMeshAlg's minMaxFaces).
%
% BUCKET: slow -- biconjQ here is ~5 min on the shared machine (one run), computed once per session.

    methods (Test)
        function theMinimumOverTWOPairsIsMESHED(testCase)
            [h, f, V] = minMeshTest.cached();
            testCase.verifyClass(h, 'AlgCrv');
            want = biconjAtoms(f);
            testCase.verifyEqual(numel(want.pairs), 2, 'the fixture must have two active pairs');
            testCase.verifyEmpty(want.tris, 'and no triangle face');
            rng(3);  W = [rand(1500,3); rand(500,3).^4];  X = (W ./ sum(W,2)) * V;
            got = h.eval(X);  ref = want.eval(X);
            testCase.verifyEqual(isinf(got), isinf(ref), 'a hole where the minimum is finite');
            fin = isfinite(ref);
            % within 1e-8: points 1e-10 from an edge read 3e-9 off, the evaluators' row tolerance
            testCase.verifyEqual(got(fin), ref(fin), 'AbsTol', 1e-8, 'RelTol', 1e-8);
        end

        function itAgreesWithTheLEGACYEnvelope(testCase)
        % The differential erratumBiconjTest runs on the other triangles (umbrella CLAUDE.md
        % section 6): convEnvCPLQ shares no code with the atom route or the mesh.
            [h, f, V] = minMeshTest.cached();
            rng(5);  W = rand(300,3);  X = (W ./ sum(W,2)) * V;
            testCase.verifyEqual(h.eval(X), f.biconj('cplq').eval(X), 'AbsTol', 1e-10);
        end
    end

    methods (Static)
        function [h, f, V] = cached()
            persistent C
            if isempty(C)
                V = [0 0; 1 0; 0 1];
                f = QuaPol(V, [1 2 1; 2 3 1; 3 1 1], [0 0 0 0 2 -3 2 0 0 0], [1 0; 1 0; 1 0]);
                C = struct('h', biconjQ(f), 'f', f, 'V', V);
            end
            h = C.h;  f = C.f;  V = C.V;
        end
    end
end
