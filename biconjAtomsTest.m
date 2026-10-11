classdef biconjAtomsTest < matlab.unittest.TestCase
% Unit tests for biconjAtoms' own contract. Its VALUES are checked in erratumBiconjTest (normal:
% erratum Example 1 and the differential against the legacy triangle envelope) and biconjQTest.
%
% BUCKET: fast (refusals only; nothing here builds a pair envelope).

    methods (Test)
        function theINPUTContractIsRefusedByName(testCase)
            testCase.verifyError(@() biconjAtoms(42), 'biconjAtoms:input');
            % CHANGED 2026-10-11: rays are in scope now (erratum Case E), so `biconjAtoms:unbounded`
            % is gone, and so (later the same night) is `unboundedConvexPiece`: a lone convex
            % unbounded piece is ONE atom with an empty D, sent back like a bounded one. What is
            % refused is an empty dom f* (co f = -infinity).
            cvxCone = QuaPol([0 0; 1 0; 0 1], [1 2 0; 1 3 0], [0 0 0 0, 2, 0, 2, 0, 0, 0], [1 0; 0 1]);
            testCase.verifyError(@() biconjAtoms(cvxCone), 'biconjAtoms:oneAtom');
            down = QuaPol([0 0; 1 0; 0 1], [1 2 0; 1 3 0], [0 0 0 0, -2, 0, -2, 0, 0, 0], [1 0; 0 1]);
            testCase.verifyError(@() biconjAtoms(down), 'biconjAtoms:emptyDomain');
        end

        function aTRIPLEFaceOnTheBoundaryOfDomFStarExtendsAlongItsNormalCone(testCase)
        % erratum Case E: a tangent plane <s,x> - c with s on the boundary of D = dom f* holds on
        % T + N_D(s), because moving x along an outward normal r of D at s costs exactly <s,r>.
        % EnvelopeMin must evaluate it there and nowhere else.
            T = [0 0; 1 0; 0 1];  s = [1; 2];  c = 0.5;
            one = EnvelopeMin({}, struct('T', T, 's', s, 'c', c, 'G', [0; 1]), {});
            X = [0.2 0.3; 0.2 5; 0 40; 2 5; 0.5 -1];
            want = [X(1:3,:)*s - c; inf; inf];
            testCase.verifyEqual(one.eval(X), want, 'AbsTol', 1e-12, 'T + cone((0,1))');
            two = EnvelopeMin({}, struct('T', T, 's', s, 'c', c, 'G', [0 1; 1 0]), {});
            X = [5 5; 0.3 0.2; -1 0.5; 0.5 -0.1];
            want = [X(1:2,:)*s - c; inf; inf];
            testCase.verifyEqual(two.eval(X), want, 'AbsTol', 1e-12, 'T + the first quadrant');
            none = EnvelopeMin({}, struct('T', T, 's', s, 'c', c, 'G', zeros(2,0)), {});
            testCase.verifyEqual(none.eval([0.2 5; 0.2 0.3]), [inf; 0.2 + 0.6 - c], 'AbsTol', 1e-12, ...
                'no generators: the triangle alone');
        end

        function aSINGLEConvexPieceIsOneAtomAndIsSentBack(testCase)
        % A convex piece is its own atom, and one atom is convex -- biconjQ's convex path answers
        % it, so this route refuses rather than build an envelope of nothing.
            f = QuaPol([0 0; 2 0; 0 2], [1 2 1; 2 3 1; 3 1 1], [0 0 0 0, 2, 0, 2, 0, 0, 0], ...
                       [1 0; 1 0; 1 0]);
            testCase.verifyError(@() biconjAtoms(f), 'biconjAtoms:oneAtom');
        end
    end
end
