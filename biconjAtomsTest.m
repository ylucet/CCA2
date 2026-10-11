classdef biconjAtomsTest < matlab.unittest.TestCase
% Unit tests for biconjAtoms' own contract. Its VALUES are checked in erratumBiconjTest (normal:
% erratum Example 1 and the differential against the legacy triangle envelope) and biconjQTest.
%
% BUCKET: fast (refusals only; nothing here builds a pair envelope).

    methods (Test)
        function theINPUTContractIsRefusedByName(testCase)
            testCase.verifyError(@() biconjAtoms(42), 'biconjAtoms:input');
            % an indefinite piece on an unbounded cone: rays are out of scope (erratum Case E)
            fin = QuaPol([0 0; 1 0; 1 1], [1 2 0; 1 3 0], [0 0 0 0, 1, 0, -1, 0, 0, 0], [1 0; 0 1]);
            testCase.verifyError(@() biconjAtoms(fin), 'biconjAtoms:unbounded');
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
