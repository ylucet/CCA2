classdef erratumBiconjEx2Test < matlab.unittest.TestCase
% Erratum Example 2's biconjugate through biconjQ -- the tests erratumBiconjTest held until
% 2026-10-17. Three indefinite pieces, four edge atoms, three active pairs and two triangle faces:
% an EnvelopeMin (a minimum over several pairs is not meshed; TODO item 1).
%
% BUCKET: slow -- biconjQ on Example 2 is ~5 min on the shared machine (one run), computed once per
% process by erratumBiconjTest.biconjCached and shared by the tests below.

    methods (Test)

        function ex2BiconjugateMatchesTheErratumCells(testCase)
        % Cell 6's test point was (32,-5.5), which lies in cell 4, not cell 6 (fixed 2026-10-11 in
        % ex2BiconjCells); with a point inside cell 6 the erratum's formula and biconjQ agree.
            h = erratumBiconjTest.biconjCached(testCase, 2);
            [X, want, names] = erratumExamplesTest.ex2BiconjCells();
            got = h.eval(X);
            for i = 1:size(X,1)
                testCase.verifyEqual(got(i), want(i), 'RelTol', 1e-9, 'AbsTol', 1e-7, ...
                    sprintf('Example 2, %s at (%g,%g)', names{i}, X(i,1), X(i,2)));
            end
        end

        function ex2BiconjugateIsAffineOnTheMiddleCell(testCase)
        % Cell 4 = conv{G1,G2,G3} carries s_bar'x - f*(s_bar), AFFINE; s_bar has no closed form,
        % so assert the property: the value at a convex combination is the combination of values.
            h = erratumBiconjTest.biconjCached(testCase, 2);
            P = [36.5 0.9; 37.5 1.4; 35.5 0.2];       % inside conv{G1,G2,G3}
            w = [0.2 0.5 0.3];
            testCase.verifyEqual(h.eval(w*P), w*h.eval(P), 'RelTol', 1e-9, 'AbsTol', 1e-7);
        end

        function ex2BiconjugateSatisfiesTheDefinition(testCase)
            h = erratumBiconjTest.biconjCached(testCase, 2);
            [f, pieces] = erratumExamplesTest.example2();
            erratumExamplesTest.checkEnvelopeDefinition(testCase, h, f, pieces, ...
                erratumExamplesTest.ex2DualSample());
        end

        function exampleTwoAssemblesOnlyItsACTIVEPairs(testCase)
        % erratum Step 1: only atom pairs whose conjugate forms tie for f* along an arc, and only
        % triangles whose three forms need three atoms. Example 2 has four edge atoms: 6 pairs and
        % 60 triangles before (2026-10-16), of which exactly pairs (1,3), (2,3), (3,4) ever attained
        % the minimum on a 4000-point sample. Values are pinned by the tests above.
            h = erratumBiconjTest.biconjCached(testCase, 2);
            testCase.verifyClass(h, 'EnvelopeMin');
            testCase.verifyEqual(numel(h.pairs), 3, 'active pairs of Example 2');
            testCase.verifyEqual(numel(h.tris), 2, 'active triangles of Example 2');
        end
    end
end
