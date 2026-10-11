classdef erratumBiconjTest < matlab.unittest.TestCase
% The BICONJUGATE of the three Karmarkar-Lucet erratum examples through biconjQ, which reaches them
% by the atom route (biconjAtoms), plus a differential check of that route against the legacy
% single-triangle envelope. Fixtures, oracles and the erratum's closed forms are erratumExamplesTest's
% statics; its header says how each oracle was validated.
%
% BUCKET: normal. biconjQ on Example 1 is ~45 s (one envelopeAssemble pair), so each example is
% computed ONCE per run and shared by its tests (`biconjCached`); the differential's 3-convex-edge
% triangle is ~100 s.
%
% QUARANTINES, by identifier only (umbrella CLAUDE.md section 8), measured 2026-10-10:
%   Example 2  bridgeFaceAlg:polynomialOverflow -- an algebraic bridge face's implicit polynomial,
%              content removed, needs 59 bits: the ANSWER is past what a double stores exactly.
%              (Until 2026-10-11 it stopped earlier at ratQ:overflow in three constructors whose
%              intermediates, not results, were large; fixed.) TODO "integer width".
%   (Example 3 was quarantined on PLQ:biconjQ:notImplemented until 2026-10-11, when the atom
%   route learned rays -- erratum Case E. It runs now.)
% Any other error, or a wrong value, fails as usual, and a test runs by itself once its refusal goes.

    methods (Test)

        function ex1BiconjugateMatchesTheErratumCells(testCase)
        % (H) fails here (P1 is strictly convex): cell 1 is bounded by a CIRCULAR arc and cell 2
        % is affine + 2 sqrt(21) |x - v|, neither of which a rational/polyhedral form can hold.
            h = erratumBiconjTest.biconjCached(testCase, 1);
            [X, want] = erratumExamplesTest.ex1BiconjCells();
            testCase.verifyEqual(h.eval(X), want, 'RelTol', 1e-9, 'AbsTol', 1e-9, ...
                'biconjQ disagrees with the erratum''s four cells on Example 1');
        end

        function ex1BiconjugateSatisfiesTheDefinition(testCase)
            h = erratumBiconjTest.biconjCached(testCase, 1);
            [f, pieces] = erratumExamplesTest.example1();
            erratumExamplesTest.checkEnvelopeDefinition(testCase, h, f, pieces, ...
                erratumExamplesTest.ex1DualSample());
        end

        function ex2BiconjugateMatchesTheErratumCells(testCase)
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

        function ex3BiconjugateIsX1Squared(testCase)
            h = erratumBiconjTest.biconjCached(testCase, 3);
            rng(20261011);
            X = [randn(200,1)*5, abs(randn(200,1))*5; 1 0; -3 0; 0 4];
            testCase.verifyEqual(h.eval(X), X(:,1).^2, 'RelTol', 1e-12, 'AbsTol', 1e-12);
            Y = [randn(50,1)*5, -abs(randn(50,1))*5 - 1e-3];
            testCase.verifyTrue(all(isinf(h.eval(Y)) & h.eval(Y) > 0), ...
                'f** must be +inf off the closed upper half-plane');
        end

        function theATOMRouteAgreesWithTheLegacyTriangleEnvelope(testCase)
        % Differential (umbrella CLAUDE.md section 6): on a single bounded triangle the legacy
        % route is convEnvCPLQ, verified against the definition to 8.9e-16 (SUPPORT_MATRIX.md
        % section 3), and shares no code with biconjAtoms. Indefinite triangles with 1, 2 and 3
        % convex edges -- vertex-against-edge cells, edge-against-edge cells (delta = 0, they
        % share a vertex) and a triple tangent plane. Before the rational corner landed
        % (cornerFaceAlg, 2026-10-10) these disagreed by up to 2.2.
            cases = { [0 0; 1 0; 0 1], [2 -3 2 0 0 0];      % 3 convex edges
                      [0 0; 2 0; 0 2], [1 2 -1 1 0 0];      % 1
                      [0 0; 3 0; 1 2], [2 1 -1 0 1 0];      % 2
                      [0 0; 2 0; 1 2], [-1 2 1 0 0 1] };    % 1
            rng(5);
            for k = 1:size(cases,1)
                V = cases{k,1};  c = cases{k,2};
                f = QuaPol(V, [1 2 1; 2 3 1; 3 1 1], [0 0 0 0 c], [1 0; 1 0; 1 0]);
                h = biconjQ(f);
                testCase.verifyClass(h, 'EnvelopeMin', 'an edge-convex indefinite piece takes the atom route');
                W = rand(300,3);  X = (W ./ sum(W,2)) * V;
                testCase.verifyEqual(h.eval(X), f.biconj().eval(X), 'AbsTol', 1e-10, ...
                    sprintf('triangle %d: atom route and legacy envelope differ', k));
            end
        end
    end

    methods (Static)
        function h = biconjCached(testCase, which)
        % biconjQ of erratum example `which`, computed once per MATLAB session, with the
        % quarantines of this file's header. A refusal is cached too, so it is not recomputed.
            persistent cache
            if isempty(cache), cache = cell(1,3); end
            if isempty(cache{which})
                switch which
                    case 1, f = erratumExamplesTest.example1();
                    case 2, f = erratumExamplesTest.example2();
                    case 3, f = erratumExamplesTest.example3();
                end
                try
                    cache{which} = struct('h', biconjQ(f), 'err', []);
                catch ME
                    cache{which} = struct('h', [], 'err', ME);
                end
            end
            ME = cache{which}.err;
            if isempty(ME), h = cache{which}.h; return, end
            if which == 2 && strcmp(ME.identifier, 'bridgeFaceAlg:polynomialOverflow')
                testCase.assumeFail(['QUARANTINED: Example 2 overflows the exact integer range ' ...
                    '(TODO: integer width -- an implicit bridge polynomial needs 59 bits). ' ME.message]);
            end
            rethrow(ME)
        end
    end
end
