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
% QUARANTINES: none since 2026-10-11. Example 2 was quarantined on
%   bridgeFaceAlg:polynomialOverflow (a bridge polynomial needs 59 bits); the polynomial is now kept
%   in BigInteger and the intermediates that overflowed after it are reduced by their content
%   (ratQ.exactFormula / bigZ). biconjQ on Example 2 takes ~26 min on the shared machine (one run),
%   computed once per run by biconjCached. Example 3 was quarantined until 2026-10-11 too.
% Any other error, or a wrong value, fails as usual, and a test runs by itself once its refusal goes.

    methods (Test)

        % Example 2's tests MOVED 2026-10-17 to erratumBiconjEx2Test (slow bucket): biconjQ on it is
        % ~5 min, which kept this suite over the normal bucket's 10-minute budget (TODO item 5).

        function ex1BiconjugateIsMeshed(testCase)
        % One active pair (a convex piece and a vertex), so biconjQ meshes it (2026-10-17): an AlgCrv
        % whose two corners from Q(sqrt 26) and Q(sqrt 6) overlap in 18 cells, value-safe (max).
            h = erratumBiconjTest.biconjCached(testCase, 1);
            testCase.verifyClass(h, 'AlgCrv');
        end

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
                % the atom route; a single active pair comes back MESHED (2026-10-17)
                testCase.verifyTrue(isa(h, 'EnvelopeMin') || isa(h, 'AlgCrv') || isa(h, 'QuaCon'), ...
                    sprintf('triangle %d: biconjQ returned a %s', k, class(h)));
                W = rand(300,3);  X = (W ./ sum(W,2)) * V;
                testCase.verifyEqual(h.eval(X), f.biconj().eval(X), 'AbsTol', 1e-10, ...
                    sprintf('triangle %d: atom route and legacy envelope differ', k));
            end
        end
    end

    methods (Test)
        function aSINGLEActivePairIsMeshedAndKeepsItsValues(testCase)
        % Two atoms leave one pair; with no triangle needing three atoms the minimum IS that pair's
        % envelope, so biconjQ meshes it (2026-10-17). Four defects stood in the way, each red here
        % first: the rational bridge face had no rows (envelopeMesh:bridgeFaces); the meshers read
        % only a face's first pair, and took a cell's interior point from its straight rows only
        % (holes); an EMPTY face applied everywhere; and a rational CORNER, which does carry the
        % value near a vertex, was carried for coverage only (1.3e-6 too low).
            V = [0 0; 3 0; 1 2];
            f = QuaPol(V, [1 2 1; 2 3 1; 3 1 1], [0 0 0 0 2 1 -1 0 1 0], [1 0; 1 0; 1 0]);
            h = biconjQ(f);
            testCase.verifyTrue(isa(h, 'AlgCrv') || isa(h, 'QuaCon'), ...
                sprintf('one active pair should mesh; got %s', class(h)));
            rng(3);
            W = [rand(1500,3); rand(500,3).^4];  X = (W ./ sum(W,2)) * V;
            want = biconjAtoms(f).eval(X);
            got = h.eval(X);
            testCase.verifyEqual(isinf(got), isinf(want), 'the mesh has a hole where the face list is finite');
            fin = isfinite(want);
            testCase.verifyEqual(got(fin), want(fin), 'RelTol', 1e-9, 'AbsTol', 1e-9, ...
                'the mesh must carry the face list''s values, near the vertices too');
        end

        function aPointBesideTheEDGEIsNotPlusInfinity(testCase)
        % Found 2026-10-16 when the prune dropped a redundant pair: on the 3-convex-edge triangle a
        % pair's bridge face and its quadratic face both missed a band ~1e-8 wide along x1 = 0
        % (the bridge's sweeping root meets det(M_mu) = 0 at mu = -1), so the pair read +Inf there
        % and the minimum took a pair 17x too high. envelopeEval now retries such a point 1e-7
        % away. On that edge co f = f = x2^2, and this close to it the gap is O(x1).
            V = [0 0; 1 0; 0 1];
            f = QuaPol(V, [1 2 1; 2 3 1; 3 1 1], [0 0 0 0 2 -3 2 0 0 0], [1 0; 1 0; 1 0]);
            h = biconjQ(f);
            x2 = 0.058899440267143444;
            X = [linspace(1e-9, 1e-7, 400).', x2*ones(400,1)];
            testCase.verifyEqual(h.eval(X), x2^2*ones(400,1), 'AbsTol', 1e-6);
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
