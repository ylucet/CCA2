classdef erratumBiconjEx2SubdivisionTest < matlab.unittest.TestCase
% Erratum Example 2's biconjugate SUBDIVIDED: three active pairs overlaid (envelopeSubdivisionAlg's
% minOverlay) and two contact triangles cut in after the overlay (cutTriangles), one of them with a
% degree-4 dual point (AlgCrv's table 3). TODO item 1 step d (2026-11-01).
%
% BUCKET: verylong -- the overlay visits up to 11567 pieces and ran for hours on the shared machine
% (DECISIONS 2026-11-01), so biconjQ declines it by default (envelopeSubdivisionAlg:overlaySize) and
% returns the EnvelopeMin; this test lifts the budget. The EnvelopeMin is the reference: it is the
% minimum over the same pairs and triangles, evaluated face list by face list.

    properties
        f
        pieces
        g       % the subdivision, as biconjQ returns it with the budget lifted
        ref     % the EnvelopeMin
    end

    methods (TestClassSetup)
        function build(t)
            [t.f, t.pieces] = erratumExamplesTest.example2();
            t.ref = biconjQ(t.f);
            t.g = biconjQ(t.f, struct('overlayBudget', Inf));
        end
    end

    methods (Test)
        function itIsSubdivided(t)
            t.verifyClass(t.ref, 'EnvelopeMin', 'the default budget declines Example 2''s overlay');
            t.verifyClass(t.g, 'AlgCrv');
            t.verifyEqual(nnz(t.g.FT), 1, 'the degree-4 triangle is a face of its own');
        end

        function itEqualsTheEnvelopeMinOn4000Points(t)
            rng(20261101);
            V = t.f.VN ./ t.f.VD;  lo = min(V);  hi = max(V);
            X = lo + rand(4000, 2) .* (hi - lo);
            want = t.ref.eval(X);  got = t.g.eval(X);
            t.verifyEqual(isinf(got), isinf(want), 'a hole, or a value where the minimum is +inf');
            fin = isfinite(want);
            t.verifyGreaterThan(nnz(fin), 1000);
            t.verifyEqual(got(fin), want(fin), 'RelTol', 1e-9, 'AbsTol', 1e-7);
        end

        function itSatisfiesTheDefinition(t)
            erratumExamplesTest.checkEnvelopeDefinition(t, t.g, t.f, t.pieces, ...
                erratumExamplesTest.ex2DualSample());
        end
    end
end
