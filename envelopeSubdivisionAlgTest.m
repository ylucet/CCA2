classdef envelopeSubdivisionAlgTest < matlab.unittest.TestCase
% [normal] The first subdivision of a CURVED-TIE envelope. `envelopeSubdivision` refuses one -- a bridge face
% carries the value, and the boundary between it and a quadratic face is degree four, which a conic
% subdivision cannot hold as an edge. `envelopeSubdivisionAlg` builds the same arrangement with such edges
% allowed, so the GEOMETRY can be settled before the output type is grown to store it.
%
% `envelopeEval` is the oracle throughout: it was verified against the definition, and a subdivision is a
% re-arrangement of the same faces, so agreement is the whole statement -- on the VALUE where the
% envelope is finite and on the DOMAIN everywhere.
%
% BUCKET: normal. The assembly of a curved-tie fixture dominates, as it does in envelopeAssembleTest.

    methods (Test)
        function theCURVEDTieFixtureIsSUBDIVIDEDAndAgreesWithTheEvaluator(t)
            env = fixture(-2);
            [subdivision, info] = envelopeSubdivisionAlg(env);
            t.verifyGreaterThan(info.faces, 0, 'the subdivision must have faces');
            t.verifyGreaterThan(info.algRows, 0, ...
                'and ALGEBRAIC rows -- without one this fixture is the case envelopeSubdivision refuses');

            rng(29);
            X = [rand(400,1), rand(400,1); rand(200,1)*4-1.5, rand(200,1)*4-1.5];
            ref = envelopeEval(env, X);
            got = subdivisionAlgEval(subdivision, X);

            t.verifyEqual(isfinite(got), isfinite(ref), ...
                          'the subdivision and the evaluator must agree on the DOMAIN, point by point');
            both = isfinite(ref) & isfinite(got);
            t.verifyGreaterThan(sum(both), 200, 'and the value check must be exercised');
            t.verifyEqual(max(abs(got(both) - ref(both))), 0, 'AbsTol', 1e-9, ...
                          'and on the VALUE wherever both are finite');
        end

        function theCONVEXFixtureNeedsNoAlgebraicRowAtAll(t)
        % The control. A convex input has no curved tie, so no bridge face, so nothing to cut
        % algebraically -- and the same routine must still reproduce it exactly. Without this, the
        % test above would pass on code that answered every fixture with an algebraic row.
            env = fixture(0);
            [subdivision, info] = envelopeSubdivisionAlg(env);
            t.verifyEqual(info.algRows, 0, 'a convex fixture needs no algebraic edge');
            rng(31);
            X = [rand(200,1), rand(200,1)];
            ref = envelopeEval(env, X);
            got = subdivisionAlgEval(subdivision, X);
            both = isfinite(ref) & isfinite(got);
            t.verifyGreaterThan(sum(both), 100);
            t.verifyEqual(max(abs(got(both) - ref(both))), 0, 'AbsTol', 1e-9);
        end

        function aCurveSplitsOnlyWhereAFaceItBoundsCanApply(t)
        % The LOCAL arrangement (2026-10-31): a curve does not split a cell on which no face it
        % bounds can apply. Then it may CROSS a final cell, so the faces' rows must leave it out --
        % agreement with the evaluator on domain and value is what that has to keep. Measured: the
        % convex fixture goes from 58 cells to 32 with 20 skips, faces and values identical.
            env = fixture(0);
            [subdivision, info] = envelopeSubdivisionAlg(env);
            t.verifyGreaterThan(info.skipped, 0, 'some curve must be skipped, or this tests nothing');
            rng(37);
            X = [rand(600,1), rand(600,1); rand(200,1)*4-1.5, rand(200,1)*4-1.5];
            ref = envelopeEval(env, X);
            got = subdivisionAlgEval(subdivision, X);
            t.verifyEqual(isfinite(got), isfinite(ref), 'a skipped curve must not cut a face short');
            both = isfinite(ref) & isfinite(got);
            t.verifyGreaterThan(sum(both), 300);
            t.verifyEqual(max(abs(got(both) - ref(both))), 0, 'AbsTol', 1e-9);
        end

        function theMINIMUMOfTwoEnvelopesIsTheOverlayOfTheirSubdivisions(t)
        % SEVERAL PAIRS (2026-10-31): each envelope is subdivided alone and the two are overlaid
        % (minOverlay). The definition is pointwise -- the minimum of the two evaluators -- so
        % agreement on domain and value is again the whole statement. Two convex fixtures keep this
        % in the normal bucket; the tie rows of two irrational faces are minSubdivisionTest's.
            e1 = fixture(0);  e2 = fixture(1);
            [subdivision, info] = envelopeSubdivisionAlg({e1, e2});
            t.verifyEqual(info.pairCells > 0, [true true], 'each pair must be subdivided');
            t.verifyGreaterThan(info.faces, 0);
            rng(41);
            X = [rand(600,1), rand(600,1); rand(200,1)*4-1.5, rand(200,1)*4-1.5];
            ref = min(envelopeEval(e1, X), envelopeEval(e2, X));
            got = subdivisionAlgEval(subdivision, X);
            t.verifyEqual(isfinite(got), isfinite(ref), 'the overlay must agree on the DOMAIN');
            both = isfinite(ref) & isfinite(got);
            t.verifyGreaterThan(sum(both), 300);
            t.verifyEqual(max(abs(got(both) - ref(both))), 0, 'AbsTol', 1e-9, 'and on the VALUE');
        end

        function theARGUMENTSAreChecked(t)
            t.verifyError(@() envelopeSubdivisionAlg(), 'envelopeSubdivisionAlg:nargin');
            t.verifyError(@() envelopeSubdivisionAlg(struct('a', 1)), 'envelopeSubdivisionAlg:shape');
            t.verifyError(@() subdivisionAlgEval(struct('kind', {})), 'subdivisionAlgEval:nargin');
            t.verifyError(@() subdivisionAlgEval(struct('kind', {}), [1 2 3]), 'subdivisionAlgEval:points');
        end
    end
end

% ==================================================================================================
function env = fixture(c)
    fA = tri([0 0; 1 0; 1 1], 0);
    fB = tri([0 0; 1 1; 0 1], c);
    env = envelopePrune(envelopeAssemble(conjQ(fA), conjQ(fB)));
end

function f = tri(V, c)
    q = [0 0 0 0, 2, 0, 2, 0, 0, 0];  q(8) = q(8) - c;  q(9) = q(9) + c;
    f = QuaPol(V, [1 2 1; 2 3 1; 3 1 1], q, [1 0; 1 0; 1 0]);
end
