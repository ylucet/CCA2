classdef envelopeMeshAlgTest < matlab.unittest.TestCase
% [normal] The first mesh of a CURVED-TIE envelope. `envelopeMesh` refuses one -- a bridge face
% carries the value, and the boundary between it and a quadratic face is degree four, which a conic
% subdivision cannot hold as an edge. `envelopeMeshAlg` builds the same arrangement with such edges
% allowed, so the GEOMETRY can be settled before the output type is grown to store it.
%
% `envelopeEval` is the oracle throughout: it was verified against the definition, and a mesh is a
% re-arrangement of the same faces, so agreement is the whole statement -- on the VALUE where the
% envelope is finite and on the DOMAIN everywhere.
%
% BUCKET: normal. The assembly of a curved-tie fixture dominates, as it does in envelopeAssembleTest.

    methods (Test)
        function theCURVEDTieFixtureMESHESAndAgreesWithTheEvaluator(t)
            env = fixture(-2);
            [mesh, info] = envelopeMeshAlg(env);
            t.verifyGreaterThan(info.faces, 0, 'the mesh must have faces');
            t.verifyGreaterThan(info.algRows, 0, ...
                'and ALGEBRAIC rows -- without one this fixture is the case envelopeMesh refuses');

            rng(29);
            X = [rand(400,1), rand(400,1); rand(200,1)*4-1.5, rand(200,1)*4-1.5];
            ref = envelopeEval(env, X);
            got = meshAlgEval(mesh, X);

            t.verifyEqual(isfinite(got), isfinite(ref), ...
                          'the mesh and the evaluator must agree on the DOMAIN, point by point');
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
            [mesh, info] = envelopeMeshAlg(env);
            t.verifyEqual(info.algRows, 0, 'a convex fixture needs no algebraic edge');
            rng(31);
            X = [rand(200,1), rand(200,1)];
            ref = envelopeEval(env, X);
            got = meshAlgEval(mesh, X);
            both = isfinite(ref) & isfinite(got);
            t.verifyGreaterThan(sum(both), 100);
            t.verifyEqual(max(abs(got(both) - ref(both))), 0, 'AbsTol', 1e-9);
        end

        function theARGUMENTSAreChecked(t)
            t.verifyError(@() envelopeMeshAlg(), 'envelopeMeshAlg:nargin');
            t.verifyError(@() envelopeMeshAlg(struct('a', 1)), 'envelopeMeshAlg:shape');
            t.verifyError(@() meshAlgEval(struct('kind', {})), 'meshAlgEval:nargin');
            t.verifyError(@() meshAlgEval(struct('kind', {}), [1 2 3]), 'meshAlgEval:points');
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
