classdef envelopeMeshCornersTest < matlab.unittest.TestCase
% [fast] A CORNER face in the mesh carries its REGION, not its value.
%
% Measured on the c = -2 fixture (DECISIONS.md 2026-09-19): wherever a corner face owns a point its
% value never exceeds the best quadratic face owning the same point, and the largest excess is
% exactly 0.0e0. So the mesh carries a corner for what its wedge DECIDES -- by `envelopeEval`'s rule
% the envelope is the max over faces and a PAIR with no face owning x is +infinity there, so a pair
% whose only candidate is a corner makes the domain end where that wedge does.
%
% The envelope here is built by hand rather than assembled, because a corner face never occurs
% without a BRIDGE face in a real assembly and a bridge's region is still a predicate -- so this is
% the only way to exercise the path until that is solved. `envelopeEval` is the oracle: it was
% verified against the definition, and the mesh is a re-arrangement of the same faces.

    methods (Test)
        function theDOMAINEndsWhereTheCornersWEDGEDoes(t)
            [env, corner] = fixture();
            [g, info] = envelopeMesh(env);
            t.verifyEqual(info.corners, 1, 'the corner must have been carried, not refused');
            t.verifyGreaterThan(info.faces, 0);

            rng(23);
            X = [rand(600,1)*6-3, rand(600,1)*6-3];
            ref = envelopeEval(env, X);
            got = g.eval(X);
            t.verifyEqual(isfinite(got), isfinite(ref), ...
                          'the mesh and the evaluator must agree on the DOMAIN, point by point');

            both = isfinite(ref) & isfinite(got);
            t.verifyGreaterThan(sum(both), 0, 'the fixture must have a nonempty domain');
            t.verifyLessThan(max(abs(got(both) - ref(both))), 1e-9, ...
                             'and on the VALUE wherever both are finite');

            % The domain really is the wedge: every point the mesh answers finitely satisfies the
            % corner's own KKT predicate, and points outside it are +infinity.
            inWedge = false(size(X,1),1);
            for i = 1:size(X,1)
                mu = corner.M \ (corner.p - X(i,:).');
                inWedge(i) = all(mu >= -1e-9);
            end
            t.verifyEqual(isfinite(got), inWedge, ...
                          'the finite set must BE the wedge, not merely agree with the oracle');
        end

        function aFaceWhoseRegionHasNORowsIsStillRefusedByName(t)
        % A corner has rows; a bridge does not, and admitting one would mesh a region nothing
        % describes. The refusal names the kind so a caller can tell the two apart.
            env = fixture();
            env.bridge{end+1} = struct('kind', 'algebraic', 'pairIndex', 1);
            t.verifyError(@() envelopeMesh(env), 'envelopeMesh:bridgeFaces');
        end
    end
end

% ==================================================================================================
function [env, corner] = fixture()
% Two pairs. Pair 1 is a quadratic face on the whole plane, so it is covered everywhere and supplies
% every value. Pair 2's ONLY candidate is a corner, so the envelope is +infinity exactly outside
% that corner's wedge -- which is the rule under test, isolated.
    n = [0 0 0 0, 2, 0, 2, 0, 0, 0];
    P = [1 0 1; -1 0 1; 0 1 1; 0 -1 1];
    tie = [2 0 2 0 0 -3];
    faces = cornerFaceAlg(n, 1, P, tie, +1);
    corner = faces(1);
    corner.pairIndex = 2;

    % x^2 + y^2 over the whole plane, with a vacuous row so the face has a region at all
    quad = struct('num', [0 0 0 0, 1, 0, 1, 0, 0, 0], 'den', 1, ...
                  'con', [0 0 0 0 0 1, 1, 0], 'pair', 1);
    env = struct('quad', quad, 'bridge', {{corner}});
end
