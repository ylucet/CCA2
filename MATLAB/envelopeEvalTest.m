classdef envelopeEvalTest < matlab.unittest.TestCase
% Unit tests for envelopeEval's own rules, on hand-built face lists so that each rule is isolated
% from the assembly that normally feeds it (envelopeAssembleTest covers the assembled envelopes).
%
% BUCKET: fast (no conjugate is computed; a face list of one or two quadratic faces).

    methods (Test)

        function aPOINTSValueDoesNotDependOnTheOtherPointsInTheCall(testCase)
        % Found 2026-10-10, the same defect QuaCon.eval had: inRegion scaled its tolerance by the
        % largest |x| in the WHOLE call, so a point just OUTSIDE a face's region was accepted --
        % and got that face's value instead of +inf -- once a far point shared the call. One face,
        % region x1 >= 0, value x1 + x2; the query is 1e-3 outside it.
            env = envelopeEvalTest.halfPlaneEnv();
            P = [-1e-3 0.5];
            alone = envelopeEval(env, P);
            batch = envelopeEval(env, [P; 1e6 0]);
            testCase.verifyEqual(alone, inf, 'outside the only face: +inf');
            testCase.verifyEqual(batch(1), alone, 'and the same with a far point in the call');
            testCase.verifyEqual(batch(2), 1e6, 'the far point itself is inside', 'RelTol', 1e-12);
        end

        function aPOINTInsideTheRegionGetsTheFacesValue(testCase)
            env = envelopeEvalTest.halfPlaneEnv();
            X = [0 0; 1 2; 3 -1; 1e-3 0.5];
            testCase.verifyEqual(envelopeEval(env, X), X(:,1) + X(:,2), 'AbsTol', 1e-12);
        end
    end

    methods (Static)
        function env = halfPlaneEnv()
            % value x1 + x2 in the weighted basis; region row: the line x1 = 0 as a conic
            % [0 0 0 1 0 0], side +1 (x1 >= 0), no pair-line branch
            env = struct('quad', struct('num', [0 0 0 0 0 0 0 1 1 0], 'den', 1, ...
                                        'con', [0 0 0 1 0 0 1 0], 'pair', 1), 'bridge', {{}});
        end
    end
end
