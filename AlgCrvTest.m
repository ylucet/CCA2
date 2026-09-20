classdef AlgCrvTest < matlab.unittest.TestCase
% [fast] The (Alg, Crv) cell of the lattice: an algebraic function on an algebraic subdivision.
% Two things are pinned here -- the PIN that keeps the new edge table off every type below this one,
% and the EVALUATOR that reads such an edge -- on a fixture small enough to check by hand.
%
% The fixture is the cone z^2 = x^2 + y^2 on its upper branch, which is |x|, cut against the constant
% face q = 1: the algebraic edge is the unit circle, and inside it the constant is the larger. That
% is a face whose value is irrational, bounded by an edge that is not a conic -- which is the whole
% reason this type exists.

    methods (Test)
        function theEDGETableIsPINNEDOffEveryTypeThatIsNotACrv(t)
        % The same rule `den`, `Ec` and `FP` follow, and the reason the property lives on AlgCon: a
        % conic subdivision has no algebraic edge, and saying so by TRAIT covers types that do not
        % exist yet.
            a = AlgCon();
            t.verifyEmpty(a.EcA, 'the default is empty');
            t.verifyError(@() setEcA(a), 'AlgCon:EcAMustBeEmpty');
            q = QuaCon(zeros(0,3), zeros(0,3), zeros(0,6), zeros(0,10), zeros(0,1), zeros(0,2), {});
            t.verifyError(@() setEcA(q), 'AlgCon:EcAMustBeEmpty');
        end

        function anALGEBRAICEdgeIsACCEPTEDOnAnAlgCrv(t)
        % The other half of the pin: the type that IS allowed to carry one, does -- without which
        % the test above would pass on a class that simply refused everything.
            g = fixture();
            t.verifyTrue(isa(g, 'Crv'));
            t.verifyTrue(isa(g, 'AlgCon'));
            t.verifyEqual(numel(g.EcA), 1);
            t.verifyEqual(g.EcA(1).FB, 2);
        end

        function theEVALUATORReadsBothTheEdgeAndTheIMPLICITValue(t)
        % Inside the unit circle the constant 1 is above |x|, so the face owns those points and its
        % value is the CONSTANT; outside, the face owns nothing and the answer is +infinity, which is
        % what a hole in a mesh means everywhere in this toolbox.
            g = fixture();
            inside  = [0 0; 0.5 0; 0 0.5; -0.25 0.25];
            outside = [2 0; 0 -3; 1.5 1.5];
            v = g.eval(inside);
            t.verifyEqual(v, ones(size(inside,1),1), 'AbsTol', 1e-12, ...
                          'the face''s value is the quadratic one, which is 1 here');
            t.verifyTrue(all(isinf(g.eval(outside))), ...
                         'and outside its edge no face owns the point');
        end

        function aFaceWhoseValueIsIMPLICITIsEvaluatedFromFPAndFB(t)
        % The other face kind: the same cone, but carrying the BRANCH as its value rather than a
        % formula. On the upper branch the value is |x|, which is what makes the type an Alg.
            [FP, num, den] = cone();
            FC = {[1 1 0 1]};                       % one algebraic row, table 1
            EcA = struct('FP', FP, 'FB', 2, 'num', num, 'den', den);
            g = AlgCrv(zeros(0,3), zeros(0,3), zeros(0,6), zeros(1,10), 1, zeros(0,2), FC, ...
                       {FP}, 2, EcA);
            X = [0.5 0; 0 0.5; 0.3 0.4];
            t.verifyEqual(g.eval(X), hypot(X(:,1), X(:,2)), 'AbsTol', 1e-12, ...
                          'the value is the upper branch of the cone, |x|');
        end

        function theARGUMENTSAreChecked(t)
            g = fixture();
            t.verifyError(@() g.eval(), 'AlgCrv:nargin');
            t.verifyError(@() g.eval([1 2 3]), 'AlgCrv:points');
            t.verifyError(@() assembleAlgCrv(), 'assembleAlgCrv:nargin');
            t.verifyError(@() assembleAlgCrv(5), 'assembleAlgCrv:shape');
        end
    end
end

% ==================================================================================================
function g = fixture()
% One face: the constant 1, owning the inside of the unit circle because that is where it beats the
% cone's upper branch.
    [FP, num, den] = cone();
    EcA = struct('FP', FP, 'FB', 2, 'num', num, 'den', den);
    FC = {[1 1 0 1]};                               % side +1: q - value_FB >= 0, i.e. inside
    g = AlgCrv(zeros(0,3), zeros(0,3), zeros(0,6), num, den, zeros(0,2), FC, {[]}, 0, EcA);
end

function [FP, num, den] = cone()
    FP = zeros(3,3,3);
    FP(1,1,3) =  1;  FP(3,1,1) = -1;  FP(1,3,1) = -1;   % z^2 = x^2 + y^2
    num = [0 0 0 0, 0, 0, 0, 0, 0, 1];                  % the constant 1, in the weighted basis
    den = 1;
end

function setEcA(obj)
    obj.EcA = struct('FP', {1}, 'FB', {1}, 'num', {1}, 'den', {1}); %#ok<NASGU>
end
