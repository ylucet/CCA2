classdef AlgConTest < matlab.unittest.TestCase
% Unit tests for AlgCon, the ROOT of the type lattice and the type `biconj` returns when a face of
% the convex envelope has no formula.
%
% BUCKET: fast (type queries and property validation; no symbolic call, no root-find beyond
% algBranch's own, which has its own suite).
%
% WHAT IS UNDER TEST HERE is the LATTICE and the PINNING, not the envelope: that a new level above
% Rat leaves every existing type's behaviour alone, that the axes stay separately queryable, and
% that no type can be made to lie about itself. The mathematics of an algebraic face is
% algBranchTest's subject; the producer that fills FP/FB is a separate item.

    methods (Test)

        % ---- the axis, and that adding a level above Rat changed nothing below it -------------

        function everyExistingTypeIsNowAlgAndThatIsTheHonestReading(testCase)
        % A rational function IS an algebraic one, so `isa(f,'Alg')` must be TRUE for every type
        % that existed before the level did. If it were false the level would be a separate
        % hierarchy rather than a generalisation, and `set.FP`'s one-site pinning could not work.
            q = QuaPol([3 1 3 0 0 0]);
            testCase.verifyTrue(isa(q,'Alg'),   'a polynomial is algebraic');
            testCase.verifyTrue(isa(q,'Rat'),   'and still rational');
            testCase.verifyTrue(isa(q,'Qua'),   'and still quadratic');
            testCase.verifyTrue(isa(q,'Conic'), 'and still on a conic subdivision');
            testCase.verifyTrue(isa(q,'AlgCon'),'and now under the new root');
            testCase.verifyTrue(isa(q,'RatCon'),'without leaving the old one');
            testCase.verifyEqual(q.kind(), 'QuaPol', 'kind() still names the concrete class');
        end

        function theTwoAxesStayINDEPENDENTLYQueryable(testCase)
        % The whole purpose of the markers: ask about ONE axis without enumerating combinations.
        % An AlgCon is at the TOP of the function axis, so it is Alg and NOT Rat -- the single
        % test that separates "has a formula" from "names a root".
            a = AlgCon();
            testCase.verifyTrue(isa(a,'Alg'));
            testCase.verifyFalse(isa(a,'Rat'), 'an AlgCon is NOT rational');
            testCase.verifyFalse(isa(a,'Qua'));
            testCase.verifyTrue(isa(a,'Conic'), 'the SUBDIVISION axis is untouched');
            testCase.verifyFalse(isa(a,'Par'),  'and is not narrowed to parabolic');
        end

        function theLATTICEIsAChainOnTheFunctionAxisAndTheOrderIsRight(testCase)
        % Qua < Rat < Alg, checked as DIRECT containment rather than asserted in a comment.
        % Reading it off the metaclass rather than from `isa` on one object is what makes it a
        % statement about the TYPES and not about a lucky instance -- and DIRECT superclasses are
        % what says the chain has no level silently inserted or skipped, which `isa` cannot.
            sup = @(c) {meta.class.fromName(c).SuperclassList.Name};
            testCase.verifyTrue(ismember('Alg', sup('Rat')),    'Rat < Alg, directly');
            testCase.verifyTrue(ismember('Rat', sup('Qua')),    'Qua < Rat, directly');
            testCase.verifyTrue(ismember('Alg', sup('AlgCon')), 'AlgCon < Alg');
            testCase.verifyTrue(ismember('Conic', sup('AlgCon')), 'AlgCon < Conic');
            testCase.verifyTrue(ismember('AlgCon', sup('RatCon')), 'RatCon < AlgCon');
            testCase.verifyTrue(ismember('Rat', sup('RatCon')),    'RatCon < Rat');
            % and the root really is the root: it has no data-carrying parent above it
            testCase.verifyEqual(sort(sup('AlgCon')), {'Alg','Conic'}, ...
                'the root sits on the two AXES and on nothing else');
        end

        % ---- the pinning: no type can be made to lie about itself -----------------------------

        function aRATIONALTypeCannotBeGivenAnImplicitFace(testCase)
        % The pin, and the reason the level is safe to add: every face of a Rat has a formula, so
        % FP must be empty on one. Enforced from the SINGLE site the property is declared, by
        % reading the object's own trait -- so it covers every rational type, including types
        % added later, without either of them naming the other.
            for f = {QuaPol([3 1 3 0 0 0]), QuaPar(), RatPol(), QuaCon()}
                obj = f{1};
                testCase.verifyError(@() setFP(obj), 'AlgCon:FPMustBeEmpty', class(obj));
                testCase.verifyError(@() setFB(obj), 'AlgCon:FBMustBeEmpty', class(obj));
            end
            function setFP(o), o.FP = {zeros(2,2,2)}; end %#ok<INUSD>
            function setFB(o), o.FB = 1; end              %#ok<INUSD>
        end

        function anALLZEROBranchListIsAcceptedEverywhere(testCase)
        % A face WITH a formula carries branch 0, so all-zero FB must be legal on a rational type:
        % it is what a producer writes while filling a mesh whose faces are not all algebraic yet,
        % and refusing it would force the producer to special-case its own bookkeeping.
            q = QuaPol([3 1 3 0 0 0]);
            q.FB = zeros(size(q.f,1),1);
            testCase.verifyEqual(q.FB, 0, 'all-zero is the "no branch" value, not a violation');
        end

        function anALGEBRAICFaceIsACCEPTEDOnAnAlgCon(testCase)
        % The other half of the pin: the type that IS allowed to carry one, does. Without this the
        % test above would pass on a class that simply refuses everything.
            a = AlgCon();
            P = zeros(3,3,3);  P(1,1,3) = 1;  P(3,1,1) = -1;  P(1,3,1) = -1;   % z^2 = x^2 + y^2
            a.FP = {P};
            a.FB = 2;
            testCase.verifyEqual(a.FB, 2);
            testCase.verifyEqual(size(a.FP{1}), [3 3 3]);
            % and the payload really is the face: the upper branch of that cone is +|x|
            testCase.verifyEqual(algBranch(a.FP{1}, [3 4], a.FB), 5, 'AbsTol', 1e-12);
        end

        function PROVENANCEIsNeverWhatAValueIsReadFrom(testCase)
        % `src`/`FR` exist for the classification table and for debugging. The contract this pins
        % is that they are OPTIONAL: a face's value comes from FP/FB alone, so an AlgCon with no
        % provenance at all still evaluates. If a value ever needed `src`, the 2026-08-23 "affine
        % cell naming a dual vertex" design would be back, and DECISIONS.md 2026-09-13 refutes it.
            a = AlgCon();
            P = zeros(3,3,3);  P(1,1,3) = 1;  P(3,1,1) = -1;  P(1,3,1) = -1;
            a.FP = {P};  a.FB = 2;
            testCase.verifyEmpty(a.src, 'no provenance recorded');
            testCase.verifyEqual(algBranch(a.FP{1}, [0 1], a.FB), 1, 'AbsTol', 1e-12, ...
                'and the value is still recoverable');
        end

        % ---- the constructor protocol survives the deeper chain -------------------------------

        function theROOTSOwnConstructorDoesNotClobberALeaf(testCase)
        % RatCon's CONSTRUCTOR PROTOCOL note records that MATLAB re-runs a shared base constructor
        % ONCE PER INHERITANCE PATH and the second run overwrites the first. Making the root
        % CONCRETE gives it a constructor of its own, so the protocol is re-checked here on the
        % real classes rather than only on the probe: a QuaPol built through the deeper chain must
        % still carry what its own leaf constructor wrote.
            q = QuaPol([3 1 3 0 0 0]);
            testCase.verifyEqual(q.f, [0 0 0 0 3 1 3 0 0 0], ...
                'the leaf''s coefficients survived every base invocation');
            testCase.verifyEqual(q.den, [0 0 1], 'and the pinned default is still filled in');
            testCase.verifyEmpty(QuaPol().f, 'the no-argument path still writes NOTHING');
            testCase.verifySize(QuaPol.empty(), [0 0], 'and .empty() still works');
        end

    end
end
