classdef RatParTest < matlab.unittest.TestCase
    % Tests for RatPar, the common parent of every type CCA2's operators return.
    %
    % These pin the CONTRACT a downstream caller depends on -- "conj returns a RatPar, and kind()
    % tells me what it actually is" -- plus the two structural decisions that contract rests on:
    % the four types are SIBLINGS (not a QuaPoly < RatPol chain), and kind() is derived from the
    % real class rather than stored. See RatPar.m and RETURN_TYPE.md.

    methods (Test)
        function everyReturnedTypeIsARatPar(testCase)
            % The whole point of the class: one static return type for callers.
            for obj = RatParTest.oneOfEach()
                testCase.verifyTrue(isa(obj{1}, 'RatPar'), ...
                    sprintf('%s must be a RatPar', class(obj{1})));
            end
        end

        function kindReportsTheConcreteType(testCase)
            samples = RatParTest.oneOfEach();
            expected = {'QuaPoly', 'RatPol', 'QuaPar', 'QuaParCPLQ'};
            for i = 1:numel(samples)
                testCase.verifyEqual(samples{i}.kind(), expected{i});
            end
        end

        function kindCannotDriftFromTheActualClass(testCase)
            % kind() is a METHOD over class(), not a stored flag, precisely so it can never
            % disagree with what the object is -- including after a copy, which is how every
            % operator returns its result (all four types are value classes).
            for obj = RatParTest.oneOfEach()
                copyOfIt = obj{1};                      % value-class copy
                testCase.verifyEqual(copyOfIt.kind(), class(copyOfIt));
                testCase.verifyEqual(copyOfIt.kind(), obj{1}.kind());
            end
        end

        function theFourTypesAreSiblingsNotAChain(testCase)
            % Deliberate: QuaPoly is mathematically a special case of RatPol, but making it a
            % SUBCLASS would make isa(aQuaPoly,'RatPol') true, and code that (correctly) reads
            % p.den after that test would then fail on a QuaPoly, which has no denominator --
            % e.g. conjPieceCPLQ's rational-piece guard. See RatPar.m's header.
            V = [0 0; 1 0; 0 1]; E = [1 2 1; 2 3 1; 3 1 1]; F = [1 0; 1 0; 1 0];
            qpoly = QuaPoly(V, E, [0 1 0 0 0 0], F);
            rpol  = RatPol(V, E, [0 1 0 0 0 0], F);
            qpar  = QuaPar(V, E, zeros(3,6), [0 1 0 0 0 0], F);

            testCase.verifyFalse(isa(qpoly, 'RatPol'));
            testCase.verifyFalse(isa(qpoly, 'QuaPar'));
            testCase.verifyFalse(isa(qpar,  'RatPol'));
            testCase.verifyFalse(isa(rpol,  'QuaPoly'));

            % ...and the guard that motivated the decision still behaves: only the RatPol has a
            % denominator to read at all.
            testCase.verifyTrue(isprop(rpol, 'den'));
            testCase.verifyFalse(isprop(qpoly, 'den'));
            testCase.verifyFalse(isprop(qpar, 'den'));
            testCase.verifyTrue(isprop(qpar, 'Ec'));
            testCase.verifyFalse(isprop(qpoly, 'Ec'));
        end

        function eachSubclassKeepsOnlyItsOwnExtraProperty(testCase)
            % The nine mesh properties moved to RatPar; the children add exactly one each (RatPol
            % `den`, QuaPar `Ec`) or none (QuaPoly). This pins that the move was a pure relocation:
            % every inherited property is still present and readable on every child.
            common = {'nv','ne','nf','V','E','f','F','P','dom'};
            for obj = RatParTest.oneOfEach()
                for c = common
                    testCase.verifyTrue(isprop(obj{1}, c{1}), ...
                        sprintf('%s lost inherited property %s', class(obj{1}), c{1}));
                end
            end
        end

        function ratParItselfCannotBeInstantiated(testCase)
            % RatPar is abstract: a bare RatPar would be a function of no particular type.
            testCase.verifyError(@() RatPar(), 'MATLAB:class:abstract');
        end

        function isMeshedSeparatesTheSymbolicForm(testCase)
            % QuaParCPLQ is the one child with no reconstructed V/E/Ec/F/P mesh -- it holds cPLQ's
            % own symbolic form. Its inherited mesh properties are empty for a real reason (the
            % geometric reconstruction is still open, see DESIGN.md II.5.1), and isMeshed() is how
            % a caller tells the two apart without dispatching on class.
            V = [0 0; 1 0; 0 1]; E = [1 2 1; 2 3 1; 3 1 1]; F = [1 0; 1 0; 1 0];
            testCase.verifyTrue(QuaPoly(V, E, [0 1 0 0 0 0], F).isMeshed());
            testCase.verifyTrue(RatPol(V, E, [0 1 0 0 0 0], F).isMeshed());
            testCase.verifyTrue(QuaPar(V, E, zeros(3,6), [0 1 0 0 0 0], F).isMeshed());
            testCase.verifyFalse(QuaParCPLQ().isMeshed());
        end

        function conjReturnsARatParAcrossEverySupportedInputShape(testCase)
            % End to end, the actual promise: whatever supported shape goes in, a RatPar comes
            % out, and kind() names it. The three shapes are conjCPLQ's own Cases A/B/C.
            fullQuad = QuaPoly([1 0 1 0 0 0]);                       % Case A: full-domain PD
            gA = fullQuad.conj();
            testCase.verifyTrue(isa(gA, 'RatPar'));
            testCase.verifyEqual(gA.kind(), 'QuaPoly');

            V = [0 0; 1 0; 0 1]; E = [1 2 1; 2 3 1; 3 1 1]; F = [1 0; 1 0; 1 0];
            gB = QuaPoly(V, E, [0 1 0 0 0 0], F).conj();             % Case B: bounded triangle
            testCase.verifyTrue(isa(gB, 'RatPar'));
            testCase.verifyEqual(gB.kind(), 'QuaPar');
        end
    end

    methods (Static)
        function s = oneOfEach()
            % One instance of each concrete RatPar subclass, as a cell row.
            V = [0 0; 1 0; 0 1]; E = [1 2 1; 2 3 1; 3 1 1]; F = [1 0; 1 0; 1 0];
            s = { QuaPoly(V, E, [0 1 0 0 0 0], F), ...
                  RatPol(V, E, [0 1 0 0 0 0], F), ...
                  QuaPar(V, E, zeros(3,6), [0 1 0 0 0 0], F), ...
                  QuaParCPLQ() };
        end
    end
end
