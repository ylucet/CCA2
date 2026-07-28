classdef RatParTest < matlab.unittest.TestCase
    % Tests for RatPar and the two-axis type lattice.
    %
    % These pin the contract a downstream caller depends on -- "conj returns a RatPar, and kind()
    % tells me what it actually is" -- plus the structural decisions it rests on: the four named
    % types are the cells of a Rat/Qua x Par/Pol grid, the inheritance is that grid's own partial
    % order (a genuine diamond), all data lives on RatPar, and each type's pinned values are
    % enforced by trait-keyed validators. See RatPar.m and RETURN_TYPE.md.

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
            expected = {'QuaPol', 'RatPol', 'QuaPar', 'QuaParCPLQ'};
            for i = 1:numel(samples)
                testCase.verifyEqual(samples{i}.kind(), expected{i});
            end
        end

        function kindCannotDriftFromTheActualClass(testCase)
            % kind() is a METHOD over class(), not a stored flag, precisely so it can never
            % disagree with what the object is -- including after a copy, which is how every
            % operator returns its result (all types here are value classes).
            for obj = RatParTest.oneOfEach()
                copyOfIt = obj{1};                      % value-class copy
                testCase.verifyEqual(copyOfIt.kind(), class(copyOfIt));
                testCase.verifyEqual(copyOfIt.kind(), obj{1}.kind());
            end
        end

        function theLatticeMatchesTheMathematics(testCase)
            % The full isa table, asserted entry by entry. This IS the design:
            %   Qua < Rat   a polynomial is a rational with unit denominator
            %   Pol < Par   a polyhedral subdivision is a parabolic one with no curvature
            %   QuaPol < RatPol & QuaPar   -- a genuine diamond, because the four named types are
            %                                  the cells of a 2x2 product lattice
            % Note the two subsumptions propagate for free: RatPol answers true to Par (via Pol),
            % and QuaPar answers true to Rat (via Qua).
            types = {'Rat','Qua','Par','Pol','RatPar','RatPol','QuaPar','QuaPol'};
            expected = [ ...
                1 1 1 1 1 1 1 1     % QuaPol    -- everything: bottom of the lattice
                1 0 1 1 1 1 0 0     % RatPol     -- rational, polyhedral
                1 1 1 0 1 0 1 0     % QuaPar     -- quadratic, parabolic
                1 1 1 0 1 0 0 0 ];  % QuaParCPLQ -- same math as QuaPar, but no mesh (see below)
            samples = RatParTest.oneOfEach();
            for i = 1:numel(samples)
                for j = 1:numel(types)
                    testCase.verifyEqual(double(isa(samples{i}, types{j})), expected(i,j), ...
                        sprintf('isa(%s, ''%s'') should be %d', ...
                            class(samples{i}), types{j}, expected(i,j)));
                end
            end
        end

        function traitsAnswerEachAxisIndependently(testCase)
            % What the markers buy over the combination names: asking about ONE axis without
            % enumerating combinations. `isa(f,'Qua')` replaces `isa(f,'QuaPol')||isa(f,'QuaPar')`,
            % and keeps working when a type is added later.
            V = [0 0; 1 0; 0 1]; E = [1 2 1; 2 3 1; 3 1 1]; F = [1 0; 1 0; 1 0];
            quaPol = QuaPol(V, E, [0 1 0 0 0 0], F);
            ratPol  = RatPol(V, E, [0 1 0 0 0 0], F);
            quaPar  = QuaPar(V, E, zeros(3,6), [0 1 0 0 0 0], F);

            % function axis
            testCase.verifyTrue(isa(quaPol,'Qua'));   testCase.verifyTrue(isa(quaPar,'Qua'));
            testCase.verifyFalse(isa(ratPol,'Qua'));
            % subdivision axis
            testCase.verifyTrue(isa(quaPol,'Pol'));   testCase.verifyTrue(isa(ratPol,'Pol'));
            testCase.verifyFalse(isa(quaPar,'Pol'));
            % the axes are independent: knowing one tells you nothing about the other
            testCase.verifyNotEqual(isa(ratPol,'Qua'), isa(ratPol,'Pol'));
            testCase.verifyNotEqual(isa(quaPar,'Qua'), isa(quaPar,'Pol'));
        end

        function pinnedValuesAreEnforcedByTraitKeyedValidators(testCase)
            % den and Ec are declared ONLY on RatPar -- they have to be, since a property defined
            % in two superclasses is fatal and unresolvable in MATLAB, which would make
            % QuaPol < RatPol & QuaPar impossible. So the pinning cannot be expressed by
            % overriding; it is enforced instead by set validators that read the object's TRAITS.
            % One definition site, whole lattice covered, and no type can lie about itself.
            V = [0 0; 1 0; 0 1]; E = [1 2 1; 2 3 1; 3 1 1]; F = [1 0; 1 0; 1 0];
            quaPol = QuaPol(V, E, [0 1 0 0 0 0], F);
            ratPol  = RatPol(V, E, [0 1 0 0 0 0], F);
            quaPar  = QuaPar(V, E, zeros(3,6), [0 1 0 0 0 0], F);

            % a Rat may carry a genuine denominator; a Qua may not
            ratPol.den = [1 2 3];                       % must not throw
            testCase.verifyError(@() setfield(quaPol, 'den', [1 2 3]), ...
                'RatPar:denMustBeTrivial'); %#ok<SFLD>
            testCase.verifyError(@() setfield(quaPar,  'den', [1 2 3]), ...
                'RatPar:denMustBeTrivial'); %#ok<SFLD>

            % a Par may carry a genuine conic; a Pol may not
            quaPar.Ec = [1 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0];   % must not throw
            testCase.verifyError(@() setfield(quaPol, 'Ec', ones(3,6)), ...
                'RatPar:EcMustBeZero'); %#ok<SFLD>
            testCase.verifyError(@() setfield(ratPol,  'Ec', ones(3,6)), ...
                'RatPar:EcMustBeZero'); %#ok<SFLD>
        end

        function everyObjectIsFullyFormedWithRespectToTheLattice(testCase)
            % A QuaPol really IS a RatPol with unit denominators and a QuaPar with zero conics, so
            % it carries them rather than leaving them empty. This matters concretely: code that
            % legitimately reads p.den after an isa(p,'RatPol') test -- e.g. conjPieceCPLQ's
            % rational-piece guard -- would index into an empty matrix otherwise, now that a
            % QuaPol genuinely answers true to that test.
            V = [0 0; 1 0; 0 1]; E = [1 2 1; 2 3 1; 3 1 1]; F = [1 0; 1 0; 1 0];
            quaPol = QuaPol(V, E, [0 1 0 0 0 0], F);
            testCase.verifyEqual(quaPol.den, [0 0 1], 'AbsTol', 0);
            testCase.verifyTrue(all(quaPol.Ec(:) == 0));
            testCase.verifySize(quaPol.Ec, [3 6]);

            % and the guard that motivated it does not trip
            testCase.verifyTrue(isa(quaPol, 'RatPol'));
            testCase.verifyFalse(any(abs(quaPol.den(1,1:2)) > sqrt(eps)));

            ratPol = RatPol(V, E, [0 1 0 0 0 0], F);
            testCase.verifyTrue(all(ratPol.Ec(:) == 0));      % Pol pins the conics
            quaPar = QuaPar(V, E, zeros(3,6), [0 1 0 0 0 0], F);
            testCase.verifyEqual(quaPar.den, [0 0 1], 'AbsTol', 0);   % Qua pins the denominator
        end

        function allDataLivesOnRatParAlone(testCase)
            % No child declares a property of its own -- forced by MATLAB's multiple-inheritance
            % rules, and checked here so a future edit does not quietly reintroduce one and make
            % the diamond unconstructible.
            own = @(c) setdiff(properties(c), properties('RatPar'));
            for c = {'QuaPol','RatPol','QuaPar'}
                testCase.verifyEmpty(own(c{1}), ...
                    sprintf('%s must declare no properties of its own', c{1}));
            end
            % QuaParCPLQ is the one exception: it adds the symbolic payload it wraps.
            testCase.verifyEqual(own('QuaParCPLQ'), {'fnd'});
        end

        function ratParItselfCannotBeInstantiated(testCase)
            % RatPar is abstract: no operator produces a bare rational-on-parabolic result yet, so
            % there is nothing to construct. Making it concrete later is a one-word change.
            testCase.verifyError(@() RatPar(), 'MATLAB:class:abstract');
        end

        function isMeshedSeparatesTheSymbolicForm(testCase)
            % QuaParCPLQ holds cPLQ's own symbolic form with no reconstructed V/E/Ec/F/P mesh, so
            % its inherited mesh properties are empty for a real reason (the reconstruction is
            % still open, DESIGN.md II.5.1). It deliberately does NOT inherit from QuaPar -- that
            % would make isa(x,'QuaPar') true and invite callers to read a mesh it does not have --
            % but it does carry the same TRAITS, since it is mathematically quadratic-on-parabolic.
            V = [0 0; 1 0; 0 1]; E = [1 2 1; 2 3 1; 3 1 1]; F = [1 0; 1 0; 1 0];
            testCase.verifyTrue(QuaPol(V, E, [0 1 0 0 0 0], F).isMeshed());
            testCase.verifyTrue(RatPol(V, E, [0 1 0 0 0 0], F).isMeshed());
            testCase.verifyTrue(QuaPar(V, E, zeros(3,6), [0 1 0 0 0 0], F).isMeshed());

            sym = QuaParCPLQ();
            testCase.verifyFalse(sym.isMeshed());
            testCase.verifyFalse(isa(sym, 'QuaPar'));
            testCase.verifyTrue(isa(sym, 'Qua') && isa(sym, 'Par'));
        end

        function noArgConstructionWritesNothing(testCase)
            % The constructor protocol the diamond depends on: a no-argument constructor must
            % return a blank object, because MATLAB re-runs a shared base constructor once per
            % inheritance path and the second run would otherwise clobber the first path's writes.
            % See RatPar.m's CONSTRUCTOR PROTOCOL note.
            for c = {'QuaPol','RatPol','QuaPar'}
                obj = feval(c{1});
                testCase.verifyEmpty(obj.f,  sprintf('%s() must write nothing', c{1}));
                testCase.verifyEmpty(obj.V,  sprintf('%s() must write nothing', c{1}));
                testCase.verifyEmpty(obj.nv, sprintf('%s() must write nothing', c{1}));
            end
        end

        function conjReturnsARatParAcrossEverySupportedInputShape(testCase)
            % End to end, the actual promise: whatever supported shape goes in, a RatPar comes
            % out, and kind() names it. The shapes are conjCPLQ's own Cases A/B.
            fullQuad = QuaPol([1 0 1 0 0 0]);                       % Case A: full-domain PD
            gA = fullQuad.conj();
            testCase.verifyTrue(isa(gA, 'RatPar'));
            testCase.verifyEqual(gA.kind(), 'QuaPol');

            V = [0 0; 1 0; 0 1]; E = [1 2 1; 2 3 1; 3 1 1]; F = [1 0; 1 0; 1 0];
            gB = QuaPol(V, E, [0 1 0 0 0 0], F).conj();             % Case B: bounded triangle
            testCase.verifyTrue(isa(gB, 'RatPar'));
            testCase.verifyEqual(gB.kind(), 'QuaPar');
        end
    end

    methods (Static)
        function s = oneOfEach()
            % One instance of each concrete type, in the row order the lattice table above uses.
            V = [0 0; 1 0; 0 1]; E = [1 2 1; 2 3 1; 3 1 1]; F = [1 0; 1 0; 1 0];
            s = { QuaPol(V, E, [0 1 0 0 0 0], F), ...
                  RatPol(V, E, [0 1 0 0 0 0], F), ...
                  QuaPar(V, E, zeros(3,6), [0 1 0 0 0 0], F), ...
                  QuaParCPLQ() };
        end
    end
end
