classdef exactQTest < matlab.unittest.TestCase
% Unit tests for exactQ, the arithmetic a sym-free CCA2 rests on.
%
% These are deliberately heavier than the type looks like it needs, because everything downstream
% inherits its guarantees: `region.merge` finds a shared facet by asking whether one constraint is
% the NEGATION of another, and Step 3's cell count depends on that answer being exact. A number
% type that is nearly right reproduces the ULP defect this whole exercise exists to remove.
%
% The oracle is the Symbolic Toolbox -- which is legitimate here and nowhere else: tests may use
% sym, the compute path may not.

    methods (Static)
        function verifyAgreesWithSym(tc, got, want, name)
        % `got` (exactQ) equals `want` (sym), exactly.
            tc.verifyTrue(isAlways(simplify(sym(got) - want) == 0, 'Unknown', 'false'), ...
                sprintf('%s: got %s, want %s', name, char(got), char(simplify(want))));
        end
    end

    methods (Test)

        % ---- construction and normalisation ---------------------------------------------

        function rationalsAreInLowestTerms(testCase)
        % Read through coeffOf rather than off the storage: the canonical FORM is the contract
        % (it is what makes `eq` and `isZero` exact), the field layout is not.
            a = exactQ(6, 8);
            testCase.verifyEqual(a.cn, int64(3));
            testCase.verifyEqual(a.cd, int64(4));
            b = exactQ(-6, -8);
            testCase.verifyEqual(b.cn, int64(3), 'a negative denominator must move to the numerator');
            testCase.verifyEqual(b.cd, int64(4));
            z = exactQ(0, 7);
            testCase.verifyTrue(isZero(z));
            testCase.verifyEmpty(z.m, 'zero is the EMPTY support, not a zero coefficient');
        end

        function aVanishingSurdPartNormalisesTheFieldAway(testCase)
        % Two values that are both rational must be COMPARABLE even if one was built in a field.
        % Without this, 0*sqrt(5) and 0 would live in different fields and `align` would refuse.
            a = exactQ(1, 2, 0, 1, 5);
            testCase.verifyTrue(a.isRational);
            testCase.verifyTrue(isZero(a.coeffOf(5)), 'no sqrt(5) term may survive');
            testCase.verifyTrue(a == exactQ(1,2));
        end

        function sqrtOfOneAndZeroFoldIntoTheRationalPart(testCase)
            testCase.verifyTrue(exactQ(1, 1, 2, 1, 1) == exactQ(3), 'sqrt(1) = 1');
            testCase.verifyTrue(exactQ(1, 1, 2, 1, 0) == exactQ(1), 'sqrt(0) = 0');
        end

        function squareFactorsAreExtractedSoFieldsMatch(testCase)
        % sqrt(8) and sqrt(2) are the SAME field. If squarefree reduction did not happen, values
        % that belong together would be refused by align -- the failure would look like a design
        % limit rather than the bug it is.
            a = exactQ.surd(8);
            testCase.verifyEqual(a.m, int64(2));
            testCase.verifyTrue(a.coeffOf(2) == exactQ(2), 'sqrt(8) = 2*sqrt(2)');
            b = exactQ.surd(2);
            c = a + b;
            exactQTest.verifyAgreesWithSym(testCase, c, 3*sqrt(sym(2)), 'sqrt(8)+sqrt(2)');
            testCase.verifyTrue(c.coeffOf(2) == exactQ(3), 'and they add as ONE term');
            testCase.verifyEqual(exactQ.surd(12).m, int64(3));
            testCase.verifyTrue(exactQ.surd(9).isRational, 'sqrt(9) is rational');
            testCase.verifyTrue(exactQ.surd(9) == exactQ(3));
        end

        function fromDoubleRefusesWhatItCannotRepresent(testCase)
            testCase.verifyTrue(exactQ(0.5) == exactQ(1, 2));
            testCase.verifyTrue(exactQ(-3) == exactQ(-3, 1));
            % pi is not a small rational and must RAISE rather than round
            testCase.verifyError(@() exactQ(pi), 'exactQ:inexact');
        end

        % ---- arithmetic against the symbolic oracle ------------------------------------

        function arithmeticMatchesSymbolicOnRandomValues(testCase)
        % The whole field, exercised: +, -, *, / on random elements of Q(sqrt(5)), each checked
        % against the Symbolic Toolbox. Random rather than hand-picked because the interesting
        % failures are in normalisation and cross-cancellation, not in any one case.
            rng(20260819);
            d = 5;
            for k = 1:60
                p = randi([-9 9], 1, 4); q = randi([1 7], 1, 4);
                x = exactQ(p(1), q(1), p(2), q(2), d);
                y = exactQ(p(3), q(3), p(4), q(4), d);
                sx = sym(p(1))/q(1) + sym(p(2))/q(2)*sqrt(sym(d));
                sy = sym(p(3))/q(3) + sym(p(4))/q(4)*sqrt(sym(d));
                exactQTest.verifyAgreesWithSym(testCase, x + y, sx + sy, 'plus');
                exactQTest.verifyAgreesWithSym(testCase, x - y, sx - sy, 'minus');
                exactQTest.verifyAgreesWithSym(testCase, x * y, sx * sy, 'times');
                if ~y.isZero
                    exactQTest.verifyAgreesWithSym(testCase, x / y, sx / sy, 'divide');
                end
            end
        end

        function divisionRationalisesTheDenominator(testCase)
            x = exactQ(1);
            y = exactQ(1, 1, 1, 1, 2);            % 1 + sqrt(2)
            z = x / y;                            % = -1 + sqrt(2)
            exactQTest.verifyAgreesWithSym(testCase, z, 1/(1 + sqrt(sym(2))), '1/(1+sqrt2)');
            testCase.verifyTrue(z * y == exactQ(1), 'z*y must be exactly 1');
        end

        function divisionByZeroRaises(testCase)
            testCase.verifyError(@() exactQ(1) / exactQ(0), 'exactQ:divideByZero');
        end

        % ---- sign, the operation everything downstream depends on -----------------------

        function signIsExactWhereFloatingPointWouldNotBe(testCase)
        % The case the type exists for. Successive convergents of sqrt(2) ALTERNATE about it, and
        % the margins shrink fast, so the sign of `p/q - sqrt(2)` is exactly the question a
        % tolerance gets wrong:
        %
        %   1393/985       = 1.41421319...  BELOW sqrt(2) by 3.7e-7   -> negative
        %   665857/470832  = 1.414213562... ABOVE sqrt(2) by 1.6e-12  -> positive
        %
        % Both are decided here by comparing a^2*bd^2 against b^2*d*ad^2 in integers -- 1940449
        % against 1940450 in the first case, one apart -- so no margin is too small.
        %
        % (This test caught its own author: the first expectation was written as positive. The
        % integers say otherwise, which is the point of having them.)
        %
        % BOTH ENTRY POINTS are pinned: `signExact` is the decision procedure, and `sign` is what
        % callers use -- it screens with a certified floating-point bound first, so a case this
        % tight is exactly where the screen must DECLINE to answer rather than guess.
            a = exactQ(1393, 985, -1, 1, 2);
            testCase.verifyEqual(sign(a), -1, '1393/985 is BELOW sqrt(2)');
            testCase.verifyEqual(signExact(a), -1);
            b = exactQ(665857, 470832, -1, 1, 2);
            testCase.verifyEqual(sign(b), 1, 'a 1.6e-12 margin is still a positive number');
            testCase.verifyEqual(signExact(b), 1);
            c = exactQ(-665857, 470832, 1, 1, 2);
            testCase.verifyEqual(sign(c), -1, 'and its negation is negative');
            testCase.verifyEqual(signExact(c), -1);

            % The screen's own contract: it answers only with a margin, and it must not answer
            % here. 1.6e-12 against terms of size 1.4 is about 1e-12 relative -- four orders
            % above eps -- so this one it CAN take; the pair below is the one it cannot.
            d = exactQ(1, 1, -1, 1, 2) + exactQ(-1, 1, 1, 1, 2);      % exactly 0
            testCase.verifyTrue(isZero(d), 'zero is decided by the representation, not the screen');
            testCase.verifyEqual(sign(d), 0);
        end

        function signExactAgreesWithTheScreenedSign(testCase)
        % `sign` screens in floating point and falls back to `signExact`. The two must never
        % disagree -- if they can, the screen's error bound is wrong, and the failure would be a
        % silent wrong verdict rather than a crash. Checked on values small enough for the exact
        % recursion to fit in int64 (it squares once per extension, so that is the constraint).
            rng(20260820);
            rad = [1 2 3 5 6 15];
            for trial = 1:40
                x = exactQ(0);
                for j = 1:3
                    x = x + exactQ(randi([-9 9]), randi([1 6])) * exactQ.surd(rad(randi(numel(rad))));
                end
                testCase.verifyEqual(sign(x), signExact(x), sprintf('trial %d: %s', trial, char(x)));
            end
        end

        function signIsZeroOnlyWhenTheValueIs(testCase)
            testCase.verifyEqual(sign(exactQ(0)), 0);
            testCase.verifyEqual(sign(exactQ.surd(2) - exactQ.surd(2)), 0);
            % 2*sqrt(2) is exactQ(0,1, 2,1, 2) -- rational part ZERO. Writing exactQ(2,1,1,1,2)
            % gives 2 + sqrt(2), which is a different number; the constructor's argument order is
            % (an, ad, bn, bd, d) and this is the easy way to misread it.
            testCase.verifyEqual(sign(exactQ.surd(8) - exactQ(0,1,2,1,2)), 0, ...
                'sqrt(8) - 2*sqrt(2) is exactly zero');
            testCase.verifyNotEqual(sign(exactQ(1, 1000000000)), 0, ...
                'a tiny rational is not zero');
        end

        function comparisonsFollowTheSign(testCase)
            a = exactQ(3, 2);                 % 1.5
            b = exactQ(0, 1, 1, 1, 2);        % sqrt(2) ~ 1.4142
            testCase.verifyTrue(a > b);
            testCase.verifyTrue(b < a);
            testCase.verifyTrue(a >= a);
            testCase.verifyTrue(a <= a);
            testCase.verifyFalse(a < b);
        end

        function equalityIsExactNotApproximate(testCase)
            a = exactQ(1, 3);
            b = exactQ(333333333, 1000000000);
            testCase.verifyFalse(a == b, '1/3 is not 0.333333333');
            testCase.verifyTrue(a == exactQ(2, 6));
        end

        % ---- the guard rails -----------------------------------------------------------

        function mixingTwoDifferentSurdsIsEXACTNotRefused(testCase)
        % THIS TEST USED TO ASSERT THE OPPOSITE, and the change is deliberate. It read
        %
        %     verifyError(@() exactQ.surd(2) + exactQ.surd(3), 'exactQ:fieldMismatch')
        %
        % on the argument that silently building Q(sqrt2,sqrt3) is how an exact type turns back
        % into a symbolic engine. The RULE was right and is kept -- nothing here promotes to an
        % arbitrary tower -- but the FIELD was too small, and that is a measured fact, not a
        % preference: the A.5 split of one triangle of the A.4/A.5 quadrilateral produces the
        % single coordinate sqrt(30)/12 - sqrt(15)/6 + 5/4, which no Q(sqrt(d)) can hold
        % (.claude/t1_multiquadratic_example.md, DECISIONS 2026-08-20). A type that refuses its
        % own target input is not carrying a design, it is failing.
        %
        % What replaces the refusal is a CLOSED family, not an open one: square roots of
        % squarefree integers, whose products stay in the family.
            x = exactQ.surd(2) + exactQ.surd(3);
            exactQTest.verifyAgreesWithSym(testCase, x, sqrt(sym(2)) + sqrt(sym(3)), 'sqrt2+sqrt3');
            testCase.verifyEqual(numel(x.m), 2, 'two independent terms, not a collapsed one');
            testCase.verifyTrue(exactQ(1) + exactQ.surd(3) == exactQ(1, 1, 1, 1, 3));
        end

        function overflowRaisesRatherThanWrapping(testCase)
        % MATLAB's int64 SATURATES silently. An exact type that saturates is a wrong type, so the
        % multiply is checked and a grown coefficient is reported as the signal it is.
            big = exactQ(intmax('int64')/3, 1);
            testCase.verifyError(@() big * big, 'exactQ:overflow');
        end

        function crossCancellationKeepsOrdinaryValuesSmall(testCase)
        % The defence that makes the overflow guard rare: cancel before multiplying. Without it a
        % chain of products of moderate fractions overflows in a handful of steps.
            x = exactQ(1);
            for k = 2:40
                x = x * exactQ(k, k+1);        % telescopes to 2/41
            end
            testCase.verifyTrue(x == exactQ(2, 41), 'the product must telescope exactly');
        end

        % ---- the values this project actually produces ----------------------------------

        function theA5CevianFootIsRepresentable(testCase)
        % 5/2 - sqrt(5)/2, recorded in ALGORITHM.md as the split foot that made attempt 4 work and
        % that T8 measured to be unavoidable. If the type cannot hold this, it is the wrong type.
            foot = exactQ(5, 2, -1, 2, 5);
            exactQTest.verifyAgreesWithSym(testCase, foot, sym(5)/2 - sqrt(sym(5))/2, 'A.5 foot');
            testCase.verifyEqual(foot.sign, 1);
            testCase.verifyEqual(double(foot), 2.5 - sqrt(5)/2, 'AbsTol', 1e-12);
        end

        function theOneUlpFacetCaseIsDecidedExactly(testCase)
        % THE defect this type exists to prevent, from DECISIONS.md: `4 - 2*sqrt(2)` arrived as two
        % doubles one ULP apart, `merge` asked "is this constraint the negation of that one",
        % got false, and Step 3's cell count grew without bound. In exact arithmetic the two are
        % the same number and the negation test answers yes.
            a = exactQ(4, 1, -2, 1, 2);
            b = exactQ(4, 1, -2, 1, 2);
            testCase.verifyTrue(a == b);
            testCase.verifyTrue(isZero(a - b));
            testCase.verifyTrue((-a) == exactQ(-4, 1, 2, 1, 2), 'negation must be exact');
            % and the doubles really are indistinguishable at the ULP level
            testCase.verifyEqual(double(a), 4 - 2*sqrt(2), 'AbsTol', 0);
        end

        function theRotationFactorIsExact(testCase)
        % 1/sqrt(2), from biconjCPLQ's 45-degree frame, where reading the coefficients back at 10
        % decimals was measured to put 9.5e-12 into co(x*y) at the origin.
            r = exactQ(1) / exactQ.surd(2);
            exactQTest.verifyAgreesWithSym(testCase, r, 1/sqrt(sym(2)), '1/sqrt(2)');
            testCase.verifyTrue(r * r == exactQ(1, 2), '(1/sqrt2)^2 must be exactly 1/2');
        end

        function powersStayExact(testCase)
            x = exactQ(1, 1, 1, 1, 2);          % 1 + sqrt(2)
            exactQTest.verifyAgreesWithSym(testCase, x^3, (1 + sqrt(sym(2)))^3, '(1+sqrt2)^3');
            testCase.verifyTrue(x^0 == exactQ(1));
            exactQTest.verifyAgreesWithSym(testCase, x^-1, 1/(1 + sqrt(sym(2))), 'negative power');
        end

        function twoDifferentSurdsCombineInsteadOfRaising(testCase)
        % THE T1 COUNTEREXAMPLE, measured 2026-08-20 and written up in
        % .claude/t1_multiquadratic_example.md. One triangle of the A.4/A.5 quadrilateral --
        % conv{(5/2,3/2),(0,0),(1/2,1)}, three convex edges, so [COAP] A.5 -- splits into
        % sub-triangles one of whose vertices has x-coordinate
        %
        %       sqrt(30)/12 - sqrt(15)/6 + 5/4
        %
        % A SINGLE NUMBER needing two extensions. The old type raised exactQ:fieldMismatch on it
        % by design, and that design cannot carry this pipeline: there is no Q(sqrt(d)) holding
        % that coordinate, so no caller can route around it by keeping cells apart.
            a = exactQ.surd(15);
            b = exactQ.surd(30);
            s = a + b;
            testCase.verifyEqual(double(s), sqrt(15) + sqrt(30), 'AbsTol', 1e-12);
            testCase.verifyEqual(sign(s), 1);
            testCase.verifyEqual(sign(a - b), -1, 'sqrt(15) < sqrt(30)');

            % The exact coordinate itself, and its sign, with no floating point in the decision.
            v = exactQ.surd(30)/exactQ(12) - exactQ.surd(15)/exactQ(6) + exactQ(5,4);
            testCase.verifyEqual(double(v), sqrt(30)/12 - sqrt(15)/6 + 5/4, 'AbsTol', 1e-12);
            testCase.verifyEqual(sign(v), 1);
        end

        function aProductOfSurdsLandsInAThirdFieldExactly(testCase)
        % sqrt(15)*sqrt(30) = sqrt(450) = 15*sqrt(2). Squarefree radicands are closed under
        % multiplication up to a rational factor, which is exactly why the field to carry is the
        % MULTIQUADRATIC one and not a tower of arbitrary algebraic numbers.
            p = exactQ.surd(15) * exactQ.surd(30);
            testCase.verifyTrue(p == exactQ(0,1,15,1,2));
            testCase.verifyEqual(double(p), 15*sqrt(2), 'AbsTol', 1e-9);

            % ... and the basis really is a basis: sqrt(2)*sqrt(3) and sqrt(6) are the SAME value,
            % so their difference is exactly zero -- not "zero to a tolerance".
            testCase.verifyTrue(isZero(exactQ.surd(2)*exactQ.surd(3) - exactQ.surd(6)));
            testCase.verifyTrue(isZero((exactQ.surd(2) + exactQ.surd(3))^2 ...
                                       - (exactQ(5) + exactQ(2)*exactQ.surd(6))));
        end

        function zeroIsDECIDEDNotApproximated(testCase)
        % The property everything downstream rests on: sqrt(m) over distinct squarefree m are
        % linearly INDEPENDENT over Q, so a value is zero exactly when every coefficient is, and
        % a combination that is merely SMALL is still decided as nonzero. Both directions matter
        % -- region's redundancy and facet tests ask "is this exactly 0" and a wrong yes merges
        % two regions that do not share a facet.
            z = exactQ.surd(2) + exactQ.surd(3) - exactQ.surd(2) - exactQ.surd(3);
            testCase.verifyTrue(isZero(z));
            testCase.verifyEqual(sign(z), 0);

            tiny = exactQ.surd(2) + exactQ.surd(3) - exactQ(3146264,1000000);   % ~2.6e-7
            testCase.verifyFalse(isZero(tiny));
            testCase.verifyEqual(sign(tiny), 1, ...
                'sqrt(2) + sqrt(3) = 3.14626436..., so this is positive by 2.6e-7');
            testCase.verifyEqual(sign(exactQ.surd(2) + exactQ.surd(3) - exactQ(3146265,1000000)), ...
                -1, 'and negative one unit further out');
        end

        function signIsExactAcrossSeveralExtensions(testCase)
        % The sign recursion: x = a + b*sqrt(p) with a, b in the field of the REMAINING primes;
        % same-sign parts answer immediately, opposite-sign parts compare a^2 against b^2*p in
        % that smaller field. Checked against vpa, which is an oracle here and not on any compute
        % path.
            cases = { exactQ.surd(2) + exactQ.surd(3) - exactQ.surd(5) - exactQ(1,2), ...
                      exactQ.surd(6) - exactQ.surd(2) - exactQ.surd(3) + exactQ(1), ...
                      exactQ.surd(15) - exactQ.surd(30)/exactQ(2) - exactQ(1,8), ...
                      exactQ(7,3)*exactQ.surd(10) - exactQ.surd(2)*exactQ.surd(5)*exactQ(7,3) };
            for k = 1:numel(cases)
                want = double(sign(double(vpa(sym(cases{k}), 40))));
                testCase.verifyEqual(sign(cases{k}), want, sprintf('case %d', k));
            end
            testCase.verifyEqual(sign(cases{4}), 0, ...
                'sqrt(10) and sqrt(2)*sqrt(5) are the same number, so this one is exactly 0');
        end

        function divisionRationalisesAcrossEveryExtension(testCase)
        % 1/x in Q(sqrt(p1),...,sqrt(pk)) needs the product of the NON-TRIVIAL conjugates, not a
        % single one: multiplying by the sqrt(p)-flip removes p and leaves the rest, so the
        % rationalisation is one pass per prime.
            x = exactQ(1) + exactQ.surd(2) + exactQ.surd(3);
            r = x / x;
            testCase.verifyTrue(r == exactQ(1));
            y = inv(x);
            testCase.verifyTrue(isZero(x*y - exactQ(1)));
            testCase.verifyEqual(double(y), 1/(1 + sqrt(2) + sqrt(3)), 'AbsTol', 1e-12);

            w = (exactQ(5,4) + exactQ.surd(30)/exactQ(12)) / (exactQ(2) - exactQ.surd(15)/exactQ(6));
            testCase.verifyEqual(double(w), (5/4 + sqrt(30)/12)/(2 - sqrt(15)/6), 'AbsTol', 1e-12);
        end

        function arithmeticMatchesSymbolicOverSeveralExtensions(testCase)
        % The same differential check the single-extension test makes, over the multiquadratic
        % field: random rational combinations of sqrt(m) for m in {1,2,3,5,6,10,15,30}, against
        % `sym` as the oracle. `sym` is the thing being REPLACED, so it is the right oracle here
        % and it appears on no compute path.
            rng(20260820);
            rad = [1 2 3 5 6 10 15 30];
            nt = 3;
            for trial = 1:20
                [x, sx] = deal(exactQ(0), sym(0));
                [y, sy] = deal(exactQ(0), sym(0));
                for j = 1:nt
                    ix = rad(randi(numel(rad))); px = randi([-9 9]); qx = randi([1 7]);
                    iy = rad(randi(numel(rad))); py = randi([-9 9]); qy = randi([1 7]);
                    x  = x  + exactQ(px, qx) * exactQ.surd(ix);
                    sx = sx + sym(px)/qx * sqrt(sym(ix));
                    y  = y  + exactQ(py, qy) * exactQ.surd(iy);
                    sy = sy + sym(py)/qy * sqrt(sym(iy));
                end
                exactQTest.verifyAgreesWithSym(testCase, x + y, sx + sy, 'mq plus');
                exactQTest.verifyAgreesWithSym(testCase, x - y, sx - sy, 'mq minus');
                exactQTest.verifyAgreesWithSym(testCase, x * y, sx * sy, 'mq times');
                testCase.verifyEqual(sign(x - y), double(sign(double(vpa(sx - sy, 40)))), ...
                    sprintf('trial %d: sign of the difference', trial));
                if ~isZero(y)
                    exactQTest.verifyAgreesWithSym(testCase, x / y, sx / sy, 'mq divide');
                end
            end
        end
    end
end
