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
            a = exactQ(6, 8);
            testCase.verifyEqual(a.an, int64(3));
            testCase.verifyEqual(a.ad, int64(4));
            b = exactQ(-6, -8);
            testCase.verifyEqual(b.an, int64(3), 'a negative denominator must move to the numerator');
            testCase.verifyEqual(b.ad, int64(4));
            z = exactQ(0, 7);
            testCase.verifyEqual(z.ad, int64(1), 'zero must normalise its denominator to 1');
        end

        function aVanishingSurdPartNormalisesTheFieldAway(testCase)
        % Two values that are both rational must be COMPARABLE even if one was built in a field.
        % Without this, 0*sqrt(5) and 0 would live in different fields and `align` would refuse.
            a = exactQ(1, 2, 0, 1, 5);
            testCase.verifyTrue(a.isRational);
            testCase.verifyEqual(a.d, int64(0));
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
            testCase.verifyEqual(a.d, int64(2));
            testCase.verifyEqual(a.bn, int64(2), 'sqrt(8) = 2*sqrt(2)');
            b = exactQ.surd(2);
            c = a + b;                                   % must not raise
            exactQTest.verifyAgreesWithSym(testCase, c, 3*sqrt(sym(2)), 'sqrt(8)+sqrt(2)');
            testCase.verifyEqual(exactQ.surd(12).d, int64(3));
            testCase.verifyEqual(exactQ.surd(9).d, int64(0), 'sqrt(9) is rational');
        end

        function fromDoubleRefusesWhatItCannotRepresent(testCase)
            testCase.verifyEqual(exactQ(0.5).an, int64(1));
            testCase.verifyEqual(exactQ(0.5).ad, int64(2));
            testCase.verifyEqual(exactQ(-3).an, int64(-3));
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
            a = exactQ(1393, 985, -1, 1, 2);
            testCase.verifyEqual(sign(a), -1, '1393/985 is BELOW sqrt(2)');
            b = exactQ(665857, 470832, -1, 1, 2);
            testCase.verifyEqual(sign(b), 1, 'a 1.6e-12 margin is still a positive number');
            c = exactQ(-665857, 470832, 1, 1, 2);
            testCase.verifyEqual(sign(c), -1, 'and its negation is negative');
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

        function mixingTwoDifferentSurdsRaises(testCase)
        % By DESIGN, and the message has to say so: silently building Q(sqrt2, sqrt3) is how an
        % exact type turns back into a symbolic engine. Whatever needs both is the thing to look at.
            testCase.verifyError(@() exactQ.surd(2) + exactQ.surd(3), 'exactQ:fieldMismatch');
            % ...but a RATIONAL combines with anything
            testCase.verifyWarningFree(@() exactQ(1) + exactQ.surd(3));
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
    end
end
