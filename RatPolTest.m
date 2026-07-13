classdef RatPolTest < matlab.unittest.TestCase
    % Tests for RatPol (rational quadratic-over-linear functions on a polyhedral subdivision).
    % Covers: equivalence to QuaPoly when the denominator is 1, full-domain and bounded
    % rational evaluation, the NaN behaviour where the denominator vanishes, the evalRational
    % helper, and the (not-yet-implemented) operator guards.

    methods (Test)
        function denominatorOneMatchesQuaPoly(testCase)
            % den = [0 0 1] everywhere => plain polynomial => same values as QuaPoly.
            pPoly = QuaPoly.oneNorm();      % |x| + |y|
            pRat  = RatPol.oneNorm();       % same data; default denominators = 1
            [X,Y] = meshgrid(linspace(-1.7,1.7,21), linspace(-1.3,1.3,17));
            P = [X(:)+0.011, Y(:)+0.006];   % offset to avoid exact edges
            testCase.verifyEqual(pRat.eval(P), pPoly.eval(P), 'AbsTol', 1e-10);
        end

        function fullDomainRationalEval(testCase)
            % r(x,y) = (x^2 + y^2) / (x + y + 3) on all of R^2.
            p = RatPol.fullRational();
            S = [1 1; 0 0; 2 1; -1 2];
            num = S(:,1).^2 + S(:,2).^2;
            den = S(:,1) + S(:,2) + 3;
            testCase.verifyEqual(p.eval(S), num./den, 'AbsTol', 1e-12);
        end

        function denominatorZeroGivesNaN(testCase)
            % On the line x + y + 3 = 0 the denominator vanishes -> NaN (continuity TODO).
            p = RatPol.fullRational();
            v = p.eval([-1.5 -1.5]);        % x + y + 3 = 0
            testCase.verifyTrue(isnan(v));
        end

        function boundedRationalEval(testCase)
            % r on the unit square, +inf outside.
            p = RatPol.squareRational();
            testCase.verifyEqual(p.eval([0.5 0.5]), 0.5/4, 'AbsTol', 1e-12); % (0.25+0.25)/(0.5+0.5+3)
            testCase.verifyEqual(p.eval([0.25 0.75]), (0.0625+0.5625)/4, 'AbsTol', 1e-12);
            testCase.verifyEqual(p.eval([2 2]), Inf);      % outside the square
            testCase.verifyEqual(p.eval([-0.1 0.5]), Inf); % outside the square
        end

        function boundedRationalVectorized(testCase)
            p = RatPol.squareRational();
            S = [0.5 0.5; 2 2; 0.25 0.75; -0.1 0.5];
            expected = [0.125; Inf; (0.0625+0.5625)/4; Inf];
            testCase.verifyEqual(p.eval(S), expected, 'AbsTol', 1e-12);
        end

        function evalRationalHelper(testCase)
            % numerator x^2+y^2 (=[2 0 2 0 0 0]) over denominator x+y+3.
            z = RatPol.evalRational([2 0 2 0 0 0], [1 1 3], [1 1; 2 1]);
            testCase.verifyEqual(z, [2/5; 5/6], 'AbsTol', 1e-12);
        end

        function denConstructorValidation(testCase)
            % den must be nf x 3.
            V = [0 0;1 0;1 1;0 1]; E = [1 2 1;2 3 1;3 4 1;4 1 1]; F = [1 0;1 0;1 0;1 0];
            f = [2 0 2 0 0 0];
            testCase.verifyError(@() RatPol(V,E,f,F,[1 1]), 'RatPol:denSize');     % wrong width
            testCase.verifyError(@() RatPol(V,E,f,F,[1 1 3; 0 0 1]), 'RatPol:denSize'); % wrong height
        end

        function conjugateNotImplemented(testCase)
            p = RatPol.fullRational();
            testCase.verifyError(@() p.conj(),   'RatPol:conj:notImplemented');
            testCase.verifyError(@() p.biconj(), 'RatPol:biconj:notImplemented');
        end

        function scalarMulAndNegateScaleNumeratorOnly(testCase)
            p = RatPol.energy();   % 0.5*(x^2+y^2), full domain, denominator 1
            S = [1 2; -3 0.5];
            q = p.scalarMul(3);
            testCase.verifyEqual(q.den, p.den);   % denominator untouched
            testCase.verifyEqual(q.eval(S), 3*p.eval(S), 'AbsTol', 1e-12);
            r = p.negate();
            testCase.verifyEqual(r.eval(S), -p.eval(S), 'AbsTol', 1e-12);
        end
    end
end
