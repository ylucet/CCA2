classdef bigZTest < matlab.unittest.TestCase
% bigZ (exact integers with MATLAB operators) and ratQ.exactFormula, which re-runs a double formula
% on bigZ inputs when the double run overflows.
%
% BUCKET: fast (small BigInteger arrays).

    methods (Test)

        function theOPERATORSAgreeWithDoublesWhereDoublesAreExact(t)
        % Differential: every supported operation on small integers, in doubles (the oracle) and in
        % bigZ, must give the same integers.
            rng(3);
            for k = 1:30
                a = randi([-50 50], 1, 4);  b = randi([-50 50], 1, 4);
                M = randi([-9 9], 2, 2);   s = randi([-7 7]);
                A = bigZ(a);  B = bigZ(b);  BM = bigZ(M);
                t.verifyEqual(double(A + B), a + b);
                t.verifyEqual(double(A - B), a - b);
                t.verifyEqual(double(A .* B), a .* b);
                t.verifyEqual(double(s * A), s * a);
                t.verifyEqual(double(A * s), a * s);
                t.verifyEqual(double(-A), -a);
                t.verifyEqual(double(BM * BM), M * M);
                t.verifyEqual(double(BM * b(1:2).'), M * b(1:2).');
                t.verifyEqual(double(conv(A, B)), conv(a, b));
                t.verifyEqual(double(A .^ 3), a .^ 3);
                t.verifyEqual(double(bigZ(s) ^ 4), s ^ 4);
                t.verifyEqual(double(A(2:3)), a(2:3));
                t.verifyEqual(double(A(end)), a(end));
                t.verifyEqual(double(BM(2,1)), M(2,1));
                t.verifyEqual(double([A, 0, B]), [a, 0, b]);
                t.verifyEqual(double([A; B]), [a; b]);
                t.verifyEqual(double(sum(A)), sum(a));
                t.verifyEqual(double(fliplr(A)), fliplr(a));
                C = A;  C(2) = 7;  c = a;  c(2) = 7;
                t.verifyEqual(double(C), c);
                C = A;  C(6) = s;  c = a;  c(6) = s;
                t.verifyEqual(double(C), c, 'growth fills with zeros');
            end
        end

        function everyREFUSALIsNamed(t)
        % what the class does not support is refused by name, never computed wrongly
            A = bigZ([1 2 3]);
            t.verifyError(@() bigZ('x'), 'bigZ:type');
            t.verifyError(@() bigZ(0.5), 'bigZ:notInteger');
            t.verifyError(@() bigZ(2^60), 'bigZ:inexactInput');
            t.verifyError(@() A + bigZ([1 2]), 'bigZ:dims');
            t.verifyError(@() A .^ -1, 'bigZ:power');
            t.verifyError(@() A ^ 2, 'bigZ:mpower');
            t.verifyError(@() subsasgn(A, substruct('{}', {1}), 5), 'bigZ:subsasgn');
        end

        function aVALUEPast2to53IsExact(t)
        % (2^40 + 1)^2 = 2^80 + 2^41 + 1 is not a double; bigZ must hold it exactly.
            v = bigZ(2^40 + 1) ^ 2;
            want = java.math.BigInteger('1208925819614629174706176').add( ...
                   java.math.BigInteger('2199023255553'));
            c = v.c;
            t.verifyTrue(c{1}.equals(want));
        end

        function exactFormulaReturnsTheREDUCEDResultWhenTheRawOneOverflows(t)
        % raw = 2^40 * [3 -5 7] * (2^30 + 7) passes 2^53; its content-free form is [3 -5 7].
            f = @(a, b, c) a * b * c;
            out = ratQ.exactFormula(f, 'test', 2^40, [3 -5 7], 2^30 + 7);
            t.verifyEqual(out, [3 -5 7]);
            t.verifyEqual(ratQ.exactFormula(f, 'test', 2, [3 -5 7], 5), [30 -50 70], ...
                'where doubles suffice the result is NOT reduced');
            t.verifyError(@() ratQ.exactFormula(@(a) a .^ 3, 'test', [2^20 + 1, 3]), ...
                'ratQ:overflow', 'a reduced result still past 2^53 is refused');
        end
    end
end
