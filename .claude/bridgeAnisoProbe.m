function bridgeAnisoProbe()
% ITEM 3, rung 1: does the elimination that worked for the ISOTROPIC pair survive DIFFERENT
% Hessians? Observe the structure before building anything.
%
% WHAT MADE THE ISOTROPIC CASE EASY, and it was not luck. With y_k(s) = a_k + inv(H_k)*s, both
%
%   x - y_1 = (x - a_1) - inv(H_1)*s        and        y_2 - y_1 = (a_2 - a_1) + (inv(H_2)-inv(H_1))*s
%
% are AFFINE in s, so their determinant D is generally QUADRATIC in s. It came out linear before
% only because H_1 = I and H_2 = 2I are both multiples of the identity: the two s-parts are then
% parallel, and s wedge s vanishes. That is special, not general.
%
% So the question is whether D keeps a genuine s1*s2 or s2^2 term once the Hessians differ -- and
% if it does, whether the elimination still fits inside cofactor expansion, since eliminating s2
% from two conics leaves degree 4 in s1 and the second elimination is then an 8x8 Sylvester.

    fprintf('%-28s %-9s %-9s %-9s %s\n', 'H1, H2', 'D:s1*s2', 'D:s2^2', 'deg_s2 D', 'shape');
    fprintf('%s\n', repmat('-', 1, 74));
    cases = { {[1 0;0 1], [2 0;0 2], 'isotropic pair'}, ...
              {[1 0;0 1], [2 0;0 4], 'anisotropic'}, ...
              {[2 0;0 2], [4 0;0 4], 'both multiples of I'}, ...
              {[2 1;1 2], [4 2;2 4], 'proportional, non-diagonal'}, ...
              {[2 1;1 2], [3 0;0 5], 'general'} };
    for i = 1:numel(cases)
        H1 = cases{i}{1};  H2 = cases{i}{2};  nm = cases{i}{3};
        [c11, c02, dg] = dStructure(H1, H2);
        if dg <= 1, shape = 'LINEAR in s2 -- substitute';
        else,       shape = 'conic -- needs the big elimination'; end
        fprintf('%-28s %-9d %-9d %-9d %s\n', nm, c11, c02, dg, shape);
    end

    fprintf(['\nThe s1*s2 and s2^2 columns are the coefficients of D as a polynomial in s. D stays\n' ...
             'linear in s2 exactly when inv(H_2) - inv(H_1) is a multiple of inv(H_1) -- i.e. when\n' ...
             'the two Hessians are PROPORTIONAL. Otherwise it is a genuine conic.\n']);

    fprintf('\n=== so what does the general elimination cost? ===\n');
    fprintf('  eliminate s2 from (C,D): both conics -> 4x4 Sylvester, degree <= 4 in s1\n');
    fprintf('  eliminate s2 from (C,Z): both conics -> 4x4 Sylvester, degree <= 4 in s1\n');
    fprintf('  eliminate s1 from those two: 8x8 Sylvester -- AT ratQ.resultantVar''s cap of 8\n');
    fprintf('  cofactor expansion of 8x8 is 8! = %d products of 8 polynomial arrays.\n', factorial(8));
    fprintf('  Timing that is the next measurement; Bareiss is the documented alternative.\n');
end

% ------------------------------------------------------------------------------------------------
function [c11, c02, dg] = dStructure(H1, H2)
% The coefficients of s1*s2 and s2^2 in D = det[x - y_1(s), y_2(s) - y_1(s)], and D's degree in s2.
% Symbolic in x is unnecessary: those two coefficients come from the s-parts alone.
%
%   x - y_1 = u0 - M1*s      with M1 = inv(H1)
%   y_2-y_1 = v0 + M21*s     with M21 = inv(H2) - inv(H1)
%
% det(p + P*s, q + Q*s) has an s-quadratic part det-like in (P*s) and (Q*s):
%   (P*s) x (Q*s) = (P11 s1 + P12 s2)(Q21 s1 + Q22 s2) - (P21 s1 + P22 s2)(Q11 s1 + Q12 s2)
    M1  = inv(H1);
    M21 = inv(H2) - inv(H1);
    P = -M1;   Q = M21;
    % coefficient of s1*s2 and of s2^2 in the wedge
    c11 = P(1,1)*Q(2,2) + P(1,2)*Q(2,1) - P(2,1)*Q(1,2) - P(2,2)*Q(1,1);
    c02 = P(1,2)*Q(2,2) - P(2,2)*Q(1,2);
    c20 = P(1,1)*Q(2,1) - P(2,1)*Q(1,1);
    tol = 1e-12;
    c11 = roundz(c11, tol);  c02 = roundz(c02, tol);  c20 = roundz(c20, tol); %#ok<NASGU>
    if c02 ~= 0, dg = 2; elseif c11 ~= 0, dg = 1; else, dg = 1; end
end

function v = roundz(v, tol)
    if abs(v) < tol, v = 0; end
end
