classdef dualPointTest < matlab.unittest.TestCase
% dualPoint / dualPointSign: the exact name of a degree-4 dual point and the sign of a polynomial
% there -- step (a) of Example 2's triangle faces (TODO item 1, DECISIONS 2026-10-30 (fifth)).
%
% BUCKET: fast. Integer polynomial arithmetic and Sturm counts; no symbolic call.

    methods (Test)

        function aRATIONALPointHasExactZeroAndNonzeroSigns(t)
        % the circle s1^2 + s2^2 = 2 meets the line s1 = s2 at (1,1) and (-1,-1)
            dp = dualPoint([1 0 1 0 0 -2], [0 0 0 1 -1 0], [1 1]);
            t.verifyEqual(dp.s, [1; 1], 'AbsTol', 1e-15);
            t.verifyEqual(dualPointSign(dp, [-1; 1]), 0, 's1 - 1 vanishes');
            t.verifyEqual(dualPointSign(dp, [-2 1; 1 0]), 0, 's1 + s2 - 2 vanishes');
            t.verifyEqual(dualPointSign(dp, [-1; 2]), 1, '2 s1 - 1 > 0');
            t.verifyEqual(dualPointSign(dp, [3 -4]), -1, '3 - 4 s2 < 0');
        end

        function aQUARTICPointIsZeroOnItsOwnIdeal(t)
        % s1^2 + s2^2 = 4 and s2 = s1^2 - 1: the eliminant u^4 - u^2 - 3 is irreducible, so s0 is of
        % degree 4. Every combination P*C1 + P'*C2 vanishes there EXACTLY, the identity a double
        % evaluation can only approximate.
            c1 = [1 0 1 0 0 -4];  c2 = [1 0 0 0 -1 -1];
            u = sqrt((1 + sqrt(13))/2);
            for s = [u, -u; u^2 - 1, u^2 - 1]
                dp = dualPoint(c1, c2, s);
                t.verifyEqual(dp.s, s, 'AbsTol', 1e-12);
                t.verifyEqual(dualPointSign(dp, conicArr(c1)), 0);
                t.verifyEqual(dualPointSign(dp, mulArr([2 -1; 3 0], conicArr(c1)) + ...
                                                mulArr([0 5; 0 0], conicArr(c2))), 0);
                t.verifyEqual(dualPointSign(dp, [0 1; 1 0]), sign(s(1) + s(2)), 's1 + s2');
            end
        end

        function twoRootsOnOneVERTICALLineAreSeparatedByAShear(t)
        % s1^2 + s2^2 = 5 and s1^2 - s2^2 = 3 meet at (+-2, +-1): (2,1) and (2,-1) share s1, so
        % s2 is not a function of s1 on the eliminant's root and the name must shear.
            c1 = [1 0 1 0 0 -5];  c2 = [1 0 -1 0 0 -3];
            up = dualPoint(c1, c2, [2 1]);  dn = dualPoint(c1, c2, [2 -1]);
            t.verifyEqual([up.s, dn.s], [2 2; 1 -1], 'AbsTol', 1e-12);
            t.verifyNotEqual(up.k, 0, 'the unsheared name cannot tell (2,1) from (2,-1)');
            t.verifyEqual(dualPointSign(up, [0 1]), 1);
            t.verifyEqual(dualPointSign(dn, [0 1]), -1);
            t.verifyEqual(dualPointSign(dn, [1 1]), 0, '1 + s2 = 0 at (2,-1)');
        end

        function theSignDoesNotDependOnTheORDEROfTheConics(t)
        % With the line first, v = -S0/S1 has S1 = -1 < 0: the sign of S1 enters every polynomial
        % of odd degree in v, and dropping it flips sign(s2) at (1,1).
            line = [0 0 0 1 -1 0];  circ = [1 0 1 0 0 -2];
            for s = [1 -1; 1 -1]
                a = dualPoint(circ, line, s);  b = dualPoint(line, circ, s);
                for Q = {[0 1], [0 1; 1 0], [-1 0 1], [0 0 0 1]}
                    t.verifyEqual(dualPointSign(b, Q{1}), dualPointSign(a, Q{1}), mat2str(Q{1}));
                    t.verifyEqual(dualPointSign(b, Q{1}), sign(evalArr(Q{1}, s)), mat2str(Q{1}));
                end
            end
        end

        function aVERTICALPairIsSeparatedWhereS1VanishesONLYAtTheRoot(t)
        % (s1-2)*s2 added to the second conic keeps (2,1) and (2,-1) on it, but makes S1 a nonzero
        % polynomial that vanishes at u = 2 -- the case the shear test decides by isolatedRootSign,
        % not by S1 being identically zero.
            c1 = [1 0 1 0 0 -5];  c2 = [1 1 -1 0 -2 -3];
            up = dualPoint(c1, c2, [2 1]);  dn = dualPoint(c1, c2, [2 -1]);
            t.verifyEqual([up.s, dn.s], [2 2; 1 -1], 'AbsTol', 1e-12);
            t.verifyEqual([dualPointSign(up, [0 1]), dualPointSign(dn, [0 1])], [1 -1]);
        end

        function signsAgreeWithDoublesOnRANDOMConics(t)
        % Property: wherever a double evaluation is clear of rounding, the exact sign agrees with
        % it; and a polynomial in the conics' ideal is exactly zero.
            rng(41);
            nPts = 0;
            for trial = 1:60
                c1 = randi([-4 4], 1, 6);  c2 = randi([-4 4], 1, 6);
                if all(c1(1:3) == 0) || all(c2(1:3) == 0), continue, end
                try
                    P = conicMeet(c1, c2);
                catch
                    continue
                end
                for i = 1:size(P,1)
                    try
                        dp = dualPoint(c1, c2, P(i,:));
                    catch ME
                        t.verifyTrue(ismember(ME.identifier, {'dualPoint:commonComponent', ...
                            'dualPoint:noMatch', 'dualPoint:noShear'}), ME.message);
                        continue
                    end
                    nPts = nPts + 1;
                    s = dp.s;
                    dq = dualPoint(c2, c1, s);          % the other order names the same point
                    t.verifyEqual(dualPointSign(dp, mulArr(randi([-3 3], 2, 2), conicArr(c1)) ...
                                                    + mulArr(randi([-3 3], 2, 2), conicArr(c2))), 0);
                    for r = 1:4
                        Q = randi([-6 6], 3, 3);
                        v = evalArr(Q, s);
                        if abs(v) > 1e-6 * sum(abs(Q(:))) * max(1, norm(s))^4
                            t.verifyEqual(dualPointSign(dp, Q), sign(v), ...
                                sprintf('trial %d point %d: Q = %s', trial, i, mat2str(Q)));
                            t.verifyEqual(dualPointSign(dq, Q), sign(v), 'swapped order');
                        end
                    end
                end
            end
            t.verifyGreaterThan(nPts, 40, 'the random draw must reach real points');
        end

        function namingIsREFUSEDByNameWhereNoPointIsDefined(t)
            t.verifyError(@() dualPoint([0 0 0 1 0 -1], [0 0 0 0 1 -1], [1 1]), 'dualPoint:twoLines');
            t.verifyError(@() dualPoint([1 0 1 0 0 -1], [2 0 2 0 0 -2], [1 0]), ...
                          'dualPoint:commonComponent');
            t.verifyError(@() dualPoint([1 0 1 0 0 -2], [0 0 0 1 -1 0], [3 0]), 'dualPoint:noMatch');
            t.verifyError(@() dualPoint([1 0 1 0 0 -2], [0 0 0 1 -1 0], [NaN 1]), 'dualPoint:approx');
            dp = dualPoint([1 0 1 0 0 -2], [0 0 0 1 -1 0], [1 1]);
            t.verifyError(@() dualPointSign(dp, ones(2,2,2)), 'dualPointSign:shape');
        end
    end
end

function A = conicArr(c)
% [a b c d e f] as the polyN array A(i+1,j+1) = coefficient of s1^i s2^j
    A = [c(6) c(5) c(3); c(4) c(2) 0; c(1) 0 0];
end

function C = mulArr(A, B)
    C = conv2(A, B);
end

function v = evalArr(Q, s)
    v = 0;
    for i = 0:size(Q,1)-1
        for j = 0:size(Q,2)-1
            v = v + Q(i+1,j+1) * s(1)^i * s(2)^j;
        end
    end
end
