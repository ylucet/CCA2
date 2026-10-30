function sg = triangleFaceSign(tf, P)
% triangleFaceSign  The EXACT sign of an integer polynomial P(s1, s2) at a triangle face's dual
%   point s0 -- rational (an exact integer evaluation) or named (dualPointSign).
%
% [input]  tf : as triangleFace returns it
%          P  : integer array, P(i+1,j+1) the coefficient of s1^i s2^j
% [output] sg : -1, 0 or +1
    if strcmp(tf.kind, 'algebraic')
        sg = dualPointSign(tf.dp, P);
        return
    end
    % s0 = s0n / s0d: sum P_ij s1n^i s2n^j s0d^(deg - i - j) has P(s0)'s sign (s0d > 0), in BigInteger
    B = @(v) java.math.BigInteger(sprintf('%.0f', v));
    P = ratQ.chk(double(P), 'polynomial');
    deg = (size(P,1) - 1) + (size(P,2) - 1);
    acc = java.math.BigInteger.ZERO;
    for i = 0:size(P,1)-1
        for j = 0:size(P,2)-1
            if P(i+1,j+1) == 0, continue, end
            t = B(P(i+1,j+1)).multiply(B(tf.s0n(1)).pow(i)).multiply(B(tf.s0n(2)).pow(j)) ...
                .multiply(B(tf.s0d).pow(deg - i - j));
            acc = acc.add(t);
        end
    end
    sg = acc.signum();
end
