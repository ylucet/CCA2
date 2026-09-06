function anisoPartsProbe()
% RUNG 1 for the anisotropic elimination's wrong answer. The pipeline RUNS -- big integers, 50 s,
% a reduced polynomial -- but the result does not vanish on the face. So one of the hand-built
% intermediates is wrong, and guessing which would be the third guess of the session.
%
% Check each one against a DIRECT numeric evaluation of what it is supposed to be. Every quantity
% below is defined twice: once as the integer polynomial the pipeline builds, once as the formula
% it is meant to encode, evaluated in doubles at random points. They must agree up to the recorded
% clearing factor.

    h1 = 1;  h2 = 1;  k1 = 2;  k2 = 4;
    pn = 13; qn = -7; ad = 10; cn = -1; cd = 2;
    a = [pn; qn]/ad;   c = cn/cd;
    H1 = diag([h1 h2]);  H2 = diag([k1 k2]);

    nv = 4;
    S1 = ratQ.monoN([1 0 0 0],1,nv);  S2 = ratQ.monoN([2 0 0 0],1,nv);
    X  = ratQ.monoN([0 1 0 0],1,nv);  Y  = ratQ.monoN([0 0 1 0],1,nv);
    one = ratQ.monoN([0 0 0 0],1,nv);

    L  = ad*h1*h2*k1*k2;
    e1 = L*(1/k1 - 1/h1);   e2 = L*(1/k2 - 1/h2);
    A1 = L/h1;              A2 = L/h2;
    a1n = pn*L/ad;          a2n = qn*L/ad;

    F1x = ratQ.polyAddN(X*L, S1*(-A1));
    As = ratQ.polyAddN(ratQ.polyMulN(F1x, one*e2), ...
                       ratQ.polyMulN(ratQ.polyAddN(one*a1n, S1*e1), one*A2));
    Bs = ratQ.polyAddN(ratQ.polyMulN(F1x, one*a2n), ...
                       ratQ.polyMulN(ratQ.polyAddN(one*a1n, S1*e1), Y*(-L)));

    rng(7);
    fprintf('%-26s %-14s %-14s %s\n', 'quantity', 'from the poly', 'direct', 'ratio');
    fprintf('%s\n', repmat('-', 1, 72));
    worstD = 0;
    for t = 1:5
        s1 = randn;  s2 = randn;  x = randn;  y = randn;
        % the polynomial's D, at scale L^2
        Dp = ev(As, s1, x, y, 0)*s2 + ev(Bs, s1, x, y, 0);
        % the direct determinant
        y1 = [s1/h1; s2/h2];
        y2 = a + [s1/k1; s2/k2];
        Dd = det([[x; y] - y1, y2 - y1]);
        r  = Dp / (Dd * L^2);
        worstD = max(worstD, abs(r - 1));
        if t <= 3
            fprintf('%-26s %-14.6g %-14.6g %.10f\n', sprintf('D at sample %d', t), Dp, Dd*L^2, r);
        end
    end
    fprintf('D: worst |ratio - 1| over 5 samples = %.3e\n\n', worstD);

    % ---- C, the tie conic -------------------------------------------------------------------------
    K  = 2*h1*k1*h2*k2*ad*cd;
    v  = ratQ.contentFree(ratQ.chk([K*(1/h1-1/k1)/2, K*(1/h2-1/k2)/2, -K*pn/ad, -K*qn/ad, K*cn/cd], 'C'));
    worstC = 0;
    for t = 1:5
        s1 = randn;  s2 = randn;
        Cp = v(1)*s1^2 + v(2)*s2^2 + v(3)*s1 + v(4)*s2 + v(5);
        s  = [s1; s2];
        Cd = (s.'*(H1\s))/2 - (s.'*a + (s.'*(H2\s))/2 - c);
        r  = Cp / Cd;
        if t == 1, ref = r; end
        worstC = max(worstC, abs(r/ref - 1));
        if t <= 3
            fprintf('%-26s %-14.6g %-14.6g %.6f\n', sprintf('C at sample %d', t), Cp, Cd, r);
        end
    end
    fprintf('C: worst deviation from a CONSTANT ratio = %.3e\n', worstC);
    fprintf(['   (C only needs to be a constant multiple of q1* - q2*, since it is set to zero;\n' ...
             '    a varying ratio would mean the conic is not that curve at all)\n']);
end

function v = ev(P, s1, x, y, z)
    v = 0;
    for i = 1:size(P,1)
        for j = 1:size(P,2)
            for k = 1:size(P,3)
                for l = 1:size(P,4)
                    cf = P(i,j,k,l);
                    if cf == 0, continue, end
                    v = v + cf * s1^(i-1) * x^(j-1) * y^(k-1) * z^(l-1);
                end
            end
        end
    end
end
