function bridgeAnisoBuild()
% ITEM 3: the elimination for pieces with DIFFERENT Hessians, built and checked end to end.
%
% .claude/bridgeAnisoProbe.m established that D = det[x - y1(s), y2(s) - y1(s)] stays LINEAR in s2
% for every Hessian pair that is simultaneously diagonal or proportional -- much wider than the
% isotropic pair the method was derived on. This builds that case and checks it against the primal
% oracle, which forms no conjugate and knows nothing about the elimination.
%
% WHAT CHANGES FROM THE ISOTROPIC CASE, and it is not free. There, D's s2-coefficient A was free of
% s1, so clearing A^2 left degree 2 in s1 and a 4x4 resultant. Here A depends on s1, so clearing
% A^2 takes the degree in s1 to 4 and the resultant is 8x8 -- exactly resultantVar's cap, and the
% reason that cap was measured rather than guessed.
%
% THE INSTANCE. q1(v) = v'H1 v/2 and q2(v) = (v-a)'H2(v-a)/2 + c, both DIAGONAL:
%   H1 = diag(h1,h2), H2 = diag(k1,k2), a = (pn,qn)/ad, c = cn/cd.
% Conjugates: qk*(s) = <s,ak> + s'inv(Hk)s/2 - ck, with a1 = 0 and c1 = 0.

    h1 = 1;  h2 = 1;          % H1 = I
    k1 = 2;  k2 = 4;          % H2 = diag(2,4) -- the anisotropic pair
    pn = 13; qn = -7; ad = 10;
    cn = -1; cd = 2;

    fprintf('H1 = diag(%d,%d), H2 = diag(%d,%d), a = (%d,%d)/%d, c = %d/%d\n\n', ...
            h1, h2, k1, k2, pn, qn, ad, cn, cd);

    [p, degS1] = buildAniso(h1, h2, k1, k2, pn, qn, ad, cn, cd);
    fprintf('degree in s1 after clearing A^2: %d  (isotropic case was 2)\n', degS1);
    fprintf('p: degrees (x,y,z) = (%d,%d,%d), %d terms, max |coeff| %.4g\n', ...
            size(p,2)-1, size(p,3)-1, size(p,4)-1, nnz(p), max(abs(p(:))));

    % ---- the ORACLE, which shares no algebra with the above ---------------------------------------
    % co f by the primal definition through the perspective inf-convolution: for two quadratics the
    % inner infimal convolution is closed form, leaving a 1-D minimisation over t. Both Hessians are
    % diagonal, so it vectorises.
    a = [pn; qn]/ad;   c = cn/cd;
    ts = linspace(0, 1, 400001).';
    m1 = (1-ts)/h1 + ts/k1;          % the inf-convolved Hessian's inverse, diagonally
    m2 = (1-ts)/h2 + ts/k2;

    fprintf('\npoints where the ORACLE says the bridge is active (0 < t < 1)\n');
    rng(20260906);
    X = 6*(rand(400,2) - 0.5);
    nb = 0;  worst = 0;
    for i = 1:size(X,1)
        d1 = X(i,1) - ts*a(1);   d2 = X(i,2) - ts*a(2);
        vals = 0.5*(d1.^2./m1 + d2.^2./m2) + ts*c;
        [z, j] = min(vals);
        if ts(j) <= 1e-6 || ts(j) >= 1-1e-6, continue, end
        nb = nb + 1;
        [v, sc] = evalRel(p, X(i,1), X(i,2), z);
        worst = max(worst, abs(v)/sc);
        if nb <= 6
            fprintf('   x=(%7.4f,%7.4f)  t=%.4f  z=%9.5f   relative |p| = %.2e\n', ...
                    X(i,1), X(i,2), ts(j), z, abs(v)/sc);
        end
    end
    fprintf('   %d bridge points, worst relative |p| = %.3e\n', nb, worst);

    % ---- and it must NOT vanish off the face ------------------------------------------------------
    fprintf('\nthe same points with z perturbed (p must NOT vanish)\n');
    shown = 0;
    for i = 1:size(X,1)
        d1 = X(i,1) - ts*a(1);   d2 = X(i,2) - ts*a(2);
        vals = 0.5*(d1.^2./m1 + d2.^2./m2) + ts*c;
        [z, j] = min(vals);
        if ts(j) <= 1e-6 || ts(j) >= 1-1e-6, continue, end
        [v, sc] = evalRel(p, X(i,1), X(i,2), z + 0.05);
        fprintf('   x=(%7.4f,%7.4f)  relative |p| = %.2e\n', X(i,1), X(i,2), abs(v)/sc);
        shown = shown + 1;
        if shown >= 4, break, end
    end
end

% ------------------------------------------------------------------------------------------------
function [p, degS1] = buildAniso(h1, h2, k1, k2, pn, qn, ad, cn, cd)
% Variables (s1, x, y, z). Everything is cleared to integers; the clearing factors are recorded in
% the comments because each one is a place an exactness claim could quietly fail.
    nv = 4;
    S1 = ratQ.monoN([1 0 0 0],1,nv);  S2 = ratQ.monoN([2 0 0 0],1,nv);
    X  = ratQ.monoN([0 1 0 0],1,nv);  Y  = ratQ.monoN([0 0 1 0],1,nv);
    Zv = ratQ.monoN([0 0 0 1],1,nv);  one = ratQ.monoN([0 0 0 0],1,nv);

    % ---- C = q1* - q2* -----------------------------------------------------------------------
    % q1*(s) = s1^2/(2h1) + s2^2/(2h2);  q2*(s) = <s,a> + s1^2/(2k1) + s2^2/(2k2) - c
    % C = s1^2*(1/h1 - 1/k1)/2 + s2^2*(1/h2 - 1/k2)/2 - <s,a> + c
    % Cleared by 2*h1*k1*h2*k2*ad*cd (positive), so each coefficient is an integer.
    K  = 2*h1*k1*h2*k2*ad*cd;
    ca = K*(1/h1 - 1/k1)/2;          % coefficient of s1^2
    cb = K*(1/h2 - 1/k2)/2;          % coefficient of s2^2
    c1 = -K*pn/ad;                   % coefficient of s1
    c2 = -K*qn/ad;                   % coefficient of s2
    c0 =  K*cn/cd;
    v  = ratQ.contentFree(ratQ.chk([ca cb c1 c2 c0], 'tie conic'));
    ca = v(1); cb = v(2); c1 = v(3); c2 = v(4); c0 = v(5);

    % ---- D = det[x - y1(s), y2(s) - y1(s)] ---------------------------------------------------
    % y1 = (s1/h1, s2/h2),  y2 = a + (s1/k1, s2/k2)
    % y2 - y1 = (a1 + s1*(1/k1 - 1/h1),  a2 + s2*(1/k2 - 1/h2))
    % D = (x - s1/h1)*(a2 + s2*e2) - (y - s2/h2)*(a1 + s1*e1),  e_i = 1/k_i - 1/h_i
    % Cleared by ad*h1*h2*k1*k2.
    L  = ad*h1*h2*k1*k2;
    e1 = L*(1/k1 - 1/h1);   e2 = L*(1/k2 - 1/h2);
    A1 = L/h1;              A2 = L/h2;              % so y1 = (A1*s1, A2*s2)/L
    % D*L^2/... -- rather than track it symbolically, build D term by term over the integers:
    %   D = (L*x - A1*s1)*(L*a2 + e2*s2*L/L) ... which is fiddly; instead scale each factor by L:
    F1x = ratQ.polyAddN(X*L, S1*(-A1));                       % L*x - A1*s1
    F1y = ratQ.polyAddN(Y*L, S2_placeholder(nv)*0);           % L*y   (s2 is handled below)
    % D is linear in s2 for diagonal Hessians, so write D = As*s2 + Bs with s2 eliminated by hand:
    %   D = F1x*(a2*L + e2*s2) - (L*y - A2*s2)*(a1*L + e1*s1)   all over L^2
    G1 = ratQ.polyAddN(one*(qn*L*L/ad), ratQ.polyMulN(one*0, one));   % a2*L^2/L handled in As/Bs
    %#ok<*NASGU>
    a1n = pn*L/ad;   a2n = qn*L/ad;                            % a*L, integers since ad | L
    % As (coefficient of s2) and Bs (the rest):
    As = ratQ.polyAddN(ratQ.polyMulN(F1x, one*e2), ratQ.polyMulN(ratQ.polyAddN(one*a1n, S1*e1), one*A2));
    Bs = ratQ.polyAddN(ratQ.polyMulN(F1x, one*a2n), ...
                       ratQ.polyMulN(ratQ.polyAddN(one*a1n, S1*e1), Y*(-L)));
    As = ratQ.contentFree(As);   Bs = ratQ.contentFree(Bs);
    degS1 = max(size(As,1), size(Bs,1)) - 1;

    % ---- substitute s2 = -Bs/As into C and Z, clearing As^2 -----------------------------------
    As2 = ratQ.polyMulN(As, As);   Bs2 = ratQ.polyMulN(Bs, Bs);   AB = ratQ.polyMulN(As, Bs);
    Ct = ratQ.polyAddN( ...
            ratQ.polyAddN(ratQ.polyMulN(S2, As2)*ca, Bs2*cb), ...
            ratQ.polyAddN(ratQ.polyMulN(S1, As2)*c1, ...
                ratQ.polyAddN(AB*(-c2), As2*c0)));
    % Z: z = <s,x> - q1*(s) = s1*x + s2*y - s1^2/(2h1) - s2^2/(2h2)
    %    => 2*h1*h2*z = 2*h1*h2*(s1 x + s2 y) - h2 s1^2 - h1 s2^2 ; clear by As^2 as well
    Zt = ratQ.polyAddN( ...
            ratQ.polyAddN(ratQ.polyMulN(S2, As2)*(-h2), Bs2*(-h1)), ...
            ratQ.polyAddN(ratQ.polyMulN(ratQ.polyMulN(S1, X), As2)*(2*h1*h2), ...
                ratQ.polyAddN(ratQ.polyMulN(ratQ.polyMulN(Y, AB), one)*(-2*h1*h2), ...
                              ratQ.polyMulN(Zv, As2)*(-2*h1*h2))));
    Ct = ratQ.contentFree(Ct);   Zt = ratQ.contentFree(Zt);
    degS1 = max(size(Ct,1), size(Zt,1)) - 1;

    p = ratQ.trimN(ratQ.contentFree(ratQ.resultantVar(Ct, Zt, 1, nv)));
end

function A = S2_placeholder(nv)
    A = ratQ.monoN(zeros(1,nv), 0, nv);
end

function [v, s] = evalRel(p, x, y, z)
    v = 0;  s = 0;
    for j = 1:size(p,2)
        for k = 1:size(p,3)
            for l = 1:size(p,4)
                cf = p(1,j,k,l);
                if cf == 0, continue, end
                t = cf * x^(j-1) * y^(k-1) * z^(l-1);
                v = v + t;  s = s + abs(t);
            end
        end
    end
    if s == 0, s = 1; end
end
