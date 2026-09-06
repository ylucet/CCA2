function bridgeImplicitProbe()
% Produce the bridge face's IMPLICIT POLYNOMIAL p(x,y,z) = 0 exactly, over the integers, using
% ratQ.resultantVar -- the step .claude/bridgeEliminationSpike.m named as missing -- and then find
% out how close the computation runs to ratQ's 2^53 exactness ceiling.
%
% THE INSTANCE. q1(v) = |v|^2/2 and q2(v) = |v - a|^2 + c, with a = (pn,qn)/ad and c = cn/cd. The
% three equations in the dual variable s are
%
%   C = (ad*cd)(s1^2 + s2^2) - 4*cd*(pn*s1 + qn*s2) + 4*ad*cn        (the tie, q1* = q2*)
%   D = A*s2 + B,   A = 2*pn - ad*x,  B = (ad*y - 2*qn)*s1 + 2*qn*x - 2*pn*y
%   Z = s1^2 + s2^2 - 2*s1*x - 2*s2*y + 2*z                          (z = <s,x> - q1*(s))
%
% ONE RESULTANT, NOT TWO, and this is the method rather than an optimisation. Eliminating s2 from
% (C,D) and from (C,Z) and then s1 needs a Sylvester matrix of size 12, far past what cofactor
% expansion can afford. But D is LINEAR in s2, so it is SOLVED instead: s2 = -B/A. Substituting and
% clearing A^2 leaves two polynomials of degree 2 in s1, so a single 4x4 resultant finishes it. An
% equation linear in one of the eliminated variables should always be substituted, never
% resultanted.
%
% WHAT p IS AND IS NOT. It VANISHES on the face -- that is what a resultant guarantees -- but it is
% not minimal: a resultant carries the extraneous factors that come from clearing A^2 and from the
% roots of the system that are not on the bridge. Reducing it would need squarefree factorisation,
% which ratQ does not have.

    fprintf('=== the implicit polynomial, and whether it vanishes on the face ===\n');
    pn = 13;  qn = -7;  ad = 10;  cn = -1;  cd = 2;
    p = buildImplicit(pn, qn, ad, cn, cd);
    fprintf('p(x,y,z): size %s, %d nonzero terms, max |coeff| %.4g\n', ...
            mat2str(size(p)), nnz(p), max(abs(p(:))));
    fprintf('degrees: x %d, y %d, z %d\n', size(p,2)-1, size(p,3)-1, size(p,4)-1);

    a = [pn; qn]/ad;   c = cn/cd;
    ctr = 2*a;         R = 2*sqrt(a.'*a - c);
    fprintf('\npoints constructed ON the bridge (p must vanish)\n');
    worst = 0;
    for th = [0.7, 1.9, 3.4, 5.0]
        s0 = ctr + R*[cos(th); sin(th)];
        x  = 0.5*s0 + 0.5*(a + s0/2);
        z  = 2*(a.'*(x - 2*a)) + R*norm(x - 2*a) + 2*c;
        [v, sc] = evalRel(p, x(1), x(2), z);
        worst = max(worst, abs(v)/sc);
        fprintf('   th=%.2f  x=(%8.5f,%8.5f)  z=%10.6f   relative |p| = %.2e\n', ...
                th, x(1), x(2), z, abs(v)/sc);
    end
    fprintf('   worst relative |p| on the bridge: %.3e\n', worst);

    fprintf('\nthe same point with z perturbed (p must NOT vanish)\n');
    s0 = ctr + R*[cos(0.7); sin(0.7)];
    x  = 0.5*s0 + 0.5*(a + s0/2);
    z0 = 2*(a.'*(x - 2*a)) + R*norm(x - 2*a) + 2*c;
    for dz = [0.05, -0.05]
        [v, sc] = evalRel(p, x(1), x(2), z0 + dz);
        fprintf('   dz=%+.2f   relative |p| = %.2e\n', dz, abs(v)/sc);
    end

    % ---- THE EXACTNESS CEILING ------------------------------------------------------------------
    % TODO.md parks the integer-overflow item until an input a user would plausibly write actually
    % overflows. This elimination is the first thing in the project to come within a factor of a few
    % of 2^53 on a small instance, so the question is worth asking directly rather than assumed.
    fprintf('\n=== how close does the elimination run to 2^53 = %.4g ? ===\n', 2^53);
    fprintf('%-18s %-14s %s\n', 'a = (13,-7)/d', 'max |coeff|', 'outcome');
    fprintf('%s\n', repmat('-', 1, 56));
    for d = [4 6 8 10 12 14 16 20]
        try
            pp = buildImplicit(13, -7, d, -1, 2);
            fprintf('%-18d %-14.4g ok (%.1f%% of the ceiling)\n', ...
                    d, max(abs(pp(:))), 100*max(abs(pp(:)))/2^53);
        catch e
            fprintf('%-18d %-14s %s\n', d, '--', e.identifier);
        end
    end
end

% ------------------------------------------------------------------------------------------------
function p = buildImplicit(pn, qn, ad, cn, cd)
% The elimination, for q1 = |v|^2/2 and q2 = |v - (pn,qn)/ad|^2 + cn/cd. Variables (s1, x, y, z).
    nv = 4;
    S1 = mono(nv, [1 0 0 0]);   S2 = mono(nv, [2 0 0 0]);
    X  = mono(nv, [0 1 0 0]);   Y  = mono(nv, [0 0 1 0]);
    Zv = mono(nv, [0 0 0 1]);   one = mono(nv, [0 0 0 0]);

    A  = add(one*(2*pn), X*(-ad));                                   % 2*pn - ad*x
    B  = add(mul(S1, add(Y*ad, one*(-2*qn))), add(X*(2*qn), Y*(-2*pn)));
    A2 = mul(A, A);   B2 = mul(B, B);   AB = mul(A, B);

    % REDUCE C TO ITS PRIMITIVE PART FIRST. A conic equation has no scale, and in a 4x4 Sylvester
    % determinant an avoidable common factor is multiplied roughly four times over -- measured, a
    % factor of 4 left in here took the elimination from 2.7e15 to an overflow at 1.3e16.
    Cv  = ratQ.contentFree([ad*cd, -4*cd*pn, -4*cd*qn, 4*ad*cn]);
    Ca  = Cv(1);  Cb1 = Cv(2);  Cb2 = Cv(3);  Cc = Cv(4);
    Ct = add(add(mul(S2, A2)*Ca, B2*Ca), ...
             add(mul(S1, A2)*Cb1, add(AB*(-Cb2), A2*Cc)));
    Zt = add(add(mul(S2, A2), B2), ...
             add(mul(X, mul(S1, A2))*(-2), add(mul(Y, AB)*2, mul(Zv, A2)*2)));
    Ct = ratQ.contentFree(Ct);   Zt = ratQ.contentFree(Zt);

    p = ratQ.trimN(ratQ.contentFree(ratQ.resultantVar(Ct, Zt, 1, nv)));
end

function A = mono(nv, e)
    sz = e + 1;  if numel(sz) < 2, sz = [sz 1]; end
    A = zeros(sz);  idx = num2cell(e + 1);  A(idx{:}) = 1;
end
function C = mul(A, B), C = ratQ.polyMulN(A, B); end
function C = add(A, B), C = ratQ.polyAddN(A, B); end

function [v, s] = evalRel(p, x, y, z)
% the value, and the sum of |term| -- a residual means nothing next to coefficients of size 1e15
% unless it is read relative to the terms that produced it.
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
