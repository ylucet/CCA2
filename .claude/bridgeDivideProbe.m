function bridgeDivideProbe()
% RUNG 1 after a wrong guess. `.claude/bridgeImplicitProbe.m` divided the eliminated polynomial p by
% the two factors the derivation suggested -- A and the leading s1^2 coefficient -- and NEITHER went
% in, 0 times each. So the extraneous content is not what it was guessed to be.
%
% Ask the question the other way round instead of guessing again: build the face's KNOWN quadric
% M explicitly, confirm it vanishes on the face, and then divide p by IT. Both vanish on the same
% 2-dimensional surface, so if M is irreducible it must divide p -- and the quotient is then the
% extraneous part, whatever it turns out to be.
%
% M comes from the closed form (z - <2a,x-2a> - 2c)^2 - R^2|x-2a|^2 with a = (13,-7)/10, c = -1/2
% and R^2 = 268/25, cleared to integers by 625:
%
%   M = (25z - 65x + 35y + 243)^2 - 268*((5x-13)^2 + (5y+7)^2)

    nv = 4;                                   % variables: s1 (singleton), x, y, z
    L  = lin(nv, [0 -65  35  25], 243);
    P1 = lin(nv, [0   5   0   0], -13);
    P2 = lin(nv, [0   0   5   0],   7);
    M  = ratQ.polyAddN(ratQ.polyMulN(L, L), ...
                       -268 * ratQ.polyAddN(ratQ.polyMulN(P1,P1), ratQ.polyMulN(P2,P2)));
    M  = ratQ.trimN(ratQ.contentFree(M));
    fprintf('M: degrees (x,y,z) = (%d,%d,%d), %d nonzero terms, max |coeff| %g\n', ...
            size(M,2)-1, size(M,3)-1, size(M,4)-1, nnz(M), max(abs(M(:))));

    a = [13; -7]/10;  c = -1/2;  ctr = 2*a;  R = 2*sqrt(a.'*a - c);
    fprintf('\nM on the face (it should vanish -- this is the sanity check on M itself)\n');
    for th = [0.7 1.9 3.4]
        [x, z] = bridgePoint(a, c, ctr, R, th);
        fprintf('   th=%.2f  relative |M| = %.2e\n', th, relAt(M, x(1), x(2), z));
    end

    fprintf('\nnow divide the eliminated p by M\n');
    p = buildImplicit(13, -7, 10, -1, 2);
    fprintf('   p degrees (x,y,z) = (%d,%d,%d)\n', size(p,2)-1, size(p,3)-1, size(p,4)-1);
    try
        q = ratQ.trimN(ratQ.polyDivExact(p, M, nv));
        fprintf('   M DIVIDES p. quotient degrees (x,y,z) = (%d,%d,%d), %d terms, max |coeff| %.4g\n', ...
                size(q,2)-1, size(q,3)-1, size(q,4)-1, nnz(q), max(abs(q(:))));
        fprintf('   -> the extraneous part is that quotient, and it is free of z: %d\n', ...
                size(q,4) == 1);
    catch e
        fprintf('   M does NOT divide p: %s\n', e.identifier);
        fprintf(['   That would mean p and M vanish on the same surface without sharing a factor,\n' ...
                 '   which cannot happen for an irreducible M -- so either M is reducible or the\n' ...
                 '   elimination is not producing what it is believed to.\n']);
    end
end

% ------------------------------------------------------------------------------------------------
function A = lin(nv, coef, k)
% the affine polynomial coef * v + k, as an nv-dimensional coefficient array
    A = 0;
    for d = 1:nv
        if coef(d) == 0, continue, end
        e = zeros(1, nv);  e(d) = 1;
        A = ratQ.polyAddN(A, ratQ.monoN(e, coef(d), nv));
    end
    A = ratQ.polyAddN(A, ratQ.monoN(zeros(1,nv), k, nv));
end

function [x, z] = bridgePoint(a, c, ctr, R, th)
    s0 = ctr + R*[cos(th); sin(th)];
    x  = 0.5*s0 + 0.5*(a + s0/2);
    z  = 2*(a.'*(x - 2*a)) + R*norm(x - 2*a) + 2*c;
end

function r = relAt(P, x, y, z)
    v = 0;  s = 0;
    for j = 1:size(P,2)
        for k = 1:size(P,3)
            for l = 1:size(P,4)
                cf = P(1,j,k,l);
                if cf == 0, continue, end
                t = cf * x^(j-1) * y^(k-1) * z^(l-1);
                v = v + t;  s = s + abs(t);
            end
        end
    end
    r = abs(v)/max(s, 1);
end

% ------------------------------------------------------------------------------------------------
function p = buildImplicit(pn, qn, ad, cn, cd)
% Identical to bridgeImplicitProbe's, repeated here so this probe stands alone.
    nv = 4;
    S1 = mono(nv, [1 0 0 0]);   S2 = mono(nv, [2 0 0 0]);
    X  = mono(nv, [0 1 0 0]);   Y  = mono(nv, [0 0 1 0]);
    Zv = mono(nv, [0 0 0 1]);   one = mono(nv, [0 0 0 0]);

    A  = add(one*(2*pn), X*(-ad));
    B  = add(mul(S1, add(Y*ad, one*(-2*qn))), add(X*(2*qn), Y*(-2*pn)));
    A2 = mul(A, A);   B2 = mul(B, B);   AB = mul(A, B);

    Cv  = ratQ.contentFree([ad*cd, -4*cd*pn, -4*cd*qn, 4*ad*cn]);
    Ca  = Cv(1);  Cb1 = Cv(2);  Cb2 = Cv(3);  Cc = Cv(4);
    Ct = add(add(mul(S2, A2)*Ca, B2*Ca), ...
             add(mul(S1, A2)*Cb1, add(AB*(-Cb2), A2*Cc)));
    Zt = add(add(mul(S2, A2), B2), ...
             add(mul(X, mul(S1, A2))*(-2), add(mul(Y, AB)*2, mul(Zv, A2)*2)));
    Ct = ratQ.contentFree(Ct);   Zt = ratQ.contentFree(Zt);

    p = ratQ.trimN(ratQ.contentFree(ratQ.resultantVar(Ct, Zt, 1, nv)));
end

function A = mono(nv, e), A = ratQ.monoN(e, 1, nv); end
function C = mul(A, B),   C = ratQ.polyMulN(A, B); end
function C = add(A, B),   C = ratQ.polyAddN(A, B); end
