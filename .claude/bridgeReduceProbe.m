function bridgeReduceProbe()
% Reduce the eliminated polynomial p to the face's MINIMAL polynomial, using only ratQ -- no
% knowledge of the closed form. This is what `.claude/bridgeDivideProbe.m` established was possible
% and what `ratQ.polyGcdBiv` was built to do.
%
% THE REDUCTION. p vanishes on the face but carries an extraneous factor that is FREE OF z:
% measured, p has degrees (x,y,z) = (10,6,2) and dividing by the closed form's quadric leaves a
% quotient of (8,4,0) with no z in it at all. So the extraneous part is exactly the CONTENT of p
% with respect to z -- the gcd of its z-coefficients as polynomials in (x,y) -- and removing it is
% a division by that gcd.
%
% Nothing below uses the closed form to FIND the answer; the closed form is used only at the end,
% to check it.

    fprintf('=== reducing p with no knowledge of the answer ===\n');
    p = buildImplicit(13, -7, 10, -1, 2);
    fprintf('p:  degrees (x,y,z) = (%d,%d,%d), %d terms, max |coeff| %.4g\n', ...
            size(p,2)-1, size(p,3)-1, size(p,4)-1, nnz(p), max(abs(p(:))));

    % the z-coefficients of p, each a bivariate polynomial in (x,y)
    nz = size(p,4);
    C  = cell(1, nz);
    for k = 1:nz
        C{k} = ratQ.trimN(reshape(p(1,:,:,k), size(p,2), size(p,3)));
    end

    % THE CONTENT, BY A SHORTCUT WITH A CHECK RATHER THAN BY A GCD.
    %
    % ratQ.polyGcdBiv is the general answer and it OVERFLOWS here, measured: p's coefficients are
    % already ~2.7e12, and polyPseudoRem's lc^(da-db+1) scaling squares them to 3.3e16, past 2^53.
    % The primitive-part discipline reduces BETWEEN steps and cannot help WITHIN the scaling.
    %
    % But the answer is tiny -- the reduced polynomial's coefficients top out at 486 -- so the size
    % is in the path, not the destination. The structure gives a way round: if the face's leading
    % z-coefficient is a CONSTANT, then c_top is exactly that constant times the content, so the
    % content is contentFree(c_top). That is a GUESS, so it is VERIFIED below -- the candidate must
    % divide every other z-coefficient exactly, and polyDivExact raises if it does not.
    g = ratQ.trimN(ratQ.contentFree(C{nz}));
    fprintf('candidate content (leading z-coefficient): (x,y) = (%d,%d), %d terms\n', ...
            size(g,1)-1, size(g,2)-1, nnz(g));
    for k = 1:nz-1
        ratQ.polyDivExact(C{k}, g, 2);            % raises unless the candidate is right
    end
    fprintf('   verified: it divides every z-coefficient exactly\n');

    % divide it out, coefficient by coefficient
    M = zeros(1, size(p,2), size(p,3), nz);
    for k = 1:nz
        q = ratQ.polyDivExact(C{k}, g, 2);
        M(1, 1:size(q,1), 1:size(q,2), k) = q;
    end
    M = ratQ.trimN(ratQ.contentFree(M));
    fprintf('reduced: degrees (x,y,z) = (%d,%d,%d), %d terms, max |coeff| %g\n', ...
            size(M,2)-1, size(M,3)-1, size(M,4)-1, nnz(M), max(abs(M(:))));

    % ---- now, and only now, compare with the closed form ----------------------------------------
    fprintf('\n=== checking against the closed form ===\n');
    a = [13; -7]/10;  c = -1/2;  ctr = 2*a;  R = 2*sqrt(a.'*a - c);
    for th = [0.7 1.9 3.4 5.0]
        [x, z] = bridgePoint(a, c, ctr, R, th);
        fprintf('   th=%.2f  relative |M| on the face = %.2e\n', th, relAt(M, x(1), x(2), z));
    end
    [x0, z0] = bridgePoint(a, c, ctr, R, 0.7);
    fprintf('   z + 0.05: relative |M| = %.2e   (must NOT vanish)\n', ...
            relAt(M, x0(1), x0(2), z0 + 0.05));

    % and it must agree with the quadric derived by hand, up to content
    L  = lin(4, [0 -65  35  25], 243);
    P1 = lin(4, [0   5   0   0], -13);
    P2 = lin(4, [0   0   5   0],   7);
    Mc = ratQ.trimN(ratQ.contentFree(ratQ.polyAddN(ratQ.polyMulN(L,L), ...
             -268*ratQ.polyAddN(ratQ.polyMulN(P1,P1), ratQ.polyMulN(P2,P2)))));
    same = isequal(size(M), size(Mc)) && isequal(M, Mc);
    fprintf('\n   identical to the hand-derived quadric (up to content): %d\n', same);
    if ~same && isequal(size(M), size(Mc))
        fprintf('   max |difference| = %g\n', max(abs(M(:) - Mc(:))));
    end
end

% ------------------------------------------------------------------------------------------------
function A = lin(nv, coef, k)
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

function p = buildImplicit(pn, qn, ad, cn, cd)
    nv = 4;
    S1 = ratQ.monoN([1 0 0 0],1,nv);   S2 = ratQ.monoN([2 0 0 0],1,nv);
    X  = ratQ.monoN([0 1 0 0],1,nv);   Y  = ratQ.monoN([0 0 1 0],1,nv);
    Zv = ratQ.monoN([0 0 0 1],1,nv);   one = ratQ.monoN([0 0 0 0],1,nv);
    ad_ = ad;
    A  = ratQ.polyAddN(one*(2*pn), X*(-ad_));
    B  = ratQ.polyAddN(ratQ.polyMulN(S1, ratQ.polyAddN(Y*ad_, one*(-2*qn))), ...
                       ratQ.polyAddN(X*(2*qn), Y*(-2*pn)));
    A2 = ratQ.polyMulN(A,A);  B2 = ratQ.polyMulN(B,B);  AB = ratQ.polyMulN(A,B);
    Cv = ratQ.contentFree([ad*cd, -4*cd*pn, -4*cd*qn, 4*ad*cn]);
    Ct = ratQ.polyAddN(ratQ.polyAddN(ratQ.polyMulN(S2,A2)*Cv(1), B2*Cv(1)), ...
             ratQ.polyAddN(ratQ.polyMulN(S1,A2)*Cv(2), ...
                 ratQ.polyAddN(AB*(-Cv(3)), A2*Cv(4))));
    Zt = ratQ.polyAddN(ratQ.polyAddN(ratQ.polyMulN(S2,A2), B2), ...
             ratQ.polyAddN(ratQ.polyMulN(X, ratQ.polyMulN(S1,A2))*(-2), ...
                 ratQ.polyAddN(ratQ.polyMulN(Y,AB)*2, ratQ.polyMulN(Zv,A2)*2)));
    Ct = ratQ.contentFree(Ct);   Zt = ratQ.contentFree(Zt);
    p  = ratQ.trimN(ratQ.contentFree(ratQ.resultantVar(Ct, Zt, 1, nv)));
end
