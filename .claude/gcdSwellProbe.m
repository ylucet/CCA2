function gcdSwellProbe()
% RUNG 1 after three failed edits at the same red. STOP editing polyGcdBiv and look at what it is
% actually being asked to do: the sizes of its real inputs, and where the swell happens.
%
% The question is not "which line overflows" -- that has been answered three times, each time in a
% different place. The question is whether a polynomial-remainder-sequence gcd can do this job AT
% ALL on inputs of this shape, or whether the algorithm is wrong for it.
    p = getP();
    nz = size(p,4);
    fprintf('the z-coefficients of p, which are what the gcd is asked about:\n');
    C = cell(1,nz);
    for k = 1:nz
        Ck = ratQ.trimN(reshape(p(1,:,:,k), size(p,2), size(p,3)));
        C{k} = ratQ.contentFree(Ck);
        fprintf('   c_%d: degrees (x,y) = (%d,%d), %d terms, max |coeff| %.4g (after contentFree %.4g)\n', ...
                k-1, size(Ck,1)-1, size(Ck,2)-1, nnz(Ck), max(abs(Ck(:))), max(abs(C{k}(:))));
    end

    fprintf('\ntrying the gcd of the top two, with the swell instrumented\n');
    a = ratQ.bivToCell(C{nz});  b = ratQ.bivToCell(C{nz-1});
    a = ratQ.bivPrimitive(a);   b = ratQ.bivPrimitive(b);
    fprintf('   after bivPrimitive: deg_x a = %d, deg_x b = %d\n', numel(a)-1, numel(b)-1);
    fprintf('   max |coeff| a = %.4g, b = %.4g\n', cellMax(a), cellMax(b));

    step = 0;
    while true
        step = step + 1;
        if numel(b) == 1, fprintf('   b became constant in x after %d steps\n', step-1); break, end
        try
            r = ratQ.bivPseudoRem(a, b);
        catch e
            fprintf('   step %d: %s -- %s\n', step, e.identifier, e.message(1:min(90,end)));
            break
        end
        if isempty(r), fprintf('   remainder vanished after %d steps\n', step); break, end
        a = b;  b = ratQ.bivPrimitive(r);
        fprintf('   step %d: deg_x = %d, max |coeff| = %.4g\n', step, numel(b)-1, cellMax(b));
        if step > 12, fprintf('   (stopping the trace)\n'); break, end
    end

    fprintf(['\nWhat this says. A PRS gcd multiplies by a leading COEFFICIENT at every elimination,\n' ...
             'and here those coefficients are polynomials in y with 12-digit integers. Reducing per\n' ...
             'step bounds the growth but does not remove it, because the leading coefficient itself\n' ...
             'is large. This is the classical intermediate expression swell, and the standard\n' ...
             'answer to it is a MODULAR gcd -- compute modulo several small primes, where every\n' ...
             'coefficient is bounded by the prime, and reconstruct -- not a better PRS.\n']);
end

function m = cellMax(c)
    m = 0;
    for i = 1:numel(c), m = max(m, max(abs(c{i}))); end
end

function p = getP()
    nv=4; pn=13; qn=-7; ad=10; cn=-1; cd=2;
    S1=ratQ.monoN([1 0 0 0],1,nv); S2=ratQ.monoN([2 0 0 0],1,nv);
    X=ratQ.monoN([0 1 0 0],1,nv);  Y=ratQ.monoN([0 0 1 0],1,nv);
    Zv=ratQ.monoN([0 0 0 1],1,nv); one=ratQ.monoN([0 0 0 0],1,nv);
    A=ratQ.polyAddN(one*(2*pn), X*(-ad));
    B=ratQ.polyAddN(ratQ.polyMulN(S1, ratQ.polyAddN(Y*ad, one*(-2*qn))), ...
                    ratQ.polyAddN(X*(2*qn), Y*(-2*pn)));
    A2=ratQ.polyMulN(A,A); B2=ratQ.polyMulN(B,B); AB=ratQ.polyMulN(A,B);
    Cv=ratQ.contentFree([ad*cd, -4*cd*pn, -4*cd*qn, 4*ad*cn]);
    Ct=ratQ.polyAddN(ratQ.polyAddN(ratQ.polyMulN(S2,A2)*Cv(1), B2*Cv(1)), ...
           ratQ.polyAddN(ratQ.polyMulN(S1,A2)*Cv(2), ratQ.polyAddN(AB*(-Cv(3)), A2*Cv(4))));
    Zt=ratQ.polyAddN(ratQ.polyAddN(ratQ.polyMulN(S2,A2), B2), ...
           ratQ.polyAddN(ratQ.polyMulN(X, ratQ.polyMulN(S1,A2))*(-2), ...
               ratQ.polyAddN(ratQ.polyMulN(Y,AB)*2, ratQ.polyMulN(Zv,A2)*2)));
    p=ratQ.trimN(ratQ.contentFree(ratQ.resultantVar(ratQ.contentFree(Ct), ratQ.contentFree(Zt), 1, nv)));
end
