function z = algFibreVanishes(P, X)
% algFibreVanishes  Is p(x,y,.) IDENTICALLY ZERO at each point -- a base point?
%
% objective: where every coefficient of the fibre p(x,y,z) in z vanishes, every z is a root and
%   algBranch returns NaN: a row of AlgCrv comparing that value is UNDEFINED there, not failed
%   (#108 residual: D-153's domain corner (0,0) was a hole). The zero test is algBranch's own: exact
%   for a polynomial kept in BigInteger, within degTol = 1e-12 of the terms' own size otherwise.
% [input]  P : algBranch's polynomial (dense polyN array in x, y, z, or a big form)
%          X : n x 2 points
% [output] z : n x 1 logical
    n = size(X,1);
    z = false(n,1);
    if ratQ.isBigPoly(P)
        Pd = ratQ.bigToDouble(P);
        sz = size(Pd);  sz(end+1:3) = 1;
        for t = 1:n
            X0 = java.math.BigDecimal(X(t,1));  Y0 = java.math.BigDecimal(X(t,2));
            allZero = true;
            for k = 1:sz(3)
                acc = java.math.BigDecimal.ZERO;
                for i = 1:sz(1)
                    for j = 1:sz(2)
                        v = P.c{sub2ind(P.sz, i, j, k)};
                        if isempty(v), continue, end
                        acc = acc.add(java.math.BigDecimal(v).multiply(X0.pow(i-1)).multiply(Y0.pow(j-1)));
                    end
                end
                if acc.signum() ~= 0, allZero = false; break, end
            end
            z(t) = allZero;
        end
        return
    end
    sz = size(P);  sz(end+1:3) = 1;
    for t = 1:n
        xs = X(t,1) .^ (0:sz(1)-1);  ys = (X(t,2) .^ (0:sz(2)-1)).';
        allZero = true;
        for k = 1:sz(3)
            c = xs * P(:,:,k) * ys;
            ca = abs(xs) * abs(P(:,:,k)) * abs(ys);
            if abs(c) > 1e-12 * ca, allZero = false; break, end
        end
        z(t) = allZero;
    end
end
