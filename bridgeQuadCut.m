function [R, info] = bridgeQuadCut(FP, num, den)
% bridgeQuadCut  The curve where a BRIDGE face's value meets a QUADRATIC face's: the first edge the
%   conic subdivision cannot hold, made exactly and over the integers.
%
% objective: a bridge face carries the envelope's value -- measured strictly above the best
%   quadratic face at 22 of 60 points -- so where the two compete the mesh has to draw the boundary
%   between them. `envelopeMesh:bridgeCompetes` refuses there today. This is that boundary.
%
% [input]  FP  : the bridge face's implicit polynomial, ratQ's dense polyN in (x, y, z), degree 2
%                in z -- P(i+1,j+1,k+1) is the coefficient of x^i y^j z^k
%          num : the quadratic face's numerator, 1 x 10 in RatCon's weighted basis
%          den : its positive integer denominator
% [output] R   : integer polynomial in (x, y), as a dense polyN in TWO variables
%          info : .degree, .zdeg
%
% ==================================================================================================
% WHY SUBSTITUTION IS THE WHOLE CONSTRUCTION, AND WHY THE RESULT IS RATIONAL
%
% The bridge's value z is not a formula -- it is the branch `FB` of `FP(x, y, .) = 0`, an irrational
% function of x. The quadratic face's value q(x, y) = num . basis / den IS a formula. So the set
% where they agree is
%
%       { x : FP(x, y, q(x,y)) = 0 } ,
%
% and although one side of that equation is irrational as a VALUE, the equation itself is rational:
% substituting a rational function into an integer polynomial and clearing the denominator leaves an
% integer polynomial. With FP of degree 2 in z and q of degree 2 in (x, y),
%
%       R = FP_2 q_n^2 + FP_1 q_n den + FP_0 den^2 ,
%
% where FP_k is FP's coefficient of z^k -- a polynomial in (x, y) alone -- and q_n is num's
% quadratic. Degree 2 in z is what makes this three terms rather than a resultant.
%
% WHAT IT COSTS. R has degree up to 2*max(deg FP_k) + 4, well above two, so it is NOT a conic: this
% is precisely the edge the `Conic` subdivision level cannot store, and the reason the subdivision
% axis has to grow (DECISIONS.md 2026-09-19). R contains the boundary but need not equal it -- the
% branch matters, and the other branch of FP meets q on the same curve -- so a consumer still has to
% say WHICH side of R the face is on, exactly as a branch row does one degree down.

    if nargin ~= 3
        error('bridgeQuadCut:nargin', 'expects (FP, num, den).');
    end
    if ratQ.isBigPoly(FP)
        % THE BIG FORM (a polynomial past 2^53, kept in BigInteger -- bridgeFaceAlg). This cut is
        % exact, so its input is read exactly: refused on z-degree first, as a double FP is, and
        % otherwise brought back by ratQ.fromBig, which raises ratQ:overflow if it does not fit.
        szB = FP.sz;  szB(end+1:3) = 1;
        if szB(3) ~= 3
            error('bridgeQuadCut:zDegree', ...
                  ['FP must have degree 2 in z -- a bridge face''s value is one of two conjugate ' ...
                   'roots; got degree %d.'], szB(3) - 1);
        end
        FP = ratQ.fromBig(FP);
    end
    if ~isnumeric(FP) || ndims(FP) > 3
        error('bridgeQuadCut:shape', 'FP must be a polyN in at most three variables.');
    end
    num = ratQ.chk(double(num(:).'), 'face numerator');
    if numel(num) ~= 10
        error('bridgeQuadCut:shape', 'num must be 1 x 10; got %d entries.', numel(num));
    end
    den = ratQ.chk(double(den), 'face denominator');
    if ~isscalar(den) || den <= 0
        error('bridgeQuadCut:denominator', 'den must be a positive scalar; got %s.', mat2str(den));
    end
    if any(num(1:4) ~= 0)
        error('bridgeQuadCut:cubic', ...
              'num must be QUADRATIC -- its cubic entries are nonzero: %s.', mat2str(num(1:4)));
    end

    % ---- FP's coefficients of z^k, each a polynomial in (x, y) -------------------------------------
    zdeg = size(FP, 3) - 1;
    if zdeg ~= 2
        error('bridgeQuadCut:zDegree', ...
              ['FP must have degree 2 in z -- a bridge face''s value is one of two conjugate ' ...
               'roots; got degree %d.'], zdeg);
    end
    F0 = FP(:,:,1);  F1 = FP(:,:,2);  F2 = FP(:,:,3);

    % ---- the quadratic face's numerator as a polyN in (x, y) ---------------------------------------
    % RatCon's weighted basis: entries 5..10 are [xx xy yy x y 1] with the SAME weights `eval` uses.
    Q = zeros(3,3);
    Q(3,1) = num(5);   Q(2,2) = num(6);   Q(1,3) = num(7);
    Q(2,1) = num(8);   Q(1,2) = num(9);   Q(1,1) = num(10);
    D = zeros(1,1);  D(1,1) = den;

    % ---- R = F2 Q^2 + F1 Q D + F0 D^2 -------------------------------------------------------------
    R = ratQ.polyAddN(ratQ.polyMulN(F2, ratQ.polyMulN(Q, Q)), ...
                      ratQ.polyAddN(ratQ.polyMulN(F1, ratQ.polyMulN(Q, D)), ...
                                    ratQ.polyMulN(F0, ratQ.polyMulN(D, D))));
    R = ratQ.trimN(ratQ.contentFree(ratQ.chk(R, 'bridge-quadratic cut')));
    info = struct('degree', max(size(R)) - 1, 'zdeg', zdeg);
end
