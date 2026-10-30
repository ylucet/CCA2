function dp = dualPoint(c1, c2, sApprox)
% dualPoint  The EXACT name of one real common root of two integer conics -- the dual point s0 of a
%   triangle face of the biconjugate, which is a degree-4 algebraic point with no closed form.
%
% objective: name s0 once, so that every later question about it -- the sign of a polynomial at s0,
%   which is all a triangle face's rows and comparisons ever ask -- is decided exactly by
%   `dualPointSign`, and its value is read in doubles from a cache. DECISIONS.md 2026-10-30 (fifth).
%
% [input]  c1, c2  : 1 x 6 conics [a b c d e f] for a s1^2 + b s1 s2 + c s2^2 + d s1 + e s2 + f = 0,
%                    integer-valued up to scale (ratQ.conic canonicalises them)
%          sApprox : 2 x 1 (or 1 x 2) double approximation of the root to name
% [output] dp : struct
%                 .c1 .c2  the canonical integer conics
%                 .k       the integer shear: s1 = u + k v, s2 = v
%                 .R       the eliminant in u, integer, highest degree first, SQUAREFREE
%                 .iv      [lo hi] dyadic, neither end a root of R, exactly one root u0 of R inside
%                 .S1 .S0  integer polynomials in u with v = -S0(u)/S1(u) on the common roots,
%                          and S1(u0) ~= 0 (decided, not sampled)
%                 .s       2 x 1 double s0 -- the CACHE of the name above
%
% ==================================================================================================
% WHY THE POINT AND NOT A BRANCH OF z. A triangle face's value <s0,x> - g(s0) could be stored like a
% bridge's, as a polynomial FP(x,y,z) of degree 4 in z plus an index among its real roots. That index
% is not constant on the triangle: the four conjugate values are AFFINE in x, two real ones cross
% along a line, and the index of ours changes there. s0 does not move with x, so it is named once.
%
% THE NAME. Shear so that the second variable is determined on the common roots, eliminate it, and
% isolate one root of the eliminant:
%
%   * in (u, v) with s1 = u + k v, each conic is A2 v^2 + A1(u) v + A0(u); the resultant in v of
%     two such is R(u) = (A2 B0 - B2 A0)^2 - (A2 B1 - B2 A1)(A1 B0 - B1 A0), of degree <= 4;
%   * B2*C1 - A2*C2 has no v^2 term: S1(u) v + S0(u) with S1 = B2 A1 - A2 B1, S0 = B2 A0 - A2 B0,
%     so v = -S0/S1 wherever S1 is nonzero -- one v per root u0, which is what makes (R, iv) a
%     NAME of a POINT rather than of a vertical line through two;
%   * a shear k is accepted only when S1(u0) ~= 0, decided by isolatedRootSign. Two common roots on one
%     vertical line make S1 vanish there; a different k separates them, and at most finitely many
%     k fail, so k = 0..12 is tried in turn.
%
% THE FORMULA FOR R IS THE RESULTANT ONLY WHEN A2 OR B2 IS NONZERO; with A2 = 0 it is B2 times the
% resultant, a nonzero constant, which does not move a root. Both zero (two curves linear in v after
% the shear) is skipped to the next k, which makes one of them nonzero unless both conics have no
% quadratic part -- two lines, whose crossing is rational and needs no name; refused by name.

    A = ratQ.conic(c1);  B = ratQ.conic(c2);
    s = double(sApprox(:));
    if numel(s) ~= 2 || any(~isfinite(s))
        error('dualPoint:approx', 'sApprox must be a finite 2-vector; got %s.', mat2str(sApprox));
    end
    if all(A(1:3) == 0) && all(B(1:3) == 0)
        error('dualPoint:twoLines', ['both curves are lines, so their crossing is RATIONAL and ' ...
              'needs no algebraic name: %s and %s.'], mat2str(A), mat2str(B));
    end

    for k = 0:12
        [A2, A1, A0] = inV(A, k);  [B2, B1, B0] = inV(B, k);
        if A2 == 0 && B2 == 0, continue, end
        R = ratQ.polyTrim(padd(pmul(padd(A2*B0, -B2*A0), padd(A2*B0, -B2*A0)), ...
                               -pmul(padd(A2*B1, -B2*A1), padd(pmul(A1, B0), -pmul(B1, A0)))));
        if isempty(R) || all(R == 0)
            error('dualPoint:commonComponent', ['the conics %s and %s share a component, so ' ...
                  'they meet in a curve and not in a point.'], mat2str(A), mat2str(B));
        end
        R = ratQ.polyPrimitive(ratQ.polySquarefree(R));
        S1 = ratQ.polyTrim(padd(B2*A1, -A2*B1));
        S0 = ratQ.polyTrim(padd(B2*A0, -A2*B0));
        if isempty(S1) || all(S1 == 0), continue, end

        u = s(1) - k*s(2);
        iv = isolateNear(R, u);
        if isempty(iv)
            error('dualPoint:noMatch', ['no real root of the eliminant of %s and %s lies near ' ...
                  'u = %.17g; the approximation is not a common point of the two conics.'], ...
                  mat2str(A), mat2str(B), u);
        end
        if isolatedRootSign(S1, R, iv) == 0, continue, end     % two roots on one vertical line
        u0 = narrow(R, iv);
        v0 = -polyval(S0, u0) / polyval(S1, u0);
        s0 = [u0 + k*v0; v0];
        if norm(s0 - s) > 1e-6 * max(1, norm(s))
            error('dualPoint:noMatch', ['no common root of %s and %s lies near %s; the nearest ' ...
                  'is %s.'], mat2str(A), mat2str(B), mat2str(s.', 17), mat2str(s0.', 17));
        end
        dp = struct('c1', A, 'c2', B, 'k', k, 'R', R, 'iv', iv, 'S1', S1, 'S0', S0, 's', s0);
        return
    end
    error('dualPoint:noShear', ['no shear k = 0..12 gives the conics %s and %s a second ' ...
          'coordinate determined at the root near %s.'], mat2str(A), mat2str(B), mat2str(s.'));
end

function [a2, a1, a0] = inV(C, k)
% The conic in (u, v), s1 = u + k v, s2 = v, as a2 v^2 + a1(u) v + a0(u), polynomials highest first.
    a = C(1); b = C(2); c = C(3); d = C(4); e = C(5); f = C(6);
    a2 = ratQ.chkB(a*k^2 + b*k + c, abs(a)*k^2 + abs(b)*k + abs(c), 'sheared conic');
    a1 = ratQ.chkB([2*a*k + b, d*k + e], [2*abs(a)*k + abs(b), abs(d)*k + abs(e)], 'sheared conic');
    a0 = [a, d, f];
end

function iv = isolateNear(R, u)
% A dyadic interval around u holding EXACTLY one root of R, endpoints not roots -- found by
% widening from double resolution, every count exact (Sturm in BigInteger). Empty when no root is
% within 2^4 * max(1,|u|) of u. Which root it isolates is checked by the caller against sApprox.
    Rb = ratQ.bigRow(R);  B = struct('sz', [1 numel(Rb)], 'c', {Rb(:)});
    sc = max(1, abs(u));
    for e = -48:4
        d = sc * 2^e;
        lo = u - d;  hi = u + d;
        if rsign(Rb, lo) == 0 || rsign(Rb, hi) == 0, continue, end
        n = ratQ.countRootsInBig(B, lo, hi);
        if n == 1, iv = [lo hi]; return, end
        if n > 1, break, end
    end
    iv = zeros(0,2);
end

function u0 = narrow(R, iv)
% The root in doubles, by bisecting a COPY of the isolating interval -- the cache's accuracy, not a
% decision: the stored interval is the one isolateNear returned.
    Rb = ratQ.bigRow(R);  B = struct('sz', [1 numel(Rb)], 'c', {Rb(:)});
    lo = iv(1);  hi = iv(2);
    for it = 1:200
        mid = (lo + hi) / 2;
        if mid <= lo || mid >= hi, break, end
        if rsign(Rb, mid) == 0, u0 = mid; return, end          % a RATIONAL root, met exactly
        if ratQ.countRootsInBig(B, lo, mid) == 1, hi = mid; else, lo = mid; end
    end
    u0 = (lo + hi) / 2;
end

function s = rsign(P, x)
% sign of the big polynomial P at the dyadic x, exactly
    [xn, xd] = ratQ.dyadicParts(x);
    n = numel(P) - 1;  acc = java.math.BigInteger.ZERO;
    for j = 1:numel(P)
        acc = acc.add(P{j}.multiply(xn.pow(n-j+1)).multiply(xd.pow(j-1)));
    end
    s = acc.signum();
end

function c = pmul(a, b)
    c = ratQ.chkB(conv(a, b), conv(abs(a), abs(b)), 'polynomial product');
end

function c = padd(a, b)
    n = max(numel(a), numel(b));
    a = [zeros(1, n - numel(a)), a];  b = [zeros(1, n - numel(b)), b];
    c = ratQ.chkB(a + b, abs(a) + abs(b), 'polynomial sum');
end
