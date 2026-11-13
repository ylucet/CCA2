function s = isolatedRootSign(N, R, iv)
% isolatedRootSign  The EXACT sign of an integer polynomial at an isolated real root of another.
%
% [input]  N, R : integer polynomials, highest degree first, as double rows or as cell rows of
%                 java.math.BigInteger; R squarefree
%          iv   : [lo hi] dyadic with exactly one root of R in (lo, hi] -- isolateRoots' convention,
%                 so hi may BE the root and lo may be another one; both are handled
% [output] s    : -1, 0 or +1, the sign of N at that root
%
% ratQ.signAt's two steps -- a ZERO is a common root of gcd(R, N) in iv, decided by Sturm; otherwise
% refine until N has no root left in the interval and read its sign there -- in BigInteger, after
% reducing N modulo R by a sign-true pseudo-remainder. ratQ.signAt runs the same steps in doubles
% and overflowed 2^53 on the first random cubic dualPointTest tried (2026-10-30).
    if ~iscell(N), N = ratQ.bigRow(N); end
    if ~iscell(R), R = ratQ.bigRow(R); end
    N = strip(N);
    if isempty(N), s = 0; return, end
    Rb = strip(R);
    lo = iv(1);  hi = iv(2);
    if bsign(Rb, hi) == 0, s = bsign(N, hi); return, end      % the root IS hi, a rational
    while bsign(Rb, lo) == 0                                  % lo is ANOTHER root: move it up
        mid = (lo + hi) / 2;
        if bsign(Rb, mid) == 0, s = bsign(N, mid); return, end
        if ratQ.countRootsInBig(asBig(Rb), mid, hi) == 1, lo = mid; else, hi = mid; end
    end
    iv = [lo hi];
    r = strip(ratQ.bigPrem(N, Rb));
    if isempty(r), s = 0; return, end
    if numel(r) == 1, s = r{1}.signum(); return, end
    g = bgcd(Rb, r);
    if numel(g) > 1 && ratQ.countRootsInBig(asBig(g), iv(1), iv(2)) > 0, s = 0; return, end
    for it = 1:400
        if countSafe(r, lo, hi) == 0
            s = bsign(r, (lo + hi) / 2);  return
        end
        mid = (lo + hi) / 2;
        if bsign(Rb, mid) == 0, s = bsign(r, mid); return, end   % a RATIONAL root, met exactly
        if ratQ.countRootsInBig(asBig(Rb), lo, mid) == 1, hi = mid; else, lo = mid; end
    end
    error('isolatedRootSign:noConverge', ['400 bisections did not separate the root from the ' ...
          'zeros of the reduced polynomial; with the zero test above that is a defect.']);
end

function n = countSafe(p, lo, hi)
% Roots of p in (lo, hi]; an endpoint that is itself a root counts as "not yet separated".
    try
        n = ratQ.countRootsInBig(asBig(p), lo, hi);
    catch ME
        if ~strcmp(ME.identifier, 'ratQ:sturmEndpointIsRoot'), rethrow(ME), end
        n = inf;
    end
end

function s = bsign(P, x)
% sign of the big polynomial P at the dyadic x, exactly
    [xn, xd] = ratQ.dyadicParts(x);
    n = numel(P) - 1;  acc = java.math.BigInteger.ZERO;
    for j = 1:numel(P)
        acc = acc.add(P{j}.multiply(xn.pow(n-j+1)).multiply(xd.pow(j-1)));
    end
    s = acc.signum();
end

function g = bgcd(a, b)
% primitive gcd of two big polynomials (Euclid on pseudo-remainders, primitive part each round)
    a = ratQ.bigPrimitive(a);  b = ratQ.bigPrimitive(b);
    while true
        r = strip(ratQ.bigPrem(a, b));
        if isempty(r), g = b; return, end
        if numel(r) == 1, g = {java.math.BigInteger.ONE}; return, end
        a = b;  b = ratQ.bigPrimitive(r);
    end
end

function B = asBig(P)
% a cell row of BigInteger in the big-form struct ratQ's Sturm routines take
    B = struct('sz', [1 numel(P)], 'c', {P(:)});
end

function P = strip(P)
    while ~isempty(P) && P{1}.signum() == 0, P(1) = []; end
end
