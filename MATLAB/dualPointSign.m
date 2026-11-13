function sg = dualPointSign(dp, Q)
% dualPointSign  The EXACT sign of an integer polynomial at a named dual point.
%
% objective: decide sign(Q(s0)) for the point `dualPoint` names, with no tolerance -- the one
%   predicate a triangle face's rows and comparisons reduce to (DECISIONS.md 2026-10-30 (fifth)).
%
% [input]  dp : as dualPoint returns it
%          Q  : integer array, Q(i+1,j+1) the coefficient of s1^i s2^j (ratQ's polyN convention)
% [output] sg : -1, 0 or +1
%
% THE REDUCTION. In dualPoint's sheared coordinates s1 = u + k v, s2 = v, Q becomes
% sum_j Q_j(u) v^j, and on the common roots v = -S0(u)/S1(u). Clearing S1^m (m the degree in v),
%
%       N(u) = sum_j Q_j(u) (-S0(u))^j S1(u)^(m-j),     sign Q(s0) = sign N(u0) * sign S1(u0)^m ,
%
% two signs of integer polynomials at the root u0 of R isolated in dp.iv -- ratQ.signAt's method, which
% decides a zero by a gcd and a nonzero by refining, both exactly. S1(u0) ~= 0 was decided when the
% name was made.
    Q = ratQ.chk(double(Q), 'polynomial');
    if ~ismatrix(Q)
        error('dualPointSign:shape', 'Q must be a 2-D coefficient array in (s1, s2).');
    end
    if all(Q(:) == 0), sg = 0; return, end
    k = dp.k;
    % Q in (u, v): s1^i s2^j = (u + k v)^i v^j = sum_t C(i,t) k^t u^(i-t) v^(t+j)
    [ni, nj] = size(Q);
    m = (ni - 1) + (nj - 1);
    Quv = zeros(ni, m + 1);                      % Quv(a+1, b+1): coefficient of u^a v^b
    Bnd = zeros(ni, m + 1);                      % ... and the sum of its terms' magnitudes
    for i = 0:ni-1
        for j = 0:nj-1
            if Q(i+1,j+1) == 0, continue, end
            for t = 0:i
                term = Q(i+1,j+1) * nchoosek(i,t) * k^t;
                Quv(i-t+1, t+j+1) = Quv(i-t+1, t+j+1) + term;
                Bnd(i-t+1, t+j+1) = Bnd(i-t+1, t+j+1) + abs(term);
            end
        end
    end
    Quv = ratQ.chkB(Quv, Bnd, 'sheared polynomial');
    while size(Quv,2) > 1 && all(Quv(:,end) == 0), Quv(:,end) = []; end
    m = size(Quv,2) - 1;

    % IN BigInteger THROUGHOUT. N's degree is m*deg(S) + deg(Q), and with a random cubic Q the
    % double path's gcd overflowed 2^53 at once (measured in dualPointTest); the name is exact only
    % if the arithmetic is.
    N = {};
    negS0 = ratQ.bigRow(-dp.S0);  S1 = ratQ.bigRow(dp.S1);
    for j = 0:m
        Qj = ratQ.polyTrim(fliplr(Quv(:,j+1).'));   % highest degree first
        if isempty(Qj) || all(Qj == 0), continue, end
        N = badd(N, bmul(bmul(ratQ.bigRow(Qj), bpow(negS0, j)), bpow(S1, m - j)));
    end
    sN = isolatedRootSign(N, dp.R, dp.iv);
    if sN == 0 || mod(m, 2) == 0, sg = sN; return, end
    sg = sN * isolatedRootSign(S1, dp.R, dp.iv);
end

function C = bmul(A, B)
    if isempty(A) || isempty(B), C = {}; return, end
    C = repmat({java.math.BigInteger.ZERO}, 1, numel(A) + numel(B) - 1);
    for i = 1:numel(A)
        for j = 1:numel(B)
            C{i+j-1} = C{i+j-1}.add(A{i}.multiply(B{j}));
        end
    end
end

function C = badd(A, B)
    n = max(numel(A), numel(B));
    z = java.math.BigInteger.ZERO;
    A = [repmat({z}, 1, n - numel(A)), A];  B = [repmat({z}, 1, n - numel(B)), B];
    C = cellfun(@(a, b) a.add(b), A, B, 'UniformOutput', false);
end

function P = bpow(A, n)
    P = {java.math.BigInteger.ONE};
    for t = 1:n, P = bmul(P, A); end
end
