function [s, info] = algCurveSide(FP, FB, num, den, X)
% algCurveSide  Which side of the boundary {bridge value = quadratic value} a point is on -- exactly,
%   for ONE named branch, without ever forming the irrational value.
%
% objective: the edge the conic subdivision cannot hold. `bridgeQuadCut` builds the rational curve
%   R = FP(x, y, q(x,y)) that contains that boundary, but R vanishes wherever q meets EITHER root of
%   FP, so its sign is the PRODUCT over both branches -- measured, the product agrees with R
%   everywhere the fibre is real while one branch agrees only 63% to 97% of the time. So R names the
%   PAIR, and the side of one branch is what a face actually needs. This decides it.
%
% [input]  FP  : the bridge face's polynomial, polyN in (x, y, z), degree 2 in z
%          FB  : its branch index -- 1-based from the BOTTOM when positive, from the TOP when
%                negative, exactly as `algBranch` reads it
%          num, den : the quadratic face, numerator in the weighted basis over a positive integer
%          X   : n x 2 RATIONAL points, as doubles (every double is a rational)
% [output] s   : n x 1 in {-1, 0, +1} -- the sign of q(x) - value_FB(x), so +1 means the QUADRATIC
%                face is larger there. NaN where the branch does not exist at that point.
%          info : .nreal, .below -- the fibre's real-root count and how many roots lie at or below q
%
% ==================================================================================================
% THE TEST IS A ROOT COUNT, AND THAT IS WHY IT IS EXACT
%
% The bridge's value at x is the FB-th real root of the integer polynomial `FP(x, y, .)`. The
% quadratic's value q(x) is RATIONAL. So
%
%       q > value_FB   <=>   at least FB of the fibre's real roots lie strictly below q,
%
% counting from the bottom, and `ratQ.countRootsIn` decides exactly that by Sturm's theorem -- a
% sign-variation count in an integer chain, with no root ever computed. Equality is the case
% `FP(x, y, q) = 0` with the count landing exactly at FB, which is `bridgeQuadCut`'s own curve, and
% it is reported as 0.
%
% This is the same move as `ratQ.pairLineSign` one degree up: the object is irrational, the NAME is
% not, and the predicate that uses it is a sign or a count over the integers.
%
% WHY THE POINT MUST BE RATIONAL. Clearing denominators is what makes the fibre an integer
% polynomial. Every double IS a rational, so a caller passing doubles is exact here; a caller that
% has an algebraic point must isolate it first, which is the degree-4 kernel and not this routine.

    if nargin ~= 5
        error('algCurveSide:nargin', 'expects (FP, FB, num, den, X).');
    end
    if ratQ.isBigPoly(FP)
        % the big form (bridgeFaceAlg): a side is a DECISION, so it is read exactly -- refused on
        % z-degree as a double FP is, otherwise brought back by ratQ.fromBig (raises if too big)
        szB = FP.sz;  szB(end+1:3) = 1;
        if szB(3) == 3, FP = ratQ.fromBig(FP); else, FP = []; end
    end
    if ~isnumeric(FP) || size(FP,3) ~= 3
        error('algCurveSide:zDegree', ...
              'FP must have degree 2 in z -- a bridge face''s value is one of two roots.');
    end
    if ~isscalar(FB) || FB == 0 || FB ~= fix(FB)
        error('algCurveSide:branch', ...
              'FB must be a nonzero integer branch index; got %s.', mat2str(FB));
    end
    if ~isnumeric(X) || size(X,2) ~= 2
        error('algCurveSide:points', 'X must be n x 2; got %s.', mat2str(size(X)));
    end
    den = ratQ.chk(double(den), 'face denominator');

    n = size(X,1);
    s = nan(n,1);
    info = struct('nreal', zeros(n,1), 'below', zeros(n,1));

    for i = 1:n
        [p, qn, ok] = fibreAndValue(FP, num, den, X(i,:));
        if ~ok, continue, end                          % the fibre degenerated: no branch here

        % MULTIPLICITY IS WHERE THE TWO CONVENTIONS MEET. `algBranch` sorts the fibre's real roots
        % WITH multiplicity -- its header says so, and a double root is two entries there -- while
        % Sturm counts DISTINCT roots, because that is what a vertex count wants. They disagree
        % exactly at a double root, which for a bridge face is the curve where its two branches
        % merge: the cone's apex in the small, the discriminant curve in general. The fibre is
        % quadratic, so the discriminant settles it outright.
        [a, bb, cc] = coeffs3(p);
        disc = ratQ.chkB(bb^2 - 4*a*cc, bb^2 + 4*abs(a)*abs(cc), 'fibre discriminant');
        if disc < 0, info.nreal(i) = 0;  continue, end % no real branch at this point

        if disc == 0
            % ONE root of multiplicity two: both branches ARE that root, so the side is the same
            % for either. r = -bb/(2a), and sign(q - r) = sign(2 a qn + bb) sign(a) over integers.
            info.nreal(i) = 2;                         % with multiplicity, as algBranch counts
            info.below(i) = NaN;
            t = ratQ.chkB(2*a*qn + bb, 2*abs(a)*abs(qn) + abs(bb), 'double root comparison');
            if t == 0, s(i) = 0; else, s(i) = sign(t) * sign(a); end
            continue
        end

        nreal = 2;                                     % disc > 0: two distinct real roots
        info.nreal(i) = nreal;
        b = FB;
        if b < 0, b = nreal + 1 + b; end               % NEGATIVE counts from the TOP, as algBranch
        if b < 1 || b > nreal, continue, end           % the branch does not exist at this point

        % The fibre is in the SCALED variable u = qd*z, so the comparison point is the INTEGER qn
        % and `countRootsIn`'s endpoint is exact -- which it would not be if q were a division.
        below = ratQ.countRootsIn(p, -inf, qn);        % roots in (-inf, qn], so qn itself counts
        info.below(i) = below;
        if polyAtInt(p, qn) == 0 && below == b
            s(i) = 0;                                  % q IS the b-th root: on the boundary
        elseif below >= b
            s(i) = 1;                                  % at least b roots at or below q
        else
            s(i) = -1;
        end
    end
end

% ==================================================================================================
function [p, qn, ok] = fibreAndValue(FP, num, den, x)
% The fibre `FP(x, y, .)` as an integer polynomial, in the SCALED variable u = qd*z, together with
% the integer qn for which u = qn means z = q(x). Scaling is what keeps the comparison exact: a
% rational endpoint would have to be a division, and `countRootsIn` wants an exact one.
    ok = false;  p = [];  qn = 0;
    % THE POINT MUST BE A SIMPLE RATIONAL, and the limit is measured rather than assumed. Clearing a
    % denominator 2^e through FP costs about |coef| (|x| 2^e)^degx 2^(e degy); with FP degree 4 in
    % each variable and coefficients of size 4 to 16 on the real fixtures, e = 4 is comfortably
    % inside 2^53 and e = 6 is over it. An arbitrary double IS a rational, but with a denominator
    % near 2^52 -- so it is refused BY NAME here rather than overflowing three calls down.
    [n1, d1] = simpleRational(x(1));
    [n2, d2] = simpleRational(x(2));

    % ---- the fibre's z^k coefficients, cleared by d1^I d2^J (positive, so the roots are unmoved) --
    I = size(FP,1) - 1;  J = size(FP,2) - 1;
    c = zeros(1,3);
    for k = 1:3
        acc = 0;
        for i = 0:I
            for j = 0:J
                a = FP(i+1, j+1, k);
                if a == 0, continue, end
                acc = acc + a * n1^i * d1^(I-i) * n2^j * d2^(J-j);
            end
        end
        c(k) = acc;
    end
    c = ratQ.chk(c, 'fibre coefficients');

    % ---- q(x) = qn/qd exactly, in the WEIGHTED basis (x^2 and y^2 carry 1/2) ----------------------
    qn = ratQ.chkB(num(5)*n1^2*d2^2 + 2*num(6)*n1*n2*d1*d2 + num(7)*n2^2*d1^2 ...
                  + 2*num(8)*n1*d1*d2^2 + 2*num(9)*n2*d2*d1^2 + 2*num(10)*d1^2*d2^2, ...
                  abs(num(5:10)) * abs([n1^2*d2^2; 2*n1*n2*d1*d2; n2^2*d1^2; ...
                                        2*n1*d1*d2^2; 2*n2*d2*d1^2; 2*d1^2*d2^2]), ...
                  'face numerator at the point');
    qd = ratQ.chk(2*den*d1^2*d2^2, 'face denominator at the point');

    % LOWEST TERMS FIRST. The scaling below squares qd, so an unreduced fraction is what puts this
    % over the ceiling: with a point of denominator 8 and a face denominator of its own, qd^2 alone
    % reached 1e17 before this line existed.
    g = gcd(abs(qn), abs(qd));
    if g > 1, qn = qn / g;  qd = qd / g; end

    % ---- u = qd z, so c2 z^2 + c1 z + c0 becomes c2 u^2 + c1 qd u + c0 qd^2 -----------------------
    p = [c(3), c(2)*qd, c(1)*qd^2];
    gp = gcdAll(p);
    if gp > 1, p = p / gp; end
    p = ratQ.polyTrim(ratQ.chk(p, 'scaled fibre'));
    if numel(p) < 2, return, end                      % constant in u: no branch to name
    ok = true;
end

% --------------------------------------------------------------------------------------------------
function [n, d] = simpleRational(v)
% A double as an exact fraction with a SMALL power-of-two denominator, or a named refusal.
    lim = 2^16;
    d = 1;
    while d <= lim && v*d ~= round(v*d)
        d = d * 2;
    end
    if d > lim
        error('algCurveSide:pointTooFine', ...
              ['the point %.17g is not a simple rational -- its exact denominator exceeds %d, and ' ...
               'clearing it through FP would leave ratQ''s 2^53 ceiling. Pass a point with a small ' ...
               'power-of-two denominator; a cell has many.'], v, lim);
    end
    n = ratQ.chk(round(v*d), 'point numerator');
end

% --------------------------------------------------------------------------------------------------
function g = gcdAll(v)
    v = abs(v(v ~= 0));
    if isempty(v), g = 0; return, end
    g = v(1);
    for k = 2:numel(v), g = gcd(g, v(k)); end
end

% --------------------------------------------------------------------------------------------------
function [a, b, c] = coeffs3(p)
% The scaled fibre as a u^2 + b u + c, whatever its stored length.
    p = p(:).';
    p = [zeros(1, 3 - numel(p)), p];
    a = p(1);  b = p(2);  c = p(3);
end

% --------------------------------------------------------------------------------------------------
function v = polyAtInt(p, t)
% p(t) for an integer polynomial at an integer, exactly.
    v = 0;  vb = 0;                    % vb bounds every partial sum: Horner on |p| at |t|
    for k = 1:numel(p)
        vb = vb * abs(t) + abs(p(k));
        v = ratQ.chkB(v * t + p(k), vb, 'fibre at the comparison point');
    end
end
