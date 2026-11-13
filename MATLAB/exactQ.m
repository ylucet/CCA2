classdef exactQ
% exactQ  Exact arithmetic in the MULTIQUADRATIC field Q(sqrt(p1), ..., sqrt(pk)): a value is a
% rational combination of sqrt(m) over squarefree integers m.
%
% WHY THIS TYPE AND NOT ANOTHER. Removing the Symbolic Toolbox from CCA2's compute path is a
% REPRESENTATION change -- `symbolicFunction.f` is a `sym` and `region.vx/vy` are `sym`, so every
% `subs`/`isAlways` is a consequence of the data type. The number type therefore has to be chosen
% first, and two candidates are already refused by measurement (DECISIONS.md):
%
%   * DOUBLES -- one ULP made a shared facet invisible to `region.merge`, and Step 3's cell count
%     then grew without bound. The failure is silent and it is quadratic in cost.
%   * RATIONAL SNAPPING -- bounding the VERTEX denominators does not bound the downstream ones,
%     because the conjugate is a rational function of those coordinates; a few squarings carried
%     1e5 to 1e25 and the run hung.
%
% And rationals ALONE do not suffice: measured 2026-08-19 (T8, DECISIONS.md), A.5's split foot is
% irrational and no rational split reduces the 3-convex-edge case to something conjugable directly.
%
% WHY MULTIQUADRATIC AND NOT ONE EXTENSION. This type carried ONE quadratic extension until
% 2026-08-20 and RAISED when two were mixed, on the argument that silently building a tower is how
% an exact type turns back into a symbolic engine. The rule was right; the field was too small.
% Measured (`.claude/t1_multiquadratic_example.md`): the A.5 split of the single triangle
% conv{(5/2,3/2),(0,0),(1/2,1)} produces a sub-triangle whose vertex is
%
%       ( sqrt(30)/12 - sqrt(15)/6 + 5/4 ,  sqrt(30)/20 - sqrt(15)/10 + 3/4 )
%
% -- ONE number needing two extensions, so no caller can route around it by keeping cells apart --
% and the neighbouring triangle of the same quadrilateral brings sqrt(5), which Step 3 then
% subtracts from a sqrt(15) cell.
%
% What makes the generalisation the right one rather than a slide into a general number field:
% square roots of SQUAREFREE integers are closed under multiplication up to a rational factor
% (sqrt(15)*sqrt(30) = 15*sqrt(2)), so the family the algorithm generates is exactly
% Q(sqrt(p1),...,sqrt(pk)) over the primes it actually meets. Both properties the pipeline depends
% on survive:
%
%   * ZERO-TESTING STAYS EXACT AND CHEAP. Those sqrt(m) are linearly INDEPENDENT over Q, so a value
%     is zero exactly when every coefficient is zero -- which is what `region`'s facet and
%     redundancy tests ask, and a wrong "yes" there merges two regions that do not share a facet.
%   * SIGN IS DECIDED, NOT ESTIMATED. `signExact` writes x = a + b*sqrt(p) with a, b in the field
%     of the remaining primes: if a and b share a sign that is the answer, otherwise compare a^2
%     against b^2*p in that smaller field. One prime leaves per step, so it terminates in the
%     rationals, with no floating point at all. `sign` runs a CERTIFIED floating-point screen
%     first and falls back to it -- see that method for why, and for why the answer is still
%     exact.
%
% THE INVARIANT, and everything below depends on it: `m` is strictly increasing, every entry is a
% squarefree positive integer (1 allowed, and it carries the rational part), and every coefficient
% cn/cd is nonzero, in lowest terms, with cd > 0. Zero is the empty support. That form is CANONICAL
% -- two values are equal exactly when their three vectors are.
%
% RATIONALS are int64 numerator/denominator. Overflow raises rather than wrapping -- a wrong answer
% that looks exact is the one outcome worse than a slow one. Sign squares its operands, so a value
% carrying several extensions costs more coefficient growth than one carrying none; that is a
% signal about the pipeline, not a case to round away.
%
% `fromDouble` REFUSES what it cannot represent exactly rather than rounding: converting at the
% boundary is the caller's job, and a type that quietly approximates in its constructor has no
% exactness to offer.

    properties (SetAccess = immutable)
        m  int64 = zeros(1,0,'int64')   % squarefree radicands, strictly increasing (1 = rational)
        cn int64 = zeros(1,0,'int64')   % coefficient numerators, all nonzero
        cd int64 = zeros(1,0,'int64')   % coefficient denominators, all positive
    end

    methods
        function o = exactQ(varargin)
        % exactQ()                      zero
        % exactQ(x)                     from a double or an integer, exactly when it is a dyadic
        %                               or small rational (see fromDouble)
        % exactQ(n, dd)                 the rational n/dd
        % exactQ(m, cn, cd)             sum_i (cn(i)/cd(i)) * sqrt(m(i)), any m >= 0
        % exactQ(an, ad, bn, bd, d)     (an/ad) + (bn/bd)*sqrt(d)   -- the single-extension form
            switch nargin
                case 0
                    return
                case 1
                    a = varargin{1};
                    if isa(a, 'exactQ'), o = a; return, end
                    [n, dd] = exactQ.fromDouble(a);
                    [mm, nn, ddv] = exactQ.canon(int64(1), n, dd);
                case 2
                    [mm, nn, ddv] = exactQ.canon(int64(1), int64(varargin{1}), int64(varargin{2}));
                case 3
                    [mm, nn, ddv] = exactQ.canon(int64(varargin{1}(:).'), ...
                                                 int64(varargin{2}(:).'), int64(varargin{3}(:).'));
                case 5
                    [mm, nn, ddv] = exactQ.canon([int64(1), int64(varargin{5})], ...
                                                 [int64(varargin{1}), int64(varargin{3})], ...
                                                 [int64(varargin{2}), int64(varargin{4})]);
                otherwise
                    error('exactQ:args', 'exactQ takes 0, 1, 2, 3 or 5 arguments.');
            end
            o.m = mm; o.cn = nn; o.cd = ddv;
        end

        % ---- predicates ---------------------------------------------------------------
        function t = isZero(o),     t = isempty(o.m); end
        function t = isRational(o), t = isempty(o.m) || (numel(o.m) == 1 && o.m(1) == 1); end

        function c = coeffOf(o, mm)
        % The coefficient of sqrt(mm), as a rational exactQ -- the observable way to read a value
        % apart, so tests do not depend on the storage layout.
            [k, ms] = exactQ.squarefree(int64(mm));
            i = find(o.m == ms, 1);
            if isempty(i)
                c = exactQ(0);
            else
                c = exactQ(o.cn(i), exactQ.mulChecked(o.cd(i), k));
            end
        end

        function t = eq(x, y)
        % EXACT equality, and it is a vector comparison because the form is canonical: the sqrt(m)
        % for distinct squarefree m are linearly independent over Q, so equal values have equal
        % coefficients term by term.
            [x, y] = exactQ.pair(x, y);
            t = numel(x.m) == numel(y.m) && all(x.m == y.m) && ...
                all(x.cn == y.cn) && all(x.cd == y.cd);
        end
        function t = ne(x, y), t = ~eq(x, y); end

        function s = sign(o)
        % The sign, exactly -- and "exactly" survives the floating-point screen below, because
        % that screen only ever answers when it has a CERTIFIED margin.
        %
        % WHY THERE IS A SCREEN AT ALL. signExact squares its operands once per extension, so a
        % value carrying k of them costs coefficients to the power 2^k. Measured: sqrt(2) +
        % sqrt(3) - 3146264/1000000 -- two extensions and a seven-digit denominator, which is an
        % ordinary size for a vertex coordinate here -- overflows int64 on the second squaring
        % (1e24 against a 9.2e18 ceiling). Refusing to answer there would make the type unusable
        % on its own target inputs.
        %
        % WHY IT IS STILL EXACT. Zero is decided by the REPRESENTATION and never by arithmetic
        % (isZero: the sqrt(m) are independent over Q). So the screen is only ever asked about a
        % value already known to be nonzero, and it answers only when the computed value exceeds a
        % rigorous bound on its own rounding error; otherwise the exact recursion runs. Doubles
        % remain refused as a REPRESENTATION -- that is the ULP defect this type exists to remove
        % -- and this is not that: nothing stored, nothing rounded, no verdict without a margin.
            if isempty(o.m), s = 0; return, end
            if numel(o.m) == 1
                s = double(sign(o.cn(1))); return          % one term, and sqrt(m) > 0
            end
            [v, err] = approxWithBound(o);
            if v > err,  s =  1; return, end
            if v < -err, s = -1; return, end
            s = signExact(o);
        end

        function [v, err] = approxWithBound(o)
        % The value in double arithmetic, with a bound on the error of THAT computation. Each term
        % costs at most a handful of correctly-rounded operations (two int64-to-double
        % conversions, a division, a square root, a multiply) and the sum costs one more per term,
        % so (numel + 6) units of eps against the sum of the terms' magnitudes is generous. The
        % magnitudes, not the total: cancellation is exactly the case this has to bound honestly.
            v = 0; mag = 0;
            for i = 1:numel(o.m)
                t = double(o.cn(i))/double(o.cd(i)) * sqrt(double(o.m(i)));
                v = v + t;
                mag = mag + abs(t);
            end
            err = (numel(o.m) + 6) * eps * mag;
            if ~isfinite(v) || ~isfinite(err)
                v = 0; err = inf;                          % out of double range: decide exactly
            end
        end

        function s = signExact(o)
        % The decision procedure, with no floating point anywhere.
        %
        % Split off ONE prime: x = a + b*sqrt(p), with a and b in the field of the remaining
        % primes. If a or b vanishes the answer is the other's sign (sqrt(p) > 0). If they share a
        % sign, that is the sign. Otherwise |a| and |b|*sqrt(p) have to be compared, which squares
        % both sides into the smaller field: sign(x) = sign(a) * sign(a^2 - b^2*p), since with
        % a > 0 > b we have x > 0 exactly when a^2 > b^2*p, and with a < 0 < b the reverse.
        %
        % One prime leaves per level, so the recursion bottoms out in the rationals. Each level
        % SQUARES, which is why `sign` screens first: the growth is what int64 cannot always hold.
            if isempty(o.m), s = 0; return, end
            if numel(o.m) == 1
                s = double(sign(o.cn(1))); return
            end
            p = exactQ.pickPrime(o.m);
            [a, b] = exactQ.splitAt(o, p);
            sa = signExact(a); sb = signExact(b);
            if sa == 0, s = sb; return, end
            if sb == 0, s = sa; return, end
            if sa == sb, s = sa; return, end
            t = a*a - (b*b)*exactQ(p, 1);
            s = sa * signExact(t);
        end

        function t = lt(x, y), t = sign(minus(x, y)) < 0; end
        function t = gt(x, y), t = sign(minus(x, y)) > 0; end
        function t = le(x, y), t = sign(minus(x, y)) <= 0; end
        function t = ge(x, y), t = sign(minus(x, y)) >= 0; end

        % ---- arithmetic ---------------------------------------------------------------
        function o = plus(x, y)
            [x, y] = exactQ.pair(x, y);
            o = exactQ([x.m, y.m], [x.cn, y.cn], [x.cd, y.cd]);
        end

        function o = uminus(x)
            o = exactQ(x.m, -x.cn, x.cd);
        end

        function o = minus(x, y)
            [x, y] = exactQ.pair(x, y);
            o = plus(x, uminus(y));
        end

        function o = mtimes(x, y)
        % sqrt(m1)*sqrt(m2) = g*sqrt(m1*m2/g^2) with g = gcd(m1,m2) -- and m1*m2/g^2 is squarefree
        % when m1 and m2 are, which is the closure property this whole field rests on.
            [x, y] = exactQ.pair(x, y);
            nx = numel(x.m); ny = numel(y.m);
            mm = zeros(1, nx*ny, 'int64'); nn = mm; ddv = mm;
            t = 0;
            for i = 1:nx
                for j = 1:ny
                    g = gcd(x.m(i), y.m(j));
                    r = exactQ.mulChecked(x.m(i)/g, y.m(j)/g);
                    [pn, pd] = exactQ.mulRat(x.cn(i), x.cd(i), y.cn(j), y.cd(j));
                    [pn, pd] = exactQ.mulRat(pn, pd, g, int64(1));
                    t = t + 1;
                    mm(t) = r; nn(t) = pn; ddv(t) = pd;
                end
            end
            o = exactQ(mm(1:t), nn(1:t), ddv(1:t));
        end

        function o = inv(x)
        % 1/x by rationalising ONE PRIME AT A TIME. Multiplying by the sqrt(p)-flip of the current
        % partial product removes p from it (a^2 - b^2*p has no sqrt(p) term), so after one pass
        % per prime the product is rational -- that product is the norm, and the accumulated
        % conjugates are the numerator. For k primes this is k multiplications, not 2^k.
            if isZero(x)
                error('exactQ:divideByZero', 'exact division by zero.');
            end
            acc = exactQ(1);
            r = x;
            ps = exactQ.primesOf(x.m);
            for i = 1:numel(ps)
                p = ps(i);
                if ~any(mod(r.m, p) == 0), continue, end     % already gone: do not square for free
                c = exactQ.conjAt(r, p);
                acc = acc * c;
                r = r * c;
            end
            if ~isRational(r)
                error('exactQ:internal', ...
                    'the norm of a multiquadratic element must be rational (got %s).', char(r));
            end
            o = acc * exactQ(r.cd(1), r.cn(1));
        end

        function o = mrdivide(x, y)
            [x, y] = exactQ.pair(x, y);
            o = mtimes(x, inv(y));
        end

        function o = power(x, k)
            k = double(k);
            if k < 0, o = inv(power(x, -k)); return, end
            o = exactQ(1);
            for i = 1:k, o = mtimes(o, x); end
        end
        function o = mpower(x, k), o = power(x, k); end

        % ---- interop ------------------------------------------------------------------
        function v = double(o)
            v = 0;
            for i = 1:numel(o.m)
                v = v + double(o.cn(i))/double(o.cd(i)) * sqrt(double(o.m(i)));
            end
        end

        function s = sym(o)
        % For tests and printing ONLY. Nothing on the compute path may call this -- that is the
        % whole point of the type.
            s = sym(0);
            for i = 1:numel(o.m)
                s = s + sym(o.cn(i))/sym(o.cd(i)) * sqrt(sym(o.m(i)));
            end
        end

        function disp(o)
            fprintf('%s\n', char(o));
        end

        function str = char(o)
            if isempty(o.m), str = '0'; return, end
            str = '';
            for i = 1:numel(o.m)
                if o.cd(i) == 1, c = sprintf('%d', o.cn(i));
                else,            c = sprintf('%d/%d', o.cn(i), o.cd(i));
                end
                if o.m(i) == 1, term = c; else, term = sprintf('%s*sqrt(%d)', c, o.m(i)); end
                if i == 1, str = term; else, str = sprintf('%s + %s', str, term); end
            end
        end
    end

    methods (Static)
        function [n, dd] = fromDouble(x, maxDen)
        % A double to an EXACT rational, refusing anything it cannot represent exactly rather
        % than rounding. Integers and dyadic rationals are exact by construction; everything else
        % goes through a continued-fraction convergent and is ACCEPTED ONLY IF it reproduces the
        % double bit for bit. A silent approximation here would defeat the whole type.
            if nargin < 2, maxDen = int64(1e9); end
            if ~isscalar(x) || ~isnumeric(x) || ~isfinite(x)
                error('exactQ:fromDouble', 'need a finite numeric scalar.');
            end
            if x == fix(x) && abs(x) < 9.2e18
                n = int64(x); dd = int64(1); return
            end
            [nn, ddd] = rat(x, 1e-15);
            if abs(ddd) <= double(maxDen) && nn/ddd == x
                [n, dd] = exactQ.norm2(int64(nn), int64(ddd)); return
            end
            error('exactQ:inexact', ...
                ['%.17g is not exactly representable as a small rational; convert at the ' ...
                 'source instead of rounding here.'], x);
        end

        function o = surd(d)
        % sqrt(d) itself, with d reduced to squarefree form times a rational factor.
            o = exactQ(int64(d), int64(1), int64(1));
        end

        function [k, dsq] = squarefree(d)
        % d = k^2 * dsq with dsq squarefree. Extracting the square factor is what keeps two values
        % that live in the same field COMPARABLE -- sqrt(8) and sqrt(2) must not be different
        % radicands, or the canonical form stops being canonical and equality stops being exact.
            k = int64(1); dsq = int64(d);
            if dsq <= 1, return, end
            f = int64(2);
            while f*f <= dsq
                while mod(dsq, f*f) == 0
                    dsq = dsq / (f*f); k = k * f;
                end
                f = f + 1;
            end
        end
    end

    methods (Static, Access = private)
        function [mm, nn, ddv] = canon(mv, nv, dv)
        % The canonical form: squarefree radicands, coefficients in lowest terms, no zero
        % coefficients, strictly increasing radicands, duplicates merged. Everything the type
        % promises -- exact equality, exact zero-testing -- is a consequence of this being the
        % ONLY form a value is ever stored in.
            mv = int64(mv(:).'); nv = int64(nv(:).'); dv = int64(dv(:).');
            if numel(nv) ~= numel(mv) || numel(dv) ~= numel(mv)
                error('exactQ:args', 'radicands and coefficients must have the same length.');
            end
            if any(mv < 0)
                error('exactQ:negativeRadicand', 'radicands must be >= 0.');
            end
            keep = (nv ~= 0) & (mv ~= 0);        % a zero coefficient, or sqrt(0), is no term
            mv = mv(keep); nv = nv(keep); dv = dv(keep);
            for i = 1:numel(mv)
                [k, ms] = exactQ.squarefree(mv(i));
                mv(i) = ms;
                if k ~= 1
                    [nv(i), dv(i)] = exactQ.mulRat(nv(i), dv(i), k, int64(1));
                else
                    [nv(i), dv(i)] = exactQ.norm2(nv(i), dv(i));
                end
            end
            [mv, ord] = sort(mv); nv = nv(ord); dv = dv(ord);
            mm = zeros(1, numel(mv), 'int64'); nn = mm; ddv = mm;
            t = 0;
            i = 1;
            while i <= numel(mv)
                n = nv(i); d = dv(i);
                j = i + 1;
                while j <= numel(mv) && mv(j) == mv(i)
                    [n, d] = exactQ.addRat(n, d, nv(j), dv(j));
                    j = j + 1;
                end
                if n ~= 0
                    t = t + 1;
                    mm(t) = mv(i); nn(t) = n; ddv(t) = d;
                end
                i = j;
            end
            mm = mm(1:t); nn = nn(1:t); ddv = ddv(1:t);
        end

        function [x, y] = pair(x, y)
            if ~isa(x, 'exactQ'), x = exactQ(x); end
            if ~isa(y, 'exactQ'), y = exactQ(y); end
        end

        function p = pickPrime(mv)
        % Any prime dividing one of the radicands. The recursion in `sign` only needs SOME prime
        % to eliminate; taking the smallest keeps the intermediate values as small as the choice
        % can make them.
            ps = exactQ.primesOf(mv);
            if isempty(ps)
                error('exactQ:internal', 'no prime to split on.');
            end
            p = ps(1);
        end

        function ps = primesOf(mv)
        % The primes dividing any radicand, ascending. The radicands are squarefree and small --
        % they are products of slopes of the input polygon -- so trial division is the right tool.
            ps = zeros(1,0,'int64');
            for i = 1:numel(mv)
                v = mv(i);
                f = int64(2);
                while f*f <= v
                    if mod(v, f) == 0
                        ps(end+1) = f; %#ok<AGROW>
                        while mod(v, f) == 0, v = v / f; end
                    end
                    f = f + 1;
                end
                if v > 1, ps(end+1) = v; end %#ok<AGROW>
            end
            ps = unique(ps);
        end

        function [a, b] = splitAt(o, p)
        % o = a + b*sqrt(p), with a and b free of sqrt(p).
            isp = mod(o.m, p) == 0;
            a = exactQ(o.m(~isp), o.cn(~isp), o.cd(~isp));
            b = exactQ(idivide(o.m(isp), p, 'fix'), o.cn(isp), o.cd(isp));
        end

        function c = conjAt(o, p)
        % The sqrt(p)-flip: a + b*sqrt(p) -> a - b*sqrt(p). Every term whose radicand is divisible
        % by p carries one factor of sqrt(p), so flipping is a sign change on those terms.
            s = ones(1, numel(o.m), 'int64');
            s(mod(o.m, p) == 0) = int64(-1);
            c = exactQ(o.m, o.cn .* s, o.cd);
        end

        function [n, dd] = norm2(n, dd)
            if dd == 0, error('exactQ:zeroDenominator', 'zero denominator.'); end
            if dd < 0, n = -n; dd = -dd; end
            g = gcd(abs(n), dd);
            if g > 1, n = n / g; dd = dd / g; end
            if n == 0, dd = int64(1); end
        end

        function [n, dd] = addRat(n1, d1, n2, d2)
            g = gcd(d1, d2);
            l = exactQ.mulChecked(d1, d2 / g);
            n = exactQ.mulChecked(n1, l / d1);
            n = n + exactQ.mulChecked(n2, l / d2);
            [n, dd] = exactQ.norm2(n, l);
        end

        function [n, dd] = mulRat(n1, d1, n2, d2)
            if n1 == 0 || n2 == 0, n = int64(0); dd = int64(1); return, end
            g1 = gcd(abs(n1), d2); g2 = gcd(abs(n2), d1);   % cross-cancel BEFORE multiplying:
            n1 = n1 / g1; d2 = d2 / g1;                     % the cheapest overflow defence there
            n2 = n2 / g2; d1 = d1 / g2;                     % is not to create the big number
            n = exactQ.mulChecked(n1, n2);
            dd = exactQ.mulChecked(d1, d2);
            [n, dd] = exactQ.norm2(n, dd);
        end

        function p = mulChecked(a, b)
        % int64 multiply that RAISES on overflow instead of saturating. MATLAB's int64 saturates
        % silently, which would turn an exact type into a wrong one -- the single worst outcome
        % available here.
            p = a * b;
            if a ~= 0 && (p / a ~= b || (abs(a) > 1 && abs(b) > 1 && abs(p) == intmax('int64')))
                error('exactQ:overflow', ...
                    ['int64 overflow multiplying %d by %d. The coefficients have grown past ' ...
                     'what this type can hold exactly -- that is a signal about the pipeline, ' ...
                     'not a case to round away.'], a, b);
            end
        end
    end
end
