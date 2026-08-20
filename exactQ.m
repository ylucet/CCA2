classdef exactQ
% exactQ  Exact arithmetic in Q(sqrt(d)): values `(a + b*sqrt(d))` with a, b rational.
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
% irrational and no rational split reduces the 3-convex-edge case to something conjugable
% directly. One quadratic extension is exactly what the algorithm generates -- A.4's foot, the
% 45-degree frame's 1/sqrt(2), `bilinearFrame`'s sqrt of the eigenvalues -- and no more.
%
% THE INVARIANT, and everything below depends on it: `d` is a squarefree POSITIVE integer, or 0
% for a purely rational value. Two exactQ values combine only when their `d` agree or one of them
% is rational. That is a real restriction and it is DELIBERATE: silently promoting to a
% multi-extension tower is how an exact type turns into a slow symbolic engine by accident. A
% mismatch raises, so the caller learns which operation needs the tower rather than paying for one
% everywhere.
%
% RATIONALS are int64 numerator/denominator, always in lowest terms with a positive denominator.
% Overflow raises rather than wrapping -- a wrong answer that looks exact is the one outcome worse
% than a slow one.

    properties (SetAccess = immutable)
        an int64 = 0      % rational part: an/ad
        ad int64 = 1
        bn int64 = 0      % surd part:     bn/bd * sqrt(d)
        bd int64 = 1
        d  int64 = 0      % squarefree > 1, or 0 when the value is rational
    end

    methods
        function o = exactQ(an, ad, bn, bd, d)
        % exactQ(x)                     from a double or an integer, exactly when it is a dyadic
        %                               or small rational (see fromDouble)
        % exactQ(n, dd)                 the rational n/dd
        % exactQ(an, ad, bn, bd, d)     (an/ad) + (bn/bd)*sqrt(d)
            if nargin == 0, return, end
            if nargin == 1
                if isa(an, 'exactQ'), o = an; return, end
                [n, dd] = exactQ.fromDouble(an);
                an = n; ad = dd; bn = 0; bd = 1; d = 0;
            elseif nargin == 2
                bn = 0; bd = 1; d = 0;
            elseif nargin ~= 5
                error('exactQ:args', 'exactQ takes 1, 2 or 5 arguments.');
            end
            [an, ad] = exactQ.norm2(int64(an), int64(ad));
            [bn, bd] = exactQ.norm2(int64(bn), int64(bd));
            d = int64(d);
            if d < 0
                error('exactQ:negativeRadicand', 'd must be >= 0 (got %d).', d);
            end
            if bn == 0
                d = int64(0); bn = int64(0); bd = int64(1);   % no surd part: normalise d away
            elseif d == 0 || d == 1
                % sqrt(0) = 0 and sqrt(1) = 1: fold the surd part into the rational part.
                if d == 1
                    [an, ad] = exactQ.addRat(an, ad, bn, bd);
                end
                bn = int64(0); bd = int64(1); d = int64(0);
            end
            o.an = an; o.ad = ad; o.bn = bn; o.bd = bd; o.d = d;
        end

        % ---- predicates ---------------------------------------------------------------
        function t = isRational(o), t = (o.d == 0); end
        function t = isZero(o),     t = (o.an == 0) && (o.bn == 0); end

        function t = eq(x, y)
        % EXACT equality. Two values in the same field are equal iff both components are, because
        % {1, sqrt(d)} is a basis over Q when d is squarefree and not a perfect square.
            [x, y] = exactQ.pair(x, y);
            t = (x.an == y.an) && (x.ad == y.ad) && (x.bn == y.bn) && (x.bd == y.bd) && (x.d == y.d);
        end
        function t = ne(x, y), t = ~eq(x, y); end

        function s = sign(o)
        % EXACT sign, with no floating point anywhere. For a + b*sqrt(d): if a and b share a sign
        % the answer is that sign; otherwise compare a^2 against b^2*d, which is exact in integers.
            if o.bn == 0
                s = double(sign(o.an)); return
            end
            if o.an == 0
                s = double(sign(o.bn)); return
            end
            sa = sign(o.an); sb = sign(o.bn);
            if sa == sb, s = double(sa); return, end
            % a and b have opposite signs: |a| vs |b|*sqrt(d)  <=>  a^2*bd^2 vs b^2*d*ad^2
            l = exactQ.mulChecked(exactQ.mulChecked(o.an, o.an), exactQ.mulChecked(o.bd, o.bd));
            r = exactQ.mulChecked(exactQ.mulChecked(o.bn, o.bn), ...
                                  exactQ.mulChecked(o.d, exactQ.mulChecked(o.ad, o.ad)));
            if l == r
                s = 0;
            elseif l > r
                s = double(sa);
            else
                s = double(sb);
            end
        end

        function t = lt(x, y), t = sign(minus(x, y)) < 0; end
        function t = gt(x, y), t = sign(minus(x, y)) > 0; end
        function t = le(x, y), t = sign(minus(x, y)) <= 0; end
        function t = ge(x, y), t = sign(minus(x, y)) >= 0; end

        % ---- arithmetic ---------------------------------------------------------------
        function o = plus(x, y)
            [x, y, dd] = exactQ.align(x, y);
            [an, ad] = exactQ.addRat(x.an, x.ad, y.an, y.ad);
            [bn, bd] = exactQ.addRat(x.bn, x.bd, y.bn, y.bd);
            o = exactQ(an, ad, bn, bd, dd);
        end

        function o = uminus(x)
            o = exactQ(-x.an, x.ad, -x.bn, x.bd, x.d);
        end

        function o = minus(x, y)
            [x, y] = exactQ.pair(x, y);
            o = plus(x, uminus(y));
        end

        function o = mtimes(x, y)
        % (a1 + b1 s)(a2 + b2 s) = (a1 a2 + b1 b2 d) + (a1 b2 + a2 b1) s
            [x, y, dd] = exactQ.align(x, y);
            [p1n, p1d] = exactQ.mulRat(x.an, x.ad, y.an, y.ad);
            [p2n, p2d] = exactQ.mulRat(x.bn, x.bd, y.bn, y.bd);
            [p2n, p2d] = exactQ.mulRat(p2n, p2d, dd, int64(1));
            [an, ad]   = exactQ.addRat(p1n, p1d, p2n, p2d);
            [q1n, q1d] = exactQ.mulRat(x.an, x.ad, y.bn, y.bd);
            [q2n, q2d] = exactQ.mulRat(y.an, y.ad, x.bn, x.bd);
            [bn, bd]   = exactQ.addRat(q1n, q1d, q2n, q2d);
            o = exactQ(an, ad, bn, bd, dd);
        end

        function o = mrdivide(x, y)
        % Division by rationalising the denominator: 1/(a + b s) = (a - b s)/(a^2 - b^2 d).
            [x, y, dd] = exactQ.align(x, y);
            if isZero(y)
                error('exactQ:divideByZero', 'exact division by zero.');
            end
            conj_y = exactQ(y.an, y.ad, -y.bn, y.bd, dd);
            num = mtimes(x, conj_y);
            den = mtimes(y, conj_y);            % rational by construction
            if ~isRational(den)
                error('exactQ:internal', 'the norm of a Q(sqrt(d)) element must be rational.');
            end
            [an, ad] = exactQ.mulRat(num.an, num.ad, den.ad, den.an);
            [bn, bd] = exactQ.mulRat(num.bn, num.bd, den.ad, den.an);
            o = exactQ(an, ad, bn, bd, dd);
        end

        function o = inv(x), o = mrdivide(exactQ(1), x); end

        function o = power(x, k)
            k = double(k);
            if k < 0, o = inv(power(x, -k)); return, end
            o = exactQ(1);
            for i = 1:k, o = mtimes(o, x); end
        end
        function o = mpower(x, k), o = power(x, k); end

        % ---- interop ------------------------------------------------------------------
        function v = double(o)
            v = double(o.an)/double(o.ad);
            if o.bn ~= 0
                v = v + double(o.bn)/double(o.bd) * sqrt(double(o.d));
            end
        end

        function s = sym(o)
        % For tests and printing ONLY. Nothing on the compute path may call this -- that is the
        % whole point of the type.
            s = sym(o.an)/sym(o.ad);
            if o.bn ~= 0
                s = s + sym(o.bn)/sym(o.bd) * sqrt(sym(o.d));
            end
        end

        function disp(o)
            fprintf('%s\n', char(o));
        end

        function str = char(o)
            if o.ad == 1, str = sprintf('%d', o.an); else, str = sprintf('%d/%d', o.an, o.ad); end
            if o.bn ~= 0
                if o.bd == 1, b = sprintf('%d', o.bn); else, b = sprintf('%d/%d', o.bn, o.bd); end
                str = sprintf('%s + %s*sqrt(%d)', str, b, o.d);
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
            [k, dsq] = exactQ.squarefree(int64(d));
            o = exactQ(0, 1, k, 1, dsq);
        end

        function [k, dsq] = squarefree(d)
        % d = k^2 * dsq with dsq squarefree. Extracting the square factor is what keeps two values
        % that live in the same field COMPARABLE -- sqrt(8) and sqrt(2) must not be different `d`.
            k = int64(1); dsq = int64(d);
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
        function [x, y] = pair(x, y)
            if ~isa(x, 'exactQ'), x = exactQ(x); end
            if ~isa(y, 'exactQ'), y = exactQ(y); end
        end

        function [x, y, dd] = align(x, y)
        % Both operands into ONE field, or an error. See the header on why this refuses to build
        % a tower.
            [x, y] = exactQ.pair(x, y);
            if x.d == y.d
                dd = x.d;
            elseif x.d == 0
                dd = y.d;
            elseif y.d == 0
                dd = x.d;
            else
                error('exactQ:fieldMismatch', ...
                    ['cannot combine sqrt(%d) with sqrt(%d) -- exactQ carries ONE quadratic ' ...
                     'extension by design. Whatever needs both is the operation to look at.'], ...
                    x.d, y.d);
            end
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
