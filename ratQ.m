classdef ratQ
% ratQ  EXACT rational arithmetic on coefficient VECTORS, for the sym-free conjugate.
%
% WHAT IT IS FOR, AND WHY IT IS NOT exactQ. CONJ_FIELD_PROOF.md Theorem 1 says exactly where the
% irrationality of `f*` is and is not: the FACE FUNCTIONS and the EDGE CONICS are always rational,
% and only the VERTICES can leave Q (degree <= 4, and S4 in the worst case, so no tower of square
% roots reaches them). QuaCon acts on that split -- it stores faces and conics exactly, over Q, and
% stores vertices only as NAMES ("the intersection of conic i and conic j") with floating-point
% coordinates as a rebuildable cache. This class is the exact half of that arrangement.
%
% So the layer exactQ was built for -- a quadratic extension for the A.5 cevian foot, which is a
% VERTEX -- is not needed here at all, and neither is a degree-4 real-algebraic kernel. See
% ratQTest.m for what is pinned, and QuaCon.m for how it is used.
%
% THE TWO CANONICAL FORMS, and the distinction is load-bearing:
%
%   * A FACE FUNCTION has a SCALE. It is stored as an integer numerator vector over one positive
%     integer denominator, reduced by gcd -- `canon`. Doubling both is the same function; doubling
%     only the numerator is a different one.
%   * A CONIC HAS NO SCALE. {c'*m(x) = 0} and {1000*c'*m(x) = 0} are the SAME CURVE, so a conic is
%     stored as a PRIMITIVE integer vector whose first nonzero entry is positive -- `conic`. Two
%     conics are the same curve exactly when those vectors are bitwise equal.
%
% That second form is the whole point. DECISIONS.md 2026-08-17 records the defect this removes: two
% doubles of one exact number, one ULP apart, made a shared facet invisible to `merge`, and Step 3's
% cell count then grew without bound. A canonical integer vector cannot have two spellings, so the
% failure is not merely less likely -- it is unrepresentable.
%
% INTEGER-VALUED DOUBLES, NOT int64, AND THE REASON IS SAFETY. MATLAB's int64 arithmetic SATURATES
% on overflow: intmax('int64')+1 is intmax('int64'), silently, with no error and no wrap. A type
% whose failure mode is a plausible wrong answer is the one outcome worse than a slow one, and
% detecting saturation after the fact costs a comparison at every operation anyway. Doubles hold
% every integer up to 2^53 exactly, mix freely with the rest of CCA2 (which is all double), and
% support `gcd` directly -- and exactness is decided by a single magnitude test, which `chk` applies
% to every result. Beyond 2^53 this class RAISES (`ratQ:overflow`) rather than returning a rounded
% value. If that ever fires in normal use the answer is a wider integer type, never a tolerance.
%
% NOTHING HERE IS APPROXIMATE. `fromDouble` refuses a value it cannot represent exactly rather than
% snapping it, for the reason exactQ.fromDouble gives: converting at the boundary is the caller's
% job, and a type that quietly approximates in its constructor has no exactness to offer.

    properties (Constant)
        % The exact-integer ceiling for a double. Every result is checked against it.
        LIMIT = 2^53
    end

    methods (Static)

        function v = chk(v, what)
        % objective: assert every entry is an exact integer within the double's exact range.
        % [input]  v : numeric array; what : (optional) name used in the error message
        % [output] v : v unchanged, so this can wrap an expression
        %
        % This is the ONLY thing standing between an exact answer and a plausible wrong one, so it
        % is applied to every result rather than only at the boundary.
            if nargin < 2, what = 'value'; end
            if ~isreal(v) || ~all(isfinite(v(:)))
                error('ratQ:notFinite', '%s must be real and finite.', what);
            end
            if any(v(:) ~= round(v(:)))
                error('ratQ:notInteger', '%s must be integral; got a non-integer entry.', what);
            end
            if any(abs(v(:)) >= ratQ.LIMIT)
                error('ratQ:overflow', ...
                    ['%s reached %.17g, at or beyond 2^53, so it is no longer an exact integer. ' ...
                     'Widen the integer type; do NOT introduce a tolerance.'], what, max(abs(v(:))));
            end
        end

        function [n, d] = canon(n, d)
        % objective: canonical form of the rational vector n/d -- gcd removed, denominator positive.
        % [input]  n : 1 x k integer numerator vector; d : nonzero integer scalar denominator
        % [output] n, d : the same value, canonically spelled
        %
        % CANONICAL means two values are equal exactly when both vectors are bitwise equal, which is
        % what `merge` needs and what a symbolic form cannot provide (no canonical form there: equal
        % quantities can compare unequal, and `isAlways` can return Unknown outright).
            ratQ.chk(n, 'numerator'); ratQ.chk(d, 'denominator');
            if ~isscalar(d) || d == 0
                error('ratQ:badDenominator', 'the denominator must be a nonzero scalar.');
            end
            if d < 0, n = -n; d = -d; end
            g = d;
            for i = 1:numel(n)
                if n(i) ~= 0, g = gcd(g, abs(n(i))); end
                if g == 1, break; end
            end
            if g > 1, n = n / g; d = d / g; end
        end

        function c = conic(c)
        % objective: canonical form of a CONIC coefficient vector -- primitive, first nonzero entry
        %            positive. A conic has no scale: {c'm = 0} and {k*c'm = 0} are one curve.
        % [input]  c : 1 x 6 integer vector [a b c d e f] for a x^2 + b xy + c y^2 + d x + e y + f
        % [output] c : the same CURVE, canonically spelled
        %
        % An all-zero input is refused: in a QuaCon every edge names its own curve explicitly, and
        % the empty conic names nothing. (QuaPar spells a straight edge as an all-zero Ec row and
        % recovers the line from the edge's two endpoint COORDINATES -- exactly the dependence on
        % stored coordinates that H-form removes, so a straight edge here carries [0 0 0 d e f].)
            ratQ.chk(c, 'conic');
            nz = find(c ~= 0, 1);
            if isempty(nz)
                error('ratQ:zeroConic', ...
                    ['an all-zero conic names no curve. In a QuaCon every edge carries its own ' ...
                     'equation, and a straight edge is [0 0 0 d e f], not all-zero.']);
            end
            g = 0;
            for i = 1:numel(c)
                if c(i) ~= 0, g = gcd(g, abs(c(i))); end
            end
            c = c / g;
            if c(nz) < 0, c = -c; end
        end

        function tf = sameConic(c1, c2)
        % objective: do two conics describe the SAME CURVE. Exact; no tolerance anywhere.
            tf = isequal(ratQ.conic(c1), ratQ.conic(c2));
        end

        function [n, d] = add(n1, d1, n2, d2)
        % objective: n1/d1 + n2/d2, canonically. The lcm is taken rather than the product, because
        %            the product is what makes coefficients grow for no reason down a fold chain.
            l = lcm(d1, d2);
            ratQ.chk(l, 'common denominator');
            n = ratQ.chk(n1 * (l/d1) + n2 * (l/d2), 'sum numerator');
            [n, d] = ratQ.canon(n, l);
        end

        function [n, d] = sub(n1, d1, n2, d2)
        % objective: n1/d1 - n2/d2, canonically.
            [n, d] = ratQ.add(n1, d1, -n2, d2);
        end

        function [n, d] = scale(n1, d1, kn, kd)
        % objective: (kn/kd) * (n1/d1), canonically. kd defaults to 1.
            if nargin < 4, kd = 1; end
            n = ratQ.chk(n1 * kn, 'scaled numerator');
            [n, d] = ratQ.canon(n, ratQ.chk(d1 * kd, 'scaled denominator'));
        end

        function tf = eqRat(n1, d1, n2, d2)
        % objective: exact equality of two rational vectors, by canonical form.
            [a, p] = ratQ.canon(n1, d1);
            [b, q] = ratQ.canon(n2, d2);
            tf = isequal(a, b) && isequal(p, q);
        end

        function c = diffConic(n1, d1, n2, d2)
        % objective: the conic {g1 = g2} where g1 = n1/d1 and g2 = n2/d2 are QUADRATIC face
        %            functions in CCA2's 10-wide cubic basis.
        % [input]  n1,d1,n2,d2 : the two face functions, exactly
        % [output] c : 1 x 6 canonical conic [a b c d e f]
        %
        % THIS IS THE ROUTINE THE WHOLE DESIGN TURNS ON. Every edge of f* is such a conic (Theorem 1
        % of CONJ_FIELD_PROOF.md), and its discriminant b^2-4ac = -det(H1-H2) is what decides whether
        % the edge is a parabola. Theorem 3 makes it vanish for two ADJACENT pieces of f; Step 3 also
        % compares non-adjacent ones, where it does not, and that is why the Conic level exists.
        %
        % Reading the weighted basis (RatCon.m's `f`): with the cubic entries zero,
        %     g = c5*x^2/2 + c6*xy + c7*y^2/2 + c8*x + c9*y + c10
        % so clearing the halves means multiplying the whole difference by 2.
            if any(n1(1:4) ~= 0) || any(n2(1:4) ~= 0)
                error('ratQ:cubicFace', ...
                    'diffConic is for QUADRATIC faces; a cubic numerator has no conic level set.');
            end
            [dn, ~] = ratQ.sub(n1, d1, n2, d2);      % the difference function, exactly
            % [a b c d e f] = [c5/2, c6, c7/2, c8, c9, c10] / dd; times 2*dd to clear, and the
            % common denominator dd then drops out because a conic has no scale.
            c = [dn(5), 2*dn(6), dn(7), 2*dn(8), 2*dn(9), 2*dn(10)];
            c = ratQ.conic(ratQ.chk(c, 'difference conic'));
        end

        function disc = discriminant(c)
        % objective: b^2 - 4ac of a conic -- 0 for a parabola or a line, >0 hyperbolic,
        %            <0 elliptic. Exact integer arithmetic, so the classification is DECIDED,
        %            not estimated.
            ratQ.chk(c, 'conic');
            disc = ratQ.chk(c(2)^2 - 4*c(1)*c(3), 'discriminant');
        end

        function k = conicKind(c)
        % objective: name the conic's type from its EXACT discriminant.
        % [output] k : 'line' | 'parabola' | 'hyperbolic' | 'elliptic'
        %
        % 'line' is the degenerate a=b=c=0 case, which QuaPar spells as an all-zero row and QuaCon
        % spells explicitly. The two nondegenerate signs are named for the discriminant rather than
        % for the point set, because a real ellipse, an imaginary one and a single point all share
        % disc < 0 and the arrangement code only ever needs the sign.
            if all(c(1:3) == 0), k = 'line'; return, end
            d = ratQ.discriminant(c);
            if d == 0,     k = 'parabola';
            elseif d > 0,  k = 'hyperbolic';
            else,          k = 'elliptic';
            end
        end

        function s = signQ(n, d)
        % objective: the exact sign of the rational n/d, elementwise on n.
        % [output] s : -1, 0 or +1 per entry
        %
        % The leaf of every predicate. Trivial by itself, and it is here rather than at the call
        % sites so that "the sign of a rational" has ONE spelling: sign(n/d) in a caller would
        % perform a floating-point division first, which is the one thing this layer exists to
        % avoid, and n*d > 0 would overflow for no reason.
            ratQ.chk(n, 'numerator'); ratQ.chk(d, 'denominator');
            if ~isscalar(d) || d == 0
                error('ratQ:badDenominator', 'the denominator must be a nonzero scalar.');
            end
            s = sign(n) * sign(d);
        end

        function [n, d] = div(n1, d1, n2, d2)
        % objective: (n1/d1) / (n2/d2), canonically. The divisor is a SCALAR rational.
        %
        % The counterpart of scale, which is already the multiplication (scale IS kn/kd times a
        % vector, so no separate mul is added -- a second spelling of one operation is how two
        % canonical forms drift apart).
            if ~isscalar(n2)
                error('ratQ:badDivisor', ...
                    'the divisor must be a scalar rational; dividing by a vector is not defined.');
            end
            if n2 == 0
                error('ratQ:divideByZero', ...
                    'division by the exact zero: a defect in the caller, not a case to regularise.');
            end
            [n, d] = ratQ.scale(n1, d1, d2, n2);
        end

        function [N, d] = combineDen(Ns, ds)
        % objective: put several rational vectors over ONE common denominator, values unchanged.
        % [input]  Ns : k x m integer numerators, one row per value; ds : k x 1 nonzero denominators
        % [output] N  : k x m over the single denominator d = lcm(|ds|)
        %
        % Deliberately NOT canonicalised row by row: the caller wants the rows COMPARABLE (a common
        % denominator is what turns a matrix of rationals into integer arithmetic), and reducing
        % each row by its own gcd would undo exactly that. Use canon when a single value is wanted
        % in canonical form; use this when several must share a scale.
        %
        % lcm, not the product -- the same reason add gives: the product is what makes
        % coefficients grow down a fold chain for no mathematical reason.
            ratQ.chk(Ns, 'numerators'); ratQ.chk(ds, 'denominators');
            ds = ds(:);
            if size(Ns,1) ~= numel(ds)
                error('ratQ:sizeMismatch', ...
                    'combineDen needs one denominator per row: %d rows, %d denominators.', ...
                    size(Ns,1), numel(ds));
            end
            if any(ds == 0)
                error('ratQ:badDenominator', 'every denominator must be nonzero.');
            end
            d = 1;
            for i = 1:numel(ds)
                d = lcm(d, abs(ds(i)));
                ratQ.chk(d, 'common denominator');
            end
            N = ratQ.chk(Ns .* (d ./ ds), 'rescaled numerator');
        end

        function [x, xd] = solve2(A, b)
        % objective: the EXACT solution of A*x = b, for a 2x2 or 3x3 integer system.
        % [input]  A : n x n integer matrix, n = 2 or 3; b : n x 1 integer vector
        % [output] x : n x 1 integer numerators; xd : positive integer denominator (x/xd solves it)
        %
        % Where it is used: a polyhedral vertex is two lines meeting (n = 2) and the KKT systems of
        % the per-piece conjugate are n = 3. Both are small and both must be DECIDED, not estimated.
        %
        % A RATIONAL system is cleared by the caller: multiplying one row of [A b] by a nonzero
        % integer changes no solution, so scaling each row by its own denominator reduces to this.
        %
        % CRAMER, AND NOT MATLAB det. det factorises in floating point: on an integer matrix it
        % returns 6.0000000000000009 for 6 and, worse, something like 4.4e-16 for a singular
        % system -- so a singularity test built on it is a tolerance test wearing a disguise. The
        % cofactor expansion below is exact integer arithmetic and the singular case is DECIDED.
        % Cramer is the wrong algorithm at size 20 and the right one at size 3, where it is a
        % handful of multiplications and needs no pivoting, hence no ordering choice to get wrong.
            ratQ.chk(A, 'system matrix'); ratQ.chk(b, 'right-hand side');
            n = size(A, 1);
            if size(A,2) ~= n || ~ismember(n, [2 3]) || ~isequal(size(b), [n 1])
                error('ratQ:badSystem', ...
                    'solve2 takes a 2x2 or 3x3 matrix with a matching n x 1 column.');
            end
            D = ratQ.detExact(A);
            if D == 0
                error('ratQ:singular', ...
                    ['the system is exactly singular (determinant 0), so it has no unique ' ...
                     'solution. The caller must handle the degenerate configuration -- parallel ' ...
                     'edges, a collapsed cell -- rather than perturbing it.']);
            end
            x = zeros(n, 1);
            for j = 1:n
                Aj = A;  Aj(:,j) = b;
                x(j) = ratQ.detExact(Aj);
            end
            [x, xd] = ratQ.canon(x, D);
        end

        function D = detExact(A)
        % objective: the determinant of a 2x2 or 3x3 integer matrix, by cofactor expansion.
        % Exact integer arithmetic; see solve2 for why MATLAB det cannot be used here.
            ratQ.chk(A, 'matrix');
            switch size(A,1)
                case 1
                    D = A(1,1);
                case 2
                    D = ratQ.chk(A(1,1)*A(2,2) - A(1,2)*A(2,1), 'determinant');
                case 3
                    D = ratQ.chk( ...
                        A(1,1)*(A(2,2)*A(3,3) - A(2,3)*A(3,2)) ...
                      - A(1,2)*(A(2,1)*A(3,3) - A(2,3)*A(3,1)) ...
                      + A(1,3)*(A(2,1)*A(3,2) - A(2,2)*A(3,1)), 'determinant');
                otherwise
                    error('ratQ:badSystem', 'detExact handles n <= 3; got %d.', size(A,1));
            end
        end

        function [H, d] = hessQ(n, dn)
        % objective: the exact Hessian of a QUADRATIC face n/dn, as H/d with H integer.
        % [output] H : 2 x 2 integer, symmetric; d : positive integer
        %
        % This is one line of arithmetic because the weighted basis was chosen to make it one:
        % RatCon.m stores f with the weights [1/6 1/2 1/2 1/6 1/2 1 1/2 1 1 1] precisely "to
        % easily manipulate Hessians", so H is [c5 c6; c6 c7] with no factor to remember.
        %
        % A CUBIC IS REFUSED rather than evaluated at a point, and that is the contract: a cubic
        % Hessian is not constant, so returning one would be a wrong answer rather than a missing
        % feature. QuaPar.evalHessian is the routine that takes a point.
            if any(n(1:4) ~= 0)
                error('ratQ:cubicFace', ...
                    ['hessQ is for QUADRATIC faces: a cubic Hessian varies with the point, so ' ...
                     'there is no constant matrix to return. Use QuaPar.evalHessian(c,x).']);
            end
            [H, d] = ratQ.canon([n(5) n(6); n(6) n(7)], dn);
        end

        function [G, d] = gradQ(n, dn)
        % objective: the exact gradient of a QUADRATIC face n/dn, as the AFFINE map G/d.
        % [output] G : 2 x 3 integer, so that grad(x,y) = G*[x; y; 1] / d
        %          d : positive integer
        %
        % Returning the MAP rather than a value at a point is what the callers want: the normal
        % cone of a face, the KKT stationarity row and the conjugate change of variables are all
        % statements about grad AS A FUNCTION, and evaluating it early forces a point to be chosen
        % before it is known.
        %
        % Refuses a cubic for the same reason hessQ does -- the gradient is then quadratic, not
        % affine, and would not fit this shape.
            if any(n(1:4) ~= 0)
                error('ratQ:cubicFace', ...
                    ['gradQ is for QUADRATIC faces: a cubic gradient is quadratic, not affine, ' ...
                     'so it does not fit the 2x3 map. Use QuaPar.evalHessian(c,x).']);
            end
            [G, d] = ratQ.canon([n(5) n(6) n(8); n(6) n(7) n(9)], dn);
        end

        function [n, d] = fromDouble(x, maxDen)
        % objective: the rational vector a caller MEANT, recovered from the doubles they passed.
        % [input]  x : 1 x k double; maxDen : largest denominator to accept (default 10^6)
        % [output] n, d : integers, canonical, with n/d ROUNDING BACK to x bit for bit
        %
        % THIS IS A RECONSTRUCTION AT THE INPUT BOUNDARY, AND IT HAS TO BE. A double is a BINARY
        % rational, so 1/3, 1/11 and 3/10 are not representable at all: a caller who writes 1/3
        % hands over 6004799503160661/2^54, and demanding a bit-exact rational reading of that
        % returns a 2^54 denominator, which is useless and then overflows. What is wanted is the
        % rational the caller meant, and it is UNIQUELY determined: among rationals with
        % denominator <= maxDen there is at most one within half an ulp of x whenever maxDen is
        % below the separation bound, which 10^6 is for any double of moderate size.
        %
        % SO THIS IS NOT THE REFUTED "RATIONAL SNAPPING", and the distinction is exactly the one
        % DECISIONS.md draws. That attempt took COMPUTED geometry -- vertices out of a
        % double-precision convEnvCPLQ -- and snapped it to the simplest rational within 1e-10;
        % it failed because the conjugate is a rational function of those coordinates, so bounding
        % the vertex denominators did not bound the downstream ones and a few squarings carried
        % 1e5 to 1e25. Here nothing computed is ever converted: this runs ONCE, on the caller's own
        % input data, and everything after it is exact integer arithmetic. Applying it to a
        % computed coordinate would be that refuted design, and QuaCon never does -- vertices are
        % NAMED, not converted (see QuaCon.m).
        %
        % IT STILL REFUSES AN IRRATIONAL, and the reason it can is a counting argument rather than
        % a tolerance: the best rational approximation to an irrational with denominator q is off
        % by about 1/q^2, so at q <= 10^6 the error is ~10^-12, four orders of magnitude larger
        % than the ~10^-16 window this accepts. sqrt(2) is rejected; 50/11 is accepted.
            if nargin < 2, maxDen = 1e6; end
            d = 1;
            for i = 1:numel(x)
                d = lcm(d, ratQ.denomOf(x(i), maxDen));
                if d > ratQ.LIMIT
                    error('ratQ:overflow', ...
                        'the common denominator exceeded 2^53 while converting %.17g.', x(i));
                end
            end
            n = round(x * d);
            ratQ.chk(n, 'converted numerator');
            bad = find(n/d ~= x, 1);
            if ~isempty(bad)
                error('ratQ:notExact', ...
                    ['%.17g does not round-trip through the common denominator %d. The vector ' ...
                     'mixes denominators whose lcm is too coarse for one of its entries; convert ' ...
                     'the entries separately, or raise maxDen.'], x(bad), d);
            end
            [n, d] = ratQ.canon(n, d);
        end

        function d = denomOf(x, maxDen)
        % objective: the smallest q <= maxDen such that some p/q rounds to the double x; raises if
        %            there is none, which is what makes an irrational input an error rather than a
        %            silent approximation. See fromDouble for why this is a reconstruction.
        %
        % Continued fractions, not a scan: a scan to 10^6 is 10^6 trials per entry. The expansion
        % is stopped by the ROUND-TRIP test (does p/q evaluate back to exactly this double), which
        % is the sharpest stopping rule available and needs no tolerance parameter of its own.
            if x == round(x), d = 1; return, end
            ax = abs(x);
            p0 = 1; q0 = 0; p1 = floor(ax); q1 = 1; r = ax - p1;
            for it = 1:64
                if r == 0, break, end
                a = floor(1/r); r = 1/r - a;
                p2 = a*p1 + p0; q2 = a*q1 + q0;
                if q2 > maxDen, break, end
                p0 = p1; q0 = q1; p1 = p2; q1 = q2;
                if p1/q1 == ax, d = q1; return, end
            end
            error('ratQ:notExact', ...
                ['%.17g is not a rational with denominator <= %d to within a rounding. It is ' ...
                 'most likely an irrational produced upstream (a VERTEX, not a coefficient) -- ' ...
                 'and by CONJ_FIELD_PROOF.md Theorem 1 no COEFFICIENT of f* can be one, so this ' ...
                 'is a defect in the caller and not a case to round away.'], x, maxDen);
        end

        function v = evalFace(n, d, X)
        % objective: evaluate the face function n/d at points X (k x 2), in double precision.
        % Evaluation is deliberately NOT exact: CONJ_FIELD_PROOF.md 8.0 lists the operations that
        % must leave Q, and evaluating f* at a point is not one of them. Exactness is for the
        % PREDICATES that decide the mesh; the mesh once built is evaluated in doubles.
            v = QuaPar.evalPoly(n, X) / d;
        end
    end
end
