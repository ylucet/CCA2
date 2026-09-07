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

        function [n2, d2] = substAffine(n, d, Mn, tn, md)
        % objective: the EXACT affine change of variables h(u) = g(M*u + t), where g = n/d is a
        %            face in RatCon's 10-wide weighted cubic basis and M = Mn/md, t = tn/md.
        % [input]  n, d : the face, exactly; Mn : 2x2 integer; tn : 2x1 integer;
        %          md   : (optional) the common denominator of the map, default 1
        % [output] n2, d2 : h, canonically, in the same basis
        %
        % This is substituteFrame.m's exact twin. That routine does the same substitution with
        % `subs` on a sym, which is where 21 of plq_1p's engine calls live; this one is integer
        % arithmetic and returns a canonical form, so two faces that ARE equal after a frame change
        % compare equal -- which `subs` cannot promise (no canonical form; equal quantities can
        % compare Unknown).
        %
        % WHY THE CUBIC TERMS ARE CARRIED. `f` is stored in the cubic basis and isConvex accepts a
        % cubic (Rat.m's closing note), so a routine that silently truncated to the quadratic part
        % would pass every quadratic test and be wrong on the one input that reaches it.
        %
        % A SINGULAR M IS REFUSED. A change of VARIABLES is invertible by definition; a singular M
        % collapses the plane onto a line, and g restricted to that line is a different function on
        % a lower-dimensional domain, not a reparametrisation of this one. Returning the collapsed
        % coefficients would be a plausible wrong answer, which is the outcome this class exists to
        % prevent.
        %
        % HOW IT STAYS INTEGRAL. The weighted basis carries sixths and halves, so the plain
        % monomial coefficients are taken over 6*d:
        %       6*g = c1 x^3 + 3 c2 x^2 y + 3 c3 x y^2 + c4 y^3
        %           + 3 c5 x^2 + 6 c6 xy + 3 c7 y^2 + 6 c8 x + 6 c9 y + 6 c10      (all over d)
        % Substituting a map of denominator md multiplies a degree-k term by md^-k, so every term
        % is lifted to the common denominator md^3 and the result sits over 6*d*md^3. Converting
        % back to the weighted basis is the inverse scaling, and `canon` removes whatever the two
        % conversions left behind.
            if nargin < 5 || isempty(md), md = 1; end
            ratQ.chk(n, 'face'); ratQ.chk(d, 'face denominator');
            ratQ.chk(Mn, 'map matrix'); ratQ.chk(tn, 'map translation');
            ratQ.chk(md, 'map denominator');
            if ~isequal(size(Mn), [2 2]) || ~isequal(size(tn), [2 1])
                error('ratQ:badMap', 'substAffine needs a 2x2 M and a 2x1 t.');
            end
            if ~isscalar(md) || md == 0
                error('ratQ:badDenominator', 'the map denominator must be a nonzero scalar.');
            end
            if ratQ.detExact(Mn) == 0
                error('ratQ:singular', ...
                    ['the change of variables is singular: M collapses the plane onto a line, ' ...
                     'so g(M*u+t) is a different function on a lower-dimensional domain, not a ' ...
                     'reparametrisation. Handle the degenerate map at the call site.']);
            end
            if md < 0, Mn = -Mn; tn = -tn; md = -md; end

            % ---- the source face, as PLAIN monomial coefficients over 6*d ---------------------
            % P is indexed P(i+1, j+1) = coefficient of x^i y^j, with i + j <= 3.
            P = zeros(4,4);
            P(4,1) = n(1);      P(3,2) = 3*n(2);    P(2,3) = 3*n(3);    P(1,4) = n(4);
            P(3,1) = 3*n(5);    P(2,2) = 6*n(6);    P(1,3) = 3*n(7);
            P(2,1) = 6*n(8);    P(1,2) = 6*n(9);
            P(1,1) = 6*n(10);
            P = ratQ.chk(P, 'plain face coefficient');

            % ---- the two substituted linear forms, as polynomials in (u,v) --------------------
            % Numerators only; the md's are accounted for by the degree lift below.
            X = zeros(4,4);  X(2,1) = Mn(1,1);  X(1,2) = Mn(1,2);  X(1,1) = tn(1);
            Y = zeros(4,4);  Y(2,1) = Mn(2,1);  Y(1,2) = Mn(2,2);  Y(1,1) = tn(2);

            % ---- accumulate P(i,j) * X^i * Y^j, each lifted to the denominator md^3 -----------
            R = zeros(4,4);
            Xp = ratQ.polyPow(X, 3);      % Xp{k+1} is X^k
            Yp = ratQ.polyPow(Y, 3);
            for i = 0:3
                for j = 0:(3-i)
                    cij = P(i+1, j+1);
                    if cij == 0, continue, end
                    term = ratQ.polyMul(Xp{i+1}, Yp{j+1});
                    lift = md^(3 - i - j);
                    R = ratQ.chk(R + cij * lift * term, 'substituted coefficient');
                end
            end

            % ---- back to the weighted basis, over 6*d*md^3 ------------------------------------
            den = ratQ.chk(6 * d * md^3, 'substituted denominator');
            % Inverting the basis weights: the entry multiplying x^i y^j is R(i+1,j+1), and the
            % weighted slot for that monomial carries 1/6 (the two cubic corners), 1/2 (the mixed
            % cubics and the pure quadratics) or 1 (xy, linear, constant), so each slot is scaled
            % by the reciprocal of its own weight.
            n2 = ratQ.chk([6*R(4,1), 2*R(3,2), 2*R(2,3), 6*R(1,4), ...
                           2*R(3,1),   R(2,2), 2*R(1,3), ...
                             R(2,1),   R(1,2), ...
                             R(1,1)] , 'substituted face');
            [n2, d2] = ratQ.canon(n2, den);
        end

        function C = polyMul(A, B)
        % objective: multiply two bivariate polynomials held as 4x4 coefficient arrays, where
        %            A(i+1,j+1) is the coefficient of x^i y^j and every term has i + j <= 3.
        % Exact integer arithmetic. A product leaving total degree 3 is a defect in the caller
        % (substAffine only ever multiplies factors whose degrees it has already summed), so it
        % raises rather than truncating.
            C = zeros(4,4);
            for i1 = 0:3
                for j1 = 0:(3-i1)
                    if A(i1+1,j1+1) == 0, continue, end
                    for i2 = 0:3
                        for j2 = 0:(3-i2)
                            if B(i2+1,j2+1) == 0, continue, end
                            i = i1 + i2;  j = j1 + j2;
                            if i + j > 3
                                error('ratQ:degreeOverflow', ...
                                    ['the product reached total degree %d, above the cubic ' ...
                                     'basis. Truncating would silently change the function.'], i+j);
                            end
                            C(i+1,j+1) = C(i+1,j+1) + A(i1+1,j1+1)*B(i2+1,j2+1);
                        end
                    end
                end
            end
            C = ratQ.chk(C, 'polynomial product');
        end

        function R = resultantCore(P, Q, k, nv, wantBig)
        % objective: the SYLVESTER RESULTANT of two nv-variate integer polynomials with respect to
        %            variable k, exactly. The general form; resultantY is the nv = 2 case.
        %
        % [input]  P, Q : arrays whose (i1+1,...,inv+1) entry is the coefficient of v1^i1...vnv^inv.
        %                 MATLAB drops trailing singleton dimensions, so nv is passed rather than
        %                 inferred -- a polynomial of degree 0 in its last variable is indexed the
        %                 same as one with no last variable at all, and guessing which would be a
        %                 silent wrong answer.
        %          k    : which variable to eliminate; nv : how many there are.
        % [output] R    : an nv-variate array whose k-th dimension is a singleton -- the resultant,
        %                 a polynomial in the REMAINING variables.
        %
        % WHY IT IS NEEDED IN THIS GENERALITY. The convex envelope's bridge face is cut out by three
        % equations in the dual variable s = (s1,s2) whose coefficients involve x, y and the value z.
        % Eliminating s takes TWO resultants, and the first one's output already has coefficients in
        % three remaining variables -- so a bivariate routine cannot chain. See DECISIONS.md
        % 2026-09-06 and .claude/bridgeEliminationSpike.m.
        %
        % The mechanics are resultantY's: the Sylvester matrix's entries are polynomials in the
        % remaining variables, its determinant is one too, and cofactor expansion uses only multiply
        % and add -- no division, so exactness is immediate.
            ratQ.chk(P, 'first polynomial');  ratQ.chk(Q, 'second polynomial');
            if ~ismember(k, 1:nv)
                error('ratQ:badVariable', 'cannot eliminate variable %d of %d.', k, nv);
            end
            ord = [setdiff(1:nv, k), k];                 % put the eliminated variable last
            Pp  = permute(P, [ord, nv+1]);
            Qp  = permute(Q, [ord, nv+1]);
            m = size(Pp, nv) - 1;
            n = size(Qp, nv) - 1;
            if m < 1 || n < 1
                error('ratQ:resultantDegree', ...
                    ['the resultant eliminates a variable, so both arguments must have degree ' ...
                     '>= 1 in it; got %d and %d. A polynomial free of it has no root to share.'], ...
                     m, n);
            end
            N = m + n;
            if N > 8
                error('ratQ:resultantTooLarge', ...
                    ['the Sylvester matrix is %dx%d, and cofactor expansion is factorial in that. ' ...
                     'Above 8 this wants Bareiss elimination instead.'], N, N);
            end

            % Slice out the coefficient of each power of the eliminated variable. Each slice is a
            % polynomial in the remaining nv-1 variables.
            Pc = cell(1, m+1);  Qc = cell(1, n+1);
            idx = repmat({':'}, 1, nv);
            for j = 0:m, ii = idx; ii{nv} = j+1; Pc{j+1} = ratQ.squeezeLast(Pp(ii{:}), nv); end
            for j = 0:n, ii = idx; ii{nv} = j+1; Qc{j+1} = ratQ.squeezeLast(Qp(ii{:}), nv); end

            S = cell(N, N);
            [S{:}] = deal(0);
            for r = 1:n
                for j = 0:m, S{r, r+j} = Pc{m-j+1}; end
            end
            for r = 1:m
                for j = 0:n, S{n+r, r+j} = Qc{n-j+1}; end
            end
            % FAST PATH FIRST, BIG INTEGERS ONLY IF IT RAISES. The determinant is where every
            % measured overflow on this line happens, and the shape is always the same: the
            % intermediates are large while the reduced ANSWER is small. So on ratQ:overflow the
            % whole determinant is redone in BigInteger, the content is taken THERE, and only then
            % does it come back to doubles -- which raises again only if the answer itself does not
            % fit, which is a different and honest failure.
            big = [];
            if wantBig
                big = ratQ.detInBig(S, N);
            else
                try
                    R = ratQ.chk(ratQ.trimN(ratQ.polyDetN(S)), 'resultant');
                catch err
                    if ~strcmp(err.identifier, 'ratQ:overflow'), rethrow(err); end
                    big = ratQ.detInBig(S, N);
                    R = ratQ.trimN(ratQ.fromBig(big));
                end
            end

            % Put the remaining variables back where they were, with a singleton in slot k.
            back = zeros(1, nv);
            rest = setdiff(1:nv, k);
            for t = 1:numel(rest), back(rest(t)) = t; end
            back(k) = nv;
            if wantBig
                R = ratQ.bigPermute(big, back, nv);
            else
                R = permute(cat(nv, R), back);
            end
        end

        function R = resultantVar(P, Q, k, nv)
        % objective: the Sylvester resultant, returned as DOUBLES. See resultantCore for the
        % mathematics and resultantVarBig for the form that does not insist the answer fit.
            R = ratQ.resultantCore(P, Q, k, nv, false);
        end

        function B = detInBig(S, N)
        % objective: the Sylvester determinant in BigInteger, content removed.
            Sb = cell(N, N);
            for r = 1:N
                for j = 1:N, Sb{r,j} = ratQ.toBig(S{r,j}); end
            end
            % COFACTOR EXPANSION, NOT BAREISS, and that is measured rather than assumed. Bareiss is
            % O(n^3) against cofactor's n!, so it looks like the obvious choice at the 8x8 this
            % needs -- but it was SLOWER on three of four cell pairs (the anisotropic one 51 s ->
            % 82 s) for two compounding reasons. A Sylvester matrix is SPARSE and banded, so
            % cofactor expansion skips most of its n! branches, while Bareiss fills in. And every
            % Bareiss step is an exact polynomial DIVISION, which profiling put at 115 s of 121 --
            % cofactor expansion divides not at all.
            B = ratQ.bigContentFree(ratQ.polyDetBig(Sb));
        end

        function B = bigPermute(A, ord, nv)
        % objective: permute a big polynomial's variables, the counterpart of MATLAB's permute.
            sa = [A.sz, ones(1, nv - numel(A.sz))];
            sb = sa(ord);
            B.sz = sb;
            B.c  = cell(prod(sb), 1);
            for i = 1:numel(A.c)
                if isempty(A.c{i}), continue, end
                e = ratQ.linToSub(i, sa);
                B.c{ratQ.subToLin(e(ord), sb)} = A.c{i};
            end
        end

        function Q = polyDivExact(P, D, nv)
        % objective: divide the integer polynomial array P by D, EXACTLY. Raises unless the
        %            division leaves no remainder and an integral quotient.
        %
        % [input]  P, D : coefficient arrays, index (i1+1,...,inv+1) = coeff of v1^i1...vnv^inv
        %          nv   : the number of variables (see resultantVar for why it is passed)
        % [output] Q    : the quotient
        %
        % WHY EXACT DIVISION AND NOT FACTORISATION. A resultant carries extraneous factors, and
        % removing them LOOKS like a factorisation problem -- but it is not, when the factors are
        % already known. The elimination multiplies by A^2 to clear a denominator, and a Sylvester
        % determinant multiplies by powers of the leading coefficients; dividing by what one
        % multiplied by is a division, not a factorisation.
        %
        % Lexicographic leading-term division, first variable most significant. It terminates
        % because the leading exponent of the remainder strictly decreases, and it is exact because
        % a non-integral quotient coefficient or a stuck remainder RAISES rather than rounding.
            ratQ.chk(P, 'dividend');  ratQ.chk(D, 'divisor');
            [eD, cD] = ratQ.leadTerm(D, nv);
            if isempty(eD)
                error('ratQ:divideByZero', 'the divisor is the zero polynomial.');
            end
            R = P;  Q = 0;  guard = 0;
            while true
                [eR, cR] = ratQ.leadTerm(R, nv);
                if isempty(eR), break, end
                guard = guard + 1;
                if guard > 100000
                    error('ratQ:divideNoConverge', 'exact division did not terminate.');
                end
                eQ = eR - eD;
                if any(eQ < 0) || mod(cR, cD) ~= 0
                    error('ratQ:notDivisible', ...
                        ['the division is not exact: the remainder''s leading term is not a ' ...
                         'multiple of the divisor''s. Nothing is truncated or rounded here.']);
                end
                T = ratQ.monoN(eQ, cR/cD, nv);
                Q = ratQ.polyAddN(Q, T);
                R = ratQ.chk(ratQ.polyAddN(R, -ratQ.polyMulN(T, D)), 'division remainder');
            end
            Q = ratQ.chk(Q, 'quotient');
        end

        function [e, c] = leadTerm(A, nv)
        % objective: the LEXICOGRAPHICALLY largest exponent vector carrying a nonzero coefficient,
        %            and that coefficient. Empty when A is the zero polynomial.
        % The first variable is most significant, which is what makes polyDivExact terminate.
            e = [];  c = 0;
            idx = find(A ~= 0);
            if isempty(idx), return, end
            sz = arrayfun(@(d) size(A,d), 1:nv);
            subs = cell(1, nv);
            [subs{:}] = ind2sub([sz 1], idx);
            E = zeros(numel(idx), nv);
            for d = 1:nv, E(:,d) = subs{d}(:) - 1; end
            [~, best] = sortrows(E, -(1:nv));
            e = E(best(1), :);
            c = A(idx(best(1)));
        end

        function A = monoN(e, c, nv) %#ok<INUSD>
        % objective: the single monomial c * v1^e1 ... vnv^env, as a coefficient array.
            sz = e + 1;
            if numel(sz) < 2, sz = [sz 1]; end
            A = zeros(sz);
            ii = num2cell(e + 1);
            A(ii{:}) = c;
        end
        function G = polyGcdBiv(A, B)
        % objective: the GCD of two bivariate integer polynomials, primitive and up to sign.
        % [input]  A, B : arrays with (i+1,j+1) the coefficient of x^i y^j
        % [output] G    : their gcd, in the same layout
        %
        % WHAT IT IS FOR. The convex envelope's bridge face comes out of resultantVar multiplied by
        % an extraneous factor that is FREE OF z -- measured on the worked instance, p has degrees
        % (10,6,2) while the face's own polynomial has (2,2,2), and the quotient (8,4,0) carries no
        % z at all. Removing it is exactly dividing by the CONTENT of p with respect to z, which is
        % the gcd of its z-coefficients as polynomials in (x,y). polyGcd answers that question for
        % UNIVARIATE polynomials only.
        %
        % HOW. Euclid in x, with the coefficients themselves polynomials in y, taking the primitive
        % part with respect to the y-content at every step -- the same discipline polyPseudoRem
        % documents one dimension down, and for the same reason.
        %
        % SCOPE, MEASURED AND NOT GUESSED. This is for inputs whose COEFFICIENTS ARE SMALL. On the
        % convex envelope's own bridge polynomial -- degrees (10,6) and (9,5) with coefficients
        % around 1e12 -- it RAISES ratQ:overflow at the second elimination (1.3e17), and that is
        % not a defect to patch: a polynomial remainder sequence multiplies by a leading
        % COEFFICIENT at every step, and here those coefficients are themselves y-polynomials with
        % 12-digit integers. Three separate fixes were tried and each moved the overflow rather
        % than removing it (.claude/gcdSwellProbe.m traces the swell step by step).
        %
        % The standard answer to intermediate expression swell is a MODULAR gcd -- work modulo
        % several small primes, where every coefficient is bounded by the prime, and reconstruct --
        % not a better PRS. That is the route if this is ever needed at scale. DECISIONS.md's 2026
        % CRT refutation does NOT apply to it: that one failed because the reconstructed VALUE
        % could not live in a double, and a gcd's coefficients are small by construction.
        %
        % What the envelope actually uses is a cheaper, VERIFIED shortcut -- see
        % .claude/bridgeReduceProbe.m -- so this routine is a general fallback, not a hot path.
            % THE INTEGER CONTENT COMES OFF FIRST, and leaving it on is what made this overflow.
            % bivPrimitive removes the y-POLYNOMIAL content, and it cannot see an integer factor:
            % ratQ.polyGcd returns a PRIMITIVE gcd, so gcd(1332e9*y, 2664e9) comes back as 1 and
            % the 1.3e12 scale survives into the pseudo-remainder, where it squares to 1.8e24.
            % A gcd is defined up to content anyway -- the result is normalised below -- so
            % dropping it here changes nothing but the size of every intermediate.
            A = ratQ.contentFree(A);  B = ratQ.contentFree(B);
            a = ratQ.bivToCell(A);  b = ratQ.bivToCell(B);
            if isempty(a), G = B; return, end
            if isempty(b), G = A; return, end
            [a, ca] = ratQ.bivPrimitive(a);
            [b, cb] = ratQ.bivPrimitive(b);
            cg = ratQ.yGcd(ca, cb);
            guard = 0;
            while true
                guard = guard + 1;
                if guard > 200
                    error('ratQ:gcdNoConverge', 'the bivariate gcd did not terminate.');
                end
                if numel(b) == 1
                    G = ratQ.contentFree(ratQ.cellToBiv({cg}));  return   % coprime in x
                end
                r = ratQ.bivPseudoRem(a, b);
                if isempty(r)
                    b = ratQ.bivScale(b, cg);
                    G = ratQ.contentFree(ratQ.cellToBiv(b));  return
                end
                a = b;
                b = ratQ.bivPrimitive(r);
            end
        end

        function c = bivToCell(A)
        % a bivariate array as a cell of y-polynomials, one per power of x, ASCENDING and trimmed.
            A = ratQ.trimN(A);
            c = cell(1, size(A,1));
            for i = 1:size(A,1), c{i} = ratQ.trim1(A(i,:)); end
            while numel(c) > 1 && all(c{end} == 0), c(end) = []; end
            if numel(c) == 1 && all(c{1} == 0), c = {}; end
        end

        function A = cellToBiv(c)
            if isempty(c), A = 0; return, end
            w = max(cellfun(@numel, c));
            A = zeros(numel(c), w);
            for i = 1:numel(c), A(i, 1:numel(c{i})) = c{i}; end
            A = ratQ.trimN(A);
        end

        function [c, cont] = bivPrimitive(c)
        % divide out the y-polynomial content -- the gcd of the x-coefficients.
            cont = [0];
            for i = 1:numel(c), cont = ratQ.yGcd(cont, c{i}); end
            if numel(cont) == 1 && abs(cont) == 1, cont = [1]; return, end
            for i = 1:numel(c), c{i} = ratQ.yDivExact(c{i}, cont); end
        end

        function c = bivScale(c, m)
            for i = 1:numel(c), c{i} = ratQ.trim1(conv(c{i}, m)); end
        end

        function r = bivPseudoRem(a, b)
        % objective: the pseudo-remainder of a by b -- both polynomials in x with y-polynomial
        %            coefficients -- reduced to its primitive part at EVERY step.
        %
        % The same discipline polyPremPrimitive applies one dimension down, and it is here for the
        % same measured reason: forming lc(b)^(deg a - deg b + 1) up front overflowed at 1.8e24 on
        % the bridge elimination's own coefficients. Multiplying by lc(b) ONCE per elimination and
        % taking the primitive part immediately keeps every intermediate bounded.
        %
        % No division is performed at all: the step is r := r*lc(b) - lc(r)*x^k*b, whose leading
        % term cancels by construction, so there is no exactness condition to check.
            db = numel(b) - 1;
            r  = a;
            guard = 0;
            while true
                while numel(r) > 1 && all(r{end} == 0), r(end) = []; end
                if isempty(r) || (numel(r) == 1 && all(r{1} == 0)), r = {}; return, end
                dr = numel(r) - 1;
                if dr < db, return, end
                guard = guard + 1;
                if guard > 10000
                    error('ratQ:premNoConverge', 'the bivariate pseudo-remainder did not terminate.');
                end
                lb = b{end};  lr = r{end};
                for i = 1:numel(r)
                    r{i} = ratQ.chk(conv(r{i}, lb), 'scaled coefficient');
                end
                sh = dr - db;
                for i = 0:db
                    t = ratQ.chk(conv(lr, b{i+1}), 'elimination term');
                    r{sh+i+1} = ratQ.trim1(ratQ.polyAdd1(r{sh+i+1}, -t));
                end
                r{end} = 0;                    % cancels exactly, by construction
                while numel(r) > 1 && all(r{end} == 0), r(end) = []; end
                if isempty(r) || (numel(r) == 1 && all(r{1} == 0)), r = {}; return, end
                r = ratQ.bivPrimitive(r);      % the growth control, and it must be INSIDE the loop
            end
        end
        function g = yGcd(a, b)
        % the primitive gcd of two ASCENDING integer y-polynomials, via ratQ.polyGcd, which works
        % in DESCENDING order -- the flip is the whole adapter and getting it wrong would silently
        % take the gcd of two reversed polynomials.
            a = ratQ.trim1(a);  b = ratQ.trim1(b);
            if all(a == 0), g = b; return, end
            if all(b == 0), g = a; return, end
            g = fliplr(ratQ.polyGcd(fliplr(a), fliplr(b)));
            g = ratQ.trim1(g);
        end

        function q = yDivExact(a, b)
        % exact division of ASCENDING integer y-polynomials; raises if it does not divide evenly.
            a = ratQ.trim1(a);  b = ratQ.trim1(b);
            if numel(b) == 1
                if b == 0, error('ratQ:divideByZero', 'y-division by zero.'); end
                if any(mod(a, b) ~= 0)
                    error('ratQ:notDivisible', 'the y-division is not exact.');
                end
                q = a / b;  return
            end
            [qq, rr] = deconv(fliplr(a), fliplr(b));
            if any(abs(ratQ.trim1(fliplr(rr))) > 0) || any(qq ~= round(qq))
                error('ratQ:notDivisible', 'the y-division is not exact.');
            end
            q = ratQ.trim1(fliplr(qq));
        end
        function g = contentOf(A)
        % objective: the gcd of an integer array's entries -- its CONTENT. Zero for an all-zero
        %            array, which is the identity for gcd and so composes correctly.
        %
        % Separate from contentFree because a caller sometimes needs the content of SEVERAL arrays
        % that must be scaled TOGETHER: the coefficients of one polynomial in a variable being
        % eliminated share a scale, and reducing them independently moves the polynomial's roots.
        % That is not hypothetical -- it made an elimination return a polynomial that ran and was
        % wrong, with the substituted root exactly half its true value.
            ratQ.chk(A, 'coefficient array');
            g = 0;
            for i = 1:numel(A)
                if A(i) ~= 0, g = gcd(g, abs(A(i))); end
                if g == 1, return, end
            end
        end

        function A = contentFree(A)
        % objective: divide an integer coefficient array by the gcd of its entries -- its PRIMITIVE
        %            PART. An all-zero array is returned unchanged.
        %
        % This is `conic`'s argument one dimension up, and it matters for the same reason plus one
        % more. A polynomial EQUATION has no scale: {p = 0} and {k*p = 0} are the same set. And in a
        % resultant the scale is not free -- a Sylvester determinant of NxN entries multiplies the
        % content roughly N times over, so carrying an avoidable factor of 4 into a 4x4 determinant
        % costs a factor of 256 against the 2^53 ceiling. Measured: leaving C unreduced took the
        % bridge elimination from 2.7e15 to an overflow at 1.3e16.
            ratQ.chk(A, 'coefficient array');
            g = 0;
            for i = 1:numel(A)
                if A(i) ~= 0, g = gcd(g, abs(A(i))); end
                if g == 1, return, end
            end
            if g > 1, A = A / g; end
        end

        function A = trimN(A)
        % objective: drop trailing all-zero slices in every dimension, so an array's size states
        %            its actual degrees. Always leaves at least a scalar.
            for d = 1:ndims(A)
                while size(A,d) > 1
                    idx = repmat({':'}, 1, ndims(A));
                    idx{d} = size(A,d);
                    if any(A(idx{:}) ~= 0, 'all'), break, end
                    idx{d} = 1:size(A,d)-1;
                    A = A(idx{:});
                end
            end
        end

        function A = squeezeLast(A, nv)
        % objective: drop the trailing singleton left by slicing dimension nv, without disturbing
        %            the leading dimensions -- `squeeze` would also collapse a genuine singleton in
        %            the middle, which changes which variable is which.
            sz = arrayfun(@(d) size(A,d), 1:nv-1);
            if isempty(sz), sz = 1; end
            if numel(sz) == 1, sz = [sz 1]; end
            A = reshape(A, sz);
        end

        % ---- THE BIG-INTEGER FALLBACK ------------------------------------------------------------
        % Entered only where the double path would RAISE, which is the arrangement DECISIONS.md
        % 2026-09-05 chose after measuring the options: java.math.BigInteger is 121x doubles but
        % needs no compiler, no licence and no network, and MATLAB ships the JVM. The hot path keeps
        % its double arithmetic and its exactness check; nothing above ratQ changes.
        %
        % WHY IT IS THE POLYNOMIAL LAYER AND NOT ratQ AS A WHOLE. Every overflow measured on the
        % convex envelope's elimination is in a polynomial PRODUCT or a determinant built from
        % them -- and in each the intermediates are large while the ANSWER is small (the isotropic
        % bridge reduces to coefficients of 486). So the useful shape is: compute big, reduce, come
        % back to doubles, and raise only if the REDUCED result is still too large.
        %
        % A big polynomial is a struct with `sz` (the size vector) and `c` (a cell of BigInteger in
        % linear order). Cells rather than a javaArray because MATLAB indexes them naturally, and
        % the arrays here are sparse enough that the loops skip most entries.

        function B = toBig(A)
        % objective: promote a double coefficient array to the big representation. Exact.
            ratQ.chk(A, 'value to promote');
            % ZERO IS THE EMPTY CELL, not a BigInteger.ZERO object. Measured: filling a result
            % array with Java zero-references cost 10 ms per add or multiply and dominated the whole
            % elimination, because a 5-variable result has tens of thousands of slots and almost all
            % of them are zero. cell(n,1) allocates them empty for free.
            B.sz = size(A);
            B.c  = cell(numel(A), 1);
            for i = 1:numel(A)
                if A(i) ~= 0, B.c{i} = java.math.BigInteger.valueOf(int64(A(i))); end
            end
        end

        function A = fromBigOrKeep(B)
        % objective: come back to doubles if the value fits, and STAY BIG if it does not.
        % For the middle of a chain of eliminations, where an intermediate may not fit but the
        % final answer will: converting eagerly would raise on a value that was never meant to be
        % returned, and staying big always would make every later step pay the 121x.
            try
                A = ratQ.fromBig(B);
            catch err
                if ~strcmp(err.identifier, 'ratQ:overflow'), rethrow(err); end
                A = B;
            end
        end

        function A = fromBig(B)
        % objective: come back to doubles, RAISING if any coefficient no longer fits exactly.
        % That raise is the honest outcome: it means the answer itself, not merely the path to it,
        % is past what a double can hold.
            A = zeros(B.sz);
            for i = 1:numel(B.c)
                v = B.c{i};
                if isempty(v) || v.signum() == 0, continue, end
                if v.abs().bitLength() > 52
                    error('ratQ:overflow', ...
                        ['a coefficient needs %d bits, past the 53 a double holds exactly. The ' ...
                         'big path computed it; the RESULT itself is too large to return.'], ...
                        v.abs().bitLength());
                end
                A(i) = double(v.doubleValue());
            end
        end

        function C = bigMul(A, B)
        % objective: multiply two big polynomials -- dense convolution, skipping zeros.
            if isempty(A.c) || isempty(B.c), C = A; return, end
            nd = max(numel(A.sz), numel(B.sz));
            sa = [A.sz, ones(1, nd - numel(A.sz))];
            sb = [B.sz, ones(1, nd - numel(B.sz))];
            sc = sa + sb - 1;
            C.sz = sc;
            C.c  = cell(prod(sc), 1);
            ia = find(~cellfun(@isempty, A.c));
            ib = find(~cellfun(@isempty, B.c));
            for p1 = ia(:).'
                ea = ratQ.linToSub(p1, sa);
                for p2 = ib(:).'
                    eb = ratQ.linToSub(p2, sb);
                    k  = ratQ.subToLin(ea + eb - 1, sc);
                    t  = A.c{p1}.multiply(B.c{p2});
                    if isempty(C.c{k}), C.c{k} = t; else, C.c{k} = C.c{k}.add(t); end
                    if C.c{k}.signum() == 0, C.c{k} = []; end
                end
            end
        end

        function C = bigAdd(A, B)
        % objective: add two big polynomials, zero-padding to the union of their shapes.
            nd = max(numel(A.sz), numel(B.sz));
            sa = [A.sz, ones(1, nd - numel(A.sz))];
            sb = [B.sz, ones(1, nd - numel(B.sz))];
            sc = max(sa, sb);
            C.sz = sc;
            C.c  = cell(prod(sc), 1);
            for i = 1:numel(A.c)
                if isempty(A.c{i}), continue, end
                k = ratQ.subToLin(ratQ.linToSub(i, sa), sc);
                if isempty(C.c{k}), C.c{k} = A.c{i}; else, C.c{k} = C.c{k}.add(A.c{i}); end
                if C.c{k}.signum() == 0, C.c{k} = []; end
            end
            for i = 1:numel(B.c)
                if isempty(B.c{i}), continue, end
                k = ratQ.subToLin(ratQ.linToSub(i, sb), sc);
                if isempty(C.c{k}), C.c{k} = B.c{i}; else, C.c{k} = C.c{k}.add(B.c{i}); end
                if C.c{k}.signum() == 0, C.c{k} = []; end
            end
        end

        function A = bigNeg(A)
            for i = 1:numel(A.c)
                if ~isempty(A.c{i}), A.c{i} = A.c{i}.negate(); end
            end
        end

        function A = bigContentFree(A)
        % objective: divide a big polynomial by the gcd of its coefficients.
            g = java.math.BigInteger.ZERO;
            for i = 1:numel(A.c)
                if ~isempty(A.c{i}), g = g.gcd(A.c{i}); end
            end
            if g.signum() == 0 || g.equals(java.math.BigInteger.ONE), return, end
            for i = 1:numel(A.c)
                if ~isempty(A.c{i}), A.c{i} = A.c{i}.divide(g); end
            end
        end

        function e = linToSub(i, sz)
            e = zeros(1, numel(sz));
            i = i - 1;
            for d = 1:numel(sz)
                e(d) = mod(i, sz(d)) + 1;
                i = floor(i / sz(d));
            end
        end

        function i = subToLin(e, sz)
            i = 0;  m = 1;
            for d = 1:numel(sz)
                i = i + (e(d) - 1) * m;
                m = m * sz(d);
            end
            i = i + 1;
        end

        function B = bigMono(e, c, nv) %#ok<INUSD>
        % objective: the single big monomial c * v1^e1 ... vnv^env.
            sz = e + 1;
            if numel(sz) < 2, sz = [sz 1]; end
            B.sz = sz;
            B.c  = cell(prod(sz), 1);
            B.c{ratQ.subToLin(e + 1, sz)} = c;
        end

        function [e, c] = bigLeadTerm(A, nv)
        % objective: the LEXICOGRAPHICALLY largest exponent carrying a nonzero coefficient, and
        %            that coefficient. Empty when A is zero. First variable most significant, so
        %            that bigDivExact terminates -- the same order leadTerm uses.
            e = [];  c = java.math.BigInteger.ZERO;
            sz = [A.sz, ones(1, nv - numel(A.sz))];
            best = [];
            % EMPTY *OR* AN EXPLICIT ZERO. bigAdd normalises a cancelling sum back to empty, and
            % this is the reader that depends on it: a BigInteger.ZERO left in a cell would be read
            % as a term, the leading exponent would stop decreasing, and bigDivExact would spin --
            % which is exactly what it did until the invariant was restored.
            for i = 1:numel(A.c)
                if isempty(A.c{i}) || A.c{i}.signum() == 0, continue, end
                ee = ratQ.linToSub(i, sz) - 1;
                if isempty(best) || ratQ.lexGreater(ee, best)
                    best = ee;  c = A.c{i};
                end
            end
            e = best;
        end

        function tf = lexGreater(a, b)
            tf = false;
            for d = 1:numel(a)
                if a(d) > b(d), tf = true; return, end
                if a(d) < b(d), return, end
            end
        end

        function Q = bigDivExact(P, D, nv)
        % objective: divide big polynomial P by D, EXACTLY, raising if it does not divide evenly.
        % The BigInteger counterpart of polyDivExact, and it exists for the same reason the big
        % determinant does: on the anisotropic bridge the intermediates need 59 bits while the
        % reduced answer needs far fewer, so the reduction has to happen BEFORE coming back to
        % doubles.
            [eD, cD] = ratQ.bigLeadTerm(D, nv);
            if isempty(eD)
                error('ratQ:divideByZero', 'the divisor is the zero polynomial.');
            end
            R = P;
            Q = struct('sz', ones(1, max(2,nv)), 'c', {{[]}});
            guard = 0;
            while true
                [eR, cR] = ratQ.bigLeadTerm(R, nv);
                if isempty(eR), break, end
                guard = guard + 1;
                if guard > 100000
                    error('ratQ:divideNoConverge', 'big exact division did not terminate.');
                end
                eQ = eR - eD;
                qr = cR.divideAndRemainder(cD);
                if any(eQ < 0) || qr(2).signum() ~= 0
                    error('ratQ:notDivisible', ...
                        'the big division is not exact; nothing is truncated or rounded here.');
                end
                T = ratQ.bigMono(eQ, qr(1), nv);
                Q = ratQ.bigAdd(Q, T);
                R = ratQ.bigAdd(R, ratQ.bigNeg(ratQ.bigMul(T, D)));
            end
        end

        function R = resultantVarBig(P, Q, k, nv)
        % objective: resultantVar returning the BIG result, without insisting it fit in a double.
        % For callers that will REDUCE the answer before they need it as doubles -- which is the
        % shape every use on the envelope line has.
            R = ratQ.resultantCore(P, Q, k, nv, true);
        end

        function tf = bigIsZero(A)
            tf = true;
            for i = 1:numel(A.c)
                if ~isempty(A.c{i}) && A.c{i}.signum() ~= 0, tf = false; return, end
            end
        end

        function d = polyDetBig(S)
        % objective: polyDetN in BigInteger. Same cofactor expansion, same recursion, no ceiling.
            n = size(S,1);
            if n == 1, d = S{1,1}; return, end
            d = struct('sz', [1 1], 'c', {{[]}});
            for j = 1:n
                if all(cellfun(@isempty, S{1,j}.c)), continue, end
                M = S(2:n, [1:j-1, j+1:n]);
                t = ratQ.bigMul(S{1,j}, ratQ.polyDetBig(M));
                if mod(j,2) == 0, t = ratQ.bigNeg(t); end
                d = ratQ.bigAdd(d, t);
            end
        end

        function d = polyDetN(S)
        % objective: the determinant of a square cell array whose entries are integer polynomial
        %            ARRAYS, by cofactor expansion along the first row. Exact.
        % The N-dimensional counterpart of polyDet; multiplication is convn, which IS dense
        % polynomial multiplication in any number of variables.
            n = size(S,1);
            if n == 1, d = S{1,1}; return, end
            d = 0;
            for j = 1:n
                if isscalar(S{1,j}) && S{1,j} == 0, continue, end
                M = S(2:n, [1:j-1, j+1:n]);
                t = ratQ.polyMulN(S{1,j}, ratQ.polyDetN(M));
                if mod(j,2) == 0, t = -t; end
                d = ratQ.polyAddN(d, t);
            end
        end

        function c = polyMulN(a, b)
        % objective: multiply two integer polynomials held as dense coefficient ARRAYS, in any
        %            number of variables. convn is exactly that product.
            if (isscalar(a) && a == 0) || (isscalar(b) && b == 0), c = 0; return, end
            c = ratQ.chk(convn(a, b), 'polynomial product');
        end

        function c = polyAddN(a, b)
        % objective: add two integer polynomial arrays of any sizes, zero-padding to the union.
            if isscalar(a) && a == 0, c = b; return, end
            if isscalar(b) && b == 0, c = a; return, end
            na = ndims(a);  nb = ndims(b);  nd = max(na, nb);
            sa = arrayfun(@(d) size(a,d), 1:nd);
            sb = arrayfun(@(d) size(b,d), 1:nd);
            sc = max(sa, sb);
            A = zeros(sc);  B = zeros(sc);
            ia = arrayfun(@(d) 1:sa(d), 1:nd, 'UniformOutput', false);
            ib = arrayfun(@(d) 1:sb(d), 1:nd, 'UniformOutput', false);
            A(ia{:}) = a;  B(ib{:}) = b;
            c = A + B;
        end

        function R = resultantY(P, Q)
        % objective: the SYLVESTER RESULTANT of two bivariate integer polynomials with respect to
        %            y, exactly. R vanishes at exactly those x for which P(x,.) and Q(x,.) have a
        %            common root in y -- which is what ELIMINATION means.
        %
        % [input]  P, Q : integer matrices with P(i+1,j+1) the coefficient of x^i y^j. The two need
        %                 not have the same size; trailing zero rows/columns are trimmed.
        % [output] R    : 1 x k integer row, ASCENDING in x (R(i+1) is the coefficient of x^i).
        %
        % WHY THIS EXISTS. The convex envelope's bridge face is defined by a system in the dual
        % variable s that has no closed-form solution: the tie conic and the collinearity conic meet
        % in up to four points, and the face is the relation between x and the value that survives
        % eliminating s. `.claude/bridgeEliminationSpike.m` establishes that both conics are exact
        % over Z; this is the routine that turns them into that relation. Nothing else in ratQ can:
        % polyGcd decides whether two UNIVARIATE polynomials share a root, and the question here is
        % for WHICH x two bivariate ones do.
        %
        % HOW. Written as polynomials in y whose coefficients are polynomials in x, the two have a
        % common root in y exactly when the (m+n) x (m+n) Sylvester matrix is singular. Its entries
        % are x-polynomials, so its determinant is an x-polynomial, and that determinant IS the
        % resultant. Expanded by cofactors, which uses only integer multiply and add -- so the
        % result is exact and no division, exact or otherwise, is ever performed.
        %
        % Cofactor expansion is factorial in the matrix size and that is deliberate: the sizes here
        % are 4 (two conics) and the alternative, fraction-free Bareiss elimination, needs exact
        % polynomial division, which is a second thing to get right for no benefit at n = 4. The
        % size is capped rather than left to grow silently.
            ratQ.chk(P, 'first polynomial');  ratQ.chk(Q, 'second polynomial');
            P = ratQ.trim2(P);  Q = ratQ.trim2(Q);
            m = size(P,2) - 1;              % degree in y
            n = size(Q,2) - 1;
            if m < 1 || n < 1
                error('ratQ:resultantDegree', ...
                    ['resultantY eliminates y, so both arguments must have degree >= 1 in it; ' ...
                     'got %d and %d. A polynomial free of y has no root to share.'], m, n);
            end
            N = m + n;
            if N > 8
                error('ratQ:resultantTooLarge', ...
                    ['the Sylvester matrix is %dx%d, and cofactor expansion is factorial in that. ' ...
                     'Above 8 this wants Bareiss elimination instead.'], N, N);
            end

            % Sylvester matrix, entries are x-polynomials held as ASCENDING rows in a cell array.
            S = cell(N, N);
            [S{:}] = deal(0);
            for r = 1:n                      % n shifted copies of P's y-coefficients, high to low
                for j = 0:m
                    S{r, r+j} = P(:, m-j+1).';
                end
            end
            for r = 1:m                      % m shifted copies of Q's
                for j = 0:n
                    S{n+r, r+j} = Q(:, n-j+1).';
                end
            end
            R = ratQ.polyDet(S);
            R = ratQ.trim1(R);
            R = ratQ.chk(R, 'resultant');
        end

        function d = polyDet(S)
        % objective: the determinant of a square cell array whose entries are integer polynomials
        %            (ASCENDING rows), by cofactor expansion along the first row. Exact.
            n = size(S,1);
            if n == 1, d = S{1,1}; return, end
            d = 0;
            for j = 1:n
                if isscalar(S{1,j}) && S{1,j} == 0, continue, end
                M = S([2:n], [1:j-1, j+1:n]); %#ok<NBRAK>
                t = ratQ.polyMul1(S{1,j}, ratQ.polyDet(M));
                if mod(j,2) == 0, t = -t; end
                d = ratQ.polyAdd1(d, t);
            end
        end

        function c = polyMul1(a, b)
        % objective: multiply two UNIVARIATE integer polynomials held as ascending rows. Exact.
        % Separate from polyMul, which is the BIVARIATE cubic-basis product and caps total degree.
            if (isscalar(a) && a == 0) || (isscalar(b) && b == 0), c = 0; return, end
            c = ratQ.chk(conv(a(:).', b(:).'), 'polynomial product');
        end

        function c = polyAdd1(a, b)
        % objective: add two UNIVARIATE integer polynomials held as ascending rows, of any lengths.
            a = a(:).';  b = b(:).';
            k = max(numel(a), numel(b));
            c = [a, zeros(1, k-numel(a))] + [b, zeros(1, k-numel(b))];
        end

        function C = polyMulBiv(A, B)
        % objective: multiply two bivariate integer polynomials with NO degree cap, in the same
        %            (i+1,j+1) = coeff of x^i y^j layout resultantY uses.
        %
        % Separate from polyMul, which is the CUBIC-BASIS product and raises above total degree 3 --
        % correct for substAffine's fixed basis, and useless for building resultant test fixtures.
            ratQ.chk(A, 'first factor');  ratQ.chk(B, 'second factor');
            C = zeros(size(A,1)+size(B,1)-1, size(A,2)+size(B,2)-1);
            for i1 = 1:size(A,1)
                for j1 = 1:size(A,2)
                    if A(i1,j1) == 0, continue, end
                    for i2 = 1:size(B,1)
                        for j2 = 1:size(B,2)
                            if B(i2,j2) == 0, continue, end
                            C(i1+i2-1, j1+j2-1) = C(i1+i2-1, j1+j2-1) + A(i1,j1)*B(i2,j2);
                        end
                    end
                end
            end
            C = ratQ.chk(C, 'bivariate product');
        end

        function A = trim2(A)
        % objective: drop trailing all-zero rows and columns of a bivariate coefficient matrix, so
        %            that its size states its actual degrees. Always leaves at least a 1x1.
            while size(A,1) > 1 && all(A(end,:) == 0), A(end,:) = []; end
            while size(A,2) > 1 && all(A(:,end) == 0), A(:,end) = []; end
        end

        function a = trim1(a)
        % objective: the same for a univariate ascending row.
            a = a(:).';
            while numel(a) > 1 && a(end) == 0, a(end) = []; end
        end

        function Ps = polyPow(A, kmax)
        % objective: the powers A^0 ... A^kmax of a bivariate polynomial, as a cell array with
        %            Ps{k+1} = A^k. A^0 is the constant 1.
            Ps = cell(1, kmax+1);
            Ps{1} = zeros(4,4);  Ps{1}(1,1) = 1;
            for k = 1:kmax
                Ps{k+1} = ratQ.polyMul(Ps{k}, A);
            end
        end

        function tf = feasible2(P, strict)
        % objective: is the 2-D polyhedron { s : P(j,1)*s1 + P(j,2)*s2 + P(j,3) >= 0 for all j }
        %            nonempty -- and with `strict`, does it have nonempty INTERIOR.
        % [input]  P      : k x 3 integer rows [p1 p2 p0]; strict : (optional) default false
        % [output] tf     : logical, DECIDED exactly
        %
        % WHY THIS IS NEEDED. A cell of a conjugate's subdivision is a list of sign conditions on
        % lines (the H-form), and some of those lists describe the empty set -- a polygon vertex
        % that never attains the max contributes a cell that is empty, not small. Emitting it
        % anyway is not wrong for `eval` (no point satisfies it) but it inflates the face count and
        % puts a face into the mesh that bounds nothing. `strict` separates the other degenerate
        % case: a cell that is a single point or a segment is nonempty and still carries no
        % two-dimensional face.
        %
        % FOURIER-MOTZKIN, and it is the right algorithm at n = 2 for the same reason Cramer is at
        % n = 3: eliminating one of two variables leaves a one-dimensional problem that is just
        % "is the largest lower bound below the smallest upper bound", and every comparison is
        % between two rationals, i.e. one integer cross-multiplication. No pivoting, no ordering
        % choice, nothing to get wrong. The pair blow-up that makes Fourier-Motzkin unusable in
        % general is bounded here by one elimination.
        %
        % ALL COMPARISONS ARE EXACT. -b/a <= -d/c with a, c > 0 is b*c >= d*a, and both sides go
        % through ratQ.chk, so an overflow raises rather than deciding the wrong way.
            if nargin < 2, strict = false; end
            ratQ.chk(P, 'half-plane');
            if isempty(P), tf = true; return, end

            % ---- eliminate s2 -------------------------------------------------------------
            % Rows with p2 = 0 already constrain s1 alone; the rest split into lower and upper
            % bounds on s2 and every (lower, upper) pair yields one more constraint on s1.
            Q = zeros(0,2);                          % rows [a b] meaning a*s1 + b >= 0
            lo = find(P(:,2) > 0);  up = find(P(:,2) < 0);  fl = find(P(:,2) == 0);
            for i = fl(:).'
                Q(end+1,:) = [P(i,1), P(i,3)]; %#ok<AGROW>
            end
            for i = lo(:).'
                for j = up(:).'
                    % P(i): p2i*s2 >= -(p1i*s1 + p0i)   with p2i > 0
                    % P(j): p2j*s2 >= -(p1j*s1 + p0j)   with p2j < 0, i.e. an upper bound
                    % Combining eliminates s2 with positive weights (-p2j) and (p2i):
                    a = ratQ.chk(-P(j,2)*P(i,1) + P(i,2)*P(j,1), 'FM coefficient');
                    b = ratQ.chk(-P(j,2)*P(i,3) + P(i,2)*P(j,3), 'FM constant');
                    Q(end+1,:) = [a, b]; %#ok<AGROW>
                end
            end

            % ---- solve the 1-D problem ------------------------------------------------------
            % a*s1 + b >= 0 is s1 >= -b/a for a > 0, s1 <= -b/a for a < 0, and b >= 0 for a = 0.
            % Bounds are kept as the pair (num, den) with den > 0 and compared by cross product.
            loN = -inf; loD = 1;  hiN = inf; hiD = 1;  haveLo = false; haveHi = false;
            for r = 1:size(Q,1)
                a = Q(r,1);  b = Q(r,2);
                if a == 0
                    if (strict && b <= 0) || (~strict && b < 0), tf = false; return, end
                elseif a > 0
                    % s1 >= -b/a, and a > 0, so the bound is (-b)/a. Keep the LARGEST.
                    if ~haveLo || ratQ.chk(-b*loD - loN*a, 'lower bound') > 0
                        loN = -b;  loD = a;  haveLo = true;
                    end
                else
                    % s1 <= -b/a with a < 0; writing it with a positive denominator gives
                    % b/(-a). Keep the SMALLEST.
                    nN = b;  nD = -a;
                    if ~haveHi || ratQ.chk(nN*hiD - hiN*nD, 'upper bound') < 0
                        hiN = nN;  hiD = nD;  haveHi = true;
                    end
                end
            end
            if haveLo && haveHi
                c = ratQ.chk(loN*hiD - hiN*loD, 'bound comparison');
                if (strict && c >= 0) || (~strict && c > 0), tf = false; return, end
            end
            tf = true;
        end

        function s = conicSign(c)
        % objective: does the conic form c take ONE sign over the whole plane, and which.
        % [input]  c : 1 x 6 integer conic [a b c d e f] for a x^2 + b xy + c y^2 + d x + e y + f
        % [output] s : +1 if c(x,y) >= 0 everywhere, -1 if <= 0 everywhere, 0 if it takes both signs
        %
        % WHY THIS IS WORTH A ROUTINE. A cell of the subdivision is a list of sign conditions, and
        % a nested fold produces conditions that are vacuous (the curve has no real points, so one
        % side is everything) or contradictory (so the cell has empty interior). Both are invisible
        % to a LINEAR feasibility test, and both are decided here exactly, in integers -- which is
        % the only handle on curved emptiness that does not need the degree-4 kernel.
        %
        % THE TEST. Writing m = [x; y; 1], the form is c(x,y) = (1/2) m' M m with
        %       M = [2a  b  d;  b  2c  e;  d  e  2f]
        % (the halves cleared by the factor 2, so M is integral). c >= 0 on the whole plane exactly
        % when M is positive SEMIdefinite -- the m = (x,y,1) slice misses only the points at
        % infinity, where a PSD form is still nonnegative by continuity.
        %
        % ALL SEVEN PRINCIPAL MINORS, not the three leading ones. Sylvester's leading-minor
        % criterion characterises positive DEFINITE; for semidefinite it is false in both
        % directions -- diag(0,-1) has both leading minors 0 and is not PSD. So every principal
        % minor is tested, which for 3x3 is three diagonal entries, three 2x2s and the determinant.
            ratQ.chk(c, 'conic');
            M = [2*c(1),   c(2),   c(4); ...
                   c(2), 2*c(3),   c(5); ...
                   c(4),   c(5), 2*c(6)];
            if ratQ.isPSD3(M)
                s = 1;
            elseif ratQ.isPSD3(-M)
                s = -1;
            else
                s = 0;
            end
        end

        function tf = isPSD2(H)
        % objective: is the 2x2 symmetric integer matrix H positive semidefinite. Exact.
        %
        % For 2x2 the principal minors are the two diagonal entries and the determinant, and all
        % three must be nonnegative. The LEADING minors alone would not do: diag(0,-1) has leading
        % minors 0 and 0 and is not PSD -- the same trap isPSD3 documents.
            ratQ.chk(H, 'matrix');
            tf = H(1,1) >= 0 && H(2,2) >= 0 && ratQ.detExact(H) >= 0;
        end

        function tf = isPSD3(M)
        % objective: is the 3x3 symmetric integer matrix M positive semidefinite. Exact.
        % Every principal minor must be nonnegative -- see conicSign for why the leading ones alone
        % are not enough.
            ratQ.chk(M, 'matrix');
            if any(diag(M) < 0), tf = false; return, end
            pairs = [1 2; 1 3; 2 3];
            for k = 1:3
                p = pairs(k,:);
                if ratQ.detExact(M(p,p)) < 0, tf = false; return, end
            end
            tf = ratQ.detExact(M) >= 0;
        end

        function chain = sturmChain(p)
        % objective: the STURM CHAIN of an integer polynomial -- p, p', and the NEGATED remainders.
        % [input]  p : 1 x (n+1) integer coefficients, highest degree first
        % [output] chain : cell array of integer polynomials
        %
        % Sturm's theorem: for a SQUAREFREE p, the number of distinct real roots in (a,b] is
        % V(a) - V(b), where V counts sign changes along the chain. That is an exact statement about
        % integers, which is the whole reason this exists -- it decides how many roots are where,
        % with no root finder and no tolerance.
        %
        % PRIMITIVE PART AT EVERY STEP, and it is not an optimisation. The pseudo-remainder squares
        % its leading coefficient each round, so a quartic with six-digit entries reaches 1e23 in a
        % few steps and ratQ.chk raises -- measured in conicMeet's own squarefreePart, which
        % records the same growth. Dividing out the content keeps the chain integral AND small, and
        % it changes no SIGN, because the content is positive by construction.
            p = ratQ.polyTrim(p);
            p = ratQ.polySquarefree(p);
            chain = {p};
            if numel(p) <= 1, return, end
            chain{2} = ratQ.polyTrim(polyder(p));
            while numel(chain{end}) > 1
                r = ratQ.polyPseudoRem(chain{end-1}, chain{end});
                r = ratQ.polyTrim(r);
                if isempty(r) || all(r == 0), break, end
                chain{end+1} = ratQ.polyPrimitive(-r); %#ok<AGROW>
            end
        end

        function n = countRootsIn(p, a, b)
        % objective: the number of DISTINCT real roots of p in the half-open interval (a, b].
        % [input]  p : integer polynomial; a, b : rational endpoints (doubles that are exact)
        % [output] n : a nonnegative integer, DECIDED
        %
        % Sturm's theorem, applied to the chain above. A repeated root counts ONCE, because the
        % chain is built from the squarefree part -- which is what the callers want: two coincident
        % intersections of two conics are one vertex, not two.
            chain = ratQ.sturmChain(p);
            n = ratQ.sturmV(chain, a) - ratQ.sturmV(chain, b);
        end

        function v = sturmV(chain, x)
        % objective: the number of SIGN CHANGES in the chain evaluated at x, zeros skipped.
        % Evaluation is in double precision on purpose and it is safe: the count only changes when
        % a chain member's sign changes, and the callers only ever evaluate at endpoints they chose
        % -- integers or dyadic rationals from bisection -- where every chain member is bounded away
        % from zero unless it vanishes exactly, which skipping handles.
            v = 0;  last = 0;
            for i = 1:numel(chain)
                s = sign(polyval(chain{i}, x));
                if s == 0, continue, end
                if last ~= 0 && s ~= last, v = v + 1; end
                last = s;
            end
        end

        function I = isolateRoots(p)
        % objective: one interval per distinct real root of p, each containing exactly that root.
        % [output] I : k x 2, ascending, with I(j,1) <= root_j <= I(j,2)
        %
        % BISECTION DRIVEN BY THE EXACT COUNT. Start from a bound that provably contains every real
        % root (Cauchy's: 1 + max|a_i|/|a_n|), then split any interval holding more than one root
        % until each holds exactly one. Sturm decides "how many are in here" at every step, so the
        % recursion is driven by an exact integer and cannot be fooled by two roots that are close.
            p = ratQ.polySquarefree(ratQ.polyTrim(p));
            if numel(p) <= 1, I = zeros(0,2); return, end
            M = 1 + max(abs(p(2:end))) / abs(p(1));
            M = ceil(M);
            chain = ratQ.sturmChain(p);
            I = bisect(chain, -M, M, ratQ.sturmV(chain,-M) - ratQ.sturmV(chain,M), 0);
            if ~isempty(I), I = sortrows(I); end
        end

        function s = signAt(q, p, iv)
        % objective: the EXACT sign of q at the root of p isolated in the interval iv.
        % [input]  q : integer polynomial; p : integer polynomial; iv : [lo hi] from isolateRoots
        % [output] s : -1, 0 or +1
        %
        % THE PREDICATE EVERYTHING ELSE IS BUILT ON. A vertex of f* is a root of an exact quartic
        % (conicMeet returns it), and every mesh question -- is this vertex inside that face, which
        % of two crossings comes first, do these two cells share a facet -- is the sign of a
        % rational polynomial AT such a number. CONJ_FIELD_PROOF.md 8.0 lists exactly these.
        %
        % TWO STEPS, and the first is what makes the second terminate:
        %
        %   1. IS q(alpha) ZERO? Exactly when p and q share the root alpha, i.e. when gcd(p,q) has a
        %      root in the interval -- an exact Sturm count on the gcd, no evaluation at all. This
        %      has to come first: without it, step 2 would bisect forever trying to separate alpha
        %      from a zero of q that IS alpha.
        %   2. OTHERWISE q has NO root at alpha, so refining the interval until it contains no root
        %      of q at all leaves q of one sign throughout -- and then evaluating at either endpoint
        %      gives that sign. Sturm again decides "does q have a root in here", so the loop is
        %      driven by exact integers and terminates because alpha is isolated from q's roots by a
        %      positive distance.
            p = ratQ.polySquarefree(ratQ.polyTrim(p));
            q = ratQ.polyTrim(q);
            if isempty(q) || all(q == 0), s = 0; return, end

            g = ratQ.polyGcd(p, q);
            if numel(g) > 1 && ratQ.countRootsIn(g, iv(1), iv(2)) > 0
                s = 0;  return                      % alpha is a common root: q(alpha) = 0 exactly
            end

            chainP = ratQ.sturmChain(p);
            lo = iv(1);  hi = iv(2);
            for it = 1:200
                if ratQ.countRootsIn(q, lo, hi) == 0
                    % q has no root in [lo,hi], so it holds one sign there; alpha is inside.
                    v = polyval(q, (lo + hi)/2);
                    if v ~= 0, s = sign(v); return, end
                    v = polyval(q, hi);
                    if v ~= 0, s = sign(v); return, end
                    s = sign(polyval(q, lo));  return
                end
                mid = (lo + hi) / 2;
                if ratQ.sturmV(chainP, lo) - ratQ.sturmV(chainP, mid) == 1
                    hi = mid;
                else
                    lo = mid;
                end
            end
            error('ratQ:signAtNoConverge', ...
                ['bisection did not separate the root from q''s zeros in 200 steps. The root and ' ...
                 'a zero of q are then closer than a dyadic rational of that depth, which for ' ...
                 'exact integer input means they are equal -- so the gcd test above should have ' ...
                 'caught it, and this is a defect rather than a precision limit.']);
        end

        function p = polyTrim(p)
        % objective: drop leading zeros, so a polynomial has one spelling. [] is the zero polynomial.
            if isempty(p), return, end
            k = find(p ~= 0, 1);
            if isempty(k), p = []; else, p = p(k:end); end
        end

        function p = polyPrimitive(p)
        % objective: divide out the integer content, keeping the SIGN of the leading coefficient.
        % The content is positive, so this changes no sign anywhere -- which is what lets the Sturm
        % chain be reduced at every step without disturbing the count.
            p = ratQ.polyTrim(p);
            if isempty(p), return, end
            g = 0;
            for i = 1:numel(p)
                if p(i) ~= 0, g = gcd(g, abs(p(i))); end
            end
            if g > 1, p = p / g; end
            ratQ.chk(p, 'primitive polynomial');
        end

        function r = polyPremPrimitive(a, b)
        % objective: the PRIMITIVE PART of the pseudo-remainder of a by b, computed WITHOUT ever
        %            forming the full lc(b)^(deg a - deg b + 1) scaling. Exact, descending order.
        %
        % WHY IT EXISTS. polyPseudoRem multiplies the whole dividend by lc(b)^(da-db+1) up front,
        % and that is what overflows: measured on the convex envelope's bridge elimination, the
        % coefficients were already 2.7e12 and the scaling squared them to 3.3e16, past 2^53.
        % Taking the primitive part BETWEEN Euclid steps -- which polyGcd already did -- cannot
        % help, because the overflow happens WITHIN the single scaling.
        %
        % The fix is to interleave: multiply by lc(b) ONCE per elimination step, cancel the leading
        % term, and take the primitive part immediately. Each step is exact for the same reason the
        % full scaling is (one factor of lc(b) is exactly what one elimination needs), and the
        % content never gets a chance to accumulate. The result differs from polyPseudoRem's only
        % by a content factor, which is why polyGcd -- whose answer is primitive by construction --
        % can use either.
            a = ratQ.polyPrimitive(a);  b = ratQ.polyPrimitive(b);
            if isempty(b)
                error('ratQ:divideByZero', 'pseudo-remainder by the zero polynomial.');
            end
            db = numel(b) - 1;
            r  = a;
            guard = 0;
            while true
                r = ratQ.polyTrim(r);
                if isempty(r) || numel(r) - 1 < db, return, end
                guard = guard + 1;
                if guard > 10000
                    error('ratQ:premNoConverge', 'the reduced pseudo-remainder did not terminate.');
                end
                lcb = b(1);  lcr = r(1);
                r = ratQ.chk(r * lcb, 'reduced pseudo-remainder step');
                r(1:numel(b)) = ratQ.chk(r(1:numel(b)) - lcr * b, 'elimination');
                r(1) = 0;                      % the leading term cancels exactly, by construction
                r = ratQ.polyPrimitive(ratQ.polyTrim(r));
            end
        end

        function r = polyPseudoRemBig(a, b)
        % objective: the pseudo-remainder of a by b, computed in java.math.BigInteger and returned
        %            as DOUBLES once its content is removed.
        % [input]  a, b : integer coefficient vectors, highest degree first
        % [output] r    : the reduced pseudo-remainder, integral and within 2^53
        %
        % WHY THE BIG PATH IS ONLY THE MIDDLE. The overflow is never in the ANSWER, it is in the
        % scaling: lc(b)^(deg a - deg b + 1) times a, which for the region-membership quartic
        % reaches 3.7e18 while the reduced remainder is small. So this computes the whole
        % elimination in BigInteger, divides out the content THERE, and only then comes back to
        % doubles -- the same "try doubles, fall back, reduce, return" shape ratQ.resultantVar
        % already uses for the Sylvester determinant, and for the same reason.
        %
        % IF THE REDUCED RESULT STILL DOES NOT FIT, this raises rather than returning a big form.
        % The Sturm chain is consumed by sturmV, which evaluates in doubles, so a chain entry that
        % is not a double is not usable by the thing that asked for it; a named refusal is the
        % honest answer and no caller is left holding a value it cannot read.
            a = ratQ.polyPrimitive(a);  b = ratQ.polyPrimitive(b);
            if isempty(b), error('ratQ:divideByZero', 'pseudo-remainder by the zero polynomial.'); end
            da = numel(a) - 1;  db = numel(b) - 1;
            if da < db, r = a; return, end

            A  = ratQ.bigVec(a);
            Bv = ratQ.bigVec(b);
            lc = Bv{1};
            f  = lc.pow(da - db + 1);
            for i = 1:numel(A), A{i} = A{i}.multiply(f); end

            while true
                A = ratQ.bigVecTrim(A);
                if isempty(A) || numel(A) - 1 < db, break, end
                qr = A{1}.divideAndRemainder(lc);
                if qr(2).signum() ~= 0
                    error('ratQ:pseudoRemNotIntegral', ...
                        ['the pseudo-remainder step did not divide evenly, which the ' ...
                         'lc(b)^(deg a - deg b + 1) scaling guarantees it must.']);
                end
                c = qr(1);
                for k = 1:numel(Bv)
                    A{k} = A{k}.subtract(c.multiply(Bv{k}));
                end
                A{1} = java.math.BigInteger.ZERO;    % force the leading term to vanish exactly
            end

            A = ratQ.bigVecTrim(A);
            if isempty(A), r = []; return, end
            g = java.math.BigInteger.ZERO;
            for i = 1:numel(A), g = g.gcd(A{i}); end
            if g.signum() ~= 0 && ~g.equals(java.math.BigInteger.ONE)
                for i = 1:numel(A), A{i} = A{i}.divide(g); end
            end

            r = zeros(1, numel(A));
            for i = 1:numel(A)
                if A{i}.bitLength() > 53
                    error('ratQ:overflow', ...
                        ['the pseudo-remainder is still %d bits after its content is removed, so ' ...
                         'it is not an exact double and the Sturm chain cannot hold it. This is ' ...
                         'the answer being large, not the intermediate.'], A{i}.bitLength());
                end
                r(i) = double(A{i}.longValue());
            end
            r = ratQ.polyTrim(ratQ.chk(r, 'big pseudo-remainder'));
        end

        function C = bigVec(v)
        % objective: a univariate integer coefficient vector as a cell of BigInteger.
        % Kept separate from ratQ.toBig, which carries a size for MULTIVARIATE arrays and stores
        % zeros as EMPTY cells -- a representation the arithmetic here would have to special-case
        % on every access, for no gain at the handful of coefficients a Sturm chain entry has.
            ratQ.chk(v, 'coefficient vector');
            C = cell(1, numel(v));
            for i = 1:numel(v), C{i} = java.math.BigInteger.valueOf(int64(v(i))); end
        end

        function C = bigVecTrim(C)
        % objective: drop leading zero coefficients from a big coefficient vector.
            k = 1;
            while k <= numel(C) && C{k}.signum() == 0, k = k + 1; end
            C = C(k:end);
        end

        function r = polyPseudoRem(a, b)
        % objective: the PSEUDO-remainder of a by b -- the remainder of lc(b)^(deg a - deg b + 1) * a,
        %            which is integral where the ordinary remainder is not.
        %
        % The content is removed by the caller at every step; without that the leading coefficients
        % square each round and a six-digit quartic reaches 1e23 in a few steps.
            % PRIMITIVE PARTS FIRST. The scaling factor is lc(b)^(deg a - deg b + 1), so a large
            % leading coefficient is raised to a power: 10^6 x^2 against its own derivative
            % 2*10^6 x already needs 4e18, past 2^53, and ratQ.chk then raises. Removing the
            % content of both operands is exact, changes no root, and keeps the factor small --
            % here it turns 10^6 and 2*10^6 into 10^6 and 2, and the scaling into 4.
            a = ratQ.polyPrimitive(a);  b = ratQ.polyPrimitive(b);
            if isempty(b), error('ratQ:divideByZero', 'pseudo-remainder by the zero polynomial.'); end
            da = numel(a) - 1;  db = numel(b) - 1;
            if da < db, r = a; return, end
            % DOUBLES FIRST, BIG ONLY WHERE THEY RUN OUT -- the same shape resultantVar uses. The
            % REDUCED remainder is usually small even when the middle is not: measured on the
            % bridge region's membership quartic, the scaling reaches 3.7e18 and an elimination
            % step 1.7e16, while the answer fits comfortably. Content reduction cannot help WITHIN
            % a step, only between steps, so the middle has to be computed wide.
            %
            % THE FALLBACK COVERS THE WHOLE DOUBLE PATH, NOT JUST THE SCALING. A first version
            % guarded only the lc^(da-db+1) product, because that is where the first overflow was
            % seen -- and the very next fixture overflowed one line further down, inside the
            % elimination loop. Guarding the site rather than the computation is how a fallback
            % ends up needing a second fallback.
            %
            % ONLY ratQ:overflow IS CAUGHT. ratQ:pseudoRemNotIntegral means the scaling identity
            % failed, which is a defect in this routine, and retrying it wider would only hide it.
            try
                r = ratQ.polyPseudoRemDouble(a, b, da, db);
            catch e
                if ~strcmp(e.identifier, 'ratQ:overflow'), rethrow(e), end
                r = ratQ.polyPseudoRemBig(a, b);
            end
        end

        function r = polyPseudoRemDouble(a, b, da, db)
        % objective: the double-precision pseudo-remainder body, split out so its overflow can be
        %            caught as a whole rather than site by site. Assumes primitive, trimmed inputs.
            r = a * b(1)^(da - db + 1);
            ratQ.chk(r, 'pseudo-remainder scaling');
            % THE SCALING IS WHAT MAKES EVERY STEP EXACT, and the loop must rely on that rather
            % than test for it. Multiplying by lc(b)^(da-db+1) is exactly enough that each
            % elimination divides evenly, so a non-integral quotient here is a DEFECT, not a case
            % to break on -- and breaking on it left the degree unreduced, which spun forever.
            % Measured on 10^6 x^2 - (2*10^6+1) against its own derivative.
            while true
                r = ratQ.polyTrim(r);
                if isempty(r) || numel(r) - 1 < db, break, end
                c = r(1) / b(1);
                if c ~= round(c)
                    error('ratQ:pseudoRemNotIntegral', ...
                        ['the pseudo-remainder step did not divide evenly, which the ' ...
                         'lc(b)^(deg a - deg b + 1) scaling guarantees it must.']);
                end
                shift = numel(r) - numel(b);
                r = r - c * [b, zeros(1, shift)];
                r = ratQ.chk(r, 'pseudo-remainder step');
                r(1) = 0;                            % force the leading term to vanish exactly
            end
            % THE PRIMITIVE PART IS RETURNED, and that is a contract rather than a convenience.
            % The BigInteger path MUST reduce before it can come back to doubles -- reduction is
            % the only reason the answer fits when the middle does not -- so if the double path
            % returned the unreduced remainder the two paths would return different vectors for
            % the same input, and the differential test that is the only reason to trust the big
            % path could not exist. It costs nothing: sturmChain is the sole caller and takes
            % polyPrimitive of the result on the very next line. A positive content changes no
            % root and no sign, which is what makes it safe here.
            r = ratQ.polyPrimitive(ratQ.polyTrim(r));
        end

        function g = polyGcd(a, b)
        % objective: the primitive integer GCD of two integer polynomials.
        % Euclid on pseudo-remainders, taking the primitive part each round so the coefficients
        % stay bounded -- see polyPseudoRem for what happens without that.
            a = ratQ.polyPrimitive(a);  b = ratQ.polyPrimitive(b);
            if isempty(a), g = b; return, end
            if isempty(b), g = a; return, end
            while true
                % polyPremPrimitive, not polyPseudoRem: the two differ only by a content factor,
                % and polyGcd's answer is primitive either way -- but the reduced form never builds
                % the lc^(da-db+1) scaling that overflows on large coefficients.
                r = ratQ.polyTrim(ratQ.polyPremPrimitive(a, b));
                if isempty(r), g = ratQ.polyPrimitive(b); return, end
                a = b;  b = ratQ.polyPrimitive(r);
                if numel(b) <= 1, g = [1]; return, end   % a nonzero constant: coprime
            end
        end

        function s = polySquarefree(p)
        % objective: p with every repeated root reduced to a simple one, i.e. p / gcd(p, p').
        %
        % Sturm's theorem needs a SQUAREFREE polynomial, and a repeated root is the normal case
        % here rather than an exotic one: two circles meeting tangentially give the resultant a
        % double root. Falls back to p unchanged when the exact division does not verify -- the
        % callers' own counts still gate every answer, and a squarefree part that cannot be
        % confirmed is not worth risking.
            p = ratQ.polyTrim(p);
            if numel(p) <= 2, s = p; return, end
            try
                g = ratQ.polyGcd(p, polyder(p));
            catch ME
                if ~strcmp(ME.identifier, 'ratQ:overflow'), rethrow(ME); end
                s = p;  return
            end
            if numel(g) <= 1, s = p; return, end
            q = deconv(p, g);
            if any(abs(conv(q, g) - p) > 0), s = p; return, end   % exact division did not verify
            s = ratQ.polyTrim(round(q));
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

        function [s, strict] = quadSignOnPoly(c, P)
        % objective: does the quadratic c take ONE SIGN on the polyhedron P, and which. Exact.
        % [input]  c : 1 x 6 integer conic [a b c d e f] for a x^2 + b xy + c y^2 + d x + e y + f
        %          P : k x 3 integer rows [p1 p2 p0] meaning p1*x + p2*y + p0 >= 0
        % [output] s      : +1 if c >= 0 on P, -1 if c <= 0 on P, 0 if it takes BOTH signs or if
        %                   this routine cannot decide
        %          strict : true only when that sign is STRICT, i.e. c never VANISHES on P
        %
        % READ THE SECOND OUTPUT. A CALLER THAT IGNORES IT WILL BE WRONG, and this is not
        % hypothetical -- the bridge assembly used `s ~= 0` as its prune test and dropped a real
        % bridge (DECISIONS 2026-09-06, cell pair (5,5) of the diagonal-kink triangle pair). `+1`
        % means `c >= 0`, which ALLOWS c = 0 somewhere: the two pieces then touch without crossing,
        % which is a TANGENTIAL TIE and carries a bridge face just as a transversal one does. Only
        % `strict` says the quadratic never reaches zero, and only that licenses a prune.
        %
        % WHY THIS EXISTS -- FILTER 2 OF THE BRIDGE ASSEMBLY. A bridge face lives where two pieces'
        % conjugates TIE, so a cell pair whose difference g1 - g2 never changes sign carries no
        % bridge and its elimination (1 to 20 s) never has to run. DECISIONS 2026-09-06 measured
        % the stakes: the candidate count is quadratic in the cell count (235 pairs on a 16-gon
        % pair) while the number of pairs that genuinely tie is small and roughly linear (11), and
        % the free whole-plane test -- ratQ.conicSign -- reaches only an eighth of that on the same
        % fixture, its share FALLING as the polygon widens. The information is in the restriction
        % to the cell, which is what this routine adds.
        %
        % 0 MEANS "NOT PRUNABLE", NEVER "PRUNE ANYWAY". A wrong +1 or -1 discards a real bridge
        % face and produces a wrong envelope with no symptom; a wrong 0 costs one elimination. So
        % every case this cannot decide exactly -- a non-pointed recession cone, a direction along
        % which the quadratic is exactly flat, an integer overflow -- returns 0. The filter is
        % SOUND and deliberately not complete, and each incompleteness is named at its own guard.
        %
        % THE METHOD. s is read off the signs of the min and the max, and sign is monotone, so
        % sign(min) is the MIN OF THE SIGNS over a candidate set -- no rational is ever compared
        % against another, only against zero. That is what keeps the whole routine inside signQ.
        %
        %   min of a quadratic on a polyhedron is attained at a point whose FACE it lies in the
        %   relative interior of, so the candidates are exactly: the unconstrained stationary point
        %   (the 2-D face), the stationary point of the restriction to each facet LINE (the 1-D
        %   faces), and the vertices (the 0-D faces) -- each kept only if it lies in P.
        %
        % Every one of those is the solution of a 2 x 2 INTEGER system, so the whole enumeration is
        % ratQ.solve2 plus one membership test plus one sign. Extra candidates are harmless and
        % missing ones are not: an extra point can only pull the computed min DOWN, which loses a
        % prune, while a missing one could push it up and prune a pair that ties. So the stationary
        % point is kept without first asking whether it is a minimum or a saddle.
        %
        % NON-ATTAINMENT IS HANDLED FIRST, because a min of -infinity would otherwise be read off a
        % finite candidate and give exactly the unsound +1 this must never return. See recCone2.
            ratQ.chk(c, 'conic'); ratQ.chk(P, 'half-plane');
            if numel(c) ~= 6 || size(P,2) ~= 3
                error('ratQ:badQuadSignInput', ...
                    'quadSignOnPoly takes a 1 x 6 conic and k x 3 half-plane rows.');
            end
            mn = ratQ.quadMinSign(c, P);
            if ~isnan(mn) && mn >= 0, s = 1;  strict = (mn > 0); return, end   % q >= 0 on P
            mx = ratQ.quadMinSign(-c, P);
            if ~isnan(mx) && mx >= 0, s = -1; strict = (mx > 0); return, end   % q <= 0 on P

            % NEVER WEAKER THAN THE FREE TEST. conicSign asks the cruder question -- one sign over
            % the WHOLE plane -- and one sign on the plane implies one sign on any subset of it, so
            % it is a sound answer wherever it fires. It is kept as the last word because the
            % polyhedral route DECLINES on degenerate configurations that the plane test happens to
            % settle, and a filter that is sometimes weaker than the test it replaces is a
            % regression however much stronger it is on average. Measured on the corpus fixtures
            % this recovers pairs the polyhedral route alone loses (DECISIONS 2026-09-06).
            % c IS PASSED RAW, NOT THROUGH ratQ.conic. That routine returns the canonical
            % spelling of a CURVE, and a curve has no scale -- it makes the first nonzero
            % coefficient positive, which NEGATES half of all inputs. For a curve that is
            % identity; for a SIGN it is the wrong answer, and it read q <= 0 as q >= 0 on two of
            % four thousand random cases before the oracle caught it. conicSign needs only
            % integrality, which chk gives it.
            if any(c ~= 0)
                s = ratQ.conicSign(c);
                % conicSign is built on POSITIVE SEMIDEFINITENESS, so its +1 permits zeros -- a
                % parabola touching the axis is PSD and does vanish. The sign is strict exactly
                % when the 3x3 form is positive (or negative) DEFINITE, which is the leading-minor
                % test rather than the all-principal-minors one semidefiniteness needs.
                strict = (s ~= 0) && ratQ.conicDefinite(c);
                return
            end
            s = 0;  strict = false;
        end

        function tf = conicDefinite(c)
        % objective: is the conic's 3x3 form DEFINITE, i.e. does the quadratic never vanish?
        % [input]  c : 1 x 6 integer conic
        % [output] tf : logical, decided exactly
        %
        % Definite is the LEADING-minor test (Sylvester), unlike semidefinite which needs every
        % principal minor -- see ratQ.conicSign, which needs the other one for the other question.
        % Either sign counts: a form that is negative definite never vanishes either.
            ratQ.chk(c, 'conic');
            M = [2*c(1),   c(2),   c(4); ...
                   c(2), 2*c(3),   c(5); ...
                   c(4),   c(5), 2*c(6)];
            tf = leadingPositive(M) || leadingPositive(-M);
        end

        function s = quadMinSign(c, P)
        % objective: the SIGN of the minimum of the quadratic c over the polyhedron P.
        % [output] s : -1 if the min is negative or -infinity, 0 if it is exactly zero, +1 if it is
        %              positive, NaN if this routine declines to decide
        %
        % NaN IS THE CONSERVATIVE ANSWER and the caller must treat it as "unknown", never as 0.
        % The declines are named where they are taken; all of them are degenerate configurations,
        % and each costs at most one elimination that did not have to run.
            H  = [2*c(1), c(2); c(2), 2*c(3)];     % the Hessian, integral because the halves clear
            g0 = [c(4); c(5)];                     % grad q = H*x + g0

            % ---- 1. IS THE MINIMUM ATTAINED AT ALL? ---------------------------------------
            [rec, ok] = ratQ.recCone2(P);
            if ~ok, s = NaN; return, end           % declined: see recCone2
            switch size(rec, 1)
                case 0
                    % P is bounded: the min is attained, and every direction is irrelevant.
                case 1
                    a = rec(1,:) * H * rec(1,:).';
                    if a < 0, s = -1; return, end  % q -> -infinity along the one recession ray
                    if a == 0
                        r = ratQ.flatRaySign(H, g0, P, rec(1,:).');
                        if r < 0 || isnan(r), s = r; return, end
                    end
                case 2
                    u1 = rec(1,:).';  u2 = rec(2,:).';
                    t = ratQ.coneFormMinSign(H, u1, u2);
                    if t < 0, s = -1; return, end  % q -> -infinity inside the recession cone
                    if t == 0
                        % The form touches zero somewhere in the cone. Where matters: on an
                        % EXTREME RAY the flat direction is known explicitly and flatRaySign
                        % settles it, and every other direction then has a strictly positive form
                        % and cannot send q to -infinity. A zero attained strictly INSIDE the cone
                        % names no single direction here and is declined -- a named incompleteness.
                        f1 = (u1.'*H*u1 == 0);  f2 = (u2.'*H*u2 == 0);
                        if ~f1 && ~f2, s = NaN; return, end
                        if f1
                            r = ratQ.flatRaySign(H, g0, P, u1);
                            if r < 0 || isnan(r), s = r; return, end
                        end
                        if f2
                            r = ratQ.flatRaySign(H, g0, P, u2);
                            if r < 0 || isnan(r), s = r; return, end
                        end
                    end
                otherwise
                    s = NaN; return
            end

            % ---- 2. THE CANDIDATE POINTS, EACH A 2 x 2 INTEGER SOLVE ----------------------
            s = 1;                                  % min of the signs, started at the top
            k = size(P, 1);

            % (a) the 2-D face: grad q = 0.
            s = min(s, ratQ.quadSignAtSolve(c, P, H, -g0));

            % (b) the 1-D faces: on facet line n'x + p0 = 0 the restriction is stationary where the
            %     gradient is orthogonal to the line's direction u = [n2; -n1]. Two linear
            %     equations, so the same solve as everything else.
            for i = 1:k
                n = P(i,1:2).';
                if all(n == 0), continue, end
                u = [n(2); -n(1)];
                s = min(s, ratQ.quadSignAtSolve(c, P, [(u.'*H); n.'], [-(u.'*g0); -P(i,3)]));
            end

            % (c) the 0-D faces: every pair of boundary lines that meets.
            for i = 1:k-1
                for j = i+1:k
                    s = min(s, ratQ.quadSignAtSolve(c, P, [P(i,1:2); P(j,1:2)], ...
                                                          [-P(i,3); -P(j,3)]));
                end
            end
        end

        function s = quadSignAtSolve(c, P, A, b)
        % objective: solve A*x = b and return the sign of the quadratic c THERE, or +1 if the
        %            system is singular or the point is outside P.
        % [output] s : -1, 0 or +1, to be folded into a running MIN of signs
        %
        % +1 IS THE NEUTRAL ELEMENT of that fold, so "no candidate here" and "a candidate whose
        % value is positive" are correctly indistinguishable to the caller.
        %
        % A SINGULAR SYSTEM IS NOT AN ERROR HERE, unlike in ratQ.solve2 where it means a caller has
        % handed over a degenerate configuration. It is the ordinary case that the facet's
        % restriction has no stationary point (the quadratic is affine along that line) or that two
        % boundary lines are parallel; the minimum is then at some other candidate, which the fold
        % already visits. Likewise an OVERFLOW declines this candidate rather than the whole
        % question -- and it declines it towards +1, which can only lose a prune.
            s = 1;
            try
                if ratQ.detExact(A) == 0, return, end
                [xn, xd] = ratQ.solve2(A, b);
                if any(P(:,1:2) * xn + P(:,3) * xd < 0), return, end   % outside P; xd > 0
                s = sign(c(1)*xn(1)^2 + c(2)*xn(1)*xn(2) + c(3)*xn(2)^2 + ...
                         c(4)*xn(1)*xd  + c(5)*xn(2)*xd  + c(6)*xd^2);
                ratQ.chk([xn; xd], 'candidate point');
            catch e
                if strcmp(e.identifier, 'ratQ:overflow'), s = 1; return, end
                rethrow(e)
            end
        end

        function r = flatRaySign(H, g0, P, u)
        % objective: along a recession direction u on which the quadratic form is EXACTLY FLAT
        %            (u'Hu = 0), does the quadratic go to -infinity somewhere in P?
        % [output] r : -1 if it does, 0 if it does not, NaN if undecided
        %
        % WHY THIS CASE IS NOT RARE, and why it had to be built rather than declined. Over a VERTEX
        % cell a conjugate is AFFINE, so a pair involving one contributes a difference whose
        % Hessian is that of the other piece alone -- singular in the vertex cell's own directions.
        % Measured on the square-vs-square fixture, every pair the polyhedral route declined but the
        % whole-plane test settled was exactly this case (DECISIONS 2026-09-06).
        %
        % THE TEST IS LINEAR, which is the whole point. With u'Hu = 0,
        %       q(x + t*u) = q(x) + t * ((H*u)'x + g0'u),
        % exactly -- the t^2 term is gone -- so the slope along u is an AFFINE function of x, and q
        % is unbounded below exactly when that slope is negative somewhere in P. A constant slope
        % (H*u = 0) is read off directly; otherwise it is one feasibility test on P plus one row.
            a = H * u;
            b = g0.' * u;
            if all(a == 0)
                if b < 0, r = -1; else, r = 0; end
                return
            end
            % { x in P : slope(x) < 0 } is an open half-plane meeting a full-dimensional P, so it
            % is nonempty exactly when the closed system has a nonempty INTERIOR -- which is what
            % feasible2's strict mode decides, in integers.
            if ratQ.feasible2([P; -a.', -b], true), r = -1; else, r = 0; end
        end

        function [U, ok] = recCone2(P)
        % objective: the extreme rays of the recession cone { u : P(:,1:2)*u >= 0 }.
        % [output] U  : m x 2 integer directions, m = 0 (the cone is {0}), 1 (a ray) or 2 (a
        %               pointed 2-D cone); ok : false if the cone is NOT POINTED and the caller
        %               must decline
        %
        % WHY POINTEDNESS IS THE DIVIDING LINE. A pointed 2-D cone is the set of nonnegative
        % combinations of its two extreme rays, which is what lets coneFormMinSign decide the sign
        % of u'Hu over the whole cone from three integers. A cone containing a LINE -- P a slab, a
        % half-plane, or a single half-plane's worth of constraints -- is not spanned that way, and
        % rather than add a second parametrisation for a configuration a dual cell does not
        % produce, this declines. It is one of the named incompleteness of quadSignOnPoly.
        %
        % THE ENUMERATION. Every extreme ray of a pointed polyhedral cone in the plane lies on the
        % boundary of one of its own constraints, so the candidates are the two directions along
        % each constraint line, and an extreme ray is exactly a candidate that satisfies every
        % constraint. Directions are compared up to POSITIVE scaling, by cross product and sign.
            N = P(:,1:2);
            N = N(any(N ~= 0, 2), :);
            U = zeros(0,2);  ok = true;
            if isempty(N), ok = false; return, end          % nothing constrains: the cone is R^2
            C = zeros(0,2);
            for i = 1:size(N,1)
                C = [C; N(i,2), -N(i,1); -N(i,2), N(i,1)]; %#ok<AGROW>
            end
            keep = false(size(C,1),1);
            for i = 1:size(C,1)
                keep(i) = all(N * C(i,:).' >= 0);
            end
            C = C(keep,:);
            for i = 1:size(C,1)
                new = true;
                for j = 1:size(U,1)
                    if U(j,1)*C(i,2) - U(j,2)*C(i,1) == 0 && U(j,:) * C(i,:).' > 0
                        new = false; break                  % same direction, different length
                    end
                end
                if new, U = [U; C(i,:)]; end %#ok<AGROW>
            end
            if size(U,1) == 2 && U(1,1)*U(2,2) - U(1,2)*U(2,1) == 0
                ok = false; U = zeros(0,2); return          % opposite rays: the cone holds a line
            end
            if size(U,1) > 2, ok = false; U = zeros(0,2); end
        end

        function s = coneFormMinSign(H, u1, u2)
        % objective: the sign of the minimum of the quadratic FORM u'Hu over the pointed cone
        %            spanned by u1 and u2, excluding the origin.
        % [output] s : -1, 0 or +1
        %
        % The form is homogeneous of degree 2, so its sign on the cone is decided on the SEGMENT
        % from u1 to u2: with u = lambda*u1 + (1-lambda)*u2 it becomes the one-variable quadratic
        %       psi(lambda) = a2*lambda^2 + a1*lambda + a0,
        %       a2 = A - 2B + C,  a1 = 2*(B - C),  a0 = C,   A = u1'Hu1, B = u1'Hu2, C = u2'Hu2,
        % and psi's minimum over [0,1] is at an endpoint -- psi(0) = C, psi(1) = A -- or, when the
        % parabola opens upwards and its vertex lies strictly inside, at that vertex, where
        % 4*a2*psi = 4*a2*a0 - a1^2. Only signs are needed, so the vertex value never has to be
        % formed as a rational: with a2 > 0 the sign of psi there is the sign of that integer.
            A = u1.' * H * u1;   B = u1.' * H * u2;   C = u2.' * H * u2;
            a2 = A - 2*B + C;    a1 = 2*(B - C);
            s = min(sign(A), sign(C));
            if a2 > 0 && a1 < 0 && -a1 < 2*a2                 % vertex -a1/(2*a2) lies in (0,1)
                s = min(s, sign(ratQ.chk(4*a2*C - a1^2, 'cone form vertex')));
            end
        end
    end
end

function tf = leadingPositive(M)
% All three LEADING principal minors of a 3x3 integer symmetric matrix strictly positive, decided
% in exact integer arithmetic. Sylvester's criterion for POSITIVE DEFINITE.
    tf = M(1,1) > 0 ...
      && ratQ.detExact(M(1:2,1:2)) > 0 ...
      && ratQ.detExact(M) > 0;
end

function I = bisect(chain, lo, hi, n, depth)
% objective: split [lo,hi], which holds n roots, until every returned interval holds exactly one.
% The count comes from Sturm at every step, so the recursion is driven by an exact integer and
% cannot be fooled by two roots that lie close together.
    if n == 0, I = zeros(0,2); return, end
    if n == 1, I = [lo hi]; return, end
    if depth > 200
        error('ratQ:isolateNoConverge', ...
            ['bisection did not separate %d roots in 200 levels. For a SQUAREFREE integer ' ...
             'polynomial the roots are distinct and separated by a positive distance, so this ' ...
             'is a defect rather than a precision limit.'], n);
    end
    mid = (lo + hi) / 2;
    vL = ratQ.sturmV(chain, lo);
    vM = ratQ.sturmV(chain, mid);
    vH = ratQ.sturmV(chain, hi);
    I = [bisect(chain, lo, mid, vL - vM, depth+1); bisect(chain, mid, hi, vM - vH, depth+1)];
end
