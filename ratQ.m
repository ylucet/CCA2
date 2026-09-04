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
