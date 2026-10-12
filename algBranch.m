function [z, info] = algBranch(P, X, b, opts)
% algBranch  Evaluate ONE BRANCH of an implicit algebraic face p(x,y,z) = 0.
%
% objective: an `Alg` face carries no formula -- it is the set {p(x,y,z) = 0} for an exact integer
%   polynomial p, plus a BRANCH INDEX saying which of p(x,y,.)'s real roots is the function's
%   value. This is the evaluator for that pair, and it is the whole of what `eval` does on such a
%   face. See Alg.m for why the representation is implicit and AlgCon.m for where P and b live.
%
% [input]  P    : the polynomial, as ratQ's dense polyN array in THREE variables: P(i+1,j+1,k+1)
%                 is the integer coefficient of x^i y^j z^k. Sizes state the degrees (ratQ.trimN).
%          X    : n x 2 evaluation points [x y].
%          b    : the branch index in the CANONICAL order defined below -- 1-based from the
%                 BOTTOM when positive, and from the TOP when negative (-1 is the largest real
%                 root). Never 0: that value means "this face carries a formula, not a branch".
%          opts : (optional) struct, fields `imagTol` and `degTol` -- see TOLERANCES.
% [output] z    : n x 1 values, NaN where the branch does not exist at that point.
%          info : struct with `nreal` (n x 1 real-root count) and `allReal` (n x 1 cell of the
%                 sorted real roots), for callers that need to see the whole fibre -- the mesh
%                 builder does, to check that the branch count is constant on a face.
%                 COUNTED WITH MULTIPLICITY: a double root is two entries, so branches b and b+1
%                 agree there. Collapsing them would need a tolerance (`roots` separates a true
%                 double root by about sqrt(eps), never exactly) and that tolerance would be
%                 deciding a STRUCTURAL question -- do two branches merge here -- which is the one
%                 kind of question this project never lets a tolerance answer. The merge locus is
%                 the discriminant curve of p in z, which a well-formed face has on its BOUNDARY.
%
% ------------------------------------------------------------------------------------------------
% THE CANONICAL ORDER, and why an index is a legitimate way to name a value.
%
% The real roots of p(x,y,.) are sorted ASCENDING and `b` selects among them, 1-based. This is
% exactly the convention QuaCon.Vname already uses one dimension down, where a vertex is named as
% [edgeA edgeB rootIdx] -- the rootIdx-th intersection of two edge conics -- rather than by its
% coordinates, because the coordinates are irrational and the name is not. The same argument
% applies here: the VALUE of a bridge face is irrational, the polynomial and the index are not.
%
% WHY THE INDEX IS WELL-DEFINED ON A FACE, which is the property that makes it a name rather than
% a guess. Two distinct real roots can only swap places by colliding, and at a collision p and
% dp/dz share a root, i.e. the discriminant of p with respect to z vanishes. That is an algebraic
% curve in (x,y), so within one open face the real roots vary continuously without crossing and
% "the b-th real root" names one continuous branch. ACROSS faces it may name a different branch,
% which is why b is stored PER FACE.
%
% COUNTING FROM WHICH END IS NOT COSMETIC, and the answer was measured rather than chosen. Two
% roots can also collide and LEAVE the real line entirely -- the count drops by two without the
% surviving branches moving at all -- and then every index ABOVE the departed pair shifts down
% while every index BELOW it is untouched. So an end matters, and the stable end is the one the
% face's own value sits at.
%
% For the convex envelope's bridge face the value is the LARGEST real root, and that is not a
% coincidence to be rechecked per fixture: the face is a SUPREMUM over the tie curve, the real
% roots of p are exactly that problem's critical values, and a supremum is the largest of them.
% Measured on three pairs over 4758 region points: the value is the largest real root at every one
% (worst 7.4e-14), while its ASCENDING index takes both 2 and 4 inside a single face. Ascending
% indexing therefore cannot name a bridge face, and a negative b -- counting from the top, -1 for
% the largest -- can. DECISIONS.md 2026-09-17.
%
% Positive b keeps its original meaning bit for bit, so nothing stored under the ascending reading
% changes; the negative spelling is an addition, for the producer that needed it.
%
% WHAT THIS ROUTINE IS NOT. It is a DOUBLE-PRECISION evaluator -- the cache, in the sense
% QuaCon.Vx is a cache of the vertex it names. The exact object is (P, b); this returns the double
% that names the same value, and deleting every number it returns changes runtime only. An exact
% route exists and is deliberately not taken here: isolate the real roots by a Sturm sequence over
% Q (ratQ already counts Sturm sequences for the bridge region's own membership test) and refine
% to a rational interval. That is the follow-up for a caller that needs a DECISION rather than a
% value; no mesh predicate calls this routine, precisely so that no decision rests on it.
%
% TOLERANCES, both of which are of the cache and not of the mathematics.
%   imagTol : a computed root counts as real when |imag| <= imagTol * max(1,|root|). `roots` returns
%             a genuinely real root with an imaginary part of order eps*|root| -- scaling by the
%             root's own magnitude is what makes the test mean the same thing at |z| = 1e-3 and at
%             1e6, and an absolute tolerance measured against the wrong scale is the single most
%             repeated measurement defect in this project's own record (DECISIONS.md 2026-09-07).
%   degTol  : a leading coefficient within degTol * (the sum of the absolute values of its own
%             terms) is a zero computed inexactly and is dropped before root-finding, because
%             `roots` on a vector with a ~0 leading entry returns a spurious root near infinity.
%             Relative to its OWN terms, not to the largest coefficient: that version dropped a
%             genuine leader (2026-10-11). For a big P the zero test is exact.
%   the DOUBLE-ROOT test (quadratic fibres only, fixed at 1e-12): a discriminant within 1e-12 of
%             zero, relative to the size of its TERMS before cancellation, is a double root -- two
%             equal real entries, the multiplicity contract above. Without it `roots` can return
%             the double root as a complex pair past imagTol, and the point has no branch at all:
%             the c = -2 AlgCrv evaluated to +Inf on parts of its domain edges x = 0 and y = 1, where
%             the bridge's two branches meet (2026-09-25). It decides nothing a mesh depends on;
%             it only stops the cache from dropping a value that exists.

    if nargin < 4, opts = struct(); end
    % A polynomial kept in BigInteger (past 2^53 -- bridgeFaceAlg) is NOT simply rounded: its
    % fibre coefficients cancel heavily, and rounding P first was measured wrong by 4 in 1167 on
    % erratum Example 2 (2026-10-11). So for such a P the fibre p(x,y,.) is summed EXACTLY, in
    % BigDecimal at the double point, and only its coefficients are rounded -- after the
    % cancellation, not before. The rounded P is kept for the shape checks and the double-root test.
    Pbig = [];
    if ratQ.isBigPoly(P), Pbig = P; end
    P = ratQ.bigToDouble(P);
    if ~isfield(opts,'imagTol'), opts.imagTol = 1e-10; end
    if ~isfield(opts,'degTol'),  opts.degTol  = 1e-12; end

    if ~isnumeric(P) || isempty(P)
        error('algBranch:polynomial', 'P must be a non-empty numeric array; got %s.', class(P));
    end
    if ndims(P) > 3
        error('algBranch:polynomial', ...
            ['P must be a polynomial in at most THREE variables (x,y,z); got %d dimensions. ' ...
             'A face of a bivariate function has no fourth variable to carry.'], ndims(P));
    end
    if ~isnumeric(X) || size(X,2) ~= 2
        error('algBranch:shape', 'X must be n x 2 [x y]; got %dx%d.', size(X,1), size(X,2));
    end
    if ~isscalar(b) || ~isnumeric(b) || b ~= fix(b) || b == 0
        error('algBranch:branch', ...
            ['the branch index must be a nonzero integer -- positive counts up from the ' ...
             'smallest real root, negative down from the largest; got %s.'], mat2str(b));
    end

    sz = size(P);  sz(end+1:3) = 1;                 % pad: a face of degree 0 in z is still 3-D here
    nx = sz(1);  ny = sz(2);  nz = sz(3);
    if nz < 2
        error('algBranch:noZ', ...
            ['P has degree 0 in z, so {p = 0} does not define a value: it is a condition on ' ...
             '(x,y) alone. A face stored this way is a producer defect, not an evaluation one.']);
    end

    n = size(X,1);
    z = nan(n,1);
    allReal = cell(n,1);
    nreal = zeros(n,1);

    for t = 1:n
        % ---- substitute (x,y), leaving a univariate polynomial in z --------------------------
        % c(k) = sum_{i,j} P(i+1,j+1,k) x^i y^j, i.e. one bilinear form per power of z. Done as a
        % vector-matrix-vector product rather than a double loop so the cost is one BLAS call per
        % power of z; for the measured degrees (at most (6,6,4)) either is free, and this one says
        % what it means.
        xs = X(t,1) .^ (0:nx-1);                     % 1 x nx
        ys = (X(t,2) .^ (0:ny-1)).';                 % ny x 1
        c = zeros(1, nz);  cabs = zeros(1, nz);
        if isempty(Pbig)
            for k = 1:nz
                c(k) = xs * P(:,:,k) * ys;
                cabs(k) = abs(xs) * abs(P(:,:,k)) * abs(ys);   % the sum's own size, for its error
            end
            isZero = abs(c) <= opts.degTol * cabs;
        else
            [c, isZero] = exactFibre(Pbig, X(t,:), nx, ny, nz);
        end

        % ---- real roots, ascending -------------------------------------------------------------
        scale = max(abs(c));
        if scale == 0
            % p(x,y,.) is identically zero: every z satisfies it, so the face names no value here.
            % Left as NaN rather than raised -- it is a legitimate question with no answer, and the
            % caller (a mesh builder checking branch counts) wants to see it, not to be stopped.
            continue
        end
        cz = c(:).';
        % DROP A LEADER ONLY IF IT IS ZERO -- exactly, for a big P; within its OWN rounding error
        % (degTol times the sum of the absolute values of its terms), for a double one. Not
        % "small next to the largest coefficient": that dropped a genuine leading 16 against
        % 1.6e17 on erratum Example 2 and returned a cubic's root (2026-10-11).
        keep = find(~isZero, 1, 'last');
        cz = cz(1:keep);
        if numel(cz) < 2
            continue                                  % degenerated to a constant: no root in z
        end
        r = roots(fliplr(cz));                        % `roots` wants highest power first

        isReal = abs(imag(r)) <= opts.imagTol * max(1, abs(r));
        rr = sort(real(r(isReal)));
        % A DOUBLE ROOT OF A QUADRATIC FIBRE -- where the two branches meet -- can come back from
        % `roots` as a complex pair with imaginary part ~ sqrt(eps), past imagTol, leaving the
        % point with no branch at all: the c = -2 AlgCrv evaluated to +Inf on its domain edge x = 0
        % that way (algBranchTest.aDOUBLERootThatRootsReturnsAsAComplexPairIsStillReal). The
        % discriminant is the right scale to ask "is this a double root": where it vanishes to
        % relative precision the answer is -b/(2a) twice, which is the multiplicity contract of
        % theBRANCHCOUNTIsReportedWITHMULTIPLICITY. Every other fibre is decided by `roots` as before.
        % The scale is the size of the TERMS before cancellation, not of the coefficients: along
        % the edge y = 1 of the same object the discriminant is 0 everywhere, but near the corner
        % x = sqrt2 - 1 the coefficients themselves are ~1e-7 (they cancel), and a tolerance
        % relative to them rejected the rounding of an exact zero.
        if numel(cz) == 3
            ca = zeros(1, 3);
            for k = 1:3, ca(k) = abs(xs) * abs(P(:,:,k)) * abs(ys); end
            dq = cz(2)^2 - 4*cz(3)*cz(1);
            if abs(dq) <= 1e-12 * (ca(2)^2 + 4*ca(3)*ca(1))
                rr = [1; 1] * (-cz(2) / (2*cz(3)));
            end
        end
        allReal{t} = rr;
        nreal(t) = numel(rr);
        idx = b;
        if idx < 0, idx = numel(rr) + 1 + b; end        % -1 is the largest, -2 the next below
        if idx >= 1 && idx <= numel(rr)
            z(t) = rr(idx);
        end
    end

    if nargout > 1
        info = struct('nreal', nreal, 'allReal', {allReal});
    end
end

function [c, isZero] = exactFibre(B, x, nx, ny, nz)
% sum_ij P(i,j,k) x^i y^j for each k, exactly (BigDecimal holds every double exactly), then rounded
    X = java.math.BigDecimal(x(1));  Y = java.math.BigDecimal(x(2));
    xp = cell(1, nx);  yp = cell(1, ny);
    xp{1} = java.math.BigDecimal.ONE;  yp{1} = java.math.BigDecimal.ONE;
    for i = 2:nx, xp{i} = xp{i-1}.multiply(X); end
    for j = 2:ny, yp{j} = yp{j-1}.multiply(Y); end
    acc = repmat({java.math.BigDecimal.ZERO}, 1, nz);
    sz = [nx ny nz];
    for t = 1:numel(B.c)
        v = B.c{t};
        if isempty(v), continue, end
        [i, j, k] = ind2sub(sz, t);
        acc{k} = acc{k}.add(java.math.BigDecimal(v).multiply(xp{i}).multiply(yp{j}));
    end
    c = cellfun(@(a) a.doubleValue(), acc);
    isZero = cellfun(@(a) a.signum() == 0, acc);
end

