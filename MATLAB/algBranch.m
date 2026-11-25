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
%                 sorted real roots), for callers that need to see the whole fibre -- the subdivision
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
% value; no subdivision predicate calls this routine, precisely so that no decision rests on it.
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
%             the bridge's two branches meet (2026-09-25). It decides nothing a subdivision depends on;
%             it only stops the cache from dropping a value that exists. Inside that window but
%             past dq's rounding (16 eps of the terms), a positive discriminant still gets its two
%             roots (stable formula) -- the mean was 6.6e-6 off a branch near the poly2 kink's merge
%             curve (2026-10-30).

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

    % VECTORISED for a double fibre of degree <= 2 in z (2026-11-01): erratum Example 2's overlay
    % called `roots` 7.1 million times, one point at a time. The same decisions, made on columns:
    % the leader dropped where it is zero to its own rounding, the double-root windows below, and
    % a complex pair counted real under imagTol. The values are the stable quadratic formula's,
    % within a few ulps of `roots`; algBranchTest checks the two paths against each other.
    if nargout < 2 && nz <= 3
        z = fibreLow(P, Pbig, X, b, opts, nx, ny, nz);
        return
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
            % caller (a subdivision builder checking branch counts) wants to see it, not to be stopped.
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
        % A NEAR-DOUBLE REAL PAIR that is real only on a CURVE is not a value (#108 residual,
        % D-153): see dropIsolatedPairs. Python and Rust also examine the complex pairs their
        % double-root noise test accepts; MATLAB drops those as before (PORT_NOTES section 4).
        rr = dropIsolatedPairs(P, Pbig, X(t,:), cz, r, rr, nx, ny, nz);
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
            T = ca(2)^2 + 4*ca(3)*ca(1);
            % ...but the window is wide (1e-12 of the terms, ~5000x dq's rounding error), and a
            % computed dq > 0 past ROUNDING is two real roots up to ~1e-5 apart: collapsing them to
            % the mean put a branch 6.6e-6 off near the poly2 kink's merge curve (2026-10-30,
            % aNEARDoubleRootIsNOTCollapsedToTheMean), so those get the stable quadratic formula.
            % Within rounding (16 eps of the terms) it IS the double root and the mean is exact --
            % splitting a true D = 0 by sqrt(rounding) put erratum Example 1's edge y = 2x 2e-8 high.
            % EXACT inside the window (2026-11-03): see WITHIN THE WINDOW, EXACTLY in fibreLow.
            if abs(dq) <= 1e-12 * T
                [cx, dx] = exactQuad(P, Pbig, X(t,:), nx, ny);
                if dx > 0
                    qq = -(cx(2) + sign(cx(2) + (cx(2) == 0)) * sqrt(dx)) / 2;
                    rr = sort([qq / cx(3); cx(1) / qq]);
                else
                    rr = [1; 1] * (-cz(2) / (2*cz(3)));
                end
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

function z = fibreLow(P, Pbig, X, b, opts, nx, ny, nz)
% algBranch's loop for fibres of degree <= 2 in z, on all points at once. A big P's fibre is still
% summed exactly point by point (exactFibre); only the roots are taken together, so the big and the
% double form of one polynomial give the same doubles (bridgeFaceAlgTest).
    n = size(X,1);
    Vx = X(:,1) .^ (0:nx-1);  Vy = X(:,2) .^ (0:ny-1);
    Ax = abs(Vx);  Ay = abs(Vy);
    c = zeros(n, 3);  ca = zeros(n, 3);
    for k = 1:nz
        ca(:,k) = sum((Ax * abs(P(:,:,k))) .* Ay, 2);
    end
    if isempty(Pbig)
        for k = 1:nz
            c(:,k) = sum((Vx * P(:,:,k)) .* Vy, 2);
        end
        isZero = abs(c) <= opts.degTol * ca;
    else
        isZero = false(n, 3);
        for t = 1:n
            [c(t,1:nz), isZero(t,1:nz)] = exactFibre(Pbig, X(t,:), nx, ny, nz);
        end
    end
    isZero(:, nz+1:3) = true;
    z = nan(n, 1);
    lead = 3 * ~isZero(:,3);
    lead(lead == 0 & ~isZero(:,2)) = 2;              % the fibre's degree after dropping zero leaders
    % ---- LINEAR: one root, -c0/c1 (what `roots` computes for a 1 x 1 companion) ----------------
    L = lead == 2 & max(abs(c), [], 2) > 0;
    if any(L) && (b == 1 || b == -1)
        z(L) = -c(L,1) ./ c(L,2);
    end
    % ---- QUADRATIC ---------------------------------------------------------------------------
    Q = lead == 3;
    if ~any(Q), return, end
    a = c(Q,3);  bb = c(Q,2);  c0 = c(Q,1);
    dq = bb.^2 - 4*a.*c0;
    T = ca(Q,2).^2 + 4*ca(Q,3).*ca(Q,1);
    lo = nan(size(a));  hi = lo;
    mean2 = -bb ./ (2*a);
    % WITHIN THE WINDOW, EXACTLY (2026-11-03). Where |dq| is within 1e-12 of its terms the doubles
    % cannot say whether the two roots are equal: on the 3-convex-edge triangle a face's fibre has
    % D ~ (k d)^2 at distance d from the hypotenuse -- two branches CROSSING there -- and up to
    % d ~ 3e-7 the rounded dq read as 0, so the evaluator returned their mean, ~1e-7 off either
    % root, and a point 1e-8 inside the domain fell out of every face. The point's coordinates and
    % P's coefficients are doubles, so the fibre's coefficients and D are EXACT in BigDecimal: D > 0
    % splits (stable formula on the exactly rounded coefficients), D <= 0 keeps the mean as before.
    win = abs(dq) <= 1e-12 * T;
    if any(win)
        iq = find(Q);  w = find(win).';
        for u = w
            [cx, dx] = exactQuad(P, Pbig, X(iq(u),:), nx, ny);
            if dx > 0, a(u) = cx(3);  bb(u) = cx(2);  c0(u) = cx(1);  dq(u) = dx;
            else,      dq(u) = 0;
            end
        end
        mean2 = -bb ./ (2*a);
    end
    dbl = win & dq <= 0;
    lo(dbl) = mean2(dbl);  hi(dbl) = mean2(dbl);
    two = ~dbl & dq > 0;
    qq = -(bb(two) + sign(bb(two) + (bb(two) == 0)) .* sqrt(dq(two))) / 2;
    r1 = qq ./ a(two);  r2 = c0(two) ./ qq;
    lo(two) = min(r1, r2);  hi(two) = max(r1, r2);
    cpx = ~dbl & dq < 0;                             % a complex pair: real only within imagTol
    im = sqrt(-dq(cpx)) ./ (2 * abs(a(cpx)));
    re = mean2(cpx);
    ok = im <= opts.imagTol * max(1, hypot(re, im));
    w = find(cpx);  w = w(ok);
    lo(w) = re(ok);  hi(w) = re(ok);
    v = nan(size(a));
    if b == 1 || b == -2, v = lo; elseif b == 2 || b == -1, v = hi; end
    zq = nan(nnz(Q), 1);  zq(:) = v;
    z(Q) = zq;
end

function [c, D] = exactQuad(P, Pbig, x, nx, ny)
% A quadratic fibre's coefficients c (z^0, z^1, z^2) and discriminant D, computed EXACTLY in
% BigDecimal from the double point x and P (double or big), then each rounded to a double: D's
% sign is exact.
    X = java.math.BigDecimal(x(1));  Y = java.math.BigDecimal(x(2));
    xp = cell(1, nx);  yp = cell(1, ny);
    xp{1} = java.math.BigDecimal.ONE;  yp{1} = java.math.BigDecimal.ONE;
    for i = 2:nx, xp{i} = xp{i-1}.multiply(X); end
    for j = 2:ny, yp{j} = yp{j-1}.multiply(Y); end
    acc = repmat({java.math.BigDecimal.ZERO}, 1, 3);
    for k = 1:3
        for i = 1:nx
            for j = 1:ny
                if isempty(Pbig)
                    v = P(i,j,k);
                    if v == 0, continue, end
                    v = java.math.BigDecimal(v);
                else
                    v = Pbig.c{sub2ind(Pbig.sz, i, j, k)};
                    if isempty(v), continue, end
                    v = java.math.BigDecimal(v);
                end
                acc{k} = acc{k}.add(v.multiply(xp{i}).multiply(yp{j}));
            end
        end
    end
    Dx = acc{2}.multiply(acc{2}).subtract(java.math.BigDecimal(4).multiply(acc{3}).multiply(acc{1}));
    c = cellfun(@(a) a.doubleValue(), acc);
    D = Dx.doubleValue();
    if D == 0 && Dx.signum() ~= 0, D = Dx.signum() * realmin; end
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

function rr = dropIsolatedPairs(P, Pbig, x, cz, r, rr, nx, ny, nz)
% Drop the NEAR-DOUBLE real roots that are not the value of any branch (#108 residual, D-153).
%
% A face is the b-th real root of its fibre, and f** is closed and convex: its value at a point is
% the limit from the points around it. So a double root that is real only on a CURVE -- a complex
% pair touching the real axis there and leaving it on both sides -- must not be counted. Such a
% pair appears where a bridge's dual point runs off to infinity, on a boundary line of dom f**. For
% a continuous f it coincides with the true value there; for a discontinuous one it does not
% (D-153: the convex piece's own value 10.2578 on its hull edge next to the jump vertex, f** =
% 9.5177), and rounding returns it as two close reals in a ~1e-5 strip inside the domain too.
%
% Only a pair the doubles cannot separate is examined: two adjacent real roots within the double-
% root noise of each other, 8 sqrt(eps sum|c_k t^k| / |p''(t)/2|) (alg_branch.py's
% _doubleRootNoise). The fibre is then evaluated EXACTLY (doubles are dyadic) at x +- h e for four
% directions e and a tiny dyadic h (~1e-12 relative), and its distinct real roots near the pair are
% counted by Sturm in BigInteger. A pair with no real root on EITHER side in at least two
% directions is dropped. A crossing of two real branches, a fold (real on one side) and a double
% root that is real along a curve all keep it: at most one of the four directions runs along it.
    cands = [];
    for i = 1:numel(rr)-1
        a = rr(i);  b = rr(i+1);
        if b > a && doubleRootNoise(cz, (a+b)/2, (b-a)/2), cands(end+1) = (a+b)/2; end %#ok<AGROW>
    end
    % a near-real complex pair (dropped above) is not examined, but it may be a multiple root's
    % rounding: then the multiple roots are decided exactly below, as in Python and Rust
    snap = ~isempty(cands) || any(arrayfun(@(v) imag(v) > 0 && ...
        doubleRootNoise(cz, real(v), imag(v)), r(:)));
    for z0 = cands
        d = sort(abs(r(:) - z0));  d = d(3:end);       % the pair itself is nearest
        if ~isempty(d) && d(1) > 0, rad = d(1)/2; else, rad = 1 + abs(z0); end
        if isolatedPair(P, Pbig, x, z0, rad, nx, ny, nz)
            for rep = 1:2
                [m, k] = min(abs(rr - z0));
                if ~isempty(k) && m <= rad, rr(k) = []; end
            end
        end
    end
    if snap
        rr = snapMultipleRoots(P, Pbig, x, cz, rr, nx, ny, nz);
    end
end

function rr = snapMultipleRoots(P, Pbig, x, cz, rr, nx, ny, nz)
% A GENUINE multiple root of the fibre, to full precision (#108 residual: the jump vertices).
% Rounding returns a k-fold real root as k roots ~eps^(1/k) apart, and at a jump vertex f** is
% such a root: D-153's vertex (-3,-2), f** = 7.5 a double root, came back as 7.5000017 (2.3e-7
% relative), the poly2 kink's (1,2) as 1.0000000256. The multiple roots are the roots of
% gcd(p, p'), where they are SIMPLE: the fibre at the exact dyadic point is split by multiplicity
% (Yun), each real root of a factor of multiplicity k >= 2 is isolated by Sturm and bisected
% exactly to adjacent doubles, and whatever rounding left within 8 (eps sum|c_j m^j| /
% |p^(k)(m)/k!|)^(1/k) of it (close reals, or a lone real whose partners came back complex)
% becomes exactly k copies of it, whatever the LAPACK in use (alg_branch.py _snapMultipleRoots).
    if isempty(rr), return, end
    [xn, xd] = ratQ.dyadicParts(x(1));
    [yn, yd] = ratQ.dyadicParts(x(2));
    D = xd.max(yd);
    row = ratQ.bigPrimitive(exactFibreRow(P, Pbig, xn.multiply(D.divide(xd)), ...
                                          yn.multiply(D.divide(yd)), D, nx, ny, nz));
    if numel(row) <= 2, return, end
    fac = yunBig(row);
    for f = 1:numel(fac)
        k = fac{f}{1};  a = fac{f}{2};
        if k < 2, continue, end
        I = ratQ.isolateRootsBig(a);
        for q = 1:size(I,1)
            m = bisectBig(a, I(q,1), I(q,2));
            w = multipleRootNoise(cz, m, k);
            near = abs(rr - m) <= w;
            if any(near)                    % rounding left something there: k copies of m
                rr = [rr(~near); repmat(m, k, 1)];
            end
        end
    end
    rr = sort(rr);
end

function w = multipleRootNoise(cz, m, k)
% 8 (eps sum|c_j m^j| / |p^(k)(m)/k!|)^(1/k): how far rounding moves a k-fold root m (cz z^0 first)
    j = 0:numel(cz)-1;
    cabs = sum(abs(cz) .* abs(m).^j);
    jj = j(j >= k);
    pk = 0;
    for t = jj, pk = pk + nchoosek(t, k) * cz(t+1) * m^(t-k); end
    if pk == 0 || cabs == 0, w = 0; return, end
    w = 8 * (eps * cabs / abs(pk))^(1/k);
end

function fac = yunBig(p)
% Yun's square-free decomposition of a BigInteger row (highest first): {k, a_k} with
% p = c prod a_k^k, each a_k square-free, primitive, of positive degree.
    g = ratQ.bigGcdRow(p, ratQ.bigDeriv(p));
    if numel(g) <= 1, fac = {{1, p}}; return, end
    fac = {};
    b = bigDivRow(ratQ.bigPrimitive(p), g);          % every distinct root once
    k = 1;  rest = g;
    while numel(b) > 1
        a = ratQ.bigGcdRow(b, rest);                 % the roots of multiplicity > k
        if numel(a) > 1, f = bigDivRow(b, a); else, f = b; end
        if numel(f) > 1, fac{end+1} = {k, ratQ.bigPrimitive(f)}; end %#ok<AGROW>
        if numel(a) <= 1, break, end
        rest = bigDivRow(rest, a);
        b = a;  k = k + 1;
    end
end

function Q = bigDivRow(P, G)
% P / G for BigInteger rows (highest first) when the division is exact over the integers.
    n = numel(P) - numel(G) + 1;
    Q = cell(1, max(n, 0));
    R = P;
    for i = 1:n
        qr = R{i}.divideAndRemainder(G{1});
        if qr(2).signum() ~= 0
            error('algBranch:internal', 'an exact polynomial division left a remainder.');
        end
        Q{i} = qr(1);
        for j = 1:numel(G), R{i+j-1} = R{i+j-1}.subtract(Q{i}.multiply(G{j})); end
    end
    for i = n+1:numel(R)
        if R{i}.signum() ~= 0
            error('algBranch:internal', 'an exact polynomial division left a remainder.');
        end
    end
end

function m = bisectBig(a, lo, hi)
% The simple root of the square-free a in (lo, hi], bisected exactly down to adjacent doubles;
% the upper end, which holds the root. A root at z = 0 is returned as 0: bisecting toward it ran
% into subnormal endpoints, which ratQ.dyadicParts cannot write (D-069: a Java exception).
    if a{end}.signum() == 0 && lo < 0 && hi >= 0, m = 0; return, end
    for it = 1:2000
        mid = (lo + hi) / 2;
        if mid <= lo || mid >= hi, break, end
        if ratQ.countRootsInHalfOpenBig(a, lo, mid) == 1, hi = mid; else, lo = mid; end
    end
    m = hi;
end

function tf = doubleRootNoise(cz, t, im)
% Is a pair at t +- i*im within the rounding of a DOUBLE root of the fibre cz (z^0 first)?
    k = 0:numel(cz)-1;
    cabs = sum(abs(cz) .* abs(t).^k);
    kk = k(k >= 2);
    p2 = sum(kk.*(kk-1)/2 .* cz(kk+1) .* t.^(kk-2));
    tf = p2 ~= 0 && cabs ~= 0 && abs(im) <= 8 * sqrt(eps * cabs / abs(p2));
end

function tf = isolatedPair(P, Pbig, x, z0, rad, nx, ny, nz)
% No real root within rad of z0 on either side of x, in at least two of four directions.
    [~, e] = log2(max([1, abs(x(1)), abs(x(2))]));      % frexp's exponent
    [xn, xd] = ratQ.dyadicParts(x(1));
    [yn, yd] = ratQ.dyadicParts(x(2));
    two = java.math.BigInteger.valueOf(2);
    if e >= 40, hn = two.pow(e - 40);  hd = java.math.BigInteger.ONE;
    else,       hn = java.math.BigInteger.ONE;  hd = two.pow(40 - e);
    end
    D = xd.max(yd).max(hd);                             % all powers of two: the lcm
    X0 = xn.multiply(D.divide(xd));  Y0 = yn.multiply(D.divide(yd));  H = hn.multiply(D.divide(hd));
    probes = [1 0; 0 1; 1 1; 1 -1];
    both = 0;
    for q = 1:4
        none = true;
        for sg = [1 -1]
            XN = X0.add(H.multiply(java.math.BigInteger.valueOf(sg*probes(q,1))));
            YN = Y0.add(H.multiply(java.math.BigInteger.valueOf(sg*probes(q,2))));
            row = exactFibreRow(P, Pbig, XN, YN, D, nx, ny, nz);
            if ratQ.countRootsInHalfOpenBig(row, z0 - rad, z0 + rad) ~= 0
                none = false;  break
            end
        end
        if none, both = both + 1; end
    end
    tf = both >= 2;
end

function row = exactFibreRow(P, Pbig, XN, YN, D, nx, ny, nz)
% The fibre at (XN/D, YN/D) times D^(nx-1+ny-1), as a BigInteger row, highest power of z first.
    xp = cell(1, nx);  yp = cell(1, ny);
    for i = 1:nx, xp{i} = XN.pow(i-1).multiply(D.pow(nx-i)); end
    for j = 1:ny, yp{j} = YN.pow(j-1).multiply(D.pow(ny-j)); end
    acc = repmat({java.math.BigInteger.ZERO}, 1, nz);
    for k = 1:nz
        for i = 1:nx
            for j = 1:ny
                if isempty(Pbig)
                    v = P(i,j,k);
                    if v == 0, continue, end
                    v = java.math.BigInteger(sprintf('%.0f', v));
                else
                    v = Pbig.c{sub2ind(Pbig.sz, i, j, k)};
                    if isempty(v), continue, end
                    v = java.math.BigInteger(v.toString());
                end
                acc{k} = acc{k}.add(v.multiply(xp{i}).multiply(yp{j}));
            end
        end
    end
    row = fliplr(acc);
end

