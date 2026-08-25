function [P, info] = conicMeet(A, B, opts)
% conicMeet  The real intersection points of two conics, from an EXACT integer quartic.
%
% objective: build the vertex layer of a QuaCon without any symbolic call and without ever
%   converting a computed coordinate to a rational. Two rational conics meet in at most four
%   points; those points have degree <= 4 over Q and are generically irrational (of twelve
%   feasible continuous three-piece configurations the vertex quartic is irreducible over Q in
%   TEN -- CONJ_FIELD_PROOF.md 8.0), so they are the one thing this design does not store exactly.
%
% [input]  A, B : 1 x 6 conics [a b c d e f] for a x^2 + b xy + c y^2 + d x + e y + f = 0.
%                 Integer-valued (they come from ratQ.conic); no scale is assumed.
%          opts : (optional) struct
%                   .tol      residual accepted when matching a root back to B (default 1e-9,
%                             scaled by the conic's own magnitude)
%                   .polish   Newton polish steps in 2D (default 8)
% [output] P    : m x 2 real intersection points, m <= 4, in CANONICAL ORDER (lexicographic by
%                 (x,y)). The order is the vertex NAME's third component, so it must depend only
%                 on the two conics and not on how they were built.
%          info : struct with the certificate and the diagnostics the filter needs
%                   .quartic  1 x 5 EXACT integer coefficients of the resultant in the sheared x
%                   .shear    the integer k used (x = x' + k*y'), 0 when none was needed
%                   .resid    m x 1 max(|A|,|B|) residual at each returned point, normalised
%                   .sep      the smallest gap between two returned points (inf when m < 2)
%                   .degenerate  true when the two conics share a component (infinitely many
%                             intersections) -- P is then empty and the caller must not treat
%                             the pair as defining a vertex
%
% WHAT IS EXACT HERE AND WHAT IS NOT, precisely. The quartic's COEFFICIENTS are exact integers: it
% is the Sylvester resultant of two integer polynomials, computed by integer polynomial arithmetic
% with ratQ.chk on every entry. So the certificate that names the vertex is exact, and a future
% degree-<=4 real-algebraic kernel would take exactly this quartic as its input. What is floating
% point is only the ROOT of that quartic, which is what the premise of this design permits:
% vertices are stored as intersections of conics with approximate coordinates as needed.
%
% WHY A RESULTANT AND NOT A PENCIL. Both work. The resultant is chosen because it yields the exact
% integer quartic as a by-product, which is the certificate; the pencil route factors a cubic
% first and its intermediate quantities are not the thing anyone wants to keep.
%
% THE SHEAR, AND WHY IT IS FREE. Eliminating y by a 4x4 Sylvester determinant needs both conics to
% genuinely have degree 2 in y, i.e. a nonzero y^2 coefficient; otherwise the formal determinant
% is not the resultant. Under x = x' + k*y' the y^2 coefficient of [a b c d e f] becomes
% a*k^2 + b*k + c, a quadratic in k that is not identically zero unless a=b=c=0 (a line, handled
% separately). Two conics rule out at most four values of k, so one of k = 0,1,2,3,4 always works,
% and it is an integer so the sheared conics stay integral and the quartic stays exact.

    if nargin < 3, opts = struct(); end
    if ~isfield(opts, 'tol'),    opts.tol = 1e-9;  end
    if ~isfield(opts, 'polish'), opts.polish = 8;  end

    A = ratQ.conic(A);  B = ratQ.conic(B);
    info = struct('quartic', [], 'shear', 0, 'resid', zeros(0,1), 'sep', inf, 'degenerate', false);

    if isequal(A, B)
        info.degenerate = true;  P = zeros(0,2);  return
    end

    % ---- a straight edge is a line, and a line meets a conic in a QUADRATIC -------------------
    % Not an optimisation: the resultant of a degree-1 and a degree-2 polynomial is a different
    % size of determinant, and going through the general path would mean padding a formal degree
    % the polynomial does not have -- which is exactly the failure the shear exists to prevent.
    if all(A(1:3) == 0) || all(B(1:3) == 0)
        if all(A(1:3) == 0), L = A; C = B; else, L = B; C = A; end
        [P, q] = lineMeetsConic(L, C, opts);
        info.quartic = [0 0 q];                    % the quadratic, padded to the quartic slot
    else
        [P, info] = generalMeet(A, B, info, opts);
    end

    % ---- canonical order, and the diagnostics the filter reads --------------------------------
    % Lexicographic by (x,y). This is the vertex NAME's root index, so it must be a function of the
    % two conics alone. Sorting a floating-point pair is safe here precisely because a genuine
    % near-tie is reported through .sep rather than silently resolved -- see QuaCon's use of it.
    if ~isempty(P)
        P = sortrows(round(P, 12));                % round only for the ORDERING key
        [P, info.resid] = polishAndScore(P, A, B, opts);
        if size(P,1) > 1
            D = pdist2rows(P);
            info.sep = min(D);
        end
    end
end

% ==============================================================================================

function [P, q] = lineMeetsConic(L, C, opts)
% objective: the real points where the line d x + e y + f = 0 meets the conic C.
% [output] P : m x 2; q : the exact integer QUADRATIC whose roots are the parameter values
%
% Parametrise the line by its own direction so that no division by a possibly-zero coefficient is
% needed: with n = (d,e) the line is {p0 + t*u} for u = (-e, d) and any p0 on it. Both p0 and u are
% rational, so substituting gives an exact integer quadratic in t.
    d = L(4); e = L(5); f = L(6);
    if abs(d) >= abs(e)
        p0 = [-f/d, 0];  scale = d;               % p0 = (-f/d, 0) is on the line
    else
        p0 = [0, -f/e];  scale = e;
    end
    % Clear p0's denominator: work with P0 = scale*p0 (integral) and remember the scale.
    P0 = ratQ.chk(round(p0 * scale), 'line base point');
    u  = [-e, d];
    % C(P0/scale + t u) = 0, multiplied by scale^2, is an integer quadratic in t.
    q2 = C(1)*u(1)^2 + C(2)*u(1)*u(2) + C(3)*u(2)^2;
    q1 = 2*C(1)*P0(1)*u(1) + C(2)*(P0(1)*u(2) + P0(2)*u(1)) + 2*C(3)*P0(2)*u(2) ...
         + scale*(C(4)*u(1) + C(5)*u(2));
    q0 = C(1)*P0(1)^2 + C(2)*P0(1)*P0(2) + C(3)*P0(2)^2 ...
         + scale*(C(4)*P0(1) + C(5)*P0(2)) + scale^2*C(6);
    q = ratQ.chk([q2 q1 q0], 'line-conic quadratic');
    t = realRootsExact(q);
    P = [P0(1)/scale + t*u(1), P0(2)/scale + t*u(2)];
    P = P(all(isfinite(P), 2), :);
    if size(P,1) == 2 && norm(P(1,:) - P(2,:)) <= opts.tol * max(1, norm(P(1,:)))
        P = P(1,:);                                % a tangency is ONE point, not two
    end
end

function [P, info] = generalMeet(A, B, info, opts)
% objective: the general conic-conic case, via the exact Sylvester resultant in y.
    k = shearFor(A, B);
    info.shear = k;
    As = shearConic(A, k);  Bs = shearConic(B, k);
    q  = resultantY(As, Bs);
    info.quartic = q;
    if all(q == 0)
        % The resultant vanishes identically: the two conics share a whole component. That is not
        % a vertex and must not be reported as four coincident ones.
        info.degenerate = true;  P = zeros(0,2);  return
    end
    xs = realRootsExact(q);
    P = zeros(0,2);
    for i = 1:numel(xs)
        ys = yAt(As, xs(i));
        for j = 1:numel(ys)
            % POLISH BEFORE ACCEPTING, not after. A candidate coming off the quartic can be a
            % long way from the true intersection -- a double root of the quartic is only
            % accurate to about sqrt(eps) even after squarefree reduction removes the worst of
            % it -- so testing the RAW candidate against B rejects genuine intersections. Newton
            % on the pair (As,Bs) is well-conditioned wherever the crossing is transversal, so the
            % refined point is the honest thing to judge. A candidate that is not an intersection
            % at all does not converge and is still rejected.
            pt = newtonOnPair(As, Bs, [xs(i), ys(j)], opts.polish);
            if abs(evalConicAt(Bs, pt)) <= opts.tol * conicScale(Bs, pt) && ...
               abs(evalConicAt(As, pt)) <= opts.tol * conicScale(As, pt)
                P(end+1, :) = pt; %#ok<AGROW>
            end
        end
    end
    if ~isempty(P)
        P = [P(:,1) + k*P(:,2), P(:,2)];           % undo the shear x = x' + k*y'
        P = dedupe(P, opts.tol);
    end
end

function k = shearFor(A, B)
% objective: the smallest k >= 0 for which BOTH sheared conics have a nonzero y^2 coefficient.
% Under x = x' + k*y' the y^2 coefficient of [a b c ...] becomes a*k^2 + b*k + c. Each conic rules
% out at most two k, so k <= 4 always succeeds; the loop is bounded rather than trusting that.
    for k = 0:4
        if (A(1)*k^2 + A(2)*k + A(3)) ~= 0 && (B(1)*k^2 + B(2)*k + B(3)) ~= 0
            return
        end
    end
    error('conicMeet:noShear', ...
        ['no integer shear in 0..4 gives both conics a y^2 term, which is impossible unless one ' ...
         'of them has a=b=c=0 and should have taken the line branch.']);
end

function C = shearConic(C, k)
% objective: substitute x = x' + k*y' into the conic, exactly.
%   a(x'+ky')^2 + b(x'+ky')y' + c y'^2 + d(x'+ky') + e y' + f
    a = C(1); b = C(2); c = C(3); d = C(4); e = C(5); f = C(6);
    C = ratQ.chk([a, 2*a*k + b, a*k^2 + b*k + c, d, d*k + e, f], 'sheared conic');
end

function q = resultantY(A, B)
% objective: Res_y(A, B) -- the EXACT integer quartic in x whose roots are the x-coordinates of the
%            intersections. Both conics must have a nonzero y^2 coefficient (the shear guarantees
%            it), so both are genuinely degree 2 in y and the 4x4 Sylvester determinant IS the
%            resultant rather than a formal expression that can vanish for the wrong reason.
%
% As a polynomial in y, [a b c d e f] is   c*y^2 + (b*x + e)*y + (a*x^2 + d*x + f).
% Each Sylvester entry is therefore a polynomial in x, stored here as a coefficient vector in
% DESCENDING powers, and the determinant is expanded with integer polynomial arithmetic -- no
% floating point and no symbolic engine anywhere in it.
    a1 = {A(3)};                 b1 = {[A(2) A(5)]};      c1 = {[A(1) A(4) A(6)]};
    a2 = {B(3)};                 b2 = {[B(2) B(5)]};      c2 = {[B(1) B(4) B(6)]};
    z  = {0};
    M = {a1{1}, b1{1}, c1{1}, z{1};
         z{1},  a1{1}, b1{1}, c1{1};
         a2{1}, b2{1}, c2{1}, z{1};
         z{1},  a2{1}, b2{1}, c2{1}};
    q = polyDet(M);
    q = trimPoly(q);
    q = ratQ.chk([zeros(1, 5 - numel(q)), q], 'resultant quartic');
end

function d = polyDet(M)
% objective: determinant of a matrix whose entries are polynomials (coefficient vectors), by
%            cofactor expansion. Exact integer arithmetic throughout; 4x4 is 24 terms, so the
%            naive expansion is the right one.
    n = size(M, 1);
    if n == 1, d = M{1,1}; return, end
    d = 0;
    for j = 1:n
        if all(M{1,j} == 0), continue, end
        sub = M(2:end, [1:j-1, j+1:n]);
        term = polyMul(M{1,j}, polyDet(sub));
        if mod(j, 2) == 0, term = -term; end
        d = polyAdd(d, term);
    end
end

function r = polyMul(p, q), r = ratQ.chk(conv(p, q), 'polynomial product'); end

function r = polyAdd(p, q)
    n = max(numel(p), numel(q));
    r = ratQ.chk([zeros(1, n-numel(p)), p] + [zeros(1, n-numel(q)), q], 'polynomial sum');
end

function p = trimPoly(p)
    i = find(p ~= 0, 1);
    if isempty(p) || isempty(i), p = 0; else, p = p(i:end); end
end

function r = realRoots(p)
% objective: the real roots of an exact integer polynomial, as doubles, deduplicated.
% This is the ONE place floating point enters the vertex layer, and it is where the premise of the
% design puts it. The polynomial itself is kept exact by the caller so that a later exact kernel
% can take over here without anything upstream changing.
%
% THE SQUAREFREE STEP IS NOT AN OPTIMISATION. A repeated root is the normal case here, not an
% exotic one: two circles meeting tangentially, or the pair (x^2+4y^2=4, x^2=y^2), give the
% resultant a double root, and `roots` then returns it with an error of about sqrt(eps) -- eight
% digits, not sixteen -- often with a spurious imaginary part of the same size. Measured: the
% four-point case above lost two of its four intersections to exactly that, because the candidate
% was rejected as not lying on the second conic.
%
% Dividing by gcd(p, p') removes every repetition EXACTLY, using integer polynomial arithmetic
% (pseudo-remainder sequence, primitive part at each step), so the roots that remain are simple and
% `roots` is back to full accuracy. It is the same move a real-algebraic kernel would make first,
% which is the point: this routine is the placeholder that kernel replaces, and the exact part of
% it will not have to change.
    p = trimPoly(p);
    if numel(p) <= 1, r = zeros(0,1); return, end
    p = p / max(abs(p));                            % scale only -- the ROOTS are unchanged
    rr = roots(p);
    rr = rr(abs(imag(rr)) <= 1e-7 * max(1, abs(real(rr))));
    r  = sort(real(rr));
    if numel(r) > 1
        keep = [true; diff(r) > 1e-10 * max(1, abs(r(2:end)))];
        r = r(keep);
    end
end

function r = realRootsExact(p)
% objective: realRoots for a polynomial whose coefficients are EXACT INTEGERS, with the repeated
%            roots removed first. Split from realRoots deliberately: the y-solve inside the sweep
%            has a double-precision coefficient (it substitutes a computed x), so the exact
%            squarefree step is meaningless there and gcd would refuse it outright. Only the two
%            certificates -- the resultant quartic and the line-conic quadratic -- are integral,
%            and they are the only two places a repeated root actually costs accuracy.
    r = realRoots(squarefreePart(p));
end

function s = squarefreePart(p)
% objective: p divided by gcd(p, p'), exactly, so every root becomes simple.
% Falls back to p unchanged if the exact division does not verify -- a squarefree part that cannot
% be confirmed is not worth risking, and the caller's residual test still gates every candidate.
    p = trimPoly(p);
    if numel(p) <= 2, s = p; return, end
    try
        g = polyGcdInt(p, polyder(p));
    catch ME
        % COEFFICIENT GROWTH IS A REAL OUTCOME HERE, not a theoretical one: the pseudo-remainder
        % chain squares its leading coefficients, and a resultant quartic with six-digit entries
        % reached 2e23 on the random sweep before the content was removed inside the loop. The
        % guard below stays because correctness must not depend on this step at all -- it buys
        % ACCURACY on repeated roots, and the caller's residual test gates every candidate either
        % way. Anything other than the exactness guard firing is a real defect and is rethrown.
        if ~strcmp(ME.identifier, 'ratQ:overflow'), rethrow(ME); end
        s = p; return
    end
    if numel(g) <= 1, s = p; return, end
    s = polyDivExact(p, g);
    if isempty(s), s = p; end
end

function g = polyGcdInt(p, q)
% objective: the primitive integer gcd of two integer polynomials, by the pseudo-remainder
%            sequence. Exact: no division that is not exact, and the content is removed each step
%            so the coefficients cannot grow the way a naive Euclidean chain makes them.
    p = primitivePart(trimPoly(p));
    q = primitivePart(trimPoly(q));
    while numel(q) > 1 || (numel(q) == 1 && q ~= 0)
        r = pseudoRem(p, q);
        p = q;
        q = primitivePart(trimPoly(r));
        if numel(q) == 1 && q(1) == 0, break, end
    end
    g = primitivePart(p);
end

function r = pseudoRem(p, q)
% objective: the pseudo-remainder of p by q -- lc(q)^k * p reduced by q, all in integers.
    r = trimPoly(p); q = trimPoly(q);
    lq = q(1);
    while numel(r) >= numel(q) && any(r ~= 0)
        shift = numel(r) - numel(q);
        r = ratQ.chk(lq * r - r(1) * [q, zeros(1, shift)], 'pseudo-remainder');
        r = trimPoly(r(2:end));                     % the leading term cancels by construction
        % REMOVE THE CONTENT EVERY STEP, not only at the end. Each reduction multiplies by lc(q),
        % so a chain of them squares the coefficients; measured, a quartic with six-digit entries
        % reached 2e23 in three steps. The gcd is unchanged by a constant factor, so dividing the
        % content out here costs nothing and is what keeps the arithmetic inside 2^53.
        r = primitivePart(r);
        if numel(r) == 1 && r == 0, break, end
    end
end

function p = primitivePart(p)
    p = trimPoly(p);
    g = 0;
    for i = 1:numel(p)
        if p(i) ~= 0, g = gcd(g, abs(p(i))); end
    end
    if g > 1, p = p / g; end
    if p(1) < 0, p = -p; end
end

function s = polyDivExact(p, g)
% objective: p / g when g divides p, verified by multiplying back. Returns [] if it does not.
% Division is done in doubles and then CHECKED exactly, rather than trusted: for the degrees here
% that is both the simplest and the only version whose failure is visible.
    [s, rem] = deconv(p, g);
    s = round(s);
    if any(abs(rem) > 1e-6 * max(1, max(abs(p)))) || ~isequal(trimPoly(conv(s, g)), trimPoly(p))
        s = [];
    end
end

function p = newtonOnPair(A, B, p, steps)
% objective: refine a candidate intersection by Newton on the pair (A,B). See the call site for
%            why this happens BEFORE the acceptance test rather than after it.
    for it = 1:steps
        J = [gradConic(A, p); gradConic(B, p)];
        r = [evalConicAt(A, p); evalConicAt(B, p)];
        if ~all(isfinite(J(:))) || rcond(J) < 1e-14, break, end
        step = (J \ r).';
        if ~all(isfinite(step)), break, end
        p = p - step;
        if norm(step) <= eps * max(1, norm(p)), break, end
    end
end

function ys = yAt(C, x)
% objective: the real y with C(x, y) = 0 for a fixed x. C has a nonzero y^2 coefficient here.
    ys = realRoots([C(3), C(2)*x + C(5), C(1)*x^2 + C(4)*x + C(6)]);
end

function v = evalConicAt(C, p)
    v = C(1)*p(1)^2 + C(2)*p(1)*p(2) + C(3)*p(2)^2 + C(4)*p(1) + C(5)*p(2) + C(6);
end

function s = conicScale(C, p)
% objective: the magnitude the residual of C at p should be compared against, so that `tol` means
%            the same thing for a conic scaled by 1000 and for one scaled by 1.
    s = max(1, max(abs(C)) * max(1, max(abs(p))^2));
end

function P = dedupe(P, tol)
    keep = true(size(P,1), 1);
    for i = 2:size(P,1)
        for j = 1:i-1
            if keep(j) && norm(P(i,:) - P(j,:)) <= tol * max(1, norm(P(j,:)))
                keep(i) = false; break
            end
        end
    end
    P = P(keep, :);
end

function [P, resid] = polishAndScore(P, A, B, opts)
% objective: refine each point with 2D Newton on (A, B) and report its normalised residual.
%
% WHY POLISH AT ALL. The quartic's roots are ill-conditioned exactly where two intersections are
% close, which is the case the filter must be able to SEE. Newton on the pair (A,B) is
% well-conditioned wherever the intersection is transversal, so polishing both improves the point
% and makes the residual a meaningful measurement rather than an artefact of the root-finder.
% Where it is NOT transversal (a tangency) the Jacobian is singular, the step is skipped, and the
% residual stays large -- which is the report the caller needs, not something to hide.
    resid = zeros(size(P,1), 1);
    for i = 1:size(P,1)
        p = P(i,:);
        for it = 1:opts.polish
            J = [gradConic(A, p); gradConic(B, p)];
            r = [evalConicAt(A, p); evalConicAt(B, p)];
            if rcond(J) < 1e-12 || ~all(isfinite(J(:))), break, end
            step = (J \ r).';
            if ~all(isfinite(step)), break, end
            p = p - step;
            if norm(step) <= eps * max(1, norm(p)), break, end
        end
        P(i,:) = p;
        resid(i) = max(abs(evalConicAt(A, p)) / conicScale(A, p), ...
                       abs(evalConicAt(B, p)) / conicScale(B, p));
    end
end

function g = gradConic(C, p)
    g = [2*C(1)*p(1) + C(2)*p(2) + C(4), C(2)*p(1) + 2*C(3)*p(2) + C(5)];
end

function d = pdist2rows(P)
% objective: all pairwise distances between rows, without the Statistics Toolbox.
    n = size(P,1); d = zeros(n*(n-1)/2, 1); t = 0;
    for i = 1:n-1
        for j = i+1:n
            t = t + 1; d(t) = norm(P(i,:) - P(j,:));
        end
    end
end
