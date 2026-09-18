function [face, info] = bridgeFaceAlg(n1, d1, n2, d2, P)
% bridgeFaceAlg  The BRIDGE face of a convex envelope, as an AlgCon payload: the implicit
%   polynomial p(x,y,z), the branch index naming which of its roots is the value, and the region.
%
% objective: fill `AlgCon.FP`/`.FB`. The type has carried that pair since 2026-09-13 and nothing
%   produced one; this is the producer. It promotes three probes into one routine --
%   `.claude/bridgeImplicitProbe.m` (the polynomial), `.claude/bridgeRegionBuild.m` (the two
%   boundary conics) and `.claude/bridgeRegionSides.m` (which side of each, exactly) -- and it
%   GENERALISES the first: that probe eliminated one hand-written instance (q1 = |v|^2/2 against
%   q2 = |v-a|^2 + c) and could not have done the general pair by its own method. See below.
%
% [input]  n1, d1 : the first dual cell's conjugate face, as an integer 1x10 numerator in the
%                   weighted basis over a positive integer denominator, i.e.
%                   g1(s) = (s'N1 s/2 + m1's + q1)/d1 with N1 = [n1(5) n1(6); n1(6) n1(7)],
%                   m1 = [n1(8); n1(9)], q1 = n1(10). The cubic entries n1(1:4) must be zero.
%          n2, d2 : the second, likewise.
%          P      : (optional) the CELL PAIR's k x 3 integer half-planes [a b c] meaning
%                   a*s1 + b*s2 + c >= 0. Given, the face is cut to the tie ARC inside the pair
%                   (bridgeArcCut); omitted, the sweep is over the WHOLE tie -- see the closing
%                   note, and do not use the uncut face inside an assembly.
% [output] face : struct -- one CELL in the sense assembleAlgConCells takes, plus the payload
%                   .FP  the integer polynomial p(x,y,z), as ratQ's dense polyN array in three
%                        variables: FP(i+1,j+1,k+1) multiplies x^i y^j z^k
%                   .FB  the branch index: which ascending real root of FP(x,y,.) is the value
%                   .cells cell array of k x 8 row blocks [a b c d e f sign branch]: the region is
%                        the UNION over the blocks, and each block is an intersection of sign
%                        conditions. One block unless a CONIC tie was cut, which splits the face --
%                        see bridgeArcCut. `branch` 0 is the whole conic, 1 or 2 ONE LINE of a
%                        degenerate one (ratQ.pairLineData)
%                   .wit [xn1 xn2 xd], the rational witness the sides and the branch were read at
%          info : the measurements -- tie conic, degrees, largest coefficient against the 2^53
%                 ceiling, and the root counts each exact decision rested on.
%
% ==================================================================================================
% THE FACE, IN ONE LINE, AND WHY IT IS NOT A FORMULA
%
% Over a cell pair the envelope's contribution is sup_s [ <s,x> - G(s) ] with G = max(g1, g2), and
% on the BRIDGE the supremum is attained on the tie curve { g1 = g2 }, under neither quadratic
% alone. So
%
%       z(x)  =  max { <s,x> - g1(s)  :  g1(s) = g2(s) } ,
%
% a quadratic objective on a CONIC constraint. Its value is algebraic in x of degree up to four --
% measured here: degree 4 in z for an anisotropic pair, 2 for an isotropic one -- which is exactly
% what the `Alg` level exists to hold. Alg.m and DECISIONS.md 2026-09-13.
%
% ==================================================================================================
% THE ELIMINATION IS A DISCRIMINANT, NOT A RESULTANT, AND THAT IS THE WHOLE ROUTINE
%
% `.claude/bridgeImplicitProbe.m` eliminated s1 and s2 from three equations with a 4x4 Sylvester
% determinant, and could only afford that because ITS instance made one of the three linear in s2.
% For a general pair the stationarity relation is genuinely quadratic in both components, the
% honest elimination is a 12x12 Sylvester matrix, and that is not affordable at any coefficient
% size. The route below needs no resultant machinery at all.
%
% LAGRANGE, READ AS A ONE-PARAMETER FAMILY OF CONJUGATES. Relax the tie with a multiplier mu:
%
%       psi(mu) = sup_s [ <s,x> - g1(s) - mu (g1(s) - g2(s)) ]
%               = sup_s [ <s,x> - G_mu(s) ]  =  G_mu*(x) ,      G_mu = (1+mu) g1 - mu g2 ,
%
% and G_mu is just another quadratic, so its conjugate is closed form:
%
%       M_mu = (1+mu) A1 - mu A2 ,  c_mu = (1+mu) b1 - mu b2 ,  k_mu = (1+mu) k1 - mu k2 ,
%       psi(mu) = (x - c_mu)' M_mu^-1 (x - c_mu) / 2  -  k_mu .
%
% The constrained value is a STATIONARY value of psi: by the envelope theorem psi'(mu) is minus the
% constraint residual g1(s(mu)) - g2(s(mu)), so "the tie constraint holds" and "psi is stationary"
% are the SAME equation, and no second relation has to be eliminated. Clearing det(M_mu),
%
%       Psi(x, z, mu) = 2 (d z + k_mu) det(M_mu) - (d x - c_mu)' adj(M_mu) (d x - c_mu)
%
% vanishes exactly when z = psi(mu), so { Psi = 0 and dPsi/dmu = 0 } says z is a stationary value.
% That pair is by definition the DISCRIMINANT of Psi in mu, and Psi is a CUBIC in mu -- so the
% elimination is the four-term discriminant of a cubic,
%
%       FP = Disc_mu(Psi) = 18abce - 4b^3e + b^2c^2 - 4ac^3 - 27a^2e^2 ,   reduced by its
%                                                                         content in z,
%
% with no Sylvester matrix anywhere. Every coefficient is an integer: with d = lcm(d1,d2) and both
% faces put over it, M_mu, c_mu and k_mu are integer-valued, and the 2 d^3 that clears Psi's own
% denominator is a constant, which a discriminant's zero set does not see.
%
% THE CUBIC CAN DROP TO A QUADRATIC, and the condition is meaningful rather than a corner case: the
% mu^3 coefficient is kappa det(Delta) - delta' adj(Delta) delta / 2, which is HALF THE DETERMINANT
% OF THE TIE CONIC'S OWN 3x3 MATRIX. It vanishes exactly when the tie conic is DEGENERATE -- a line
% pair, a point, a double line -- and then Psi is a quadratic in mu and the discriminant is
% c^2 - 4be. Both degrees are built; a Psi linear or constant in mu means there is no tie curve at
% all, and that is refused by name.
%
% ==================================================================================================
% THE REGION, AND WHY IT IS TWO SIGN CONDITIONS
%
% At a tie point the subdifferential of a max is the segment between the two gradients, so the set
% of x this face owns is the union over the tie of conv{ grad g1(s), grad g2(s) }: a family of
% segments swept along the tie. Its two extreme members are the IMAGES of the tie conic under the
% affine maps grad g1 and grad g2, each a conic reached by exact substitution (ratQ.substAffine),
% and `.claude/bridgeRegionCross.m` measured that the family has NO envelope arc -- every traced
% boundary point lies on one of the two constructed conics -- so the region really is an
% intersection of two sign conditions and an H-form cell can hold it.
%
% WHICH SIDE is read at a witness, and the witness is PROPOSED in floating point and VERIFIED
% EXACTLY, this project's standing idiom for a decision that must not rest on a tolerance. The
% verification is a root count: x is in the region exactly when the sweeping segment through it
% exists, i.e. when the tie residual along the family, Phi(x, .) -- a quartic in mu -- has a root
% with mu in [-1, 0] (mu = 0 is grad g1's end of the segment, mu = -1 is grad g2's). Sturm decides
% that over the integers, so no proposal is ever trusted.
%
% THE CUT IS NOT OPTIONAL INSIDE AN ASSEMBLY, and the numbers say why. Without `P` the region is the
% sweep over the WHOLE tie, which is bigger than the face: measured, the uncut region equals the
% envelope at every point whose maximiser lies in the pair (14/14) and is BEATEN at every point
% whose maximiser lies outside it (72/72). Pass `P` and `bridgeArcCut` adds the two boundaries the
% arc's endpoints generate; it answers for a LINE tie -- the common one, since adjacent pieces share
% a Hessian -- and refuses a genuine conic tie by name, where the cut segments are irrational.

    if nargin < 4 || nargin > 5
        error('bridgeFaceAlg:nargin', ...
            ['expects (n1, d1, n2, d2) or (n1, d1, n2, d2, P): two faces, each a numerator and ' ...
             'a denominator, and optionally the cell pair to cut the tie to.']);
    end
    if nargin < 5, P = zeros(0,3); end
    if ~isnumeric(P) || (~isempty(P) && size(P,2) ~= 3)
        error('bridgeFaceAlg:pairShape', 'P must be k x 3 integer half-planes; got %s.', ...
              mat2str(size(P)));
    end
    n1 = checkFace(n1, 'n1');  n2 = checkFace(n2, 'n2');
    d1 = checkDen(d1, 'd1');   d2 = checkDen(d2, 'd2');

    % ---- both faces over ONE integer denominator ----------------------------------------------
    d  = ratQ.chk(lcm(d1, d2), 'common denominator');
    N1 = faceHess(n1) * (d/d1);  m1 = faceLin(n1) * (d/d1);  q1 = n1(10) * (d/d1);
    N2 = faceHess(n2) * (d/d2);  m2 = faceLin(n2) * (d/d2);  q2 = n2(10) * (d/d2);
    ratQ.chk([N1(:); N2(:); m1(:); m2(:); q1; q2].', 'rescaled faces');

    Dl = N1 - N2;  dl = m1 - m2;  kp = q1 - q2;       % the tie's own data: Delta, delta, kappa
    tie = ratQ.diffConic(n1, d1, n2, d2);             % the tie conic, canonically spelled

    % ---- 1. the implicit polynomial -----------------------------------------------------------
    [FP, ipInfo] = implicitPoly(d, N1, m1, q1, N2, m2, q2);

    % ---- 2. the two boundary conics -----------------------------------------------------------
    C1 = imageConic(tie, N1, m1, d, 'first');
    C2 = imageConic(tie, N2, m2, d, 'second');

    % ---- 3. the witnesses, proposed then verified ------------------------------------------------
    [W, wInfo] = bridgeWitness(d, N1, m1, N2, m2, Dl, dl, kp, C1, C2, P);

    % ---- 4. the sides, and the branch index -----------------------------------------------------
    % A WITNESS MUST DECIDE EVERYTHING, so the search does not stop at the first point that is
    % merely IN the region: the branch certification substitutes it into FP and clears xd to the
    % sixth power, which can exceed 2^53 where the membership test did not. A candidate that cannot
    % certify is passed over rather than trusted, and only an empty list is a refusal.
    FB = [];  lastErr = [];
    for w = 1:size(W,1)
        xn = W(w,1:2).';  xd = W(w,3);
        try
            s1 = signAtRational(C1, xn, xd);
            s2 = signAtRational(C2, xn, xd);
            [FB, bInfo] = branchIndex(FP, d, N1, m1, q1, N2, m2, q2, Dl, dl, kp, xn, xd);
            break
        catch ME
            if ~any(strcmp(ME.identifier, {'ratQ:overflow', 'bridgeFaceAlg:branchUndecided'}))
                rethrow(ME)
            end
            lastErr = ME;  FB = [];
        end
    end
    if isempty(FB)
        error('bridgeFaceAlg:branchOverflow', ...
            ['none of the %d verified witnesses could CERTIFY the branch index; the last one ' ...
             'reported "%s". Substituting a witness into FP clears its denominator to the total ' ...
             'degree, and the Sturm chain then squares that, so this is the 2^53 ceiling rather ' ...
             'than a malformed face. The general answer is the modular route DECISIONS.md ' ...
             'records for ratQ.'], size(W,1), lastErr.message);
    end

    base = [C1, s1, 0; C2, s2, 0];
    cutInfo = struct('kind', 'none', 'nCut', 0, 'nCells', 1);
    cells = {base};
    if ~isempty(P)
        [blocks, cutInfo] = bridgeArcCut(tie, P, N1, m1, N2, m2, d, [xn(:).', xd], C1, C2, kp);
        cells = cellfun(@(b) [base; b], blocks, 'UniformOutput', false);
    end
    face = struct('FP', {FP}, 'FB', FB, 'cells', {cells}, 'wit', [xn(:).', xd]);
    info = struct('tie', tie, 'degrees', ipInfo.degrees, 'maxCoeff', ipInfo.maxCoeff, ...
                  'ceiling', ipInfo.maxCoeff / 2^53, 'muDegree', ipInfo.muDegree, ...
                  'zContentRemoved', ipInfo.zContentRemoved, ...
                  'witnessRoots', wInfo.nRoots, 'witnessTries', wInfo.tries, ...
                  'value', bInfo.value, 'nRealZ', bInfo.nRealZ, 'bracket', bInfo.bracket, ...
                  'above', bInfo.above, 'cut', cutInfo.kind, 'nCut', cutInfo.nCut, 'nCells', numel(cells));
end

% ==================================================================================================
function [p, info] = implicitPoly(d, N1, m1, q1, N2, m2, q2)
% Disc_mu(Psi), over the integers -- the header's mathematics, transcribed.
%
% Psi is held as a cell array over powers of mu, ASCENDING, whose entries are polyN arrays in
% (x,y,z). M_mu, c_mu and k_mu carry no x, so they stay plain integer vectors in mu and only the
% two components of (d x - c_mu) ever become polynomials. That is what keeps the arithmetic small.
    X = ratQ.monoN([1 0 0], 1, 3);  Y = ratQ.monoN([0 1 0], 1, 3);
    Z = ratQ.monoN([0 0 1], 1, 3);  U = ratQ.monoN([0 0 0], 1, 3);

    Dl  = N1 - N2;
    e11 = [N1(1,1), Dl(1,1)];  e12 = [N1(1,2), Dl(1,2)];  e22 = [N1(2,2), Dl(2,2)];
    Dt  = ratQ.chk(vadd(conv(e11, e22), -conv(e12, e12)), 'det M_mu');   % degree 2 in mu
    a11 = e22;  a12 = -e12;  a22 = e11;                                  % adj(M_mu), symmetric

    c1 = [m1(1), m1(1)-m2(1)];  c2 = [m1(2), m1(2)-m2(2)];
    kk = [q1, q1-q2];

    W1 = {ratQ.polyAddN(d*X, -c1(1)*U), -c1(2)*U};       % (d x - c_mu), component by component
    W2 = {ratQ.polyAddN(d*Y, -c2(1)*U), -c2(2)*U};

    term1 = muScale(muMul({ratQ.polyAddN(d*Z, kk(1)*U), kk(2)*U}, numCell(Dt, U)), 2);
    term2 = muAdd(muAdd(muMul(numCell(a11, U), muMul(W1, W1)), ...
                        muScale(muMul(numCell(a12, U), muMul(W1, W2)), 2)), ...
                  muMul(numCell(a22, U), muMul(W2, W2)));

    Psi = muTrim(muAdd(term1, muScale(term2, -1)));

    switch numel(Psi) - 1
        case 3
            A = Psi{4};  B = Psi{3};  C = Psi{2};  E = Psi{1};
            p =       pAdd(pMul(18, pMul(A, pMul(B, pMul(C, E)))), ...
                           pMul(-4, pMul(B, pMul(B, pMul(B, E)))));
            p = pAdd(p, pMul(pMul(B, B), pMul(C, C)));
            p = pAdd(p, pMul(-4, pMul(A, pMul(C, pMul(C, C)))));
            p = pAdd(p, pMul(-27, pMul(pMul(A, A), pMul(E, E))));
        case 2
            % The tie conic is DEGENERATE -- see the header. Psi = B mu^2 + C mu + E.
            B = Psi{3};  C = Psi{2};  E = Psi{1};
            p = pAdd(pMul(C, C), pMul(-4, pMul(B, E)));
        otherwise
            error('bridgeFaceAlg:degreeInMu', ...
                ['Psi has degree %d in the multiplier, and a stationary value needs at least ' ...
                 'two. Both faces then share a Hessian AND a linear part, so g1 - g2 is a ' ...
                 'constant: { g1 = g2 } is empty or the whole plane, there is no tie curve, and ' ...
                 'so there is no bridge face.'], numel(Psi) - 1);
    end
    muDegree = numel(Psi) - 1;

    p = ratQ.trimN(ratQ.contentFree(ratQ.chk(p, 'discriminant')));
    if size(p, 3) < 2
        error('bridgeFaceAlg:noZ', ...
            ['the eliminated polynomial has degree 0 in z, so it states a condition on x alone ' ...
             'and names no value. That is a producer defect, not an input one.']);
    end

    % ---- remove the content in z ---------------------------------------------------------------
    % The discriminant vanishes on the face AND on an extraneous locus free of z -- measured on the
    % worked instance, degrees (6,6,2) where the face's own quadric has (2,2,2). That extraneous
    % part is exactly the content of p with respect to z, and the candidate for it is the leading
    % z-coefficient made primitive. The candidate is a GUESS, and polyDivExact VERIFIES it: it
    % raises unless the candidate divides every other z-coefficient exactly. This shortcut, and why
    % it is used instead of ratQ.polyGcdBiv (which overflows at these coefficient sizes), is
    % `.claude/bridgeReduceProbe.m`.
    nz = size(p, 3);
    Cz = cell(1, nz);
    for k = 1:nz, Cz{k} = ratQ.trimN(reshape(p(:,:,k), size(p,1), size(p,2))); end
    g = ratQ.trimN(ratQ.contentFree(Cz{nz}));
    removed = false;
    if ~(isscalar(g) && abs(g) == 1)
        Q = cell(1, nz);  ok = true;
        for k = 1:nz
            if isscalar(Cz{k}) && Cz{k} == 0, Q{k} = 0; continue, end
            try
                Q{k} = ratQ.polyDivExact(Cz{k}, g, 2);
            catch ME
                if ~startsWith(ME.identifier, 'ratQ:'), rethrow(ME), end
                ok = false;  break
            end
        end
        if ok
            p = ratQ.trimN(ratQ.contentFree(fromZSlices(Q)));
            removed = true;
        end
    end

    sz = size(p);  sz(end+1:3) = 1;
    info = struct('degrees', sz - 1, 'maxCoeff', max(abs(p(:))), 'muDegree', muDegree, ...
                  'zContentRemoved', removed);
end

% ==================================================================================================
function C = imageConic(tie, Nk, mk, d, which)
% The conic through grad g_k(tie), by exact substitution.
%
% grad g_k(s) = (Nk s + mk)/d is an invertible affine map when Nk is nonsingular, and the image of
% a conic under one is a conic whose coefficients are a substitution away: inverting
% x = (Nk s + mk)/d gives s = adj(Nk)(d x - mk)/det(Nk), which is exactly ratQ.substAffine's
% (Mn, tn, md). A conic goes through a routine that takes a FACE because the 10-wide weighted basis
% spells the same function -- c5 = 2a, c6 = b, c7 = 2c -- and a function's zero set is its zero set.
    if ratQ.detExact(Nk) == 0
        if all(Nk(:) == 0)
            error('bridgeFaceAlg:vertexCell', ...
                ['the %s face is AFFINE, so its dual cell is a VERTEX cell: grad g is constant, ' ...
                 'its end of the sweep collapses to a POINT, and the region is a fan rather ' ...
                 'than a band between two conics. A legitimate cell, and not this routine''s ' ...
                 'shape.'], which);
        end
        error('bridgeFaceAlg:edgeCell', ...
            ['the %s face has a SINGULAR Hessian, so its dual cell is an EDGE cell: grad g ' ...
             'collapses the plane onto a LINE and that end of the sweep is not a conic. A ' ...
             'legitimate cell, and not this routine''s shape.'], which);
    end
    adjN = [Nk(2,2), -Nk(1,2); -Nk(2,1), Nk(1,1)];
    Mn = ratQ.chk(d * adjN,  'image map');
    tn = ratQ.chk(-adjN*mk,  'image translation');
    md = ratQ.detExact(Nk);
    [n2, ~] = ratQ.substAffine(faceOfConic(tie), 1, Mn, tn, md);
    C = ratQ.conic(conicOfFace(n2));
end


% ==================================================================================================
function [b, info] = branchIndex(FP, d, N1, m1, q1, N2, m2, q2, Dl, dl, kp, xn, xd)
% Which ascending real root of FP(x,y,.) is this face's value at the witness.
%
% THE VALUE IS COMPUTED FROM THE OTHER SIDE OF THE CONSTRUCTION, which is what makes this a check
% rather than a restatement: the mu-roots of Phi in [-1,0] are the stationary multipliers, each
% gives z = psi(mu) directly from the closed form, and the face's value is the largest of them --
% the face IS a supremum. FP is never used to produce the value, only to be indexed by it.
%
% THE ANSWER IS ALWAYS -1, THE TOP BRANCH, AND IT IS CERTIFIED RATHER THAN ASSUMED. The face is a
% supremum over the tie and the real roots of FP are that problem's critical values, so the value
% is the largest of them. That is an argument, so the routine checks it: the proposed value is
% bracketed by a DYADIC interval (exact in binary) and two Sturm counts finish it -- exactly one
% root inside the bracket, and NONE above it. Multiplicity does not enter: Sturm counts distinct
% roots, and "nothing above" is the same statement either way, which is why the top index needs no
% squarefree fibre where an ascending one would.
%
% The ascending index is NOT used, and the reason is measured: over 4758 region points of three
% pairs the value is the largest real root every time, while its ascending index takes both 2 and 4
% WITHIN ONE FACE, because a pair of lower roots leaves the real line in the face's interior.
% algBranch.m and DECISIONS.md 2026-09-17.
    Phi = ratQ.polyTrim(fliplr(bridgePhi(N1, m1, N2, m2, Dl, dl, kp, d, xn, xd)));
    x = xn(:) / xd;
    mus = roots(Phi);
    mus = real(mus(abs(imag(mus)) <= 1e-10*max(1, abs(mus))));
    mus = mus(mus >= -1 - 1e-12 & mus <= 1e-12);
    if isempty(mus)
        error('bridgeFaceAlg:noStationary', ...
            ['the witness was VERIFIED to be in the region, so Phi has a root in [-1,0], yet ' ...
             'none survived the double-precision root finder. The exact test and the numeric ' ...
             'one disagree, which is a defect here and not an input.']);
    end
    zs = arrayfun(@(mu) psiValue(mu, d, N1, m1, q1, N2, m2, q2, x), mus);
    zs = zs(isfinite(zs));
    if isempty(zs)
        error('bridgeFaceAlg:noStationary', ...
            'every stationary multiplier at the witness gave a non-finite value.');
    end
    zStar = max(zs);

    % ---- the fibre FP(x,y,.), exactly, as an integer polynomial in z ---------------------------
    [Pz, ~] = fibre(FP, xn, xd);
    Pz = ratQ.polyTrim(Pz);
    if ~isempty(Pz)
        % The integer content is free to remove and is not free to KEEP: it multiplies every
        % pseudo-remainder the Sturm chain forms, and the chain squares its leading coefficient at
        % every step. Clearing xd to the total degree leaves a content of xd^k for exactly that
        % reason, so this is the difference between deciding and overflowing.
        Pz = Pz / gcdAll(Pz);
    end
    if numel(Pz) <= 1
        error('bridgeFaceAlg:emptyFibre', ...
            'FP degenerates to a constant at the witness, so it names no value there.');
    end
    rz = roots(Pz);
    rz = rz(abs(imag(rz)) <= 1e-10*max(1, abs(rz)));
    rr = sort(real(rz));
    [gap, ~] = min(abs(rr - zStar));
    scale = max(1, abs(zStar));
    if gap > 1e-6 * scale
        error('bridgeFaceAlg:valueNotARoot', ...
            ['the value %.12g computed from the multiplier is not a root of FP at the witness ' ...
             '(the nearest root is %.3e away). The elimination and the closed form disagree, ' ...
             'which is a defect in this routine and not in its input.'], zStar, gap);
    end

    % ---- and now certify the TOP branch exactly --------------------------------------------------
    b  = -1;
    w  = max(1e-7, 1e-7*scale);
    lo = dyadic(zStar - w, 'floor');  hi = dyadic(zStar + w, 'ceil');
    M  = ceil(1 + max(abs(Pz(2:end)))/abs(Pz(1))) + 1;        % Cauchy: every real root is inside
    inBracket = ratQ.countRootsIn(Pz, lo, hi);
    above     = ratQ.countRootsIn(Pz, hi, M);
    if inBracket ~= 1 || above ~= 0
        error('bridgeFaceAlg:branchUndecided', ...
            ['the top branch could not be CERTIFIED at the witness: %d root(s) in the bracket ' ...
             '(needs 1) and %d above it (needs 0). The double-precision value and the exact ' ...
             'count disagree, so no branch index is claimed.'], inBracket, above);
    end

    info = struct('value', zStar, 'nRealZ', numel(rr), 'bracket', [lo hi], 'above', above);
end

% --------------------------------------------------------------------------------------------------
function z = psiValue(mu, d, N1, m1, q1, N2, m2, q2, x)
% psi(mu) = (x - c_mu)' M_mu^-1 (x - c_mu)/2 - k_mu, in the scaled integers: with everything over
% d, this is [ (d x - c_mu)' adj(M_mu) (d x - c_mu) / (2 det M_mu) - k_mu ] / d.
    M = (1+mu)*N1 - mu*N2;
    c = (1+mu)*m1 - mu*m2;
    k = (1+mu)*q1 - mu*q2;
    dt = M(1,1)*M(2,2) - M(1,2)*M(2,1);
    if dt == 0, z = inf; return, end
    adjM = [M(2,2), -M(1,2); -M(2,1), M(1,1)];
    w = d*x - c;
    z = ((w.'*adjM*w)/(2*dt) - k) / d;
end

% --------------------------------------------------------------------------------------------------
function [Pz, deg] = fibre(FP, xn, xd)
% FP(xn/xd, .) as an integer polynomial in z, HIGHEST DEGREE FIRST, cleared by xd^deg where deg is
% the largest total degree in (x,y) that appears. Clearing by one power for the whole polynomial
% rather than per term is what keeps it a polynomial in z with integer coefficients.
    sz = size(FP);  sz(end+1:3) = 1;
    deg = 0;
    for i = 1:sz(1)
        for j = 1:sz(2)
            if any(FP(i,j,:) ~= 0), deg = max(deg, (i-1) + (j-1)); end
        end
    end
    c = zeros(1, sz(3));
    for k = 1:sz(3)
        acc = 0;
        for i = 1:sz(1)
            for j = 1:sz(2)
                v = FP(i,j,k);
                if v == 0, continue, end
                acc = acc + v * xn(1)^(i-1) * xn(2)^(j-1) * xd^(deg-(i-1)-(j-1));
            end
        end
        c(k) = acc;
    end
    Pz = ratQ.chk(fliplr(c), 'fibre in z');
end

% --------------------------------------------------------------------------------------------------
function v = dyadic(x, mode)
% The nearest dyadic rational with 2^-20 spacing, below or above. Dyadic so that the value is
% EXACT as a double: Sturm evaluates its chain at these endpoints, and an endpoint that is not
% exactly representable would be a tolerance smuggled into an exact count.
    s = 2^20;
    if strcmp(mode, 'floor'), v = floor(x*s)/s; else, v = ceil(x*s)/s; end
end

% ==================================================================================================
% small helpers: mu-polynomials with polyN coefficients, and the two bases
% ==================================================================================================
function c = muMul(a, b)
    c = repmat({0}, 1, numel(a) + numel(b) - 1);
    for i = 1:numel(a)
        if isscalar(a{i}) && a{i} == 0, continue, end
        for j = 1:numel(b)
            if isscalar(b{j}) && b{j} == 0, continue, end
            c{i+j-1} = ratQ.polyAddN(c{i+j-1}, ratQ.polyMulN(a{i}, b{j}));
        end
    end
end
function c = muAdd(a, b)
    n = max(numel(a), numel(b));
    c = repmat({0}, 1, n);
    for i = 1:n
        if i <= numel(a), c{i} = ratQ.polyAddN(c{i}, a{i}); end
        if i <= numel(b), c{i} = ratQ.polyAddN(c{i}, b{i}); end
    end
end
function c = muScale(a, k)
    c = a;
    for i = 1:numel(a), c{i} = k * a{i}; end
end
function c = muTrim(a)
    last = 0;
    for i = 1:numel(a)
        if ~(isscalar(a{i}) && a{i} == 0) && any(a{i}(:) ~= 0), last = i; end
    end
    c = a(1:max(last,1));
end
function c = numCell(v, U)
    c = cell(1, numel(v));
    for i = 1:numel(v), c{i} = v(i)*U; end
end
function g = gcdAll(p)
    g = 0;
    for v = p, g = gcd(g, abs(round(v))); end
    if g == 0, g = 1; end
end
function c = vadd(a, b)
    n = max(numel(a), numel(b));
    c = [a, zeros(1, n-numel(a))] + [b, zeros(1, n-numel(b))];
end
function c = pMul(a, b), c = ratQ.polyMulN(a, b); end
function c = pAdd(a, b), c = ratQ.polyAddN(a, b); end

function A = fromZSlices(Q)
% Reassemble a 3-D polyN from its z-slices, zero-padding to the union of their (x,y) sizes.
    nx = 1;  ny = 1;
    for k = 1:numel(Q)
        if isscalar(Q{k}) && Q{k} == 0, continue, end
        nx = max(nx, size(Q{k},1));  ny = max(ny, size(Q{k},2));
    end
    A = zeros(nx, ny, numel(Q));
    for k = 1:numel(Q)
        if isscalar(Q{k}) && Q{k} == 0, continue, end
        A(1:size(Q{k},1), 1:size(Q{k},2), k) = Q{k};
    end
end

function n = faceOfConic(q)
% a conic [a b c d e f] as a FACE in the 10-wide weighted basis: c5 = 2a, c6 = b, c7 = 2c.
    n = [0 0 0 0, 2*q(1), q(2), 2*q(3), q(4), q(5), q(6)];
end
function q = conicOfFace(n)
% TWICE the face, so the result is integral: the weighted basis carries halves on the squared
% terms, and a conic has NO SCALE, so doubling is free where halving is not.
    q = [n(5), 2*n(6), n(7), 2*n(8), 2*n(9), 2*n(10)];
end

function s = signAtRational(C, xn, xd)
% The sign of the integer conic C at the rational point xn/xd, decided in integers.
    v = C(1)*xn(1)^2 + C(2)*xn(1)*xn(2) + C(3)*xn(2)^2 + ...
        C(4)*xn(1)*xd + C(5)*xn(2)*xd + C(6)*xd^2;
    s = sign(ratQ.chk(v, 'conic value at the witness'));
end

function H = faceHess(n), H = [n(5) n(6); n(6) n(7)]; end
function m = faceLin(n),  m = [n(8); n(9)];           end

function n = checkFace(n, name)
    if ~isnumeric(n) || ~isvector(n) || numel(n) ~= 10
        error('bridgeFaceAlg:shape', '%s must be a 1x10 face numerator; got %s.', ...
              name, mat2str(size(n)));
    end
    n = ratQ.chk(double(n(:).'), name);
    if any(n(1:4) ~= 0)
        error('bridgeFaceAlg:cubicFace', ...
            ['%s has a CUBIC part. A conjugate face is quadratic -- the sup of an affine ' ...
             'function minus a quadratic is one -- so a cubic here means the caller passed ' ...
             'something other than a conjugate cell.'], name);
    end
end
function d = checkDen(d, name)
    if ~isscalar(d) || ~isnumeric(d) || d ~= fix(d) || d <= 0
        error('bridgeFaceAlg:denominator', '%s must be a positive integer; got %s.', ...
              name, mat2str(d));
    end
    d = double(d);
end
