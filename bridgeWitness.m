function [W, info] = bridgeWitness(d, N1, m1, N2, m2, Dl, dl, kp, C1, C2, P)
% bridgeWitness  Rational points of a bridge face's region, PROPOSED in doubles and VERIFIED exactly.
%
% objective: the one piece of machinery `bridgeFaceAlg` and `bridgeArcCut` share. The producer needs
%   a witness to read the region's sides and certify the branch index at; the cut needs one for the
%   COMPLEMENTARY arc, to learn which side of its own boundary the face is NOT on. Same search, same
%   verification, different arc -- so it lives in one place rather than twice.
%
% [input]  d, N1, m1, N2, m2 : the two faces over one integer denominator
%          Dl, dl, kp        : the tie's own data (Delta, delta, kappa), also over that denominator
%          C1, C2            : the region's two boundary conics, so a proposal ON one is rejected
%          P                 : k x 3 half-planes the base point must satisfy; EMPTY for the whole tie
% [output] W    : m x 3 rows [xn1 xn2 xd], verified, coarsest denominator first
%          info : .nRoots, .tries, .overflows, .nWitnesses
%
% WHY SEVERAL AND NOT ONE. A witness has to decide everything asked of it, and the branch
% certification can overflow 2^53 where the membership test did not -- so the caller needs a list to
% walk, not a single point to trust. Ordered coarsest-denominator first, because every exact test
% downstream clears xd to a power of the total degree.
%
% WHAT "VERIFIED" MEANS, and it is the whole point of the file. x lies in the region exactly when the
% sweeping segment through it exists: Phi(x, .) -- the tie residual along the family, a quartic in
% the multiplier -- must have a root mu in [-1, 0]. With a pair P given, that root must also put the
% BASE POINT s(mu) inside P, and s(mu) is algebraic, so each row of P is the SIGN OF A POLYNOMIAL AT
% A ROOT, which `ratQ.signAt` decides in integers. Nothing here is ever trusted because it looked
% right in floating point.

    % THE WITNESS MUST BE SWEPT FROM THE ARC THE FACE ACTUALLY OWNS. With a pair given, a tie point
    % outside it produces a point of the UNCUT sweep, and the cut would then read its own sides at a
    % point on the wrong side of itself.
    S = tieSample(Dl, dl, kp, P);
    tries = 0;  overflows = 0;  W = zeros(0,3);  nRoots = [];
    for i = 1:size(S, 1)
        s = S(i,:).';
        x0 = ((N1*s + m1) + (N2*s + m2)) / (2*d);
        if any(~isfinite(x0)), continue, end
        for xdTry = [1 2 4 8 16 32 64]
            tries = tries + 1;
            xnTry = round(x0 * xdTry);
            % A CANDIDATE MAY BE UNDECIDABLE RATHER THAN OUTSIDE, and the two must not be
            % confused. Phi's coefficients carry xd^2 and the Sturm chain squares its leading
            % coefficient at every step, so a large witness can exceed 2^53 while a smaller one
            % decides instantly. An overflow therefore moves to the next candidate and is COUNTED,
            % so that "no witness" can say which of the two reasons it means.
            try
                if signAtRational(C1, xnTry, xdTry) == 0 || signAtRational(C2, xnTry, xdTry) == 0
                    continue                          % ON a boundary: names no side
                end
                if isempty(P)
                    [inR, nR] = memberExact(N1, m1, N2, m2, Dl, dl, kp, d, xnTry, xdTry);
                else
                    [inR, nR] = memberExactArc(N1, m1, N2, m2, Dl, dl, kp, d, xnTry, xdTry, P);
                end
            catch ME
                if ~strcmp(ME.identifier, 'ratQ:overflow'), rethrow(ME), end
                overflows = overflows + 1;  continue
            end
            if inR
                if ~any(all(W == [xnTry(:).', xdTry], 2))
                    W(end+1,:) = [xnTry(:).', xdTry]; %#ok<AGROW>
                    nRoots(end+1) = nR;               %#ok<AGROW>
                end
                break                                 % this tie point has answered; take the next
            end
        end
        if size(W,1) >= 12, break, end
    end
    if ~isempty(W)
        [~, ord] = sortrows([W(:,3), max(abs(W(:,1:2)), [], 2)]);   % coarsest, then smallest
        W = W(ord,:);
        info = struct('nRoots', nRoots(ord(1)), 'tries', tries, 'overflows', overflows, ...
                      'nWitnesses', size(W,1));
        return
    end
    if overflows > 0
        error('bridgeFaceAlg:witnessOverflow', ...
            ['no interior point of the bridge region could be verified: %d of %d candidates ' ...
             'were UNDECIDABLE at 2^53 (the Sturm chain on Phi overflowed) and the rest were ' ...
             'outside. The face is not refused because it is ill-formed but because its ' ...
             'membership test does not fit in a double -- the general answer is the modular ' ...
             'route DECISIONS.md records for ratQ, not a different witness.'], overflows, tries);
    end
    error('bridgeFaceAlg:noWitness', ...
        ['no interior point of the bridge region could be VERIFIED, after %d attempts over %d ' ...
         'tie points. Either the tie curve has no real points -- the two quadratics never tie, ' ...
         'so the pair carries no bridge and should not have reached here -- or the region is ' ...
         'thinner than the coarsest rational grid this routine will use, which the caller ' ...
         'should be told rather than have guessed at.'], tries, size(S,1));
end

% --------------------------------------------------------------------------------------------------
function S = tieSample(Dl, dl, kp, P)
% Real points of the tie { s'Dl s + 2 dl's + 2 kp = 0 }, in doubles, and INSIDE the pair when one is
% given. They only ever PROPOSE, so doubles are the right tool.
%
% SAMPLING THEN FILTERING IS NOT ENOUGH, and the failure is silent rather than loud: a fixed grid
% over a long tie can put a SINGLE point inside a short arc, and that point is then a corner of the
% pair -- whose swept midpoint lies exactly ON a cut boundary, where no side can be read. Measured:
% a 41-point grid left exactly one sample in the arc and the routine reported "no witness" for a
% face with a perfectly good parallelogram of interior. So the arc's own extent is found first and
% the sampling happens INSIDE it.
    if all(Dl(:) == 0)
        % A LINE TIE, parameterised directly -- which also removes the solve-for-s2 blind spot: a
        % VERTICAL tie has no s2 to solve for, and the grid returned nothing at all for one.
        D = dl(1);  E = dl(2);
        if D == 0 && E == 0, S = zeros(0,2); return, end
        w  = [E; -D];
        s0 = -kp * [D; E] / (D^2 + E^2);
        [lo, hi] = arcExtent(s0, w, P);
        if isnan(lo), S = zeros(0,2); return, end
        taus = linspace(lo, hi, 43);
        taus = taus(2:end-1);                         % the ENDS are the cut boundaries; stay inside
        S = (s0 + w*taus).';
        S = S(all(isfinite(S), 2), :);
        return
    end

    if ratQ.detExact(Dl) ~= 0
        c0 = -Dl \ dl;
        v0 = (c0.'*Dl*c0)/2 + dl.'*c0 + kp;
        ev = max(abs(eig(Dl)));
        L  = 4*sqrt(abs(2*v0)/max(ev, eps)) + 4;
        ctr = c0(1);
    else
        L = 20;  ctr = 0;
    end
    for refine = 0:2
        S = conicPoints(Dl, dl, kp, ctr, L, 41 * 10^refine);
        if isempty(P) || isempty(S), break, end
        S = S(all(S*P(:,1:2).' + P(:,3).' >= 0, 2), :);
        if size(S,1) >= 5, break, end
    end
    if ~isempty(S)                                   % middle first: the least extreme proposals
        [~, ord] = sort(vecnorm(S - mean(S,1), 2, 2));
        S = S(ord, :);
    end
end

% --------------------------------------------------------------------------------------------------
function [lo, hi] = arcExtent(s0, w, P)
% The parameter interval of { s0 + tau*w } inside P. NaN when the pair excludes the tie entirely.
%
% THE WINDOW FOLLOWS THE ARC rather than the origin. A fixed +-12 window looks reasonable and is
% wrong: an arc that P places far from the foot of the line falls outside it, the sampler returns
% nothing, and the face is refused for "no witness" when it has a perfectly good region. Only an
% end that P leaves INFINITE needs a window at all, and then it is measured from the other end.
    lo = -inf;  hi = inf;
    for i = 1:size(P,1)
        a = P(i,1:2).';  c = P(i,3);
        A = a.'*w;  B = a.'*s0 + c;
        if abs(A) <= 0
            if B < 0, lo = NaN;  hi = NaN;  return, end
            continue
        end
        if A > 0, lo = max(lo, -B/A); else, hi = min(hi, -B/A); end
    end
    if lo >= hi, lo = NaN;  hi = NaN;  return, end
    span = 24;
    if isinf(lo), lo = ternary(isinf(hi), -12, hi - span); end
    if isinf(hi), hi = ternary(isinf(lo),  12, lo + span); end
end

function v = ternary(c, a, b)
    if c, v = a; else, v = b; end
end

% --------------------------------------------------------------------------------------------------
function S = conicPoints(Dl, dl, kp, ctr, L, n)
% Solving the conic in s2 for each s1 on a grid beats walking rays, which miss a branch no ray
% crosses (DECISIONS.md 2026-09-06 records that costing two wrong verdicts).
    S = zeros(0,2);
    for s1 = ctr + linspace(-L, L, n)
        a = Dl(2,2);
        b = 2*Dl(1,2)*s1 + 2*dl(2);
        c = Dl(1,1)*s1^2 + 2*dl(1)*s1 + 2*kp;
        if abs(a) < eps*max(1, abs(b))
            if abs(b) > 0, S(end+1,:) = [s1, -c/b]; end %#ok<AGROW>
            continue
        end
        disc = b^2 - 4*a*c;
        if disc < 0, continue, end
        r = sqrt(disc);
        S(end+1,:) = [s1, (-b + r)/(2*a)]; %#ok<AGROW>
        if r > 0, S(end+1,:) = [s1, (-b - r)/(2*a)]; end %#ok<AGROW>
    end
    S = S(all(isfinite(S), 2) & all(abs(S) < 1e6, 2), :);
end

% --------------------------------------------------------------------------------------------------
function [inR, n] = memberExactArc(N1, m1, N2, m2, Dl, dl, kp, d, xn, xd, P)
% Is x = xn/xd in the region of the face CUT to the arc inside P? EXACTLY.
%
% Membership in the uncut sweep is a root count; membership in the CUT one asks more of the same
% root: not only must Phi have a root mu in [-1,0], the base point s(mu) it names must lie in P.
% s(mu) is algebraic -- it is a rational function of a root of a quartic -- so every row of P is
% the SIGN OF A POLYNOMIAL AT A ROOT, which `ratQ.signAt` decides in integers. That is the same
% predicate the vertex names rest on (CONJ_FIELD_PROOF.md 8.0), one dimension up.
%
% WHY THIS CANNOT BE LEFT TO THE CUT ROWS. The cut reads its sides AT the witness, so a witness
% outside the arc's sweep would have the cut describe a region the face does not occupy -- measured,
% a witness rounded to the nearest integer landed outside a narrow strip and both cut rows then came
% back as the SAME side, leaving an unbounded region where the face is a parallelogram.
    [Phi, Dt, U1, U2] = bridgePhi(N1, m1, N2, m2, Dl, dl, kp, d, xn, xd);
    Ph = ratQ.polyTrim(fliplr(Phi));
    if numel(Ph) <= 1, inR = false;  n = 0;  return, end
    I = ratQ.isolateRoots(Ph);
    n = 0;  inR = false;
    for r = 1:size(I,1)
        lo = max(I(r,1), -1);  hi = min(I(r,2), 0);
        if lo > hi, continue, end
        if ratQ.countRootsIn(Ph, lo, hi) ~= 1 && ~(polyval(Ph, lo) == 0)
            continue                                  % that root is not in [-1,0]
        end
        n = n + 1;
        iv = [I(r,1), I(r,2)];
        sD = ratQ.signAt(fliplr(Dt), Ph, iv);
        if sD == 0, continue, end                     % det M vanishes there: s(mu) is not a point
        ok = true;
        for i = 1:size(P,1)
            q = vadd(vadd(P(i,1)*U1, P(i,2)*U2), P(i,3)*xd*Dt);
            sq = ratQ.signAt(fliplr(q), Ph, iv);
            % STRICT. A witness ON a facet of the pair is swept to a point ON one of the cut
            % boundaries, where no side can be read -- so the interior is what is demanded, and a
            % candidate that only touches is passed over for the next one.
            if sq * sD <= 0, ok = false; break, end
        end
        if ok, inR = true; return, end
    end
end

% --------------------------------------------------------------------------------------------------
function [inR, n] = memberExact(N1, m1, N2, m2, Dl, dl, kp, d, xn, xd)
% Is x = xn/xd in the region? EXACTLY, by Sturm.
%
% x is swept at parameter mu when s(mu) = adj(M_mu)(d x - c_mu)/det(M_mu) lies on the tie, and mu
% in [-1,0] is the segment (mu = 0 at grad g1, mu = -1 at grad g2). Writing U = xd det(M_mu) s(mu),
% which is integral, the tie equation s'Dl s + 2 dl's + 2 kp = 0 clears to
%
%       Phi(mu) = U' Dl U + 2 xd det(M_mu) dl'U + 2 kp xd^2 det(M_mu)^2 ,
%
% a quartic in mu with integer coefficients, and the question is whether it has a root in [-1,0].
    [Phi, ~] = bridgePhi(N1, m1, N2, m2, Dl, dl, kp, d, xn, xd);
    Ph = fliplr(Phi);                                 % ratQ's univariate convention: highest first
    Ph = ratQ.polyTrim(Ph);
    if numel(Ph) <= 1, inR = false;  n = 0;  return, end
    n = ratQ.countRootsIn(Ph, -1, 0);                 % half-open (-1, 0]
    inR = n > 0 || polyval(Ph, -1) == 0;
end

% --------------------------------------------------------------------------------------------------
function c = vadd(a, b)
    n = max(numel(a), numel(b));
    c = [a, zeros(1, n-numel(a))] + [b, zeros(1, n-numel(b))];
end

function s = signAtRational(C, xn, xd)
% The sign of an integer conic at the rational point xn/xd, decided in integers.
    v = C(1)*xn(1)^2 + C(2)*xn(1)*xn(2) + C(3)*xn(2)^2 + ...
        C(4)*xn(1)*xd + C(5)*xn(2)*xd + C(6)*xd^2;
    s = sign(ratQ.chk(v, 'conic value at the witness'));
end
