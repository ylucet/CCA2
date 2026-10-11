function faces = cornerFaceAlg(n, d, P, tie, side)
% cornerFaceAlg  The candidates a CURVED-tie cell needs beyond the straight ones: the CORNERS, where
%   the tie meets a facet of the cell.
%
% objective: complete the candidate set of `sup { <s,x> - g(s) : s in P, side*q(s) >= 0 }`. The
%   straight candidates -- stationary, facet, vertex -- come from `envelopeCellRegions`, and the
%   curved facet's interior is the bridge face. What is left is where those two boundaries MEET, and
%   it is not a detail: measured on the c = -2 fixture, EVERY uncovered point of both curved pairs
%   (91 of 200 and 109 of 200) has its maximiser at exactly such a corner.
%
% [input]  n, d : the cell's face, as an integer 1x10 numerator over a positive denominator
%          P    : k x 3 integer half-planes of the cell
%          tie  : the tie conic, integer [a b c d e f]
%          side : +1 for the side where tie >= 0, -1 for the other
% [output] faces : struct array, one per corner that can carry a candidate, with
%                    .kind  'corner'
%                    .FP    the implicit polynomial p(x,y,z), integer, degree 2 in z
%                    .FB    the branch index (negative: counted from the top)
%                    .s     the corner in doubles -- the CACHE of the exact name below
%                    .sigma the exact name: the facet's parameter quadratic and which root
%                    .M, .p the membership data: the corner wins where M \\ (p - x) >= 0
%
% ==================================================================================================
% THE VALUE IS AFFINE IN x AND IRRATIONAL IN ITS COEFFICIENTS
%
% The corner `s_e` is a FIXED point -- it does not move with x -- so the candidate's value is simply
%
%       z(x) = <s_e, x> - g(s_e) ,
%
% an AFFINE function of x. What makes it more than an ordinary face is that `s_e` is a quadratic
% irrational: it is the tie conic restricted to a facet line, i.e. a root of an integer quadratic in
% the facet's own parameter. So the coefficients live in Q(sigma), and the face is stored the way
% every irrational value in this toolbox is stored -- implicitly:
%
%       with K z = E0(x) + E1(x) sigma  and  alpha sigma^2 + beta sigma + gamma = 0 ,
%       [ 2 alpha K z - 2 alpha E0 + beta E1 ]^2  =  E1^2 (beta^2 - 4 alpha gamma) ,
%
% a polynomial in (x, y, z) over the INTEGERS, of degree 2 in z, whose two roots are the two
% conjugate corners. `FB` names which. No square root is ever formed.
%
% THE REGION IS A WEDGE, AND MEMBERSHIP IS A PREDICATE. `s_e` is the maximiser exactly when the KKT
% multipliers of the two active constraints -- the facet and the tie -- are both nonnegative:
%
%       x - grad g(s_e) + mu1 a + mu2 side grad q(s_e) = 0 ,   mu >= 0 ,
%
% which is `M \\ (grad g(s_e) - x) >= 0` for the fixed matrix `M = [a, side grad q(s_e)]`. Two linear
% conditions with irrational coefficients: a WEDGE whose apex is irrational. It is carried as a
% predicate for the same reason a bridge face's region is (DECISIONS.md 2026-09-17) -- only the
% conjugate PAIR of such wedges is rational, and a pair does not say which one.
%
% CARRYING IT AS ROWS NOW WORKS, and `.con` holds them: `cornerWedgeRows` spells each edge as one
% BRANCH of its rational conjugate pair, measured against the predicate above at 6400 of 6400 points
% of the c = -2 fixture. The paragraph below is the FIRST attempt at it, kept because the measurement
% is the reason the second one is written the way it is.
%
% CARRYING IT AS ROWS WAS TRIED ONCE AND FAILED. `ratQ.pairLineData` and the BRANCH column
% on `QuaCon.FC` name one line of a rational pair exactly, including the PARALLEL pair this wedge's
% facet edge always is, so the naming was not the obstacle. The rows built from it disagreed with the
% predicate above at 450 of 3200 points of the c = -2 fixture and could not name a branch at all for
% 8 of 16 corners. The step to redo is the reduction of the KKT multipliers modulo the facet's
% parameter quadratic, not the naming. DECISIONS.md 2026-09-18.

    if nargin ~= 5
        error('cornerFaceAlg:nargin', 'expects (n, d, P, tie, side).');
    end
    n = ratQ.chk(double(n(:).'), 'face');
    if numel(n) ~= 10 || any(n(1:4) ~= 0)
        error('cornerFaceAlg:shape', 'n must be a 1x10 QUADRATIC face numerator.');
    end
    if ~ismember(side, [-1 1])
        error('cornerFaceAlg:side', 'side must be +1 or -1; got %s.', mat2str(side));
    end
    tie = ratQ.chk(tie(:).', 'tie conic');

    N = [n(5) n(6); n(6) n(7)];  m = [n(8); n(9)];  qq = n(10);
    % `src` is what the face was BUILT from, kept so a consumer can redo the construction exactly
    % rather than re-deriving it from doubles -- `cornerWedgeRows` needs precisely this.
    faces = struct('kind', {}, 'FP', {}, 'FB', {}, 's', {}, 'sigma', {}, 'M', {}, 'p', {}, ...
                   'src', {}, 'con', {});

    for i = 1:size(P,1)
        a = P(i,1:2).';  c = P(i,3);
        na = ratQ.chk(a.'*a, 'facet normalisation');
        if na == 0, continue, end
        w   = [a(2); -a(1)];
        s0n = [-c*a(1); -c*a(2)];
        [al, be, ga] = restrictTie(tie, s0n, w, na);
        if al == 0
            % LINEAR along the facet: one crossing, at a RATIONAL point. This is the GENERIC case
            % for a vertex against an edge -- the parabola's axis is parallel to the strip's
            % facets -- and skipping it lost every such corner (2026-10-10). be = 0 too means the
            % facet runs along the tie or misses it: no corner.
            if be == 0, continue, end
            f1 = rationalCorner(n, d, N, m, qq, P, a, s0n, w, na, be, ga, tie, side, i);
            faces = [faces, f1]; %#ok<AGROW>
            continue
        end
        D = ratQ.chkB(be^2 - 4*al*ga, be^2 + 4*abs(al)*abs(ga), 'corner discriminant');
        if D <= 0, continue, end                   % the facet misses the tie, or touches it

        % ---- the value, implicitly ------------------------------------------------------------
        [E0, E1, K] = valueForm(N, m, qq, d, s0n, w, na, al, be, ga);
        FP = implicitValue(E0, E1, K, al, be, D);

        for root = [1 -1]                          % sigma = (-be + root*sqrt(D)) / (2 al)
            sg = (-be + root*sqrt(D)) / (2*al);
            s  = (s0n + sg*w) / na;
            if ~all(P(:,1:2)*s + P(:,3) >= -1e-9 * max(1, max(abs(s))))
                continue                            % the corner is not in the cell at all
            end
            gq = [2*tie(1)*s(1) + tie(2)*s(2) + tie(4); tie(2)*s(1) + 2*tie(3)*s(2) + tie(5)];
            M  = [a, side*gq];
            if abs(det(M)) < 1e-12 * max(1, max(abs(M(:))))
                continue                            % the two active constraints are parallel there
            end
            p  = (N*s + m) / d;                     % grad g at the corner
            [FB, ok] = branchOf(FP, s, p, d, N, m, qq);
            if ~ok, continue, end
            face = struct('kind', 'corner', 'FP', FP, 'FB', FB, 's', s, ...
                          'sigma', struct('quad', [al be ga], 'root', root, ...
                                          'facet', i, 's0n', s0n, 'w', w, 'na', na), ...
                          'M', M, 'p', p, ...
                          'src', struct('n', n, 'd', d, 'tie', tie, 'side', side), 'con', []);
            % AND THE SAME WEDGE AS ROWS, when it can be spelled: that is what lets a corner enter an
            % arrangement instead of being carried as a predicate. A refusal leaves `.con` empty and
            % costs nothing -- the predicate above is unchanged and remains the definition.
            face.con = cornerWedgeRows(face, n, d, tie, side);
            faces(end+1) = face; %#ok<AGROW>
        end
    end
end

% ==================================================================================================
function faces = rationalCorner(n, d, N, m, qq, P, a, s0n, w, na, be, ga, tie, side, i)
% The corner of a facet the tie crosses ONCE: sigma = -ga/be, so s = sn/sd is rational, the face is
% the AFFINE function <s,x> - g(s) with rational coefficients, and its wedge is two straight rows.
%
%   s = (s0n*be - ga*w) / (na*be),   and   2 d sd^2 z - 2 d sd <sn,x> + sn'N sn + 2 sd m'sn
%   + 2 sd^2 qq = 0   is the value, cleared to integers (the implicit form every face carries).
    faces = struct('kind', {}, 'FP', {}, 'FB', {}, 's', {}, 'sigma', {}, 'M', {}, 'p', {}, ...
                   'src', {}, 'con', {});
    % Every sum of products goes through ratQ.chkB with its absolute-value bound: a cancelling
    % sum past 2^53 can come out small and wrong while passing a check on itself.
    A = @abs;
    sn = ratQ.chkB(s0n*be - ga*w, A(s0n)*A(be) + A(ga)*A(w), 'rational corner');
    sd = ratQ.chk(na*be, 'rational corner');
    if sd < 0, sn = -sn;  sd = -sd; end
    gs = gcd(gcd(sn(1), sn(2)), sd);                           % LOWEST TERMS before anything is
    sn = sn / gs;  sd = sd / gs;                               % squared: erratum Example 2, 2026-10-11
    inCell = ratQ.chkB(P(:,1:2)*sn + P(:,3)*sd, A(P(:,1:2))*A(sn) + A(P(:,3))*sd, ...
                       'rational corner in the cell');
    if any(inCell < 0), return, end                            % exact: the corner is not in the cell
    FP = zeros(2,2,2);
    FP(1,1,2) = ratQ.chk(2*d*sd^2, 'rational corner value');
    FP(2,1,1) = ratQ.chk(-2*d*sd*sn(1), 'rational corner value');
    FP(1,2,1) = ratQ.chk(-2*d*sd*sn(2), 'rational corner value');
    FP(1,1,1) = ratQ.chkB(sn.'*N*sn + 2*sd*(m.'*sn) + 2*sd^2*qq, ...
                          A(sn).'*A(N)*A(sn) + 2*sd*(A(m).'*A(sn)) + 2*sd^2*A(qq), ...
                          'rational corner value');
    FP = ratQ.trimN(ratQ.contentFree(FP));
    % the tie's gradient at s, times sd (a positive scaling of M's second column changes no sign)
    T = [2*tie(1) tie(2) tie(4); tie(2) 2*tie(3) tie(5)];
    gqn = ratQ.chkB(T*[sn; sd], A(T)*[A(sn); sd], 'rational corner gradient');
    Mn = [a, side*gqn];
    dM = ratQ.chkB(Mn(1,1)*Mn(2,2) - Mn(1,2)*Mn(2,1), A(Mn(1,1)*Mn(2,2)) + A(Mn(1,2)*Mn(2,1)), ...
                   'rational corner det');
    if dM == 0, return, end                                    % the two constraints are parallel
    pn = ratQ.chkB(N*sn + m*sd, A(N)*A(sn) + A(m)*sd, 'rational corner gradient of g');
    pd = ratQ.chk(d*sd, 'rational corner');                    % grad g(s) = pn/pd
    % mu = Mn \ (p - x) >= 0  <=>  sign(dM) * adj(Mn) (pn - pd x) >= 0, two rows LINEAR in x
    adjM = [Mn(2,2), -Mn(1,2); -Mn(2,1), Mn(1,1)];
    con = zeros(2,8);
    for r = 1:2
        R = [ratQ.chk(-sign(dM)*pd*adjM(r,1:2), 'rational corner row'), ...
             ratQ.chkB(sign(dM)*(adjM(r,:)*pn), A(adjM(r,:))*A(pn), 'rational corner row')];
        con(r,1:6) = ratQ.conic([0 0 0, R]);
        k = find([0 0 0, R] ~= 0, 1);
        con(r,7) = sign(R(k-3));
    end
    s = sn / sd;  p = pn / pd;
    M = [a, side*gqn/sd];
    face = struct('kind', 'corner', 'FP', FP, 'FB', -1, 's', s, ...
                  'sigma', struct('quad', [0 be ga], 'root', 0, 'facet', i, 's0n', s0n, 'w', w, ...
                                  'na', na), ...
                  'M', M, 'p', p, 'src', struct('n', n, 'd', d, 'tie', tie, 'side', side), ...
                  'con', con);
    faces(end+1) = face;
end

% --------------------------------------------------------------------------------------------------
function [al, be, ga] = restrictTie(tie, s0n, w, na)
% The tie along s = (s0n + sigma*w)/na, cleared by na^2: an integer quadratic in sigma.
    A = tie(1); B = tie(2); C = tie(3); Dc = tie(4); E = tie(5); F = tie(6);
    al = A*w(1)^2 + B*w(1)*w(2) + C*w(2)^2;
    be = 2*A*s0n(1)*w(1) + B*(s0n(1)*w(2) + s0n(2)*w(1)) + 2*C*s0n(2)*w(2) + na*(Dc*w(1) + E*w(2));
    ga = A*s0n(1)^2 + B*s0n(1)*s0n(2) + C*s0n(2)^2 + na*(Dc*s0n(1) + E*s0n(2)) + F*na^2;
    al = ratQ.chk(al, 'restricted tie');  be = ratQ.chk(be, 'restricted tie');
    ga = ratQ.chk(ga, 'restricted tie');
end

% --------------------------------------------------------------------------------------------------
function [E0, E1, K] = valueForm(N, m, qq, d, s0n, w, na, al, be, ga)
% K*z = E0(x) + E1(x)*sigma, with E0 and E1 AFFINE in x and everything integral.
%
% z = <s,x> - g(s) at s = (s0n + sigma w)/na. The only obstacle is sigma^2, and the tie's own
% quadratic removes it: alpha sigma^2 = -beta sigma - gamma, so multiplying through by alpha leaves
% a form that is LINEAR in sigma.
    T0 = al*(s0n.'*N*s0n) - ga*(w.'*N*w) + 2*na*al*(m.'*s0n) + 2*al*na^2*qq;
    T1 = 2*al*(s0n.'*N*w) - be*(w.'*N*w) + 2*na*al*(m.'*w);
    E0 = [2*al*na*d*s0n(1), 2*al*na*d*s0n(2), -T0];      % [gx gy const], affine in x
    E1 = [2*al*na*d*w(1),   2*al*na*d*w(2),   -T1];
    K  = ratQ.chk(2*al*na^2*d, 'value scaling');
    E0 = ratQ.chk(E0, 'value form');  E1 = ratQ.chk(E1, 'value form');
end

% --------------------------------------------------------------------------------------------------
function FP = implicitValue(E0, E1, K, al, be, D)
% [2 al K z - 2 al E0 + be E1]^2 - E1^2 D = 0, as ratQ's dense polyN array in (x, y, z).
    L = affPlusZ(-2*al*E0 + be*E1, 2*al*K);              % the bracket: affine in x, linear in z
    try
        FP = ratQ.polyAddN(sq(L), -D * sq(affPlusZ(E1, 0)));
        FP = ratQ.trimN(ratQ.contentFree(ratQ.chk(FP, 'corner polynomial')));
    catch ME
        % THE SQUARE CAN PASS 2^53 WHILE THE POLYNOMIAL, CONTENT REMOVED, DOES NOT -- biconjQ of
        % the poly2 kink at c = -2 reached it (2026-09-25). Then the square is taken in BigInteger,
        % the content removed THERE, and only the result converted back; fromBig still raises if
        % the reduced polynomial itself does not fit, which is an honest overflow.
        if ~strcmp(ME.identifier, 'ratQ:overflow'), rethrow(ME), end
        % THE BRACKET ITSELF CAN PASS 2^53 TOO -- erratum Example 2 reached 1.8e17 (2026-10-11) --
        % so it is built in BigInteger from its factors, every one of which already fits.
        big = @(v) ratQ.toBig(v);
        Lb = ratQ.bigAdd(ratQ.bigAdd( ...
                ratQ.bigMul(big(-2*ratQ.chk(al, 'corner tie')), big(affPlusZ(E0, 0))), ...
                ratQ.bigMul(big(be), big(affPlusZ(E1, 0)))), ...
                ratQ.bigMul(big(2*al), big(affPlusZ([0 0 0], K))));
        Eb = ratQ.toBig(affPlusZ(ratQ.chk(E1, 'corner slope'), 0));
        Db = ratQ.toBig(ratQ.chk(D, 'corner discriminant'));
        FPb = ratQ.bigAdd(ratQ.bigMul(Lb, Lb), ratQ.bigNeg(ratQ.bigMul(Db, ratQ.bigMul(Eb, Eb))));
        FP = ratQ.trimN(ratQ.fromBig(ratQ.bigContentFree(FPb)));
    end
end

function A = affPlusZ(e, kz)
% (e1 x + e2 y + e3) + kz*z, as a polyN in three variables.
    A = zeros(2,2,2);
    A(2,1,1) = e(1);  A(1,2,1) = e(2);  A(1,1,1) = e(3);  A(1,1,2) = kz;
end

function B = sq(A)
    B = ratQ.polyMulN(A, A);
end

% --------------------------------------------------------------------------------------------------
function [FB, ok] = branchOf(FP, s, p, d, N, m, qq) %#ok<INUSD>
% WHICH root of FP(x,y,.) this corner is. The two roots are the two CONJUGATE corners, so they are
% told apart by evaluating at a point where the face's own value is known: the corner's value at
% x = grad g(s) is <s, grad g(s)> - g(s), which is elementary.
    x = p(:).';
    z = s.'*p - ((s.'*N*s)/2 + m.'*s + qq)/d;
    [~, info] = algBranch(FP, x, -1);
    rr = info.allReal{1};
    if isempty(rr), FB = 0;  ok = false;  return, end
    [gap, j] = min(abs(rr - z));
    if gap > 1e-6 * max(1, abs(z)), FB = 0;  ok = false;  return, end
    FB = j - numel(rr) - 1;                 % negative: counted from the TOP, as AlgCon.FB expects
    ok = true;
end
