function faces = envelopeCellRegions(N, m, q, dd, P)
% envelopeCellRegions  The PRIMAL subdivision of one cell's contribution to the convex envelope,
%   constructed exactly -- each candidate's region as integer half-planes in x, and the quadratic
%   it carries there.
%
% objective: the second step of the envelope's assembly, and the one that was only ever MEASURED
%   before. `.claude/cellFaceExactProbe.m` names the candidate that wins AT a point; an assembly
%   needs the REGION each candidate owns, because that region is the face's domain.
%   `.claude/cellFaceRegionProbe.m` traced those regions along rays and found them POLYHEDRAL
%   (25 arcs, worst line residual 1.3e-3 against a trace step of 4.5e-3, nothing above tolerance).
%   This CONSTRUCTS them instead, so the answer is exact and no longer rests on a fit.
%
% [input]  N  : 2x2 integer, the face's quadratic part in the weighted basis ([n5 n6; n6 n7])
%          m  : 2x1 integer linear part ([n8; n9])
%          q  : integer constant (n10)
%          dd : positive integer denominator -- together g(s) = (s'Ns/2 + m's + q)/dd
%          P  : k x 3 integer half-planes of the cell, rows [a b c] meaning a*s1 + b*s2 + c >= 0
% [output] faces : struct array, one per candidate that owns a 2-D region:
%                    .who  struct('kind', 'stationary'|'facet'|'vertex', 'facet', i, 'facet2', j)
%                    .R    r x 3 integer rows [a b c] meaning a*x1 + b*x2 + c >= 0 -- the region
%                    .fN   1x10 integer numerator in the 10-wide weighted basis
%                    .fD   positive integer denominator, so the face is ratQ.evalFace(fN, fD, x)
%                    .S    the MAXIMISER that produced the value, as the integer affine map
%                          s*(x) = (S.Mn*x + S.tn)/S.d with S.d > 0
%
% WHY THE MAXIMISER IS RETURNED AND NOT JUST THE VALUE. On a BRIDGE pair the sup is taken over one
% SIDE of the tie conic, not over the whole cell, and a candidate is admissible exactly where its
% own maximiser is on that side -- i.e. where q(s*(x)) has the right sign. Since s*(x) is affine,
% that pulls back to a CONIC condition in x (ratQ.substAffine), which an H-form cell can carry and
% which is decided in integers. Without the map the caller would have to re-derive it per candidate
% kind, and the KKT case analysis that produced it lives here.
%
% ------------------------------------------------------------------------------------------------
% THE MATHEMATICS, because the whole routine is one KKT condition written three times.
%
% Over the cell P the envelope's contribution is the concave program
%
%       h(x) = sup { <s,x> - g(s) : s in P } ,     g(s) = ( s'Ns/2 + m's + q ) / dd .
%
% Clearing dd and writing Y = dd*x - m -- INTEGER-affine in x, which is what keeps everything below
% inside the integers -- this is
%
%       dd * h(x) = sup { <s,Y> - s'Ns/2 - q : s in P } ,        phi(s) := <s,Y> - s'Ns/2 - q .
%
% phi is concave with grad phi(s) = Y - Ns, so s0 is the maximiser exactly when
%
%       Y - N s0  lies in the NORMAL CONE of P at s0 .                                    (KKT)
%
% A maximiser lies in the relative interior of exactly one face of P, and that face's normal cone
% is fixed, so there are three cases and in each one BOTH SIDES OF (KKT) ARE AFFINE IN x:
%
%   STATIONARY (the 2-D face). The cone is {0}, so s0 = N^-1 Y and the condition is just
%       s0 in P -- affine in x once cleared by det(N) > 0.
%
%   FACET i (a 1-D face, normal cone = the ray spanned by the OUTWARD normal -a_i). Parametrise the
%       facet line as s = (P0 + tau*w)/na with a_i'P0 = -c_i*na, w = [a_i2; -a_i1] and na = a_i'a_i,
%       all integers -- the scaling by na is what avoids a denominator. Along it phi is a concave
%       1-D quadratic whose stationary point is tau* = U/DND with DND = w'Nw > 0 and
%       U = na*(w'Y) - w'N P0, affine in x. The conditions are then
%         (a) tau* lies within the facet's extent -- one affine row per OTHER constraint, and
%         (b) the multiplier is nonnegative, i.e. (Y - N s0)'a_i <= 0,
%       each cleared by DND > 0 into an integer row.
%
%   VERTEX v (a 0-D face, normal cone spanned by the two active outward normals). With
%       [mu_i; mu_j] = -Acol^-1 (Y - N v) and Acol = [a_i a_j] by COLUMNS, the condition mu >= 0 is
%       two affine rows once cleared by det(Acol), whose SIGN must be carried -- clearing by a
%       negative determinant without flipping the inequality turns a region inside out.
%       TWO MATRICES APPEAR HERE AND THEY ARE TRANSPOSES: v itself solves the ROW system
%       [a_i'; a_j'] v = [-c_i; -c_j], while the cone combines the COLUMNS. Using one for the other
%       still type-checks and still yields a plausible vertex -- it put a cone's apex at (2,0)
%       instead of (2,2), silently dropping that face's whole x2 term.

% WHERE THE REGIONS DO *NOT* COVER THE PLANE, and it is the answer rather than a gap. A cell of a
% conjugate is generally an unbounded cone, so for many x the sup diverges and h(x) = +infinity --
% there is no maximiser, no candidate and correctly no region. The union of the regions is exactly
% the set where the sup is finite, which is why a coverage count below 100% is expected and why the
% test asserts the COMPLEMENT: off the regions, the objective is unbounded above on the cell.
%
% WHY IT MATTERS THAT THESE ARE EXACT. The value at a point may be a double -- CONJ_FIELD_PROOF.md
% 8.0 lists the operations that must leave Q and evaluating a finished function is not among them.
% But a REGION is a mesh decision: it says which formula holds where, and a tolerance there returns
% a plausible number from the wrong face. So every row here is an integer row and every comparison
% that built it was a comparison against zero.
%
% WHAT IS NOT HERE. The three cases above cover a cell of a conjugate of a CONVEX piece, which is
% what the envelope's assembly hands in. A cell whose face is affine along a facet (w'Nw = 0) has
% no facet-interior maximiser -- the max sits at an endpoint, and the vertex case already visits it
% -- so that facet is skipped rather than special-cased. Regions with empty interior are dropped by
% `ratQ.feasible2`, so a candidate that never wins contributes no face.

    if nargin ~= 5
        error('envelopeCellRegions:nargin', 'expects (N, m, q, dd, P).');
    end
    if ~isequal(size(N), [2 2]) || N(1,2) ~= N(2,1)
        error('envelopeCellRegions:shape', 'N must be a symmetric 2x2; got %dx%d.', size(N,1), size(N,2));
    end
    if dd <= 0
        error('envelopeCellRegions:denominator', 'dd must be positive; got %g.', dd);
    end
    m = m(:);
    faces = struct('who', {}, 'R', {}, 'fN', {}, 'fD', {}, 'S', {});
    k = size(P,1);
    if k == 0, return, end

    detN = N(1,1)*N(2,2) - N(1,2)*N(2,1);
    K    = [N(2,2), -N(1,2); -N(1,2), N(1,1)];        % adj(N), so N*K = detN*I

    % ---- 1. THE STATIONARY CANDIDATE ----------------------------------------------------------
    % s0 = N^-1 Y = K*Y/detN, feasible iff P(:,1:2)*s0 + P(:,3) >= 0. Cleared by detN > 0.
    if detN > 0
        A = P(:,1:2) * K;                              % k x 2
        R = [dd * A, -A*m + P(:,3)*detN];
        if ratQ.feasible2(R, true)
            % dd*h = Y'KY/(2 detN) - q, and Y = dd*x - m
            Km = K*m;
            fN = zeros(1,10);
            fN(5)  = 2 * dd^2 * K(1,1);
            fN(6)  = 2 * dd^2 * K(1,2);
            fN(7)  = 2 * dd^2 * K(2,2);
            fN(8)  = -2 * dd * Km(1);
            fN(9)  = -2 * dd * Km(2);
            fN(10) = m.'*Km - 2*q*detN;
            % s*(x) = K*Y/detN with Y = dd*x - m
            S = mkMap(dd*K, -Km, detN);
            faces(end+1) = mkFace('stationary', 0, 0, R, fN, 2*detN*dd, S); %#ok<AGROW>
        end
    end

    % ---- 2. ONE CANDIDATE PER FACET -----------------------------------------------------------
    for i = 1:k
        a = P(i,1:2).';  c = P(i,3);
        na = a.'*a;
        if na == 0, continue, end
        w  = [a(2); -a(1)];                            % along the facet
        DND = w.'*N*w;
        if DND <= 0, continue, end                     % phi affine along the facet: see the header
        P0 = -c*a;                                     % a'P0 = -c*na, so P0/na is on the line

        % U = na*(w'Y) - w'N P0 = alpha'x + beta,  affine in x
        alpha = na * dd * w;
        beta  = -na*(w.'*m) - w.'*N*P0;

        R = zeros(0,3);
        dead = false;
        % (a) the facet's extent: for every OTHER row j, tau*Aj + Bj >= 0, cleared by DND > 0
        for j = 1:k
            if j == i, continue, end
            aj = P(j,1:2).';  cj = P(j,3);
            Aj = aj.'*w;  Bj = aj.'*P0 + cj*na;
            if Aj == 0
                % The row does not vary along the facet, so it holds everywhere on that line or
                % nowhere. NOWHERE means the facet is not part of P's boundary at all and carries
                % no candidate -- flagged rather than left to an empty R, because `zeros(0,3)` is
                % itself `isempty` and testing for that would also discard the legitimate case
                % where no row constrains the facet.
                if Bj < 0, dead = true; break, end
                continue
            end
            R(end+1,:) = Aj*[alpha.', beta] + [0 0 Bj*DND]; %#ok<AGROW>
        end
        if dead, continue, end
        % (b) the multiplier is nonnegative: (Y - N s0)'a <= 0, cleared by DND > 0 and negated
        %     into the >= 0 convention every row here uses.
        NaV = N*a;
        lhsX = DND * na * dd * a.';                    % DND*na*(a'Y) part, in x
        lhsC = -DND * na * (a.'*m) - DND * (a.'*(N*P0));
        row  = -( [lhsX, lhsC] - (w.'*NaV) * [alpha.', beta] );
        R(end+1,:) = row; %#ok<AGROW>

        if ~ratQ.feasible2(R, true), continue, end

        % dd*h = P0'Y/na - P0'NP0/(2na^2) - q + U^2/(2 na^2 DND); cleared by 2 na^2 DND.
        fN = zeros(1,10);
        fN(5)  = 2 * alpha(1)^2;
        fN(6)  = 2 * alpha(1) * alpha(2);
        fN(7)  = 2 * alpha(2)^2;
        fN(8)  = 2*na*DND*dd*P0(1) + 2*beta*alpha(1);
        fN(9)  = 2*na*DND*dd*P0(2) + 2*beta*alpha(2);
        fN(10) = -2*na*DND*(P0.'*m) - DND*(P0.'*N*P0) - 2*na^2*DND*q + beta^2;
        % s*(x) = (P0 + tau* w)/na with tau* = U/DND and U = alpha'x + beta
        S = mkMap(w*alpha.', P0*DND + beta*w, na*DND);
        faces(end+1) = mkFace('facet', i, 0, R, fN, 2*na^2*DND*dd, S); %#ok<AGROW>
    end

    % ---- 3. ONE CANDIDATE PER VERTEX ----------------------------------------------------------
    for i = 1:k
        for j = i+1:k
            ai = P(i,1:2).';  aj = P(j,1:2).';
            % TWO MATRICES, AND THEY ARE TRANSPOSES OF EACH OTHER -- keeping them apart by name is
            % the point. The vertex SOLVE stacks the rows (ai'v + ci = 0 is a row equation), while
            % the multiplier condition Y - Nv = -mu_i ai - mu_j aj combines the COLUMNS. Using one
            % where the other belongs still type-checks, still produces a plausible vertex, and put
            % the apex of a cone at (2,0) instead of (2,2) -- dropping the whole x2 term of that
            % face's value. Caught by the no-overshoot check, which needs no oracle at all.
            Mrow = [ai.'; aj.'];                       % the SOLVE
            Acol = [ai,   aj  ];                       % the normal CONE
            D2 = Mrow(1,1)*Mrow(2,2) - Mrow(1,2)*Mrow(2,1);
            if D2 == 0, continue, end                  % parallel rows name no vertex
            % v = V0/vd solves Mrow * v = [-ci; -cj]
            rhs = [-P(i,3); -P(j,3)];
            adjM = [Mrow(2,2), -Mrow(1,2); -Mrow(2,1), Mrow(1,1)];
            V0  = adjM * rhs;
            vd  = D2;
            if vd < 0, V0 = -V0;  vd = -vd; end
            if ~all(P(:,1:2)*V0 + P(:,3)*vd >= 0), continue, end   % not a vertex OF P

            % mu = -Acol^-1 (Y - N v) >= 0, cleared by det(Acol) with the SIGN carried -- clearing
            % by a negative determinant without flipping the inequality turns a region inside out.
            adjA2 = [Acol(2,2), -Acol(1,2); -Acol(2,1), Acol(1,1)];
            sgn = sign(D2);
            % Y - N v = (vd*Y - N V0)/vd, and vd > 0, so clear it away
            MX = vd * dd * eye(2);                     % the x part of (vd*Y)
            MC = -vd*m - N*V0;                         % the constant part
            R = -sgn * [adjA2 * MX, adjA2 * MC];
            if ~ratQ.feasible2(R, true), continue, end

            % dd*h = V0'Y/vd - V0'NV0/(2 vd^2) - q; cleared by 2 vd^2
            fN = zeros(1,10);
            fN(8)  = 2 * vd * dd * V0(1);
            fN(9)  = 2 * vd * dd * V0(2);
            fN(10) = -2*vd*(V0.'*m) - V0.'*N*V0 - 2*q*vd^2;
            % s*(x) = v, a CONSTANT: the vertex does not move with x
            S = mkMap(zeros(2), V0, vd);
            faces(end+1) = mkFace('vertex', i, j, R, fN, 2*vd^2*dd, S); %#ok<AGROW>
        end
    end
end

% ------------------------------------------------------------------------------------------------
function f = mkFace(kind, i, j, R, fN, fD, S)
    if fD < 0, fN = -fN;  fD = -fD; end
    [fN, fD] = ratQ.canon(fN, fD);
    f = struct('who', struct('kind', kind, 'facet', i, 'facet2', j), ...
               'R', R, 'fN', fN, 'fD', fD, 'S', S);
end

% ------------------------------------------------------------------------------------------------
function S = mkMap(Mn, tn, d)
% The maximiser as s*(x) = (Mn*x + tn)/d, with the denominator kept POSITIVE so a caller can read
% a sign off the numerator without carrying one.
    if d < 0, Mn = -Mn;  tn = -tn;  d = -d; end
    S = struct('Mn', ratQ.chk(Mn, 'maximiser map'), ...
               'tn', ratQ.chk(tn(:), 'maximiser translation'), ...
               'd',  ratQ.chk(d, 'maximiser denominator'));
end
