function rays = pieceRecessionRays(piece)
% Extreme rays of a maxQuaPar PIECE's constraint recession cone -- the directions in which the region
% QuaPar.eval will actually admit (every bounding conic, sign-oriented, <= tol) extends to infinity.
%
% This is a CLOSED-FORM decision, not a numerical probe: a direction d recesses the region iff it
% recesses every one of its bounding edges, and that is decided exactly.
%   * a straight edge V(i)->V(j) (the piece is CCW, interior on the LEFT): recessive iff
%     cross(V(j)-V(i), d) >= 0;
%   * a ray edge (dirIn arrives at V(1) from infinity, so its CCW traversal is -dirIn; dirOut leaves
%     V(nv), so +dirOut): same left-of-travel test on that direction;
%   * the one parabolic edge (curveEc, normalised by facePoly/assignSide to evalConic > 0 on the
%     piece's OWN interior, so the constraint is evalConic >= 0): recessive iff the leading
%     coefficient A(d) = [dx dy] Q [dx dy]' of curveEc's quadratic part is >= 0. Q of a parabola is
%     semidefinite, so A(d) keeps one sign except along its null direction, where the LINEAR term
%     breaks the tie (handled below).
%
% Extreme rays of the cone occur where a single constraint is tight, so the candidate directions are
% the edge directions and the arc's null direction; each is kept iff it recesses every constraint.
% All sign tests are exact (sym over the pieces' rational coefficients) -- no sampling, no tolerance
% on a probed value.
%
% Returns a k-by-2 array of unit representative directions. k == 0 means the recession cone is {0},
% i.e. the piece is a COMPACT face. For a piece WITH rays the ray edges are constraints too, so the
% result is the recession cone the representation actually encodes -- the caller compares it against
% the intended cone(dirIn,dirOut). See FARFIELD_FIX_PLAN.md.
    V   = sym(piece.V);
    nv  = size(V,1);
    unb = ~isempty(piece.dirIn);
    % curveAfter ~= 0 does NOT mean "this edge is curved": boundedPiece tags every piece it builds
    % with the closing edge's index, including the straight-splitting-curve case where curveEc is
    % all zeros (maxQuaPar's pieceIsCurved says why the tag has to stay). Treating such an edge as
    % an arc skips it in the edge-direction loops below AND replaces its half-plane with a conic
    % that is identically zero, so the piece comes out with one fewer constraint than it has edges
    % -- unbounded where it is not.
    ca  = piece.curveAfter;
    if ca ~= 0 && (isempty(piece.curveEc) || all(piece.curveEc == 0)), ca = 0; end

    % CCW straight-edge direction vectors (the arc edge is excluded; it enters via Q below).
    E = {};
    if ~unb
        for i = 1:nv
            if i == ca, continue; end
            j = mod(i,nv)+1;
            E{end+1} = V(j,:) - V(i,:); %#ok<AGROW>
        end
    else
        din = sym(piece.dirIn); dout = sym(piece.dirOut);
        E{end+1} = -din;                                  % incoming ray, walked toward V(1)
        for i = 1:nv-1
            if i == ca, continue; end                     % arc on segment V(i)->V(i+1)
            E{end+1} = V(i+1,:) - V(i,:); %#ok<AGROW>
        end
        E{end+1} = dout;                                  % outgoing ray, walked away from V(nv)
    end

    hasArc = (ca ~= 0);
    Q = sym(zeros(2)); cc = sym(zeros(1,6));
    if hasArc
        cc = sym(piece.curveEc);
        Q  = [cc(1), cc(2)/2; cc(2)/2, cc(3)];            % A(d) = d*Q*d'
        % THE ARC'S CHORD IS A CONSTRAINT TOO, and leaving it out is what made a perfectly bounded
        % arc-piece look non-compact. A curved edge is a bounded ARC; its conic is not, and on the
        % parabola's concave side the kept side of the conic wraps around past the arc's own
        % endpoints. The chord cuts that off. It is not an EDGE of the piece, so it is derived --
        % and only when every other vertex of the piece lies on one side of it, which is exactly
        % when the piece is contained in that half-plane and the constraint is redundant for it.
        % QuaPar.chordCuts applies the same rule to an assembled face.
        %
        % WHICH SIDE, AND WHEN. Not decidable from the other vertices, which is what this routine
        % used to do: a LENS (arc plus its chord as a real edge) has both of its vertices ON the
        % chord and no others, so they say nothing at all and the side came out of whichever
        % branch was written first. That is the same mistake QuaPar.chordCuts records in
        % DECISIONS.md, where reading the side off the vertices killed two green tests.
        %
        % Decide it from the conic instead, which is exact and needs no interior point. Along the
        % chord X0 + t*ch the conic restricts to a quadratic q(t) with q(0) = q(1) = 0, so
        % q(t) = A*t*(t-1) with A = ch*Q*ch'. Two consequences, and they settle both questions:
        %   * A <= 0 means q >= 0 on the open chord, i.e. the chord's interior is INSIDE the kept
        %     side {evalConic >= 0}. The piece then straddles the chord line and no chord may be
        %     emitted -- this is exactly the lens, and also every face on the parabola's CONVEX
        %     side, where the conic is already the constraint that closes the region.
        %   * A > 0 means the chord's interior is excluded, so the piece meets the chord line only
        %     at X0 and X1 and lies on ONE side of it: the side the arc's own interior points are
        %     on. Along the parabola's axis direction n (the null direction of Q, where the conic
        %     is LINEAR) from the chord midpoint M, the arc is reached at s* = -conic(M)/(grad.n),
        %     and its side of the chord is the sign of s*.cross(ch,n) -- no point construction and
        %     no tolerance.
        % The vertex test below is then kept only as a VETO, which is what makes the constraint
        % provably redundant for the piece and so unable to shrink it.
        X0 = V(ca,:); X1 = V(mod(ca,nv)+1,:);
        ch = X1 - X0;
        Achord = ch*Q*ch.';
        chDir = sym([]);
        if logical(Achord > 0)
            % THE AXIS DIRECTION MUST COME FROM AN EIGENVECTOR, NOT FROM THE DISCRIMINANT.
            % arcNullDirs solves d*Q*d' = 0 exactly and returns NOTHING when b^2-4ac comes out
            % negative -- which is what a floating-point parabola's Q does about half the time, since
            % it is only semidefinite up to rounding (measured: -2.78e-17 on the first arc of the
            % arcVsArcRefusesAnUnboundedTwoArcSplit fixture). The chord was then silently never
            % emitted, the piece's constraint region stayed a slab open at BOTH ends, and
            % reccConeViolation reported it receding in both senses. parabolaArcFrame has always
            % taken the smallest-magnitude eigenvector instead, with a tolerance; do the same.
            % Only the SIGN of sideArc is used below and its magnitude is bounded away from zero
            % (|cross(ch,n0)| ~ 1.3 and |sStar| ~ 0.04 on that fixture), so deriving n0 numerically
            % does not weaken the exact sign tests that follow.
            Qd = double(Q);
            [Vq0, Dq0] = eig(Qd);
            [~, i00] = min(abs(diag(Dq0)));
            n0d = Vq0(:,i00)';
            nulls = {sym(n0d/norm(n0d))};
            for q0 = nulls
                n0 = q0{1};
                M  = (X0 + X1)/2;
                gM = [2*cc(1)*M(1) + cc(2)*M(2) + cc(4), ...
                      cc(2)*M(1) + 2*cc(3)*M(2) + cc(5)];
                den = gM(1)*n0(1) + gM(2)*n0(2);
                if logical(den == 0), continue, end
                sStar   = -(cc(1)*M(1)^2 + cc(2)*M(1)*M(2) + cc(3)*M(2)^2 ...
                            + cc(4)*M(1) + cc(5)*M(2) + cc(6)) / den;
                sideArc = sStar * (ch(1)*n0(2) - ch(2)*n0(1));
                if logical(sideArc > 0),     chDir = ch;  break
                elseif logical(sideArc < 0), chDir = -ch; break
                end
            end
        end
        if ~isempty(chDir)
            % VETO: every other vertex and ray direction must lie on that same side, so that the
            % half-plane contains the whole piece and the constraint cannot shrink it.
            side = sym([]);
            for i = 1:nv
                if i == ca || i == mod(ca,nv)+1, continue, end
                side(end+1) = chDir(1)*(V(i,2)-X0(2)) - chDir(2)*(V(i,1)-X0(1)); %#ok<AGROW>
            end
            if unb
                for d0 = {sym(piece.dirIn), sym(piece.dirOut)}
                    side(end+1) = chDir(1)*d0{1}(2) - chDir(2)*d0{1}(1); %#ok<AGROW>
                end
            end
            if isempty(side) || all(logical(side >= 0))
                E{end+1} = chDir;                          % interior on the left of the chord
            end
        end
    end

    % Candidate extreme directions: every straight edge direction (both signs) and the arc null dir.
    cand = {};
    for k = 1:numel(E), cand{end+1} = E{k}; cand{end+1} = -E{k}; end %#ok<AGROW>
    if hasArc, cand = [cand, arcNullDirs(Q)]; end

    rays = zeros(0,2);
    for k = 1:numel(cand)
        d = cand{k};
        if all(logical(d == 0)), continue; end
        if recessesAll(d, E, hasArc, Q, cc, V, ca, unb)
            rays = appendDir(rays, d);
        end
    end
end

function tf = recessesAll(d, E, hasArc, Q, cc, V, ca, unb)
% EXACT ARITHMETIC ON INEXACT DATA IS NOT AN EXACT DECISION, and this routine spent a long time
% pretending otherwise. Its coefficients arrive as DOUBLES and are then lifted with sym(), so the
% comparisons below are exact -- on numbers whose last bits are rounding noise. Two of them decide
% the answer, and both mathematically compare against ZERO:
%
%   * A(d) = d*Q*d' is exactly 0 along a parabola's null direction, which for a half-strip piece is
%     the very direction being tested. Measured on the G4 fixture, the SAME direction stored twice
%     (dirIn and dirOut differ in the last bit) gave A = +4.4e-20 and A = -2.5e-17: one passed the
%     `A < 0` test outright, the other was rejected, and because +4.4e-20 is not `== 0` the
%     linear-term TIE-BREAK below -- the branch written precisely for this case -- never ran.
%   * cross(E,d) is exactly 0 when d is along a ray edge, which again is this piece. Measured
%     -5.6e-17, and `< 0` rejected the direction that every actual constraint admits.
%
% The result was a half-strip reporting recession along the NEGATIVE of its own declared ray --
% the whole decision made by noise of order 1e-17. Compare against a tolerance scaled by the
% operands instead, so "zero" means zero at the precision the data actually carries. That widens
% the TIE case; it does not weaken it, because the tie is then settled by the linear term, which
% is bounded away from zero here (measured |grad.d| = 8.1e-02).
%
% The header's claim of exactness survives where it is true: the tie-break, the chord side and the
% veto all compare quantities that are genuinely bounded away from zero.
    tf = false;
    dd = double(d);
    for m = 1:numel(E)
        Em = double(E{m});
        tolC = 1e-12 * max(1, norm(Em)) * max(1, norm(dd));
        if logical(E{m}(1)*d(2) - E{m}(2)*d(1) < -tolC), return; end   % cross(E,d) < 0
    end
    if hasArc
        A = d*Q*d.';
        tolA = 1e-12 * max(1, norm(double(Q))) * max(1, norm(dd))^2;
        if logical(A < -tolA), return; end                        % leading term < 0 -> not recessive
        if logical(abs(A) <= tolA)
            % Tie: along the arc's null direction evalConic is linear in t; recessive iff its slope
            % from an interior point is >= 0. Interior point = a vertex not on the arc (or the arc
            % chord midpoint), stepped slightly inward is unnecessary -- the SLOPE is B = grad.d and
            % is independent of which interior point, so evaluate the gradient of evalConic at any
            % arc endpoint dotted with d.
            if unb, i0 = ca; else, i0 = ca; end
            p0 = V(i0,:);                                          % an arc endpoint
            g  = [2*cc(1)*p0(1) + cc(2)*p0(2) + cc(4), ...
                  cc(2)*p0(1) + 2*cc(3)*p0(2) + cc(5)];            % grad evalConic at p0
            if logical(g(1)*d(1) + g(2)*d(2) < 0), return; end     % slope < 0 -> not recessive
        end
    end
    tf = true;
end

function dirs = arcNullDirs(Q)
% Directions d with d*Q*d' = 0. For a parabola Q is rank-1, giving one (rational) null direction.
    a = Q(1,1); b = 2*Q(1,2); c = Q(2,2);
    dirs = {};
    if logical(a ~= 0)
        disc = b^2 - 4*a*c;
        if logical(disc == 0)
            dirs{end+1} = [-b/(2*a), sym(1)];
        elseif logical(disc > 0)
            s = sqrt(disc);
            dirs{end+1} = [(-b+s)/(2*a), sym(1)];
            dirs{end+1} = [(-b-s)/(2*a), sym(1)];
        end
    elseif logical(c ~= 0)
        dirs{end+1} = [sym(1), sym(0)];
        dirs{end+1} = [-c, b];
    else
        dirs{end+1} = [sym(1), sym(0)];
        dirs{end+1} = [sym(0), sym(1)];
    end
end

function rays = appendDir(rays, d)
    d = double(d);
    n = norm(d);
    if n < 1e-15, return; end
    d = d/n;
    for r = 1:size(rays,1)
        if norm(rays(r,:) - d) < 1e-9, return; end
    end
    rays(end+1,:) = d; %#ok<AGROW>
end
