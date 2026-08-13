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
    ca  = piece.curveAfter;

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
    tf = false;
    for m = 1:numel(E)
        if logical(E{m}(1)*d(2) - E{m}(2)*d(1) < 0), return; end   % cross(E,d) < 0 -> not recessive
    end
    if hasArc
        A = d*Q*d.';
        if logical(A < 0), return; end                            % leading term < 0 -> not recessive
        if logical(A == 0)
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
