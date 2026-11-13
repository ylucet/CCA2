function lb = conjEdgeLowerBound(q, S)
% conjEdgeLowerBound  A CLOSED-FORM lower bound on f*(s), from the boundary of dom f alone.
%
% objective: a check that works on EVERY conjugate route, including the single-piece one. The fold
%   cross-check in conjCPLQ verifies an assembled result against the identity
%   f* = max_k (q_k + I_P_k)*, which needs at least two pieces; a single triangle has no such
%   identity and, until this, no check at all. G4 -- conj returning a MINORANT on some indefinite
%   triangles -- is exactly that hole.
%
% [input]  q : the QuaPol whose conjugate is being checked
%          S : m x 2 dual points
% [output] lb : m x 1 lower bound, +Inf where the sup along an unbounded edge diverges
%
% THE BOUND, and why it is valid. For ANY subset A of dom f,
%       f*(s) = sup_{x in dom f} <s,x> - f(x)  >=  sup_{x in A} <s,x> - f(x).
% Taking A to be one EDGE of the domain makes the right-hand side a one-dimensional quadratic
% maximisation with a closed form: along v_i + t(v_j - v_i) the objective is A2 t^2 + A1 t + A0,
% so the max over t in [0,1] is the better endpoint, or the stationary point when the quadratic is
% concave and that point is interior. Taking the best over all edges of all faces gives `lb`.
%
% MEASURED, and this is why the EDGE and not the VERTEX version. The vertex-only bound -- the same
% argument with A a single point -- is violated by NOTHING: 0 of 24 random cases, including the one
% that is wrong by 2.7e-2, because that sup is not attained at a vertex. The edge bound fires on
% exactly that case and on nothing else (every other case at ~1e-15). DECISIONS.md 2026-08-25.
%
% WHAT IT CANNOT SEE. A sup attained strictly INSIDE the domain. For a concave or indefinite piece
% the maximiser of <s,x> - q(x) is on the boundary, so the bound is tight there; for a convex piece
% it can be interior, and then the bound is genuinely slack. That is the route conjConvexPolygon
% handles in closed form and which the random sweep found exact in every case -- but the limit is
% real and is stated rather than glossed.

    S = double(S);
    lb = -inf(size(S,1), 1);
    for i = 1:q.nf
        f10 = q.f(i,:);
        ej = find(any(q.F == i, 2));
        for t = 1:numel(ej)
            j = ej(t);
            a = q.V(q.E(j,1), :);
            b = q.V(q.E(j,2), :);
            isSeg = q.E(j,3) == 1;
            for k = 1:size(S,1)
                lb(k) = max(lb(k), maxOnEdge(S(k,:), a, b, f10, isSeg));
            end
        end
    end
end

function m = maxOnEdge(s, a, b, f10, isSeg)
% max of <s,x> - q(x) over the segment a->b, or over the RAY from a in the direction b - a.
%
% A ray's second stored point is a direction marker, not a point of the domain (RatCon.m's `E`), so
% the direction is b - a either way and only the parameter RANGE differs: [0,1] for a segment,
% [0,inf) for a ray. Reading a ray's marker as a domain point would make the bound wrong -- too
% large on the marker side, which would turn a correct conjugate into a reported defect.
    d = b - a;
    Q = [f10(5), f10(6); f10(6), f10(7)];
    L = [f10(8); f10(9)];
    A2 = -0.5 * (d * Q * d.');
    A1 = s * d.' - (a * Q * d.' + L.' * d.');
    A0 = s * a.' - (0.5*(a*Q*a.') + L.'*a.' + f10(10));
    m = A0;
    if isSeg
        m = max(m, A0 + A1 + A2);
        if A2 < 0
            tt = -A1/(2*A2);
            if tt > 0 && tt < 1, m = max(m, A2*tt^2 + A1*tt + A0); end
        end
        return
    end
    % A RAY: t in [0, inf).
    if A2 > 0 || (A2 == 0 && A1 > 0)
        m = inf; return                       % the sup diverges: f* is +Inf at this s
    end
    if A2 < 0
        tt = -A1/(2*A2);
        if tt > 0, m = max(m, A2*tt^2 + A1*tt + A0); end
    elseif A1 > 0
        m = inf;
    end
end
