function [tf, why] = isConvexQ(obj)
% isConvexQ  Is a QuaPol CONVEX. Decided EXACTLY, by integer sign tests only.
%
% objective: the multi-piece convexity test that `biconjQ`'s `co f = f` short-circuit needs.
%   Per-piece PSD is free and already checked there; what was missing is the condition ACROSS a
%   shared edge, and that is the whole content of this file.
%
% [input]  obj : QuaPol, operable (degree <= 2) and EXACT (fN/fD/VN/VD present).
% [output] tf  : logical.
%          why : '' when tf is true, otherwise a short string naming the FIRST condition that
%                failed and where -- so a caller (and a test) learns which of the four it was
%                rather than only that the answer is false.
%
% ------------------------------------------------------------------------------------------------
% THE THEOREM. Bauschke, Lucet and Phan (ESAIM COCV 22 (2016) 728-742, Theorem 5.1) reduce
% convexity of a piecewise-defined function to the convex subdifferential being nonempty ON THE
% EDGES -- local, but still infinitely many points. Singh and Lucet (SIAM J. Optim. 31 (2021)
% 972-990, "Linear-time convexity test for low-order piecewise polynomials", Theorem 3.7) turn
% that into finitely many tests for PLQ. In their notation, with the facet e = P_i cap P_j
% contained in {x : a'x = beta}, P_i in {a'x >= beta} and P_j in {a'x <= beta}, f is convex iff
%
%   (a) every Q_k is positive semidefinite,
%   (b) f is continuous,
%   (c) at the vertices x_k of e there are delta_k >= 0 with (Q_i-Q_j)x_k + l_i-l_j = delta_k a,
%   (d) and delta_i <= delta_j along every unbounded direction (x_j - x_i) of e.
%
% Both references are in reference/ (BAUSCHKE-16.pdf, SINGH-21.pdf). This follows the PAPER's
% Theorem 3.7. It does NOT follow the thesis Algorithm 3 that QuaPar.isEdgeConvex cites -- that
% thesis is not in reference/ and has not been read here, so it is not cited.
%
% WHY THE DOT PRODUCT AND NOT A DIVISION. Condition (c) is an equation between vectors, and the
% obvious reading of it -- divide (Q_i - Q_j)x_k + l_i - l_j componentwise by a and check the two
% quotients agree -- is what QuaPar.isEdgeConvex does, which is why that routine has to divide by
% a zero component, filter the resulting NaN with rmmissing, and compare against sqrt(eps). None
% of that survives here. Write d = f_i - f_j. Continuity (b) says d vanishes on the whole edge
% LINE, so the derivative of d ALONG the edge is zero and grad d is automatically parallel to a;
% the only content left in (c) is the SIGN, and
%
%       delta(x) = <grad d(x), a> / |a|^2 ,
%
% so delta(x) >= 0 is exactly <grad d(x), a> >= 0 -- one integer dot product, no division at all.
% |a|^2 > 0 is a positive scale, so it cannot change a sign; and every other scaling below
% (clearing vertex denominators, the factor 2 out of the weighted basis) is positive for the same
% reason. Nothing here needs ratQ's degree-4 sign kernel: grad d is affine, hence delta is affine
% along the edge, and an affine function on a segment is nonnegative iff it is at both endpoints.
%
% ORIENTATION IS LOAD-BEARING and is fixed by the mesh, not chosen: F(j,1) is the face on the LEFT
% of edge j (QuaPol.orderEdges), and a = (-(v2y - v1y), v2x - v1x) is the edge direction rotated a
% quarter turn counterclockwise, so a points into the LEFT face. Hence i = F(j,1), j = F(j,2) in
% the theorem's labelling and d = f_{F(j,1)} - f_{F(j,2)}. Getting this backwards would report
% every convex f as concave and every concave f as convex, so it is pinned by the |x| / -|x| pair
% in isConvexQTest.
%
% THE DOMAIN. f = min_k (f_k + iota_{P_k}) can only be convex if dom f is, and that is tested here
% EXACTLY rather than read from obj.dom.isConvex, which orderEdges computes in doubles. A
% polyhedral domain is convex iff it is the intersection of the closed half-planes of its own
% boundary edges, and it is contained in that intersection iff every vertex and every recession
% direction is: if some x satisfied all the half-planes but lay outside, the segment from an
% interior point to x would leave through one boundary edge and stay on that edge's negative side.
% O(ne * nv) rather than the paper's O(ne), and deliberately: nv is small here, the predicate is a
% 2x2 integer determinant, and correctness outranks speed (umbrella CLAUDE.md 7).
%
% NOT COVERED, by assertion rather than by omission: a domain of dimension < 2. `biconjQ` routes
% needles, segments and edge chains to thinEnvelope BEFORE it reaches any convexity question, so
% there is no caller; convexity of a quadratic along a lone segment is a univariate test and would
% be a different routine.

    if ~isa(obj, 'QuaPol')
        error('PLQ:isConvexQ:input', ...
            'isConvexQ takes a QuaPol (quadratic on polyhedral); got a %s.', class(obj));
    end
    obj.assertOperable();
    obj.assertExact();

    why = '';
    nf = size(obj.fN, 1);

    % ---- (a) every piece's Hessian is PSD -- free, and necessary -------------------------------
    % First because it is the cheapest and because the rest assumes it: Theorem 3.7(ii) takes all
    % Q_k PSD as its hypothesis, so a delta test on top of an indefinite piece would be answering
    % a question the theorem does not ask.
    for k = 1:nf
        H = ratQ.hessQ(obj.fN(k,:), obj.fD(k));
        if ~ratQ.isPSD2(H)
            tf = false;
            why = sprintf('face %d has a Hessian that is not positive semidefinite', k);
            return
        end
    end

    % ---- the whole plane, one face: PSD is the whole story --------------------------------------
    if obj.nv == 0
        tf = true;
        return
    end
    if obj.dom.dim < 2
        error('PLQ:isConvexQ:thinDomain', ...
            ['isConvexQ is for a 2-dimensional domain; this one has dimension %d. See the ' ...
             'header: biconjQ routes thin domains to thinEnvelope before asking.'], obj.dom.dim);
    end

    % ---- (a') the DOMAIN is convex ---------------------------------------------------------------
    [ok, jBad] = domainIsConvex(obj);
    if ~ok
        tf = false;
        why = sprintf('the domain is not convex: it leaves the half-plane of boundary edge %d', jBad);
        return
    end

    % ---- (b), (c), (d) one shared edge at a time -------------------------------------------------
    for j = 1:size(obj.F, 1)
        iL = obj.F(j,1);  iR = obj.F(j,2);
        if iL == 0 || iR == 0
            % A boundary edge: the other side is +infinity, whose only convexity content is that
            % the domain is convex, and that has just been established.
            continue
        end
        if iL == iR
            % The same quadratic on both sides -- d is identically zero, so continuity holds and
            % delta = 0. Nothing to test.
            continue
        end

        [dn, dd] = ratQ.sub(obj.fN(iL,:), obj.fD(iL), obj.fN(iR,:), obj.fD(iR));

        [P1, D1] = vertexOf(obj, obj.E(j,1));
        [P2, D2] = vertexOf(obj, obj.E(j,2));

        % a and the edge, both cleared of denominators by the POSITIVE factor D1*D2.
        a = ratQ.chk([-(P2(2)*D1 - P1(2)*D2), P2(1)*D1 - P1(1)*D2], 'edge normal');
        if all(a == 0)
            error('PLQ:isConvexQ:degenerateEdge', ...
                'edge %d has both endpoints at the same vertex, so it has no normal.', j);
        end

        % ---- (b) continuity: d vanishes on the edge LINE. d restricted to the line is a
        % univariate quadratic, so three distinct points decide it (Singh-Lucet Proposition 3.1;
        % Remark 3.2 counts the three parameters). The midpoint is the third.
        Pm = ratQ.chk([P1(1)*D2 + P2(1)*D1, P1(2)*D2 + P2(2)*D1], 'midpoint numerator');
        Dm = ratQ.chk(2*D1*D2, 'midpoint denominator');
        if evalNum(dn, P1, D1) ~= 0 || evalNum(dn, P2, D2) ~= 0 || evalNum(dn, Pm, Dm) ~= 0
            tf = false;
            why = sprintf('f is discontinuous across edge %d', j);
            return
        end

        % ---- (c) delta >= 0 at both endpoints, up to the positive factors D_k and |a|^2 --------
        G = ratQ.gradQ(dn, dd);               % grad d (x,y) = G*[x;y;1] / dd, with dd > 0
        s1 = deltaNum(G, a, P1, D1);
        s2 = deltaNum(G, a, P2, D2);
        if s1 < 0 || s2 < 0
            tf = false;
            why = sprintf(['the subdifferential is empty on edge %d: the gradient jump points ' ...
                           'the wrong way at an endpoint'], j);
            return
        end

        % ---- (d) along a ray, delta must not DECREASE ------------------------------------------
        % E(j,3) == 0 marks a ray, whose base is E(j,1) and whose direction is E(j,2) - E(j,1);
        % Theorem 3.7 asks delta_i <= delta_j along it. s1 and s2 carry different positive scales
        % (D1 and D2), so they are compared cross-multiplied rather than directly.
        if obj.E(j,3) == 0
            if ratQ.chk(s1*D2, 'delta comparison') > ratQ.chk(s2*D1, 'delta comparison')
                tf = false;
                why = sprintf(['the subdifferential empties out along the unbounded edge %d: ' ...
                               'delta decreases'], j);
                return
            end
        end
    end

    tf = true;
end

% ------------------------------------------------------------------------------------------------
function [P, D] = vertexOf(obj, i)
% The i-th vertex as the exact rational P/D, with P integer 1x2 and D > 0.
    P = obj.VN(i,:);  D = obj.VD(i);
    if D <= 0
        error('PLQ:isConvexQ:badVertex', 'vertex %d has a nonpositive denominator.', i);
    end
end

% ------------------------------------------------------------------------------------------------
function v = evalNum(n, P, D)
% objective: a number with the SIGN of the quadratic n (RatCon's weighted basis, any positive
%   denominator) at the rational point P/D -- and in particular zero exactly when the value is.
%
% The weighted basis is [x^3 x^2y xy^2 y^3 x^2 xy y^2 x y 1] with weights
% [1/6 1/2 1/2 1/6 1/2 1 1/2 1 1 1] (QuaPol.evalPoly), so for a quadratic the value is
% (n5 x^2)/2 + n6 xy + (n7 y^2)/2 + n8 x + n9 y + n10. Multiplying by the POSITIVE 2*D^2 clears
% both the halves and the point's denominator and leaves an integer.
    v = ratQ.chk(n(5)*P(1)^2 + 2*n(6)*P(1)*P(2) + n(7)*P(2)^2 ...
               + 2*n(8)*P(1)*D + 2*n(9)*P(2)*D + 2*n(10)*D^2, 'jump value');
end

% ------------------------------------------------------------------------------------------------
function s = deltaNum(G, a, P, D)
% objective: a number with the SIGN of delta at the rational point P/D, and comparable BETWEEN the
%   two endpoints of one edge once divided by that point's D.
%
% delta = <grad d, a> / |a|^2 with grad d (x) = G*[x;y;1]/dd. Multiplying by the positive
% D*dd*|a|^2 leaves <G*[P; D], a>, an integer.
    g = ratQ.chk(G * [P(1); P(2); D], 'gradient jump');
    s = ratQ.chk(g(1)*a(1) + g(2)*a(2), 'delta numerator');
end

% ------------------------------------------------------------------------------------------------
function [ok, jBad] = domainIsConvex(obj)
% objective: is dom f convex, exactly. See the header for why this is the intersection-of-its-own-
%   boundary-half-planes test, and why the vertices plus the recession directions decide it.
    ok = true;  jBad = 0;

    % Every recession direction of the domain, cleared of denominators by a positive factor.
    rays = zeros(0,2);
    for k = 1:size(obj.E,1)
        if obj.E(k,3) == 0
            [Q1, C1] = vertexOf(obj, obj.E(k,1));
            [Q2, C2] = vertexOf(obj, obj.E(k,2));
            rays(end+1,:) = ratQ.chk(Q2*C1 - Q1*C2, 'ray direction'); %#ok<AGROW>
        end
    end

    for j = 1:size(obj.F,1)
        iL = obj.F(j,1);  iR = obj.F(j,2);
        if iL ~= 0 && iR ~= 0, continue, end          % an interior edge is not a domain constraint
        if iL == 0 && iR == 0
            error('PLQ:isConvexQ:danglingEdge', ...
                'edge %d has +infinity on both sides, so the mesh is not a subdivision.', j);
        end
        [P1, D1] = vertexOf(obj, obj.E(j,1));
        [P2, D2] = vertexOf(obj, obj.E(j,2));
        % l(x) = a'x + c, scaled by the positive D1*D2, is > 0 strictly to the LEFT of the edge.
        a = ratQ.chk([-(P2(2)*D1 - P1(2)*D2), P2(1)*D1 - P1(1)*D2], 'boundary normal');
        c = ratQ.chk(P1(1)*P2(2) - P2(1)*P1(2), 'boundary offset');
        sgn = 1;
        if iL == 0, sgn = -1; end                     % the finite side is the RIGHT one
        % The scales: l at the vertex Q/C is (a'Q + c*C)/(C*D1*D2), and D1*D2 is already inside a
        % and c, so a'Q + c*C carries the sign.
        for i = 1:size(obj.VN,1)
            [Q, C] = vertexOf(obj, i);
            if sgn * ratQ.chk(a(1)*Q(1) + a(2)*Q(2) + c*C, 'domain half-plane') < 0
                ok = false;  jBad = j;  return
            end
        end
        for r = 1:size(rays,1)
            if sgn * ratQ.chk(a(1)*rays(r,1) + a(2)*rays(r,2), 'domain recession') < 0
                ok = false;  jBad = j;  return
            end
        end
    end
end
