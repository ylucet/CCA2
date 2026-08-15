function [expr, why] = convEnvUnbounded(r, q, vars)
% convEnvUnbounded  The convex envelope of a quadratic over one piece of an unbounded face,
%   in the case where that envelope is AFFINE -- Step 1's nCE==0 branch, generalized off the
%   bounded triangle and, unlike cPLQ's own closed forms, actually a function of q.
%
% [input]  r     : region, one piece produced by fanUnboundedFace (triangle, half-strip or
%                  wedge) or a bounded triangle.
%          q     : sym, the quadratic 1/2 x'Qx + L'x + c in vars.
%          vars  : 1x2 sym, [x y].
% [output] expr  : sym, the affine envelope.
%          why   : char, which shape decided it (for diagnostics and tests).
%
% WHAT THIS COMPUTES, AND THE ARGUMENT THAT IT IS THE ENVELOPE. co q = q** is the supremum of
% the affine minorants of q on P. That sup is in general PIECEWISE affine, not affine -- so
% "here is an affine minorant that touches q" is NOT on its own a proof that it is the
% envelope. What is a proof: parametrize the affine minorants, minimize the gap at an
% arbitrary point of P, and observe that the minimizer comes out INDEPENDENT of that point.
% One affine function is then best everywhere at once, so the sup is affine and equals it.
% Each of the three shapes fanUnboundedFace emits supplies exactly the 3 conditions an affine
% function on the plane has degrees of freedom:
%
%   WEDGE v + cone(d1,d2), with d1'Qd1 = d2'Qd2 = 0. Write a minorant by its gap parameters
%     e0 = q(v) - l(v) and bi = <grad q(v) - a, di>. At z = v + alpha*d1 + beta*d2 the gap is
%         e0 + alpha*b1 + beta*b2 + alpha*beta*d1'Qd2 .
%     Letting alpha or beta run to infinity forces b1 >= 0 and b2 >= 0, and e0 >= 0 is the gap
%     at v itself. Since alpha,beta >= 0, the gap at ANY z is minimized by e0 = b1 = b2 = 0 --
%     the same minorant for every z. Hence
%         l(z) = q(v) + <grad q(v), z - v>,
%     the tangent plane at the apex. It is a minorant exactly when d1'Qd2 >= 0; if d1'Qd2 < 0
%     then (d1+d2)'Q(d1+d2) < 0 is a recession direction of negative curvature, q is unbounded
%     below, and the envelope is -inf rather than affine.
%
%   HALF-STRIP conv{v1,v2} + cone(d), with d'Qd = 0 and w'Qw <= 0 for w = v2 - v1 (q concave
%     along the base edge; were it convex the envelope would bend along that edge). With
%     e0 = q(v1) - l(v1), A = <grad q(v1) - a, w> and B = <grad q(v1) - a, d>, the gap at
%     z = v1 + s*w + t*d, s in [0,1], t >= 0, is
%         e0 + s*A + t*B + s^2*w'Qw/2 + s*t*w'Qd .
%     t -> infinity forces B + s*w'Qd >= 0 for every s in [0,1], i.e. B >= max(0, -w'Qd);
%     s = 1, t = 0 forces A >= -e0 - w'Qw/2; and e0 >= 0. Minimizing e0 + s0*A + t0*B at a
%     point (s0,t0) drives A to its bound, leaving e0*(1 - s0) - s0*w'Qw/2, which for s0 in
%     [0,1] is minimized at e0 = 0 -- again independently of (s0,t0). So l(v1) = q(v1),
%     l(v2) = q(v2), and
%         <a,d> = min( <grad q(v1),d>, <grad q(v2),d> ) ,
%     the two agreeing exactly when w'Qd = 0.
%
%   TRIANGLE conv{v1,v2,v3}: three vertex values, l(vi) = q(vi) -- the classical interpolant,
%     valid when q is CONCAVE along every edge (d'Qd <= 0). This is cPLQ's own nCE == 0 branch
%     ([COAP]), with one difference that matters: cPLQ's closed form is a function of the
%     vertex COORDINATES only and silently computes the envelope of x*y whatever q is (pinned
%     by test: it returns the same envelope for x*y, x^2-y^2, (x^2+y^2)/2 and 3xy+7x-2y+5).
%     This one interpolates the actual values of q.
%
% WHAT IT REFUSES, and why refusing beats guessing. A recession direction with d'Qd > 0 means q
% is strictly convex along a ray, so the envelope is curved (cPLQ's nCE == 1/2 shapes) and no
% affine answer is correct -- 'convEnvUnbounded:convexAlongRay'. A face on which q is unbounded
% below has co q = -inf and hence conjugate identically +inf, which is a different kind of
% answer than a function -- 'convEnvUnbounded:minusInfinity'. Both are reported, never
% silently approximated, because both would otherwise surface as a plausible-looking wrong
% envelope far downstream.

    x = vars(1); y = vars(2);
    Q = double(hessian(q, [x y]));
    L = double(subs(gradient(q, [x y]), [x y], [0 0]));
    L = L(:);
    gradAt = @(v) Q*v(:) + L;
    qAt    = @(v) double(subs(q, [x y], double(v(:)')));

    tolQ = 1e-9 * max(1, max(abs(Q(:))));

    [D, kind] = r.recessionRays;
    [nP, px, py] = r.finiteVertices;
    V = uniqueRows(double([px(:), py(:)]));

    % The -inf gate comes first: it is the only outcome that is not a function, and every
    % branch below would otherwise produce a finite-looking affine answer for it.
    if ~strcmp(kind, 'bounded')
        if r.quadUnboundedBelow(Q, L)
            error('convEnvUnbounded:minusInfinity', ...
                ['q is unbounded below on this piece, so its convex envelope is -inf and its ' ...
                 'conjugate is identically +inf. That is not an affine envelope and must not ' ...
                 'be reported as one.']);
        end
    end

    % A CONVEX q needs no Step 1 at all: co(q|P) = q on a convex P. Worth taking before the
    % shape branches, not after, because for a positive definite Q every direction has
    % d'Qd > 0 and every branch below would refuse it as "convex along a ray" -- correctly, in
    % that the envelope is genuinely not affine, but pointlessly, since there is nothing to
    % compute.
    %
    % UNBOUNDED PIECES ONLY, deliberately. On a BOUNDED triangle this routine is called from
    % convexEnvelope1's nCE==0 branch, whose contract is "no convex edge, so the envelope is the
    % affine interpolant" -- and returning a quadratic there changes what Step 2 is handed, which
    % its nCE==0 branch cannot conjugate. Taking the short-circuit for bounded pieces turned
    % conjCPLQTest's own convex Case C from a (wrong, see below) answer into an error. Note what
    % that reveals: nCE tests edge SLOPES, which only detects convex edges for f = x*y, so a
    % convex q on a triangle can be misclassified as nCE==0 and get an interpolant where the
    % true envelope is q itself. That misclassification is pre-existing and recorded in
    % SUPPORT_MATRIX.md section 7; it is not made worse here, and fixing it means teaching Step 2
    % the convex-over-a-triangle case, which is a separate job.
    if ~strcmp(kind, 'bounded') && all(eig((Q + Q')/2) >= -tolQ)
        expr = q;
        why  = 'convex (envelope is q itself)';
        return
    end

    switch kind
    case 'bounded'
        if size(V,1) ~= 3
            error('convEnvUnbounded:notTriangle', ...
                'a bounded piece must be a triangle here; this one has %d vertices.', size(V,1));
        end
        for e = [1 2; 2 3; 3 1]'
            d = V(e(2),:) - V(e(1),:);
            if d * Q * d' > tolQ
                error('convEnvUnbounded:convexAlongRay', ...
                    ['q is strictly convex along a triangle edge (a cPLQ "convex edge"), so ' ...
                     'the envelope is not affine; this is the nCE >= 1 case.']);
            end
        end
        M = [V, ones(3,1)];
        if abs(det(M)) < 1e-12
            error('convEnvUnbounded:degenerate', 'the three vertices are collinear.');
        end
        abc = M \ [qAt(V(1,:)); qAt(V(2,:)); qAt(V(3,:))];
        expr = abc(1)*x + abc(2)*y + abc(3);
        why  = 'triangle';

    case 'wedge'
        if size(V,1) ~= 1
            error('convEnvUnbounded:notWedge', ...
                'a wedge piece must carry exactly one finite vertex; this one has %d.', size(V,1));
        end
        v = V(1,:);
        nCvx = 0; iCvx = 0;
        for i = 1:size(D,1)
            if D(i,:) * Q * D(i,:)' > tolQ, nCvx = nCvx + 1; iCvx = i; end
        end
        if nCvx == 1 && size(D,1) == 2
            % ONE ray convex, the other flat: the envelope is CURVED, and it is q with its cross
            % term deleted. Derivation, in the same style as the affine cases above.
            %
            % Write z = v + alpha*d1 + beta*d2 with d1 the FLAT ray (A11 = d1'Q d1 = 0) and d2 the
            % convex one (A22 > 0), alpha,beta >= 0, and let gi = <grad q(v), di>, Aij = di'Q dj.
            % Then
            %       q(z) = q(v) + alpha*g1 + beta*g2 + alpha*beta*A12 + beta^2*A22/2 .
            % Parametrise the affine minorants by their gap parameters e0 = q(v) - l(v) and
            % bi = <grad q(v) - a, di>; the gap at z is
            %       e0 + alpha*b1 + beta*b2 + alpha*beta*A12 + beta^2*A22/2 .
            % alpha -> infinity at fixed beta forces b1 + beta*A12 >= 0 for every beta >= 0, i.e.
            % b1 >= 0 AND A12 >= 0 -- if A12 < 0 then d1 + t*d2 is a recession direction of
            % negative curvature, q is unbounded below, and the envelope is -inf, which the guard
            % below reports. Taking alpha = 0 and minimising over beta >= 0 forces
            % b2^2 <= 2*A22*e0. Minimising the gap at a target (alpha0,beta0) then gives b1 = 0,
            % b2 = -sqrt(2*A22*e0) and e0 = beta0^2*A22/2, whose value is exactly
            %       alpha0*beta0*A12 .
            % So the envelope is q minus its cross term:
            %       co q(z) = q(v) + alpha*g1 + beta*g2 + beta^2*A22/2 ,
            % which is convex (linear in alpha, convex quadratic in beta, composed with a linear
            % map) and minorises q because alpha*beta*A12 >= 0 on the wedge. Unlike the affine
            % cases the minimiser DEPENDS on the target point (through beta0), which is precisely
            % why the envelope here is not affine.
            %
            % CHECKED against unboundedFaceTest/curvedEnvelopeOverAWedgeIsExact: q = x*y on
            % K = {0 <= y <= x} has v = (0,0), d1 = (1,0), d2 = (1,1), A11 = 0, A22 = 2, A12 = 1,
            % grad q(v) = 0, and beta = y -- giving co q = y^2, which is the envelope that test
            % derives by hand (its best affine minorant at (1,1) is worth 1, and the apex tangent
            % plane 0).
            i2 = iCvx; i1 = 3 - iCvx;
            d1 = D(i1,:); d2 = D(i2,:);
            A12 = d1 * Q * d2';
            if A12 < -tolQ
                error('convEnvUnbounded:minusInfinity', ...
                    ['q recedes to -infinity along d1 + t*d2 (the cross term d1''Q d2 = %g is ' ...
                     'negative while d1 is flat), so the envelope is -inf, not a function.'], A12);
            end
            if abs(det([d1; d2])) < 1e-12
                error('convEnvUnbounded:degenerate', 'the wedge''s two rays are parallel.');
            end
            gv = gradAt(v);
            g1 = gv' * d1(:);  g2 = gv' * d2(:);
            ab = [d1; d2]' \ [x - v(1); y - v(2)];      % [alpha; beta] as affine forms in x,y
            expr = expand(qAt(v) + g1*ab(1) + g2*ab(2) + (d2*Q*d2')*ab(2)^2/2);
            why  = 'wedge (curved: one convex ray)';
            return
        end
        for i = 1:size(D,1)
            d = D(i,:);
            if d * Q * d' > tolQ
                error('convEnvUnbounded:convexAlongRay', ...
                    ['q is strictly convex along a wedge ray, so the envelope is curved, not ' ...
                     'affine.']);
            end
        end
        a = gradAt(v);
        expr = qAt(v) + a(1)*(x - v(1)) + a(2)*(y - v(2));
        why  = 'wedge';

    case 'ray'
        if size(V,1) ~= 2
            error('convEnvUnbounded:notHalfStrip', ...
                'a half-strip piece must carry exactly two finite vertices; this one has %d.', ...
                size(V,1));
        end
        d = D(1,:);
        if d * Q * d' > tolQ
            % CONVEX along the ray. The envelope is curved, and when the base edge and the ray
            % are Q-ORTHOGONAL (w'Q d = 0) it separates: q is then a function of s plus a
            % function of t in the frame z = v1 + s*w + t*d, so
            %       co q = (envelope of the s-part on [0,1]) + (the t-part, already convex).
            % The s-part is concave (w'Q w <= 0, which the guard below requires for the affine
            % case and which is what makes the base edge need convexifying at all), so its
            % envelope on [0,1] is the AFFINE INTERPOLANT between its endpoint values -- and
            % those endpoints are q(v1) and q(v2). Hence
            %       co q(z) = q(v1) + s*(q(v2) - q(v1)) + t*<grad q(v1), d> + t^2*(d'Q d)/2 .
            %
            % CHECKED against
            % unboundedFaceTest/nonconvexQuadraticWithACurvedEnvelopeOverAHalfStripIsExact:
            % q = -x^2 + y^2 on {0 <= x <= 1, y >= 0} has v1 = (0,0), v2 = (1,0), w = (1,0),
            % d = (0,1), w'Q d = 0, q(v1) = 0, q(v2) = -1, <grad q(v1), d> = 0 and d'Q d = 2,
            % giving co q = -x + y^2 -- the envelope that test states.
            %
            % w'Q d ~= 0 is NOT handled: the two directions then interact and the minimiser of the
            % gap moves with the target point in both coordinates. Refused loudly rather than
            % approximated.
            w0 = V(2,:) - V(1,:);
            if abs(w0 * Q * d') > tolQ
                error('convEnvUnbounded:convexAlongRay', ...
                    ['q is strictly convex along the half-strip''s direction AND its base edge ' ...
                     'is not Q-orthogonal to it (w''Q d = %g), so the envelope is curved and ' ...
                     'does not separate.'], w0 * Q * d');
            end
            if w0 * Q * w0' > tolQ
                error('convEnvUnbounded:convexAlongRay', ...
                    ['q is strictly convex along BOTH the half-strip''s base edge and its ' ...
                     'direction; the envelope is q itself only if the whole piece is convex, ' ...
                     'which the short-circuit above already handles.']);
            end
            if abs(det([w0; d])) < 1e-12
                error('convEnvUnbounded:degenerate', ...
                    'the half-strip''s base edge is parallel to its rays.');
            end
            st = [w0; d]' \ [x - V(1,1); y - V(1,2)];   % [s; t] as affine forms in x,y
            expr = expand(qAt(V(1,:)) + (qAt(V(2,:)) - qAt(V(1,:)))*st(1) ...
                          + (gradAt(V(1,:))' * d(:))*st(2) + (d*Q*d')*st(2)^2/2);
            why  = 'half-strip (curved: convex along the ray, separable)';
            return
        end
        v1 = V(1,:); v2 = V(2,:); w = v2 - v1;
        if w * Q * w' > tolQ
            error('convEnvUnbounded:convexAlongRay', ...
                ['q is strictly convex along the half-strip''s base edge, so the envelope bends ' ...
                 'along that edge and is not affine.']);
        end
        if abs(det([w; d])) < 1e-12
            error('convEnvUnbounded:degenerate', ...
                'the half-strip''s base edge is parallel to its rays.');
        end
        slope = min(gradAt(v1)' * d(:), gradAt(v2)' * d(:));
        % Solve <a,w> = q(v2)-q(v1) and <a,d> = slope, then fix the constant at v1.
        a = [w; d] \ [qAt(v2) - qAt(v1); slope];
        c0 = qAt(v1) - a' * v1(:);
        expr = a(1)*x + a(2)*y + c0;
        why  = 'half-strip';

    otherwise
        error('convEnvUnbounded:notPointed', ...
            ['recessionRays reports ''%s'': this piece has no apex, so it is not one of the ' ...
             'three shapes an affine envelope is rigid on.'], kind);
    end

    expr = simplify(expr);
    if nP < 1 && ~strcmp(kind,'bounded')
        error('convEnvUnbounded:noVertex', 'piece has no finite vertex to anchor the envelope.');
    end
end

% ------------------------------------------------------------------------------------------
function V = uniqueRows(V)
    keep = true(size(V,1),1);
    for i = 2:size(V,1)
        for j = 1:i-1
            if keep(j) && norm(V(i,:) - V(j,:)) < 1e-9, keep(i) = false; break, end
        end
    end
    V = V(keep,:);
end
