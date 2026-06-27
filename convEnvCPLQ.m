function r = convEnvCPLQ(obj)
% convEnvCPLQ  Step 1 of the 'cplq' conjugate pipeline: convex envelope of ONE quadratic piece.
%
% objective: compute conv(q + I_P), the convex envelope of a single quadratic q over its
%   (convex) polyhedral domain P, returned as a RatPol. This is Step 1 of the algorithm in
%     [COAP] Karmarkar & Lucet, Comput. Optim. Appl. 94 (2026) 747-780 (Sect. 3, Appendix A);
%     [JOGO] Karmarkar & Lucet, J. Glob. Optim. 94 (2026) 3-34.
%
% [input]  obj : QuaPoly with a SINGLE piece (nf==1: either the full domain nv==0, or one face).
%                Must be operable (quadratic numerator; cubic is rejected).
% [output] r   : RatPol = conv(q + I_P).
%
% Classification by the (constant) Hessian Q of q = 1/2 x'Q x + L'x + kappa:
%   * Q positive semidefinite  -> q is convex on the convex set P, so conv = q itself
%                                 (valid on the full domain or any single face).
%   * Q negative semidefinite over a TRIANGLE -> conv is the affine function interpolating the
%                                 three vertex values (the "no convex edges" case, Appendix A.2).
%   * Q indefinite, BILINEAR part beta*x*y (beta>0, no x^2/y^2 terms) over a triangle:
%                                 classify the 0/1/2 strictly convex edges (slope > 0) and apply
%                                 [COAP] Appendix A.2 (0 -> affine), A.3 eq.16 (1 -> rational
%                                 quadratic/linear), A.4 (2 -> quadratic, harmonic-mean form),
%                                 then add the linear part. Verified against the Appendix A
%                                 examples (e.g. conv(xy) over conv{(1,1),(0,0),(2,0)} =
%                                 2 y^2 / (y - x + 2)).
%   * Q indefinite, GENERAL (x^2/y^2 terms, beta<=0, or 3 convex edges) -> NOT IMPLEMENTED YET:
%                                 needs the eigen-rotation to x*y (and back) and/or the
%                                 3-convex-edge triangle split (Appendix A.5). See DESIGN.md
%                                 II.5.1 and codeOld/cPLQ (plq_1p.convexEnvelope).
%
% Multi-piece inputs (nf>1) are handled by the pipeline that loops over pieces (later); here a
% single piece is required.

    obj.assertOperable();                 % quadratic numerator (cubic rejected)
    if obj.nf ~= 1
        error('convEnvCPLQ:notImplemented', ...
            ['convEnvCPLQ handles a single quadratic piece (nf==1); got nf=%d. ' ...
             'Multi-piece convex envelopes are assembled by the pipeline (not implemented yet).'], obj.nf);
    end

    [L,Q,C] = QuaPoly.matrixForm(obj.f(1,:)); %#ok<ASGLU>
    ev   = eig(Q);
    isPSD = all(ev >= -sqrt(eps));
    isNSD = all(ev <=  sqrt(eps));
    f6 = obj.f(1,5:10);                   % quadratic numerator in 6-coeff [x^2 xy y^2 x y const]

    % ---- Convex case: envelope = q itself (on the full domain or any single convex face) ----
    if isPSD
        if obj.nv == 0
            r = RatPol(f6);                       % full-domain quadratic
        else
            r = RatPol(obj.V, obj.E, f6, obj.F);  % same quadratic on the same face
        end
        return
    end

    % ---- Non-convex on the full domain: convex envelope is -inf (not a finite RatPol) -------
    if obj.nv == 0
        error('convEnvCPLQ:notImplemented', ...
            'Convex envelope of a non-convex full-domain quadratic is not finite; not implemented.');
    end

    % ---- Non-convex requires a bounded triangle for the closed-form cases below -------------
    isTriangle = (obj.nv == 3) && (obj.ne == 3) && obj.isDomBounded;
    if ~isTriangle
        error('convEnvCPLQ:notImplemented', ...
            ['Convex envelope of a non-convex quadratic is currently implemented only over a ' ...
             'bounded triangle (got nv=%d, ne=%d, bounded=%d).'], obj.nv, obj.ne, obj.isDomBounded);
    end

    % ---- Concave over a triangle: affine interpolation of the three vertex values -----------
    if isNSD
        z   = QuaPoly.evalPoly(obj.f(1,:), obj.V);   % 3x1 values q(v_i) at the triangle vertices
        abc = [obj.V, ones(3,1)] \ z;                % [alpha; beta; gamma]: z = alpha x + beta y + gamma
        flin = [0 0 0 abc(1) abc(2) abc(3)];         % linear function in 6-coeff format
        r = RatPol(obj.V, obj.E, flin, obj.F);
        return
    end

    % ---- Indefinite over a triangle ---------------------------------------------------------
    % General indefinite q has the form (quadratic part) = x^2/y^2/xy terms with one positive
    % and one negative eigenvalue. Reducing a general indefinite quadratic to the bilinear form
    % x*y requires an eigen-rotation (and tracking the change of variables back) -- not done yet.
    % Here we handle the already-bilinear case: quadratic part = beta*x*y (no x^2, y^2 terms)
    % with beta>0, plus an arbitrary linear part. This covers [COAP] Appendix A's worked
    % examples and produces real rational/quadratic/linear envelope pieces.
    if abs(f6(1)) > sqrt(eps) || abs(f6(3)) > sqrt(eps)
        error('convEnvCPLQ:notImplemented', ...
            ['Convex envelope of a general indefinite quadratic (with x^2 or y^2 terms) needs a ' ...
             'rotation to the bilinear form x*y, which is not implemented yet. Only the pure ' ...
             'bilinear case (quadratic part = beta*x*y) is supported. See DESIGN.md II.5.1.']);
    end
    beta = f6(2);                          % actual coefficient of x*y
    if beta <= sqrt(eps)
        error('convEnvCPLQ:notImplemented', ...
            ['Only a positive bilinear coefficient (beta>0) is supported; beta<=0 needs the ' ...
             'general reduction (reflection/rotation), not implemented yet. See DESIGN.md II.5.1.']);
    end
    [num6, den3] = bilinearTriangleEnvelope(obj.V, beta, f6(4:6));
    r = RatPol(obj.V, obj.E, num6, obj.F, den3);
end

% ===================== local functions ======================================================
function [num6, den3] = bilinearTriangleEnvelope(V, beta, lin)
% conv(beta*x*y + d*x + e*y + f0) over the triangle with vertices V (3x2), beta>0.
% Returns the envelope as numerator coefficients num6 (STORED weighted basis [x^2 xy y^2 x y 1])
% and a linear denominator den3 = [g h k] (den3 = [0 0 1] when the envelope is a polynomial).
% Implements [COAP] Appendix A.2 (0 convex edges -> linear), A.3 eq.16 (1 -> rational),
% A.4 (2 -> quadratic). 3 convex edges (needs splitting) is not implemented.
    tol = sqrt(eps);
    d = lin(1); e = lin(2); f0 = lin(3);
    % classify each of the 3 edges: convex (for beta>0) iff its slope m > 0
    ed = [1 2; 2 3; 3 1];
    cm = []; cq = []; copp = [];          % slope, intercept, opposite-vertex index of convex edges
    for t = 1:3
        i = ed(t,1); j = ed(t,2); opp = 6 - i - j;
        vi = V(i,:); vj = V(j,:); dx = vj(1)-vi(1); dy = vj(2)-vi(2);
        if abs(dx) < tol, continue; end    % vertical edge: x*y is linear along it (not convex)
        m = dy/dx;
        if m > tol
            cm(end+1) = m; cq(end+1) = vi(2) - m*vi(1); copp(end+1) = opp; %#ok<AGROW>
        end
    end
    nCE = numel(cm);
    switch nCE
        case 0   % Appendix A.2: affine interpolation of q at the 3 vertices
            qv  = beta*(V(:,1).*V(:,2)) + d*V(:,1) + e*V(:,2) + f0;
            abc = [V, ones(3,1)] \ qv;     % q ~ abc(1) x + abc(2) y + abc(3)
            numPlain = [0 0 0 abc(1) abc(2) abc(3)]; den3 = [0 0 1];
        case 1   % Appendix A.3 eq.16: rational (quadratic over linear)
            m = cm(1); q = cq(1); x1 = V(copp(1),1); y1 = V(copp(1),2);
            a = -m*y1; b = q; c = x1;
            dd = -q*y1 + m*x1*y1; ee = -q*x1 - x1*y1; ff = q*x1*y1;
            g = -m; h = 1; k = -y1 + m*x1;
            numB = [a b c dd ee ff];       % conv(x*y) numerator (plain), over denB
            denB = [g h k];
            % conv(beta*xy + lin) = (beta*numB + lin*denB) / denB
            linDen = [d*g, d*h+e*g, e*h, d*k+f0*g, e*k+f0*h, f0*k];
            numPlain = beta*numB + linDen; den3 = denB;
        case 2   % Appendix A.4: quadratic (harmonic-mean form)
            q2 = twoEdgeQuadPlain(cm(1), cq(1), cm(2), cq(2), V);  % conv(x*y), plain coeffs
            numPlain = beta*q2 + [0 0 0 d e f0]; den3 = [0 0 1];
        otherwise % nCE == 3
            error('convEnvCPLQ:notImplemented', ...
                ['Bilinear triangle with 3 convex edges must be split into two triangles ' ...
                 '(Appendix A.5); not implemented yet.']);
    end
    num6 = [2*numPlain(1), numPlain(2), 2*numPlain(3), numPlain(4), numPlain(5), numPlain(6)];
end

function q = twoEdgeQuadPlain(mh, qh, mw, qw, V)
% conv(x*y) over a triangle with two convex edges y=mh x+qh and y=mw x+qw (Appendix A.4),
% as plain coefficients [a b c d e f] of a x^2+b xy+c y^2+d x+e y+f. Of the two (+/-) solutions
% we pick the one that touches f=x*y along a convex edge.
    for s = [1 -1]
        cand = buildTwoEdge(mh, qh, mw, qw, s);
        xt = mean(V(:,1)); yt = mh*xt + qh;             % a point on convex edge h
        if abs(evalPlain(cand, xt, yt) - xt*yt) < 1e-7*(1+abs(xt*yt))
            q = cand; return;
        end
    end
    q = buildTwoEdge(mh, qh, mw, qw, 1);                 % fallback (matches [COAP] example)
end

function c = buildTwoEdge(mh, qh, mw, qw, s)
    r = sqrt(mh*mw); denom = mh + mw + s*2*r;
    c = [ mh*mw/denom, s*2*r/denom, 1/denom, ...
          (mh*qw + mw*qh)/denom, -(qh+qw)/denom, qh*qw/denom ];
end

function v = evalPlain(c, x, y)
    v = c(1)*x.^2 + c(2)*x.*y + c(3)*y.^2 + c(4)*x + c(5)*y + c(6);
end
