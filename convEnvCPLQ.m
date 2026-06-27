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
%   * Q indefinite (bilinear)  -> NOT IMPLEMENTED YET. The envelope is piecewise linear /
%                                 quadratic / rational (quadratic over linear) on a polyhedral
%                                 subdivision: rotate q to the bilinear form x*y, classify the
%                                 0/1/2 strictly convex edges, and apply the Appendix A.2-A.5
%                                 formulae (linear; rational eq.16; quadratic harmonic-mean).
%                                 See DESIGN.md II.5.1 and codeOld/cPLQ (plq_1p.convexEnvelope).
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

    % ---- Indefinite (bilinear) over a triangle: Appendix A.2-A.5 (rational/quadratic) -------
    error('convEnvCPLQ:notImplemented', ...
        ['Convex envelope of an indefinite (bilinear) quadratic over a triangle is not ' ...
         'implemented yet: rotate to x*y, classify the 0/1/2 strictly convex edges, and apply ' ...
         'the Appendix A.2-A.5 formulae. See DESIGN.md II.5.1 and codeOld/cPLQ.']);
end
