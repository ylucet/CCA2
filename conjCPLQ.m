function g = conjCPLQ(obj, idx)
% conjCPLQ  Fenchel conjugate via the symbolic per-piece ('cplq') engine.
%
% objective: Compute the Legendre-Fenchel conjugate f*(s) = sup_x <s,x> - f(x) of a 2D
%   piecewise linear-quadratic function, following the three-step algorithm of
%     [JOGO] Karmarkar & Lucet, J. Glob. Optim. 94 (2026) 3-34, and
%     [COAP] Karmarkar & Lucet, Comput. Optim. Appl. 94 (2026) 747-780:
%       Step 1  convex envelope of each quadratic piece  -> rational (quad/linear) on polyhedral
%       Step 2  conjugate of each rational piece (Lagrange multipliers) -> quadratic on parabolic
%       Step 3  maximum of the conjugates                -> f* (quadratic on parabolic)
%
% [input]  obj : QuaPoly, operable (degree<=2)
%          idx : (optional) variable index 1 or 2 for the PARTIAL conjugate (not yet implemented)
% [output] g   : the conjugate. Returned as a QuaPoly when it is itself quadratic-on-polyhedral
%                (e.g. a full-domain strictly convex quadratic); otherwise it is a QuaPar
%                (quadratic on a parabolic subdivision) -- see "not implemented" notes below.
%
% STATUS (incremental implementation -- see DESIGN.md II.5.1):
%   * IMPLEMENTED (exact): full-domain strictly convex quadratic -> full-domain quadratic.
%   * IMPLEMENTED: a single bounded-triangle piece (any classification -- affine, convex
%     (PD or rank-1 PSD), concave, or genuinely indefinite with 0 or 1 convex edge) -> QuaPar.
%     Step 2 (conjPieceCPLQ) conjugates the ORIGINAL piece directly whenever it can (since
%     f*=(conv f)*, this needs no Step 1 envelope at all); when it can't (a concave piece, or
%     an indefinite piece with 2 convex edges), Step 1 (convEnvCPLQ) is computed first and its
%     (single-face) envelope is conjugated instead -- see conjSingleTriangle below.
%   * TODO: Step 3 (maximum of several conjugates), needed whenever a domain is genuinely
%     covered by MORE THAN ONE piece: a multi-face input (nf>1), a single non-triangular face,
%     or a single triangle whose indefinite quadratic has 3 convex edges (Step 1 splits it into
%     2 sub-triangles, each convexified -- see convEnvCPLQ -- with no single-piece shortcut).
%
% NOTE on arithmetic: the design target is exact symbolic + rational arithmetic
%   ([COAP]/[JOGO]); this first version uses double precision for the closed-form quadratic
%   case. Upgrading conjCPLQ to symbolic/rational coefficients is a follow-up.

    if nargin >= 2 && ~isempty(idx)
        error('PLQ:conjCPLQ:partialNotImplemented', ...
            'Partial conjugate (''cplq'' engine) is not implemented yet.');
    end
    obj.assertOperable();   % degree<=2 (cubic rejected; cubic is for isConvex only)

    % ---- Case A: full-domain quadratic (no vertices, single face) -----------------------
    % f(x) = 1/2 x'Q x + L'x + kappa over all of R^2.
    if obj.nv == 0 && obj.nf == 1
        [L, Q, C] = QuaPoly.matrixForm(obj.f);   %#ok<ASGLU> (C is empty for quadratics)
        kappa = obj.f(end);
        ev = eig(Q);
        if all(ev > sqrt(eps))
            % Strictly convex: the conjugate is again a full-domain quadratic
            %   f*(s) = 1/2 (s-L)' inv(Q) (s-L) - kappa.
            M    = inv(Q);                 %#ok<MINV> small 2x2, explicit inverse is fine/clear
            grad = -M * L;                 % linear part of f*
            d    = 0.5 * (L' * M * L) - kappa;
            % Store in the 6-coefficient quadratic format [x^2 xy y^2 x y const]:
            %   matrixForm reads Q=[c5 c6; c6 c7], L=[c8;c9], const=c10.
            f6 = [M(1,1), M(1,2), M(2,2), grad(1), grad(2), d];
            g  = QuaPoly(f6);
            return
        end
        % Not strictly convex: the conjugate is an indicator / parabolic object (a QuaPar),
        % not a full-domain quadratic. Handled by the general pipeline (not yet implemented).
        error('PLQ:conjCPLQ:notImplemented', ...
            ['Conjugate of a non-strictly-convex full-domain quadratic yields an indicator ' ...
             'or parabolic result (QuaPar), which is not implemented yet. See DESIGN.md II.5.1.']);
    end

    % ---- Case B: a single bounded-triangle piece -----------------------------------------
    if obj.nf == 1 && obj.nv == 3 && obj.ne == 3 && obj.isDomBounded
        g = conjSingleTriangle(obj);
        return
    end

    % ---- General case: [JOGO] Step 3 (max of conjugates) is still TODO ------------------
    error('PLQ:conjCPLQ:notImplemented', ...
        ['General ''cplq'' conjugate over more than one piece is not implemented yet: it needs ' ...
         'Step 3 (pointwise maximum of the per-piece conjugates), which is not implemented. ' ...
         'See DESIGN.md II.5.1.']);
end

% ================================================================================================
function g = conjSingleTriangle(obj)
% Step 2 (+ Step 1 when needed) for a single bounded-triangle piece, following the "sidestep"
% noted in conjPieceCPLQ's header: since f*=(conv f)*, try to conjugate the ORIGINAL piece
% directly first (conjPieceCPLQ handles affine, PD/rank-1-PSD convex, and indefinite with 0 or 1
% convex edge this way). That fails only for a concave piece or an indefinite piece with >=2
% convex edges (both raise 'conjPieceCPLQ:notImplemented' -- the only error conjPieceCPLQ can
% raise here, since obj already satisfies its shape/type/degree preconditions); in that case
% compute Step 1's convex envelope (convEnvCPLQ) and conjugate that instead: a concave piece's
% envelope is a single affine piece, and an indefinite piece with exactly 2 convex edges envelopes
% to a single rank-1 PSD quadratic (COAP Appendix A.4) -- both directly conjugable. An indefinite
% piece with 3 convex edges envelopes to TWO sub-triangle pieces (convEnvCPLQ splits it): that
% case genuinely needs Step 3 (max of the two sub-conjugates), not yet implemented.
    try
        g = conjPieceCPLQ(obj);
        return
    catch ME
        if ~strcmp(ME.identifier, 'conjPieceCPLQ:notImplemented')
            rethrow(ME);
        end
    end
    env = convEnvCPLQ(obj);
    if env.nf == 1
        g = conjPieceCPLQ(env);
        return
    end
    error('PLQ:conjCPLQ:notImplemented', ...
        ['Indefinite triangle with 3 convex edges: Step 1''s envelope splits into %d sub-triangle ' ...
         'pieces, whose conjugates must be combined by Step 3 (max of conjugates), which is not ' ...
         'implemented yet. See DESIGN.md II.5.1.'], env.nf);
end
