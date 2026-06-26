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
%   * TODO: the general piecewise / nonconvex pipeline (Steps 1-3) and the QuaPar output class.
%           Port the per-piece symbolic envelope+conjugate from codeOld/cPLQ
%           (conjugateExpr.m, plq_1p.m, functionNDomain.m, region.m, quadQuad.m).
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

    % ---- General case: [JOGO] 3-step pipeline (TODO) ------------------------------------
    error('PLQ:conjCPLQ:notImplemented', ...
        ['General ''cplq'' conjugate (piecewise / nonconvex domains) is not implemented yet. ' ...
         'It requires the QuaPar (quadratic-on-parabolic) class and the per-piece symbolic ' ...
         'envelope + conjugate ported from codeOld/cPLQ. See DESIGN.md II.5.1.']);
end
