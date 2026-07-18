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
% [output] g   : the conjugate. QuaPoly when it is itself quadratic-on-polyhedral (e.g. a
%                full-domain strictly convex quadratic); QuaPar for a single bounded-triangle
%                piece (Case B); for a genuinely multi-face or non-triangular bounded domain
%                (Case C, below), g is instead a cPLQ `functionNDomain` array (NOT a QuaPoly/
%                QuaPar) -- evaluate it with `evalFunctionNDomain(g, s)`, not `g.eval(s)`. Case C
%                does not (yet) support further composition (biconj/infConv/moreau/...): those
%                call `.conj()` again on `g`, which only QuaPoly/QuaPar provide.
%
% STATUS (incremental implementation -- see DESIGN.md II.5.1):
%   * IMPLEMENTED (exact): full-domain strictly convex quadratic -> full-domain quadratic.
%   * IMPLEMENTED: a single bounded-triangle piece (any classification -- affine, convex
%     (PD or rank-1 PSD), concave, genuinely indefinite with 0 or 1 convex edge, or genuinely
%     indefinite with 3 convex edges) -> QuaPar. Step 2 (conjPieceCPLQ) conjugates the ORIGINAL
%     piece directly whenever it can (since f*=(conv f)*, this needs no Step 1 envelope at all);
%     when it can't (a concave piece, or an indefinite piece with 2 or 3 convex edges), Step 1
%     (convEnvCPLQ) is computed first. A 2-convex-edge envelope is a single rank-1-PSD face,
%     conjugated directly; a 3-convex-edge envelope splits into TWO 2-convex-edge sub-triangles
%     (COAP Appendix A.5), each conjugated separately and combined via Step 3 (maxQuaPar,
%     pointwise maximum) -- see conjSingleTriangle/conjMaxOfSubTriangles below.
%   * IMPLEMENTED (Phase 1, this session -- see DESIGN.md II.5.1 / .claude/SESSION_HANDOFF.md):
%     Step 3 for a bounded domain genuinely covered by more than one ORIGINAL piece (nf>1) or a
%     single non-triangular face, via the integrated cPLQ symbolic pipeline (quaPolyToPlq.m ->
%     plq.triangulate -> plq.maximum), NOT the numeric conjPieceCPLQ/maxQuaPar path (which still
%     cannot combine curved-edge QuaPars from independent triangles -- see maxQuaPar.m's own TODO).
%     This is Case C below. NOTE the return-type caveat above: g is a functionNDomain array here,
%     not a QuaPoly/QuaPar, so composition (biconj/infConv/moreau/...) is not supported for this
%     case yet -- only direct conjugate + evalFunctionNDomain. `cPLQ`'s own biconjugateF has known
%     open bugs (region.getNormalConeVertexQ) not needed for Case C itself.
%
% NOTE on arithmetic: the design target is exact symbolic + rational arithmetic
%   ([COAP]/[JOGO]); Case A/B use double precision for the closed-form quadratic case, while
%   Case C uses cPLQ's own exact symbolic + rational arithmetic directly. Upgrading Case A/B to
%   symbolic/rational coefficients (or Case C to closed-form numerics, Phase 2) is a follow-up.

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

    % ---- Case C: general bounded domain (nf>1 and/or a non-triangular face) -------------
    % [JOGO] Step 3 (max of conjugates) via the integrated cPLQ symbolic pipeline (Phase 1;
    % DESIGN.md II.5.1 / .claude/SESSION_HANDOFF.md) -- the case Case B's own numeric path
    % (conjPieceCPLQ + maxQuaPar) cannot do, since maxQuaPar refuses curved-edge QuaPar inputs
    % from independent triangles. g is a cPLQ functionNDomain array here, not a QuaPoly/QuaPar --
    % see this function's own header for the return-type caveat.
    if obj.isDomBounded
        p = quaPolyToPlq(obj);
        p = p.triangulate;
        p = p.maximum;
        g = p.maxConjugate;
        return
    end

    error('PLQ:conjCPLQ:notImplemented', ...
        ['General ''cplq'' conjugate over an unbounded multi-face domain is not implemented yet. ' ...
         'See DESIGN.md II.5.1.']);
end

% ================================================================================================
function g = conjSingleTriangle(obj)
% Step 2 (+ Step 1/3 when needed) for a single bounded-triangle piece, following the "sidestep"
% noted in conjPieceCPLQ's header: since f*=(conv f)*, try to conjugate the ORIGINAL piece
% directly first (conjPieceCPLQ handles affine, PD/rank-1-PSD convex, and indefinite with 0 or 1
% convex edge this way). That fails only for a concave piece or an indefinite piece with >=2
% convex edges (both raise 'conjPieceCPLQ:notImplemented' -- the only error conjPieceCPLQ can
% raise here, since obj already satisfies its shape/type/degree preconditions); in that case
% compute Step 1's convex envelope (convEnvCPLQ) and conjugate that instead: a concave piece's
% envelope is a single affine piece, and an indefinite piece with exactly 2 convex edges envelopes
% to a single rank-1 PSD quadratic (COAP Appendix A.4) -- both directly conjugable. An indefinite
% piece with 3 convex edges envelopes to TWO sub-triangle pieces (COAP Appendix A.5, convEnvCPLQ's
% splitThreeConvex): Step 2 conjugates each sub-triangle separately, and Step 3 (maxQuaPar)
% combines them into the true conjugate of the original piece -- see conjMaxOfSubTriangles.
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
    g = conjMaxOfSubTriangles(env);
end

% ================================================================================================
function g = conjMaxOfSubTriangles(env)
% Step 3: conjugate every triangle face of a multi-face convex-envelope RatPol (Step 1's output)
% independently (Step 2, conjPieceCPLQ), then combine them via repeated pairwise maxQuaPar. Only
% reached from conjSingleTriangle's 3-convex-edge case, where env is always the 2-sub-triangle
% split of COAP Appendix A.5 -- each sub-triangle is provably 2-convex-edge (a rank-1 PSD
% quadratic, no curved domain edges), so both conjPieceCPLQ (Step 2) and maxQuaPar's polyhedral-
% domain requirement (Step 3) are satisfied exactly; see maxQuaParTest.m for the same construction
% validated against ground truth.
    g = toQuaPar(conjPieceCPLQ(extractTriangleFace(env, 1)));
    for k = 2:env.nf
        gk = toQuaPar(conjPieceCPLQ(extractTriangleFace(env, k)));
        g = maxQuaPar(g, gk);
    end
end

function p = extractTriangleFace(r, k)
% Extract face k of a multi-face RatPol as a standalone single-triangle RatPol -- the shape
% conjPieceCPLQ requires (nf=1, nv=3, ne=3, bounded). conjPieceCPLQ reads only p.V/p.f/p.den (it
% fixes V's CCW order itself via its own triSignedArea check), so the vertex order returned by
% orderEdges' walk (documented as clockwise) needs no correction here.
    iVs = faceVertexIndices(r, k);
    if numel(iVs) ~= 3
        error('PLQ:conjCPLQ:notImplemented', ...
            ['conjCPLQ: Step 3 currently supports only triangular envelope pieces (face %d has ' ...
             '%d vertices).'], k, numel(iVs));
    end
    V3 = r.V(iVs, :);
    E3 = [1 2 1; 2 3 1; 3 1 1];
    F3 = [1 0; 1 0; 1 0];
    p = RatPol(V3, E3, r.f(k,:), F3, r.den(k,:));
end

function iVs = faceVertexIndices(obj, k)
% Vertex indices around face k, in the order of its ordered edge list obj.P{k} (same convention
% as convEnvCPLQ.m's own file-local helper of the same name, duplicated here since it is
% file-local there).
    face = obj.P{k}; iVs = zeros(1, numel(face));
    for i = 1:numel(face)
        j = face(i);
        if j > 0, iVs(i) = obj.E(j,1); else, iVs(i) = obj.E(-j,2); end
    end
end
