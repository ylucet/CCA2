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
% [input]  obj : QuaPol, operable (degree<=2)
%          idx : (optional) variable index 1 or 2 for the PARTIAL conjugate (not yet implemented)
% [output] g   : the conjugate. QuaPol when it is itself quadratic-on-polyhedral (e.g. a
%                full-domain strictly convex quadratic); QuaPar for a single bounded-triangle
%                piece (Case B); for a genuinely multi-face or non-triangular bounded domain
%                (Case C, below), g is a QuaParCPLQ -- a thin wrapper around cPLQ's own
%                `functionNDomain` array (NOT a QuaPol/QuaPar: see QuaParCPLQ.m's own header for
%                why) that provides the same operator surface (conj/add/scalarMul/addQuadratic/
%                eval), so biconj/infConv/moreau/proxAverage compose with it exactly as they do
%                with QuaPol/QuaPar -- see QuaParCPLQ.m for the STATUS/limitations of that
%                composition (in particular biconjugateF's known open mergeL/removeTangent exact-
%                tie-point gap -- DESIGN.md II.5.1 -- is inherited by QuaParCPLQ.conj).
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
%     single non-triangular face, via the integrated cPLQ symbolic pipeline (quaPolToPlq.m ->
%     plq.triangulate -> plq.maximum), NOT the numeric conjPieceCPLQ/maxQuaPar path (which still
%     cannot combine curved-edge QuaPars from independent triangles -- see maxQuaPar.m's own TODO).
%     This is Case C below. g is a QuaParCPLQ (see the return-type note above and QuaParCPLQ.m):
%     composition (biconj/infConv/moreau/proxAverage) IS supported, by driving the same per-piece
%     machinery `plq.biconjugateF` already used (functionNDomain.conjugateOfPiecePoly/mergeL/
%     addEq) -- inheriting that recipe's own known open mergeL/removeTangent exact-tie-point gap
%     (DESIGN.md II.5.1) for domains it affects.
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
        [L, Q, C] = QuaPol.matrixForm(obj.f);   %#ok<ASGLU> (C is empty for quadratics)
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
            g  = QuaPol(f6);
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
    % from independent triangles. g is a QuaParCPLQ here, not a QuaPol/QuaPar -- see this
    % function's own header and QuaParCPLQ.m for the return-type rationale.
    % UNBOUNDED domains take the same route. The isDomBounded gate that used to stand here was
    % load-bearing while quaPolToPlq threw a face's ray direction away and plq_1p read region's
    % intmax direction markers as ordinary coordinates; both are fixed (see quaPolToPlq's header
    % and fanUnboundedFace/convEnvUnbounded). What is NOT yet done is the conjugate of a CURVED
    % convex envelope over an unbounded face, and that raises
    % plq_1p:conjugateFunction:unboundedNonAffine from inside this call rather than silently
    % returning a number -- which is why the gate can come out here rather than being replaced
    % by a narrower one that would have to duplicate the same test.
    p = quaPolToPlq(obj);
    p = p.triangulate;
    p = p.maximum;
    % Case C gets the same Step 3 cross-check conjEnvelopeViaCPLQ already applies. It was not
    % wired here before, and that was a real hole: the per-piece conjugates can each be right
    % while the cross-piece maximum drops cells, which returns plausible numbers rather than
    % erroring. Measured on the 4-cone fan with convex faces -- each of the 4 faces produces a
    % correct 4-cell conjugate, the assembled maximum keeps only 4 of the 16, and f*(-0.5,2)
    % comes back 1.125 where the truth is 2.
    assertStep3MatchesPieces(p, obj);
    g = QuaParCPLQ(p.maxConjugate);
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
    try
        if env.nf == 1
            g = conjPieceCPLQ(env);
        else
            g = conjMaxOfSubTriangles(env);
        end
        return
    catch ME
        if ~strcmp(ME.identifier, 'conjPieceCPLQ:notImplemented')
            rethrow(ME);
        end
        % The only way to get here is Step 1 having produced a RATIONAL face: Step 2 has no
        % rational branch (conjPieceCPLQ.m:107). cPLQ's Step 2 does -- fall back to it.
    end
    g = conjEnvelopeViaCPLQ(env);
end

% ================================================================================================
function g = conjEnvelopeViaCPLQ(env)
% Conjugate Step 1's envelope through cPLQ's OWN symbolic Step 2/3 instead of the numeric
% conjPieceCPLQ + maxQuaPar pair, for the envelopes the numeric path cannot take: those with a
% genuinely rational face, which Step 1 emits whenever it has to split an indefinite triangle
% (2- or 3-convex-edge). Returns a QuaParCPLQ, exactly as Case C does -- the same wrapper, for the
% same reason (the V/E/Ec/F/P mesh is not reconstructed from cPLQ's symbolic regions; see
% QuaParCPLQ.m).
%
% This is "start from cPLQ" applied where it belongs. Note carefully WHICH half is reused: cPLQ's
% Step 2/3, not its Step 1. Running cPLQ end to end on the ORIGINAL triangle does not work, and
% that is a statement about cPLQ's Step 1 only:
%   * for a 2-convex-edge triangle cPLQ's own envelope is the single Appendix A.4 formula, which
%     CCA2 established is not always tight (convEnvCPLQTest.bilinearTwoConvexEdgesSplitIsTight);
%     its conjugate is correspondingly wrong, and on conv{(2,1),(0,0),(1,0)} it leaves the paper's
%     own flagged dual point s=(-0.008727,-0.999962) covered by NO region at all (evaluates NaN);
%   * cPLQ's Step 1 has no nCE==3 branch (neither convexEnvelope1 nor conjugateFunction), so fed
%     a 3-convex-edge triangle directly it silently returns an EMPTY envelope.
% Neither is a reason to teach anything a 3-convex-edge case: [COAP] Appendix A.5's split turns
% such a triangle into 2-convex-edge sub-triangles, and CCA2's Step 1 already does that. By the
% time cPLQ's Step 2 is called here, every face it sees has come through that split. So the
% working combination is CCA2 Step 1 + cPLQ Step 2/3, which is what this does.
%
% VALIDATION (2026-07-28): on conv{(2,1),(0,0),(1,0)} with f=xy -- the 2-convex-edge split, whose
% envelope has a rational face -- the result matches sup_{x in T} <s,x> - f(x) to <= 8.9e-16 at
% all 10 sampled dual points, including the flagged point above where cPLQ's own end-to-end path
% returns NaN. See conjCPLQTest.indefiniteTriangleTwoConvexEdgesSplitViaCPLQStep2.
    p = ratPolToPlq(env);
    try
        p = p.maximum;
    catch ME
        where = '';
        if ~isempty(ME.stack)
            where = sprintf(' at %s:%d', ME.stack(1).name, ME.stack(1).line);
        end
        error('PLQ:conjCPLQ:cplqFailed', ...
            ['Step 2 fell back to cPLQ for a rational envelope face, and cPLQ''s symbolic ' ...
             'pipeline failed%s (%s: %s). Step 2 itself (the per-piece conjugate) is known to ' ...
             'complete on these envelopes; what is not yet reliable is Step 3, cPLQ''s ' ...
             'cross-piece maximum (plq.maximumConjugate -> functionNDomain.maximumP -> ' ...
             'region.maximum). See SUPPORT_MATRIX.md section 1.2.'], ...
            where, ME.identifier, ME.message);
    end
    if isempty(p.maxConjugate)
        error('PLQ:conjCPLQ:cplqFailed', ...
            'cPLQ''s Step 2/3 returned an empty conjugate for Step 1''s envelope.');
    end
    assertStep3MatchesPieces(p, env);
    g = QuaParCPLQ(p.maxConjugate);
end

% ================================================================================================
function assertStep3MatchesPieces(p, env)
% Cross-check cPLQ's Step 3 (the assembled cross-piece maximum, p.maxConjugate) against the
% SAME quantity computed the other way: the pointwise max of the per-piece Step 2 conjugates,
% p.pieces(k).maxConjugate. Both are f*, since f* = (conv f)* = max_k (env_k + I_{face k})*, so
% any disagreement is a bug -- and Step 3 is the half that has one.
%
% WHY THIS GATE EXISTS. It was added while Step 3's assembly was known wrong: region.merge
% unioned two regions by deleting their shared facet and intersecting what remained (an
% over-approximation unless the union happens to be convex), and simplifyUnboundedRegion deleted
% any constraint not passing through a finite vertex (not a redundancy test at all for an
% unbounded or curved region). Both let a region claim territory that was never its own, carrying
% the wrong value there -- 57 of 289 points of a 17x17 dual grid on f=xy over conv{(0,0),(3,3),
% (1,2)}. Both are now FIXED, by LP certificates: region.redundantSubset decides redundancy
% exactly, and region.unionIsExact decides union-exactness exactly (see region.m's own
% LP-certificate header). That case now reports 0 of 289 wrong at every round of the fold.
%
% The gate STAYS anyway, and not out of caution: it is an independent cross-check of the whole
% assembly against the same f* computed the other way, it is cheap next to the symbolic pipeline
% it follows, and the failure it guards against is silent by nature -- a wrong partition returns
% plausible numbers rather than erroring. See .claude/SESSION_HANDOFF.md and SUPPORT_MATRIX.md
% section 1.2.
%
% Sampling, not proof: the grid below is a screen, not a certificate. It is sized from the
% primal domain (the dual scale that matters is set by the gradients over it), and it is cheap
% next to the symbolic pipeline it follows.
    R = 2*max(abs(env.V), [], 'all') + 1;
    t = linspace(-R, R, 11);
    for i = 1:numel(t)
        for j = 1:numel(t)
            s = [t(i), t(j)];
            got = evalFunctionNDomain(p.maxConjugate, s);
            % A piece that covers s with no region is +inf THERE, and max(anything, +inf) = +inf.
            % This loop used to skip such a piece and take the max of the finite ones instead,
            % which is only right when every piece's conjugate is finite everywhere -- true when
            % every PRIMAL face is bounded, and false as soon as one is unbounded, where each
            % piece's conjugate is +inf off a proper cone. Skipping then invents a finite ref and
            % the gate fires on a correct answer: measured on f = max(0,x,y) as three wedges,
            % whose conjugate really is the indicator of the simplex {s>=0, s1+s2<=1}, exact at
            % 9 probes, yet reported as a Step 3 disagreement everywhere outside it.
            % ...but "no region covers s" means +inf only for a piece whose PRIMAL face is
            % unbounded. A bounded face has conjugate finite EVERYWHERE, so an uncovered point
            % there is a coverage gap in Step 2's own output, not a genuine +inf, and demanding
            % that the assembled maximum be +inf on it would fire the gate on cases that are
            % right (conjCPLQTest's indefiniteTriangleThreeConvexEdgesUsesStep3 is one). The
            % piece's own domain decides which reading applies.
            refIsInf = false;
            ref = -inf;
            for k = 1:p.nPieces
                v = evalFunctionNDomain(p.pieces(k).maxConjugate, s);
                if isnan(v)
                    if pieceIsUnbounded(p.pieces(k))
                        refIsInf = true;    % genuinely +inf here, so the max is
                    end
                elseif isfinite(v)
                    ref = max(ref, v);
                end
            end
            if refIsInf
                if isfinite(got)
                    error('PLQ:conjCPLQ:cplqFailed', ...
                        ['cPLQ''s Step 3 disagrees with its own Step 2 at s=(%.6f,%.6f): the ' ...
                         'assembled maximum gives the finite value %.6f, but at least one piece ' ...
                         'has f* = +inf there and the max of anything with +inf is +inf.'], ...
                        s(1), s(2), got);
                end
                continue
            end
            if ~isfinite(ref)
                continue        % no piece covers s: nothing to check against
            end
            if ~isfinite(got) || abs(got - ref) > 1e-7 * max(1, abs(ref))
                error('PLQ:conjCPLQ:cplqFailed', ...
                    ['cPLQ''s Step 3 disagrees with its own Step 2 at the dual point ' ...
                     's=(%.6f,%.6f): the assembled maximum gives %.6f but the pointwise max ' ...
                     'of the per-piece conjugates -- the same f*, computed the other way -- ' ...
                     'gives %.6f. Since Step 1 and Step 2 agree with ground truth on the cases ' ...
                     'this covers, suspect the ASSEMBLY: a region claiming territory that was ' ...
                     'never its own. The two historical causes are fixed and certified by LP ' ...
                     '(region.redundantSubset, region.unionIsExact), so start by checking ' ...
                     'whether either certificate returned UNDECIDED and fell back to its ' ...
                     'conservative branch. See SUPPORT_MATRIX.md section 1.2.'], ...
                    s(1), s(2), got, ref);
            end
        end
    end
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

% ================================================================================================
function tf = pieceIsUnbounded(pc)
% Is this plq_1p piece's PRIMAL face unbounded? Decides whether an uncovered dual point means
% f* = +inf (unbounded face) or a coverage gap (bounded face -- f* is finite everywhere there).
    tf = false;
    try
        [~, kind] = pc.d.polygon.recessionRays;
        tf = ~strcmp(kind, 'bounded');
    catch
        tf = false;     % non-affine facet: not a case this check can decide, so assume bounded
    end
end
