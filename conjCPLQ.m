function g = conjCPLQ(obj, idx, route)
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
    % ROUTE. 'auto' (the default, and what every caller but one uses) tries the numeric path first
    % for a bounded domain and falls back to Case C's symbolic pipeline. 'symbolic' skips straight
    % to Case C.
    %
    % The two produce the SAME function in different representations, so choosing between them is
    % a representation decision, not a mathematical one -- and biconjCPLQ has to make it. Since
    % 2026-08-13 the numeric path completes on a bounded multi-face domain and returns a MESHED
    % QuaPar with parabolic edges, which is what the SCIP bridge wants; but f** = (f*)* then hands
    % that curved QuaPar to a second conjugation, and quaPolToPlq requires a polyhedral domain, so
    % biconj died at quaPolToPlq:curvedEdge on every such input. The symbolic form conjugates
    % again through cPLQ's own machinery, which handles a curved domain. So conj keeps returning
    % the mesh and biconj asks for the symbolic form on purpose.
    if nargin < 3 || isempty(route), route = 'auto'; end
    if ~ismember(lower(route), {'auto','symbolic','numeric'})
        error('PLQ:conjCPLQ:route', ...
            'route must be ''auto'', ''symbolic'' or ''numeric''; got ''%s''.', route);
    end
    forceSymbolic = strcmpi(route, 'symbolic');
    % 'numeric' REFUSES to fall back. It exists because "conj is sym-free except as a fallback" is
    % a claim that has to be testable: with 'auto' a numeric refusal is swallowed and the answer
    % still comes back, so a regression that pushes a case onto the symbolic path is invisible in
    % the values. Under 'numeric' the refusal is rethrown with its own identifier, which is what
    % lets a test pin the ROUTE and a measurement (checkConjSymFree) name the reason.
    forceNumeric = strcmpi(route, 'numeric');
    obj.assertOperable();   % degree<=2 (cubic rejected; cubic is for isConvex only)

    % ---- STEP 0: NORMALISE the subdivision ---------------------------------------------------
    % An interior edge whose two sides carry the same quadratic is a line the caller drew, not a
    % property of f, and every case below dispatches on the mesh: Case A needs one face and no
    % vertices, Case B one bounded triangle, Case C triangulates and then pays Step 3 to maximise
    % ACROSS the pieces -- pieces that this normalisation may have removed entirely. Deleting such
    % an edge changes the mesh and not the function. See mergeSameQuadFaces.m, ALGORITHM.md Step 0.
    obj = mergeSameQuadFaces(obj);

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
        % NOT STRICTLY CONVEX. This used to be one lumped `notImplemented`, and that hid the fact
        % that it is THREE different objects, only one of which is even the same KIND of gap. The
        % mathematics is three lines in each case and is written out below; what none of them has
        % is somewhere to be PUT. A QuaPol/QuaPar face is two-dimensional, and while the mesh does
        % carry a `dim < 2` domain (nv=1 needle, nv=2/ne=1 segment-or-ray, both with nf=0), QuaPar
        % `eval` locates a point on such a domain with `belongToEdge`, which is membership of a
        % SEGMENT between two vertices -- so a whole LINE has no representation and a needle
        % evaluates to +inf even at its own point. The route that could hold these is cPLQ's
        % `region`, which takes equality constraints (`functionNDomain.addEq`), i.e. a
        % QuaParCPLQ -- see `TODO.md` and `RETURN_TYPE.md`.
        %
        % Classified rather than implemented, so the next person starts from the answer.
        evs  = sort(ev);
        tolE = sqrt(eps) * max(1, norm(Q));
        if evs(1) < -tolE
            % A negative eigenvalue: f decreases without bound along that eigendirection, so no
            % affine function minorises it, conv f = -inf, and f*(s) = +inf for EVERY s. dom f*
            % is EMPTY -- not a lower-dimensional domain but no domain at all, which the mesh has
            % no encoding for either (nf=0 already means "dim < 2", not "empty").
            error('PLQ:conjCPLQ:conjugateHasEmptyDomain', ...
                ['f has a negative curvature direction (eig(Q) = [%g %g]), so conv f = -inf and ' ...
                 'f* = +inf everywhere. dom f* is EMPTY, which no QuaPol/QuaPar mesh encodes.'], ...
                evs(1), evs(2));
        elseif evs(2) <= tolE
            % Q = 0: f is AFFINE, f(x) = <L,x> + kappa, so
            %       f*(s) = sup_x <s - L, x> - kappa
            % which is +inf unless s = L, and -kappa there. dom f* is a single POINT, and that is
            % a NEEDLE -- nv=1, ne=0, nf=0, which QuaPar's constructor has always anticipated.
            % It used to raise `PLQ:conjCPLQ:conjugateIsAPoint` because QuaPar.eval returned +inf
            % on a needle, including at the needle's own vertex; eval now has that branch, so the
            % answer can simply be returned. See DECISIONS.md 2026-08-25 (overnight, B3).
            g = QuaPar(L(:).', zeros(0,3), [0 0 0 0 0 -kappa], zeros(0,2));
            return
        else
            % Q is PSD of rank 1: along the null direction n of Q, f is affine with slope <L,n>,
            % so the sup is finite only where <s-L,n> = 0. dom f* is that LINE, and on it
            % f*(s) = 1/2 (s-L)' pinv(Q) (s-L) - kappa.
            n = null((Q+Q.')/2);
            n = n(:,1).' / norm(n(:,1));
            error('PLQ:conjCPLQ:conjugateIsALine', ...
                ['Q is positive semidefinite of rank 1 (eig(Q) = [%g %g]), so f* is finite only ' ...
                 'on the LINE <s - L, n> = 0 with n = (%g,%g), where ' ...
                 'f*(s) = 1/2 (s-L)'' pinv(Q) (s-L) - %g. dom f* is a line: QuaPar''s dim<2 ' ...
                 'domain is a SEGMENT or ray between two vertices, not a line.'], ...
                evs(1), evs(2), n(1), n(2), kappa);
        end
    end

    % ---- Case B: a single bounded-triangle piece -----------------------------------------
    if obj.nf == 1 && obj.nv == 3 && obj.ne == 3 && obj.isDomBounded
        if forceSymbolic
            % route='symbolic' for a single triangle: cPLQ's own per-triangle form, which is what
            % `conjSingleTriangle` already falls back to, NOT Case C. Case C does not cover a
            % single triangle -- its header scopes it to nf>1 and/or a non-triangular face, and
            % sending one there raises cplqFailed after 102 s (measured, DECISIONS.md 2026-08-25
            % G12). `conjEnvelopeViaCPLQ` is the symbolic form for exactly this shape, so the
            % route now has somewhere to go and `biconjCPLQ`'s curved-mesh escape stops being a
            % no-op here: it asks for the symbolic form precisely because the second conjugation
            % cannot take a curved mesh, and it used to get the same curved mesh back.
            g = conjEnvelopeViaCPLQ(convEnvCPLQ(obj));
            return
        end
        g = conjSingleTriangle(obj, forceNumeric);
        verifyEdgeBound(g, obj);
        return
    end

    % ---- Case B2: a BOUNDED polygon, taken through Case B's closed-form path -------------
    % A bounded polygon is a union of triangles, and Case B conjugates a triangle in closed
    % form, so f* = max_k (q + I_T_k)* needs no symbolic engine -- the very decomposition
    % Case C's own `triangulate` step performs, only executed by the numeric machinery instead.
    % Nothing was missing mathematically; the polygon was simply never offered to the code that
    % can already do it, and every one went to Case C's symbolic pipeline.
    %
    % Attempted, not assumed. Step 3 (maxQuaPar) accepts at most ONE curved operand, and an
    % INDEFINITE triangle conjugates to a parabolic QuaPar -- so two or more of those refuse and
    % this falls through to Case C exactly as before. Convex, affine and concave triangles
    % conjugate to POLYHEDRAL pieces, so any number of them combine. That is the real boundary
    % of what this branch buys, and it is the arc-vs-arc clipping gap, not a gap here.
    %
    % Measured before wiring, on the four polygons of checkPolygonNumericPath: convex, affine
    % and concave over a unit square, a general quadrilateral, a parallelogram and a pentagon
    % all come back exact -- 4.4e-16 against an EXACT QP reference, not a sampled one. (A
    % sampled reference reported ~3e-5 on the convex cases; refining it 160 -> 320 -> 640 drove
    % that to 1e-6, which is how it was identified as the reference's error, not CCA2's.)
    % The isDomBounded gate came off on 2026-08-24. It was there because everything below this
    % point needed a bounded TRIANGLE; conjConvexPolygon needs neither, so an unbounded domain
    % whose faces carry CONVEX quadratics now has a numeric route. An unbounded face carrying a
    % non-convex one still declines and still falls through to Case C, by name.
    if ~forceSymbolic
        try
            g = conjPolygonalDomain(obj);
            verifyEdgeBound(g, obj);
            return
        catch ME
            if forceNumeric, rethrow(ME); end
            if ~strcmp(ME.identifier, 'conjPieceCPLQ:notImplemented') && ...
               ~strncmp(ME.identifier, 'maxQuaPar:', 10) && ...
               ~strcmp(ME.identifier, 'PLQ:conjCPLQ:notImplemented') && ...
               ~strcmp(ME.identifier, 'convEnvCPLQ:notImplemented')
                rethrow(ME);
            end
            % Step 2 or Step 3 cannot take some triangle of this domain; Case C's symbolic
            % pipeline can, so fall through to it rather than failing.
            recordFallback(ME.identifier);
        end
    end
    recordFallback('unbounded-or-not-attempted');
    if forceNumeric
        error('PLQ:conjCPLQ:numericRouteDeclined', ...
            ['route=''numeric'' was requested but the numeric path does not cover this input ' ...
             '(nf=%d, nv=%d, bounded=%d). The symbolic Case C would have answered it.'], ...
            obj.nf, obj.nv, obj.isDomBounded);
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
function verifyEdgeBound(g, obj)
% objective: check the returned conjugate against a lower bound that needs NO second computation --
%   `f*(s) >= max over the boundary of dom f of [<s,x> - f(x)]`, closed form per edge -- and raise
%   when it is violated. OPT-IN: set the global CCA2_CONJ_VERIFY.
%
% ON BY DEFAULT, and the measurement that justifies that is worth keeping. Set the global
% CCA2_CONJ_VERIFY to 0 to turn it off.
%
%   * over EVERY fast and slow suite with it on: 363 pass, 0 fail -- it flags no correct answer;
%   * over 24 random polygons and quadratics: it fires on exactly ONE, the case that is genuinely
%     wrong by 2.7e-2, and on nothing else (every other case sits at ~1e-15).
%
% So it converts a class of SILENT wrong answer into a named refusal without touching anything that
% works, which is the trade this file makes everywhere else.
%
% IT RAISES RATHER THAN FALLING BACK, and that is measured too: on the known-bad case the SYMBOLIC
% route returns the same wrong value to six digits. The defect is in the shared Step 1/Step 2
% closed form, not in a route, so there is nothing to fall back TO -- and a fallback that silently
% produced the same minorant would be worse than an error.
%
% WHY IT EARNS ITS KEEP. `conjPolygonalDomain`'s fold cross-check needs at least two pieces (it
% verifies against `max_k (q_k + I_P_k)*`). A SINGLE triangle has no such identity, and that is
% exactly where G4 lives -- `conj` of `xy` on some triangles returns a MINORANT. This bound is the
% only check that covers that route.
    global CCA2_CONJ_VERIFY %#ok<GVMIS>
    if ~isempty(CCA2_CONJ_VERIFY) && ~CCA2_CONJ_VERIFY, return, end
    if ~isa(g, 'RatPar') || ~g.isMeshed(), return, end
    S = edgeBoundProbes(g);
    lb = conjEdgeLowerBound(obj, S);
    for i = 1:size(S,1)
        if ~isfinite(lb(i)), continue, end
        got = g.eval(S(i,:));
        if got < lb(i) - 1e-7*(1 + abs(lb(i)))
            error('PLQ:conjCPLQ:belowEdgeBound', ...
                ['conj returned %.10g at (%.6g,%.6g), but the boundary of dom f alone already ' ...
                 'gives %.10g. f*(s) >= <s,x> - f(x) for every x in dom f, so the answer is a ' ...
                 'MINORANT, not the conjugate.'], got, S(i,1), S(i,2), lb(i));
        end
    end
end

function S = edgeBoundProbes(g)
% Directions x magnitudes, plus the result's own vertices -- a minorant shows up on a REGION, so a
% spread finds it, and a cell boundary is where it shows first.
    R = [0.25 1 3 10];
    th = (0:15) * (2*pi/16);
    S = zeros(0,2);
    for r = R, S = [S; r*[cos(th).' sin(th).']]; end %#ok<AGROW>
    S = [S; 0 0];
    if ~isempty(g.V), S = [S; g.V]; end
    S = unique(round(S, 10), 'rows');
end

% ================================================================================================
function recordFallback(why)
% objective: record WHY this call is about to use the symbolic Case C, so the rate and the reasons
%            can be measured rather than guessed.
%
% This exists because "conj is sym-free except as a fallback" is a claim about a RATE, and a rate
% has to be counted. Off unless the caller asks: set the global to an empty cell and it fills with
% one identifier per fallback. Costs one global read on a path that is already about to start
% MuPAD, so it is free where it matters.
    global CCA2_CONJ_FALLBACK %#ok<GVMIS>
    if iscell(CCA2_CONJ_FALLBACK)
        CCA2_CONJ_FALLBACK{end+1} = why;
    end
end

% ================================================================================================
function g = conjSingleTriangle(obj, forceNumeric)
% forceNumeric: rethrow instead of falling through to the symbolic Case C. See the ROUTE note at
% the top of this file -- `route='numeric'` exists so that "conj is sym-free except as a fallback"
% is TESTABLE, and this function ignored it, so Case B took the symbolic fallback under every
% route. Measured 2026-08-25: a `route='numeric'` call on TODO.md G4's triangle spent 420 s in
% cPLQ before failing there, when the numeric route had already refused it in about 3 s. Callers
% that are themselves inside the numeric path (conjEnvelopeViaCPLQ's per-triangle loop) pass
% false, because for them the fallback is the point.
    if nargin < 2, forceNumeric = false; end
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
            g = conjFaceOrOriginal(env, 1, obj);
        else
            g = conjMaxOfSubTriangles(env, obj);
        end
        return
    catch ME
        if forceNumeric, rethrow(ME); end
        if ~strcmp(ME.identifier, 'conjPieceCPLQ:notImplemented') && ...
           ~strncmp(ME.identifier, 'maxQuaPar:', 10)
            rethrow(ME);
        end
        % Step 2 could not take some face even via the original piece (conjPieceCPLQ), or Step 3
        % could not assemble the faces it did produce (maxQuaPar). cPLQ's symbolic Step 2/3 can --
        % fall back to it.
    end
    g = conjEnvelopeViaCPLQ(env);
end

% ================================================================================================
function tf = allFacesAffine(obj)
% Every face carries a function with no quadratic part. The 6-coefficient row is
% [x^2 xy y^2 x y const], so the first three decide it.
    tf = true;
    for i = 1:obj.nf
        q = obj.f(i, 5:7);
        if any(abs(q) > sqrt(eps) * max(1, max(abs(obj.f(i, 5:10)))))
            tf = false; return
        end
    end
end

function faces = affineFaceList(obj)
% Each face as conjAffinePLQ wants it: the same boundary description conjConvexPolygon takes, plus
% the affine function's gradient and constant.
    faces = struct('W', {}, 'dFirst', {}, 'dLast', {}, 'a', {}, 'b', {});
    for i = 1:obj.nf
        [W, dF, dL] = faceBoundary(obj, i);
        f6 = obj.f(i, 5:10);
        faces(end+1) = struct('W', W, 'dFirst', dF, 'dLast', dL, ...
                              'a', [f6(4), f6(5)], 'b', f6(6)); %#ok<AGROW>
    end
end

% ================================================================================================
function g = conjPolygonalDomain(obj)
% Case B2: conjugate every face of the domain by the cheapest closed form that fits it, and
% combine with Step 3's numeric max. Bounded or unbounded.
%
% Sound because a sup over a union is the max of the sups: f* = max_k (q_k + I_{P_k})*. That holds
% for any COVER, which is also why a fan triangulation -- a cover, not a minimal triangulation --
% is enough where one is still used.
%
% TWO ROUTES PER FACE, and the first is the one that matters:
%
%   * a CONVEX (positive definite) face goes to conjConvexPolygon WHOLE. No triangulation, any
%     number of facets, bounded or not, and the answer is a QuaPol -- polyhedral, so nothing
%     curved reaches Step 3 from it. ALGORITHM.md: "a convex piece never needed a triangle;
%     splitting it only forces Step 3 to glue back together what was never broken."
%   * anything else is fan-triangulated and each triangle goes through Case B's own path, exactly
%     as before. That path needs BOUNDED triangles, so a non-convex UNBOUNDED face declines here
%     and the caller routes the whole domain to Case C.
%
% Any per-piece result that is not a mesh is refused rather than worked around: it means that
% piece's own Step 2 fell back to cPLQ's symbolic form, and mixing the two representations inside
% Step 3 is the confusion this branch exists to avoid.
    % ALL-AFFINE INPUT GOES TO conjAffinePLQ, WHOLE (TODO.md G2). Every face affine means f* is a
    % max of SUPPORT functions, each +inf off a cone -- exactly the shape maxQuaPar refuses, which
    % is why `max(0,x,y)` and its family always went to the symbolic Case C. conjAffinePLQ builds
    % the answer directly from the definition and never enters Step 3 at all. It declines by name
    % when dom f* comes out unbounded, which is the case Step 3 can already do, so falling through
    % on that identifier costs nothing and keeps this a strict addition.
    if allFacesAffine(obj)
        try
            g = conjAffinePLQ(affineFaceList(obj));
            return
        catch ME
            if ~strcmp(ME.identifier, 'PLQ:conjAffinePLQ:unboundedDual'), rethrow(ME); end
        end
    end

    E3 = [1 2 1; 2 3 1; 3 1 1]; F3 = [1 0; 1 0; 1 0];
    gs = {};
    for i = 1:obj.nf
        [Lf, Qf, Cf] = QuaPol.matrixForm(obj.f(i,:));
        f6 = obj.f(i, 5:10);
        isConvexFace = isempty(Cf) && all(eig((Qf+Qf.')/2) > sqrt(eps)*max(1, norm(Qf)));
        if isConvexFace
            [W, dF, dL] = faceBoundary(obj, i);
            gs{end+1} = toQuaPar(conjConvexPolygon(W, dF, dL, Qf, Lf, obj.f(i,end))); %#ok<AGROW>
            continue
        end
        tris = faceTrianglesCCW(obj, i);
        for t = 1:numel(tris)
            gt = conjSingleTriangle(QuaPol(tris{t}, E3, f6, F3), false);
            if ~(isa(gt, 'QuaPar') || isa(gt, 'QuaPol') || isa(gt, 'RatPol'))
                error('PLQ:conjCPLQ:notImplemented', ...
                    ['triangle %d of face %d conjugated to a %s (its own Step 2 fell back to ' ...
                     'the symbolic path), which Step 3 numeric max cannot take.'], ...
                    t, i, class(gt));
            end
            gs{end+1} = toQuaPar(gt); %#ok<AGROW>
        end
    end
    if isempty(gs)
        error('PLQ:conjCPLQ:notImplemented', 'no face to conjugate.');
    end
    g = gs{1};
    for k = 2:numel(gs)
        g = maxQuaPar(g, gs{k});
    end
    assertFoldMatchesPieces(g, gs, obj);
end

function assertFoldMatchesPieces(g, gs, obj)
% objective: verify the assembled Step 3 result against the identity it was built from --
%   f* = max_k (q_k + I_{P_k})* -- and DECLINE if they disagree, so the caller falls back.
%
% WHY THIS IS NOT OPTIONAL, measured. Removing the isDomBounded gate let a 4-cone fan through the
% numeric route, and at s = (-2,-3) it returned 2.0 where the definition sup is 4.5: the fold had
% dropped the cell carrying face 4's strip, leaving 4 cells for what needs many more. Every probe
% point of the OTHER orientation of the same fixture was exact, which is precisely the signature of
% a dropped region -- right almost everywhere, wrong on a set, and silent. `conjCPLQ` already
% applies this same cross-check to Case C (`assertStep3MatchesPieces`) for the same reason; the
% numeric route had no equivalent and now has one.
%
% THE ORACLE IS AN IDENTITY, NOT A SAMPLE OF THE TRUTH. `gs` are the per-face conjugates that were
% folded, so their pointwise max IS f* exactly -- no reference implementation, no quadrature, no
% tolerance on the mathematics. Only the SAMPLING is incomplete, which makes this a one-sided
% check: it can miss a defect, it cannot invent one.
%
% WHERE IT SAMPLES. A dropped cell is a REGION, so a spread of directions and magnitudes finds it;
% the result's own vertices are added because a cell boundary is where a drop shows first, and the
% dual points of the input's vertices because those are where f*'s own cells meet.
    if numel(gs) < 2, return, end
    S = probePoints(g, obj);
    for i = 1:size(S,1)
        s = S(i,:);
        best = -inf;
        for k = 1:numel(gs)
            v = gs{k}.eval(s);
            if isfinite(v), best = max(best, v); end
        end
        if ~isfinite(best), continue, end
        got = g.eval(s);
        tol = 1e-7 * (1 + abs(best));
        if ~isfinite(got) || abs(got - best) > tol
            % A maxQuaPar: identifier on purpose -- it is what conjCPLQ's own fallback catch tests,
            % so a disagreement routes the whole domain to the symbolic Case C instead of returning
            % a number nobody checked.
            error('maxQuaPar:assemblyDisagreesWithPieces', ...
                ['the assembled Step 3 result is %.10g at (%.6g,%.6g) but the max of the pieces ' ...
                 'it was built from is %.10g. f* = max_k (q_k + I_P_k)* is an identity, so the ' ...
                 'assembly has dropped or over-extended a cell.'], got, s(1), s(2), best);
        end
    end
end

function S = probePoints(g, obj)
% objective: points at which to check the fold. Directions x magnitudes, plus the structural points
%   where a dropped cell shows up first.
    R = [0.25 1 3 10];
    th = (0:15) * (2*pi/16);
    S = zeros(0,2);
    for r = R
        S = [S; r*[cos(th).' sin(th).']]; %#ok<AGROW>
    end
    S = [S; 0 0];
    if ~isempty(g.V), S = [S; g.V]; end
    for i = 1:obj.nf                       % the dual points of the input's own vertices
        [L, Q, C] = QuaPol.matrixForm(obj.f(i,:));
        if ~isempty(C), continue, end
        for v = 1:size(obj.V,1)
            S = [S; (Q*obj.V(v,:).' + L).']; %#ok<AGROW>
        end
    end
    S = unique(round(S, 10), 'rows');
end

function [W, dF, dL] = faceBoundary(obj, i)
% objective: face i's boundary as a CCW vertex list plus, when it is unbounded, the two recession
%            directions -- the input conjConvexPolygon takes.
%
% Built from the EDGE INCIDENCE rather than from obj.P{i}. That list is documented clockwise and
% starting from a particular ray, and re-deriving the convention here is exactly the kind of
% off-by-one that produces a well-formed mesh on the wrong side. A convex face's boundary is a
% chain, so walking the incidence is unambiguous, and the orientation is then MEASURED.
    ej = find(any(obj.F == i, 2));
    if isempty(ej)
        error('PLQ:conjCPLQ:notImplemented', 'face %d has no edges.', i);
    end
    isR = obj.E(ej,3) == 0;
    rays = ej(isR); segs = ej(~isR);
    if ~ismember(numel(rays), [0 2])
        error('PLQ:conjCPLQ:notImplemented', ...
            'face %d has %d rays; a convex face has 0 or 2.', i, numel(rays));
    end

    if isempty(rays)
        iVs = unique([obj.E(segs,1); obj.E(segs,2)], 'stable');
        W = obj.V(iVs,:);
        ctr = mean(W, 1);
        [~, ord] = sort(atan2(W(:,2)-ctr(2), W(:,1)-ctr(1)));   % convex face: angle order IS CCW
        W = W(ord,:);
        dF = []; dL = [];
        return
    end

    % Unbounded: the two rays' BASE vertices are the ends of the vertex chain; walk the segments
    % between them.
    ends = obj.E(rays,1);
    W = obj.V(chainOrder(obj, segs, ends(1), ends(2)), :);
    dF = obj.V(obj.E(rays(1),2),:) - obj.V(obj.E(rays(1),1),:);
    dL = obj.V(obj.E(rays(2),2),:) - obj.V(obj.E(rays(2),1),:);
    dF = dF / norm(dF);  dL = dL / norm(dL);
    if orientedArea(W, dF, dL) < 0
        W = flipud(W);  tmp = dF; dF = dL; dL = tmp;
    end
end

function iVs = chainOrder(obj, segs, vStart, vEnd)
% objective: the vertices of a chain of segment edges, from vStart to vEnd.
    iVs = vStart;
    used = false(numel(segs),1);
    cur = vStart;
    guard = 0;
    while cur ~= vEnd
        guard = guard + 1;
        if guard > numel(segs) + 2
            error('PLQ:conjCPLQ:notImplemented', 'the boundary chain does not terminate.');
        end
        nxt = 0;
        for t = 1:numel(segs)
            if used(t), continue, end
            e = obj.E(segs(t),1:2);
            if e(1) == cur, nxt = e(2); used(t) = true; break
            elseif e(2) == cur, nxt = e(1); used(t) = true; break
            end
        end
        if nxt == 0
            error('PLQ:conjCPLQ:notImplemented', ...
                'the boundary chain of an unbounded face is not connected.');
        end
        iVs(end+1) = nxt; %#ok<AGROW>
        cur = nxt;
    end
end

function a = orientedArea(W, dF, dL)
% objective: the sign that says whether the walk [in along dF] -> W -> [out along dL] runs
%   counter-clockwise. Closing the unbounded set with two far points is enough: for a convex set
%   and R past every vertex, the sign of the closed polygon signed area is the set orientation.
    R = 1e3 * (1 + max(abs(W(:))));
    Q = [W(1,:) + R*dF; W; W(end,:) + R*dL];
    a = signedArea(Q);
end

function tris = faceTrianglesCCW(obj, i)
% Vertices of bounded face i in boundary order -> CCW polygon -> fan triangulation. The face of a
% QuaPol is convex, which is what makes a fan from one vertex valid. Same construction as
% convEnvCPLQ's file-local extractFaceTrianglesCCW, duplicated here because it is file-local
% there (as faceVertexIndices already is).
    ej = find(any(obj.F == i, 2));
    if any(obj.E(ej,3) == 0)
        error('PLQ:conjCPLQ:notImplemented', ...
            ['the fan-triangulation route needs a BOUNDED face (face %d is unbounded). A convex ' ...
             'unbounded face does not come here at all -- conjConvexPolygon takes it whole -- so ' ...
             'reaching this means the face is unbounded AND not positive definite.'], i);
    end
    W = obj.V(faceVertexIndices(obj, i), :);
    if signedArea(W) < 0, W = flipud(W); end
    n = size(W,1);
    if n < 3
        error('PLQ:conjCPLQ:notImplemented', 'face %d has fewer than 3 vertices.', i);
    end
    tris = cell(1, n-2);
    for t = 1:n-2, tris{t} = [W(1,:); W(t+1,:); W(t+2,:)]; end
end

function a = signedArea(W)
    x = W(:,1); y = W(:,2); n = size(W,1); a = 0;
    for i = 1:n, j = mod(i,n)+1; a = a + (x(i)*y(j) - x(j)*y(i)); end
    a = a/2;
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
function g = conjMaxOfSubTriangles(env, obj)
% Step 3: conjugate every triangle face of a multi-face convex-envelope RatPol (Step 1's output)
% independently (Step 2, conjFaceOrOriginal), then combine them via repeated pairwise maxQuaPar.
% Reached from conjSingleTriangle whenever Step 1 split the triangle (2- or 3-convex-edge).
%
% `obj` is the ORIGINAL piece Step 1 was given -- needed because a split leaves RATIONAL faces
% that Step 2 cannot conjugate directly; see conjFaceOrOriginal for why the original is an exact
% substitute for them.
% THE FOLD IS CROSS-CHECKED, for the reason assertFoldMatchesPieces gives at length: f* = max_k
% (q_k + I_{T_k})* is an IDENTITY, so the per-face conjugates are their own oracle and a fold that
% disagrees with them has dropped or over-extended a cell. conjPolygonalDomain has had this check
% since the 4-cone fan incident; this route did not, and that is how G4 stayed silent.
%
% MEASURED on the G4 fixture (sweep case 21, a triangle carrying x*y plus an affine part).
% Step 1 splits it into FOUR faces -- the nCE==3 cevian split, each half re-split -- two of them
% slivers of area 8.7e-05 against 2.7e-02. Every face's own Step 2 conjugate is EXACT at the bad
% point (1.032507658472, to 12 digits, four times over), and folding faces 2 and 3 keeps it; the
% fourth fold returns 1.005089907622. So the defect is not the per-piece closed form -- which is
% what `TODO.md` recorded -- but the assembly, and this is the check that says so.
    gs = cell(1, env.nf);
    for k = 1:env.nf
        gs{k} = toQuaPar(conjFaceOrOriginal(env, k, obj));
    end
    g = gs{1};
    for k = 2:env.nf
        g = maxQuaPar(g, gs{k});
    end
    % IT REFUSES, IT DOES NOT FALL BACK -- the same trade G6 made for the edge lower bound, and
    % for the same measured reason. A `maxQuaPar:` identifier is what conjSingleTriangle's own
    % catch tests, so raising one here would route this triangle to the symbolic Case C; on the
    % G4 fixture that did not finish in 25 MINUTES, against 0.4 s for the refusal, and nothing
    % says Case C's answer on it is the better one. Re-identify so the refusal propagates.
%
% UNDER THE SAME SWITCH AS THE EDGE BOUND. `CCA2_CONJ_VERIFY = 0` is the documented escape hatch
% for "conj checks its answer against the definition and refuses when it is violated", and this is
% exactly that kind of check, so it must answer to the same global -- otherwise turning the
% verification off no longer gives back the old answer, which is the one thing the hatch is for.
    global CCA2_CONJ_VERIFY %#ok<GVMIS>
    if ~isempty(CCA2_CONJ_VERIFY) && ~CCA2_CONJ_VERIFY, return, end
    try
        assertFoldMatchesPieces(g, gs, obj);
    catch ME
        if ~strcmp(ME.identifier, 'maxQuaPar:assemblyDisagreesWithPieces'), rethrow(ME); end
        error('PLQ:conjCPLQ:foldDroppedACell', '%s', ME.message);
    end
end

% ================================================================================================
function g = conjFaceOrOriginal(env, k, obj)
% Step 2 for ONE face of Step 1's envelope, conjugating the face itself when it can and the
% ORIGINAL quadratic over that face's domain when it cannot.
%
% WHY THE ORIGINAL IS AN EXACT SUBSTITUTE, not an approximation. Step 1 splits the input triangle
% T into sub-triangles T_k and puts on each the envelope of the ORIGINAL q over THAT sub-triangle:
% R_k = conv(q + I_{T_k}) (convEnvCPLQ's solveTriangleBF classifies T_k's own edges, not T's).
% Conjugation cannot see the difference between a function and its convex envelope, so
%
%       (R_k + I_{T_k})* = (conv(q + I_{T_k}))* = (q + I_{T_k})*
%
% exactly. Step 3 then maxes these, which is right for the same reason it is right at the top
% level: the T_k cover T, so conv(q+I_T) = min_k (R_k + I_{T_k}) and the conjugate of a min is the
% max of the conjugates.
%
% This matters because a split ALWAYS leaves at least one RATIONAL face (the 1-convex-edge
% sub-triangle's envelope is [COAP] A.3's quadratic-over-linear), and Step 2 has no rational
% branch -- that gap is what used to send the whole 2-/3-convex-edge case to cPLQ's symbolic
% Step 2/3. It does not need a rational-conjugate formula at all: the rational face is a
% RE-EXPRESSION of a piece whose conjugate Step 2 already computes in closed form (an indefinite
% quadratic over a 1-convex-edge triangle, conjIndefiniteQuadTriangle).
%
% Measured on V = conv{(2,1),(0,0),(1,0)}, f = xy -- the 2-convex-edge tightness split, whose
% face 2 is rational: the original-q route reproduces sup_{x in T_2} <s,x> - f(x) to 0 (exactly)
% at 6 dual points, and the assembled Step 3 result matches ground truth over the WHOLE triangle
% to 4.4e-16 at 8 dual points, including the paper's own flagged s = (-0.008727, -0.999962).
% That case previously required the symbolic engine end to end.
    face = extractTriangleFace(env, k);
    try
        g = conjPieceCPLQ(face);
        return
    catch ME
        if ~strcmp(ME.identifier, 'conjPieceCPLQ:notImplemented')
            rethrow(ME);
        end
    end
    if ~isa(obj, 'QuaPol')
        rethrow(ME);   % no original quadratic to fall back on (e.g. a direct RatPol conjugate)
    end
    g = conjPieceCPLQ(QuaPol(face.V, face.E, obj.f(1,5:10), face.F));
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
