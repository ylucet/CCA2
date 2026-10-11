function h = biconjQ(obj)
% biconjQ  The biconjugate f** = cl co f, computed EXACTLY over the rationals.
%
% objective: the closed convex envelope of a QuaPol, with every coefficient an exact rational and
%   every decision made by exact integer arithmetic.
%
% [input]  obj : QuaPol, operable (degree <= 2) and EXACT (it must carry fN/fD -- see
%                QuaPol.assertExact for why computing exactly on inexact data is worse than not
%                computing at all).
% [output] h   : QuaCon.
%
% ------------------------------------------------------------------------------------------------
% DIRECT, NOT conj(conj(f)). ALGORITHM.md states the principle and the measurement behind it:
% CCA2 has a direct algorithm for the envelope, and computing it as a double conjugation is a
% detour that cost 436 s on an input whose answer was "return f unchanged". The double conjugation
% would also need the conjugate of a QuaCon, whose return type is an open question (the plan's
% Phase 6), so the direct route is both faster and better defined.
%
% WHAT IS IMPLEMENTED
%
%   * f CONVEX on its own domain -> f itself. Nothing to compute: co f = f. Decided exactly, by
%     Sylvester on integers per piece plus `isConvexQ` for the multi-piece part (the gradient jump
%     across every shared edge, and the convexity of the domain).
%
%   * an EDGE-CONCAVE piece -> the affine interpolant over the LOWER HULL of the lifted vertices.
%     This is the case the convex envelope is cleanest for, and the reason is that such a function
%     on a polytope has its whole envelope determined by the VERTEX VALUES: the graph of co f is
%     the lower convex hull of {(v, q(v))}, so the computation is a hull of m points in R^3 with
%     rational coordinates -- no calculus, no calls to an optimiser, and every predicate an integer
%     sign test.
%
%   * UNBOUNDED pieces too, since 2026-09-05, by homogenising that hull: a lifted RECESSION RAY
%     joins the lifted vertices as a generator, and conv(vertices) + cone(rays) is the answer. It
%     is one code path rather than two -- with no rays every test below reduces to the bounded one
%     verbatim. The condition is that each recession direction d has d'Hd = 0 AND H*d = 0, which
%     makes the lifted set recede in the FIXED direction (d, L'd); H*d = 0 is automatic whenever H
%     is semidefinite, so only an indefinite piece can fail it, and that is refused by name.
%
% WHY THE ANSWER IS A QuaCon, AND WHAT THE MISSING TYPE ACTUALLY IS. Both implemented cases have
% RATIONAL face coefficients: the convex case returns f, and the concave case returns affine
% functions interpolating rational vertex values.
%
% CORRECTED 2026-09-05. This header used to say the missing faces were AFFINE cells
% `<p,x> - f*(p)` naming a dual vertex p of degree up to 4, and that the planned AlgAlg would hold
% them. That is true only for the case where THREE pieces tie at one isolated dual point. The
% commoner case is TWO pieces tying along a dual CURVE, and there p varies with x, so the face is
% not affine at all. Derived in closed form and checked against the primal definition by
% `.claude/bridgeFaceProbe.m`: for q1 = |v|^2/2 and q2 = |v-a|^2 + c the bridge face is
%
%       co f(x) = <2a, x - 2a> + 2*sqrt(|a|^2 - c) * |x - 2a| + 2c ,
%
% a NORM -- the square root of a quadratic. Measured there: the best affine fit to that face is off
% by 2.9 and the best quadratic fit by 0.36, against an oracle agreement of 1.5e-11. So `AlgAlg` as
% specified would not hold it.
%
% BUT THE FACE IS STILL ALGEBRAIC, AND THAT IS THE USABLE PART -- `.claude/bridgeAlgebraicProbe.m`.
% "Not a polynomial in x" is not "not algebraic": sqrt(Q) satisfies z^2 = Q. Written IMPLICITLY the
% face above is
%
%       (z - <2a, x-2a> - 2c)^2 - R^2*|x - 2a|^2 = 0 ,     R^2 = 4(|a|^2 - c),
%
% a QUADRIC CONE whose coefficients are RATIONAL, since R^2 is rational even though R is not; it
% holds to 1.8e-10 there. That is the same shape QuaCon already stores for its edge CONICS, one
% dimension up -- so the isotropic case needs an implicit quadric face plus a branch choice, not a
% new arithmetic.
%
% The general case is bigger but still implicit and still rational. With DIFFERENT Hessians
% (H2 = diag(2,4)) no quadric fits at all, yet a polynomial relation does exist: degree <= 4 in z
% -- which is what eliminating the minimiser t predicts, since dval/dt = 0 clears to a quartic in t
% -- with coefficients of degree 6 in (x,y). Measured as a COLLAPSE, the smallest singular value
% dropping from 1.2e-9 at degree 5 to 5.4e-14 at degree 6 and staying flat after, which is an
% identity rather than a better fit.
%
% THE RATIONAL/IRRATIONAL BOUNDARY IS SHARP, and the same script pins the other side of it: the
% bridge is quadratic exactly when the two tying cells have the SAME Hessian, because then the two
% conjugates differ by an affine function, the dual tie set is a LINE rather than a conic, and the
% sup of a concave quadratic over a line is again quadratic (measured residual 4.6e-12). Domain
% clipping breaks that condition even between pieces that share a Hessian -- an interior cell tying
% against a vertex cell has effective Hessians H and 0 -- which is why no sub-case of the coupling
% problem is rational enough to implement on its own. See DECISIONS.md 2026-09-05.

    if ~isa(obj, 'QuaPol')
        error('PLQ:biconjQ:input', ...
            'biconjQ takes a QuaPol (quadratic on polyhedral); got a %s.', class(obj));
    end

    % A CLOSED conic edge (E(:,3) = 2) has no endpoints, and the edge handling below is
    % written as "ray or segment". Refuse it by name rather than reading flag 2 as a
    % segment and indexing vertex 0. See RatCon.assertEdgeFlags. Nothing produces one
    % today, so this fires on nothing; it is here for when a bridge face is first stored.
    %
    % AFTER the input-type check, deliberately: a non-QuaPol has no mesh to inspect, and
    % reading obj.E first replaced this routine's named type refusal with an access error.
    RatCon.assertNoClosedEdges(obj.E, 'biconjQ');
    obj.assertOperable();
    obj.assertExact();

    nf = size(obj.fN, 1);

    % ---- a domain of dimension < 2, BEFORE the convex short-circuit -------------------------------------------------------------
    % co f is convex, so its domain is conv(dom f). For a needle or a single segment that is the
    % domain itself -- a point or a segment -- and both are THIN, which the H-form now expresses
    % with an EQUALITY side.
    %
    % BEFORE the convex short-circuit, not after, and the ordering is load-bearing: an affine or
    % PSD q on such a domain IS convex, so that branch would claim it and then read a face through
    % pieceGeometryB -- and these meshes have nf = 0 and an empty F, so there is no face to read.
    % Measured: the needle and the convex segment both raised noFace while the CONCAVE segment,
    % which falls past the short-circuit, worked.
    if obj.dom.dim < 2
        % A CHAIN of edges is NOT thin, and this used to refuse it as `ambiguousChain` on the
        % grounds that several function rows with no F left it undecided which edge took which.
        % That premise was simply wrong: for a chain, ROW j IS EDGE j's function -- QuaPol.eval
        % reads exactly obj.f(i,:) for edge i, and a chain must carry one row per edge for eval to
        % work at all. So the convention exists and is already relied on.
        %
        % Two non-collinear segments have a 2-DIMENSIONAL convex hull, so co f is an ordinary
        % lower-hull answer read from the whole mesh. That is concaveEnvelope, unchanged; only the
        % GATHER differs, iterating edges instead of faces.
        if isChainMesh(obj)
            assertChainEdgeConcave(obj);
            if isChain(obj)
                h = concaveEnvelope(obj);       % non-collinear: a 2-D hull
            else
                h = collinearChainEnvelope(obj); % collinear: a 1-D hull along its own span
            end
            return
        end
        h = thinEnvelope(obj);
        return
    end

    % ---- f CONVEX: there is nothing to compute -------------------------------------------------
    % The necessary-and-here-sufficient test for a SINGLE piece is that its Hessian is PSD. For a
    % multi-piece f that is necessary but not sufficient: the gradient jump across every shared
    % edge has to be consistent too, and the domain has to be convex. `isConvexQ` DECIDES that,
    % exactly (Singh-Lucet SIOPT 2021 Theorem 3.7), so the multi-piece case is now established
    % rather than asserted.
    %
    % fIsConvex STAYS, and only its role changes. It is a genuine short-circuit -- a caller that
    % already knows f is convex should not pay for the edge loop -- but it is no longer the only
    % route, so an honest input no longer has to be trusted to be answered. It is checked FIRST
    % for exactly that reason: `||` short-circuits, so isConvexQ runs only when nothing cheaper
    % has settled the question.
    allPSD = true;
    for k = 1:nf
        Hn = [obj.fN(k,5) obj.fN(k,6); obj.fN(k,6) obj.fN(k,7)];
        if ~ratQ.isPSD2(Hn), allPSD = false; break, end
    end
    flagged = isequal(obj.fIsConvex, true);
    if flagged, assertFlagNotRefuted(obj); end
    if allPSD && (nf == 1 || flagged || isConvexQ(obj))
        if obj.nv == 0
            % The whole plane, one face: co f = f with no constraints at all. QuaPol.energy is
            % exactly this, and it used to raise noFace because quaPolAsQuaCon looks for a FACE and
            % a full-domain mesh has none.
            h = assembleQuaConCells(struct('num', obj.fN(1,:), 'den', obj.fD(1), ...
                                           'con', zeros(0,7)));
            return
        end
        h = quaPolAsQuaCon(obj);
        return
    end
    if isequal(obj.fIsConvex, true) && ~allPSD
        % The free necessary condition contradicts the caller's assertion. Refused LOUDLY rather
        % than trusted -- ALGORITHM.md's rule, because the failure it prevents is silent: biconj
        % would return a non-convex f as its own convex envelope.
        error('PLQ:biconjQ:notConvex', ...
            ['fIsConvex is asserted true but some piece has a Hessian that is not positive ' ...
             'semidefinite, which is a free necessary condition for convexity.']);
    end

    % ---- THE ASSEMBLY ROUTE, for a TWO-PIECE input the edge-concavity path cannot take ---------
    % Where every piece is convex but f is not, the envelope must BRIDGE the gradient jump, and the
    % lower-hull machinery below cannot: a bridge is not determined by the vertex values. That case
    % used to be refused outright (`convexPieceCoupling`). It is now COMPUTED, by assembling the
    % envelope from the two pieces' conjugates and meshing the result -- `envelopeAssemble`,
    % `envelopePrune`, `envelopeMesh`, each verified in its own suite.
    %
    % It is attempted only for TWO pieces, because that is what the assembly takes, and it is
    % attempted BEFORE the per-piece checks so that their refusal is not raised for a case that now
    % has an answer. A failure falls through to those checks, so nothing that used to work changes
    % and nothing that used to refuse refuses for a WORSE reason.
    assemblyWhy = '';
    if allPSD && nf == 2
        try
            h = assembledEnvelope(obj);
            return
        catch ME
            if ~startsWith(ME.identifier, 'envelope') && ...
               ~startsWith(ME.identifier, 'bridgeFace') && ...
               ~startsWith(ME.identifier, 'bridgeArcCut') && ...
               ~startsWith(ME.identifier, 'cornerFace')
                rethrow(ME)
            end
            % fall through: the checks below name the case in the vocabulary a caller expects --
            % but the refusal CARRIES this reason, because "the assembly declined" and "there is no
            % route" are different news and only the first says what to fix.
            assemblyWhy = sprintf(' The assembly route declined first: %s -- %s', ME.identifier, ...
                                  ME.message(1:min(end, 200)));
        end
    end

    % ---- EDGE-CONCAVITY, per face, which is weaker than concavity and is the real condition ----
    % The envelope is the lower hull of the lifted vertices whenever q is concave along every EDGE
    % DIRECTION of its face -- H itself need not be negative semidefinite. Why the weaker condition
    % suffices, on a triangle and hence on each cell of the hull: q - L vanishes at the three
    % corners and is concave along each edge, so q >= L on the BOUNDARY; and with H not positive
    % semidefinite q - L has no interior minimum, so its minimum over the closed triangle is on
    % that boundary. (The classical vertex-polyhedral / edge-concave criterion, verified here
    % rather than taken on trust -- .claude/edgeconcave-sweep.m checks the envelope's four defining
    % properties on 124 cases, 48 of them genuinely indefinite.)
    %
    % IT IS WHERE McCORMICK LIVES. q = xy on the unit square has H = [0 1; 1 0], genuinely
    % indefinite, yet d'Hd = 0 along both edge directions -- so its envelope is the lower hull of
    % the four corner values, which is exactly max(0, x+y-1).
    %
    % MULTI-PIECE WORKS FOR THE SAME REASON, and this is why the envelope needs no fold. Each
    % face's graph has its lower hull spanned by that face's lifted vertices, and the convex hull
    % of a union is the hull of the union of the spanning sets -- so the whole envelope is the
    % lower hull of ALL the lifted vertices at once. The envelope COUPLES the pieces, which is why
    % it is not a max like the conjugate; taking one hull over everything is precisely that
    % coupling, done in one step.
    for k = 1:nf
        Hk = [obj.fN(k,5) obj.fN(k,6); obj.fN(k,6) obj.fN(k,7)];
        [vsk, rays, bnd] = pieceGeometryB(obj, k);
        if ~bnd
            % TWO DIFFERENT SITUATIONS, and saying which is which is the point of this branch.
            % Along a recession direction d, q behaves as t*<grad q, d> + (t^2/2) d'Hd, so the
            % envelope runs to -infinity exactly when some recession direction has d'Hd < 0, or has
            % d'Hd = 0 with the linear rate negative somewhere on the piece. Otherwise q is bounded
            % below on its own recession cone and the envelope is FINITE -- and since 2026-09-05
            % it is COMPUTED, by concaveEnvelope's homogenised hull, provided the lifted recession
            % ray is a fixed direction (the shear test below).
            Lk = [obj.fN(k,8); obj.fN(k,9)];
            % THE CONE, NOT ONLY ITS GENERATORS. On cone(d1,d2) the curvature is
            % a*l1^2 + 2b*l1*l2 + c*l2^2, and it is nonnegative for every l >= 0 exactly when
            % a >= 0, c >= 0 and (b >= 0 or b^2 <= a*c). Checking the generators alone would miss a
            % direction strictly INSIDE the cone that runs downward -- the same flaw conjQ's
            % recessionConditions had to fix.
            %
            % AND THE RATE IS READ ON THE PIECE, NOT AT THE ORIGIN. With d'Hd = 0 the function
            % along the ray is AFFINE in t with rate (H x + L)'d, which depends on x as soon as
            % H*d is nonzero -- and H*d can be nonzero here precisely because this branch is the
            % NON-convex one, so H may be indefinite (for a semidefinite H, d'Hd = 0 forces
            % H*d = 0 and the rate is the constant L'd). This test used to read L'd, the rate at
            % the ORIGIN, which is the right number only for a piece whose apex sits there.
            % Measured witness: q = x*y on {x >= 0, y >= -5} has f(t,-5) = -5t, so co f is
            % -infinity, and every other test here passes it -- .claude/recessionRateProbe.m and
            % biconjQTest/theDIVERGENCETestReadsTheRateONThePieceNotAtTheOrigin.
            %
            % The rate is affine in x, so its infimum over the piece is attained at a VERTEX unless
            % it also falls along a recession direction -- and that second possibility is exactly
            % what the cone test below decides, since d'H e is the rate's own slope along e.
            [Vi, vd] = ratQ.combineDen(obj.VN, obj.VD);
            diverges = false;
            for r = 1:size(rays,1)
                d = rays(r,:).';
                a = ratQ.chkB(d.' * Hk * d, abs(d).' * abs(Hk) * abs(d), 'recession curvature');
                if a < 0, diverges = true; break, end
                if a == 0
                    % (H*v + L)'d at every vertex v of the piece, cleared by the POSITIVE vd.
                    for i = 1:numel(vsk)
                        v = Vi(vsk(i), :).';
                        rate = ratQ.chkB((Hk*v + vd*Lk).' * d, (abs(Hk)*abs(v) + abs(vd)*abs(Lk)).' * abs(d), 'recession slope at a vertex');
                        if rate < 0, diverges = true; break, end
                    end
                    if diverges, break, end
                end
            end
            if ~diverges && size(rays,1) >= 2
                aR = abs(rays);  aH = abs(Hk);
                a = ratQ.chkB(rays(1,:)*Hk*rays(1,:).', aR(1,:)*aH*aR(1,:).', 'cone a');
                b = ratQ.chkB(rays(1,:)*Hk*rays(2,:).', aR(1,:)*aH*aR(2,:).', 'cone b');
                c = ratQ.chkB(rays(2,:)*Hk*rays(2,:).', aR(2,:)*aH*aR(2,:).', 'cone c');
                if b < 0 && ratQ.chkB(b^2 - a*c, b^2 + abs(a)*abs(c), 'cone discriminant') > 0
                    diverges = true;
                end
            end
            if diverges
                error('PLQ:biconjQ:minusInfinity', ...
                    ['piece %d is unbounded and q runs DOWNWARD along one of its recession ' ...
                     'directions, so co f is -infinity there. That is the right answer and there ' ...
                     'is nowhere to put it -- the same representational gap conjCPLQ records as ' ...
                     'convEnvUnbounded:minusInfinity.'], k);
            end
            % NOTHING MORE IS NEEDED HERE. The lifted recession ray is (d, rho) with rho the
            % MINIMUM of q's rate over the piece, and concaveEnvelope computes that minimum -- see
            % the derivation there. This used to raise PLQ:biconjQ:shearedRecession whenever H*d
            % was nonzero, i.e. whenever the rate varied; that refusal was unnecessary, since a
            % varying rate needs a minimum rather than a new algorithm.
            %
            % A recession direction with d'Hd > 0 is NOT handled here and must not be: q then grows
            % quadratically along that ray, the piece is not edge-concave along it, and the
            % edge-curvature test below is the one that should speak. It names the obstruction,
            % where a boundedness-first refusal would only name the order two checks happen to run
            % in.
        end
        own = find(any(obj.F == k, 2));
        for j = own.'
            a = obj.E(j,1);  b = obj.E(j,2);
            d = obj.VN(b,:)/obj.VD(b) - obj.VN(a,:)/obj.VD(a);
            [dn, ~] = ratQ.fromDouble(d);
            if ratQ.chkB(dn * Hk * dn.', abs(dn) * abs(Hk) * abs(dn).', 'edge curvature') > 0
                % TWO DIFFERENT SITUATIONS TRIP THIS ONE TEST, and until 2026-09-05 they shared a
                % message that was right about only one of them. Which it is turns on whether THIS
                % piece is convex, so that is what selects the identifier.
                if allPSD
                    % EVERY piece is convex, so no piece needs [COAP] at all -- each one's own
                    % envelope is itself. What is missing is the hull ACROSS pieces: f is not
                    % convex (isConvexQ said so above), so co f has to bridge the gradient jump,
                    % and a bridge is not determined by the vertex values either.
                    %
                    % The bridge face is NOT rational, and it is not affine either -- see
                    % DECISIONS.md 2026-09-05 and .claude/bridgeFaceProbe.m, which derive it in
                    % closed form for two convex pieces and check it against the primal
                    % definition. It is affine PLUS THE NORM of an affine image of x, so the type
                    % that would hold it does not exist yet and is not the planned AlgAlg.
                    error('PLQ:biconjQ:convexPieceCoupling', ...
                        ['every piece is convex but f is not, so co f must BRIDGE the gradient ' ...
                         'jump across a shared edge (piece %d curves upward along one of its ' ...
                         'own edges, so its vertex values do not determine the bridge). The ' ...
                         'bridge face carries the square root of a quadratic, which no CCA2 ' ...
                         'type can hold -- see DECISIONS.md 2026-09-05.%s'], k, assemblyWhy);
                end
                % THE ATOM ROUTE (2026-10-10). A non-convex piece's conjugate is that of its
                % vertices and convex edges (erratum Theorem 1), each of which has a polyhedral
                % conjugate, so the envelope is the minimum over the atoms' pairwise envelopes and
                % triple tangent planes -- biconjAtoms' header. Bounded domains only; the result is
                % an EnvelopeMin, which answers eval but is not meshed into one AlgCon.
                if all(obj.E(:,3) ~= 0)
                    h = biconjAtoms(obj);
                    return
                end
                error('PLQ:biconjQ:notImplemented', ...
                    ['piece %d is not convex AND curves UPWARD along one of its edges, and the ' ...
                     'domain is UNBOUNDED -- the atom route covers bounded domains only (an ' ...
                     'unbounded one needs the erratum''s Case E).'], k);
            end
        end
        % A 2-D face needs three affinely independent GENERATORS, and for an unbounded piece two of
        % them may be recession directions -- a cone has one vertex and two rays and is perfectly
        % 2-dimensional, so counting vertices alone would refuse it.
        if bnd
            if numel(vsk) < 3
                error('PLQ:biconjQ:degenerate', ...
                    'piece %d has %d vertices; a 2-D bounded face needs at least three.', ...
                    k, numel(vsk));
            end
        elseif numel(vsk) + size(rays,1) < 3
            error('PLQ:biconjQ:degenerate', ...
                'piece %d has %d vertices and %d rays; a 2-D face needs three generators.', ...
                k, numel(vsk), size(rays,1));
        end
    end

    h = concaveEnvelope(obj);
end

% ==================================================================================================

% ------------------------------------------------------------------------------------------------
function h = assembledEnvelope(obj)
% co f for a TWO-PIECE input, assembled from the pieces' conjugates.
%
% f** = sup_s [<s,x> - f*(s)] and f* = max over the pieces of each piece's conjugate, so the envelope
% splits over the overlay of the two conjugates' cells. `envelopeAssemble` walks those cell pairs,
% `envelopePrune` drops the faces that can never be the maximum, and `envelopeMesh` builds the
% arrangement -- which is the only route that does not form a partial max, and partial maxima are
% what made every fold order explode (DECISIONS.md 2026-09-18).
%
% TWO MESHES, and the cheaper one is tried first. `envelopeMesh` returns a QuaCon, which is what a
% caller wants whenever the answer fits: every face quadratic, every edge a conic. It refuses --
% `envelopeMesh:bridgeCompetes` -- exactly when a BRIDGE face's value competes with a quadratic
% one, because that boundary has degree four and no conic subdivision holds it as an edge.
%
% For that case `envelopeMeshAlg` builds the same arrangement with such an edge allowed, and
% `assembleAlgCrv` puts it in an `AlgCrv` -- algebraic faces on an algebraic subdivision, the type
% the subdivision axis grew for (Crv.m). So the RETURN TYPE says which case the input was, and both
% are verified against `envelopeEval` in their own suites.
    pieces = splitPieces(obj);
    gA = conjQ(pieces{1});
    gB = conjQ(pieces{2});
    env = envelopePrune(envelopeAssemble(gA, gB));
    try
        h = envelopeMesh(env);
    catch ME
        switch ME.identifier
            case 'envelopeMesh:bridgeCompetes'
                try
                    h = assembleAlgCrv(envelopeMeshAlg(env));
                catch ME2
                    if ~strcmp(ME2.identifier, 'envelopeMeshAlg:faceKind'), rethrow(ME2), end
                    h = unmeshed(env);
                end
            case {'envelopeMesh:bridgeFaces', 'envelopeMesh:arrangementSize'}
                % A corner face whose wedge has no rows cannot be placed in any arrangement, and an
                % arrangement past the curve budget does not finish (TODO). The face list is still
                % the envelope -- envelopeEval answers for it, verified against the definition --
                % so return it unmeshed rather than refuse.
                h = unmeshed(env);
            otherwise
                rethrow(ME)
        end
    end
end

function h = unmeshed(env)
% The pair's face list as an EnvelopeMin of one pair: eval only, no mesh (see EnvelopeMin.m).
    h = EnvelopeMin({env}, struct('T', {}, 's', {}, 'c', {}), {});
end

% ------------------------------------------------------------------------------------------------
function P = splitPieces(obj)
% The two pieces as single-face QuaPols, each carrying its own face and that face's own polygon.
    P = cell(1,2);
    for k = 1:2
        own = find(any(obj.F == k, 2));
        vs = unique(reshape(obj.E(own,1:2), [], 1));
        vs = vs(vs > 0);
        V = obj.V(vs,:);
        idx = zeros(max(vs),1);
        idx(vs) = 1:numel(vs);
        E = obj.E(own,:);
        E(:,1:2) = idx(E(:,1:2));
        F = zeros(size(E,1), 2);
        for j = 1:size(E,1)
            F(j,:) = (obj.F(own(j),:) == k);      % 1 where this piece is, 0 for the outside
        end
        P{k} = QuaPol(V, E, obj.f(k,:), F);
    end
end

% ------------------------------------------------------------------------------------------------
function assertFlagNotRefuted(obj)
% Under the global BICONJ_ASSERT, REFUTE a false fIsConvex outright instead of trusting it.
% Off by default, in the established MAXQP_ASSERT / QUAPAR_VALIDATE style: tests and investigations
% turn it on, production runs pay nothing.
%
% WHY IT IS OPT-IN RATHER THAN ALWAYS ON. The flag exists precisely to SKIP the edge loop that
% isConvexQ runs, so verifying it on every call would delete the only thing it buys. What it costs
% to leave off is exactly one thing, and it is measured: QuaPol.examples{16} has a non-convex
% domain, the legacy dom.isConvex says otherwise, and with the flag asserted biconjQ returned that
% f as its own convex envelope -- silently. DECISIONS.md 2026-09-05 carries the witness.
%
% The free necessary condition (every Hessian PSD) is checked on EVERY call regardless, above and
% below this; this raises the same identifier because to a caller it is the same mistake, only
% detected by a stronger test.
    global BICONJ_ASSERT %#ok<GVMIS>
    if isempty(BICONJ_ASSERT) || BICONJ_ASSERT < 1, return, end
    [tf, why] = isConvexQ(obj);
    if ~tf
        error('PLQ:biconjQ:notConvex', ...
            ['fIsConvex is asserted true but the exact test refutes it: %s. ' ...
             '(BICONJ_ASSERT is on; unset it to trust the flag again.)'], why);
    end
end

function h = concaveEnvelope(obj)
% objective: co(q + I_P) for q CONCAVE or AFFINE on a bounded polygon P, exactly.
%
% THE MATHEMATICS IS ONE SENTENCE. A concave function on a polytope lies above every chord, so the
% lower convex hull of its graph is the lower convex hull of the finitely many LIFTED VERTICES
% {(v_i, q(v_i))} -- and that hull is the graph of co f. A convex piecewise-affine function is the
% MAX of the affine functions of its facets, so
%       co f (x) = max over lower-hull facets of  l_F(x),      x in P,
% and +infinity off P. Both halves are exactly the structure conjQ already builds: a max of affine
% functions, intersected with the polygon's own half-planes.
%
% EVERYTHING IS INTEGER. The lifted points are put over one denominator, a candidate facet is the
% plane through three of them (an integer cross product), and "is this a LOWER hull facet" is the
% sign of an integer for every other point. No hull library, no orientation predicate that could be
% read backwards, and no tolerance.
    [Vi, vd] = ratQ.combineDen(obj.VN, obj.VD);
    nf = size(obj.fN, 1);
    % For a CHAIN the loop below runs over EDGES rather than faces -- see the gather. nf is a count
    % of function ROWS either way, and a chain carries one per edge, so the loop bound is the same.
    chain = isChain(obj);

    % ---- the lifted vertices of EVERY face, over ONE denominator --------------------------------
    % The faces may carry different denominators, so the common one is 2*vd^2*lcm(fD). A vertex
    % shared by two faces gets the SMALLER of the two values: the envelope is of cl f, and the
    % lower hull is what the smaller lifted point contributes. For a continuous f they agree, and
    % where they do not, f is discontinuous there and the lower value is the closure's.
    Lden = 1;
    for k = 1:nf, Lden = lcm(Lden, obj.fD(k)); end
    Dl = ratQ.chk(2 * Lden * vd^2, 'lift denominator');

    lifted = containers.Map('KeyType', 'double', 'ValueType', 'double');
    R3 = zeros(0,3);
    for k = 1:nf
        [fN, fD] = obj.faceQ(k);
        Hn = [fN(5) fN(6); fN(6) fN(7)];
        Ln = [fN(8); fN(9)];
        kn = fN(10);
        scale = ratQ.chk(Lden / fD, 'denominator scaling');
        if chain
            % A CHAIN's units are EDGES, not faces: row k of fN is edge k's function (QuaPol.eval
            % reads it that way), the edge's own endpoints are its vertices, and a ray edge
            % contributes its direction. pieceGeometryB cannot be asked, since it looks up F and a
            % chain's F is all zeros.
            if obj.E(k,3) ~= 0
                vsk = [obj.E(k,1); obj.E(k,2)];  raysk = zeros(0,2);
            else
                vsk = obj.E(k,1);
                raysk = ratQ.chk(Vi(obj.E(k,2),:) - Vi(obj.E(k,1),:), 'chain ray direction');
            end
        else
            [vsk, raysk, ~] = pieceGeometryB(obj, k);
        end
        % THE LIFTED RECESSION RAYS, AND THE HEIGHT IS A MINIMUM. Along a recession direction d,
        % q(x + t d) = q(x) + t*r(x,d) with r(x,d) = (H x + L)'d, since d'Hd = 0 is already forced
        % here. The lifted RAY of the epigraph is then (d, rho) for any rho >= r(x,d) at SOME x of
        % the piece -- and a direction that recedes from one point of a convex set recedes from all
        % of them, so the smallest admissible height is
        %
        %       rho = min over the piece of r(.,d),
        %
        % which is a MINIMUM OVER THE PIECE'S VERTICES because r is affine in x and (the cone test
        % in the piece loop) does not fall along either recession direction.
        %
        % WHY THIS IS THE WHOLE HULL, and not merely a minorant. Writing x = sum(lam_i v_i) +
        % sum(mu_j d_j),
        %
        %   q(x) = q(sum lam v) + sum_j mu_j * r(sum lam v, d_j) + sum_{j,k} mu_j mu_k d_j'H d_k / 2
        %
        % using d_j'H d_j = 0. The last sum is >= 0 because the cone test established d_1'H d_2 >= 0;
        % r is affine so r(sum lam v, d) >= min_v r(v,d) = rho; and q(sum lam v) is above the lower
        % hull of the lifted vertices by the edge-concavity argument. So q >= the hull everywhere,
        % with equality attained -- which makes the hull the envelope rather than a bound.
        %
        % REPLACED the constant L'd on 2026-09-05. That was the special case H*d = 0, where r does
        % not vary; the general one used to be refused as PLQ:biconjQ:shearedRecession, and the
        % refusal was unnecessary -- the varying rate needs a minimum, not a new algorithm.
        % q = x*y on the first quadrant is the smallest instance: r((0,0),(1,0)) = 0 gives the
        % lifted rays (1,0,0) and (0,1,0) and hence co f = 0, which is right.
        %
        % The scaling matches Q3's below: both are the true (x, y, q) multiplied by the positive
        % 2*Lden*vd^2, and for a DIRECTION only the positive multiple matters.
        for r = 1:size(raysk,1)
            d = raysk(r,:).';
            rho = [];
            for i = 1:numel(vsk)
                v = Vi(vsk(i), :).';
                ri = ratQ.chkB((Hn*v + vd*Ln).' * d, (abs(Hn)*abs(v) + abs(vd)*abs(Ln)).' * abs(d), 'lifted ray height');
                if isempty(rho) || ri < rho, rho = ri; end
            end
            R3(end+1,:) = [ratQ.chk(2*Lden*vd*d(1), 'lift ray x'), ...
                           ratQ.chk(2*Lden*vd*d(2), 'lift ray y'), ...
                           ratQ.chk(2*scale*rho, 'lift ray z')]; %#ok<AGROW>
        end
        for i = 1:numel(vsk)
            v  = Vi(vsk(i), :).';
            qv = ratQ.chkB(scale * (v.'*Hn*v + 2*vd*(Ln.'*v) + 2*vd^2*kn), abs(scale) * (abs(v).'*abs(Hn)*abs(v) + 2*abs(vd)*(abs(Ln).'*abs(v)) + 2*vd^2*abs(kn)), 'vertex value');
            key = vsk(i);
            if ~isKey(lifted, key) || qv < lifted(key)
                lifted(key) = qv;
            end
        end
    end

    vs = cell2mat(keys(lifted)).';
    m  = numel(vs);
    R3 = unique(R3, 'rows');
    nr = size(R3, 1);
    if m + nr < 3 || m < 1
        error('PLQ:biconjQ:degenerate', ...
            'the domain has %d vertices and %d recession rays; a 2-D one needs three generators.', ...
            m, nr);
    end
    Q3 = zeros(m, 3);
    for i = 1:m
        v = Vi(vs(i), :).';
        Q3(i,:) = [ratQ.chk(2*Lden*vd*v(1), 'lift x'), ...
                   ratQ.chk(2*Lden*vd*v(2), 'lift y'), lifted(vs(i))];
    end

    % ---- the LOWER hull facets, by brute force over triples of GENERATORS -----------------------
    % A generator is a lifted VERTEX or a lifted RECESSION RAY, and the hull is
    % conv(vertices) + cone(rays). Homogenising that way is the whole extension to unbounded
    % domains: a facet is spanned by three generators of which at least one must be a vertex (a
    % plane needs a point to fix its offset), and the two IN-PLANE DIRECTIONS come from differences
    % of vertices and from ray directions interchangeably.
    %
    % The validity test splits the same way. A supporting plane must have every lifted vertex on or
    % above it -- n.p >= d -- and must not be undercut at infinity, which is the HOMOGENEOUS form
    % of the same condition: n.r >= 0 for every lifted ray. With no rays present both reduce to the
    % bounded test that was here before, so this is one code path rather than two.
    %
    % m + nr is a handful, so O((m+nr)^4) is a few hundred integer operations and needs no
    % incremental hull. Every test is exact.
    G   = [Q3; R3];
    isP = [true(m,1); false(nr,1)];
    ng  = m + nr;
    facets = zeros(0, 4);                       % rows [a b c d] with a X + b Y + c Z = d, c > 0
    for i = 1:ng
        for j = i+1:ng
            for k = j+1:ng
                idx = [i j k];
                pts = idx(isP(idx));
                if isempty(pts), continue, end   % three rays fix no offset: not a facet
                dirs = zeros(0,3);
                base = pts(1);
                for t = idx
                    if t == base, continue, end
                    if isP(t)
                        dirs(end+1,:) = G(t,:) - G(base,:); %#ok<AGROW>
                    else
                        dirs(end+1,:) = G(t,:);             %#ok<AGROW>
                    end
                end
                n = ratQ.chk(cross(dirs(1,:), dirs(2,:)), 'facet normal');
                if n(3) == 0, continue, end     % vertical plane: not the graph of a function
                if n(3) < 0, n = -n; end        % orient upward, so "above" is n.p >= d
                d = ratQ.chkB(n * G(base,:).', abs(n) * abs(G(base,:)).', 'facet offset');
                below = false;
                for p = 1:m
                    if ratQ.chkB(n * Q3(p,:).' - d, abs(n) * abs(Q3(p,:)).' + abs(d), 'facet side') < 0, below = true; break, end
                end
                if ~below
                    for p = 1:nr
                        if ratQ.chkB(n * R3(p,:).', abs(n) * abs(R3(p,:)).', 'facet side at infinity') < 0
                            below = true;  break
                        end
                    end
                end
                if below, continue, end          % undercut by a vertex or at infinity
                % Reduce to a primitive integer row so that one plane found through several
                % different triples is ONE row and `unique` collapses it -- the same canonical-form
                % argument ratQ.conic makes for curves.
                g = gcd(gcd(abs(n(1)), abs(n(2))), gcd(abs(n(3)), abs(d)));
                if g == 0, g = 1; end
                facets(end+1, :) = [n, d] / g; %#ok<AGROW>
            end
        end
    end
    facets = unique(facets, 'rows');
    if isempty(facets)
        error('PLQ:biconjQ:noFacet', ...
            'no lower-hull facet was found, which cannot happen for three non-collinear vertices.');
    end

    % ---- each facet's affine function, and the cell where it is the largest ---------------------
    % From a X + b Y + c Z = d in the lifted coordinates, the function is
    %       l(x,y) = (d - a*Dl*x - b*Dl*y) / (Dl*c)  ... written over the denominator Dl*c.
    nfac = size(facets,1);
    % THE ENVELOPE'S DOMAIN IS conv(P), NOT P, and for a non-convex face those differ. co f must
    % be CONVEX, so its domain is convex: the closed convex hull of dom f. For a convex face the
    % two coincide and nothing changes; for a non-convex one, using the face's own edges produced a
    % region that was not even the face -- measured on an L, where biconjQ returned +Inf at (2,0),
    % (2,1), (1,2) and (1.5,0.5), all of which are IN the L and where the envelope is finite.
    %
    % The FACE FUNCTIONS needed no change: the lower hull of the lifted vertices is already the
    % right answer there, reflex vertices included. Checked against an LP over a dense sample of
    % the true L -- the hull and the LP agree at every probe, inside the L and in conv(L) beyond
    % it, so only the domain was wrong.
    % conv(dom f) INCLUDING ITS RECESSION CONE. R3(:,1:2) are the lifted rays' shadows, which are
    % exactly the primal recession directions, so the same list serves both halves of the hull.
    dom = hullHalfPlanes(Vi, vd, vs, R3(:,1:2));
    cells = struct('num', {}, 'den', {}, 'con', {});
    for t = 1:nfac
        a = facets(t,1);  b = facets(t,2);  c = facets(t,3);  d = facets(t,4);
        num = [0 0 0 0, 0, 0, 0, ratQ.chk(-a*Dl,'c8'), ratQ.chk(-b*Dl,'c9'), d];
        den = ratQ.chk(Dl * c, 'facet denominator');

        rows = dom;
        for s = 1:nfac
            if s == t, continue, end
            [dn, ~] = ratQ.sub(num, den, ...
                [0 0 0 0, 0, 0, 0, ratQ.chk(-facets(s,1)*Dl,'c8'), ...
                 ratQ.chk(-facets(s,2)*Dl,'c9'), facets(s,4)], ...
                ratQ.chk(Dl * facets(s,3), 'facet denominator'));
            if all(dn == 0), continue, end
            rows(end+1,:) = [ratQ.chk(2*dn(8),'d1'), ratQ.chk(2*dn(9),'d2'), ...
                             ratQ.chk(2*dn(10),'d0')]; %#ok<AGROW>
        end
        if ~ratQ.feasible2(rows, true), continue, end   % this facet never wins inside P
        con = zeros(size(rows,1), 7);
        for r = 1:size(rows,1)
            con(r,:) = [ratQ.conic([0 0 0, rows(r,:)]), sgnOfB(rows(r,:), +1)];
        end
        cells(end+1) = struct('num', num, 'den', den, 'con', con); %#ok<AGROW>
    end

    h = assembleQuaConCells(cells);
end

function rows = hullHalfPlanes(Vi, vd, vs, dirs)
% objective: the half-planes of conv(vertices) + cone(dirs), exactly.
% [input]  dirs : (optional) r x 2 integer recession directions over the same denominator vd.
%                 Empty or omitted gives the bounded hull.
% [output] rows : h x 3, [a b c] meaning a*x + b*y + c >= 0 on the hull, in ACTUAL coordinates
%
% A candidate line is spanned either by TWO VERTICES or by one vertex and a RECESSION DIRECTION,
% and it bounds the hull exactly when every vertex lies weakly on one side of it AND every
% recession direction points weakly to that same side. The direction test is the HOMOGENEOUS form
% of the vertex test -- the offset drops out -- which is the whole of the extension to unbounded
% domains; with no directions it reduces to the bounded test verbatim.
%
% Omitting the direction test would keep only the lines through vertex PAIRS and would accept ones
% the cone crosses, so an unbounded domain would come back as its bounded hull: a domain strictly
% too small, and +infinity returned where the envelope is finite.
%
% Every test is an integer sign, so the hull is DECIDED. O((m+r)*m^2) for a handful of corners,
% which needs no incremental hull algorithm and has no orientation convention to read backwards.
%
% Numerator coordinates decide the side; the emitted row is in actual coordinates (linear part
% scaled by vd). See conjQ's pieceShape for why mixing the two yields a parallel line.
    if nargin < 4, dirs = zeros(0,2); end
    rows = zeros(0,3);
    m = numel(vs);
    r = size(dirs,1);
    for a = 1:m
        for b = 1:(m + r)
            p = Vi(vs(a), :);
            if b <= m
                if a == b, continue, end
                n = [Vi(vs(b),2)-p(2), -(Vi(vs(b),1)-p(1))];
            else
                e = dirs(b-m, :);
                n = [e(2), -e(1)];
            end
            if all(n == 0), continue, end                 % coincident vertices, or a null direction
            c0 = ratQ.chkB(-(n * p.'), (abs(n) * abs(p).'), 'offset');
            sgn = 0;  ok = true;
            for i = 1:m
                t = ratQ.chkB(n(1)*Vi(vs(i),1) + n(2)*Vi(vs(i),2) + c0, abs(n(1)*Vi(vs(i),1)) + abs(n(2)*Vi(vs(i),2)) + abs(c0), 'hull side');
                if t == 0, continue, end
                if sgn == 0, sgn = sign(t);
                elseif sign(t) ~= sgn, ok = false; break
                end
            end
            if ok
                for i = 1:r
                    t = ratQ.chkB(n * dirs(i,:).', abs(n) * abs(dirs(i,:)).', 'hull side at infinity');
                    if t == 0, continue, end
                    if sgn == 0, sgn = sign(t);
                    elseif sign(t) ~= sgn, ok = false; break
                    end
                end
            end
            if ~ok || sgn == 0, continue, end             % generators on both sides: not a facet
            row = sgn * [ratQ.chk(vd*n(1), 'hull normal'), ...
                         ratQ.chk(vd*n(2), 'hull normal'), c0];
            row = row / gcdRow(row);
            if ~any(all(rows == row, 2)), rows(end+1,:) = row; end %#ok<AGROW>
        end
    end
end

function g = gcdRow(row)
    g = 0;
    for i = 1:numel(row)
        if row(i) ~= 0, g = gcd(g, abs(row(i))); end
    end
    if g == 0, g = 1; end
end

function rows = polygonHalfPlanes(obj, Vi, vd, vs, k, rays)
% objective: face k's own constraints, as exact half-planes [a b c] meaning a*x + b*y + c >= 0 at
%            the ACTUAL coordinates (x,y).
%
% TWO COORDINATE SYSTEMS MEET HERE AND MUST NOT BE MIXED. `Vi` holds vertex NUMERATORS over the
% shared denominator vd, so a point's actual coordinates are X/vd. With n a normal to the edge and
% c0 = -n.p taken at numerators, one and the same line reads
%
%       n.X + c0 = 0        in numerator coordinates
%       vd*n.x + c0 = 0     in actual coordinates,   since X = vd*x
%
% so the SIDE test (done at numerators, where the vertices live) uses n, while the emitted ROW
% (read by eval at actual points) uses vd*n. Using one where the other belongs rescales the offset
% and yields a PARALLEL line: measured on a triangle with half-integer vertices, whose domain came
% out as 2x + 5y >= 3 where the edge is 2x + 5y = 3/2, so the function excluded two of its own
% three vertices. Invisible on integer vertices, which is what every hand-built fixture had.
%
% The side is fixed by requiring the piece's OWN vertices to satisfy it, so no orientation
% convention is consulted and none can be read backwards.
    rows = zeros(0,3);
    for j = find(any(obj.F == k, 2)).'
        p = Vi(obj.E(j,1), :);  q = Vi(obj.E(j,2), :);
        n  = [q(2)-p(2), -(q(1)-p(1))];
        c0 = ratQ.chkB(-(n * p.'), (abs(n) * abs(p).'), 'offset');
        sgn = 0;
        for i = 1:numel(vs)
            v = Vi(vs(i), :);
            t = ratQ.chkB(n(1)*v(1) + n(2)*v(2) + c0, abs(n(1)*v(1)) + abs(n(2)*v(2)) + abs(c0), 'side');
            if t ~= 0, sgn = sign(t); break, end
        end
        if sgn == 0 && nargin >= 6
            % Every vertex is ON this edge's line. A RECESSION direction still distinguishes the
            % sides: the constraint must survive travelling to infinity along it, a condition on
            % the LINEAR part alone.
            for q = 1:size(rays,1)
                t = ratQ.chkB(n * rays(q,:).', abs(n) * abs(rays(q,:)).', 'recession side');
                if t ~= 0, sgn = sign(t); break, end
            end
        end
        if sgn == 0
            % Nothing in the piece's own geometry can say, so it is bounded BY A LINE -- a
            % half-plane. Only F knows, and n = (d2,-d1) is the RIGHT normal, so the inward normal
            % is +n when the face sits on the right and -n when on the left.
            if obj.F(j,2) == k
                sgn = 1;
            elseif obj.F(j,1) == k
                sgn = -1;
            else
                continue                            % genuinely degenerate: skip the edge
            end
        end
        rows(end+1,:) = sgn * [ratQ.chk(vd*n(1), 'edge normal'), ...
                               ratQ.chk(vd*n(2), 'edge normal'), c0]; %#ok<AGROW>
    end
end

function s = sgnOfB(row, want)
% The H-form side, after ratQ.conic may have negated the row -- the same guard conjQ's sgnOf is,
% duplicated because MATLAB gives a function file no way to share a local with another file.
    r = [0 0 0, row(1), row(2), row(3)];
    nz = find(r ~= 0, 1);
    if isempty(nz)
        error('ratQ:zeroConic', 'a constraint with all-zero coefficients names no half-plane.');
    end
    if r(nz) < 0, s = -want; else, s = want; end
end

function h = quaPolAsQuaCon(obj)
% objective: a convex QuaPol re-expressed as a QuaCon, its faces unchanged.
% co f = f when f is convex, so this is the whole computation for that branch.
    [Vi, vd] = ratQ.combineDen(obj.VN, obj.VD);
    cells = struct('num', {}, 'den', {}, 'con', {});
    for k = 1:size(obj.fN,1)
        [vs, rays, ~] = pieceGeometryB(obj, k);
        % An UNBOUNDED piece is re-expressible after all: co f = f when f is convex, whatever the
        % shape of its domain. All an unbounded piece changes is how each edge's SIDE is decided --
        % a recession direction must satisfy the constraint too, and where every vertex and every
        % ray lies ON the edge's line (a half-plane piece) only F knows. Same reasoning as conjQ's
        % pieceShape. Before this, 11 of the 13 fixtures the legacy biconj answers and this one did
        % not were refused here.
        rows = polygonHalfPlanes(obj, Vi, vd, vs, k, rays);
        % A convex f is returned UNCHANGED, so its pieces must be re-expressed exactly -- and a
        % face described by half-planes IS convex. A non-convex face cannot be, so re-expressing it
        % would silently replace the region by something else. Detected by comparing against the
        % face's own convex hull: for a convex face the two descriptions have the same rows.
        % TESTED DIRECTLY, not by comparing against the vertices' convex hull. That comparison is
        % meaningless for an UNBOUNDED piece -- a cone has one vertex and no hull edges at all, so
        % every unbounded piece looked non-convex and QuaPol.oneNorm, which IS convex, was refused.
        % The direct test is conjQ's: every vertex satisfies every half-plane, and every recession
        % direction satisfies their linear parts.
        for r = 1:size(rows,1)
            for i = 1:numel(vs)
                v = Vi(vs(i), :);
                if ratQ.chkB(rows(r,1)*v(1) + rows(r,2)*v(2) + vd*rows(r,3), abs(rows(r,1)*v(1)) + abs(rows(r,2)*v(2)) + abs(vd*rows(r,3)), 'convexity') < 0
                    error('PLQ:biconjQ:nonConvexPiece', ...
                        ['piece %d is NOT CONVEX, and co f = f here, so the answer must carry ' ...
                         'that same non-convex region -- which an intersection of half-planes ' ...
                         'cannot express. Unlike the conjugate this cannot be fixed by ' ...
                         'splitting: the convex hull of a union is not the union of hulls.'], k);
                end
            end
            for q = 1:size(rays,1)
                if ratQ.chkB(rows(r,1)*rays(q,1) + rows(r,2)*rays(q,2), abs(rows(r,1)*rays(q,1)) + abs(rows(r,2)*rays(q,2)), 'convexity') < 0
                    error('PLQ:biconjQ:nonConvexPiece', ...
                        'piece %d recedes along (%d,%d), leaving one of its own half-planes.', ...
                        k, rays(q,1), rays(q,2));
                end
            end
        end
        con = zeros(size(rows,1), 7);
        for r = 1:size(rows,1)
            con(r,:) = [ratQ.conic([0 0 0, rows(r,:)]), sgnOfB(rows(r,:), +1)];
        end
        cells(end+1) = struct('num', obj.fN(k,:), 'den', obj.fD(k), 'con', con); %#ok<AGROW>
    end
    h = assembleQuaConCells(cells);
end

function tf = isChainMesh(obj)
% objective: is this mesh a CHAIN at all -- two or more edges and no faces -- collinear or not?
    tf = ~isempty(obj.F) && ~any(obj.F(:) ~= 0) && ~isempty(obj.E) && size(obj.E,1) >= 2 ...
         && obj.nv >= 3;
end

% ------------------------------------------------------------------------------------------------
function tf = isChain(obj)
% objective: is this mesh a CHAIN -- a domain of dimension < 2 made of two or more edges whose
%   convex hull is 2-dimensional?
%
% Collinear chains are excluded deliberately: their hull is still 1-dimensional, so the answer is a
% THIN one and belongs to thinEnvelope's machinery (an H-form equality side), not to a lower hull
% that assumes a 2-D domain.
    tf = false;
    % NO FACES is the defining property and must be tested FIRST. Without it every ordinary 2-D
    % mesh with two edges and three vertices answers yes, and concaveEnvelope then reads its
    % geometry from edges instead of faces -- which is how leaving this out broke eleven tests at
    % once rather than none.
    if isempty(obj.F) || any(obj.F(:) ~= 0), return, end
    if isempty(obj.E) || size(obj.E,1) < 2, return, end
    if obj.nv < 3, return, end
    tf = ~QuaPol.isCollinear(obj.V, sqrt(eps));
end

% ------------------------------------------------------------------------------------------------
function assertChainEdgeConcave(obj)
% objective: refuse a chain edge along which q curves UPWARD, by the same rule and for the same
%   reason as a face edge.
%
% The lifted graph over one edge is a parabolic ARC. Curving downward it lies above its own chord,
% so it contributes only its two ENDPOINTS to the lower hull and the vertex values decide the
% envelope. Curving upward it bulges BELOW the chord, and then the hull needs the arc itself.
    [Vi, ~] = ratQ.combineDen(obj.VN, obj.VD);
    for j = 1:size(obj.E,1)
        [fN, ~] = obj.faceQ(j);
        Hn = [fN(5) fN(6); fN(6) fN(7)];
        d = ratQ.chk(Vi(obj.E(j,2),:) - Vi(obj.E(j,1),:), 'chain edge direction').';
        if ratQ.chkB(d.' * Hn * d, abs(d).' * abs(Hn) * abs(d), 'chain edge curvature') > 0
            error('PLQ:biconjQ:notImplemented', ...
                ['chain edge %d curves UPWARD along itself, so its lifted arc dips below its ' ...
                 'own chord and the envelope is not determined by the vertex values.'], j);
        end
    end
end

% ------------------------------------------------------------------------------------------------
function [vs, rays, bounded] = pieceGeometryB(obj, k)
% face k's true vertices and recession directions -- see conjQ's pieceGeometry, which this mirrors.
    own = find(any(obj.F == k, 2));
    if isempty(own)
        error('PLQ:biconjQ:noFace', 'face %d has no edges.', k);
    end
    E = obj.E(own, :);
    [Vi, ~] = ratQ.combineDen(obj.VN, obj.VD);
    vs = [];  rays = zeros(0,2);
    for j = 1:size(E,1)
        if E(j,3) ~= 0
            vs = [vs; E(j,1); E(j,2)]; %#ok<AGROW>
        else
            vs = [vs; E(j,1)]; %#ok<AGROW>
            rays(end+1,:) = ratQ.chk(Vi(E(j,2),:) - Vi(E(j,1),:), 'ray direction'); %#ok<AGROW>
        end
    end
    vs = unique(vs);
    bounded = isempty(rays);
end

% ------------------------------------------------------------------------------------------------
function h = collinearChainEnvelope(obj)
% objective: co f when the domain is a COLLINEAR chain -- several segments on one line, each
%   carrying its own quadratic.
%
% The hull of such a domain is still a SEGMENT, so the answer is thin and the problem is
% one-dimensional in the line's parameter: co f is the 1-D lower convex hull of the lifted
% vertices, which is piecewise AFFINE. Every edge has already been checked edge-concave, so each
% sub-segment's lifted arc lies above its own chord and contributes only its endpoints -- the same
% argument the 2-D case makes, one dimension down.
%
% EVERY CELL CARRIES THE LINE AS AN EQUALITY, exactly as thinEnvelope's single-segment case does (a
% side of 0 means ON the curve), plus two inequalities cutting it to that hull piece's own parameter
% range. So a k-piece hull is k cells, and thinEnvelope's segment is the k = 1 case.
    [Vi, vd] = ratQ.combineDen(obj.VN, obj.VD);
    ne = size(obj.E,1);

    % ---- the line, and a parameter along it ----------------------------------------------------
    % Any two distinct vertices fix the direction; the chain is collinear, so the rest lie on it.
    % The parameter t(x) = <x - A, d> is taken at NUMERATOR coordinates, so it is an integer and
    % monotone along the line.
    A = Vi(1,:).';
    d = [];
    for i = 2:size(Vi,1)
        cand = ratQ.chk(Vi(i,:).' - A, 'chain direction');
        if any(cand ~= 0), d = cand; break, end
    end
    if isempty(d)
        error('PLQ:biconjQ:degenerate', 'every vertex of the chain is the same point.');
    end
    g = gcd(abs(d(1)), abs(d(2)));  if g > 1, d = d / g; end

    % ---- the lifted vertices, over ONE denominator ---------------------------------------------
    % A vertex shared by two edges takes the SMALLER of the two values, for the reason
    % concaveEnvelope gives: the envelope is of cl f, and the lower hull is what the smaller lifted
    % point contributes.
    Lden = 1;
    for j = 1:ne, Lden = lcm(Lden, obj.fD(j)); end
    Dl = ratQ.chk(2 * Lden * vd^2, 'lift denominator');

    lifted = containers.Map('KeyType', 'double', 'ValueType', 'double');
    for j = 1:ne
        [fN, fD] = obj.faceQ(j);
        Hn = [fN(5) fN(6); fN(6) fN(7)];
        Ln = [fN(8); fN(9)];  kn = fN(10);
        scale = ratQ.chk(Lden / fD, 'denominator scaling');
        for e = 1:2
            iv = obj.E(j,e);
            v  = Vi(iv,:).';
            qv = ratQ.chkB(scale * (v.'*Hn*v + 2*vd*(Ln.'*v) + 2*vd^2*kn), abs(scale) * (abs(v).'*abs(Hn)*abs(v) + 2*abs(vd)*(abs(Ln).'*abs(v)) + 2*vd^2*abs(kn)), 'vertex value');
            if ~isKey(lifted, iv) || qv < lifted(iv), lifted(iv) = qv; end
        end
    end

    vs = cell2mat(keys(lifted)).';
    m  = numel(vs);
    T  = zeros(m,1);  Z = zeros(m,1);
    for i = 1:m
        T(i) = ratQ.chkB((Vi(vs(i),:).' - A).' * d, (abs(Vi(vs(i),:)).' + abs(A)).' * abs(d), 'chain parameter');
        Z(i) = lifted(vs(i));
    end
    [T, ord] = sort(T);  Z = Z(ord);

    % ---- the 1-D LOWER hull of (T, Z), by integer cross products --------------------------------
    % Monotone chain: walk in increasing T and pop while the last two and the new point do not turn
    % upward. Every predicate is an integer sign, so the hull is DECIDED.
    st = zeros(0,1);
    for i = 1:m
        while numel(st) >= 2
            o = st(end-1);  a1 = st(end);
            cr = ratQ.chkB((T(a1)-T(o))*(Z(i)-Z(o)) - (Z(a1)-Z(o))*(T(i)-T(o)), abs(T(a1)-T(o))*abs(Z(i)-Z(o)) + abs(Z(a1)-Z(o))*abs(T(i)-T(o)), 'hull turn');
            if cr <= 0, st(end) = []; else, break, end
        end
        st(end+1) = i; %#ok<AGROW>
    end
    if numel(st) < 2
        error('PLQ:biconjQ:degenerate', 'the collinear chain produced fewer than two hull points.');
    end

    % ---- one cell per hull piece ----------------------------------------------------------------
    n0 = [d(2); -d(1)];                                    % normal to the line
    c0 = ratQ.chkB(-(n0.' * A), (abs(n0).' * abs(A)), 'line offset');
    onLine = [ratQ.chk(vd*n0(1),'l1'), ratQ.chk(vd*n0(2),'l2'), c0];
    Ad = ratQ.chkB(A.' * d, abs(A).' * abs(d), 'anchor parameter');

    cells = struct('num', {}, 'den', {}, 'con', {});
    for k = 1:numel(st)-1
        p = st(k);  q = st(k+1);
        DT = ratQ.chk(T(q) - T(p), 'hull run');
        DZ = ratQ.chk(Z(q) - Z(p), 'hull rise');
        % Z(x)*DT = DZ*vd*<d,x> + (Z_p*DT - DZ*(A.d + T_p)),  over Dl*DT.
        num = [0 0 0 0, 0, 0, 0, ...
               ratQ.chk(DZ*vd*d(1), 'c8'), ratQ.chk(DZ*vd*d(2), 'c9'), ...
               ratQ.chkB(Z(p)*DT - DZ*(Ad + T(p)), abs(Z(p)*DT) + abs(DZ)*(abs(Ad) + abs(T(p))), 'c10')];
        den = ratQ.chk(Dl * DT, 'cell denominator');
        % t(x) >= T_p and t(x) <= T_q, at ACTUAL coordinates (linear part scaled by vd).
        lo = [ratQ.chk(vd*d(1),'s1'),  ratQ.chk(vd*d(2),'s2'),  ratQ.chk(-(Ad + T(p)),'s0')];
        hi = [ratQ.chk(-vd*d(1),'t1'), ratQ.chk(-vd*d(2),'t2'), ratQ.chk(Ad + T(q),'t0')];
        con = [ratQ.conic([0 0 0, onLine]), 0; ...
               ratQ.conic([0 0 0, lo]), sgnOfB(lo, +1); ...
               ratQ.conic([0 0 0, hi]), sgnOfB(hi, +1)];
        cells(end+1) = struct('num', num, 'den', den, 'con', con); %#ok<AGROW>
    end
    h = assembleQuaConCells(cells);
end

function h = thinEnvelope(obj)
% objective: co f when the DOMAIN itself has dimension < 2 -- a needle or a single segment.
%
% co f is convex, so its domain is conv(dom f), and for these inputs that is the domain itself: a
% point or a segment. Both are THIN, which the H-form expresses with an EQUALITY side (a side of 0
% on a curve means ON it) -- the same machinery conj uses for a point- or line-supported conjugate.
%
%   NEEDLE. dom f is one point p, and a single point is convex, so co f = f: the value q(p) at p
%   and +infinity elsewhere. Two equalities and a constant face.
%
%   SEGMENT. The problem is ONE-dimensional in the segment parameter. Where q curves UP along the
%   segment (d'Hd >= 0) it is already convex there and co f = f; where it curves DOWN, co f is the
%   CHORD -- the affine interpolant of the two endpoint values, which is the 1-D lower hull. One
%   equality for the line, two inequalities for the ends.
%
% A CHAIN of edges is not handled here and is not thin: the convex hull of two non-collinear
% segments is a two-dimensional polygon, so that case belongs with the ordinary envelope, reading
% its vertices from the whole mesh rather than from a face.
    [fN, fD] = obj.faceQ(1);
    Hn = [fN(5) fN(6); fN(6) fN(7)];
    Ln = [fN(8); fN(9)];
    kn = fN(10);
    [Vi, vd] = ratQ.combineDen(obj.VN, obj.VD);

    if isempty(obj.E)
        % ---- a NEEDLE ------------------------------------------------------------------------
        p  = Vi(1,:).';
        qp = ratQ.chkB(p.'*Hn*p + 2*vd*(Ln.'*p) + 2*vd^2*kn, abs(p).'*abs(Hn)*abs(p) + 2*abs(vd)*(abs(Ln).'*abs(p)) + 2*vd^2*abs(kn), 'needle value');   % over 2*fD*vd^2
        r1 = [ratQ.chk(vd,'e1'), 0, ratQ.chk(-p(1),'e0')];
        r2 = [0, ratQ.chk(vd,'e1'), ratQ.chk(-p(2),'e0')];
        con = [ratQ.conic([0 0 0, r1]), 0; ratQ.conic([0 0 0, r2]), 0];
        h = assembleQuaConCells(struct('num', [0 0 0 0 0 0 0 0 0 qp], ...
                                       'den', ratQ.chk(2*fD*vd^2,'needle denominator'), ...
                                       'con', con));
        return
    end

    if size(obj.E,1) ~= 1
        % UNREACHABLE, and kept as an assertion rather than deleted. Every chain -- collinear or
        % not -- is routed before thinEnvelope is called: a non-collinear one to concaveEnvelope's
        % 2-D hull, a collinear one to collinearChainEnvelope's 1-D hull along its own span. So a
        % multi-edge mesh arriving HERE means the routing above let one through.
        error('PLQ:biconjQ:chainNotImplemented', ...
            ['a chain of %d edges reached thinEnvelope, which handles only a single segment. ' ...
             'Chains are routed before this point, so this is a defect in that routing.'], ...
            size(obj.E,1));
    end
    if obj.E(1,3) == 0
        error('PLQ:biconjQ:unbounded', ...
            ['the domain is a RAY. Where q curves down along it the envelope runs to -infinity, ' ...
             'which is a correct answer with nowhere to be stored -- the gap conjCPLQ records as ' ...
             'convEnvUnbounded:minusInfinity.']);
    end

    % ---- a SEGMENT --------------------------------------------------------------------------
    a = Vi(obj.E(1,1),:).';  b = Vi(obj.E(1,2),:).';
    d = ratQ.chk(b - a, 'segment direction');
    alpha = ratQ.chkB(d.' * Hn * d, abs(d).' * abs(Hn) * abs(d), 'segment curvature');

    n  = [d(2); -d(1)];
    c0 = ratQ.chkB(-(n.' * a), (abs(n).' * abs(a)), 'line offset');
    rows = [ratQ.chk(vd*n(1),'l1'), ratQ.chk(vd*n(2),'l2'), c0];          % ON the line
    lo   = [ratQ.chk(vd*d(1),'s1'), ratQ.chk(vd*d(2),'s2'), ratQ.chkB(-(d.'*a), abs(d).'*abs(a), 's0')];
    hi   = [ratQ.chk(-vd*d(1),'t1'), ratQ.chk(-vd*d(2),'t2'), ratQ.chkB(d.'*b, abs(d).'*abs(b), 't0')];
    con  = [ratQ.conic([0 0 0, rows]), 0; ...
            ratQ.conic([0 0 0, lo]), sgnOfB(lo, +1); ...
            ratQ.conic([0 0 0, hi]), sgnOfB(hi, +1)];

    if alpha >= 0
        % q is convex along the segment, so it is already its own envelope there.
        num = fN;  den = fD;
    else
        % q curves DOWN along the segment: the envelope is the CHORD, the affine interpolant of the
        % two endpoint values. Written over 2*fD*vd^2*(d.d), and checked to reproduce q(a) and q(b)
        % exactly at the ends.
        qa = ratQ.chkB(a.'*Hn*a + 2*vd*(Ln.'*a) + 2*vd^2*kn, abs(a).'*abs(Hn)*abs(a) + 2*abs(vd)*(abs(Ln).'*abs(a)) + 2*vd^2*abs(kn), 'endpoint value');
        qb = ratQ.chkB(b.'*Hn*b + 2*vd*(Ln.'*b) + 2*vd^2*kn, abs(b).'*abs(Hn)*abs(b) + 2*abs(vd)*(abs(Ln).'*abs(b)) + 2*vd^2*abs(kn), 'endpoint value');
        dd = ratQ.chk(d.'*d, 'segment length squared');
        num = [0 0 0 0, 0, 0, 0, ...
               ratQ.chk((qb-qa)*vd*d(1), 'c8'), ratQ.chk((qb-qa)*vd*d(2), 'c9'), ...
               ratQ.chkB(qa*dd - (qb-qa)*(a.'*d), abs(qa*dd) + abs(qb-qa)*(abs(a).'*abs(d)), 'c10')];
        den = ratQ.chk(2*fD*vd^2*dd, 'chord denominator');
    end
    h = assembleQuaConCells(struct('num', num, 'den', den, 'con', con));
end
