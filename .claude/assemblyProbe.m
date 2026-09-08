function assemblyProbe(varargin)
% THE ASSEMBLY: walk the cell pairs, build every face, and stitch the bridge and non-bridge faces
% into ONE mesh whose evaluation is the biconjugate.
%
% WHAT THIS JOINS. Every ingredient was built and verified separately:
%
%   conjQ per piece                     ->  gA, gB   (single-piece input, so all cell edges are lines)
%   FILTER 1, cells meet in 2-D         ->  ratQ.feasible2(..., strict)
%   FILTER 2, the tie changes sign      ->  ratQ.quadSignOnPoly       (.claude/cellPairFilterProbe.m)
%   the NON-BRIDGE face over a cell     ->  .claude/cellFaceProbe.m, .claude/cellFaceExactProbe.m
%   the BRIDGE region and its sides     ->  .claude/bridgeRegionBuild.m, .claude/bridgeRegionSides.m
%   the bridge pipeline end to end      ->  .claude/bridgeEnvelopeProbe.m
%
% This is the join, and the join is a PARTITION plus a MAX:
%
%       f**(x) = sup_s <s,x> - G(s),   G = max(gA, gB)
%              = MAX over pieces R of   sup_{s in R} <s,x> - G|_R(s) ,
%
% for any cover {R} of the s-plane on which G is a single quadratic. The cover is the cell pairs:
% over P_ij = cell_i(gA) ^ cell_j(gB), G is max(gA_i, gB_j), and the tie conic q = gA_i - gB_j
% either has one strict sign on P_ij -- one quadratic, a NON-BRIDGE piece -- or changes sign, and
% then P_ij splits into {q >= 0} carrying gA_i and {q <= 0} carrying gB_j, whose common boundary is
% the tie arc that carries the BRIDGE face.
%
% WHY A MAX OF SUPS NEEDS NO GLUING ARGUMENT, which is what makes the assembly this short. Every
% candidate this evaluator forms is <s,x> - G(s) at a point s that is FEASIBLE for its own piece,
% so no candidate can exceed f**(x): the assembly cannot overshoot, whatever the mesh looks like.
% It can only UNDERSHOOT, and only by missing the global maximiser -- so the whole correctness
% question is "is every maximiser some piece's candidate", and that is the enumeration below. The
% two are checked separately: the coverage check says the pieces cover the plane, the value check
% says the max reproduces the definition.
%
% THE CANDIDATE ENUMERATION, per piece, with the reason each kind is needed:
%
%   * the stationary point of the piece's own quadratic, when it lies in the piece
%   * per facet of P_ij, the 1-D concave quadratic on the facet -- its vertex and its endpoints --
%     restricted to the part of the facet on the piece's side of the tie conic. The restriction is
%     NOT cosmetic: a candidate on the wrong side would be evaluating gA_i where G = gB_j > gA_i,
%     which is the one way this evaluator could overshoot.
%   * on a bridge pair, the TIE ARC. At a maximiser that ties, x is in the subdifferential of a max,
%     i.e. x = (1-t) grad gA_i(s) + t grad gB_j(s) for some t in [0,1] -- exactly the sweeping
%     segment of .claude/bridgeRegionBuild.m. Eliminating s gives the QUARTIC in t of
%     .claude/bridgeRegionSides.m, so the arc candidates are its roots in [0,1]. Nothing new is
%     derived here; the bridge machinery is called.
%
% WHAT IS EXACT AND WHAT IS NOT, said plainly. The MESH decisions are exact: which pairs meet
% (ratQ.feasible2) and which are bridges (ratQ.quadSignOnPoly) are integer predicates. The
% per-piece EVALUATION here is in doubles -- the candidate enumeration is the one
% .claude/cellFaceExactProbe.m already implements over Q, and swapping it in is a substitution, not
% new mathematics. The quartic's roots are found numerically for the same reason. So this probe
% answers "does the stitched mesh reproduce f**", not "is every stitch decided in Q".
%
%   assemblyProbe             % three fixtures, 120 sample points each
%   assemblyProbe(400)

    n = 120;
    if nargin >= 1, n = varargin{1}; end

    % ONE KNOB, THREE REGIMES. c is the coefficient of the linear form vanishing on the shared
    % diagonal, so the two pieces agree there for every c and the input is always a function:
    % c = 0 is CONVEX (and then f** = f, a second oracle that shares nothing with the first),
    % c = -1 and c = -2 are non-convex kinks of increasing depth -- the case the bridge exists for.
    %
    % A FOURTH FIXTURE THAT IS NOT THE UNIT SQUARE, because three runs on one geometry answer one
    % geometry. The rectangle splits into two triangles the same way but is anisotropic, so its
    % dual mesh is not symmetric under swapping the coordinates.
    F = { 'unit square, c = 0 (CONVEX)',  [0 0; 1 0; 1 1], [0 0; 1 1; 0 1],  0, [0 1; 0 1]; ...
          'unit square, c = -1',          [0 0; 1 0; 1 1], [0 0; 1 1; 0 1], -1, [0 1; 0 1]; ...
          'unit square, c = -2',          [0 0; 1 0; 1 1], [0 0; 1 1; 0 1], -2, [0 1; 0 1]; ...
          'rectangle 2x1, c = -2',        [0 0; 2 0; 2 1], [0 0; 2 1; 0 1], -2, [0 2; 0 1] };
    for k = 1:size(F,1)
        c = F{k,4};  box = F{k,5};
        fprintf('\n================ %s ================\n', F{k,1});
        fA = tri(F{k,2}, 0);
        fB = tri(F{k,3}, c);
        gA = conjQ(fA);   gB = conjQ(fB);
        [pieces, stats] = buildMesh(gA, gB);
        report(stats, pieces);
        coverage(pieces, gA, gB);
        domainCheck(pieces, box);
        facetHardPoints(pieces, gA, gB);
        verify(pieces, gA, gB, fA, fB, c, box, n);
    end
end

% ================================================================================================
function [pieces, stats] = buildMesh(gA, gB)
% WALK THE CELL PAIRS AND BUILD THE PIECES. Both filters are the exact ones; nothing here is a
% tolerance.
    pieces = struct('kind', {}, 'P', {}, 'g', {}, 'g2', {}, 'q', {}, 'side', {}, 'ij', {});
    stats = struct('nPairs', 0, 'nMeet', 0, 'nSame', 0, 'nNonBridge', 0, 'nBridge', 0, ...
                   'nSkipCurved', 0);
    for i = 1:gA.nf
        [RA, okA] = cellRows(gA, i);
        if ~okA, stats.nSkipCurved = stats.nSkipCurved + 1; continue, end
        if isempty(RA), continue, end
        for j = 1:gB.nf
            [RB, okB] = cellRows(gB, j);
            if ~okB || isempty(RB), continue, end
            stats.nPairs = stats.nPairs + 1;
            P = [RA; RB];
            if ~ratQ.feasible2(P, true), continue, end       % FILTER 1: meet in 2-D
            stats.nMeet = stats.nMeet + 1;

            q = tieConic(gA, i, gB, j);
            if all(q == 0)
                % IDENTICAL FACES. The two conjugates agree on the whole cell -- both primal
                % pieces share that vertex and its value -- so there is no tie curve and no bridge:
                % one quadratic over the whole piece. Counted apart because calling it a bridge is
                % the mistake .claude/bridgeEnvelopeProbe.m made on its first run.
                stats.nSame = stats.nSame + 1;
                pieces(end+1) = mkPiece('poly', P, face(gA,i), [], [0 0 0 0 0 0], 0, [i j]); %#ok<AGROW>
                continue
            end

            [sgn, strict] = ratQ.quadSignOnPoly(q, P);       % FILTER 2: does the tie change sign?
            if sgn ~= 0 && strict
                stats.nNonBridge = stats.nNonBridge + 1;
                if sgn > 0, pieces(end+1) = mkPiece('poly', P, face(gA,i), [], q, 0, [i j]); %#ok<AGROW>
                else,       pieces(end+1) = mkPiece('poly', P, face(gB,j), [], q, 0, [i j]); %#ok<AGROW>
                end
            else
                % A BRIDGE PAIR: the piece splits along the tie arc, and the arc itself carries a
                % face. Both sides are stored on one piece because the arc belongs to neither.
                stats.nBridge = stats.nBridge + 1;
                pieces(end+1) = mkPiece('bridge', P, face(gA,i), face(gB,j), q, 1, [i j]); %#ok<AGROW>
            end
        end
    end
end

function p = mkPiece(kind, P, g, g2, q, side, ij)
    p = struct('kind', kind, 'P', P, 'g', g, 'g2', g2, 'q', q, 'side', side, 'ij', ij);
end

% ================================================================================================
function [v, who, vFacet] = evalMesh(pieces, x)
% THE STITCHED MESH, EVALUATED. The max over every piece's candidates, with the winner's name.
%
% `vFacet` is the best value attained by a FACET candidate alone -- a facet vertex, a facet
% endpoint or a constant facet. It is reported because the facet candidates never come out as the
% strict winner, and an unexplained absence is worth measuring rather than asserting.
%
% THE EXPLANATION, AND WHAT MEASURING IT CORRECTED. gA and gB are each the conjugate of a STRICTLY
% convex quadratic on one triangle, so the primal maximiser is unique and each conjugate is C^1
% ACROSS ITS OWN CELL EDGES: at a point s on such an edge the subdifferential of gA is the single
% gradient. G = max(gA, gB) is therefore non-smooth only on the TIE, so a facet maximiser occurs
% only for x on the IMAGE OF A CURVE -- a measure-zero set a random draw never hits. The first
% version of this comment said a facet candidate "ties but never wins"; measured, it does not even
% tie -- 0 of 150 random points in every fixture. The tie is checked where it can actually happen,
% at constructed points (facetHardPoints below), and the facet loop's measured load-bearing role is
% a different one: it is what detects an UNBOUNDED direction, i.e. what makes dom f** come out right.
    v = -inf;  who = 'none';  vFacet = -inf;
    for k = 1:numel(pieces)
        p = pieces(k);
        if strcmp(p.kind, 'poly')
            [w, nm, wf] = polyCands(p.g, p.P, x, [0 0 0 0 0 0], 0);
            [v, who] = keepBest(v, who, w, sprintf('%s@(%d,%d)', nm, p.ij(1), p.ij(2)));
            vFacet = max(vFacet, wf);
        else
            % the {q >= 0} side carries gA_i, the {q <= 0} side carries gB_j, and the arc both
            [w, nm, wf] = polyCands(p.g,  p.P, x, p.q, +1);
            [v, who] = keepBest(v, who, w, sprintf('%s+@(%d,%d)', nm, p.ij(1), p.ij(2)));
            vFacet = max(vFacet, wf);
            [w, nm, wf] = polyCands(p.g2, p.P, x, p.q, -1);
            [v, who] = keepBest(v, who, w, sprintf('%s-@(%d,%d)', nm, p.ij(1), p.ij(2)));
            vFacet = max(vFacet, wf);
            [w, nm] = arcCands(p.g, p.g2, p.q, p.P, x);
            [v, who] = keepBest(v, who, w, sprintf('%s@(%d,%d)', nm, p.ij(1), p.ij(2)));
        end
    end
end

function [v, who] = keepBest(v, who, w, nm)
    if w > v, v = w; who = nm; end
end

% ---- the candidates over one polyhedral piece, optionally cut by the tie conic -------------------
function [best, who, bestFacet] = polyCands(g, P, x, q, want)
% sup over { s in P : want*q(s) >= 0 } of <s,x> - g(s), by candidate enumeration. want = 0 means no
% conic restriction. Every candidate returned is FEASIBLE, which is what forbids an overshoot.
% `bestFacet` is the same max taken over the FACET candidates only; see evalMesh for why.
    H = g.A;  y = x - g.b;  cc = g.k;
    best = -inf;  who = 'none';  bestFacet = -inf;

    if abs(det(H)) > 1e-12
        s = H \ y;
        if inP(P, s) && qOK(q, want, s)
            best = s.'*y - 0.5*(s.'*H*s) - cc;  who = 'stationary';
        end
    end

    for i = 1:size(P,1)
        Ni = P(i,1:2);
        if all(Ni == 0), continue, end
        d = [Ni(2); -Ni(1)];
        p = -P(i,3) * Ni(:) / (Ni * Ni.');
        [lo, hi, empty] = lineRange(P, p, d);
        if empty, continue, end

        a2 = -0.5*(d.'*H*d);  a1 = d.'*y - p.'*H*d;  a0 = p.'*y - 0.5*(p.'*H*p) - cc;
        for iv = sideIntervals(q, want, p, d, lo, hi)
            u = iv{1}(1);  w = iv{1}(2);
            if a2 < -1e-14
                ts = -a1/(2*a2);
                if ts > u && ts < w
                    val = a0 + a1*ts + a2*ts^2;
                    bestFacet = max(bestFacet, val);
                    if val > best, best = val; who = 'vertex'; end
                end
            elseif abs(a1) <= 1e-12
                bestFacet = max(bestFacet, a0);
                if a0 > best, best = a0; who = 'constant'; end
            elseif (a1 > 0 && isinf(w)) || (a1 < 0 && isinf(u))
                best = inf;  who = 'unbounded';  bestFacet = inf;  return
            end
            for t = [u w]
                if isinf(t), continue, end
                val = a0 + a1*t + a2*t^2;
                bestFacet = max(bestFacet, val);
                if val > best, best = val; who = 'endpoint'; end
            end
        end
    end
end

function IV = sideIntervals(q, want, p, d, lo, hi)
% The parts of the facet's parameter interval [lo,hi] on the requested side of the tie conic.
% q restricted to s = p + t*d is a QUADRATIC in t, so the side is at most two intervals, and its
% endpoints -- the arc crossings -- are candidates in their own right, which is why they are
% returned rather than discarded.
    if want == 0, IV = {[lo hi]}; return, end
    al =  q(1)*d(1)^2 + q(2)*d(1)*d(2) + q(3)*d(2)^2;
    be =  2*q(1)*p(1)*d(1) + q(2)*(p(1)*d(2) + p(2)*d(1)) + 2*q(3)*p(2)*d(2) + q(4)*d(1) + q(5)*d(2);
    ga =  q(1)*p(1)^2 + q(2)*p(1)*p(2) + q(3)*p(2)^2 + q(4)*p(1) + q(5)*p(2) + q(6);
    al = want*al;  be = want*be;  ga = want*ga;          % now: keep { t : al t^2 + be t + ga >= 0 }

    sc = max([abs(al), abs(be), abs(ga), 1]);
    if abs(al) < 1e-12*sc
        if abs(be) < 1e-12*sc
            if ga >= -1e-12*sc, IV = {[lo hi]}; else, IV = {}; end
            return
        end
        r = -ga/be;
        if be > 0, IV = clipOne(max(lo, r), hi); else, IV = clipOne(lo, min(hi, r)); end
        return
    end
    disc = be^2 - 4*al*ga;
    if disc < 0
        if al > 0, IV = {[lo hi]}; else, IV = {}; end
        return
    end
    rr = sort([(-be - sqrt(disc))/(2*al), (-be + sqrt(disc))/(2*al)]);
    if al > 0
        IV = [clipOne(lo, min(hi, rr(1))), clipOne(max(lo, rr(2)), hi)];
    else
        IV = clipOne(max(lo, rr(1)), min(hi, rr(2)));
    end
end

function IV = clipOne(u, w)
    if u > w, IV = {}; else, IV = {[u w]}; end
end

% ---- the candidates ON the tie arc: the bridge face ---------------------------------------------
function [best, who] = arcCands(g1, g2, q, P, x)
% A maximiser that TIES has x in the subdifferential of the max there, i.e.
%       x = (1-t) grad g1(s) + t grad g2(s) = M_t s + v_t,   t in [0,1],
% so s = M_t^-1 (x - v_t) and q(s) = 0 becomes a QUARTIC in t once det(M_t) is cleared -- the Phi
% of .claude/bridgeRegionSides.m, formed here in doubles and rooted numerically.
    best = -inf;  who = 'none';
    A1 = g1.A;  A2 = g2.A;  b1 = g1.b;  b2 = g2.b;

    m11 = [A2(1,1)-A1(1,1), A1(1,1)];   m12 = [A2(1,2)-A1(1,2), A1(1,2)];
    m21 = [A2(2,1)-A1(2,1), A1(2,1)];   m22 = [A2(2,2)-A1(2,2), A1(2,2)];
    adj11 = m22;  adj12 = -m12;  adj21 = -m21;  adj22 = m11;
    D = padd(conv(m11, m22), -conv(m12, m21));
    w1 = padd([0, x(1)], -[b2(1)-b1(1), b1(1)]);
    w2 = padd([0, x(2)], -[b2(2)-b1(2), b1(2)]);
    u1 = padd(conv(adj11, w1), conv(adj12, w2));
    u2 = padd(conv(adj21, w1), conv(adj22, w2));
    Phi = padd(padd(padd(q(1)*conv(u1,u1), q(2)*conv(u1,u2)), q(3)*conv(u2,u2)), ...
               padd(conv(D, padd(q(4)*u1, q(5)*u2)), q(6)*conv(D, D)));

    sc = max(abs(Phi));
    if sc == 0 || ~all(isfinite(Phi)), return, end
    Phi = Phi / sc;
    r = roots(Phi);
    for t = r(:).'
        if abs(imag(t)) > 1e-8, continue, end
        t = real(t);
        if t < -1e-9 || t > 1 + 1e-9, continue, end
        t = min(max(t, 0), 1);
        Mt = (1-t)*A1 + t*A2;   vt = (1-t)*b1 + t*b2;
        if abs(det(Mt)) < 1e-12, continue, end            % the segment degenerates; not a witness
        s = Mt \ (x - vt);
        if ~inP(P, s), continue, end
        if abs(evalConic(q, s)) > 1e-6 * max(1, s.'*s), continue, end   % must really be a tie
        val = s.'*x - gval(g1, s);
        if val > best, best = val; who = 'arc'; end
    end
end

% ================================================================================================
function report(stats, pieces)
    nB = sum(strcmp({pieces.kind}, 'bridge'));
    fprintf('cell pairs %d | meet in 2-D %d | identical faces %d | non-bridge %d | BRIDGE %d\n', ...
            stats.nPairs, stats.nMeet, stats.nSame, stats.nNonBridge, stats.nBridge);
    fprintf('mesh: %d pieces (%d polyhedral, %d bridge pairs carrying 3 faces each)\n', ...
            numel(pieces), numel(pieces) - nB, nB);
    if stats.nSkipCurved > 0
        fprintf('  NOTE: %d cell(s) skipped -- a curved cell edge, which this assembly does not cover\n', ...
                stats.nSkipCurved);
    end
end

function coverage(pieces, gA, gB)
% DO THE PIECES COVER THE s-PLANE? A hole in the cover is the assembly's only way to be wrong that
% the value check can miss -- it would show up only at the x whose maximiser falls in the hole. So
% it is asked directly, and asked about the CELL PAIR a random s actually belongs to.
    rng(7);
    nHole = 0;  nOK = 0;
    for t = 1:400
        s = (rand(2,1)*8 - 4);
        [~, iA] = gA.eval(s.');   [~, jB] = gB.eval(s.');
        if iA == 0 || jB == 0, continue, end
        hit = false;
        for k = 1:numel(pieces)
            if inP(pieces(k).P, s), hit = true; break, end
        end
        if hit, nOK = nOK + 1; else, nHole = nHole + 1; end
    end
    fprintf('coverage: %d random s in some piece, %d in NONE\n', nOK, nHole);
end

function domainCheck(pieces, box)
% IS THE MESH FINITE EXACTLY ON dom f**? The value check above only samples INSIDE the domain, so
% on its own it says nothing about the boundary of the answer -- and "the mesh is finite where the
% envelope is infinite" is a face too many, while the reverse is a face missing.
%
% dom f** = conv(dom f) = the primal box, so this is a decidable statement with an independent
% right answer. Points within a margin of the box's boundary are SKIPPED rather than counted:
% there the sup is finite but unbounded-in-the-limit and a double evaluation legitimately goes
% either way, so counting them would measure the tolerance rather than the mesh.
    marg = 0.02;
    nIn = 0;  nOut = 0;  badIn = 0;  badOut = 0;
    for a = box(1,1)-1 : 0.1 : box(1,2)+1
        for b = box(2,1)-1 : 0.1 : box(2,2)+1
            x = [a; b];
            d = max([box(:,1) - x; x - box(:,2)]);   % signed distance to the box, in the inf-norm
            if abs(d) < marg, continue, end
            v = evalMesh(pieces, x);
            if d < 0
                nIn = nIn + 1;
                if ~isfinite(v), badIn = badIn + 1; end
            else
                nOut = nOut + 1;
                if isfinite(v), badOut = badOut + 1; end
            end
        end
    end
    fprintf(['domain: %d points inside dom f** (%d wrongly INFINITE), ' ...
             '%d outside (%d wrongly FINITE)\n'], nIn, badIn, nOut, badOut);
end

function facetHardPoints(pieces, gA, gB)
% THE POINTS A RANDOM DRAW NEVER REACHES: x whose maximiser sits ON A CELL EDGE, so the winning
% candidate has to be a facet one. Constructed, not sampled -- the same reason
% .claude/cellFaceExactProbe.m constructs its hard points.
%
% THE CONSTRUCTION. Take s0 in the relative interior of a facet of a piece, and set x = grad(s0)
% for that piece's own face. Then s0 is a maximiser of <s,x> - G(s) provided the piece's face IS
% the max at s0, which is checked rather than assumed. The mesh must return that value AND a facet
% candidate must attain it -- if only the value matched, the mesh could be reading it off some
% other piece and the facet loop would still be dead code.
    nH = 0;  nVal = 0;  nFacet = 0;
    for k = 1:numel(pieces)
        p = pieces(k);
        if ~strcmp(p.kind, 'poly'), continue, end
        H = p.g.A;  b = p.g.b;
        for i = 1:size(p.P,1)
            Ni = p.P(i,1:2);
            if all(Ni == 0), continue, end
            d = [Ni(2); -Ni(1)];
            pp = -p.P(i,3) * Ni(:) / (Ni * Ni.');
            [lo, hi, empty] = lineRange(p.P, pp, d);
            if empty, continue, end
            if isinf(lo) && isinf(hi), t0 = 0;
            elseif isinf(lo),          t0 = hi - 1;
            elseif isinf(hi),          t0 = lo + 1;
            else,                      t0 = (lo + hi) / 2;
            end
            s0 = pp + t0*d;
            if ~inP(p.P, s0), continue, end
            % the piece's own face must be the max at s0, or x = grad(s0) names a different problem
            gv = gval(p.g, s0);
            if abs(max(gA.eval(s0.'), gB.eval(s0.')) - gv) > 1e-9 * max(1, abs(gv)), continue, end
            x = H*s0 + b;
            want = s0.'*x - gv;
            [v, ~, vF] = evalMesh(pieces, x);
            nH = nH + 1;
            if isfinite(v) && abs(v - want) <= 1e-8 * max(1, abs(want)), nVal = nVal + 1; end
            if isfinite(vF) && abs(vF - want) <= 1e-8 * max(1, abs(want)), nFacet = nFacet + 1; end
        end
    end
    fprintf(['facets: %d CONSTRUCTED points whose maximiser is on a cell edge -- %d right, ' ...
             '%d with a facet candidate attaining it\n'], nH, nVal, nFacet);
end

function verify(pieces, gA, gB, fA, fB, c, box, n)
% THE VALUE CHECK, against two oracles that share nothing with the assembly: the definition's own
% sup, and -- when the input is convex -- the primal function itself, since f** = f there.
%
% THE SAMPLE IS THE PRIMAL DOMAIN, AND THAT IS NOT A CONVENIENCE. dom f** = conv(dom f) = the
% primal box, so f** is +infinity everywhere outside it and a point drawn outside measures nothing.
% A first version of this probe drew x from [-1.5, 2.5]^2 and reported a worst gap of 1.1e-16 --
% from TWO surviving points out of 40, the other 38 skipped as infinite on both sides. The tell was
% the one MORNING.md names: the same worst-case point in all three fixtures is measuring the
% harness. The domain itself is checked separately, above.
    rng(2026);
    X = rand(n,2) .* (box(:,2) - box(:,1)).' + box(:,1).';
    worst = 0;  worstAt = [0 0];  nOver = 0;  nUnder = 0;  nSkip = 0;  kinds = {};  nFacetTie = 0;
    for t = 1:n
        x = X(t,:).';
        [vA, who, vF] = evalMesh(pieces, x);
        if isfinite(vA) && isfinite(vF) && abs(vF - vA) <= 1e-9 * max(1, abs(vA))
            nFacetTie = nFacetTie + 1;
        end
        vD = envelopeByDefinition(gA, gB, x);
        if ~isfinite(vA) || ~isfinite(vD), nSkip = nSkip + 1; continue, end
        rel = (vA - vD) / max(1, abs(vD));
        if rel >  1e-6, nOver  = nOver  + 1; end
        if rel < -1e-6, nUnder = nUnder + 1; end
        if abs(rel) > worst, worst = abs(rel); worstAt = x.'; end
        kinds{end+1} = regexprep(who, '[-+]?@.*$', ''); %#ok<AGROW>
    end
    fprintf('value: %d of %d sample points measured (%d infinite on one side)\n', ...
            n - nSkip, n, nSkip);
    fprintf('       worst relative gap to the definition %.3e at (%.3f, %.3f)\n', ...
            worst, worstAt(1), worstAt(2));
    fprintf('       %d points where the mesh OVERSHOT, %d where it UNDERSHOT (tol 1e-6)\n', ...
            nOver, nUnder);
    u = unique(kinds);
    fprintf('       winning face kinds:');
    for i = 1:numel(u), fprintf('  %s x%d', u{i}, sum(strcmp(kinds, u{i}))); end
    fprintf('\n');
    fprintf('       a FACET candidate attains the winning value at %d of %d random points\n', ...
            nFacetTie, n - nSkip);

    if c == 0
        % f is CONVEX here, so f** = f ON THE PRIMAL DOMAIN. A second oracle with no fminsearch in
        % it: if the assembly were reproducing the definition's *search* rather than the envelope,
        % this is what would catch it.
        rng(5);  worstP = 0;  nP = 0;
        for t = 1:200
            x = rand(2,1) .* (box(:,2) - box(:,1)) + box(:,1);
            [v1, r1] = fA.eval(x.');   [v2, r2] = fB.eval(x.');
            if r1 > 0, fv = v1; elseif r2 > 0, fv = v2; else, continue, end
            va = evalMesh(pieces, x);
            if ~isfinite(va), continue, end
            nP = nP + 1;
            worstP = max(worstP, abs(va - fv) / max(1, abs(fv)));
        end
        fprintf('       convex fixture: f** vs f itself on %d primal points, worst %.3e\n', nP, worstP);
    end
end

% ================================================================================================
% ---- oracle -------------------------------------------------------------------------------------
function v = envelopeByDefinition(gA, gB, x)
% f**(x) = sup_s <s,x> - max(gA(s), gB(s)), by direct search. Shares no arithmetic with the mesh.
    h = @(s) -( s(:).' * x - max(gA.eval(s(:).'), gB.eval(s(:).')) );
    best = inf;
    for s0 = [x, [0;0], -x, 2*x]
        s = fminsearch(h, s0, optimset('Display', 'off', 'TolX', 1e-12, 'TolFun', 1e-12));
        if h(s) < best, best = h(s); end
    end
    v = -best;
end

% ---- small helpers ------------------------------------------------------------------------------
function tf = inP(P, s), tf = all(P(:,1:2)*s + P(:,3) >= -1e-9); end

function tf = qOK(q, want, s)
    if want == 0, tf = true; return, end
    tf = want * evalConic(q, s) >= -1e-9;
end

function v = evalConic(q, s)
    v = q(1)*s(1)^2 + q(2)*s(1)*s(2) + q(3)*s(2)^2 + q(4)*s(1) + q(5)*s(2) + q(6);
end

function v = gval(g, s), v = 0.5*(s.'*g.A*s) + g.b.'*s + g.k; end

function [lo, hi, empty] = lineRange(P, p, d)
    lo = -inf;  hi = inf;  empty = false;
    for j = 1:size(P,1)
        A = P(j,1:2)*d;  B = P(j,1:2)*p + P(j,3);
        if abs(A) < 1e-14
            if B < -1e-9, empty = true; return, end
        elseif A > 0, lo = max(lo, -B/A);
        else, hi = min(hi, -B/A);
        end
    end
    if lo > hi + 1e-9, empty = true; end
end

function c = padd(a, b)
    n = max(numel(a), numel(b));
    c = [zeros(1, n-numel(a)), a] + [zeros(1, n-numel(b)), b];
end

function q = tieConic(gA, i, gB, j)
% The locus where the two cells' conjugates agree, cleared to integers. The weighted basis carries
% halves on the squared terms, so the conic is [n5, 2*n6, n7, 2*n8, 2*n9, 2*n10] -- getting this
% wrong named a different curve and moved real numbers (DECISIONS 2026-09-06).
    D = gA.fN(i,:) * gB.fD(j) - gB.fN(j,:) * gA.fD(i);
    q = [D(5), 2*D(6), D(7), 2*D(8), 2*D(9), 2*D(10)];
end

function g = face(gc, k)
% Cell k's face as g(s) = s'A s/2 + b's + k, which is what the 10-wide weighted basis stores.
    n = gc.fN(k,:);  d = gc.fD(k);
    g = struct('A', [n(5), n(6); n(6), n(7)] / d, 'b', [n(8); n(9)] / d, 'k', n(10) / d);
end

function [rows, ok] = cellRows(g, k)
% The cell's half-planes as integer rows. `ok` is false when an edge is a genuine CONIC -- this
% assembly covers polyhedral cells, and a curved cell edge is reported rather than dropped in
% silence.
    rows = zeros(0,3);  ok = true;
    FC = g.FC{k};
    for t = 1:size(FC,1)
        cc = g.EcQ(FC(t,1),:);
        if any(cc(1:3) ~= 0), rows = zeros(0,3); ok = false; return, end
        if FC(t,2) == 0, continue, end
        rows(end+1,:) = FC(t,2) * cc(4:6); %#ok<AGROW>
    end
end

function f = tri(V, c)
% A positive definite quadratic on one triangle, perturbed by c times the linear form vanishing on
% the shared diagonal, so the two pieces agree there for every c -- DECISIONS 2026-09-05.
    q = [0 0 0 0, 2, 0, 2, 0, 0, 0];
    q(8) = q(8) - c;
    q(9) = q(9) + c;
    f = QuaPol(V, [1 2 1; 2 3 1; 3 1 1], q, [1 0; 1 0; 1 0]);
end
