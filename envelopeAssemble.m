function [env, info] = envelopeAssemble(gA, gB)
% envelopeAssemble  The convex envelope of a two-piece function, assembled from the pieces'
%   conjugates: the third step, and the one that turns a pile of verified parts into a function.
%
% objective: `f** = sup_s [<s,x> - G(s)]` with `G = max(gA, gB)`, and that sup splits over the
%   overlay of the two conjugates' cells. Every ingredient has been built and verified separately;
%   this walks the pairs and folds them.
%
% [input]  gA, gB : QuaCon, the exact conjugates of the two pieces (conjQ of each)
% [output] env  : the envelope as its FACES -- struct with
%                   .quad   struct array of quadratic candidates: .num, .den (the face function),
%                           .con (its region, as conic rows), .pair (which cell pair it came from)
%                   .bridge cell array of the CURVED-tie faces, as `bridgeFaceAlg` returns them
%                 `envelopeEval` evaluates it; the value at x is the MAX over the faces whose
%                 region contains x.
%          info : counts per pair class
%
% ==================================================================================================
% ONE RULE PER PAIR, AND THE CENSUS SAYS WHICH RULE APPLIES
%
% Over the pair `P_ij` the max is `max(gA_i, gB_j)` and the tie is `q = gA_i - gB_j`:
%
%   NO SIGN CHANGE  -- one quadratic wins on the whole pair. Its candidates over `P_ij` are
%                      `envelopeCellRegions`' three kinds, and that is the pair's contribution.
%
%   TANGENTIAL      -- the tie TOUCHES the pair without crossing it. One quadratic still wins
%                      everywhere (with equality on the touching set), so this is the case above
%                      with the winner read off the non-strict sign. Measured, 10 of 29 bridge
%                      pairs on the standard fixtures (DECISIONS.md 2026-09-17).
%
%   STRAIGHT TIE    -- the tie is a LINE, a double line or a parallel pair. Each side of it is then
%                      a POLYGON, so the pair splits into two polyhedral cells and the ordinary
%                      candidates cover both. There is NO separate bridge face: a straight tie is
%                      an ordinary facet. 14 of 29.
%
%   CURVED TIE      -- and only here is a bridge face needed, because the facet is curved. The
%                      sides' straight candidates still come from `envelopeCellRegions` on the
%                      polyhedral pair, each restricted to the side its own MAXIMISER falls on --
%                      which is a conic condition in x, because that maximiser is affine
%                      (`envelopeCellRegions` returns it as `.S`). The curved facet itself is
%                      `bridgeFaceAlg`'s face. 5 of 29.
%
% THE FOLD ACROSS PAIRS IS A MAX, and it needs no gluing argument: every candidate is feasible for
% its own pair, so no pair can overshoot, and the pairs cover the plane -- so the max over the faces
% IS the envelope. `.claude/assemblyProbe.m` measured both halves.
%
% WHY THE FACES ARE RETURNED RATHER THAN A SINGLE MESH. Folding the pairs pairwise with `maxQ` is
% exact and it EXPLODES: each fold overlays two meshes, and measured on the convex fixture the
% running mesh went 5, 6, 38, 35, 40, 125, 610 faces over seven of eighteen pairs -- while the
% envelope itself has a handful. The growth is real rather than a defect: the intermediate maxima
% carry pieces that later pairs dominate away, and nothing in the fold knows that yet. Controlling
% it needs a domination prune (does some candidate beat this one on the whole of its region -- an
% exact quadratic sign test on a polygon, `ratQ.quadSignOnPoly`), and that is a separate item.
%
% So this returns the FACES, which are each exact, and `envelopeEval` takes the max. That is a
% complete and checkable envelope; what it is not is a single `QuaCon`, and the two things that
% stand between it and one are the prune above and a container for a predicate-decided face.

    if nargin ~= 2
        error('envelopeAssemble:nargin', 'expects (gA, gB), two QuaCon conjugates.');
    end
    for h = {gA, gB}
        if ~isa(h{1}, 'QuaCon')
            error('envelopeAssemble:type', 'both arguments must be QuaCon; got %s.', class(h{1}));
        end
    end

    [pairs, stats] = envelopeCellPairs(gA, gB);
    info = struct('nPairs', numel(pairs), 'nPoly', 0, 'nTangential', 0, 'nStraight', 0, ...
                  'nCurved', 0, 'stats', stats);
    env = struct('quad', struct('num', {}, 'den', {}, 'con', {}, 'pair', {}), 'bridge', {{}});

    for k = 1:numel(pairs)
        pr = pairs(k);
        [nA, dA] = faceOf(gA, pr.ij(1));
        [nB, dB] = faceOf(gB, pr.ij(2));

        if strcmp(pr.kind, 'poly')
            info.nPoly = info.nPoly + 1;
            if pr.win == 'A', cells = cellsFor(nA, dA, pr.P);
            else,             cells = cellsFor(nB, dB, pr.P);
            end
            env = collect(env, cells, k);
            continue
        end

        % ---- a bridge pair: which of the three rules? ------------------------------------------
        [sgn, strict] = ratQ.quadSignOnPoly(pr.q, pr.P);
        if sgn ~= 0 && ~strict
            % TANGENTIAL: one face wins everywhere, with equality only where the tie touches.
            info.nTangential = info.nTangential + 1;
            if sgn > 0, cells = cellsFor(nA, dA, pr.P); else, cells = cellsFor(nB, dB, pr.P); end
            env = collect(env, cells, k);
            continue
        end

        if all(pr.q(1:3) == 0)
            % STRAIGHT TIE (a line): each side is a polygon, so split and take the ordinary
            % candidates. The tie becomes one more facet, which is exactly what it is.
            info.nStraight = info.nStraight + 1;
            env = collect(env, cellsFor(nA, dA, [pr.P;  pr.q(4:6)]), k);
            env = collect(env, cellsFor(nB, dB, [pr.P; -pr.q(4:6)]), k);
            continue
        end

        % ---- the tie is a CONIC: a bridge face, unless it degenerates -------------------------
        try
            face = bridgeFaceAlg(nA, dA, nB, dB, pr.P);
        catch ME
            switch ME.identifier
                case {'bridgeFaceAlg:tangentialArc'}
                    % The tie crosses the pair's LINES but meets the pair itself at a point.
                    info.nTangential = info.nTangential + 1;
                    s2 = sideWinner(nA, dA, nB, dB, pr.P);
                    if s2 > 0, cells = cellsFor(nA, dA, pr.P); else, cells = cellsFor(nB, dB, pr.P); end
                    env = collect(env, cells, k);
                    continue
                case {'bridgeFaceAlg:straightTie', 'bridgeFaceRat:noSeparateFace'}
                    % A DOUBLE LINE or a parallel pair: straight after all, so the same rule as
                    % above -- but the tie is not one half-plane, so each side is taken by the
                    % maximiser test rather than by adding a row.
                    info.nStraight = info.nStraight + 1;
                    env = collect(env, sideCells(nA, dA, pr.P, pr.q, +1), k);
                    env = collect(env, sideCells(nB, dB, pr.P, pr.q, -1), k);
                    continue
                otherwise
                    rethrow(ME)
            end
        end
        info.nCurved = info.nCurved + 1;
        face.pairIndex = k;      % a bridge face belongs to a PAIR, and eval groups by pair
        env.bridge{end+1} = face; %#ok<AGROW>
        env = collect(env, sideCells(nA, dA, pr.P, pr.q, +1), k);
        env = collect(env, sideCells(nB, dB, pr.P, pr.q, -1), k);
    end
    info.nQuadFaces = numel(env.quad);
    info.nBridgeFaces = numel(env.bridge);
end

% --------------------------------------------------------------------------------------------------
function env = collect(env, cells, k)
% Add one pair's candidates to the face list, each remembering which pair it came from.
    for i = 1:numel(cells)
        env.quad(end+1) = struct('num', cells(i).num, 'den', cells(i).den, ...
                                 'con', cells(i).con, 'pair', k); %#ok<AGROW>
    end
end

% ==================================================================================================
function [n, d] = faceOf(g, k)
    n = g.fN(k,:);  d = g.fD(k);
end

% --------------------------------------------------------------------------------------------------
function cells = cellsFor(n, d, P)
% One cell per candidate of `sup { <s,x> - g(s) : s in P }`, in the form assembleQuaConCells takes.
    N = [n(5) n(6); n(6) n(7)];  m = [n(8); n(9)];  q = n(10);
    faces = envelopeCellRegions(N, m, q, d, P);
    cells = struct('num', {}, 'den', {}, 'con', {});
    for i = 1:numel(faces)
        cells(end+1) = struct('num', faces(i).fN, 'den', faces(i).fD, ...
                              'con', rowsToConics(faces(i).R)); %#ok<AGROW>
    end
end

% --------------------------------------------------------------------------------------------------
function cells = sideCells(n, d, P, tie, side)
% The candidates of one SIDE of a tie that is not a half-plane: take the candidates over the whole
% polyhedral pair and keep each one where its OWN MAXIMISER falls on that side.
%
% That test is exact and it is a CONIC condition in x, because the maximiser is an affine map of x
% -- `envelopeCellRegions` returns it as `.S` for exactly this purpose. Substituting an affine map
% into the tie conic is `ratQ.substAffine`, the same routine the bridge's own boundary uses.
    N = [n(5) n(6); n(6) n(7)];  m = [n(8); n(9)];  q = n(10);
    faces = envelopeCellRegions(N, m, q, d, P);
    cells = struct('num', {}, 'den', {}, 'con', {});
    for i = 1:numel(faces)
        S = faces(i).S;
        c = pullback(tie, S.Mn, S.tn, S.d);
        if all(c == 0)
            % the tie is constant along this candidate's maximiser: it is entirely on one side, and
            % which one cannot be read from a zero row -- so the candidate is kept whole.
            con = rowsToConics(faces(i).R);
        else
            con = [rowsToConics(faces(i).R); ratQ.conic(c), side*sign(firstNonZero(c)), 0];
        end
        cells(end+1) = struct('num', faces(i).fN, 'den', faces(i).fD, 'con', con); %#ok<AGROW>
    end
end

% --------------------------------------------------------------------------------------------------
function s = sideWinner(nA, dA, nB, dB, P)
% Which face wins on a pair whose tie only touches it: the sign of the difference, read where it is
% not zero. `quadSignOnPoly` already decided it is one sign; this recovers which.
    c = ratQ.diffConic(nA, dA, nB, dB);
    [s, ~] = ratQ.quadSignOnPoly(c, P);
    if s == 0, s = 1; end          % identical faces: either is the winner, and they agree
end

% --------------------------------------------------------------------------------------------------
function con = rowsToConics(R)
% Half-planes [a b c] as conic rows [a b c d e f sign branch]: a line IS a conic with no quadratic
% part, and its side is +1 because the row was written as ">= 0".
    if isempty(R), con = zeros(0,8); return, end
    con = zeros(size(R,1), 8);
    for i = 1:size(R,1)
        con(i,1:6) = ratQ.conic([0 0 0, R(i,:)]);
        con(i,7) = sign(firstNonZero([0 0 0, R(i,:)]));
        con(i,8) = 0;
    end
    % ratQ.conic normalises the sign of the leading entry, so the side follows that normalisation
    % rather than the original spelling -- which is why it is recomputed here instead of assumed.
end

function v = firstNonZero(c)
    i = find(c ~= 0, 1);
    if isempty(i), v = 1; else, v = sign(c(i)); end
end

function c = pullback(q, Mn, tn, dd)
% The tie conic COMPOSED with an affine map s = (Mn x + tn)/dd, as an integer conic in x.
%
% NOT ratQ.substAffine, and the difference matters: that routine inverts the map, so it needs an
% invertible one -- and the maximiser maps here are rank-deficient BY CONSTRUCTION. A VERTEX
% candidate's maximiser is CONSTANT (rank 0) and a FACET candidate's has rank 1, which is most of
% them. Composing forward needs no inverse and is elementary: substitute and clear dd^2.
    L1 = [Mn(1,1), Mn(1,2), tn(1)];
    L2 = [Mn(2,1), Mn(2,2), tn(2)];
    c = q(1)*mulAff(L1,L1) + q(2)*mulAff(L1,L2) + q(3)*mulAff(L2,L2) ...
        + dd*(q(4)*[0 0 0, L1] + q(5)*[0 0 0, L2]) + [0 0 0 0 0, q(6)*dd^2];
    c = ratQ.chk(c, 'tie pulled back through a maximiser');
end

function c = mulAff(f, g)
% (f1 x1 + f2 x2 + f3)(g1 x1 + g2 x2 + g3), as a conic [a b c d e f].
    c = [f(1)*g(1), f(1)*g(2) + f(2)*g(1), f(2)*g(2), ...
         f(1)*g(3) + f(3)*g(1), f(2)*g(3) + f(3)*g(2), f(3)*g(3)];
end
