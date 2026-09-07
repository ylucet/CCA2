function bridgeEnvelopeProbe(varargin)
% ITEM 3: THE BRIDGE PIPELINE RUN END TO END ON A REAL TRIANGLE PAIR, checked against the
% definition.
%
% WHAT THIS PUTS IN SEQUENCE. Every stage has been verified on its own; nothing had run them
% together on a fixture a user could write:
%
%   conjQ each triangle              ->  gA, gB, 7 dual cells each
%   FILTER 1  cells meet in 2-D      ->  ratQ.feasible2(..., strict)
%   FILTER 2  g1 - g2 changes sign   ->  ratQ.quadSignOnPoly
%   per surviving pair, the REGION   ->  the two boundary curves, via ratQ.substAffine
%   then CHECK against sup_s <s,x> - max(gA(s), gB(s))
%
% THE FIXTURE IS THE CORPUS ONE. The unit square split by its diagonal into two positive definite
% quadratic pieces that agree on the shared edge, with a NON-CONVEX kink -- exactly the input
% biconjQ refuses with convexPieceCoupling, and the reason the bridge exists at all.
%
% WHAT IS ASKED, AFTER A FIRST DRAFT ASKED THE WRONG THING. The first version tried to rebuild the
% WHOLE envelope, bridge faces and non-bridge faces alike, and compare it to the definition. It
% disagreed by up to 1.94, and the cause was not the bridges: a NON-bridge face needs the Legendre
% value of one cell RESTRICTED TO THAT CELL, and over a vertex or edge cell the conjugate is affine
% or rank-deficient, so that is a constrained problem this pipeline does not solve yet. Taking
% unconstrained maximisers and skipping the singular cells produced numbers from formulas outside
% their own validity -- a wrong answer from the checker, not from the thing being checked.
%
% So the question is narrowed to the one the bridge machinery is responsible for, which is also the
% one that decides whether the filters are safe:
%
%       whenever the definition's own maximiser TIES the two pieces, is that cell pair in the
%       filtered list, and does the tie-restricted sup reproduce the envelope's value?
%
% A tie the filters dropped would be a bridge face silently missing from the envelope -- exactly
% the unsound direction filter 2 was built never to take.
%
%   bridgeEnvelopeProbe            % 150 sample points
%   bridgeEnvelopeProbe(1000)      % longer

    n = 150;
    if nargin >= 1, n = varargin{1}; end

    c = -2;                                   % a NON-CONVEX kink: the case that needs a bridge
    fA = tri([0 0; 1 0; 1 1], 0);
    fB = tri([0 0; 1 1; 0 1], c);
    gA = conjQ(fA);   gB = conjQ(fB);
    fprintf('\npieces: %d dual cells and %d dual cells\n', gA.nf, gB.nf);

    % ---- the two filters -----------------------------------------------------------------------
    pairs = struct('i', {}, 'j', {}, 'P', {});
    nMeet = 0;
    for i = 1:gA.nf
        RA = cellRows(gA, i);   if isempty(RA), continue, end
        for j = 1:gB.nf
            RB = cellRows(gB, j);  if isempty(RB), continue, end
            P = [RA; RB];
            if ~ratQ.feasible2(P, true), continue, end
            nMeet = nMeet + 1;
            q = tieConic(gA, i, gB, j);
            if all(q == 0), continue, end
            % STRICTLY one sign, and the second output is what says so. A +1 from
            % quadSignOnPoly only means q >= 0, which permits q = 0 somewhere -- a TANGENTIAL tie,
            % which is still a tie and still carries a bridge. Pruning on `s ~= 0` dropped cell
            % pair (5,5) of this very fixture, where the definition ties; see DECISIONS 2026-09-06.
            [sgn, strict] = ratQ.quadSignOnPoly(q, P);
            if sgn ~= 0 && strict, continue, end
            pairs(end+1) = struct('i', i, 'j', j, 'P', P); %#ok<AGROW>
        end
    end
    fprintf('cell pairs: %d total, %d meet in 2-D, %d survive filter 2\n', ...
            gA.nf * gB.nf, nMeet, numel(pairs));

    % ---- each surviving pair's region boundaries -------------------------------------------------
    nConic = 0;  nPoint = 0;
    for k = 1:numel(pairs)
        q  = tieConic(gA, pairs(k).i, gB, pairs(k).j);
        C1 = imageConic(q, cellGrad(gA, pairs(k).i));
        C2 = imageConic(q, cellGrad(gB, pairs(k).j));
        if isempty(C1) || isempty(C2), nPoint = nPoint + 1; else, nConic = nConic + 1; end
    end
    fprintf('surviving pairs:');
    for k = 1:numel(pairs), fprintf(' (%d,%d)', pairs(k).i, pairs(k).j); end
    fprintf('\n');
    fprintf(['regions: %d bounded by two conics, %d with a DEGENERATE boundary -- a vertex or ' ...
             'edge cell,\n'], nConic, nPoint);
    fprintf(['         whose conjugate is affine or rank deficient, so that side of the region ' ...
             'collapses\n         to a point or a line\n']);

    % ---- the check -------------------------------------------------------------------------------
    rng(2026);
    X = [rand(n,1)*4 - 1.5, rand(n,1)*4 - 1.5];
    nTie = 0;  nInList = 0;  nSame = 0;  worst = 0;  missed = {};
    for t = 1:n
        x = X(t,:).';
        [vDef, sStar] = envelopeByDefinition(gA, gB, x);
        [a, iA] = gA.eval(sStar.');   [b, jB] = gB.eval(sStar.');
        if iA == 0 || jB == 0, continue, end
        if abs(a - b) > 1e-6 * max(1, max(abs([a b]))), continue, end    % not a tie: not a bridge

        % A TIE IS NOT YET A BRIDGE, and this fixture is where that became clear. Cells 1 of the
        % two conjugates have IDENTICAL face functions -- both pieces share the primal vertex
        % (0,0) and take the same value there, so their conjugates agree on a whole 2-D cell
        % rather than crossing along a curve. The maximiser then "ties" everywhere in that cell.
        % No new face arises: the common function is ALREADY a face of each conjugate, so the
        % envelope there is an ordinary vertex face and there is nothing to eliminate. Counted
        % separately rather than reported as a bridge the filters missed, which is what the first
        % run of this probe called it -- 149 of 150 points, all in pair (1,1).
        if all(tieConic(gA, iA, gB, jB) == 0), nSame = nSame + 1; continue, end

        nTie = nTie + 1;
        k = find([pairs.i] == iA & [pairs.j] == jB, 1);
        if isempty(k)
            if numel(missed) < 8, missed{end+1} = sprintf('(%d,%d)', iA, jB); end %#ok<AGROW>
            continue
        end
        nInList = nInList + 1;
        vTie = supOnTie(gA, pairs(k).i, gB, pairs(k).j, pairs(k).P, x);
        worst = max(worst, abs(vDef - vTie) / max(1, abs(vDef)));
    end

    fprintf('\n%d sample points: %d had a TIED maximiser, i.e. are genuine bridge points\n', n, nTie);
    fprintf('  of those, %d were in the filtered pair list and %d were NOT\n', nInList, nTie - nInList);
    if ~isempty(missed)
        fprintf('  pairs a filter DROPPED though the definition ties there: %s\n', strjoin(missed, ' '));
    end
    fprintf('  worst relative gap, definition vs the tie-restricted sup: %.3e\n', worst);
    if nTie > 0 && nInList == nTie && worst < 1e-6
        fprintf(['\nEvery bridge point the definition found was claimed by a pair the two filters\n' ...
                 'kept, and the tie-restricted sup reproduced the envelope value at all of them.\n' ...
                 'On this fixture the filters drop no real bridge and the bridge value is right.\n']);
    end
end

% ------------------------------------------------------------------------------------------------
function q = tieConic(gA, i, gB, j)
% The locus where the two cells' conjugates agree, cleared to integers over the two denominators.
%
% THE WEIGHTED BASIS CARRIES HALVES, and forgetting that mis-weights the curve. A face is stored
% in RatCon's 10-wide basis, whose value is
%       n5 x^2/2 + n6 xy + n7 y^2/2 + n8 x + n9 y + n10 ,
% so the raw slice n(5:10) is NOT the conic [a b c d e f]: the squared terms are half what that
% spelling means. The conic of the same function, cleared to integers by doubling -- free, since a
% conic has no scale -- is [n5, 2*n6, n7, 2*n8, 2*n9, 2*n10].
%
% This was wrong here until 2026-09-06 and it moved real numbers: the difference conic passed to
% quadSignOnPoly named a DIFFERENT CURVE from the one where the two faces actually agree, so both
% filter-2 columns were measured against the wrong tie. Found by the end-to-end run, where a cell
% pair the definition said was NOT tied came back tied. DECISIONS 2026-09-06.
    D = gA.fN(i,:) * gB.fD(j) - gB.fN(j,:) * gA.fD(i);
    q = [D(5), 2*D(6), D(7), 2*D(8), 2*D(9), 2*D(10)];
end

function g = cellGrad(gc, k)
% The gradient of face k as an affine map, returned as a struct with grad(s) = A s + b. The face is
% stored in the 10-wide weighted basis, whose value is
%     c5 s1^2/2 + c6 s1 s2 + c7 s2^2/2 + c8 s1 + c9 s2 + c10,
% so the gradient is [c5 c6; c6 c7] s + [c8; c9], over the face denominator.
    n = gc.fN(k,:);  d = gc.fD(k);
    g = struct('A', [n(5), n(6); n(6), n(7)] / d, 'b', [n(8); n(9)] / d);
end

function cimg = imageConic(q, g)
% The image of the tie conic under an invertible affine map, exactly. Empty when the map is
% singular -- a vertex cell (A = 0) or an edge cell (rank 1) -- where the image is a point or a
% line and substAffine correctly refuses, since that is not a change of variables.
    A = g.A;  b = g.b;
    if abs(det(A)) < 1e-12, cimg = []; return, end
    [An, ad] = intify(A);  [bn, bd] = intify(b);
    L = lcm(ad, bd);  An = An * (L/ad);  bn = bn * (L/bd);
    adjA = [An(2,2), -An(1,2); -An(2,1), An(1,1)];
    dAn  = ratQ.detExact(An);
    try
        [n2, ~] = ratQ.substAffine([0 0 0 0, 2*q(1), q(2), 2*q(3), q(4), q(5), q(6)], 1, ...
                                   ratQ.chk(L*adjA), ratQ.chk(-adjA*bn), dAn);
    catch
        cimg = [];  return
    end
    cimg = ratQ.conic([n2(5), 2*n2(6), n2(7), 2*n2(8), 2*n2(9), 2*n2(10)]);
end

% ---- the oracle ---------------------------------------------------------------------------------
function [v, sBest] = envelopeByDefinition(gA, gB, x)
% f**(x) = sup_s <s,x> - max(gA(s), gB(s)), and WHERE it is attained. The objective is concave, so
% a few starts and a simplex search settle it. Shares nothing with the filters or the region
% construction, which is what makes it an oracle rather than a restatement.
    h = @(s) -( s(:).' * x - maxConj(gA, gB, s(:).') );
    best = inf;  sBest = [0;0];
    for s0 = [x, [0;0], -x, 2*x]
        s = fminsearch(h, s0, optimset('Display', 'off', 'TolX', 1e-12, 'TolFun', 1e-12));
        if h(s) < best, best = h(s); sBest = s; end
    end
    v = -best;
end

function v = maxConj(gA, gB, s)
    v = max(gA.eval(s), gB.eval(s));
end

function v = supOnTie(gA, i, gB, j, P, x)
% The sup of <s,x> - gA_i(s) over the pair's tie INSIDE ITS OWN CELL. Restricting to the cell is
% not a detail: off the cell that face formula is not the conjugate, so an unrestricted sup would
% maximise the wrong function and quietly beat the true value.
    S = tieInCell(gA, i, gB, j, P);
    v = -inf;
    for t = 1:size(S,1)
        s = S(t,:).';
        v = max(v, s.' * x - faceVal(gA, i, s));
    end
end

function S = tieInCell(gA, i, gB, j, P)
% Points of the tie that lie in the cell intersection. SOLVED in the second coordinate rather than
% bisected along a ray: a ray-based search misses a branch it does not cross, which DECISIONS
% 2026-09-06 records costing two wrong conclusions earlier in this same line of work.
    q = tieConic(gA, i, gB, j);
    S = zeros(0,2);
    for u = linspace(-30, 30, 60001)
        a = q(3);  b = q(2)*u + q(5);  cc = q(1)*u^2 + q(4)*u + q(6);
        if abs(a) < 1e-14
            if abs(b) > 1e-14, S = addIfIn(S, [u, -cc/b], P); end
            continue
        end
        disc = b^2 - 4*a*cc;
        if disc < 0, continue, end
        r = sqrt(disc);
        S = addIfIn(S, [u, (-b + r)/(2*a)], P);
        S = addIfIn(S, [u, (-b - r)/(2*a)], P);
    end
end

function S = addIfIn(S, p, P)
    if all(P(:,1:2) * p(:) + P(:,3) >= -1e-9), S(end+1,:) = p; end
end

function v = faceVal(g, k, s)
    n = g.fN(k,:);
    v = (n(5)*s(1)^2/2 + n(6)*s(1)*s(2) + n(7)*s(2)^2/2 + ...
         n(8)*s(1) + n(9)*s(2) + n(10)) / g.fD(k);
end

% ---- fixture and small helpers -------------------------------------------------------------------
function f = tri(V, c)
% A positive definite quadratic on one triangle, perturbed by c times the linear form that vanishes
% on the shared diagonal y = x. The two pieces therefore agree there for EVERY c, so one knob
% sweeps convex, flat and non-convex without ever producing a mesh that is not a function --
% DECISIONS 2026-09-05's construction.
    q = [0 0 0 0, 2, 0, 2, 0, 0, 0];
    q(8) = q(8) - c;
    q(9) = q(9) + c;
    f = QuaPol(V, [1 2 1; 2 3 1; 3 1 1], q, [1 0; 1 0; 1 0]);
end

function rows = cellRows(g, k)
    rows = zeros(0,3);
    FC = g.FC{k};
    for r = 1:size(FC,1)
        e = FC(r,1);  sgn = FC(r,2);
        cc = g.EcQ(e,:);
        if any(cc(1:3) ~= 0), rows = zeros(0,3); return, end
        if sgn == 0, continue, end
        rows(end+1,:) = sgn * cc(4:6); %#ok<AGROW>
    end
end

function [N, d] = intify(X)
    d = 1;
    for v = X(:).'
        [~, dv] = rat(v, 1e-12);
        d = lcm(d, abs(dv));
    end
    N = round(X * d);
end
