function ok = proveStageD(seed, nPairs)
% proveStageD  STAGE D: the MAX across triangles -- how many arcs f* can carry, and whether an
%   arc-vs-arc face pair is ever actually reachable.
%
% f* = max_T (q + I_T)*, so beyond each triangle's own dual boundaries (Stage C) there are the
% CROSS boundaries, where two different triangles' candidates tie: e_A^{T1} - e_B^{T2}. Stage D
% classifies those by the same two discriminants, and then answers the question that actually
% governs maxQuaPar:
%
%   D3  Can a curved face of g1 and a curved face of g2 OVERLAP?
%
% That is the whole content of the arc-vs-arc gap. clipByFace refuses a pair in which both faces
% carry an arc -- but if such a pair never has a two-dimensional intersection in the supported
% regime (adjacent sub-triangles of one domain, which is the only regime maxQuaPar claims), then
% the refusal is guarding an EMPTY case and the surgery never needs writing: the fix is to test
% emptiness first, not to implement conic-vs-conic clipping.
%
% [output] ok : true iff every obligation verified.
%
% Instances are the same configuration maxQuaParTest's own curved fixtures use -- a convex
% quadrilateral split by a diagonal -- since that is what "two adjacent sub-pieces of the same
% originally-nonconvex domain" means, and it is the scoping caveat maxQuaPar's header states.

    if nargin < 1 || isempty(seed),   seed = 20260802; end
    if nargin < 2 || isempty(nPairs), nPairs = 25;     end
    rng(seed);
    fprintf('\n=== STAGE D: the max across triangles ===\n');

    E3 = [1 2 1; 2 3 1; 3 1 1]; F3 = [1 0; 1 0; 1 0];
    nBothCurved = 0; nBothCurvedOverlapping = 0; nUsable = 0;
    nCross = 0; nCrossLine = 0; nCrossPair = 0; nCrossPara = 0; nCrossOther = 0;
    arcCounts = [];

    % The main loop TARGETS the both-curved configuration. Drawing pairs at random and taking
    % whatever conjPieceCPLQ accepts gave 3 usable instances out of 25, none of them in the
    % (1CE,1CE) bucket -- so it measured 0 both-curved pairs and said nothing. D3a below shows
    % that bucket is 14% of all splits, so it is filtered for here explicitly.
    Qbil = [0 1; 1 0];
    drawn = 0;
    for it = 1:20000
        if drawn >= nPairs, break, end
        [T1, T2] = randomAdjacentPair();
        if isempty(T1), continue, end
        if numConvexEdgesOf(T1, Qbil) ~= 1 || numConvexEdgesOf(T2, Qbil) ~= 1, continue, end
        drawn = drawn + 1;
        try
            g1 = conjPieceCPLQ(QuaPol(T1, E3, [0 1 0 0 0 0], F3));
            g2 = conjPieceCPLQ(QuaPol(T2, E3, [0 1 0 0 0 0], F3));
        catch
            continue                          % a sub-triangle Step 2 cannot take directly
        end
        nUsable = nUsable + 1;

        % ---- D1: classify every CROSS face-pair difference ---------------------------------
        for k = 1:g1.nf
            for l = 1:g2.nf
                [d2, D3] = conicOfDifference(g1.f(k,:), g2.f(l,:));
                sc = max(1, norm(g1.f(k,5:10)) + norm(g2.f(l,5:10)));
                nCross = nCross + 1;
                if abs(d2) <= 1e-7*sc^2
                    if abs(D3) <= 1e-9*sc^3, nCrossLine = nCrossLine + 1;
                    else,                    nCrossPara = nCrossPara + 1;
                    end
                elseif abs(D3) <= 1e-9*sc^3, nCrossPair = nCrossPair + 1;
                else,                        nCrossOther = nCrossOther + 1;
                end
            end
        end

        % ---- D3: do a curved face of g1 and a curved face of g2 ever overlap? --------------
        c1 = curvedFaces(g1); c2 = curvedFaces(g2);
        if ~isempty(c1) && ~isempty(c2)
            nBothCurved = nBothCurved + numel(c1)*numel(c2);
            if facesOverlap(g1, c1, g2, c2)
                nBothCurvedOverlapping = nBothCurvedOverlapping + 1;
            end
        end

        % ---- D2: arcs in the assembled max, where it assembles ------------------------------
        try
            g = maxQuaPar(g1, g2);
            arcCounts(end+1) = sum(any(g.Ec ~= 0, 2)); %#ok<AGROW>
        catch
            % refused or failed to assemble; counted by D3 instead
        end
    end

    % ---- D3a: is the both-curved configuration REACHABLE in the supported regime at all? -----
    % The loop above measures only what conjPieceCPLQ accepts directly, which turned out to be a
    % small and biased subsample. The reachability question needs no conjugation: a sub-triangle's
    % conjugate is curved exactly when it has ONE convex edge (Stage C), so count how often BOTH
    % halves of a diagonal-split quadrilateral are 1CE.
    cnt = zeros(4,4); nQuad = 0;
    for it = 1:2000
        [A, B] = randomAdjacentPair();
        if isempty(A), continue, end
        nQuad = nQuad + 1;
        Q = [0 1; 1 0];                       % f = x*y, already the bilinear frame
        n1 = numConvexEdgesOf(A, Q); n2 = numConvexEdgesOf(B, Q);
        cnt(n1+1, n2+1) = cnt(n1+1, n2+1) + 1;
    end
    fprintf('  D3a nCE of the two halves over %d quadrilaterals (rows = T1, cols = T2):\n', nQuad);
    fprintf('        %6s %6s %6s %6s\n', '0CE','1CE','2CE','3CE');
    for i = 1:4
        fprintf('    %sCE %6d %6d %6d %6d\n', num2str(i-1), cnt(i,1), cnt(i,2), cnt(i,3), cnt(i,4));
    end
    fprintf('     BOTH halves 1CE (both conjugates curved): %d of %d\n', cnt(2,2), nQuad);

    fprintf('  usable adjacent pairs                     : %d\n', nUsable);
    fprintf('  D1 cross-triangle candidate differences   : %d total\n', nCross);
    fprintf('       line %d, line pair %d, PARABOLA %d, other %d\n', ...
            nCrossLine, nCrossPair, nCrossPara, nCrossOther);
    fprintf('  D3 both-curved face pairs                 : %d\n', nBothCurved);
    fprintf('     instances where two curved faces OVERLAP: %d\n', nBothCurvedOverlapping);
    if ~isempty(arcCounts)
        fprintf('  D2 arcs in the assembled max              : min %d, max %d over %d results\n', ...
                min(arcCounts), max(arcCounts), numel(arcCounts));
    end

    ok = (nCrossOther == 0);
    if nCrossOther > 0
        fprintf('  FAIL: %d cross differences are non-degenerate non-parabolas\n', nCrossOther);
    end
    fprintf('=== STAGE D %s ===\n', tf(ok));
end

% ============================================================================================
function tf_ = facesOverlap(g1, c1, g2, c2)
% Does some curved face of g1 share two-dimensional area with some curved face of g2? Decided by
% locating grid points in BOTH meshes -- QuaPar.eval returns the face index as its second output,
% so this uses the same point location the pipeline itself uses rather than a reimplementation.
    tf_ = false;
    R = 3*max(1, max(abs([g1.V(:); g2.V(:)])));
    t = linspace(-R, R, 41);
    for i = 1:numel(t)
        for j = 1:numel(t)
            s = [t(i), t(j)];
            [~, r1] = g1.eval(s);
            if isempty(r1) || ~ismember(r1, c1), continue, end
            [~, r2] = g2.eval(s);
            if isempty(r2) || ~ismember(r2, c2), continue, end
            tf_ = true; return
        end
    end
end

function idx = curvedFaces(g)
% Faces incident to an edge carrying a nonzero conic.
    idx = [];
    if isempty(g.Ec), return, end
    for j = 1:size(g.Ec,1)
        if any(g.Ec(j,:) ~= 0)
            idx = [idx, g.F(j,1), g.F(j,2)]; %#ok<AGROW>
        end
    end
    idx = unique(idx(idx >= 1));
end

function [d2, D3] = conicOfDifference(frowA, frowB)
    QA = [frowA(5) frowA(6); frowA(6) frowA(7)]; LA = [frowA(8); frowA(9)]; kA = frowA(10);
    QB = [frowB(5) frowB(6); frowB(6) frowB(7)]; LB = [frowB(8); frowB(9)]; kB = frowB(10);
    Qd = QA - QB; Ld = LA - LB; kd = kA - kB;
    a = Qd(1,1)/2; bb = Qd(1,2); c = Qd(2,2)/2; d = Ld(1); e = Ld(2); f = kd;
    d2 = bb^2 - 4*a*c;
    D3 = det([a bb/2 d/2; bb/2 c e/2; d/2 e/2 f]);
end

function [T1, T2] = randomAdjacentPair()
% A convex quadrilateral split by a diagonal -- maxQuaPar's supported regime. Returns [] if the
% draw is not convex.
    T1 = []; T2 = [];
    for k = 1:40
        P = 4*randn(4,2);
        c = mean(P,1);
        [~, ord] = sort(atan2(P(:,2)-c(2), P(:,1)-c(1)));
        P = P(ord,:);
        if ~isConvexQuad(P), continue, end
        T1 = toCCW(P([1 2 3],:));
        T2 = toCCW(P([1 3 4],:));
        return
    end
end

function tf_ = isConvexQuad(P)
    tf_ = true; n = 4; sgn = 0;
    for i = 1:n
        a = P(i,:); b = P(mod(i,n)+1,:); c = P(mod(i+1,n)+1,:);
        cr = (b(1)-a(1))*(c(2)-b(2)) - (b(2)-a(2))*(c(1)-b(1));
        if abs(cr) < 1e-9, tf_ = false; return, end
        if sgn == 0, sgn = sign(cr); elseif sign(cr) ~= sgn, tf_ = false; return, end
    end
end

function T = toCCW(T)
    a = (T(2,1)-T(1,1))*(T(3,2)-T(1,2)) - (T(2,2)-T(1,2))*(T(3,1)-T(1,1));
    if a < 0, T = T([1 3 2],:); end
end

function n = numConvexEdgesOf(V, Q)
% Convex edges of triangle V for the quadratic with Hessian Q, counted as Stage A3 established:
% edge direction d is convex iff d'Qd > 0. For Q = [0 1; 1 0] (f = x*y) this is the positive-slope
% test classifyConvexEdges uses, but written in the frame-free form so it needs no bilinear map.
    ed = [1 2; 2 3; 3 1]; n = 0;
    for t = 1:3
        d = (V(ed(t,2),:) - V(ed(t,1),:))';
        if d'*Q*d > 1e-9*max(1,norm(d)^2), n = n + 1; end
    end
end

function s = tf(v), if v, s = 'PASS'; else, s = 'FAIL'; end, end
