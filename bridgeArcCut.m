function [rows, info] = bridgeArcCut(tie, P, N1, m1, N2, m2, d, wit)
% bridgeArcCut  Cut a bridge face's sweep to the tie ARC that lies inside its own cell pair.
%
% objective: the one step between `bridgeFaceAlg` and an assembled envelope. The producer sweeps
%   the WHOLE tie curve; a cell pair owns only the part of it inside the pair, and the difference is
%   not cosmetic -- measured (DECISIONS.md 2026-09-17), on the uncut region the face equals the
%   envelope at every point whose maximiser lies in the pair (14/14) and is BEATEN at every point
%   whose maximiser lies outside it (72/72). The uncut face claims ground it loses on.
%
% [input]  tie          : the tie conic, integer [a b c d e f] (ratQ.conic form)
%          P            : k x 3 integer half-planes of the cell pair, [a b c] meaning a s1+b s2+c>=0
%          N1,m1,N2,m2  : the two faces' Hessians and linear parts, ALREADY over the common d
%          d            : that common positive integer denominator
%          wit          : [xn1 xn2 xd], a rational point of the face's region, for reading signs
% [output] rows : r x 7 integer rows [a b c d e f sign] to ADD to the face's region
%          info : .kind ('line'), .tauLo/.tauHi (the arc's extent), .nCut (rows added)
%
% ==================================================================================================
% WHAT THE CUT IS, AND WHY A LINE TIE IS THE EASY CASE AND THE COMMON ONE
%
% The sweep sends a tie point s to the segment [grad g1(s), grad g2(s)]. Restricting s to the arc
% inside P therefore adds exactly two boundaries -- the segments at the arc's two ENDPOINTS -- and
% each endpoint is the tie meeting a facet of P.
%
% When the tie is a genuine CONIC those endpoints are quadratic irrationals, the segment through one
% of them has coefficients in Q(s_e) rather than Q, and only the PAIR of conjugate segments is
% rational -- a degenerate line-pair conic PC, the product over the two conjugate roots, which is a
% symmetric function of them and so has integer coefficients. A line pair's sign condition names
% TWO OPPOSITE WEDGES, and picking the right one is the difficulty. It is REFUSED here, and the
% route that looks like it works is recorded so that it is not built a second time:
%
%   THE WEDGE ARGUMENT, AND WHY IT IS FALSE. A line cuts an ellipse into exactly two arcs, so the
%   set A swept from one of them is CONNECTED and its boundary contains the two segments. It is
%   tempting to conclude that A crosses neither LINE and therefore lies in ONE of the four wedges,
%   which would make PC's sign a valid classifier -- with the wedges' apex rational (the singular
%   point of PC) and a rational separator through it available by an integer sign test on PC's
%   quadratic part. All of that was built, and it is wrong at the first step: A does not cross the
%   SEGMENTS, but it may perfectly well cross a line's EXTENSION beyond its segment, and then A
%   spans several wedges. MEASURED on the isotropic pair over the box {-1<=s1<=3, -2<=s2<=3}: 145
%   of 296 points whose base lies strictly inside the pair were EXCLUDED by the wedge cell, with
%   the witness at PC < 0 and the excluded points at PC > 0 -- the same region, two wedges.
%   DECISIONS.md 2026-09-17.
%
%   WHAT THE FAILURE POINTS AT. The boundary of A really is the two SEGMENTS, not the two lines, and
%   a segment is not the zero set of a polynomial. So no single H-form cell over CONICS holds A:
%   either the cell list gains a subdivision that separates A's pieces from B's within each sign of
%   PC, or the subdivision axis needs a level that can carry a bounded arc. Which of the two is
%   right is not settled, and a guess here puts a face on the wrong side of its own boundary.
%
% When the tie is a LINE -- tie(1:3) all zero -- every one of those objects is RATIONAL: the arc is
% an interval of the line with rational endpoints, each endpoint is a rational point, and the
% segment through it is a rational LINE, i.e. an ordinary half-plane with an unambiguous side. No
% line pair, no ambiguity, nothing to disambiguate.
%
% AND THE LINE TIE IS NOT A SPECIAL CASE, it is the one the pipeline produces most. The tie of two
% cells is the conic { g1 = g2 }, whose quadratic part is (H1-H2)/2; [COAP]/[JOGO] Theorem 6 --
% CONJ_FIELD_PROOF.md Theorem 3 -- makes H1 = H2 exactly when the two pieces compared are ADJACENT
% pieces of f, and two cells of the SAME piece's conjugate share their Hessian outright. So the
% degenerate tie is the rule among neighbours and the genuine conic is what non-adjacent pairs give.

    if nargin ~= 8
        error('bridgeArcCut:nargin', 'expects (tie, P, N1, m1, N2, m2, d, wit).');
    end
    tie = ratQ.chk(tie(:).', 'tie conic');
    if any(tie(1:3) ~= 0)
        error('bridgeArcCut:conicTie', ...
            ['the tie is a genuine CONIC, so each cut segment lives in Q(s_e) and only the ' ...
             'conjugate PAIR is rational -- a line pair, whose sign condition names two OPPOSITE ' ...
             'wedges. The wedge-and-separator construction that would resolve that was built and ' ...
             'MEASURED WRONG (145 of 296 interior points excluded): the swept region crosses a ' ...
             'cut line''s EXTENSION beyond its own segment, so it spans several wedges. See this ' ...
             'file''s header and DECISIONS.md 2026-09-17 before attempting it again.']);
    end
    D = tie(4);  E = tie(5);  F = tie(6);
    if D == 0 && E == 0
        error('bridgeArcCut:degenerateTie', ...
            'the tie has no linear part either, so it names no curve at all.');
    end

    % ---- the arc, as an interval of the line's own parameter ------------------------------------
    % s(tau) = (s0n + tau*w)/nd with w along the line and s0n its foot, both integer: then
    % D*s1 + E*s2 + F = 0 identically, so tau parameterises the tie and nothing else.
    w   = [E; -D];
    s0n = [-F*D; -F*E];
    nd  = ratQ.chk(D^2 + E^2, 'line normalisation');

    % Row a's + c >= 0 becomes A*tau + B >= 0 with A, B integers.
    lo = [];  hi = [];                                   % rational bounds as [num den], den > 0
    for i = 1:size(P,1)
        a = P(i,1:2).';  c = P(i,3);
        A = ratQ.chk(a.'*w, 'arc bound slope');
        B = ratQ.chk(a.'*s0n + c*nd, 'arc bound constant');
        if A == 0
            % The row does not vary along the tie: it holds on all of it or on none of it.
            if B < 0
                error('bridgeArcCut:emptyArc', ...
                    ['facet %d excludes the ENTIRE tie line, so this pair contains no tie point ' ...
                     'and carries no bridge face. The pair classification and this routine ' ...
                     'disagree, which is a caller defect.'], i);
            end
            continue
        end
        % A*tau + B >= 0  <=>  tau >= -B/A (A > 0)  or  tau <= -B/A (A < 0)
        b = [-B, A];
        if A < 0, b = [B, -A]; end                       % keep the denominator positive
        if A > 0, lo = tighter(lo, b, @gt); else, hi = tighter(hi, b, @lt); end
    end
    if ~isempty(lo) && ~isempty(hi) && ratLess(hi, lo)
        error('bridgeArcCut:emptyArc', ...
            'the facets leave no tie point at all: the arc''s bounds cross.');
    end

    % ---- one half-plane per FINITE endpoint -----------------------------------------------------
    rows = zeros(0,7);
    for e = {lo, hi}
        tau = e{1};
        if isempty(tau), continue, end                   % the arc runs to infinity on that side
        row = endpointRow(tau, w, s0n, nd, N1, m1, N2, m2, d, wit);
        if isempty(row), continue, end                   % the segment degenerates to a point
        rows(end+1,:) = row; %#ok<AGROW>
    end
    info = struct('kind', 'line', 'tauLo', {lo}, 'tauHi', {hi}, 'nCut', size(rows,1));
end

% ==================================================================================================
function row = endpointRow(tau, w, s0n, nd, N1, m1, N2, m2, d, wit)
% The half-plane whose boundary is the SEGMENT swept from the arc endpoint s(tau).
%
% Every quantity here is integral because tau is: s_e = (s0n*td + tn*w)/(nd*td), and the two
% gradients over that common denominator are (Nk s_e + m_k)/d, so the segment's DIRECTION is an
% integer difference and its normal is that direction rotated. The line through the first gradient
% with that normal is the boundary, and the side is read at the witness -- which is a point of the
% face, so it names the side the face is on rather than a side someone chose.
    tn = tau(1);  td = tau(2);
    seN = ratQ.chk(s0n*td + tn*w, 'arc endpoint');       % s_e = seN / seD
    seD = ratQ.chk(nd*td,         'arc endpoint denominator');

    P1n = ratQ.chk(N1*seN + seD*m1, 'grad g1 at the endpoint');   % grad g1 = P1n/(d*seD)
    P2n = ratQ.chk(N2*seN + seD*m2, 'grad g2 at the endpoint');
    dir = P2n - P1n;
    if all(dir == 0), row = []; return, end              % the two gradients agree: no segment

    n = [-dir(2); dir(1)];
    % n'x - n'grad g1 = 0, cleared by d*seD
    c = ratQ.conic([0 0 0, ratQ.chk(d*seD*n(1), 'cut row'), ratQ.chk(d*seD*n(2), 'cut row'), ...
                    ratQ.chk(-n.'*P1n, 'cut row')]);
    s = signAtRational(c, wit(1:2).', wit(3));
    if s == 0
        error('bridgeArcCut:witnessOnCut', ...
            ['the witness lies exactly ON a cut boundary, so it names no side. A witness in the ' ...
             'INTERIOR of the arc''s sweep decides it; this one was chosen before the cut was ' ...
             'known.']);
    end
    row = [c, s];
end

% ==================================================================================================
function b = tighter(b, cand, cmp)
% Keep the tighter of two rational bounds, or the new one when there was none.
    if isempty(b), b = cand; return, end
    if cmp(cand(1)*b(2), b(1)*cand(2)), b = cand; end     % both denominators are positive
end
function tf = ratLess(p, q)
    tf = p(1)*q(2) < q(1)*p(2);
end
function s = signAtRational(C, xn, xd)
    v = C(1)*xn(1)^2 + C(2)*xn(1)*xn(2) + C(3)*xn(2)^2 + ...
        C(4)*xn(1)*xd + C(5)*xn(2)*xd + C(6)*xd^2;
    s = sign(ratQ.chk(v, 'cut row at the witness'));
end
