function [cells, info] = bridgeArcCut(tie, P, N1, m1, N2, m2, d, wit, C1, C2, kp)
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
%          C1, C2       : the region's two boundary conics -- needed to verify a witness of the
%                         COMPLEMENTARY arc, which is how a conic cut learns which side it is NOT on
%          kp           : the tie's CONSTANT q1 - q2 over the same d. The tie conic has lost it to
%                         its own scale, and the witness search needs the tie as a FUNCTION
% [output] cells : cell array of r x 8 integer row blocks [a b c d e f sign branch] to ADD to the
%                  face's region. ONE block for a line tie; a conic cut splits the face, so the
%                  face is the UNION over the blocks -- see below.
%          info  : .kind ('line' | 'conic'), .nCut (rows in the first block), .nCells
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
% rational -- a degenerate line-pair conic PC. That case is REFUSED, and TWO constructions for it
% have now been built and measured wrong. Both are written out because both are plausible enough to
% be re-derived by anyone who meets the problem.
%
%   ATTEMPT 1, THE WEDGE. A line cuts a connected tie into two arcs, so the swept set A of one arc
%   is CONNECTED and its boundary contains the two segments -- so, it seemed, A crosses neither LINE
%   and lies in one of the four WEDGES, making PC's sign a classifier. FALSE: A does not cross the
%   SEGMENTS, but it may cross a line's EXTENSION beyond its segment. Measured: 145 of 296 interior
%   points excluded.
%
%   ATTEMPT 2, NAMING THE LINE. `ratQ.pairLineData` decomposes the pair and a branch index names one
%   of the two conjugate lines exactly, so the four SIGN PATTERNS of (L1, L2) are available as
%   ordinary H-form conditions. On the first fixtures measured, every pattern was purely A or purely
%   B, which would make A the complement of one pattern -- a union of two cells. FALSE in general:
%   on an anisotropic pair one pattern holds BOTH, 1399 points of A and 2218 of B. The pattern cell
%   is DISCONNECTED there, with one component reachable from A and the other from B, and no further
%   condition on C1, C2, L1, L2 separates them.
%
% WHAT THE TWO FAILURES HAVE IN COMMON, which is the thing to take to the next attempt: the face's
% boundary is a SEGMENT -- a bounded piece of a line -- and every representation tried describes it
% with the unbounded curve the segment sits on. `ratQ.pairLineData` and the `branch` column stay,
% because naming one line of a conjugate pair is needed by any construction at all and is now exact;
% what is missing is a way to bound it. DECISIONS.md 2026-09-17.
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

    if nargin < 11
        error('bridgeArcCut:nargin', ...
            'expects (tie, P, N1, m1, N2, m2, d, wit, C1, C2, kp).');
    end
    tie = ratQ.chk(tie(:).', 'tie conic');
    if any(tie(1:3) ~= 0)
        error('bridgeArcCut:conicTie', ...
            ['the tie is a genuine CONIC, so each cut boundary is a SEGMENT of an irrational ' ...
             'line and only the conjugate PAIR is rational. Two constructions for that have been ' ...
             'built and MEASURED WRONG -- the wedge (145 of 296 interior points excluded) and the ' ...
             'sign pattern of the two named lines (one pattern holding 1399 points of the face ' ...
             'and 2218 of its complement). Read this file''s header and DECISIONS.md 2026-09-17 ' ...
             'before attempting a third.']);
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
    nd  = ratQ.chkB(D^2 + E^2, D^2 + E^2, 'line normalisation');

    % Row a's + c >= 0 becomes A*tau + B >= 0 with A, B integers.
    lo = [];  hi = [];                                   % rational bounds as [num den], den > 0
    for i = 1:size(P,1)
        a = P(i,1:2).';  c = P(i,3);
        A = ratQ.chkB(a.'*w, abs(a).'*abs(w), 'arc bound slope');
        B = ratQ.chkB(a.'*s0n + c*nd, abs(a).'*abs(s0n) + abs(c*nd), 'arc bound constant');
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
    rows = zeros(0,8);
    for e = {lo, hi}
        tau = e{1};
        if isempty(tau), continue, end                   % the arc runs to infinity on that side
        row = endpointRow(tau, w, s0n, nd, N1, m1, N2, m2, d, wit);
        if isempty(row), continue, end                   % the segment degenerates to a point
        rows(end+1,:) = row; %#ok<AGROW>
    end
    cells = {rows};
    info = struct('kind', 'line', 'nCut', size(rows,1), 'nCells', 1);
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
    seN = ratQ.chkB(s0n*td + tn*w, abs(s0n)*abs(td) + abs(tn)*abs(w), 'arc endpoint');   % s_e = seN / seD
    seD = ratQ.chk(nd*td,         'arc endpoint denominator');

    P1n = ratQ.chkB(N1*seN + seD*m1, abs(N1)*abs(seN) + abs(seD)*abs(m1), 'grad g1 at the endpoint');   % grad g1 = P1n/(d*seD)
    P2n = ratQ.chkB(N2*seN + seD*m2, abs(N2)*abs(seN) + abs(seD)*abs(m2), 'grad g2 at the endpoint');
    dir = P2n - P1n;
    if all(dir == 0), row = []; return, end              % the two gradients agree: no segment

    n = [-dir(2); dir(1)];
    % n'x - n'grad g1 = 0, cleared by d*seD
    c = ratQ.conic([0 0 0, ratQ.chk(d*seD*n(1), 'cut row'), ratQ.chk(d*seD*n(2), 'cut row'), ...
                    ratQ.chkB(-n.'*P1n, abs(n).'*abs(P1n), 'cut row')]);
    s = signAtRational(c, wit(1:2).', wit(3));
    if s == 0
        error('bridgeArcCut:witnessOnCut', ...
            ['the witness lies exactly ON a cut boundary, so it names no side. A witness in the ' ...
             'INTERIOR of the arc''s sweep decides it; this one was chosen before the cut was ' ...
             'known.']);
    end
    row = [c, s, 0];      % branch 0: a straight cut line IS the whole conic
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
