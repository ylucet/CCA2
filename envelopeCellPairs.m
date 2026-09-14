function [pairs, stats] = envelopeCellPairs(gA, gB)
% envelopeCellPairs  Enumerate and CLASSIFY the cell pairs of two conjugates, exactly.
%
% objective: the first step of the convex envelope's assembly. `co f` over two pieces is
%
%       f**(x) = sup_s [ <s,x> - G(s) ] ,        G = max(gA, gB) ,
%
% and that sup splits over any cover of the s-plane on which G is a SINGLE quadratic. The cover is
% the overlay of the two conjugates' own cells: over P_ij = cell_i(gA) ^ cell_j(gB) the max is
% max(gA_i, gB_j), and the tie conic q = gA_i - gB_j either keeps one strict sign on P_ij -- one
% quadratic, a NON-BRIDGE pair -- or changes sign, and then P_ij splits into {q >= 0} carrying
% gA_i and {q <= 0} carrying gB_j, whose common boundary is the arc that carries the BRIDGE face.
%
% This routine answers which of those each pair is, and it answers it with INTEGER PREDICATES: the
% classification decides the MESH, and `CONJ_FIELD_PROOF.md` 8.0 puts mesh decisions on the exact
% side of the line. Nothing here is a tolerance.
%
% [input]  gA, gB : QuaCon, the exact conjugates of the two pieces (conjQ of each).
% [output] pairs  : struct array, one entry per pair that MEETS in 2-D:
%                     .kind  'poly' (one quadratic wins) | 'bridge' (the tie changes sign)
%                     .ij    [i j], the two cell indices
%                     .P     the pair's half-planes, integer rows [a b c] meaning a*s1+b*s2+c >= 0
%                     .q     the tie conic as integer [a b c d e f], all-zero when the two faces
%                            are IDENTICAL on the pair
%                     .win   for 'poly': 'A' or 'B', which conjugate's face wins on the whole pair;
%                            for 'bridge': '' (both, one on each side of q)
%          stats  : counts, for the measurements the probes report -- nPairs (examined), nMeet
%                   (meet in 2-D), nSame (identical faces), nPoly, nBridge, nSkipCurved.
%
% ------------------------------------------------------------------------------------------------
% THE TWO FILTERS, AND WHY BOTH ARE NEEDED.
%
%   FILTER 1, `ratQ.feasible2(P, true)` -- do the two cells MEET IN 2-D? A pair meeting only along
%     an edge or at a point contributes no region and must not produce a face. Measured
%     (DECISIONS.md 2026-09-06): of 49 pairs for two triangles, 24 meet at all and only 13 meet in
%     2-D, so this is where most of the quadratic pair count goes.
%
%   FILTER 2, `ratQ.quadSignOnPoly(q, P)` -- does the tie conic change sign on the pair? This is
%     the one that decides bridge from non-bridge, and it is the whole classification: no sign
%     change means one of the two quadratics dominates on the entire pair and the envelope needs no
%     bridge there.
%
% THE IDENTICAL-FACE CASE IS COUNTED APART, and that is a correction rather than bookkeeping. When
% both primal pieces share the vertex a pair's cells are built from, the two conjugates agree on
% the whole pair: q is identically zero, there is no tie curve, and there is no bridge. Calling it
% a bridge is the mistake `.claude/bridgeEnvelopeProbe.m` made on its first run, and an all-zero q
% reaching `quadSignOnPoly` would be asked which side of nothing a point is on.
%
% WHAT IS OUT OF SCOPE, BY NAME. A cell whose boundary contains a genuine CONIC edge is refused
% with `envelopeCellPairs:curvedCell` rather than dropped. For the conjugate of a CONVEX piece
% every cell is polyhedral (`CONJ_FIELD_PROOF.md` 7.3: the subdivision is the KKT active set, so
% every boundary is the affine condition that a multiplier changes sign), so a curved cell here
% means the caller passed something other than a convex piece's conjugate -- a caller defect, and
% one that would otherwise show up much later as a wrong region.

    if nargin ~= 2
        error('envelopeCellPairs:nargin', 'expects (gA, gB), two QuaCon conjugates.');
    end
    for g = {gA, gB}
        if ~isa(g{1}, 'QuaCon')
            error('envelopeCellPairs:type', ...
                ['both arguments must be QuaCon (the exact conjugate of one piece); got %s. ' ...
                 'The classification below is integer arithmetic on fN/fD and EcQ, which only ' ...
                 'an exact conjugate carries.'], class(g{1}));
        end
    end

    pairs = struct('kind', {}, 'ij', {}, 'P', {}, 'q', {}, 'win', {});
    stats = struct('nPairs', 0, 'nMeet', 0, 'nSame', 0, 'nPoly', 0, 'nBridge', 0, 'nSkipCurved', 0);

    for i = 1:gA.nf
        RA = cellHalfPlanes(gA, i);
        if isempty(RA), continue, end
        for j = 1:gB.nf
            RB = cellHalfPlanes(gB, j);
            if isempty(RB), continue, end
            stats.nPairs = stats.nPairs + 1;

            P = [RA; RB];
            if ~ratQ.feasible2(P, true), continue, end        % FILTER 1: meet in 2-D
            stats.nMeet = stats.nMeet + 1;

            q = tieConic(gA, i, gB, j);
            if all(q == 0)
                % IDENTICAL faces on this pair -- see the header. One quadratic, no tie, no bridge.
                stats.nSame = stats.nSame + 1;
                stats.nPoly = stats.nPoly + 1;
                pairs(end+1) = mk('poly', [i j], P, q, 'A'); %#ok<AGROW>
                continue
            end

            [sgn, strict] = ratQ.quadSignOnPoly(q, P);        % FILTER 2: does the tie change sign?
            if sgn ~= 0 && strict
                stats.nPoly = stats.nPoly + 1;
                if sgn > 0, w = 'A'; else, w = 'B'; end
                pairs(end+1) = mk('poly', [i j], P, q, w); %#ok<AGROW>
            else
                stats.nBridge = stats.nBridge + 1;
                pairs(end+1) = mk('bridge', [i j], P, q, ''); %#ok<AGROW>
            end
        end
    end
end

% ------------------------------------------------------------------------------------------------
function p = mk(kind, ij, P, q, win)
    p = struct('kind', kind, 'ij', ij, 'P', P, 'q', q, 'win', win);
end

function rows = cellHalfPlanes(g, k)
% Cell k's boundary as INTEGER half-planes [a b c] meaning a*s1 + b*s2 + c >= 0.
%
% A QuaCon stores a face as a list of [edge, side] pairs against its exact conics `EcQ`, and for a
% polyhedral cell each such conic is a line -- its quadratic part is zero and rows 4:6 are the line.
% A side of 0 means the face is ON the curve, i.e. THIN: it bounds no 2-D region, so it carries no
% half-plane and is skipped here rather than turned into a degenerate row.
    FC = g.FC{k};
    rows = zeros(0,3);
    for t = 1:size(FC,1)
        cc = g.EcQ(FC(t,1), :);
        if any(cc(1:3) ~= 0)
            error('envelopeCellPairs:curvedCell', ...
                ['cell %d of a %s has a genuinely CONIC edge (edge %d, quadratic part [%d %d %d]). ' ...
                 'The conjugate of a CONVEX piece is polyhedral -- every cell boundary is the ' ...
                 'affine condition that a KKT multiplier changes sign -- so a curved cell means ' ...
                 'this is not one.'], k, class(g), FC(t,1), cc(1), cc(2), cc(3));
        end
        if FC(t,2) == 0, continue, end          % a THIN face: on the curve, bounding nothing
        rows(end+1,:) = FC(t,2) * cc(4:6); %#ok<AGROW>
    end
end

function q = tieConic(gA, i, gB, j)
% The locus where the two cells' faces agree, cleared to integers.
%
% THE WEIGHTED BASIS IS THE TRAP HERE. A face is stored in the 10-wide cubic basis whose squared
% terms carry halves, so the conic of the DIFFERENCE is [n5, 2*n6, n7, 2*n8, 2*n9, 2*n10] and not
% the raw entries. Getting it wrong names a different curve and moves real numbers -- measured and
% recorded, DECISIONS.md 2026-09-06.
    D = gA.fN(i,:) * gB.fD(j) - gB.fN(j,:) * gA.fD(i);
    q = [D(5), 2*D(6), D(7), 2*D(8), 2*D(9), 2*D(10)];
end
