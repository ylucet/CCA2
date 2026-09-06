function cellPairFilterProbe()
% ASSEMBLY, and the measurement that decides its architecture. A bridge lives on a PAIR of dual
% cells, one from each piece, so the natural loop is over all such pairs -- and that loop is the
% whole cost, because each pair costs an elimination (1 s to 20 s, measured in
% .claude/bridgeCellProbe.m).
%
% A single convex quadratic on an n-gon has about 2n+1 dual cells: one per vertex, one per edge, one
% interior. So two triangles give 7 x 7 = 49 pairs, and the wider fixtures in the corpus reach 34
% cells apiece, which is 1156. At 20 s each that is hours.
%
% THE QUESTION IS THEREFORE NOT HOW FAST ONE ELIMINATION IS -- item 3 already took the worst case
% from 64 s to 19 s -- BUT HOW MANY CAN BE SKIPPED WITHOUT RUNNING THEM. Two cheap exact filters
% are available before any algebra:
%
%   1. the two cells must MEET in dual space. Both are polyhedral (a convex piece's conjugate has
%      only straight cell boundaries -- conjQ's Case B), so this is one integer feasibility test.
%   2. on that intersection g1 - g2 must change sign, or there is no tie and no bridge.
%
% This measures filter 1, which is the cheap one, on real fixtures.

    fprintf('%-38s %7s %7s %9s %9s\n', 'piece pair', 'cellsA', 'cellsB', 'pairs', 'meeting');
    fprintf('%s\n', repmat('-', 1, 76));

    cases = {
        {'two unit triangles (shifted)', tri([0 0; 1 0; 0 1]),      tri([2 0; 3 0; 2 1])}, ...
        {'two unit triangles (overlapping cones)', tri([0 0; 1 0; 0 1]), tri([0 0; -1 0; 0 -1])}, ...
        {'triangle vs square',           tri([0 0; 1 0; 0 1]),      sq()}, ...
        {'square vs square (shifted)',   sq(),                      sqAt(3)} };

    for i = 1:numel(cases)
        nm = cases{i}{1};  fA = cases{i}{2};  fB = cases{i}{3};
        try
            gA = conjQ(fA);   gB = conjQ(fB);
        catch e
            fprintf('%-38s   %s\n', nm, e.identifier);
            continue
        end
        [np, nm2] = countMeeting(gA, gB);
        fprintf('%-38s %7d %7d %9d %9d\n', nm, gA.nf, gB.nf, np, nm2);
    end

    fprintf(['\nThe "meeting" column is what an assembly would actually eliminate over. Whatever it\n' ...
             'is, the loop must be written filter-first: the elimination is the expensive step and\n' ...
             'the feasibility test is integer arithmetic on a handful of half-planes.\n']);
end

% ------------------------------------------------------------------------------------------------
function [nPairs, nMeet] = countMeeting(gA, gB)
% How many cell pairs have a nonempty intersection in dual space?
    nPairs = gA.nf * gB.nf;
    nMeet  = 0;
    for i = 1:gA.nf
        RA = cellRows(gA, i);
        if isempty(RA), continue, end
        for j = 1:gB.nf
            RB = cellRows(gB, j);
            if isempty(RB), continue, end
            if ratQ.feasible2([RA; RB], false), nMeet = nMeet + 1; end
        end
    end
end

function rows = cellRows(g, k)
% A cell's half-planes as [a b c] meaning a*x + b*y + c >= 0. Every boundary of a convex piece's
% conjugate is a LINE (conjQ Case B), so the conic's quadratic part is zero and the linear part is
% the half-plane; a side of 0 is an equality and is skipped, since a thin cell meets nothing in the
% sense this filter is asking about.
    rows = zeros(0,3);
    FC = g.FC{k};
    for r = 1:size(FC,1)
        e = FC(r,1);  sgn = FC(r,2);
        c = g.EcQ(e,:);
        if any(c(1:3) ~= 0), rows = zeros(0,3); return, end   % a curved boundary: not this case
        if sgn == 0, continue, end
        rows(end+1,:) = sgn * c(4:6); %#ok<AGROW>
    end
end

% ------------------------------------------------------------------------------------------------
function f = tri(V)
    f = QuaPol(V, [1 2 1; 2 3 1; 3 1 1], [0 0 0 0, 2, 0, 2, 0, 0, 0], [1 0; 1 0; 1 0]);
end

function f = sq()
    f = sqAt(0);
end

function f = sqAt(dx)
    V = [0 0; 1 0; 1 1; 0 1];  V(:,1) = V(:,1) + dx;
    f = QuaPol(V, [1 2 1; 2 3 1; 3 4 1; 4 1 1], [0 0 0 0, 2, 0, 2, 0, 0, 0], ...
               [1 0; 1 0; 1 0; 1 0]);
end
