function inc_check()
% F(j,:) = [left, right] of edge j. Check it against what the CELLS say: the cell named on the left
% must carry a +1 side on that curve, the one on the right a -1, and no cell may be named on a side
% it does not claim.
    E3 = [1 2 1; 2 3 1; 3 1 1];  F3 = [1 0; 1 0; 1 0];
    sq = [0 0; 1 0; 1 1; 0 1];  Esq = [1 2 1; 2 3 1; 3 4 1; 4 1 1];  Fsq = [1 0;1 0;1 0;1 0];
    cases = {};
    cases{end+1} = {'fold', maxQ(conjQ(QuaPol([0 0; 1 0; 0 1], E3, [0 0 0 0 1 0 1 0 0 0], F3)), ...
                                 conjQ(QuaPol([1 0; 1 1; 0 1], E3, [0 0 0 0 4 1 3 -2 1 0], F3)))};
    cases{end+1} = {'square PD',  conjQ(QuaPol(sq, Esq, [0 0 0 0 2 0 2 1 -1 2], Fsq))};
    cases{end+1} = {'square indef', conjQ(QuaPol(sq, Esq, [0 0 0 0 1 0 -1 1 -1 2], Fsq))};
    cases{end+1} = {'oneNorm', conjQ(QuaPol.oneNorm())};

    for c = 1:numel(cases)
        nm = cases{c}{1};  g = cases{c}{2};
        bad = 0;  checked = 0;
        for j = 1:g.ne
            for col = 1:2
                k = g.F(j,col);
                if k == 0, continue, end
                checked = checked + 1;
                r = find(g.FC{k}(:,1) == j, 1);
                if isempty(r), bad = bad + 1; continue, end
                want = 3 - 2*col;                       % col 1 -> +1 (left), col 2 -> -1 (right)
                if g.FC{k}(r,2) ~= want, bad = bad + 1; end
            end
        end
        % and every vertex an edge names must actually name that edge back
        vbad = 0;
        for j = 1:g.ne
            for v = g.E(j,1:2)
                if v == 0, continue, end
                if ~any(g.Vname(v,1:2) == j), vbad = vbad + 1; end
            end
        end
        fprintf('%-14s ne=%-4d F entries %-4d inconsistent %d | E endpoints inconsistent %d\n', ...
            nm, g.ne, checked, bad, vbad);
    end
end
