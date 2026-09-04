function [fN, fD, FC] = mergeAdjacentCells(fN, fD, FC)
% objective: merge cells that carry the SAME function and are separated only by one curve.
% [input]/[output] fN, fD : the face functions; FC : cell arrays of [curveIndex, sign] rows
%
% THE RULE, and it is an exact set identity rather than a heuristic. If two cells carry the same
% function and their constraint lists are identical except for ONE curve, on which they take
% OPPOSITE sides, then
%       {A and c >= 0}  union  {A and c <= 0}  =  {A}
% because every point satisfies c >= 0 or c <= 0. So the two cells are replaced by one carrying the
% list without that curve. Nothing is approximated and no convexity is required of the union: this
% is the same identity that makes Quine-McCluskey's consensus step sound.
%
% WHY IT IS NEEDED. Step 3 folds pairwise, so each fold re-splits cells that the previous fold had
% already separated, and the same function comes out spread over many cells. Measured 2026-09-04 on
% a PSD-singular fixture: 75 occupied cells carrying 10 distinct functions. That is a merge
% problem, and it is distinct from the emptiness problem -- these cells are all NONEMPTY, they are
% just cut up for no reason.
%
% WHY IT IS NOT THE WHOLE ANSWER. Two cells can carry one function, be genuinely adjacent, and
% still differ in more than one constraint -- their union is then not describable by dropping a
% single row, and merging them needs redundancy elimination against a CURVED constraint set, which
% is Phase 2c. This handles the case that the fold actually generates, and stops there.
%
% The loop runs to a fixpoint: one merge shortens a list, which can expose another.
    for k = 1:numel(FC)
        FC{k} = normalizeRows(FC{k});
    end

    changed = true;
    while changed
        changed = false;
        % Group by exact function, so only cells that could possibly merge are compared.
        [~, ~, grp] = unique([fN, fD], 'rows');
        for g = 1:max(grp)
            idx = find(grp == g);
            if numel(idx) < 2, continue, end
            done = false;
            for a = 1:numel(idx)
                for b = a+1:numel(idx)
                    ka = idx(a);  kb = idx(b);
                    c = mergeableOn(FC{ka}, FC{kb});
                    if isempty(c), continue, end
                    rows = FC{ka};
                    FC{ka} = normalizeRows(rows(rows(:,1) ~= c, :));
                    FC{kb} = [];                       % marked dead; removed below
                    fD(kb) = NaN;
                    changed = true;  done = true;  break
                end
                if done, break, end
            end
        end
        alive = ~isnan(fD);
        fN = fN(alive, :);  fD = fD(alive);  FC = FC(alive);
    end
end

function c = mergeableOn(A, B)
% objective: the single curve index on which A and B differ by SIDE ALONE, or [] if there is none.
% Both lists arrive normalized (unique rows, sorted), so this is a comparison of sorted integers.
    c = [];
    if size(A,1) ~= size(B,1), return, end
    d = find(any(A ~= B, 2));
    if numel(d) ~= 1, return, end
    if A(d,1) ~= B(d,1) || A(d,2) == B(d,2), return, end
    c = A(d,1);
end

function R = normalizeRows(R)
% objective: one canonical spelling of a constraint list -- duplicate rows removed, rows sorted.
% Without this two lists that describe the same cell could differ only in the ORDER they were
% built in, and the merge test above, which compares row by row, would miss every such pair.
    if isempty(R), R = zeros(0,2); return, end
    R = unique(R, 'rows');
end
