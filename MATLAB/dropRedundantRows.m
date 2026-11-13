function FC = dropRedundantRows(EcQ, FC)
% dropRedundantRows  Remove constraints a cell does not need, exactly.
%
% objective: delete rows that the cell's OTHER rows already imply, so that two cells describing
%   genuinely adjacent regions end up with comparable lists.
%
% [input]  EcQ : the shared curve list; FC : cell arrays of [curveIndex, sign] rows
% [output] FC  : the same cells, with implied rows removed
%
% WHY THIS RUNS BEFORE THE MERGE, and it is the reason the merge found nothing without it.
% Measured 2026-09-04 on a PSD-singular fixture: two cells carrying ONE function, one with 15
% constraints and one with 14, whose curve sets differed by a single extra row and whose shared
% rows differed in exactly one SIDE. That is a mergeable pair wearing a disguise -- the extra row
% is implied by the others, so the two cells really are "the same list, opposite sides on one
% curve". A pairwise fold generates these in bulk: every fold step appends a comparison against a
% function that some earlier constraint has already settled.
%
% THE TEST IS EXACT AND SOUND, NOT COMPLETE. Row i is dropped when the other LINEAR rows alone make
% it impossible to violate -- i.e. {other linear rows} together with the strict negation of row i
% has empty interior, which ratQ.feasible2 decides exactly. Using only the linear rows is what
% makes it sound in the presence of curved ones: if the linear rows alone already imply row i, then
% so does the full list, whatever the curves do. The converse is not tested, so a row implied only
% WITH the help of a curve survives -- that needs Phase 2c's kernel, like everything else curved.
%
% Only LINEAR rows are candidates for removal. A curved row could also be redundant, and deciding
% that is the same missing kernel.

    for k = 1:numel(FC)
        rows = FC{k};
        if size(rows,1) < 2, continue, end

        % A THIN cell -- one carrying an EQUALITY side -- is left alone. This pass reasons about
        % HALF-PLANES: it asks whether the others make a row impossible to violate, and an equality
        % is not a half-plane. Passing one through anyway turns it into the all-zero row 0*c, which
        % feasible2 correctly reads as having no interior, so every equality was deleted and the
        % thin face silently became the whole plane. Measured while building the point-domain case:
        % ne came out 0 and the conjugate of an affine function evaluated to -kappa everywhere.
        if any(rows(:,2) == 0), continue, end

        keep = true(size(rows,1), 1);
        for i = 1:size(rows,1)
            if any(EcQ(rows(i,1), 1:3) ~= 0), continue, end     % curved: not a candidate

            % the other rows that are still kept AND are linear, as half-planes
            others = zeros(0,3);
            for j = 1:size(rows,1)
                if j == i || ~keep(j), continue, end
                if any(EcQ(rows(j,1), 1:3) ~= 0), continue, end
                others(end+1,:) = rows(j,2) * EcQ(rows(j,1), 4:6); %#ok<AGROW>
            end

            % Can row i be strictly violated while the others hold? If not, it is implied.
            violated = -rows(i,2) * EcQ(rows(i,1), 4:6);
            if ~ratQ.feasible2([others; violated], true)
                keep(i) = false;
            end
        end
        FC{k} = rows(keep, :);
    end
end
