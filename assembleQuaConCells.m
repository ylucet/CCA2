function g = assembleQuaConCells(cells)
% assembleQuaConCells  Turn a list of H-form cells into a QuaCon.
%
% objective: deduplicate the bounding curves into one canonical edge list, drop the cells that
%   describe the empty set, name the corners, and reduce each face function -- the last step every
%   exact producer (conjQ's cases, maxQ) shares.
%
% [input]  cells : struct array with fields
%                    num, den : the face function, as an integer 1x10 numerator over a denominator
%                    con      : k x 7 rows [a b c d e f  sign], meaning the cell lies where
%                               sign * (a x^2 + b xy + c y^2 + d x + e y + f) >= 0
% [output] g     : QuaCon
%
% DEDUPLICATION IS BITWISE, and that is the whole reason the conic rows are canonical. Two cells
% that share a facet reached it by different routes; in doubles those two spellings differ by an
% ULP and the shared facet becomes invisible to `merge`. DECISIONS.md 2026-08-17 measured the
% consequence -- 57 cells carrying 10 distinct functions, 4 merges out of 612 attempts, and a cell
% count that then grew without bound. Here `find(all(EcQ == c, 2), 1)` is the entire comparison.
%
% ------------------------------------------------------------------------------------------------
% TWO THINGS ARE DELIBERATELY INCOMPLETE, AND BOTH ARE STATED RATHER THAN PAPERED OVER.
%
% 1. THE EMPTINESS TEST IS SOUND, NOT COMPLETE. A cell is dropped when its LINEAR constraints alone
%    are already infeasible, which ratQ.feasible2 decides exactly. A cell that is empty only
%    because of a CURVED constraint survives. That is an over-approximation of the face set, and it
%    is harmless where it matters -- `eval` never reports such a face, because no point satisfies
%    its constraints -- but `nf` is then an upper bound rather than the count of faces that carry
%    area. Deciding it properly means asking whether a conic meets a polygon's interior, which is
%    Phase 2c's filtered predicate with the exact degree-4 kernel behind it, and that is not built.
%
% 2. ONLY LINE-LINE CORNERS ARE NAMED. Two lines meet in a rational point, so membership in the
%    cell is an exact integer sign test and the corner can be named with certainty. A corner
%    involving a curved edge is generically irrational (degree up to 4 over Q -- CONJ_FIELD_PROOF.md
%    Theorem 1), so deciding whether it lies in the cell needs the same exact kernel as (1).
%    Naming it without that test would put points into the vertex list that bound nothing, which is
%    worse than a vertex list that is honestly partial -- E and F are empty for the same reason.
%
% 3. THE MERGE IS PARTIAL TOO. mergeAdjacentCells joins two cells that carry one function and are
%    separated by a single curve, which is an exact set identity and is the pattern a pairwise fold
%    generates. Two cells carrying one function that differ in more than one constraint are left
%    apart; joining those needs redundancy elimination against a curved constraint set, i.e. (1)'s
%    kernel again.
%
% Neither gap can produce a wrong VALUE: `eval` reads the face functions and the sign conditions,
% and both are exact.

    % NO CELLS IS AN ANSWER, NOT A FAILURE. A QuaCon with zero faces evaluates to +infinity at
    % every point, which IS the function f = +inf -- the correct conjugate whenever the sup
    % diverges for every s. Raising here made the routines above decline to return something this
    % type can already store; verified directly that the empty object builds and evaluates.
    if isempty(cells)
        g = QuaCon(zeros(0,3), zeros(0,3), zeros(0,6), zeros(0,10), zeros(0,1), zeros(0,2), {});
        return
    end

    % ---- drop the cells that provably carry no two-dimensional face -------------------------
    % Two sound tests, both exact, neither complete (see note 1 in the header):
    %
    %  (a) CONTRADICTORY SIDES. A cell asking for both c >= 0 and c <= 0 on the SAME canonical
    %      curve lies inside {c = 0}, which has empty interior whatever c is -- curve or line. A
    %      nested fold produces these in quantity, because the same pair of functions can be
    %      compared again at a later fold step, and this is the only filter that sees a CURVED
    %      constraint at all. Cheap: the rows are canonical integers, so it is an equality test.
    %
    %  (b) INFEASIBLE LINEAR PART. Everything straight, decided exactly by Fourier-Motzkin.
    %
    %  (c) A CONSTANT-SIGN CONIC. A conic form that is nonnegative (or nonpositive) on the whole
    %      plane -- ratQ.conicSign decides this exactly, by testing the 3x3 form for semi-
    %      definiteness -- makes its condition either VACUOUS, in which case the row is dropped and
    %      the remaining tests get sharper, or satisfiable only ON the curve, in which case the
    %      cell has empty interior. This is the only handle on CURVED emptiness that does not need
    %      the degree-4 kernel, and it is what a nested fold generates most of.
    live = true(1, numel(cells));
    for k = 1:numel(cells)
        rows = cells(k).con;

        % A THIN cell -- one carrying an EQUALITY side -- has no interior BY CONSTRUCTION, so the
        % nonempty-interior tests below would delete exactly the faces that case exists to build.
        % Only the contradictory-sides test still applies to it, and that one is about +1 against
        % -1 on one curve, which a 0 does not participate in.
        if any(rows(:,7) == 0)
            [~, ~, ict] = unique(rows(:,1:6), 'rows');
            for u = 1:max(ict)
                sd = rows(ict == u, 7);
                if any(sd > 0) && any(sd < 0), live(k) = false; break, end
            end
            continue
        end

        % (c) first, because it can DELETE rows and so make (a) and (b) sharper. A conic that
        %     takes one sign everywhere gives a condition that is either vacuous (drop the row,
        %     it constrains nothing) or unsatisfiable off the curve itself (drop the cell).
        keep = true(size(rows,1), 1);
        for r = 1:size(rows,1)
            cs = ratQ.conicSign(rows(r,1:6));
            if cs == 0, continue, end
            if cs == rows(r,7)
                keep(r) = false;                               % vacuous
            else
                live(k) = false;  break                        % empty interior
            end
        end
        if ~live(k), continue, end
        rows = rows(keep,:);
        cells(k).con = rows;
        if isempty(rows), continue, end                        % no conditions: the whole plane

        [~, ~, ic] = unique(rows(:,1:6), 'rows');
        for u = 1:max(ic)
            sides = rows(ic == u, 7);
            if any(sides > 0) && any(sides < 0)
                live(k) = false;  break                        % (a)
            end
        end
        if ~live(k), continue, end
        isLine = all(rows(:,1:3) == 0, 2);
        live(k) = ratQ.feasible2(rows(isLine,7) .* rows(isLine,4:6), true);   % (b)
    end
    cells = cells(live);
    if isempty(cells)
        % Everything filtered away, and the filter demands a two-dimensional INTERIOR. The
        % zero-face object -- which evaluates to +infinity everywhere -- is returned, and that is
        % right whenever the domain really is empty.
        %
        % KNOWN LIMITATION, with a fixture: a conjugate whose domain has EMPTY INTERIOR is reported
        % as empty too. QuaPol.examples{19} is nine AFFINE pieces whose dual domain is the single
        % point s = (0,0), where f*(0,0) = -inf f = 0; conjQ answers +infinity there. Relaxing the
        % filter to "nonempty as a set" was tried and is WORSE: it kept 982 degenerate cells, and
        % eval's tolerance then admitted points that are genuinely outside, turning one wrong point
        % into two. The real fix is to detect a thin dual domain and emit it with EQUALITY sides,
        % as caseAFullDomain already does for the point and line cases -- see TODO.md 2026-09-05.
        g = QuaCon(zeros(0,3), zeros(0,3), zeros(0,6), zeros(0,10), zeros(0,1), zeros(0,2), {});
        return
    end

    % ---- one canonical curve list, shared by every cell that borders it ----------------------
    EcQ = zeros(0,6);
    FC  = cell(numel(cells), 1);
    for k = 1:numel(cells)
        rows = cells(k).con;
        idx  = zeros(size(rows,1), 1);
        for r = 1:size(rows,1)
            c = rows(r, 1:6);
            hit = find(all(EcQ == c, 2), 1);
            if isempty(hit)
                EcQ(end+1, :) = c; %#ok<AGROW>
                hit = size(EcQ, 1);
            end
            idx(r) = hit;
        end
        FC{k} = [idx, rows(:,7)];
    end

    fN = zeros(numel(cells), 10);  fD = ones(numel(cells), 1);
    for k = 1:numel(cells)
        [fN(k,:), fD(k)] = ratQ.canon(cells(k).num, cells(k).den);
    end

    % ---- shrink each cell to the constraints it actually needs, then merge --------------------
    % The order matters and was measured: without the redundancy pass the merge finds NOTHING,
    % because the fold leaves cells differing by one side AND one implied row, which does not match
    % the merge's "same list, opposite side" pattern until the implied row is gone.
    FC = dropRedundantRows(EcQ, FC);
    [fN, fD, FC] = mergeAdjacentCells(fN, fD, FC);

    % Some curves may now be referenced by nothing; drop them so `ne` counts the curves that
    % actually bound something, and renumber what is left.
    used = false(size(EcQ,1), 1);
    for k = 1:numel(FC), used(FC{k}(:,1)) = true; end
    renum = cumsum(used);
    EcQ = EcQ(used, :);
    for k = 1:numel(FC)
        if ~isempty(FC{k}), FC{k}(:,1) = renum(FC{k}(:,1)); end
    end

    % ---- name the corners that can be decided exactly ----------------------------------------
    ne = size(EcQ, 1);
    seen = false(ne, ne);
    Vname = zeros(0,3);
    % numel(FC), NOT numel(cells): mergeAdjacentCells above REMOVES cells, so the two counts differ
    % as soon as a merge actually fires. It never did until Case D was restructured, and then this
    % ran off the end of FC -- MATLAB:badsubscript on a four-piece input, an unnamed crash rather
    % than a named refusal.
    for k = 1:numel(FC)
        rows = FC{k};
        for p1 = 1:size(rows,1)
            for p2 = p1+1:size(rows,1)
                i = rows(p1,1);  j = rows(p2,1);
                if i == j, continue, end
                if any(EcQ(i,1:3) ~= 0) || any(EcQ(j,1:3) ~= 0), continue, end   % see note 2
                A = [EcQ(i,4) EcQ(i,5); EcQ(j,4) EcQ(j,5)];
                if ratQ.detExact(A) == 0, continue, end       % parallel: they meet nowhere
                [xn, xd] = ratQ.solve2(A, [-EcQ(i,6); -EcQ(j,6)]);
                if ~cellHoldsAt(EcQ, rows, xn, xd), continue, end
                lo = min(i,j);  hi = max(i,j);
                if ~seen(lo,hi)
                    seen(lo,hi) = true;
                    Vname(end+1, :) = [lo hi 1]; %#ok<AGROW>
                end
            end
        end
    end

    % E and F stay EMPTY: recovering which cell borders which along which segment is an
    % arrangement computation nothing consumes yet, and a plausible-looking invented incidence
    % would be worse than an honest gap.
    g = QuaCon(Vname, zeros(ne,3), EcQ, fN, fD, zeros(ne,2), FC);
end

function tf = cellHoldsAt(EcQ, rows, xn, xd)
% objective: does the point xn/xd satisfy EVERY constraint of this cell. Exact.
%
% xd > 0 (ratQ.canon normalises it), so multiplying a constraint through by xd^2 preserves its
% sign, and the quadratic terms then need xn's own square -- which keeps this exact for a CURVED
% constraint too, even though the corner it is testing must itself be a line-line one.
    tf = true;
    for r = 1:size(rows,1)
        c = EcQ(rows(r,1), :);
        val = ratQ.chk(c(1)*xn(1)^2 + c(2)*xn(1)*xn(2) + c(3)*xn(2)^2 ...
                     + c(4)*xn(1)*xd + c(5)*xn(2)*xd + c(6)*xd^2, 'constraint at corner');
        if rows(r,2) == 0
            if val ~= 0, tf = false; return, end      % an EQUALITY row: the corner must be ON it
        elseif rows(r,2) * val < 0
            tf = false; return
        end
    end
end
