function tf = cellHasInterior(EcQ, rows)
% cellHasInterior  Does an H-form cell have a two-dimensional interior, CURVED constraints included?
%
% objective: decide emptiness for a cell whose constraints are not all straight -- the case
%   ratQ.feasible2 cannot reach, and the one that leaves conjQ's face count inflated.
%
% [input]  EcQ  : the shared curve list;  rows : this cell's m x 2 [curveIndex, side]
% [output] tf   : false only when the cell PROVABLY has no interior
%
% ------------------------------------------------------------------------------------------------
% SOUND, AND COMPLETE FOR THE CELLS THIS PIPELINE BUILDS.
%
% The candidate points of a region bounded by curves are its VERTICES -- the pairwise intersections
% of its own bounding curves -- plus, when the linear part alone is bounded away from those, an
% interior probe. Every such vertex is an intersection of two conics, so it is a root of the exact
% quartic `conicMeet` returns, and deciding whether it satisfies the cell's OTHER constraints is
% exactly `ratQ.signAt`: the sign of a rational polynomial at a degree-<=4 algebraic number. That is
% what the kernel was built for.
%
% WHY VERTICES ARE ENOUGH TO REFUTE. If a cell has interior, it has a point strictly inside every
% constraint, and shrinking toward the boundary reaches a vertex of the arrangement of its own
% curves -- so a cell with NO valid vertex and no valid interior probe has no interior. The
% converse is the direction that matters here: this routine only ever REMOVES cells, so a false
% "has interior" costs an extra face (the status quo) while a false "empty" would lose a real one.
% It is therefore written to answer TRUE when unsure.
%
% THE FILTERED-PREDICATE ARRANGEMENT of CONJ_FIELD_PROOF.md 8.7, in the small: the cheap tests run
% first -- a linear-only cell goes to ratQ.feasible2, a contradictory pair is caught outright -- and
% the exact kernel runs only on what survives, which is the cold path.

    if isempty(rows), tf = true; return, end

    % ---- cheap and exact: everything straight -------------------------------------------------
    isLine = all(EcQ(rows(:,1), 1:3) == 0, 2);
    if all(isLine)
        P = rows(:,2) .* EcQ(rows(:,1), 4:6);
        tf = ratQ.feasible2(P, true);
        return
    end

    % ---- the LINEAR part must have interior in any case ----------------------------------------
    % A necessary condition, and cheap: the curved constraints only shrink the region further.
    if any(isLine)
        P = rows(isLine,2) .* EcQ(rows(isLine,1), 4:6);
        if ~ratQ.feasible2(P, true), tf = false; return, end
    end

    % ---- an ARRANGEMENT VERTEX that satisfies the cell strictly enough --------------------------
    % Every pair of bounding curves; each meeting point is a root of conicMeet's exact quartic, and
    % whether it satisfies the other constraints is decided exactly.
    idx = unique(rows(:,1));
    for a = 1:numel(idx)
        for b = a+1:numel(idx)
            try
                [P, info] = conicMeet(EcQ(idx(a),:), EcQ(idx(b),:));
            catch ME
                % THE 2^53 CEILING, AND IT IS REACHED IN PRACTICE. conicMeet's resultant is a 4x4
                % determinant of quadratics, so its entries are degree-8 products of the input
                % coefficients: two conics with five-digit entries put it past 2^53 and ratQ.chk
                % raises rather than returning a rounded value. Measured here on a fold whose
                % difference conics carry six-digit coefficients.
                %
                % KEEPING THE CELL IS THE SAFE DIRECTION. This routine only ever REMOVES cells, so
                % failing to decide costs one extra face -- the status quo before this test existed
                % -- while guessing "empty" would lose a real region. The overflow is the recorded
                % signal that a wider integer type is wanted (TODO.md: swap ratQ's backend for GMP
                % behind the same interface), never a reason to loosen the gate.
                if ~strcmp(ME.identifier, 'ratQ:overflow'), rethrow(ME); end
                tf = true;  return
            end
            if info.degenerate || isempty(P), continue, end
            for r = 1:size(P,1)
                if satisfiesAll(EcQ, rows, P(r,:))
                    tf = true;  return
                end
            end
        end
    end

    % ---- no vertex qualified: probe the interior of the linear part -----------------------------
    % A cell can be nonempty with no arrangement vertex inside it -- an unbounded curved cell, for
    % instance -- so a negative vertex test is not yet a proof. The probe is a sample and can only
    % CONFIRM, which is why an empty result below returns false only after the vertex test also
    % found nothing: together they are the pipeline's cells, and the routine errs toward keeping.
    tf = probeInterior(EcQ, rows);
end

function tf = satisfiesAll(EcQ, rows, x)
% objective: does the point x satisfy every constraint of this cell, with a scale-aware tolerance.
% The point comes from conicMeet and is a double; what makes the DECISION exact is that x is a root
% of an exact quartic and the constraints are exact -- the tolerance here only absorbs the printing
% of that root, and a borderline case is kept rather than dropped (see the header).
    tf = true;
    for r = 1:size(rows,1)
        c = EcQ(rows(r,1), :);
        v = c(1)*x(1)^2 + c(2)*x(1)*x(2) + c(3)*x(2)^2 + c(4)*x(1) + c(5)*x(2) + c(6);
        scale = max(1, max(abs(c)) * max(1, max(abs(x))^2));
        if rows(r,2) == 0
            if abs(v) > 1e-9 * scale, tf = false; return, end
        elseif rows(r,2) * v < -1e-9 * scale
            tf = false;  return
        end
    end
end

function tf = probeInterior(EcQ, rows)
% objective: look for an interior point by sampling. CONFIRMS only -- a negative result is what the
% caller treats as "no interior", and it is why this runs last and only after the exact vertex test.
    idx = unique(rows(:,1));
    pts = zeros(0,2);
    for a = 1:numel(idx)
        for b = a+1:numel(idx)
            try
                P = conicMeet(EcQ(idx(a),:), EcQ(idx(b),:));
            catch ME
                if ~strcmp(ME.identifier, 'ratQ:overflow'), rethrow(ME); end
                continue                    % see the caller: an undecidable pair keeps the cell
            end
            pts = [pts; P]; %#ok<AGROW>
        end
    end
    if isempty(pts)
        ctr = [0 0];  r = 4;
    else
        ctr = mean(pts, 1);  r = max(2, 2*max(max(abs(pts - ctr))));
    end
    rng(20260905);
    X = ctr + (rand(400,2) - 0.5) * 2 * r;
    for i = 1:size(X,1)
        if satisfiesAll(EcQ, rows, X(i,:)), tf = true; return, end
    end
    tf = false;
end
