function [env, info] = envelopePrune(env)
% envelopePrune  Drop the envelope faces that can never be the maximum, EXACTLY.
%
% objective: `envelopeAssemble` returns one candidate per cell pair per candidate kind, and most of
%   them are beaten everywhere they apply -- 120 faces for a convex fixture whose envelope has a
%   handful. This removes that surplus, exactly: measured 120 -> 52, 142 -> 94 and 110 -> 38 on the
%   three standard fixtures, with every value and every finiteness verdict unchanged.
%
% IT DOES NOT MAKE `maxQ` THE WAY TO BUILD A MESH, and the reason is worth stating precisely because
% half of it was got wrong once already (DECISIONS.md 2026-09-18, and its correction):
%
%   * folding individual FACES is SEMANTICALLY WRONG. A `QuaCon` is +infinity off its faces, while a
%     candidate that does not apply contributes nothing -- -infinity. Measured: the fold collapses
%     to an empty mesh after about ten faces.
%   * folding a whole PAIR is semantically RIGHT -- a pair's contribution really is +infinity
%     outside its candidate regions -- and it EXPLODES: 5, 8, 15, 24, 48, 63, 252, 270, 675, 1323,
%     1917, 3726 faces over twelve pairs, against an envelope with a handful. The growth is
%     intrinsic: a max over the FIRST k pairs is a partial max, which genuinely has many pieces that
%     only the full max deletes.
%
% So the mesh means building the ARRANGEMENT once -- overlay the regions, decide the winner per cell
% -- and this routine's job is to make that input as small as it can exactly be.
%
% [input]  env : as `envelopeAssemble` returns it
% [output] env : the same, with the dominated quadratic faces removed
%          info : .before, .after, .duplicates, .dominated
%
% ==================================================================================================
% WHAT "DOMINATED" MEANS HERE, AND WHY THE TEST IS SOUND RATHER THAN COMPLETE
%
% Face i may be dropped when some other face j is at least as large everywhere face i applies AND
% applies everywhere face i does:
%
%       R_i  subset of  R_j       and       v_j >= v_i  on R_i .
%
% Then every value i could contribute to the max, j already contributes. Both halves are decided
% EXACTLY, and both are asked of the LINEAR RELAXATION of R_i -- its half-plane rows, with any conic
% rows dropped:
%
%   * dropping rows only makes the region BIGGER, so a property that holds on the relaxation holds
%     on R_i. The test is therefore SOUND in the direction that matters: it never drops a face that
%     could win. It is not complete -- a face dominated only because of a curved row survives -- and
%     that is the right way round, since a surviving face costs size and a dropped one costs
%     correctness.
%   * `ratQ.feasible2` decides the inclusion (is R_i with one of j's rows NEGATED still feasible)
%     and `ratQ.quadSignOnPoly` decides the comparison, both over the integers.
%
% THE VALUE DIFFERENCE IS A QUADRATIC, not a rational function, even though the faces carry
% denominators: v_j - v_i = (num_j den_i - num_i den_j)/(den_i den_j) and both denominators are
% POSITIVE integers, so the sign is the numerator's and the numerator is a quadratic form.
%
% BRIDGE AND CORNER FACES ARE NOT PRUNED. Their regions are predicates rather than row lists, so
% neither test applies to them; they are few (one bridge and a handful of corners per curved pair)
% and they are exactly the faces no quadratic can dominate anyway -- on their own regions they ARE
% the envelope (measured 2026-09-17).

    if nargin ~= 1
        error('envelopePrune:nargin', 'expects (env).');
    end
    if ~isstruct(env) || ~all(isfield(env, {'quad', 'bridge'}))
        error('envelopePrune:shape', 'env must be what envelopeAssemble returns.');
    end

    n0 = numel(env.quad);
    keep = true(1, n0);

    % ---- 1. identical faces ---------------------------------------------------------------------
    % The same candidate reached from two pairs is the same face twice. Bitwise, because the rows
    % are canonical integers -- the same argument assembleQuaConCells' deduplication rests on.
    dup = 0;
    for i = 1:n0
        if ~keep(i), continue, end
        for j = i+1:n0
            if ~keep(j), continue, end
            if isequal(env.quad(i).num, env.quad(j).num) && ...
               isequal(env.quad(i).den, env.quad(j).den) && ...
               isequal(env.quad(i).con, env.quad(j).con)
                env.quad(i).pair = union(env.quad(i).pair, env.quad(j).pair);
                keep(j) = false;  dup = dup + 1;
            end
        end
    end

    % ---- 2. domination --------------------------------------------------------------------------
    dom = 0;
    idx = find(keep);
    for a = 1:numel(idx)
        i = idx(a);
        if ~keep(i), continue, end
        Li = linearRows(env.quad(i).con);
        if isempty(Li), continue, end                % no linear rows: nothing can be decided here
        for b = 1:numel(idx)
            j = idx(b);
            if j == i || ~keep(j) || ~keep(i), continue, end
            if ~containsRegion(Li, env.quad(j).con), continue, end
            if ~atLeastAsLarge(env.quad(j), env.quad(i), Li), continue, end
            % THE DOMINATING FACE INHERITS THE PAIR, and this is not bookkeeping. A pair's
            % coverage is what says where its sup is FINITE, so dropping a face shrinks the
            % coverage and the envelope then claims +infinity on ground it owns -- measured, 60 of
            % 300 sampled points flipped from finite to infinite before this line existed. It is
            % sound precisely because the test above established R_i subset of R_j: everything i
            % covered for its pair, j covers.
            env.quad(j).pair = union(env.quad(j).pair, env.quad(i).pair);
            keep(i) = false;  dom = dom + 1;
        end
    end

    env.quad = env.quad(keep);
    info = struct('before', n0, 'after', numel(env.quad), 'duplicates', dup, 'dominated', dom);
end

% ==================================================================================================
function L = linearRows(con)
% The half-plane rows of a region, as [a b c] meaning a*x + b*y + c >= 0. A conic row is DROPPED,
% which enlarges the region -- see the header for why that is the safe direction.
    L = zeros(0,3);
    for r = 1:size(con,1)
        if any(con(r,1:3) ~= 0), continue, end
        if size(con,2) >= 8 && con(r,8) ~= 0, continue, end   % one line of a pair: not a half-plane
        L(end+1,:) = con(r,7) * con(r,4:6); %#ok<AGROW>
    end
end

% --------------------------------------------------------------------------------------------------
function tf = containsRegion(Li, conj)
% Does every row of face j hold on the region Li? Row by row: Li together with the row NEGATED must
% be infeasible. `ratQ.feasible2` decides that exactly, by Fourier-Motzkin.
    tf = true;
    for r = 1:size(conj,1)
        if any(conj(r,1:3) ~= 0), tf = false; return, end          % a curved row: not decided here
        if size(conj,2) >= 8 && conj(r,8) ~= 0, tf = false; return, end
        row = conj(r,7) * conj(r,4:6);
        if ratQ.feasible2([Li; -row], true)
            tf = false;  return                                    % j's row can fail inside R_i
        end
    end
end

% --------------------------------------------------------------------------------------------------
function tf = atLeastAsLarge(fj, fi, Li)
% Is v_j >= v_i on the region Li? The difference is a quadratic once both denominators are cleared,
% and `ratQ.quadSignOnPoly` decides its sign on a polygon exactly.
    num = ratQ.chkB(fj.num * fi.den - fi.num * fj.den, abs(fj.num) * abs(fi.den) + abs(fi.num) * abs(fj.den), 'value difference');
    if any(num(1:4) ~= 0)
        tf = false;  return                       % a cubic difference: outside this test's reach
    end
    c = [num(5), 2*num(6), num(7), 2*num(8), 2*num(9), 2*num(10)];
    if all(c == 0), tf = true; return, end        % the two faces are the same function
    % `quadSignOnPoly` answers for the CANONICAL spelling of the conic, whose leading entry it makes
    % positive -- so the answer is read back through that flip, exactly as a region's side is.
    nz = find(c ~= 0, 1);
    [s, ~] = ratQ.quadSignOnPoly(ratQ.conic(c), Li);
    tf = (s * sign(c(nz))) >= 1;                 % +1 means v_j - v_i >= 0 on the whole region
end
