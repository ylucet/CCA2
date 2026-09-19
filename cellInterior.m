function [tf, info] = cellInterior(C, rows)
% cellInterior  Does an arrangement cell have a two-dimensional interior, when some of its edges are
%   BRANCH rows -- one named line of a rational conjugate pair?
%
% objective: an arrangement is only as small as its emptiness test. Measured on the c = -2 fixture
%   (DECISIONS.md 2026-09-18), the same 23 curves give 31104 live cells if the branch rows are
%   ignored, 8072 if the undecided cells are kept, and 486 if the zero-margin cells are decided --
%   and the 8 corner curves then add no cells at all. The cost of a weak test is not its per-call
%   price but the cells it keeps, because each one is split again by every later curve.
%
% [input]  C    : the curve table, k x 7 -- [conic(1:6), branch], branch 0 for an ordinary curve
%          rows : this cell's m x 2 [curveIndex, side]
% [output] tf   : false only when the cell PROVABLY has no interior
%          info : .path ('rational' | 'margin' | 'exactTie' | 'unsureTie'), .margin
%
% ==================================================================================================
% A FILTERED PREDICATE, in three stages, cheapest first.
%
% 1. THE RATIONAL ROWS, by the tests that already exist: `ratQ.feasible2` when they are all straight,
%    `cellHasInterior`'s exact kernel when a conic is among them. A branch row cannot go to either --
%    its line is irrational, and the conic it names is the PAIR, whose sign is the product of the two
%    lines and not the one wanted.
%
% 2. THE MARGIN, in doubles. Every straight row, branch or not, is a half-plane, and the largest
%    inscribed margin of a set of half-planes is positive exactly when they have interior. In two
%    dimensions the optimum makes three constraints tight, so the triples are enumerated and no
%    solver is needed. A margin clear of zero decides the cell, and that is almost every cell.
%
% 3. THE TIE, exactly. A margin of zero means no interior -- but in doubles it cannot be told from a
%    genuinely thin cell, so it is decided by an exact Farkas certificate: the cell has NO interior
%    iff some lambda >= 0, not all zero, has sum(lambda_i a_i) = 0 and sum(lambda_i c_i) <= 0. Every
%    PAIR and TRIPLE of the straight rows is tried, because in the plane an infeasible system always
%    has an infeasible subsystem of at most three -- and because the LP's tight set alone is not
%    enough: two rows always have interior between them, so a tight PAIR reads as "live" while a
%    third row contradicts it.
%
%    The certificate's arithmetic lives in Q(t) for the branch rows' slope or level t, so when the
%    branch rows of a subsystem name the SAME conic there is one field, `ratQ.qsMul` multiplies in it
%    and `ratQ.signAt` decides the signs -- exactly, with no square root formed. A subsystem mixing
%    two different pairs would need their compositum; that is not built, so the cell is KEPT and the
%    path reported as `unsureTie`. Keeping costs a cell; dropping one wrongly loses a face.

    if nargin ~= 2
        error('cellInterior:nargin', 'expects (C, rows).');
    end
    if size(C, 2) ~= 7
        error('cellInterior:shape', ...
              'C must be k x 7 -- [conic(1:6), branch]; got %d columns.', size(C,2));
    end
    info = struct('path', 'rational', 'margin', inf);
    if isempty(rows), tf = true; return, end

    % ---- 1. the rational rows, exactly ------------------------------------------------------------
    lin = zeros(0,3);  curved = zeros(0,6);  sides = [];  H = zeros(0,3);  src = [];
    for t = 1:size(rows,1)
        c = C(rows(t,1), 1:6);  br = C(rows(t,1), 7);  sd = rows(t,2);
        if all(c(1:5) == 0)
            error('cellInterior:constant', ...
                  ['row %d is a CONSTANT, not a curve: %s. Such a row is vacuous when it holds ' ...
                   'and empties the cell when it does not, and must be resolved before the ' ...
                   'arrangement is built.'], t, mat2str(c));
        end
        if br ~= 0
            H(end+1,:) = sd * branchLine(c, br); %#ok<AGROW>
            src(end+1) = t; %#ok<AGROW>
        elseif all(c(1:3) == 0)
            lin(end+1,:) = sd * c(4:6); %#ok<AGROW>
            H(end+1,:) = sd * c(4:6); %#ok<AGROW>
            src(end+1) = t; %#ok<AGROW>
        else
            curved(end+1,:) = c; %#ok<AGROW>
            sides(end+1) = sd; %#ok<AGROW>
        end
    end
    if ~isempty(lin) && ~ratQ.feasible2(lin, true), tf = false; return, end
    if ~isempty(curved)
        tf = cellHasInterior([curved; zeros(size(lin,1),3), lin], ...
                             [(1:size(curved,1)).', sides(:); ...
                              (size(curved,1)+(1:size(lin,1))).', ones(size(lin,1),1)]);
        if ~tf, return, end
    end

    % ---- 2. the margin, in doubles ----------------------------------------------------------------
    tf = true;
    if size(H,1) < 2, return, end
    m = marginLP(H);
    info.margin = m;  info.path = 'margin';
    tol = 1e-9 * max(1, max(abs(H(:))));
    if m > tol, return, end
    if m < -tol, tf = false; return, end

    % ---- 3. the tie, exactly ----------------------------------------------------------------------
    info.path = 'exactTie';
    [empty, skipped] = certificate(C, rows, src);
    if empty, tf = false; return, end
    if skipped
        info.path = 'unsureTie';                    % kept, because keeping only costs a cell
    end
end

% ==================================================================================================
function [empty, skipped] = certificate(C, rows, src)
% Is the cell's straight part EMPTY of interior, by an exact Farkas certificate?
%
%       no interior  <=>  exists lambda >= 0, not all zero, with
%                         sum(lambda_i a_i) = 0  and  sum(lambda_i c_i) <= 0 .
%
% In the plane such a certificate never needs more than THREE constraints -- an infeasible system of
% half-planes has an infeasible subsystem of at most dimension plus one -- so every pair and triple
% is tried. The margin LP's tight set is NOT enough on its own: two rows of a cell always have
% interior between them, so reading only a tight PAIR concludes "live" while a third row contradicts
% it. That was the defect this replaced, and the six-sector test is what caught it.
%
% A subsystem whose branch rows name two DIFFERENT pairs would need the compositum of two quadratic
% fields, which is not built here; those are skipped, and the caller then KEEPS the cell -- keeping
% costs a cell, dropping one wrongly loses a face.
    empty = false;  skipped = false;
    m = numel(src);
    if m < 2, return, end
    for i = 1:m-1
        for j = i+1:m
            [e, sk] = tryOne(C, rows, src([i j]));
            empty = empty || e;  skipped = skipped || sk;
            if empty, return, end
            for k = j+1:m
                [e, sk] = tryOne(C, rows, src([i j k]));
                empty = empty || e;  skipped = skipped || sk;
                if empty, return, end
            end
        end
    end
end

% --------------------------------------------------------------------------------------------------
function [empty, skipped] = tryOne(C, rows, idx)
% One subsystem, in the single quadratic field its branch rows share.
    empty = false;  skipped = false;
    conics = zeros(0,6);
    for t = idx(:).'
        if C(rows(t,1), 7) ~= 0, conics(end+1,:) = C(rows(t,1), 1:6); end %#ok<AGROW>
    end
    quad = [];  iv = [-1 1];  ref = 0;
    if ~isempty(conics)
        if size(unique(conics, 'rows'), 1) > 1, skipped = true;  return, end
        pd = ratQ.pairLineData(conics(1,:));
        quad = pd.quad;
        % ONE generator for the whole subsystem. Two branch rows of the same pair name the TWO
        % DIFFERENT roots of this quadratic, so writing both as "t" makes them the same line and the
        % certificate can never be found -- which is precisely what the six-sector test caught.
        % The second root is rational in the first: t + t' = -be/al, so t' = (-be - al t)/al.
        for t = idx(:).'
            if C(rows(t,1), 7) ~= 0, ref = C(rows(t,1), 7);  break, end
        end
        if numel(quad) >= 3, iv = pd.iv(ref,:); end
    end
    A1 = {};  A2 = {};  Cc = {};
    for t = idx(:).'
        [a, c, ok] = exactHalfPlane(C, rows, t, ref);
        if ~ok, skipped = true;  return, end
        A1{end+1} = a{1};  A2{end+1} = a{2};  Cc{end+1} = c; %#ok<AGROW>
    end
    if isempty(quad), quad = [1 0]; end                  % no branch row: t never occurs
    if numel(quad) == 2, quad = [1 0]; end               % a degenerate pair: every slope is rational
    if numel(A1) == 2
        % On two rows a certificate means ANTI-PARALLEL normals. They are parallel when the cross
        % product vanishes and opposed when the dot product is negative, and then lambda is
        % (|a2|^2, -a1.a2), both positive -- so nothing has to be divided to find it.
        cr = cross2(A1{1}, A2{1}, A1{2}, A2{2}, quad);
        if ratQ.qsSign(cr, quad, iv) ~= 0, return, end
        dt = ratQ.qsAdd(ratQ.qsMul(A1{1}, A1{2}, quad), ratQ.qsMul(A2{1}, A2{2}, quad));
        if ratQ.qsSign(dt, quad, iv) >= 0, return, end
        n2 = ratQ.qsAdd(ratQ.qsMul(A1{2}, A1{2}, quad), ratQ.qsMul(A2{2}, A2{2}, quad));
        lam = {n2, qneg(dt)};
    else
        lam = {cross2(A1{2}, A2{2}, A1{3}, A2{3}, quad), ...
               cross2(A1{3}, A2{3}, A1{1}, A2{1}, quad), ...
               cross2(A1{1}, A2{1}, A1{2}, A2{2}, quad)};
        sg = arrayfun(@(k) ratQ.qsSign(lam{k}, quad, iv), 1:3);
        if all(sg <= 0)
            for k = 1:3, lam{k} = qneg(lam{k}); end
            sg = -sg;
        end
        if ~all(sg >= 0), return, end                    % no nonnegative combination annihilates a
        if all(sg == 0), return, end                     % the three normals are all parallel
    end
    tot = [0 0 1];
    for k = 1:numel(Cc), tot = ratQ.qsAdd(tot, ratQ.qsMul(lam{k}, Cc{k}, quad)); end
    empty = ratQ.qsSign(tot, quad, iv) <= 0;
end

% --------------------------------------------------------------------------------------------------
function [a, c, ok] = exactHalfPlane(C, rows, t, ref)
% Row t as `a.x + c >= 0` with every coefficient an element of Q(t_ref), exactly, where t_ref is the
% root named by branch `ref` of the subsystem's one conjugate pair.
    a = {[0 0 1], [0 0 1]};  c = [0 0 1];  ok = false;
    con = C(rows(t,1), 1:6);  br = C(rows(t,1), 7);  sd = rows(t,2);
    if br == 0
        if ~all(con(1:3) == 0), return, end                 % a conic is not a half-plane
        a = {[sd*con(4) 0 1], [sd*con(5) 0 1]};
        c = [sd*con(6) 0 1];
        ok = true;
        return
    end
    pd = ratQ.pairLineData(con);
    if strcmp(pd.kind, 'crossing') && pd.vertical && br == 2
        p = pd.apex;                                        % the VERTICAL line x1 = p1, rational
        a = {[sd*p(3) 0 1], [0 0 1]};
        c = [-sd*p(1) 0 1];
        ok = true;
        return
    end
    T = rootElement(pd.quad, br, ref);                      % the slope or level, in Q(t_ref)
    if isempty(T), return, end
    if strcmp(pd.kind, 'parallel')
        % { form.x = T }, so the half-plane is sd*(form.x - T) >= 0
        a = {[sd*pd.form(1) 0 1], [sd*pd.form(2) 0 1]};
        c = ratQ.qsMul([-sd 0 1], T, pd.quad);
    else
        % (x2 - p2) - T (x1 - p1) = 0, cleared by the apex's denominator and oriented as
        % `ratQ.pairLineSign` writes it.
        p = pd.apex;
        a = {ratQ.qsMul([-sd*p(3) 0 1], T, pd.quad), [sd*p(3) 0 1]};
        c = ratQ.qsAdd([-sd*p(2) 0 1], ratQ.qsMul([sd*p(1) 0 1], T, pd.quad));
    end
    ok = true;
end

% --------------------------------------------------------------------------------------------------
function T = rootElement(quad, br, ref)
% The br-th root of `quad`, written in the field of its ref-th root.
%
% The two roots of al t^2 + be t + ga sum to -be/al, so the one that is not the generator is
% `(-be - al t)/al` -- rational in the other, which is why one quadratic field carries both lines of
% a pair and no compositum is needed for a subsystem that stays inside one pair.
    T = [];
    if numel(quad) == 2
        T = [-quad(2), 0, quad(1)];                         % degenerate: the slope is rational
        return
    end
    if numel(quad) ~= 3, return, end
    if br == ref
        T = [0 1 1];
    else
        T = [-quad(2), -quad(1), quad(1)];
    end
    T = ratQ.qsRed(T);
end

% --------------------------------------------------------------------------------------------------
function w = cross2(u1, u2, v1, v2, Q)
    w = ratQ.qsAdd(ratQ.qsMul(u1, v2, Q), qneg(ratQ.qsMul(u2, v1, Q)));
end

function u = qneg(u)
    u = [-u(1), -u(2), u(3)];
end

% --------------------------------------------------------------------------------------------------
function L = branchLine(c, b)
% One line of a degenerate conic in DOUBLES, oriented exactly as `ratQ.pairLineEval` evaluates it --
% so a side taken from one agrees with the other. Read off that routine itself rather than case-split
% on crossing and parallel pairs: the value is affine in x, so three points determine it.
    v = ratQ.pairLineEval(c, b, [0 0; 1 0; 0 1]);
    L = [v(2) - v(1), v(3) - v(1), v(1)];
    n = hypot(L(1), L(2));
    if n > 0, L = L / n; end
end

% --------------------------------------------------------------------------------------------------
function m = marginLP(H)
% max t subject to a_i.x + c_i >= t for unit a_i, inside a box. Positive iff the cell has interior.
% The optimum of a 3-variable LP sits where three constraints are tight, so they are enumerated.
    R = 1e4;                                            % far outside any cell of this pipeline
    for i = 1:size(H,1)
        n = hypot(H(i,1), H(i,2));
        if n > 0, H(i,:) = H(i,:) / n; end
    end
    M = [H(:,1:2), -ones(size(H,1),1); 1 0 0; -1 0 0; 0 1 0; 0 -1 0];
    b = [-H(:,3); -R; -R; -R; -R];
    m = -inf;  k = size(M,1);
    for i = 1:k-2
        for j = i+1:k-1
            for l = j+1:k
                S = M([i j l], :);
                if abs(det(S)) < 1e-12, continue, end
                z = S \ b([i j l]);
                if all(M*z >= b - 1e-9) && z(3) > m, m = z(3); end
            end
        end
    end
    if isinf(m), m = -1; end                            % no vertex satisfies all: empty
end
