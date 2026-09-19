function branchRowPruning()
% How much does EXACT emptiness pruning buy over RELAXED, once the corner faces' BRANCH rows are in
% the arrangement? The question is not the cost of one test -- exact is obviously the dearer call --
% but the cost of the cells a weak test keeps, because every later curve splits them again.
%
% The measure is the LIVE CELL COUNT after each curve, which is machine-independent (AI/CLAUDE.md
% section 3); wall-clock is reported as the MINIMUM of three repetitions, on a shared machine.
%
% RELAXED  : a branch row is IGNORED in the emptiness test (sound -- it can only keep cells).
% CERTIFIED: the branch rows and the straight rational rows are additionally tested TOGETHER as
%            half-planes, by maximising the inscribed margin; the cell is dropped only when that
%            margin is negative by more than the coefficient error can explain. Undecided cells are
%            KEPT, so this is sound too, and the undecided count is reported rather than hidden.
    cacheFile = fullfile(tempdir, 'branchRowPruning_env.mat');
    if exist(cacheFile, 'file')
        L = load(cacheFile);  env = L.env;                  % the assembly is ~7 minutes; reruns are
    else                                                    % about the arrangement, not about it
        fA = tri([0 0; 1 0; 1 1], 0);  fB = tri([0 0; 1 1; 0 1], -2);
        env = envelopePrune(envelopeAssemble(conjQ(fA), conjQ(fB)));
        save(cacheFile, 'env');
    end

    [C, nQ, nB] = curveTable(env);
    fprintf('curves: %d total (%d from quadratic faces, %d branch rows from corners)\n', ...
            size(C,1), nQ, nB);

    % "aggressive" is NOT a candidate implementation: a zero margin means no INTERIOR, so treating
    % it as empty is what an EXACT test concludes -- but in doubles it would also drop a genuinely
    % thin cell. It bounded what exactness could buy, and `cellInterior` is the implementation that
    % reaches that bound soundly, deciding each tie by an exact Farkas certificate instead.
    for mode = ["relaxed", "certified", "aggressive", "exact"]
        best = inf;  rec = [];  stats = [];
        for rep = 1:3
            t0 = tic;
            [cells, st] = arrange(C, mode);
            el = toc(t0);
            if el < best, best = el;  rec = cells;  stats = st; end
        end
        fprintf(['%-9s | live cells %5d | emptiness calls %6d | undecided %4d | ' ...
                 'per-curve %s | min of 3 runs %.1f s\n'], ...
                mode, numel(rec), stats.calls, stats.undecided, mat2str(stats.perCurve), best);
    end
end

% ==================================================================================================
function [cells, st] = arrange(C, mode)
% The incremental arrangement, with the emptiness test as the only difference between the two runs.
    cells = {zeros(0,2)};
    st = struct('calls', 0, 'undecided', 0, 'perCurve', []);
    for j = 1:size(C,1)
        next = {};
        for k = 1:numel(cells)
            for sd = [1 -1]
                cand = [cells{k}; j, sd];
                st.calls = st.calls + 1;
                if strcmp(mode, "exact")
                    [live, info] = cellInterior(C, cand);
                    und = double(strcmp(info.path, 'unsureTie'));
                else
                    [live, und] = isLive(C, cand, mode);
                end
                st.undecided = st.undecided + und;
                if live, next{end+1} = cand; end %#ok<AGROW>
            end
        end
        cells = next;
        st.perCurve(end+1) = numel(cells); %#ok<AGROW>
        if numel(cells) > 20000
            fprintf('  (%s: over 20000 live cells after %d of %d curves -- stopping)\n', ...
                    mode, j, size(C,1));
            return
        end
    end
end

% --------------------------------------------------------------------------------------------------
function [tf, und] = isLive(C, rows, mode)
% The SAME exact test on the rational rows in both modes; the branch rows are what differ.
    und = 0;
    lin = zeros(0,3);  curved = zeros(0,6);  sides = [];  half = zeros(0,3);
    for t = 1:size(rows,1)
        c = C(rows(t,1), 1:6);  br = C(rows(t,1), 7);  sd = rows(t,2);
        if br ~= 0
            if ~strcmp(mode, "relaxed")
                L = branchLine(c, br);                 % the named line, in doubles
                half(end+1,:) = sd * L; %#ok<AGROW>
            end
            continue
        end
        if all(c(1:3) == 0)
            lin(end+1,:) = sd * c(4:6); %#ok<AGROW>
            half(end+1,:) = sd * c(4:6); %#ok<AGROW>
        else
            curved(end+1,:) = c; %#ok<AGROW>
            sides(end+1) = sd; %#ok<AGROW>
        end
    end
    if ~isempty(lin) && ~ratQ.feasible2(lin, true), tf = false; return, end
    if ~isempty(curved)
        try
            tf = cellHasInterior([curved; zeros(size(lin,1),3), lin], ...
                                 [(1:size(curved,1)).', sides(:); ...
                                  (size(curved,1)+(1:size(lin,1))).', ones(size(lin,1),1)]);
        catch
            % The exact kernel can overflow on this fixture's curves (`ratQ.chk` on a line base
            % point). Its own contract is to answer TRUE when unsure, so that is what happens here,
            % identically in both modes so the comparison stays fair.
            tf = true;
        end
        if ~tf, return, end
    end
    tf = true;
    if ~strcmp(mode, "relaxed") && size(half,1) >= 2
        [m, scale] = margin(half);
        if m < -1e-9 * scale
            tf = false;                                % empty, by more than the error can explain
        elseif abs(m) <= 1e-9 * scale
            und = 1;                                   % too close to call
            if strcmp(mode, "aggressive"), tf = false; end   % what an EXACT test would decide here
        end
    end
end

% --------------------------------------------------------------------------------------------------
function L = branchLine(c, b)
% One line of a degenerate conic, as [u v w] with u x + v y + w, in doubles -- the same line
% `ratQ.pairLineEval` evaluates, and oriented the same way, so a side read from one matches the other.
    persistent cache
    if isempty(cache), cache = containers.Map('KeyType', 'char', 'ValueType', 'any'); end
    key = sprintf('%d,', c, b);
    if cache.isKey(key), L = cache(key); return, end
    % Two points of the line and its orientation, read off `pairLineEval` itself: fit u x + v y + w
    % to the values it returns, which is exact for a line and needs no case split on crossing or
    % parallel pairs.
    X = [0 0; 1 0; 0 1];
    v = ratQ.pairLineEval(c, b, X);
    L = [v(2) - v(1), v(3) - v(1), v(1)];
    n = hypot(L(1), L(2));
    if n > 0, L = L / n; end
    cache(key) = L;
end

% --------------------------------------------------------------------------------------------------
function [m, scale] = margin(H)
% The largest inscribed margin of { a_i.x + c_i >= 0 } inside a box: max t with a_i.x + c_i >= t for
% unit a_i. Positive means the cell has interior, negative means it is empty. In two dimensions the
% optimum sits where three of the constraints are tight, so the triples are enumerated -- there are
% at most a few hundred of them and no solver is needed.
    R = 50;                                             % the box every cell of this fixture lies in
    A = [H(:,1:2), ones(size(H,1),1)];                  % a.x + c >= t  <=>  [a, -1].[x; t] >= -c
    for i = 1:size(A,1)
        n = hypot(H(i,1), H(i,2));
        if n > 0, A(i,:) = [H(i,1:2)/n, 1];  H(i,3) = H(i,3)/n; end
    end
    A = [A; 1 0 0; -1 0 0; 0 1 0; 0 -1 0];
    b = [-H(:,3); -R; -R; -R; -R];                      % A(:,1:2) x - A(:,3) t >= b
    M = [A(:,1:2), -A(:,3)];
    scale = max(1, max(abs(b)));
    m = -inf;  k = size(M,1);
    for i = 1:k-2
        for j = i+1:k-1
            for l = j+1:k
                S = M([i j l], :);
                if abs(det(S)) < 1e-12, continue, end
                z = S \ b([i j l]);
                if all(M*z >= b - 1e-9*scale) && z(3) > m, m = z(3); end
            end
        end
    end
    if isinf(m), m = -scale; end                        % no vertex satisfies all: empty
end

% --------------------------------------------------------------------------------------------------
function [C, nQ, nB] = curveTable(env)
% Every distinct (conic, branch) the faces' regions are built from -- the quadratic faces' rational
% rows, and the corner faces' branch rows.
    C = zeros(0,7);  nQ = 0;  nB = 0;
    for i = 1:numel(env.quad)
        con = env.quad(i).con;
        for t = 1:size(con,1)
            [C, isNew] = add(C, [con(t,1:6), 0]);
            nQ = nQ + isNew;
        end
    end
    for k = 1:numel(env.bridge)
        f = env.bridge{k};
        if ~strcmp(f.kind, 'corner') || isempty(f.con), continue, end
        for t = 1:size(f.con,1)
            [C, isNew] = add(C, [f.con(t,1:6), f.con(t,8)]);
            nB = nB + isNew;
        end
    end
end

function [C, isNew] = add(C, r)
% A row with no x or y in it is not a CURVE: it is a constant, either vacuously true or a face that
% is empty, and it must not enter the arrangement. Letting one through is what made the exact kernel
% throw `ratQ:notFinite` -- it was asked to meet the "curve" 1 = 0 -- and a throw there is caught as
% "unsure", which silently turned the pruning off from curve 10 onwards.
    isNew = 0;
    if all(r(1:5) == 0), return, end
    if isempty(find(all(C == r, 2), 1)), C(end+1,:) = r;  isNew = 1; end
end

function f = tri(V, c)
    q = [0 0 0 0, 2, 0, 2, 0, 0, 0];  q(8) = q(8) - c;  q(9) = q(9) + c;
    f = QuaPol(V, [1 2 1; 2 3 1; 3 1 1], q, [1 0; 1 0; 1 0]);
end
