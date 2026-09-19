function branchRowPruning2()
% WHY does the arrangement stop pruning at curve 11? Counted, not guessed: how many emptiness calls
% reach the exact conic kernel, how many it answers, and how many it THROWS on. Plus a standalone
% check of the margin LP the certified mode rests on, because a test that is always undecided
% measures nothing.
    L = load(fullfile(tempdir, 'branchRowPruning_env.mat'));  env = L.env;
    C = curveTable(env);

    % ---- the margin LP, against cases whose answer is known --------------------------------------
    fprintf('margin: wedge %+0.3f (want > 0), opposite %+0.3f (want < 0), slab %+0.3f (want > 0)\n', ...
            margin([1 0 0; 0 1 0]), margin([1 0 -1; -1 0 -1]), margin([1 0 1; -1 0 1]));

    % ---- where the pruning goes ------------------------------------------------------------------
    cells = {zeros(0,2)};
    tot = struct('calls',0,'byLin',0,'byConic',0,'threw',0,'kept',0);
    shown = 0;
    for j = 1:min(14, size(C,1))
        next = {};
        for k = 1:numel(cells)
            for sd = [1 -1]
                cand = [cells{k}; j, sd];
                tot.calls = tot.calls + 1;
                [lin, curved, sides] = split(C, cand);
                if ~isempty(lin) && ~ratQ.feasible2(lin, true)
                    tot.byLin = tot.byLin + 1;  continue
                end
                live = true;
                if ~isempty(curved)
                    try
                        live = cellHasInterior([curved; zeros(size(lin,1),3), lin], ...
                                 [(1:size(curved,1)).', sides(:); ...
                                  (size(curved,1)+(1:size(lin,1))).', ones(size(lin,1),1)]);
                        if ~live, tot.byConic = tot.byConic + 1; end
                    catch err
                        tot.threw = tot.threw + 1;
                        if shown < 2
                            shown = shown + 1;
                            fprintf('threw (%s) with %d conics, max |coef| %.3g\n', ...
                                    err.identifier, size(curved,1), max(abs(curved(:))));
                        end
                        live = true;
                    end
                end
                if live, next{end+1} = cand;  tot.kept = tot.kept + 1; end %#ok<AGROW>
            end
        end
        cells = next;
        fprintf('curve %2d: %5d live | calls %6d linDrop %5d conicDrop %5d threw %5d\n', ...
                j, numel(cells), tot.calls, tot.byLin, tot.byConic, tot.threw);
    end
end

function [lin, curved, sides] = split(C, rows)
    lin = zeros(0,3);  curved = zeros(0,6);  sides = [];
    for t = 1:size(rows,1)
        c = C(rows(t,1), 1:6);  br = C(rows(t,1), 7);  sd = rows(t,2);
        if br ~= 0, continue, end
        if all(c(1:3) == 0), lin(end+1,:) = sd * c(4:6); %#ok<AGROW>
        else, curved(end+1,:) = c;  sides(end+1) = sd; %#ok<AGROW>
        end
    end
end

function m = margin(H)
    R = 50;
    for i = 1:size(H,1)
        n = hypot(H(i,1), H(i,2));
        if n > 0, H(i,:) = H(i,:) / n; end
    end
    M = [H(:,1:2), -ones(size(H,1),1)];
    b = -H(:,3);
    M = [M; 1 0 0; -1 0 0; 0 1 0; 0 -1 0];
    b = [b; -R; -R; -R; -R];
    m = -inf;  k = size(M,1);
    for i = 1:k-2, for j = i+1:k-1, for l = j+1:k
        S = M([i j l], :);
        if abs(det(S)) < 1e-12, continue, end
        z = S \ b([i j l]);
        if all(M*z >= b - 1e-9) && z(3) > m, m = z(3); end
    end, end, end
    if isinf(m), m = -R; end
end

function C = curveTable(env)
    C = zeros(0,7);
    for i = 1:numel(env.quad)
        con = env.quad(i).con;
        for t = 1:size(con,1), C = add(C, [con(t,1:6), 0]); end
    end
    for k = 1:numel(env.bridge)
        f = env.bridge{k};
        if ~strcmp(f.kind, 'corner') || isempty(f.con), continue, end
        for t = 1:size(f.con,1), C = add(C, [f.con(t,1:6), f.con(t,8)]); end
    end
end

function C = add(C, r)
    if isempty(find(all(C == r, 2), 1)), C(end+1,:) = r; end
end
