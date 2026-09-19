function branchRowPruning3()
% The smallest input that makes the exact kernel throw. Shrink first, fix nothing yet.
    L = load(fullfile(tempdir, 'branchRowPruning_env.mat'));  env = L.env;
    C = curveTable(env);
    seen = containers.Map('KeyType','char','ValueType','logical');
    for i = 1:size(C,1)
        for j = i+1:size(C,1)
            if C(i,7) ~= 0 || C(j,7) ~= 0, continue, end
            a = C(i,1:6);  b = C(j,1:6);
            try
                conicMeet(a, b);
            catch err
                key = sprintf('%d-%d', i, j);
                if ~seen.isKey(err.identifier)
                    seen(err.identifier) = true;
                    fprintf('%s : conicMeet(%s, %s)\n', err.identifier, mat2str(a), mat2str(b));
                end
            end
        end
    end
    if seen.Count == 0, fprintf('no pair of curves throws on its own\n'); end
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
