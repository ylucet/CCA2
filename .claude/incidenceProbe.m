function incidenceProbe()
% Is the incidence CONSISTENT with the tables it was read from? Three checks, each against the
% definition rather than against a stored number.
    L = load(fullfile(tempdir, 'branchRowPruning_env.mat'));  env = L.env;
    g = assembleAlgCrv(envelopeMeshAlg(env));
    V = algCrvVertices(g);
    [E, F, info] = algCrvIncidence(g, V);
    fprintf(['%d edges | segments %d | rays %d | shared by two faces %d | ' ...
             'more than two vertices %d | %d vertices\n'], ...
            info.edges, info.segments, info.rays, info.shared, info.oversubdivided, size(V,1));

    % 1. every endpoint named on an edge must be a vertex whose NAME mentions that edge
    bad = 0;
    for e = 1:size(E,1)
        for c = [3 4]
            v = E(e,c);
            if v == 0, continue, end
            named = (V(v,1) == E(e,1) && V(v,2) == E(e,2)) || ...
                    (V(v,3) == E(e,1) && V(v,4) == E(e,2));
            if ~named, bad = bad + 1; end
        end
    end
    fprintf('endpoints that do not mention their edge: %d\n', bad);

    % 2. every face listed in F must really reference that edge in its own rows
    bad2 = 0;
    for e = 1:size(E,1)
        for c = 1:2
            k = F(e,c);
            if k == 0, continue, end
            rows = g.FC{k};
            hit = false;
            for t = 1:size(rows,1)
                rt = 0;
                if size(rows,2) >= 4, rt = rows(t,4); end
                if rows(t,1) == E(e,1) && rt == E(e,2), hit = true; break, end
            end
            if ~hit, bad2 = bad2 + 1; end
        end
    end
    fprintf('faces listed for an edge they do not use: %d\n', bad2);

    % 3. every edge of every face must appear in E
    bad3 = 0;
    for k = 1:g.nf
        rows = g.FC{k};
        for t = 1:size(rows,1)
            rt = 0;
            if size(rows,2) >= 4, rt = rows(t,4); end
            if isempty(find(E(:,1) == rows(t,1) & E(:,2) == rt, 1)), bad3 = bad3 + 1; end
        end
    end
    fprintf('face rows with no edge entry: %d\n', bad3);
end
