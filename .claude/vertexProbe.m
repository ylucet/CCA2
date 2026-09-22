function vertexProbe()
% Do the named vertices exist, lie on BOTH their edges, and satisfy the face they belong to? The
% definition is the check: a vertex named [edgeA tableA edgeB tableB rootIdx] must be a point where
% both curves vanish and every row of its face holds.
    L = load(fullfile(tempdir, 'branchRowPruning_env.mat'));  env = L.env;
    g = assembleAlgCrv(envelopeMeshAlg(env));
    [V, info] = algCrvVertices(g);
    fprintf('%d faces | candidates %d | NAMED %d | refused %d | %.0f s\n', ...
            g.nf, info.candidates, info.named, info.refused, info.seconds);
    if isempty(V), fprintf('no vertices named\n'); return, end
    fprintf('distinct edge pairs used: %d | tables seen: %s\n', ...
            size(unique(V(:,1:4), 'rows'),1), mat2str(unique(V(:,[2 4]))'));
end
