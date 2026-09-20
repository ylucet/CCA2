function meshWithCorners()
% Does the mesh -- now carrying the corner faces' REGIONS for coverage -- agree with `envelopeEval`,
% VALUE and DOMAIN, on the fixture whose ties are curved? `envelopeEval` is the oracle: it is what
% was verified against the definition, and the mesh is only a re-arrangement of the same faces.
    L = load(fullfile(tempdir, 'branchRowPruning_env.mat'));  env = L.env;
    t0 = tic;
    [g, info] = envelopeMesh(env);
    fprintf('mesh: %d curves, %d cells, %d refined, %d faces, %d uncovered, %d corners | %.0f s\n', ...
            info.curves, info.cells, info.refined, info.faces, info.uncovered, info.corners, toc(t0));

    rng(19);
    X = [rand(400,1)*5-2.5, rand(400,1)*5-2.5];
    ref = envelopeEval(env, X);
    got = g.eval(X);
    fin = isfinite(ref);
    fprintf('domain: reference finite %d of %d | mesh finite %d | disagreements %d\n', ...
            sum(fin), numel(fin), sum(isfinite(got)), sum(fin ~= isfinite(got)));
    both = fin & isfinite(got);
    if any(both)
        fprintf('value: max |mesh - reference| over %d common points = %.3e\n', ...
                sum(both), max(abs(got(both) - ref(both))));
    end
end
