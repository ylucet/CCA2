function meshAlgProbe()
% Does the algebraic-edged mesh agree with `envelopeEval` on the CURVED-TIE fixture -- the one no
% mesh has ever been built for? `envelopeEval` is the oracle: it was verified against the definition
% and the mesh is a re-arrangement of the same faces.
    L = load(fullfile(tempdir, 'branchRowPruning_env.mat'));  env = L.env;
    t0 = tic;
    [mesh, info] = envelopeMeshAlg(env);
    fprintf('c = -2 fixture: %d cells -> %d faces, %d algebraic rows, %d cells left out | %.0f s\n', ...
            info.cells, info.faces, info.algRows, info.refused, toc(t0));
    rng(29);
    % The envelope of this fixture is finite on the UNIT SQUARE and +infinity outside it, so a wide
    % box spends its points proving the domain and almost none on the values. Both are wanted, so
    % both are sampled: the square for the VALUES, a ring around it for the DOMAIN.
    X = [rand(400,1), rand(400,1); rand(200,1)*4-1.5, rand(200,1)*4-1.5];
    ref = envelopeEval(env, X);
    got = meshAlgEval(mesh, X);
    dom = sum(isfinite(ref) ~= isfinite(got));
    both = isfinite(ref) & isfinite(got);
    fprintf('domain mismatches %d of %d | max |mesh - reference| %.3e over %d common points\n', ...
            dom, numel(ref), max(abs(got(both) - ref(both))), sum(both));
end
