function algCrvProbe()
% Does the TYPED object evaluate the same as the subdivision it was assembled from, and as the evaluator
% that was verified against the definition? Three answers, one function, on the curved-tie fixture.
    L = load(fullfile(tempdir, 'branchRowPruning_env.mat'));  env = L.env;
    subdivision = envelopeSubdivisionAlg(env);
    g = assembleAlgCrv(subdivision);
    fprintf('AlgCrv: %d faces, %d conic edges, %d algebraic edges | isa Crv %d, isa AlgCon %d\n', ...
            g.nf, size(g.EcQ,1), numel(g.EcA), isa(g,'Crv'), isa(g,'AlgCon'));
    rng(29);
    X = [rand(400,1), rand(400,1); rand(200,1)*4-1.5, rand(200,1)*4-1.5];
    ref = envelopeEval(env, X);
    m   = subdivisionAlgEval(subdivision, X);
    t   = g.eval(X);
    both = isfinite(ref) & isfinite(t);
    fprintf(['domain: subdivision vs reference %d mismatches, TYPE vs reference %d | ' ...
             'value: type vs reference %.3e over %d points\n'], ...
            sum(isfinite(m) ~= isfinite(ref)), sum(isfinite(t) ~= isfinite(ref)), ...
            max(abs(t(both) - ref(both))), sum(both));
end
