function compactProbe()
% The documented 321-face mesh, re-measured with implied rows dropped. Same fixture as
% envelopeAssembleTest's convex case: two triangles on the unit square, c = 0.
    fA = tri([0 0; 1 0; 1 1], 0);  fB = tri([0 0; 1 1; 0 1], 0);
    env = envelopePrune(envelopeAssemble(conjQ(fA), conjQ(fB)));
    t0 = tic;
    [g, info] = envelopeMesh(env);
    fprintf('convex fixture: %d curves, %d cells, %d refined, %d FACES, %d edges | %.0f s\n', ...
            info.curves, info.cells, info.refined, info.faces, g.ne, toc(t0));
    rng(13);
    X = [rand(300,1), rand(300,1)];
    want = min(fA.eval(X), fB.eval(X));
    got = g.eval(X);
    fin = isfinite(want);
    fprintf('against f itself: %d finite points, max |mesh - f| = %.3e, domain mismatches %d\n', ...
            sum(fin), max(abs(got(fin) - want(fin))), sum(isfinite(got) ~= fin));
end

function f = tri(V, c)
    q = [0 0 0 0, 2, 0, 2, 0, 0, 0];  q(8) = q(8) - c;  q(9) = q(9) + c;
    f = QuaPol(V, [1 2 1; 2 3 1; 3 1 1], q, [1 0; 1 0; 1 0]);
end
