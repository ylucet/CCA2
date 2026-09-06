function recessionRateProbe()
% RUNG 1: OBSERVE. biconjQ's divergence test says it checks whether q's linear rate along a
% recession direction is "negative SOMEWHERE ON THE PIECE", but the code evaluates `L'*d` -- the
% rate at the ORIGIN, not on the piece. Those differ whenever H*d is nonzero, because the true
% rate along d at a point x is (H*x + L)'*d, which then varies with x.
%
% If they differ on a reachable input, biconjQ is giving the WRONG DIAGNOSIS today (it says
% "finite, not computed yet" where the answer is -infinity). Both are refusals, so nothing is
% wrong-answered right now -- but the unbounded envelope is exactly the case that turns this
% branch into a computation, so it has to be right before that lands.
%
% THE FIXTURE. q = x*y, whose H = [0 1; 1 0] is indefinite, on the translated quadrant
% {x >= 0, y >= -5} = (0,-5) + cone((1,0), (0,1)).
%   * both recession directions have d'H d = 0, so neither is refused for curvature;
%   * L = 0, so the code's L'*d is 0 for both -- not negative, so it concludes "not diverging";
%   * the cone test has a = c = 0 and b = (1,0)*H*(0,1)' = 1 >= 0, so that does not fire either;
%   * but the TRUE rate along d = (1,0) at the vertex (0,-5) is (H*x)'d = y = -5 < 0.
%
% So q runs downward along that ray and co f is -infinity. No envelope machinery is needed to see
% it: f itself is evaluated below, marching out along the ray.
    V = [0 -5; 1 -5; 0 -4];
    E = [1 2 0; 1 3 0];
    F = [1 0; 0 1];
    q = [0 0 0 0, 0, 1, 0, 0, 0, 0];      % weighted basis: H = [0 1; 1 0], i.e. q = x*y
    f = QuaPol(V, E, q, F);

    fprintf('H = [0 1; 1 0] (q = x*y) on {x >= 0, y >= -5}\n\n');
    fprintf('f marched out along the ray (t, -5):\n');
    for t = [0 1 10 100 1000 10000]
        fprintf('   f(%6g, -5) = %12g\n', t, f.eval([t -5]));
    end
    fprintf('\n-> f is unbounded BELOW on its own domain, so co f = -infinity.\n\n');

    try
        h = biconjQ(f);
        fprintf('biconjQ ANSWERED with nf=%d -- that would be a wrong answer.\n', h.nf);
    catch e
        fprintf('biconjQ raises: %s\n', e.identifier);
        if strcmp(e.identifier, 'PLQ:biconjQ:minusInfinity')
            fprintf('   correct diagnosis.\n');
        else
            fprintf(['   WRONG DIAGNOSIS: it should be minusInfinity. This branch claims the\n' ...
                     '   envelope is FINITE and merely uncomputed, which is what the unbounded\n' ...
                     '   envelope item would then go on to compute.\n']);
        end
    end

    fprintf('\n--- the control: the SAME cone at the origin, where L''d and the true rate agree ---\n');
    f0 = QuaPol([0 0; 1 0; 0 1], [1 2 0; 1 3 0], q, [1 0; 0 1]);
    fprintf('f0 along (t, 0): ');
    fprintf('%g ', f0.eval([(0:3).' zeros(4,1)]));
    fprintf('\n-> x*y >= 0 on the first quadrant, so THIS one really is bounded below.\n');
    try
        biconjQ(f0);
        fprintf('biconjQ answered.\n');
    catch e
        fprintf('biconjQ raises: %s\n', e.identifier);
    end
end
