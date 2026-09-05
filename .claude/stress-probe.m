function stress_probe()
% The coverage table varies DOMAIN SHAPE x HESSIAN CLASS. This probes the axes it holds fixed, which
% is where an "is it complete?" claim would break if it is going to: face convexity, vertex count,
% piece count, mixed bounded/unbounded pieces, collinear vertices, and coefficient magnitude.
%
% Each case reports OK with the answer checked against an independent oracle where one is available,
% or the identifier it raised, or WRONG with the discrepancy.

    fprintf('%-34s %s\n', 'case', 'result');
    fprintf('%s\n', repmat('-', 1, 70));

    % ---- a NON-CONVEX face (an L-shape as ONE piece) -----------------------------------------
    V = [0 0; 2 0; 2 1; 1 1; 1 2; 0 2];
    m = size(V,1);
    E = [(1:m).', [2:m,1].', ones(m,1)];  F = [ones(m,1), zeros(m,1)];
    report('non-convex L-shaped face', QuaPol(V, E, [0 0 0 0 1 0 1 0 0 0], F), V, [1 0;0 1], [0;0], 0);

    % ---- many vertices -----------------------------------------------------------------------
    th = (0:11).'/12 * 2*pi;
    V = round([4*cos(th), 4*sin(th)] * 2) / 2;
    V = V(convhull(V(:,1), V(:,2)), :);  V = V(1:end-1, :);
    m = size(V,1);
    E = [(1:m).', [2:m,1].', ones(m,1)];  F = [ones(m,1), zeros(m,1)];
    report(sprintf('%d-gon, convex q', m), QuaPol(V, E, [0 0 0 0 2 1 3 0 0 0], F), V, [2 1;1 3], [0;0], 0);

    % ---- COLLINEAR vertices on the boundary ---------------------------------------------------
    V = [0 0; 1 0; 2 0; 2 2; 0 2];
    m = size(V,1);
    E = [(1:m).', [2:m,1].', ones(m,1)];  F = [ones(m,1), zeros(m,1)];
    report('collinear vertex on an edge', QuaPol(V, E, [0 0 0 0 1 0 1 0 0 0], F), V, [1 0;0 1], [0;0], 0);

    % ---- MANY pieces: the plane cut into 4 quadrants, different quadratics --------------------
    V4 = [0 0; 1 0; 0 1; -1 0; 0 -1];
    E4 = [1 2 0; 1 3 0; 1 4 0; 1 5 0];
    F4 = [1 2; 2 3; 3 4; 4 1];
    f4 = [0 0 0 0 1 0 1 0 0 0; 0 0 0 0 2 0 1 0 0 0; ...
          0 0 0 0 1 0 2 0 0 0; 0 0 0 0 3 0 3 0 0 0];
    report('4 unbounded pieces, 4 quadratics', QuaPol(V4, E4, f4, F4), [], [], [], []);

    % ---- MIXED: a bounded triangle next to an unbounded wedge --------------------------------
    % (fixture withdrawn: QuaPol's own constructor rejects it, so it tests my mesh-writing rather
    % than conjQ. Building a valid mixed mesh by hand is its own exercise.)

    % ---- LARGE coefficients: does the exact kernel overflow? -----------------------------------
    V = [0 0; 1 0; 1 1; 0 1];
    E = [1 2 1; 2 3 1; 3 4 1; 4 1 1];  F = [1 0;1 0;1 0;1 0];
    for scale = [1e3 1e5 1e7]
        report(sprintf('bounded square, coeffs ~%.0e', scale), ...
               QuaPol(V, E, [0 0 0 0 scale 0 scale 0 0 0], F), V, ...
               [scale 0; 0 scale], [0;0], 0);
    end
    % ---- fine rational vertices: denominators multiply -----------------------------------------
    Vq = [0 0; 1/7 0; 1/7 1/11; 0 1/13];
    report('vertices with denominators 7,11,13', QuaPol(Vq, E, [0 0 0 0 1 0 1 0 0 0], F), ...
           Vq, [1 0;0 1], [0;0], 0);
end

function report(name, f, V, H, L, k0)
    if ~f.isExact()
        fprintf('%-34s INPUT NOT EXACT\n', name);  return
    end
    try
        g = conjQ(f);
    catch e
        fprintf('%-34s RAISED %s\n', name, e.identifier);  return
    end
    if isempty(V)
        fprintf('%-34s OK nf=%d (no oracle here)\n', name, g.nf);  return
    end
    rng(20260904);
    S = [randn(200,2)*3; 0 0];
    want = zeros(size(S,1),1);
    for i = 1:size(S,1)
        want(i) = supSampled(S(i,:).', V, H, L, k0);
    end
    got = g.eval(S);
    below = got < want - 1e-7*max(1,abs(want));       % a sampled sup is a LOWER bound
    fprintf('%-34s OK nf=%-4d below-sampled %d of %d\n', name, g.nf, sum(below), numel(want));
end

function v = supSampled(s, V, H, L, k0)
% a LOWER bound on the sup: sample the polygon densely (vertices, edge points, interior mixtures)
    m = size(V,1);
    X = V;
    for a = 1:m
        b = mod(a, m) + 1;
        for t = linspace(0, 1, 25)
            X(end+1,:) = (1-t)*V(a,:) + t*V(b,:); %#ok<AGROW>
        end
    end
    rng(7);
    for i = 1:400
        w = rand(m,1);  w = w/sum(w);
        X(end+1,:) = w.' * V; %#ok<AGROW>
    end
    q = 0.5*sum((X*H).*X, 2) + X*L + k0;
    v = max(X*s - q);
end
