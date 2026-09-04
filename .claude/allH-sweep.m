function allH_sweep(nCase, nPts)
% Sweep conjQ over random polygons and EVERY Hessian class -- positive definite, indefinite,
% negative definite, PSD-singular, NSD-singular and zero -- against a general oracle.
%
% THE ORACLE enumerates the candidate maximisers of a quadratic over a polytope: every vertex,
% every edge's clamped 1-D stationary point (whatever the sign of that edge's curvature), and the
% unconstrained stationary point when it lies inside. For a concave objective one of those IS the
% maximiser; for an indefinite one the maximiser is on the boundary, hence a vertex or an edge
% stationary point, so the list is still complete. Every candidate is a genuine point of P, so a
% spurious one can only ever be dominated.
    if nargin < 1, nCase = 60; end
    if nargin < 2, nPts  = 150; end
    rng(20260904);
    tot = 0; bad = 0; unc = 0; worst = 0;
    byClass = containers.Map();
    for c = 1:nCase
        [V, E, F] = randomPolygon();
        [H, cls] = randomH(mod(c, 6));
        L = randi([-4 4], 2, 1);  k0 = randi([-4 4]);
        f = QuaPol(V, E, [0 0 0 0, H(1,1), H(1,2), H(2,2), L(1), L(2), k0], F);
        if ~f.isExact(), continue, end
        try
            g = conjQ(f);
        catch e
            fprintf('case %2d [%s] ERROR %s\n', c, cls, e.identifier);
            continue
        end
        S = [randn(nPts,2)*3; 0 0];
        [got, idx] = g.eval(S);
        unc = unc + sum(idx == 0);
        want = zeros(size(S,1),1);
        for i = 1:size(S,1)
            want(i) = supQ(S(i,:).', V, H, L, k0);
        end
        d = abs(got-want)./max(1,abs(want));  d(~isfinite(got)) = inf;
        n = sum(d > 1e-9);
        bad = bad + n;  tot = tot + numel(d);  worst = max(worst, max(d));
        if ~isKey(byClass, cls), byClass(cls) = [0 0 0]; end
        v = byClass(cls);  byClass(cls) = [v(1)+1, v(2)+n, v(3)+g.nf];
        if n > 0
            fprintf('case %2d [%s] %d wrong, nf=%d\n', c, cls, n, g.nf);
        end
    end
    ks = keys(byClass);
    for i = 1:numel(ks)
        v = byClass(ks{i});
        fprintf('  %-16s cases %2d  wrong %d  mean nf %.1f\n', ks{i}, v(1), v(2), v(3)/v(1));
    end
    fprintf('ALLH points %d  wrong %d  uncovered %d  worst %.3e\n', tot, bad, unc, worst);
end

function v = supQ(s, V, H, L, k0)
    q   = @(x) 0.5 * x.' * H * x + L.' * x + k0;
    obj = @(x) s.' * x - q(x);
    cand = V.';
    if abs(det(H)) > 1e-12
        xs = H \ (s - L);
        if inPoly(xs.', V), cand = [cand, xs]; end
    end
    m = size(V,1);
    for j = 1:m
        a = V(j,:).';  b = V(mod(j,m)+1,:).';  d = b - a;
        al = d.' * H * d;
        if abs(al) < 1e-14, continue, end          % affine along this edge: endpoints suffice
        t = min(1, max(0, (s - (H*a + L)).' * d / al));
        cand = [cand, a + t*d]; %#ok<AGROW>
    end
    v = -inf;
    for i = 1:size(cand,2), v = max(v, obj(cand(:,i))); end
end

function tf = inPoly(x, V)
    m = size(V,1);  ctr = mean(V,1);  tf = true;
    for j = 1:m
        a = V(j,:);  b = V(mod(j,m)+1,:);  d = b - a;  n = [d(2), -d(1)];
        if n * (ctr - a).' > 0, n = -n; end
        if n * (x - a).' > 1e-12, tf = false; return, end
    end
end

function [H, cls] = randomH(which)
    switch which
        case 0
            while true
                H = randi([-4 4],2,2); H = H+H.';
                if H(1,1)>0 && det(H)>0.5, break, end
            end
            cls = 'PD';
        case 1
            while true
                H = randi([-4 4],2,2); H = H+H.';
                if det(H) < -0.5, break, end
            end
            cls = 'indefinite';
        case 2
            while true
                H = randi([-4 4],2,2); H = H+H.';
                if H(1,1)<0 && det(H)>0.5, break, end
            end
            cls = 'ND';
        case 3
            u = randi([-3 3],2,1);  if all(u==0), u = [1;0]; end
            H = u*u.';  cls = 'PSD-singular';
        case 4
            u = randi([-3 3],2,1);  if all(u==0), u = [1;0]; end
            H = -u*u.';  cls = 'NSD-singular';
        otherwise
            H = zeros(2);  cls = 'affine';
    end
end

function [V, E, F] = randomPolygon()
    m  = randi([3 6]);
    th = sort(rand(1,m)*2*pi);
    r  = randi([1 6],1,m);
    V  = round([r'.*cos(th'), r'.*sin(th')]*2)/2;
    V  = uniquetol(V, 1e-9, 'ByRows', true);
    k  = convhull(V(:,1), V(:,2));
    V  = V(k(1:end-1), :);
    m  = size(V,1);
    E  = [(1:m).', [2:m,1].', ones(m,1)];
    F  = [ones(m,1), zeros(m,1)];
end
