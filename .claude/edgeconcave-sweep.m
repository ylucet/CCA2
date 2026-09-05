function edgeconc_sweep(nCase)
% Validate the EDGE-CONCAVE criterion that biconjQ now relies on: if q is concave along every edge
% DIRECTION of the polygon, its convex envelope is the lower hull of the lifted vertices -- even
% when H itself is indefinite.
%
% The criterion is a claim about mathematics, so it is checked against the DEFINITION of a convex
% envelope rather than against another implementation:
%   (1) co f <= f on the domain          -- the failure mode if the criterion is FALSE: the lower
%                                           hull would rise above f somewhere
%   (2) equality at the extreme points   -- the hull interpolates the vertex data
%   (3) co f is convex
%   (4) no affine minorant of f exceeds co f  -- it really is the LARGEST convex minorant
    if nargin < 1, nCase = 200; end
    rng(20260904);
    tried = 0; indef = 0; b1 = 0; b2 = 0; b3 = 0; b4 = 0; errs = 0;
    for c = 1:nCase
        [V, E, F] = randomPolygon();
        H = randi([-4 4], 2, 2);  H = H + H.';
        L = randi([-3 3], 2, 1);  k0 = randi([-3 3]);

        % keep only EDGE-CONCAVE cases, and skip the ones where H is PSD (then q is convex and
        % biconjQ returns f, which is a different branch)
        m = size(V,1);
        curv = zeros(m,1);
        for j = 1:m
            d = V(mod(j,m)+1,:) - V(j,:);
            curv(j) = d * H * d.';
        end
        if any(curv > 1e-12), continue, end
        isPSD = H(1,1) >= 0 && H(2,2) >= 0 && det(H) >= -1e-12;
        if isPSD, continue, end
        tried = tried + 1;
        if det(H) < -1e-12, indef = indef + 1; end

        f = QuaPol(V, E, [0 0 0 0, H(1,1), H(1,2), H(2,2), L(1), L(2), k0], F);
        if ~f.isExact(), continue, end
        try
            h = biconjQ(f);
        catch e
            errs = errs + 1;
            fprintf('case %d ERROR %s\n', c, e.identifier);
            continue
        end

        X = inHull(V, 400);
        qx = 0.5*sum((X*H).*X, 2) + X*L + k0;
        hx = h.eval(X);
        if any(hx > qx + 1e-9*max(1,abs(qx))), b1 = b1 + 1; end          % (1)

        qv = 0.5*sum((V*H).*V, 2) + V*L + k0;
        ext = extremePoints(V);
        if any(abs(h.eval(V(ext,:)) - qv(ext)) > 1e-9*max(1,abs(qv(ext)))), b2 = b2 + 1; end   % (2)

        A = inHull(V, 200);  B = inHull(V, 200);
        if any(h.eval((A+B)/2) > (h.eval(A) + h.eval(B))/2 + 1e-9), b3 = b3 + 1; end           % (3)

        bad4 = false;
        for t = 1:15
            g = randn(1,2);
            b = min(qv - V*g.');
            if any(X*g.' + b > hx + 1e-9*max(1,abs(hx))), bad4 = true; break, end
        end
        if bad4, b4 = b4 + 1; end                                                              % (4)
    end
    fprintf(['EDGE-CONCAVE cases %d (indefinite %d)  errors %d | ' ...
             'above-f %d  vertex %d  nonconvex %d  not-largest %d\n'], ...
            tried, indef, errs, b1, b2, b3, b4);
end

function ext = extremePoints(V)
    m = size(V,1);  ext = true(m,1);
    for i = 1:m
        a = V(mod(i-2,m)+1,:);  b = V(mod(i,m)+1,:);
        if abs(det([b-a; V(i,:)-a])) < 1e-9, ext(i) = false; end
    end
end

function X = inHull(V, n)
    m = size(V,1);  X = zeros(n,2);
    for i = 1:n
        w = rand(m,1);  w = w/sum(w);
        X(i,:) = w.' * V;
    end
end

function [V, E, F] = randomPolygon()
    while true
        m  = randi([3 5]);
        th = sort(rand(1,m)*2*pi);
        r  = randi([1 5],1,m);
        V  = round([r'.*cos(th'), r'.*sin(th')]);
        V  = uniquetol(V, 1e-9, 'ByRows', true);
        if size(V,1) >= 3 && rank([V - V(1,:); 0 0]) >= 2, break, end
    end
    k = convhull(V(:,1), V(:,2));
    V = V(k(1:end-1), :);
    m = size(V,1);
    E = [(1:m).', [2:m,1].', ones(m,1)];
    F = [ones(m,1), zeros(m,1)];
end
