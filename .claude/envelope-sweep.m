function env_sweep(nCase)
% biconjQ against the DEFINITION of a convex envelope, which is what plqCheck.m checks too:
%
%   (1) co f <= f  on the domain, with EQUALITY at every vertex of the domain. For a concave q the
%       equality at vertices is not optional -- the envelope is the affine interpolant of exactly
%       those values, so a facet chosen wrongly shows up there first.
%   (2) co f is CONVEX (midpoint test on random pairs).
%   (3) co f is the LARGEST convex minorant: no affine function that stays below f on the vertices
%       may exceed co f anywhere. Checked by sampling affine minorants built from the vertex data.
    if nargin < 1, nCase = 40; end
    rng(20260904);
    bad1 = 0; bad2 = 0; bad3 = 0; tot = 0; errs = 0;
    for c = 1:nCase
        [V, E, F] = randomPolygon();
        H = randi([0 4],2,2);  H = -(H + H.');
        if mod(c,4) == 0, H = zeros(2); end
        if ~(H(1,1) <= 0 && H(2,2) <= 0 && (H(1,1)*H(2,2)-H(1,2)^2) >= 0), continue, end
        L = randi([-4 4],2,1);  k0 = randi([-4 4]);
        f = QuaPol(V, E, [0 0 0 0, H(1,1), H(1,2), H(2,2), L(1), L(2), k0], F);
        if ~f.isExact(), continue, end
        try
            h = biconjQ(f);
        catch e
            errs = errs + 1;  fprintf('case %2d ERROR %s\n', c, e.identifier);  continue
        end

        % (1) co f <= f inside, equality at the vertices
        X = samplePolygon(V, 300);
        qv = arrayfun(@(i) 0.5*V(i,:)*H*V(i,:).' + L.'*V(i,:).' + k0, 1:size(V,1)).';
        % Equality is required at the EXTREME points, not at every listed vertex. convhull can
        % return a point that lies ON the segment between its neighbours, and for a concave q such
        % a vertex is DOMINATED -- its lifted point is strictly above the chord, so it is not on
        % the lower hull and co f is strictly below q there. That is the envelope being right, and
        % it is the same phenomenon conjQ's aDominatedVertexContributesNoCellAtAll pins.
        hv = h.eval(V);
        ext = true(size(V,1),1);
        for i = 1:size(V,1)
            a = V(mod(i-2,size(V,1))+1,:);  b = V(mod(i,size(V,1))+1,:);
            if abs(det([b-a; V(i,:)-a])) < 1e-9, ext(i) = false; end
        end
        nVert = sum(abs(hv(ext) - qv(ext)) > 1e-9*max(1,abs(qv(ext))));
        hx = h.eval(X);
        qx = 0.5*sum((X*H).*X, 2) + X*L + k0;
        nAbove = sum(hx > qx + 1e-9*max(1,abs(qx)));
        nInf   = sum(isinf(hx));
        if nVert > 0 || nAbove > 0
            bad1 = bad1 + 1;
            fprintf('case %2d: %d/%d vert off, %d/%d above, %d Inf | V=%s H=%s L=%s k0=%d\n', c, nVert, numel(qv), nAbove, numel(qx), nInf, mat2str(V), mat2str(H), mat2str(L), k0);
        end

        % (2) convexity
        A = samplePolygon(V, 200);  B = samplePolygon(V, 200);
        M = (A+B)/2;
        if any(h.eval(M) > (h.eval(A) + h.eval(B))/2 + 1e-9), bad2 = bad2 + 1; end

        % (3) largest convex minorant: an affine function below f AT THE VERTICES is below co f
        for t = 1:20
            g = randn(1,2);
            b = min(qv.' - V*g.');                 % the highest such affine minorant
            aff = X*g.' + b;
            if any(aff > hx + 1e-9*max(1,abs(hx))), bad3 = bad3 + 1; break, end
        end
        tot = tot + 1;
    end
    fprintf('ENV cases %d errors %d   minorant-or-vertex %d   nonconvex %d   not-largest %d\n', ...
        tot, errs, bad1, bad2, bad3);
end

function X = samplePolygon(V, n)
    m = size(V,1);  X = zeros(n,2);
    for i = 1:n
        w = rand(m,1);  w = w/sum(w);
        X(i,:) = w.' * V;
    end
end

function [V, E, F] = randomPolygon()
    while true
        m  = randi([3 6]);
        th = sort(rand(1,m)*2*pi);
        r  = randi([1 6],1,m);
        V  = round([r'.*cos(th'), r'.*sin(th')]*2)/2;
        V  = uniquetol(V, 1e-9, 'ByRows', true);
        % the rounding can collapse points onto each other or onto one line, and convhull then
        % raises rather than returning a degenerate hull -- redraw instead
        if size(V,1) >= 3 && rank([V - V(1,:); 0 0]) >= 2, break, end
    end
    k  = convhull(V(:,1), V(:,2));
    V  = V(k(1:end-1), :);
    m  = size(V,1);
    E  = [(1:m).', [2:m,1].', ones(m,1)];
    F  = [ones(m,1), zeros(m,1)];
end
