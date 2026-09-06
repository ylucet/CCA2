function biconjAllHSweep(nRep)
% biconjQ over EVERY Hessian class x domain shape, single- and multi-piece, against the definition.
% The `biconjQ` counterpart of allH-sweep.m, which does this for `conjQ` -- TODO.md 2026-09-05
% named its absence as the biggest hole in biconjQ's evidence.
%
% WHAT envelope-sweep.m ALREADY DOES, and why this is not it: that one is single-piece, bounded,
% and NEGATIVE SEMIDEFINITE only -- the case whose envelope is the lifted hull. It says nothing
% about a positive definite piece, an indefinite one, an unbounded domain, or two pieces at all.
%
% THE MULTI-PIECE FIXTURES ARE CONTINUOUS BY CONSTRUCTION, which is what makes them usable. On a
% mesh split by one edge, put q0 on one side and q0 + c*l on the other, where l is the shared
% edge's own linear form. l vanishes on the edge, so f is continuous there for every c -- and the
% SIGN of c is exactly what decides convexity across the edge (isConvexQ's delta). So one knob
% sweeps convex, affine-jump and non-convex without ever producing a mesh that is not a function.
%
% THE ORACLE is the definition of the convex envelope, in three parts, none of which mentions a
% lower hull:
%   (1) co f <= f on the domain.
%   (2) co f is CONVEX: the midpoint inequality on random pairs of domain points.
%   (3) co f is a MINORANT OF EVERY CONVEX COMBINATION: for points x1..xk in dom and weights lam,
%       co f(sum lam_i x_i) <= sum lam_i f(x_i). This holds for any convex h <= f, so it holds for
%       co f, and it is the part that catches an envelope that is too HIGH -- the direction (1)
%       and (2) together cannot see. It needs no vertex enumeration, so unlike envelope-sweep's
%       vertex-equality check it is valid for an unbounded domain and for a curved piece.
%
% Refusals are TABULATED, not treated as failures: a named refusal is the routine declining to
% answer, which is the honest outcome for the cases still unimplemented. What would be a defect is
% a wrong ANSWER, and that is what the three checks look for.
    if nargin < 1, nRep = 3; end
    rng(20260905);

    % strip1/strip2 exist because every OTHER unbounded shape here recedes along a direction that
    % an indefinite H shears, so before they were added the sweep exercised the unbounded envelope
    % only through its refusals. A half-strip recedes along (0,1) alone, which is in the kernel of
    % a diagonal H, so the lifted recession ray is a fixed direction and the hull is computed.
    shapes = {'square1', 'square2', 'poly1', 'poly2', 'cone1', 'cone2', 'halfplane', ...
              'strip1', 'strip2'};
    % indefXY is x*y, and it is NOT redundant with indef. indef is diag(2,-2), whose recession
    % directions along the axes have d'Hd = +-2, so it never produces a ray of ZERO curvature --
    % the case whose lifted height is a genuine minimum rather than a constant. x*y has d'Hd = 0
    % along both axes, so it is the class that exercises it.
    classes = {'PD', 'PSDsing', 'zero', 'indef', 'indefXY', 'NSDsing', 'ND'};
    tally = containers.Map();          % "shape|class" -> [nOK nRefused nBad]
    refusals = containers.Map();       % identifier -> count
    tot = 0; bad = 0; worst = 0; nOK = 0;

    for rep = 1:nRep
    for si = 1:numel(shapes)
    for ci = 1:numel(classes)
        for c = [-2 0 2]
            H = hessOfClass(classes{ci});
            L = randi([-3 3], 2, 1);  k0 = randi([-3 3]);
            f = buildFixture(shapes{si}, H, L, k0, c);
            if isempty(f) || ~f.isExact(), continue, end
            key = sprintf('%s|%s', shapes{si}, classes{ci});
            if ~isKey(tally, key), tally(key) = [0 0 0]; end
            v = tally(key);

            try
                h = biconjQ(f);
            catch e
                id = e.identifier;
                if ~isKey(refusals, id), refusals(id) = 0; end
                refusals(id) = refusals(id) + 1;
                v(2) = v(2) + 1;  tally(key) = v;
                continue
            end
            nOK = nOK + 1;  v(1) = v(1) + 1;

            X = sampleDomain(f, 200);
            if size(X,1) < 4, tally(key) = v; continue, end
            fx = f.eval(X);  hx = h.eval(X);
            ok = isfinite(fx) & ~isnan(fx) & isfinite(hx);
            X = X(ok,:);  fx = fx(ok);  hx = hx(ok);
            if numel(fx) < 4, tally(key) = v; continue, end
            sc = max(1, abs(fx));

            % (1) co f <= f
            d1 = max([0; (hx - fx)./sc]);

            % (2) convexity of co f
            i1 = randi(numel(fx), 200, 1);  i2 = randi(numel(fx), 200, 1);
            M = (X(i1,:) + X(i2,:))/2;
            hm = h.eval(M);
            d2 = max([0; (hm - (hx(i1)+hx(i2))/2)./max(1,abs(hm))]);

            % (3) co f (sum lam x) <= sum lam f(x)
            d3 = 0;
            for t = 1:60
                k = randi([2 4]);
                idx = randi(numel(fx), k, 1);
                lam = rand(k,1);  lam = lam/sum(lam);
                z = lam.' * X(idx,:);
                hz = h.eval(z);
                rhs = lam.' * fx(idx);
                d3 = max(d3, (hz - rhs)/max(1,abs(rhs)));
            end

            dmax = max([d1 d2 d3]);
            worst = max(worst, dmax);  tot = tot + 1;
            if dmax > 1e-7
                bad = bad + 1;  v(3) = v(3) + 1;
                fprintf(['BAD %-10s %-8s c=%+d  (1) %.2e  (2) %.2e  (3) %.2e  ' ...
                         'H=%s L=%s k0=%d\n'], shapes{si}, classes{ci}, c, d1, d2, d3, ...
                         mat2str(H), mat2str(L.'), k0);
            end
            tally(key) = v;
        end
    end
    end
    end

    fprintf('\n%-12s %-9s %6s %9s %5s\n', 'shape', 'class', 'answered', 'refused', 'bad');
    fprintf('%s\n', repmat('-', 1, 48));
    ks = sort(keys(tally));
    for i = 1:numel(ks)
        v = tally(ks{i});
        parts = split(string(ks{i}), '|');
        fprintf('%-12s %-9s %6d %9d %5d\n', parts(1), parts(2), v(1), v(2), v(3));
    end
    fprintf('\nrefusals:\n');
    rk = sort(keys(refusals));
    for i = 1:numel(rk)
        fprintf('  %-40s %d\n', rk{i}, refusals(rk{i}));
    end
    fprintf('\nBICONJ-ALLH answered %d  checked %d  bad %d  worst %.3e\n', nOK, tot, bad, worst);
end

% ------------------------------------------------------------------------------------------------
function H = hessOfClass(cls)
    switch cls
        case 'PD',      H = [2 0; 0 2];
        case 'PSDsing', H = [2 0; 0 0];
        case 'zero',    H = [0 0; 0 0];
        case 'indef',   H = [2 0; 0 -2];
        case 'indefXY', H = [0 1; 1 0];
        case 'NSDsing', H = [-2 0; 0 0];
        case 'ND',      H = [-2 0; 0 -2];
    end
end

% ------------------------------------------------------------------------------------------------
function f = buildFixture(shape, H, L, k0, c)
% Every multi-piece shape puts q0 on one side of the splitting edge and q0 + c*l on the other,
% with l the edge's own linear form -- continuous for every c, see the header.
    q0 = [0 0 0 0, H(1,1), H(1,2), H(2,2), L(1), L(2), k0];
    f = [];
    switch shape
        case 'square1'      % one piece on the unit square
            V = [0 0; 1 0; 1 1; 0 1];
            f = QuaPol(V, [1 2 1; 2 3 1; 3 4 1; 4 1 1], q0, [1 0; 1 0; 1 0; 1 0]);
        case 'square2'      % split by the diagonal (0,0)-(1,1); l = y - x
            V = [0 0; 1 0; 1 1; 0 1];
            E = [1 2 1; 2 3 1; 3 4 1; 4 1 1; 1 3 1];
            F = [1 0; 1 0; 2 0; 2 0; 2 1];
            f = QuaPol(V, E, [q0; addLin(q0, c, [-1 1 0])], F);
        case 'poly1'        % one piece on a pentagon
            V = [0 0; 2 0; 3 2; 1 3; -1 1];
            n = size(V,1);
            f = QuaPol(V, [(1:n)', [(2:n)';1], ones(n,1)], q0, [ones(n,1), zeros(n,1)]);
        case 'poly2'        % a quadrilateral split by a chord; l = y
            V = [-2 0; 2 0; 1 2; -1 2];
            E = [1 2 1; 2 3 1; 3 4 1; 4 1 1; 1 3 1];
            F = [1 0; 1 0; 2 0; 2 0; 2 1];
            f = QuaPol(V, E, [q0; addLin(q0, c, [0 1 0])], F);
        case 'cone1'        % one piece on the first quadrant
            f = QuaPol([0 0; 1 0; 0 1], [1 2 0; 1 3 0], q0, [1 0; 0 1]);
        case 'cone2'        % the first quadrant split by the diagonal ray; l = y - x
            V = [0 0; 1 0; 1 1; 0 1];
            E = [1 2 0; 1 3 0; 1 4 0];
            F = [1 0; 2 1; 0 2];
            f = QuaPol(V, E, [q0; addLin(q0, c, [-1 1 0])], F);
        case 'strip1'       % one piece on the half-strip {0 <= x <= 1, y >= 0}
            f = QuaPol([0 0; 1 0; 0 1; 1 1], [1 2 1; 1 3 0; 2 4 0], q0, [1 0; 0 1; 1 0]);
        case 'strip2'       % the same half-strip split by the horizontal y = 1; l = y - 1
            V = [0 0; 1 0; 0 1; 1 1; 0 2; 1 2];
            E = [1 2 1; 1 3 1; 2 4 1; 3 4 1; 3 5 0; 4 6 0];
            % Edge 4 runs (0,1) -> (1,1), so its direction is (1,0) and its LEFT is +y -- the
            % UPPER face. F(j,1) is the left face, so face 2 goes in column 1. Written [1 2] at
            % first, which silently swapped the two quadratics and made the sweep report 9 bad
            % answers against a mesh that was not the intended function at all
            % (.claude/strip2MeshProbe.m).
            F = [1 0; 0 1; 1 0; 2 1; 0 2; 2 0];
            f = QuaPol(V, E, [q0; addLin(q0, c, [0 1 -1])], F);
        case 'halfplane'    % the closed upper half-plane split by the ray x=0, y>=0; l = x
            V = [0 0; 0 1; -1 0; 1 0];
            E = [1 2 0; 1 3 0; 1 4 0];
            F = [1 2; 0 1; 2 0];
            f = QuaPol(V, E, [addLin(q0, c, [1 0 0]); q0], F);
    end
end

% ------------------------------------------------------------------------------------------------
function q = addLin(q0, c, lin)
% q0 + c*l with l = lin(1)*x + lin(2)*y + lin(3), in the 10-wide weighted basis (x, y and the
% constant all carry weight 1, so the coefficients go straight in).
    q = q0;
    q(8) = q(8) + c*lin(1);
    q(9) = q(9) + c*lin(2);
    q(10) = q(10) + c*lin(3);
end

% ------------------------------------------------------------------------------------------------
function X = sampleDomain(f, n)
% Points the mesh itself says are in the domain. Rejection sampling over a box, which is the only
% method that works for the unbounded shapes too.
    B = 4;
    P = B*(2*rand(6*n, 2) - 1);
    v = f.eval(P);
    X = P(isfinite(v) & ~isnan(v), :);
    if size(X,1) > n, X = X(1:n,:); end
end
