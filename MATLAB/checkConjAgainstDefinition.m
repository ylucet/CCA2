function R = checkConjAgainstDefinition(nCase, seed, verbose)
% checkConjAgainstDefinition  Randomized end-to-end check of `conj` against its DEFINITION.
%
% objective: one number that says whether the conjugate is right, over a family nobody chose by
%   hand. The per-routine tests pin specific shapes; this asks the only question that matters --
%   is `conj(f)(s)` equal to `sup_{x in dom f} <s,x> - f(x)` -- on random polygons and random
%   quadratics, and reports the worst disagreement.
%
% [input]  nCase   : number of random cases (default 40)
%          seed    : rng seed (default 20260824), so a run is reproducible
%          verbose : print a row per case (default true)
% [output] R : struct array with name, kind, symbolic, worstErr, secs
%
% WHY IT IS A CHECK AND NOT A TEST. It is slow (the reference is a scan plus a pattern search per
% probe point) and its cost scales with the family, so it belongs in a run someone starts, not in a
% bucket. Promote a case to `conjCPLQTest` when it finds something.
%
% THE REFERENCE, and its one-sidedness. `sup` is approached from BELOW by any finite sample, so
% `conj < reference - tol` is a definite defect while `conj > reference` within the sampler's own
% accuracy is expected. The sampler is a coarse scan over the domain followed by a pattern search
% whose direction set includes the polygon's own EDGE directions -- without those a maximiser that
% sits ON an edge cannot be improved by sliding along it, and the reference stalls at the scan's
% accuracy while reporting convergence (measured at 3.2e-3, which is far above any tolerance worth
% asserting). The polygon's VERTICES are always in the sample: for a sup a vertex attains, a grid
% is not a reference at all (DECISIONS.md 2026-08-18).

    if nargin < 1 || isempty(nCase),  nCase = 40; end
    if nargin < 2 || isempty(seed),   seed = 20260824; end
    if nargin < 3 || isempty(verbose), verbose = true; end
    rng(seed);
    R = struct('name', {}, 'kind', {}, 'symbolic', {}, 'worstErr', {}, 'secs', {});
    S = probeGrid();

    for c = 1:nCase
        [W, f6, name] = randomCase(c);
        q = polygonQuaPol(W, f6);
        t0 = tic;
        kd = ''; worst = NaN;
        try
            g = q.conj('cplq');
            kd = g.kind();
            worst = 0;
            for i = 1:size(S,1)
                got = evalAny(g, S(i,:));
                ref = supOverPolygon(S(i,:), W, f6);
                worst = max(worst, got - ref);          % signed: BELOW the sup is the defect
                if got < ref - 1e-7
                    worst = got - ref;
                    break
                end
            end
        catch ME
            kd = ['ERROR:' ME.identifier];
        end
        R(end+1) = struct('name', name, 'kind', kd, ...
            'symbolic', strcmp(kd, 'QuaParCPLQ'), 'worstErr', worst, 'secs', toc(t0)); %#ok<AGROW>
        if verbose
            fprintf('%-34s %-14s sym=%d worst=%+10.3e %7.2fs\n', ...
                name, kd, R(end).symbolic, worst, R(end).secs);
        end
    end

    if verbose
        bad = [R.worstErr] < -1e-7;
        err = arrayfun(@(r) ~isempty(strfind(r.kind, 'ERROR')), R); %#ok<STREMP>
        fprintf('\n%d cases: %d BELOW the sup (defects), %d errored, %d symbolic\n', ...
            numel(R), sum(bad), sum(err), sum([R.symbolic]));
        if any(bad)
            fprintf('DEFECTS:\n');
            for i = find(bad), fprintf('   %s  worst %+.3e\n', R(i).name, R(i).worstErr); end
        end
    end
end

% ================================================================================================

function [W, f6, name] = randomCase(c)
% One random convex polygon and one random quadratic on it. The polygon is built from points on a
% perturbed circle so that it is convex by construction rather than by rejection.
    n = randi([3 6]);
    th = sort(rand(1,n) * 2*pi);
    r  = 0.6 + 0.8*rand(1,n);
    W  = [r'.*cos(th)', r'.*sin(th)'];
    W  = W(convhull(W(:,1), W(:,2)), :);
    W  = W(1:end-1, :);                       % convhull repeats the first point
    if signedArea(W) < 0, W = flipud(W); end
    kind = mod(c, 4);
    switch kind
        case 0, f6 = [2 0 2 0 0 0];                 nm = 'convex';
        case 1, f6 = [0 1 0 0 0 0];                 nm = 'indefinite xy';
        case 2, f6 = [-2 0 -2 0 0 0];               nm = 'concave';
        case 3, f6 = [0 0 0 1 -2 0.5];              nm = 'affine';
    end
    f6(4:6) = f6(4:6) + [randn randn randn];        % a random affine part on every one
    name = sprintf('%d: %d-gon, %s', c, size(W,1), nm);
end

function q = polygonQuaPol(W, f6)
    n = size(W,1);
    E = [(1:n)', mod((1:n),n)'+1, ones(n,1)];
    F = [ones(n,1), zeros(n,1)];
    q = QuaPol(W, E, f6, F);
end

function v = evalAny(g, s)
    if isa(g, 'QuaParCPLQ'), v = evalFunctionNDomain(g.fnd, s); else, v = g.eval(s); end
end

function S = probeGrid()
    R = [0.3 1 2.5 6];
    th = (0:11) * (2*pi/12);
    S = zeros(0,2);
    for r = R, S = [S; r*[cos(th).' sin(th).']]; end %#ok<AGROW>
    S = [S; 0 0];
end

function v = supOverPolygon(s, W, f6)
% See this file's header for why the reference is a scan followed by a pattern search, and why the
% direction set has to contain the polygon's own edge directions.
    n = size(W,1);
    cand = W;
    for i = 1:n
        j = mod(i,n) + 1;
        t = linspace(0, 1, 120).';
        cand = [cand; (1-t)*W(i,:) + t*W(j,:)]; %#ok<AGROW>
    end
    for a = linspace(0, 1, 18)
        for i = 2:n-1
            b = linspace(0, 1, 18).';
            cand = [cand; W(1,:) + a*(W(i,:) - W(1,:)) + a*b*(W(i+1,:) - W(i,:))]; %#ok<AGROW>
        end
    end
    vals = arrayfun(@(k) objAt(s, cand(k,:), f6), (1:size(cand,1))');
    [v, bi] = max(vals);

    D = [1 0; -1 0; 0 1; 0 -1; 1 1; 1 -1; -1 1; -1 -1];
    D = D ./ vecnorm(D, 2, 2);
    for i = 1:n
        e = W(mod(i,n)+1,:) - W(i,:);
        D = [D; e/norm(e); -e/norm(e)]; %#ok<AGROW>
    end
    p = cand(bi,:); step = 1;
    for it = 1:80
        improved = false;
        for k = 1:size(D,1)
            qp = p + step*D(k,:);
            if inHull(qp, W)
                vq = objAt(s, qp, f6);
                if vq > v, v = vq; p = qp; improved = true; end
            end
        end
        if ~improved
            step = step/2;
            if step < 1e-14, break, end
        end
    end
end

function v = objAt(s, x, f6)
    v = s(1)*x(1) + s(2)*x(2) - (f6(1)*x(1)^2/2 + f6(2)*x(1)*x(2) + f6(3)*x(2)^2/2 ...
                                 + f6(4)*x(1) + f6(5)*x(2) + f6(6));
end

function tf = inHull(p, W)
    n = size(W,1); tf = true;
    tol = 1e-12 * (1 + max(abs(W(:))));
    for i = 1:n
        j = mod(i,n) + 1;
        e = W(j,:) - W(i,:);
        if (p - W(i,:)) * [e(2); -e(1)] > tol, tf = false; return, end
    end
end

function a = signedArea(W)
    x = W(:,1); y = W(:,2); n = size(W,1); a = 0;
    for i = 1:n, j = mod(i,n)+1; a = a + (x(i)*y(j) - x(j)*y(i)); end
    a = a/2;
end
