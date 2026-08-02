function checkBiconjDomainCoverage()
% checkBiconjDomainCoverage  Does biconj('cplq') work for ANY domain?  Run it and read the table.
%
% RESULT as of 2026-08-02 (see SUPPORT_MATRIX.md section 8 for the standing summary):
%   OK     triangle x {convex, indefinite, concave, affine}
%   OK     box as ONE face, f = x*y (McCormick) and f = 0 (indicator)
%   OK     AXIS-ALIGNED box [0,2]x[0,3], one face
%   OK     unbounded: full-domain quadratic, 4 cones (|x|+|y|), 3 wedges (max(0,x,y))
%   WRONG  box as TWO faces sharing a diagonal
%   ERROR  general convex quadrilateral (one face)
%   ERROR  parallelogram (one face)
%
% READ WHICH CONJUGATION FAILS before concluding anything: the general quadrilateral fails in
% the FIRST, the parallelogram and the two-face box in the SECOND. They are different defects.
%
% The quadrilateral one is a WIRING gap, not a missing algorithm. Two Step 1s exist here:
% convEnvCPLQ.m (CCA2's own) HAS the 3-convex-edge case -- splitThreeConvex cuts the triangle
% through the middle vertex into two 2-convex-edge sub-triangles, [COAP] A.5, which is the
% re-triangulation a general polyhedral set needs -- while the vendored plq_1p.convexEnvelope1
% branches on nCE == 0,1,2 and then falls off the end, leaving `envelope` EMPTY at nCE == 3.
% conjCPLQ's Case C drives Step 1 through the VENDORED one, so CCA2's own split is unreachable
% from conj/biconj on a multi-vertex domain. nCE counts edges of positive finite SLOPE, which is
% why an axis-aligned box (all pieces nCE = 0) sails through and a sheared one does not -- the
% symptom is orientation-dependent, the missing branch is not.
%  One case per domain FAMILY, each checked against a
% ground truth that owes nothing to the conjugate pipeline.
%
%   BOUNDED domains: f** = conv f, computed as the lower convex hull of the sampled graph
%     {(x, y, f(x,y))} via convhulln -- the lower envelope is the max of its lower facets' planes.
%   UNBOUNDED domains: every case here is CONVEX, so f** = f and the check is an identity.
%
% Reports, per case: OK / WRONG (with the worst error and where) / ERROR (identifier) / and
% whether the +inf-outside-the-domain behaviour is right for bounded domains.

    cases = buildCases();
    fprintf('\n%-42s %-10s %-9s %-11s %s\n', 'case', 'domain', 'verdict', 'worst err', 'note');
    fprintf('%s\n', repmat('-', 1, 100));
    for k = 1:numel(cases)
        c = cases(k);
        t0 = tic;
        try
            h = c.p.biconj('cplq');
        catch ME
            fprintf('%-42s %-10s %-9s %-11s %s\n', c.name, c.kind, 'ERROR', '-', ME.identifier);
            continue
        end
        el = toc(t0);
        [verdict, worst, at, note] = checkCase(c, h);
        fprintf('%-42s %-10s %-9s %-11.3g %s (%.0fs, %s)\n', c.name, c.kind, verdict, worst, ...
            note, el, class(h));
        if ~strcmp(verdict, 'OK')
            fprintf('%-42s %s\n', '', sprintf('worst at (%g,%g)', at(1), at(2)));
        end
    end
    fprintf('\n');
end

% ================================================================================================
function [verdict, worst, at, note] = checkCase(c, h)
    worst = 0; at = [0 0]; note = '';
    if strcmp(c.kind, 'unbounded')
        % convex by construction: f** = f
        bad = 0;
        for i = 1:size(c.X,1)
            got = evalAny(h, c.X(i,:));
            want = c.f(c.X(i,1), c.X(i,2));
            e = abs(got - want);
            if ~isfinite(e), e = inf; end
            if e > worst, worst = e; at = c.X(i,:); end
            if e > 1e-7, bad = bad + 1; end
        end
        if bad == 0, verdict = 'OK'; else, verdict = 'WRONG'; end
        note = sprintf('%d/%d pts', size(c.X,1)-bad, size(c.X,1));
        return
    end

    % bounded: ground truth = lower convex envelope of the sampled graph
    env = lowerEnvelope(c.S, c.V);
    bad = 0;
    for i = 1:size(c.X,1)
        got = evalAny(h, c.X(i,:));
        want = env(c.X(i,:));
        e = abs(got - want);
        if ~isfinite(e), e = inf; end
        if e > worst, worst = e; at = c.X(i,:); end
        if e > 2e-3, bad = bad + 1; end
    end
    % and +inf strictly outside the domain
    nOut = 0; badOut = 0;
    if isfield(c, 'Xout') && ~isempty(c.Xout)
        for i = 1:size(c.Xout,1)
            nOut = nOut + 1;
            got = evalAny(h, c.Xout(i,:));
            if isfinite(got), badOut = badOut + 1; end
        end
    end
    if bad == 0 && badOut == 0
        verdict = 'OK';
    else
        verdict = 'WRONG';
    end
    note = sprintf('%d/%d in, %d/%d out', size(c.X,1)-bad, size(c.X,1), nOut-badOut, nOut);
end

% ================================================================================================
function v = evalAny(h, pt)
    if isa(h, 'QuaParCPLQ')
        v = evalFunctionNDomain(h.fnd, pt);
    else
        v = h.eval(pt);
    end
    if isnan(v), v = inf; end
end

% ================================================================================================
function fn = lowerEnvelope(S, F)
% conv of the sampled graph: the lower convex envelope is the pointwise MAX of the planes of the
% hull facets that lie weakly below every sample. S is n-by-2, F is n-by-1.
    P = [S, F(:)];
    % An AFFINE (or constant) f has a planar graph, which convhulln rejects as degenerate --
    % and for such an f the convex envelope is f itself, so fit the plane and be done.
    A = [S, ones(size(S,1),1)];
    cf = A \ F(:);
    if max(abs(A*cf - F(:))) < 1e-12 * max(1, max(abs(F)))
        fn = @(x) cf(1)*x(1) + cf(2)*x(2) + cf(3);
        return
    end
    K = convhulln(P);
    planes = zeros(0,3);
    for t = 1:size(K,1)
        A = P(K(t,1),:); B = P(K(t,2),:); C = P(K(t,3),:);
        n = cross(B-A, C-A);
        if abs(n(3)) < 1e-10, continue; end          % vertical facet: not a graph of a plane
        % z = a x + b y + c
        a = -n(1)/n(3); b = -n(2)/n(3); cc = A(3) - a*A(1) - b*A(2);
        vals = a*S(:,1) + b*S(:,2) + cc;
        if all(F(:) >= vals - 1e-9*max(1, max(abs(F))))
            planes(end+1,:) = [a b cc]; %#ok<AGROW>
        end
    end
    fn = @(x) max(planes(:,1)*x(1) + planes(:,2)*x(2) + planes(:,3));
end

% ================================================================================================
function cases = buildCases()
    cases = struct('name', {}, 'kind', {}, 'p', {}, 'X', {}, 'Xout', {}, 'S', {}, 'V', {}, 'f', {});

    E3 = [1 2 1; 2 3 1; 3 1 1]; F3 = [1 0; 1 0; 1 0];

    % --- bounded TRIANGLE, convex quadratic -----------------------------------------------
    T = [0 0; 1 0; 0 1];
    cases(end+1) = mkBounded('triangle, convex q=(x^2+y^2)/2', ...
        QuaPol(T, E3, [1 0 1 0 0 0], F3), T, @(x,y) 0.5*(x.^2+y.^2), ...
        [0.2 0.2; 0.5 0.25; 1/3 1/3; 0 0], [1.5 0.5; -0.5 0.2]);

    % --- bounded TRIANGLE, indefinite (the classic) ---------------------------------------
    cases(end+1) = mkBounded('triangle, indefinite f=x*y', ...
        QuaPol(T, E3, [0 1 0 0 0 0], F3), T, @(x,y) x.*y, ...
        [0.2 0.2; 0.3 0.4; 0.25 0.5; 0.1 0.6], [1.5 0.5; -0.5 0.2]);

    % --- bounded TRIANGLE, concave --------------------------------------------------------
    cases(end+1) = mkBounded('triangle, concave q=-(x^2+y^2)/2', ...
        QuaPol(T, E3, [-1 0 -1 0 0 0], F3), T, @(x,y) -0.5*(x.^2+y.^2), ...
        [0.2 0.2; 0.3 0.4; 1/3 1/3], [1.5 0.5]);

    % --- bounded TRIANGLE, affine ---------------------------------------------------------
    cases(end+1) = mkBounded('triangle, affine f=2x-3y+1', ...
        QuaPol(T, E3, [0 0 0 2 -3 1], F3), T, @(x,y) 2*x-3*y+1, ...
        [0.2 0.2; 0.3 0.4; 0 0], [1.5 0.5]);

    % --- BOX as ONE face, bilinear (McCormick) --------------------------------------------
    B = [0 0; 1 0; 1 1; 0 1];
    Eb = [1 2 1; 2 3 1; 3 4 1; 4 1 1]; Fb = [1 0; 1 0; 1 0; 1 0];
    cases(end+1) = mkBounded('box (1 face), f=x*y  [McCormick]', ...
        QuaPol(B, Eb, [0 1 0 0 0 0], Fb), B, @(x,y) x.*y, ...
        [0.75 0.25; 0.5 0.5; 0.25 0.25; 0.9 0.6; 0.2 0.8], [1.5 0.5; -0.5 0.5]);

    % --- BOX as ONE face, indicator -------------------------------------------------------
    cases(end+1) = mkBounded('box (1 face), indicator f=0', ...
        QuaPol(B, Eb, [0 0 0 0 0 0], Fb), B, @(x,y) 0*x, ...
        [0.5 0.5; 0 0; 1 1; 0.25 0.75], [1.5 0.5; -0.5 0.5; 2 2]);

    % --- BOX as TWO faces (the known open case) -------------------------------------------
    E2 = [1 2 1; 2 3 1; 1 3 1; 3 4 1; 4 1 1]; F2 = [1 0; 1 0; 2 1; 2 0; 2 0];
    cases(end+1) = mkBounded('box (2 faces, diagonal), f=x*y', ...
        QuaPol([0 0; 1 0; 1 1; 0 1], E2, [0 1 0 0 0 0; 0 1 0 0 0 0], F2), B, @(x,y) x.*y, ...
        [0.75 0.25; 0.5 0.5; 0.25 0.25; 0.9 0.6], [1.5 0.5]);

    % --- non-square QUADRILATERAL, one face -----------------------------------------------
    Q = [0 0; 2 0; 2.5 1.5; 0.5 1];
    Eq = [1 2 1; 2 3 1; 3 4 1; 4 1 1]; Fq = [1 0; 1 0; 1 0; 1 0];
    cases(end+1) = mkBounded('quadrilateral (1 face), f=x*y', ...
        QuaPol(Q, Eq, [0 1 0 0 0 0], Fq), Q, @(x,y) x.*y, ...
        [1 0.5; 1.5 0.8; 0.8 0.5], [3 0.2; -0.5 0.2]);

    % --- UNBOUNDED: full-domain convex quadratic (Case A) ---------------------------------
    cases(end+1) = mkUnbounded('full domain, q=(x^2+y^2)/2 (Case A)', QuaPol.energy(), ...
        @(x,y) 0.5*(x.^2+y.^2), [0 0; 1 1; -2 3; 0.5 -1.5; 4 4]);

    % --- UNBOUNDED: 4 cones, f = |x|+|y| ---------------------------------------------------
    cases(end+1) = mkUnbounded('unbounded 4 cones, f=|x|+|y|', QuaPol.oneNorm(), ...
        @(x,y) abs(x)+abs(y), [0 0; 1 0; 0.5 0.5; -1 2; 2 -3; 0.25 0.75; -2 -2]);

    % --- UNBOUNDED: 3 wedges, f = max(0,x,y) ----------------------------------------------
    Vw = [0 0; -1 0; 1 1; 0 -1];
    Ew = [1 2 0; 1 3 0; 1 4 0]; Fw = [1 2; 2 3; 3 1];
    fw = [0 0 0 0 0 0; 0 0 0 0 1 0; 0 0 0 1 0 0];
    cases(end+1) = mkUnbounded('unbounded 3 wedges, f=max(0,x,y)', QuaPol(Vw, Ew, fw, Fw), ...
        @(x,y) max(max(0,x),y), [0 0; 2 1; -1 3; -2 -3; 0.5 0.25]);
end

% ================================================================================================
function s = mkBounded(name, p, poly, f, X, Xout)
    n = 46;
    [uu,vv] = meshgrid(linspace(min(poly(:,1))-0.001, max(poly(:,1))+0.001, n), ...
                       linspace(min(poly(:,2))-0.001, max(poly(:,2))+0.001, n));
    in = inpolygon(uu(:), vv(:), poly(:,1), poly(:,2));
    S = [uu(in), vv(in)];
    S = [S; poly];                                   % vertices matter for the envelope
    s = struct('name', name, 'kind', 'bounded', 'p', p, 'X', X, 'Xout', Xout, ...
               'S', S, 'V', f(S(:,1), S(:,2)), 'f', f);
end

function s = mkUnbounded(name, p, f, X)
    s = struct('name', name, 'kind', 'unbounded', 'p', p, 'X', X, 'Xout', [], ...
               'S', [], 'V', [], 'f', f);
end
