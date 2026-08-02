function out = sweepQuaParEvalAtVertices(seed, nCases, verbose)
% sweepQuaParEvalAtVertices  SEEDED, COMMITTED sweep for the "QuaPar.eval exactly AT a vertex"
%   defect (SUPPORT_MATRIX.md section 7) and for the relative-tolerance fix in QuaPar.eval.
%
% objective: measure how often QuaPar.eval fails to locate a point that lies EXACTLY at a vertex
%   of its own subdivision, under BOTH the current tolerant test and the exact one it replaced,
%   in a single run -- so the fix's effect is a difference between two columns of one table
%   rather than a comparison between two versions of the source.
%
% [input]  seed    : rng seed (REQUIRED -- this exists to be reproducible; see section 0.1)
%          nCases  : number of random subdivisions (default 200)
%          verbose : print a per-case line (default false)
% [output] out     : struct with the counts and rates, and the failing vertices
%
%   >> out = sweepQuaParEvalAtVertices(20260802, 200)
%
% WHY THIS FILE EXISTS. SUPPORT_MATRIX.md section 0.1 requires every quoted measurement to name a
% committed, seeded generator. Three figures in this repository failed that bar, and the most
% consequential was "QuaPar.eval is wrong at ~1.4% of result vertices (5/356)": the sweep behind
% it was a throwaway script with no recorded seed, an attempt to rebuild it assembled no
% qualifying case, and so the fix for the defect could not be checked against the number it was
% meant to fix. This replaces that number with one anybody can re-derive.
%
% SCOPE, STATED PLAINLY. The original figure was measured on the vertices of maxQuaPar RESULTS,
% assembled from random operand pairs. That population cannot be reconstructed. This sweep
% measures the SAME MECHANISM on a population that can be: the vertices of randomly generated
% polyhedral subdivisions, evaluated directly. The mechanism is identical and does not involve
% maxQuaPar at all -- a vertex is a point at which several edge conics should evaluate to exactly
% zero, they come out at +-1e-16 in floating point, and an exact `<= 0` test then admits the point
% to NO face and leaves eval's return at its Inf initialization. So the rate here is comparable in
% kind, not in value, to the retired 1.4%. The CURVED half of the retired claim (~0.8%, 9/1105)
% is NOT covered: generating valid random parabolic subdivisions (orientation invariant,
% b^2-4ac = 0, arcs that are genuine arcs) is a separate piece of work, and manufacturing them by
% running maxQuaPar is exactly the fragile route that failed before.
%
% THE FUNCTION IS THE SAME QUADRATIC ON EVERY FACE, deliberately. That makes the piecewise
% function globally smooth, so every face agrees at every shared vertex and "located by the
% neighbouring face" cannot register as an error. What remains measurable is precisely the defect
% under test -- located by NO face -- with nothing else mixed in.

    if nargin < 1 || isempty(seed)
        error('sweepQuaParEvalAtVertices:seed', ...
            'A seed is REQUIRED: an unseeded sweep is not a reproducible measurement.');
    end
    if nargin < 2 || isempty(nCases), nCases = 200; end
    if nargin < 3, verbose = false; end

    rng(seed);

    nVert = 0; failNew = 0; failOld = 0; badRing = 0;
    failedNew = zeros(0,2); failedOld = zeros(0,2);
    nBuilt = 0;
    for c = 1:nCases
        [q, coef] = randomFan(4 + randi(5));
        if isempty(q), continue; end
        nBuilt = nBuilt + 1;

        for iv = 1:q.nv
            v = q.V(iv,:);
            want = QuaPar.evalPoly(coef, v);
            nVert = nVert + 1;

            vNew = q.eval(v);
            okNew = isfinite(vNew) && abs(vNew - want) <= 1e-9 * max(1, abs(want));
            if ~okNew
                failNew = failNew + 1;
                failedNew(end+1,:) = [c iv]; %#ok<AGROW>
            end

            if ~locatedByExactTest(q, v)
                failOld = failOld + 1;
                failedOld(end+1,:) = [c iv]; %#ok<AGROW>
            end

            % The recorded diagnostic: the geometry is right, only the exactly-at-a-corner
            % tie-break is at issue -- so a ring of radius 1e-8 around the same vertex must
            % evaluate correctly whatever happens AT it.
            th = (0:5)' * (pi/3);
            ring = v + 1e-8 * [cos(th), sin(th)];
            for k = 1:size(ring,1)
                vr = q.eval(ring(k,:));
                if isfinite(vr) && abs(vr - QuaPar.evalPoly(coef, ring(k,:))) > 1e-9
                    badRing = badRing + 1;
                end
            end
        end
        if verbose
            fprintf('case %3d: nv=%d  cumulative fail new=%d old=%d\n', c, q.nv, failNew, failOld);
        end
    end

    out = struct();
    out.seed        = seed;
    out.nCases      = nCases;
    out.nBuilt      = nBuilt;
    out.nVertices   = nVert;
    out.failExact   = failOld;                       % the pre-fix `all(vals <= 0, 2)` test
    out.failTolerant= failNew;                       % QuaPar.eval as it stands today
    out.rateExact   = failOld / max(1, nVert);
    out.rateTolerant= failNew / max(1, nVert);
    out.ringMismatch= badRing;
    out.failedExact = failedOld;
    out.failedTolerant = failedNew;

    fprintf(['sweepQuaParEvalAtVertices(seed=%d, nCases=%d): %d subdivisions, %d vertices\n' ...
             '  located by NO face, EXACT test (pre-fix)    : %d  (%.2f%%)\n' ...
             '  located by NO face, TOLERANT test (current) : %d  (%.2f%%)\n' ...
             '  ring-of-1e-8 mismatches (geometry check)    : %d\n'], ...
        seed, nCases, nBuilt, nVert, out.failExact, 100*out.rateExact, ...
        out.failTolerant, 100*out.rateTolerant, out.ringMismatch);
end

% ================================================================================================
function [q, coef] = randomFan(m)
% A random convex m-gon, fan-triangulated from vertex 1, carrying ONE random quadratic on every
% face. Returns [] if the random point set did not yield at least 4 hull vertices.
    Praw = randn(3*m, 2);
    k = convhull(Praw(:,1), Praw(:,2));
    k(end) = [];                              % convhull repeats the first point
    V = Praw(k,:);
    if size(V,1) < 4
        q = []; coef = []; return
    end
    V = V(1:min(size(V,1), m), :);
    n = size(V,1);
    if n < 4
        q = []; coef = []; return
    end
    % convhull returns the hull in COUNTER-CLOCKWISE order, which is what the F convention below
    % assumes (interior on the left of each directed boundary edge).
    if signedArea(V) < 0
        V = flipud(V);
    end

    E = zeros(0,3); F = zeros(0,2);
    for j = 1:n                                % boundary edges, CCW: interior on the left
        j2 = mod(j, n) + 1;
        if j == 1
            fk = 1;
        elseif j == n
            fk = n - 2;
        else
            fk = j - 1;                        % edge (j, j+1) bounds triangle (1, j, j+1)
        end
        E(end+1,:) = [j j2 1]; %#ok<AGROW>
        F(end+1,:) = [fk 0];   %#ok<AGROW>
    end
    for j = 3:n-1                              % diagonals 1->j, between faces j-2 and j-1
        E(end+1,:) = [1 j 1];  %#ok<AGROW>
        F(end+1,:) = [j-1 j-2]; %#ok<AGROW>
    end

    coef = randn(1,6);
    f = repmat(coef, n-2, 1);
    q = QuaPar(V, E, f, F);
end

% ================================================================================================
function tf = locatedByExactTest(q, x)
% The point location QuaPar.eval performed BEFORE the relative tolerance went in: `all(vals<=0,2)`
% with no tolerance at all. Kept here, not in QuaPar, so that measuring the fix needs no second
% copy of the library and no editing of the source under test.
    tf = false;
    EC = q.edgeConics();
    for i = 1:size(q.P,1)
        Pe = q.P{i};
        if isempty(Pe), continue; end
        vals = zeros(1, numel(Pe));
        for t = 1:numel(Pe)
            cvec = EC(abs(Pe(t)),:);
            sc = max(1, max(abs(cvec)));
            vals(t) = QuaPar.evalConic(cvec, x) * sign(Pe(t)) / sc;
        end
        if all(vals <= 0)
            tf = true;
            return
        end
    end
end

% ================================================================================================
function a = signedArea(V)
    n = size(V,1); a = 0;
    for i = 1:n
        j = mod(i,n)+1;
        a = a + V(i,1)*V(j,2) - V(j,1)*V(i,2);
    end
    a = a/2;
end
