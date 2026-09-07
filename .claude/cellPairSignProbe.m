function cellPairSignProbe(varargin)
% THE CEILING ON FILTER 2, measured before it is built.
%
% .claude/cellPairFilterProbe.m settles filter 1 -- do the two dual cells MEET, and meet in two
% dimensions -- and its answer is that the surviving count is a fifth to a third of the product,
% quadratic in the cell count. That is not a change of order, so the order has to come from filter
% 2: on the cell intersection, does g1 - g2 actually CHANGE SIGN? No sign change, no tie, no
% bridge, and the elimination never runs.
%
% BUILDING FILTER 2 EXACTLY IS REAL WORK. Deciding whether a quadratic changes sign on a 2-D
% polyhedron needs the polyhedron's vertices, edges and recession directions from its half-plane
% form, plus the exact extremum of a quadratic on each. That is a bounded job but not a small one,
% and it is only worth doing if the prune rate is worth having. So this measures the CEILING first.
%
% WHAT THIS MEASURES, AND WHY THE NUMBER IS OPTIMISTIC. The conjugates are sampled on a grid;
% QuaCon.eval reports both the value and WHICH cell supplied it, so every sample lands in one
% (cellA, cellB) bucket and contributes the sign of g1 - g2 there. A bucket that shows both signs
% definitely ties. A bucket that shows one sign MIGHT still tie somewhere the grid missed -- so the
% "one sign only" count is an UPPER bound on what filter 2 can remove, never a promise. If even the
% upper bound is small, filter 2 is not worth building; only if it is large does the exact version
% earn its cost. Stated this way round on purpose: an optimistic estimate that says "do not build"
% is a sound conclusion, and one that says "build" still has to be confirmed by the real thing.
%
% Bucket counts are also reported against filter 1's full-dimensional count, because a bucket the
% grid never reaches is not evidence of anything and must not be read as a prune.
%
%   cellPairSignProbe            % default 401 x 401 grid
%   cellPairSignProbe(801)       % finer, to see whether "one sign only" survives refinement

    n = 401;
    if nargin >= 1, n = varargin{1}; end

    fprintf('grid %d x %d over the dual box each pair needs\n\n', n, n);
    fprintf('%-34s %8s %8s %8s %8s %8s\n', 'piece pair', 'full-dim', 'sampled', 'both', ...
            'one-sign', 'ceiling');
    fprintf('%s\n', repmat('-', 1, 82));

    cases = fixtures();
    for i = 1:numel(cases)
        nm = cases{i}{1};  fA = cases{i}{2};  fB = cases{i}{3};
        gA = conjQ(fA);  gB = conjQ(fB);
        nFull = countFullDim(gA, gB);
        [nSamp, nBoth, nOne] = signBuckets(gA, gB, n);
        if nSamp == 0
            ceil_ = NaN;
        else
            ceil_ = 100 * nOne / nSamp;
        end
        fprintf('%-34s %8d %8d %8d %8d %7.0f%%\n', nm, nFull, nSamp, nBoth, nOne, ceil_);
    end

    fprintf(['\n"sampled" is how many of filter 1''s pairs the grid actually reached; a pair it\n' ...
             'missed is not evidence. "ceiling" is one-sign / sampled -- the MOST filter 2 could\n' ...
             'remove, since a finer grid can only turn one-sign buckets into both-sign ones.\n' ...
             'Re-run at a finer grid: a ceiling that drops is a ceiling that was grid noise.\n']);
end

% ------------------------------------------------------------------------------------------------
function [nSampled, nBoth, nOneSign] = signBuckets(gA, gB, n)
% Sample the plane, bucket by (cellA, cellB), and ask which buckets see BOTH signs of gA - gB.
%
% The box is chosen from the two conjugates' own named vertices, padded, so the grid covers every
% bounded cell and reaches into the unbounded ones rather than being a fixed guess.
    [lo, hi] = dualBox(gA, gB);
    t1 = linspace(lo(1), hi(1), n);
    t2 = linspace(lo(2), hi(2), n);
    [S1, S2] = meshgrid(t1, t2);
    X = [S1(:), S2(:)];

    [vA, iA] = gA.eval(X);
    [vB, iB] = gB.eval(X);
    ok = (iA > 0) & (iB > 0) & isfinite(vA) & isfinite(vB);
    d  = vA(ok) - vB(ok);
    key = sub2ind([gA.nf, gB.nf], iA(ok), iB(ok));

    % A sign is only counted when it is clear of the evaluation's own noise: gA and gB are computed
    % in doubles (QuaCon.eval is deliberately not exact), so a difference at rounding scale is not a
    % sign. Calling it one would report a tie where there is none, which is the expensive error --
    % it puts an elimination back that filter 2 was meant to remove, never removes a real one.
    % POINTWISE, not global. A single scale taken over the whole grid is the box's largest value,
    % which near the box corners is enormous; using it as the noise floor everywhere would suppress
    % genuine sign changes in the small-value cells and INFLATE the ceiling -- i.e. it would
    % manufacture the answer this probe exists to test. The floor is therefore each point's own.
    tol = sqrt(eps) * max(1, max(abs(vA(ok)), abs(vB(ok))));
    pos = accumarray(key, double(d >  tol), [gA.nf * gB.nf, 1], @max, 0);
    neg = accumarray(key, double(d < -tol), [gA.nf * gB.nf, 1], @max, 0);
    seen = accumarray(key, 1, [gA.nf * gB.nf, 1], @max, 0) > 0;

    nSampled = sum(seen);
    nBoth    = sum(seen & pos > 0 & neg > 0);
    nOneSign = nSampled - nBoth;
end

function [lo, hi] = dualBox(gA, gB)
% A box containing both conjugates' named vertices, padded by half its own width so the unbounded
% cells are entered rather than merely touched.
    V = [vertsOf(gA); vertsOf(gB)];
    if isempty(V), lo = [-4 -4]; hi = [4 4]; return, end
    lo = min(V, [], 1);  hi = max(V, [], 1);
    w  = max(hi - lo, 1);
    lo = lo - w;  hi = hi + w;
end

function V = vertsOf(g)
% QuaCon.vertexCoords realises the NAMED vertices as doubles on demand -- the class deliberately
% does not store them, so this is the one call that turns the exact mesh into a bounding box.
    V = g.vertexCoords();
    V = V(all(isfinite(V), 2), :);
end

function nFull = countFullDim(gA, gB)
% Filter 1's answer, recomputed here so the two columns are from the same run.
    nFull = 0;
    for i = 1:gA.nf
        RA = cellRows(gA, i);
        if isempty(RA), continue, end
        for j = 1:gB.nf
            RB = cellRows(gB, j);
            if isempty(RB), continue, end
            if ratQ.feasible2([RA; RB], true), nFull = nFull + 1; end
        end
    end
end

function rows = cellRows(g, k)
    rows = zeros(0,3);
    FC = g.FC{k};
    for r = 1:size(FC,1)
        e = FC(r,1);  sgn = FC(r,2);
        c = g.EcQ(e,:);
        if any(c(1:3) ~= 0), rows = zeros(0,3); return, end
        if sgn == 0, continue, end
        rows(end+1,:) = sgn * c(4:6); %#ok<AGROW>
    end
end

% ------------------------------------------------------------------------------------------------
function cases = fixtures()
% The same ladder cellPairFilterProbe uses, so the two probes' columns line up.
    cases = { ...
        {'two triangles (shifted)',      tri([0 0; 1 0; 0 1]), tri([2 0; 3 0; 2 1])}, ...
        {'two triangles (cones overlap)', tri([0 0; 1 0; 0 1]), tri([0 0; -1 0; 0 -1])}, ...
        {'triangle vs square',           tri([0 0; 1 0; 0 1]), poly([0 0; 1 0; 1 1; 0 1], 0)}, ...
        {'square vs square',             poly([0 0; 1 0; 1 1; 0 1], 0), ...
                                         poly([0 0; 1 0; 1 1; 0 1], 3)}, ...
        {'hexagon vs hexagon',           poly(hexV(), 0), poly(hexV(), 6)}, ...
        {'octagon vs octagon',           poly(latV(1), 0), poly(latV(1), 8)}, ...
        {'16-gon vs 16-gon',             poly(latV(2), 0), poly(latV(2), 20)} };
end

function V = latV(k)
    D = [];
    for a = -k:k
        for b = -k:k
            if (a ~= 0 || b ~= 0) && gcd(abs(a), abs(b)) == 1
                D(end+1,:) = [a b]; %#ok<AGROW>
            end
        end
    end
    [~, ord] = sort(atan2(D(:,2), D(:,1)));
    V = cumsum([0 0; D(ord,:)]);
    V(end,:) = [];
end

function V = hexV()
    V = [0 0; 2 0; 3 2; 2 4; 0 4; -1 2];
end

function f = tri(V)
    f = QuaPol(V, [1 2 1; 2 3 1; 3 1 1], [0 0 0 0, 2, 0, 2, 0, 0, 0], [1 0; 1 0; 1 0]);
end

function f = poly(V, dx)
    n = size(V,1);
    V(:,1) = V(:,1) + dx;
    E = [(1:n).', [2:n, 1].', ones(n,1)];
    f = QuaPol(V, E, [0 0 0 0, 2, 0, 2, 0, 0, 0], repmat([1 0], n, 1));
end
