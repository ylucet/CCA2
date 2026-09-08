function cellFaceRegionProbe(varargin)
% WHICH PRIMAL REGION DOES EACH NON-BRIDGE CANDIDATE OWN? The non-bridge analogue of the question
% asked for the bridge faces, and the answer is the opposite one.
%
% WHY IT IS THE SAME QUESTION. `.claude/cellFaceExactProbe.m` names the candidate that wins AT a
% point. An assembly needs the REGION each candidate owns, because that region is the face's
% domain. For the BRIDGE faces the answer was measured on 2026-09-06 and it was CONIC -- line fits
% off by 0.30 to 0.71 where conic fits sat at 1e-14, which is what forced `QuaCon` rather than
% `QuaPar`. The same question has to be asked here rather than assumed to have the same answer.
%
% THE PREDICTION, WHICH IS WHY THE ANSWER SHOULD DIFFER. The winning candidate is decided by where
% the unconstrained maximiser `s*(x) = H^-1 y` falls relative to the cell: the KKT conditions say
% candidate `s0` wins exactly when `y - H s0` lies in the normal cone of `P` at `s0`. Both `y - H s0`
% and that cone are AFFINE in `x`, so each region is an intersection of half-planes -- POLYHEDRAL.
%
% Contrast with the bridge: there the region was a union over a TIE CONIC of segments between two
% affine images, so a curve entered through the tie. Here there is no tie and no curve, so the
% regions should be polyhedral -- and if they are, the non-bridge faces need no new subdivision
% machinery at all, which is a real saving and the reason to measure it.
%
% STATUS 2026-09-07: RESOLVED, THE REGIONS ARE POLYHEDRAL, and the three measurement defects that
% stood between the prediction and that answer are all fixed here (DECISIONS 2026-09-07):
%
%   1. POOLING every crossing of a cell into one fit -- fixed earlier; see the note in
%      traceBoundaries.
%   2. A THRESHOLD BELOW THE MEASUREMENT'S OWN RESOLUTION. The straightness test was a fixed 1e-6
%      while the radial trace step is 4.5e-3, so a perfectly straight boundary located to +-dr/2
%      lands three orders ABOVE the threshold and is reported curved. The threshold is now DERIVED
%      from the trace step, and both numbers are printed, so the comparison can be read.
%   3. ONE POINT WITH TWO NAMES, which is what the "remaining pooled label pair" actually was. A
%      cell VERTEX is an endpoint of the two facets meeting there, and the label spelled the facet
%      it came from, so the same point -- the same affine face -- carried two labels and floating
%      point decided which appeared. The "boundary" between them is not a curve but NOISE OVER A
%      2-D REGION: those were the arcs with 2000 and 4400 points, fitted by neither a line (0.20,
%      0.71) nor a conic (1.6e-2, 5.2e-2), and a conic residual that large should have been the
%      tell. An endpoint is now named by its COORDINATES. The label pair was also made ORDERED and
%      each key split into connected components, so two distinct boundaries can never share a fit.
%
% With all three fixed: 25 arcs, worst line residual 1.3e-3 against a trace step of 4.5e-3, and
% NOTHING above the tolerance. The regions are polyhedral.
%
% HOW IT IS MEASURED, WITHOUT ASSUMING THE ANSWER. Sample the plane, label each point by the
% winning candidate, and for each label trace the boundary of its region along rays; then fit both
% a LINE and a general CONIC to each traced boundary. That is exactly the procedure
% `.claude/bridgeRegionProbe.m` used, so the two answers are comparable -- and it can return
% "curved" if that is what is there.
%
%   cellFaceRegionProbe            % 201 x 201 grid
%   cellFaceRegionProbe(401)

    n = 201;
    if nargin >= 1, n = varargin{1}; end
    verbose = nargin >= 2;                 % per-arc listing: key, points, extent, residual

    % TWO FIXTURES, NOT THREE, AND THE REASON IS THE COST rather than the answer. The hexagon has
    % 13 cells and its widest has 18 candidate regions; tracing 361 rays x 1200 radii through the
    % full candidate enumeration at each point does not finish inside a sensible bound. The
    % triangle and the square already exercise every candidate KIND -- stationary, vertex,
    % endpoint and constant all win somewhere on them -- so the geometric question is answered by
    % them. Widening to the hexagon is a runtime problem, not a coverage one, and it is named here
    % rather than left as a silent truncation.
    fixtures = { 'triangle', tri([0 0; 1 0; 1 1], 0), ...
                 'square',   poly([0 0; 1 0; 1 1; 0 1]) };

    % THE THRESHOLD IS THE MEASUREMENT'S OWN RESOLUTION, not a constant. A crossing is detected
    % between two consecutive radii and placed at the midpoint, so its position carries an error of
    % up to dr/2 whatever the boundary really is; a straight boundary therefore fits a line only to
    % about that. `strTol = dr` is that bound with one factor of two of room, and it is printed
    % beside the worst residual so the comparison is legible rather than asserted. A curved
    % boundary is not near this scale: the bridge regions' line fits missed by 0.30 to 0.71.
    dr = traceStep();
    strTol = dr;

    nB = 0;  worstLine = 0;  bestConicWin = 0;  nCurved = 0;
    for i = 1:2:numel(fixtures)
        nm = fixtures{i};  g = conjQ(fixtures{i+1});
        fprintf('\n=== %s ===\n', nm);
        for k = 1:g.nf
            P = cellRows(g, k);
            if isempty(P), continue, end
            [H, b, c] = cellData(g, k);
            B = traceBoundaries(H, b, c, P, n);
            if isempty(B), continue, end
            worstHere = 0;
            for t = 1:numel(B)
                [eL, eC, sc] = fitBoth(B{t}.pts);
                eAbs = eL * sc;                       % back in the plane's own units
                if verbose
                    fprintf(['    cell %2d arc %2d  key %-8s  %4d pts  extent %.3f  ' ...
                             'line %.2e  conic %.2e\n'], ...
                            k, t, B{t}.pair, size(B{t}.pts,1), sc, eAbs, eC*sc);
                end
                nB = nB + 1;
                worstLine = max(worstLine, eAbs);
                worstHere = max(worstHere, eAbs);
                if eC > 0, bestConicWin = max(bestConicWin, eL / max(eC, 1e-18)); end
                if eAbs > strTol, nCurved = nCurved + 1; end
            end
            fprintf('  cell %2d: %d candidate region(s), %d traced boundary arc(s), worst line %.2e\n', ...
                    k, B{1}.nLab, numel(B), worstHere);
        end
    end

    fprintf('\n%d traced boundary arcs over two fixtures\n', nB);
    fprintf('radial trace step dr = %.3e, so the straightness tolerance is %.3e\n', dr, strTol);
    fprintf('worst LINE residual over all of them: %.3e\n', worstLine);
    fprintf('%d arcs where the line fit was WORSE than the tolerance (i.e. genuinely curved)\n', ...
            nCurved);
    if nCurved == 0
        fprintf(['\nEvery candidate region is bounded by STRAIGHT lines. So the non-bridge faces\n' ...
                 'need no curved subdivision at all -- the opposite of the bridge faces, whose\n' ...
                 'boundaries are conics (DECISIONS 2026-09-06), and the reason is that a bridge\n' ...
                 'region is a union over a tie CONIC while these come from affine KKT conditions.\n']);
    end
end

% ------------------------------------------------------------------------------------------------
function B = traceBoundaries(H, b, c, P, n)
% Label the plane by winning candidate, then trace the boundaries BETWEEN labels along rays from
% the centroid of the labelled set. Each traced arc is returned with its points.
    lo = -6;  hi = 6;
    t  = linspace(lo, hi, n);
    [X1, X2] = meshgrid(t, t);
    L = zeros(numel(X1), 1);
    keys = containers.Map('KeyType','char','ValueType','double');
    for i = 1:numel(X1)
        w = winner(H, b, c, P, [X1(i); X2(i)]);
        if isempty(w), continue, end
        if ~isKey(keys, w), keys(w) = keys.Count + 1; end
        L(i) = keys(w);
    end
    nLab = keys.Count;
    B = {};
    if nLab < 2, return, end

    ctr = [mean(X1(L > 0)); mean(X2(L > 0))];

    % Crossings are grouped BY LABEL PAIR, and that grouping is the whole measurement.
    %
    % A first version POOLED every crossing point of a cell into one set and fitted a single line
    % to it. That was wrong, and wrong in the direction that manufactures a finding: a cell has
    % several DISTINCT boundaries, each between a different pair of candidates, so no single line
    % fits their union and the residual came out 0.36 on every cell -- reported as "22 of 22 arcs
    % genuinely curved" when nothing had been shown about any individual boundary at all. The
    % comment justifying the pooling claimed it was conservative; it was the opposite.
    %
    % Keying on the ORDERED pair (label before, label after) separates the boundaries, since two
    % different boundaries of the same region have different neighbours across them.
    %
    % THE KEY WAS UNORDERED, WHICH IS THE SECOND HALF OF THE SAME DEFECT. `min|max` merges the
    % crossing INTO a region with the crossing OUT of it, so two regions that touch on OPPOSITE
    % sides of a cell shared one key and one fit -- the pooling signature all over again, one label
    % pair at a time, and the residual it produced was the 2218-point outlier of 2026-09-07. The
    % key is now the ordered pair, which is what the paragraph above always claimed it was.
    [rLo, rHi, nR] = traceRadii();
    dr = traceStep();
    byPair = containers.Map('KeyType','char','ValueType','any');
    for a = linspace(0, 2*pi, 361)
        u = [cos(a); sin(a)];
        prev = 0;
        for r = linspace(rLo, rHi, nR)
            x = ctr + r*u;
            w = winner(H, b, c, P, x);
            if isempty(w), cur = 0; elseif ~isKey(keys, w), cur = -1; else, cur = keys(w); end
            if prev > 0 && cur > 0 && cur ~= prev
                key = sprintf('%d>%d', prev, cur);
                pt  = (ctr + (r - dr/2)*u).';       % the crossing is in [r-dr, r]; take the middle
                if isKey(byPair, key), byPair(key) = [byPair(key); pt];
                else,                  byPair(key) = pt;
                end
            end
            prev = cur;
        end
    end

    % AND THEN SPLIT EACH KEY INTO CONNECTED COMPONENTS, because an ordered key is not a proof of
    % one boundary either: a boundary that curls round the centre can be crossed in both senses, and
    % two far-apart straight boundaries would still share a key. Single linkage with a gap of
    % `gapTol` separates them. The gap is not a free parameter: consecutive rays are 2*pi/360 apart,
    % so two crossings of the SAME boundary are at most about rHi*dA apart, and anything beyond a
    % few times that is a different boundary. Splitting too eagerly is safe here -- a subset of a
    % straight boundary is still straight, and a component too small to fit is dropped and counted.
    gapTol = 4 * rHi * (2*pi/360);

    B = {};
    ks = keys_of(byPair);
    for i = 1:numel(ks)
        for comp = components(byPair(ks{i}), gapTol)
            pts = comp{1};
            if size(pts,1) < 8, continue, end      % too few to fit anything; not evidence
            B{end+1} = struct('pts', pts, 'nLab', nLab, 'pair', ks{i}); %#ok<AGROW>
        end
    end
end

function [rLo, rHi, nR] = traceRadii()
% The radial trace, in ONE place, so the straightness tolerance and the crossing midpoint are both
% derived from the same numbers rather than restated as constants.
    rLo = 0.05;  rHi = 5.5;  nR = 1200;
end

function dr = traceStep()
    [rLo, rHi, nR] = traceRadii();
    dr = (rHi - rLo) / (nR - 1);
end

function C = components(pts, gapTol)
% Single-linkage components of a point set at scale `gapTol`, returned as a cell row. Written out
% rather than taken from a toolbox: the sets are a few hundred points, and one fewer dependency in
% a probe is worth more than the speed.
    m = size(pts,1);
    lab = zeros(m,1);  nc = 0;
    for i = 1:m
        if lab(i) > 0, continue, end
        nc = nc + 1;  lab(i) = nc;  stack = i;
        while ~isempty(stack)
            j = stack(end);  stack(end) = [];
            d = vecnorm(pts - pts(j,:), 2, 2);
            nb = find(d <= gapTol & lab == 0);
            lab(nb) = nc;  stack = [stack; nb]; %#ok<AGROW>
        end
    end
    C = cell(1, nc);
    for k = 1:nc, C{k} = pts(lab == k, :); end
end

function ks = keys_of(m)
    ks = keys(m);
end

function w = winner(H, b, c, P, x)
% The winning candidate's LABEL at x, as a string. Empty when the sup is unbounded.
%
% THE LABEL NAMES THE FACE, NOT THE SPELLING, and that distinction was the third measurement defect
% of this probe. An ENDPOINT candidate is a FIXED point of the plane -- a vertex of the cell -- so
% the face it carries is the affine function `<v,x> - g(v)`, and it does not matter which of the two
% facets meeting at `v` produced it. Labelling it `end<facet>_<t>` gave that one point TWO names,
% and where it won, floating point decided which name came out. The "boundary" between those two
% names is then not a curve at all but NOISE SCATTERED OVER A 2-D REGION -- which is what the arcs
% with 2000 and 4400 points were, and why neither a line nor a conic fitted them. Naming an
% endpoint by its COORDINATES merges the two spellings, because they are the same face.
    y = x - b;  best = -inf;  w = '';
    if abs(det(H)) > 1e-12
        s = H \ y;
        if all(P(:,1:2)*s + P(:,3) >= -1e-9)
            best = s.'*y - 0.5*(s.'*H*s) - c;  w = 'stationary';
        end
    end
    for i = 1:size(P,1)
        Ni = P(i,1:2);
        if all(Ni == 0), continue, end
        d = [Ni(2); -Ni(1)];
        p = -P(i,3) * Ni(:) / (Ni * Ni.');
        [lo, hi, empty] = lineRange(P, p, d);
        if empty, continue, end
        a2 = -0.5*(d.'*H*d);  a1 = d.'*y - p.'*H*d;  a0 = p.'*y - 0.5*(p.'*H*p) - c;
        if a2 < -1e-14
            ts = -a1/(2*a2);
            if ts > lo && ts < hi
                v = a0 + a1*ts + a2*ts^2;
                if v > best, best = v; w = sprintf('vertex%d', i); end
            end
        elseif abs(a1) > 1e-12
            if (a1 > 0 && isinf(hi)) || (a1 < 0 && isinf(lo)), w = ''; return, end
        else
            if a0 > best, best = a0; w = sprintf('const%d', i); end
        end
        for q = [lo hi]
            if isinf(q), continue, end
            v = a0 + a1*q + a2*q^2;
            sq = p + q*d;                              % the endpoint itself, which IS the face
            if v > best, best = v; w = sprintf('pt(%.6f,%.6f)', sq(1), sq(2)); end
        end
    end
    if isinf(best), w = ''; end
end

function [lo, hi, empty] = lineRange(P, p, d)
    lo = -inf;  hi = inf;  empty = false;
    for j = 1:size(P,1)
        A = P(j,1:2)*d;  B = P(j,1:2)*p + P(j,3);
        if abs(A) < 1e-14
            if B < -1e-9, empty = true; return, end
        elseif A > 0, lo = max(lo, -B/A);
        else, hi = min(hi, -B/A);
        end
    end
    if lo > hi + 1e-9, empty = true; end
end

function [eL, eC, sc] = fitBoth(Pts)
% Least-squares residuals of the best LINE and the best general CONIC through Pts, normalised by
% the extent so the two are comparable and scale-free -- the same measure bridgeRegionProbe used,
% so the two answers can be read side by side. `sc` is that extent, and it is RETURNED because the
% straightness test needs the residual in the ORIGINAL units: the tolerance comes from the trace
% step, which is a distance in the plane, and comparing it against a scale-free number is how this
% probe reported straight boundaries as curved.
    if size(Pts,1) < 6, eL = 0; eC = 0; sc = 1; return, end
    Pts = Pts - mean(Pts,1);
    sc = max(vecnorm(Pts, 2, 2));
    if sc == 0, eL = 0; eC = 0; sc = 1; return, end
    Pts = Pts / sc;
    ML = [Pts, ones(size(Pts,1),1)];
    eL = min(svd(ML)) / sqrt(size(Pts,1));
    MC = [Pts(:,1).^2, Pts(:,1).*Pts(:,2), Pts(:,2).^2, Pts(:,1), Pts(:,2), ones(size(Pts,1),1)];
    eC = min(svd(MC)) / sqrt(size(Pts,1));
end

% ---- fixtures and helpers ------------------------------------------------------------------------
function [H, b, c] = cellData(g, k)
    n = g.fN(k,:);  d = g.fD(k);
    H = [n(5), n(6); n(6), n(7)] / d;
    b = [n(8); n(9)] / d;
    c = n(10) / d;
end

function rows = cellRows(g, k)
    rows = zeros(0,3);
    FC = g.FC{k};
    for t = 1:size(FC,1)
        cc = g.EcQ(FC(t,1),:);
        if any(cc(1:3) ~= 0), rows = zeros(0,3); return, end
        if FC(t,2) == 0, continue, end
        rows(end+1,:) = FC(t,2) * cc(4:6); %#ok<AGROW>
    end
end

function f = tri(V, c)
    q = [0 0 0 0, 2, 0, 2, 0, 0, 0];
    q(8) = q(8) - c;   q(9) = q(9) + c;
    f = QuaPol(V, [1 2 1; 2 3 1; 3 1 1], q, [1 0; 1 0; 1 0]);
end

function f = poly(V)
    n = size(V,1);
    E = [(1:n).', [2:n, 1].', ones(n,1)];
    f = QuaPol(V, E, [0 0 0 0, 2, 0, 2, 0, 0, 0], repmat([1 0], n, 1));
end
