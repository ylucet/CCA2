% T8, second idea: AXIS-PARALLEL splits in the bilinear frame.
%
% WHY THE CENTROID FAILED. For f = u1*u2 an edge is CONVEX exactly when its direction d has
% d1*d2 > 0 -- the restriction of the bilinear form to the edge is a quadratic with that leading
% coefficient. A cevian from a vertex has a slope BETWEEN the two edges it separates, so on a
% triangle whose edges all have positive slope every cevian has positive slope too, and the
% convex-edge count does not drop. Measured: all three sub-triangles still refused.
%
% WHAT SHOULD WORK. An AXIS-PARALLEL edge has d1 = 0 or d2 = 0, so d1*d2 = 0 and the restriction
% is AFFINE -- not a convex edge. Cutting horizontally through the middle vertex, then vertically
% through each resulting middle vertex, gives sub-triangles carrying one horizontal and one
% vertical edge, leaving at most ONE original (possibly convex) edge each. And the cut
% coordinates are the vertices' OWN coordinates, so the split is rational whenever the input is.
%
% CAVEAT, stated up front: this is rational only when the bilinear FRAME is rational, i.e. when f
% is already of the form b*x*y + linear. A general indefinite quadratic needs M = bilinearFrame(Q),
% whose entries carry sqrt(lambda). So this route removes the surd for the x*y family -- which is
% every fixture in the slow bucket and every QPLIB term -- and not for the general case.
cd(getenv('CCA2DIR'));
warning('off','symbolic:sym:isAlways:TruthUnknown');

E3 = [1 2 1; 2 3 1; 3 1 1]; F3 = [1 0; 1 0; 1 0];
S  = [0 0; 1 0; 0 1; 1 1; -1 1; 1 -1; 2 0.5; -0.5 2; -1 -1];

cases = { ...
    'A.5 3-convex-edge (COAP)', [0 0; 1 1; 3 2]; ...
    'A.5 variant b',            [0 0; 2 1; 1 3]; ...
    'A.5 variant c',            [0 0; 3 1; 1 2]; ...
    'A.4 2-convex-edge',        [2 1; 0 0; 1 0]; ...
    '1-convex-edge',            [0 0; 2 0; 1 1] };
f6 = [0 1 0 0 0 0];        % f = x*y, identity bilinear frame

for c = 1:size(cases,1)
    nm = cases{c,1}; V = cases{c,2};
    tris = axisSplit(V);
    fprintf('\n=== %s ===  split into %d sub-triangles\n', nm, numel(tris));
    g = []; ok = true; why = '';
    for t = 1:numel(tris)
        W = tris{t};
        if abs(det([W(2,:)-W(1,:); W(3,:)-W(1,:)])) < 1e-12, continue, end
        if det([W(2,:)-W(1,:); W(3,:)-W(1,:)]) < 0, W = W([1 3 2],:); end
        try
            gt = conjPieceCPLQ(QuaPol(W, E3, f6, F3));
        catch ME
            ok = false; why = sprintf('sub %d %s: %s', t, mat2str(W), ME.identifier); break
        end
        if isempty(g), g = gt; else, g = maxQuaPar(g, gt); end
    end
    if ~ok
        fprintf('  REFUSED (%s)\n', why);
        continue
    end
    fprintf('  every sub-triangle conjugated DIRECTLY -- no envelope, no A.5, no surd\n');
    worst = 0; bad = 0;
    for i = 1:size(S,1)
        s = S(i,:);
        want = supOverTriangle(f6, V, s);
        got = NaN;
        try, got = g.eval(s); catch, end
        if isnan(got)
            fprintf('    s=(%5.2f,%5.2f) UNCOVERED (want %.6g)\n', s(1), s(2), want); bad=bad+1; continue
        end
        e = abs(got - want); worst = max(worst, e);
        if e > 1e-6 * max(1, abs(want))
            fprintf('    s=(%5.2f,%5.2f) got %.9g want %.9g ERR %.3g\n', s(1), s(2), got, want, e);
            bad = bad + 1;
        end
    end
    fprintf('  vs sup over the ORIGINAL triangle: worst %.3g, %d of %d wrong\n', worst, bad, size(S,1));
end

function tris = axisSplit(V)
% Cut horizontally through the middle-y vertex, then vertically through each piece's middle-x
% vertex. All cut coordinates are existing vertex coordinates, so the result is rational.
    tris = {};
    for T = horizontalCut(V)
        for U = verticalCut(T{1})
            tris{end+1} = U{1}; %#ok<AGROW>
        end
    end
end

function out = horizontalCut(V)
    out = {};
    [~, ord] = sort(V(:,2));
    A = V(ord(1),:); B = V(ord(2),:); C = V(ord(3),:);
    if abs(C(2) - A(2)) < 1e-12, out = {V}; return, end
    t = (B(2) - A(2)) / (C(2) - A(2));
    D = A + t*(C - A);                       % the point on AC at B's height
    out = {[A; B; D], [B; C; D]};
end

function out = verticalCut(V)
    out = {};
    [~, ord] = sort(V(:,1));
    A = V(ord(1),:); B = V(ord(2),:); C = V(ord(3),:);
    if abs(C(1) - A(1)) < 1e-12, out = {V}; return, end
    t = (B(1) - A(1)) / (C(1) - A(1));
    D = A + t*(C - A);
    out = {[A; B; D], [B; C; D]};
end

function v = supOverTriangle(f6, V, s)
    n = 120;
    [a, b] = meshgrid(linspace(0,1,n), linspace(0,1,n));
    keep = (a(:) + b(:)) <= 1;
    L = [a(keep), b(keep)];
    X = V(1,:) + L(:,1).*(V(2,:)-V(1,:)) + L(:,2).*(V(3,:)-V(1,:));
    fv = QuaPol.evalPoly(f6, X);
    v = max(X * s(:) - fv(:));
end
