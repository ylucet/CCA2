classdef conjDefinitionAuditTest < matlab.unittest.TestCase
% [normal, ~1-2 min] conj (default engine) against the DEFINITION f*(s) = sup_{x in dom f} <s,x> - f(x), on
% every QuaPol example the engine answers (conjQExamplesSweepTest says which), the erratum
% examples, and two unbounded fixtures.
%
% THE ORACLE does not use QuaPol.eval (which mis-reads rays and non-convex boundaries -- TODO).
% Each piece's boundary is walked from V/E/F, rays truncated at radius R, the gap closed by an arc;
% the sup over that polygon is exact up to rounding (every edge in closed form, plus the interior
% critical point of a strictly convex q). Two radii decide +infinity: a sup that still grows from
% R1 to R2 is +infinity. The arc can overshoot a face by O(|vertex|/R), so R is large (measured on
% examples2{9}: 1.5e-2 off at R = 400, exact at 4000).
%
% THE POINTS include an integer half-grid, not only random s: a thin dom f* (a line or a point,
% e.g. f = 0 on a half-plane) is never hit by a random point, and missing it is how two wrong
% answers went unseen until 2026-10-11.

    methods (Test)
        function everyAnswerMatchesTheDefinition(t)
            [F, names] = conjDefinitionAuditTest.inputs();
            R1 = 4000;  R2 = 8000;  nChecked = 0;
            for i = 1:numel(F)
                try, g = F{i}.conj(); catch, continue, end      % refusals: conjQExamplesSweepTest
                nChecked = nChecked + 1;
                rng(i);  [ga, gb] = meshgrid(-3:3);
                S = [randn(30,2)*2; randn(10,2)*8; ga(:) gb(:)];
                v = g.eval(S);  want = zeros(size(S,1),1);
                for j = 1:size(S,1)
                    a = supDef(F{i}, S(j,:), R1);  b = supDef(F{i}, S(j,:), R2);
                    if b > a + 1e-7*max(1,abs(a)), want(j) = inf; else, want(j) = b; end
                end
                t.verifyEqual(isinf(v), isinf(want), sprintf('%s: where f* is finite', names{i}));
                fin = isfinite(v) & isfinite(want);
                t.verifyEqual(v(fin), want(fin), 'AbsTol', 1e-7, 'RelTol', 1e-7, ...
                    sprintf('%s: the value', names{i}));
            end
            t.verifyGreaterThan(nChecked, 30, 'the audit must reach the inputs');
        end
    end

    methods (Static)
        function [F, names] = inputs()
            sets = {'examples','examples2','examplesNonconvex','examplesDiscontinuous'};
            F = {};  names = {};
            for si = 1:numel(sets)
                P = QuaPol.(sets{si})();
                for k = 1:numel(P), F{end+1} = P{k}; names{end+1} = sprintf('%s{%d}', sets{si}, k); end %#ok<AGROW>
            end
            F = [F, {erratumExamplesTest.example1(), erratumExamplesTest.example2(), ...
                     erratumExamplesTest.example3()}];
            names = [names, {'erratum1', 'erratum2', 'erratum3'}];
            F{end+1} = QuaPol([0 -2; 1 -2; 1 -1], [1 2 0; 1 3 0], [0 0 0 0, 1, 0, -1, 0, -2, -2], [1 0; 0 1]);
            names{end+1} = 'shiftedCone';
            F{end+1} = QuaPol([0 0; 1 0; 0 1; -1 0], [1 2 0; 1 3 0; 1 4 0], ...
                              [0 0 0 0, 2, 0, 2, 0, 0, 0; 0 0 0 0, 0, -2, 2, 0, 0, 0], [1 0; 2 1; 0 2]);
            names{end+1} = 'quadrants';
        end
    end
end

% ==================================================================================================
function v = supDef(f, s, R)
% sup over dom f of <s,x> - f(x), with every unbounded piece truncated at radius ~R.
  V = f.V; E = f.E; v = -inf;
  if f.nv == 0                                   % one quadratic on the whole plane
    v = supPoly(f.f(1,:), circ(R), s); return
  end
  if all(f.F(:) == 0) || f.nf == 0               % dimension < 2: the edges are the domain
    if isempty(E), v = s*V(1,:).' - qv(f.f(1,:), V(1,:)); return, end
    for e = 1:size(E,1)
      row = f.f(min(e, size(f.f,1)),:);
      [a, b] = segOf(V, E(e,:), R);
      v = max(v, supSeg(row, a, b, s));
    end
    return
  end
  for k = 1:f.nf
    v = max(v, supPoly(f.f(k,:), faceLoop(f, k, R), s));
  end
end

function [a, b] = segOf(V, e, R)
  a = V(e(1),:);
  if e(3) == 0, d = V(e(2),:) - a; b = a + d/norm(d)*R; else, b = V(e(2),:); end
end

function L = faceLoop(f, k, R)
% Face k's boundary as a closed polygon: rays truncated at length R, the gap closed by an arc
% of radius R (counter-clockwise, since the face is on the left of its oriented boundary).
  V = f.V; own = find(any(f.F == k, 2)).';
  seg = zeros(0,4); outEnd = []; inStart = [];
  for e = own
    a = V(f.E(e,1),:);
    if f.E(e,3) == 0
      d = V(f.E(e,2),:) - a; far = a + d/norm(d)*R;
      if f.F(e,1) == k, seg(end+1,:) = [a far]; outEnd = far; else, seg(end+1,:) = [far a]; inStart = far; end
    else
      b = V(f.E(e,2),:);
      if f.F(e,1) == k, seg(end+1,:) = [a b]; else, seg(end+1,:) = [b a]; end
    end
  end
  if ~isempty(outEnd)
    t1 = atan2(outEnd(2), outEnd(1)); t2 = atan2(inStart(2), inStart(1));
    if t2 <= t1, t2 = t2 + 2*pi; end
    t = linspace(t1, t2, 400).'; arc = [cos(t) sin(t)]*norm(outEnd);
    arc(1,:) = outEnd; arc(end,:) = inStart;
    seg = [seg; arc(1:end-1,:) arc(2:end,:)];
  end
  % chain the oriented segments into a loop
  L = seg(1,1:2); cur = seg(1,3:4); used = false(size(seg,1),1); used(1) = true;
  for it = 1:size(seg,1)
    L(end+1,:) = cur; %#ok<AGROW>
    j = find(~used & vecnorm(seg(:,1:2) - cur, 2, 2) < 1e-9*max(1,norm(cur)), 1);
    if isempty(j), break, end
    used(j) = true; cur = seg(j,3:4);
  end
end

function C = circ(R), t = linspace(0, 2*pi, 2000).'; C = R*[cos(t) sin(t)]; end

function v = supPoly(row, L, s)
% max of <s,x> - q(x) over the closed polygon L (vertices in order): every edge in closed form,
% plus the interior critical point when q is strictly convex and the point is inside.
  v = -inf; n = size(L,1);
  for i = 1:n
    v = max(v, supSeg(row, L(i,:), L(mod(i,n)+1,:), s));
  end
  H = [row(5) row(6); row(6) row(7)]; Lc = [row(8); row(9)];
  if H(1,1) > 0 && det(H) > 0
    x = (H \ (s(:) - Lc)).';
    if inpolygon(x(1), x(2), L(:,1), L(:,2)), v = max(v, s*x.' - qv(row, x)); end
  end
end

function v = supSeg(row, a, b, s)
  d = b - a; H = [row(5) row(6); row(6) row(7)]; Lc = [row(8); row(9)];
  A = -0.5*d*H*d.'; B = s*d.' - (a*H + Lc.')*d.'; C = s*a.' - qv(row, a);
  t = [0 1];
  if A < 0, t(end+1) = min(max(-B/(2*A), 0), 1); end
  v = max(A*t.^2 + B*t + C);
end

function y = qv(row, x)
  H = [row(5) row(6); row(6) row(7)];
  y = 0.5*x*H*x.' + x*[row(8); row(9)] + row(10);
end
