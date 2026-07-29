% Direct unit checks of the new LP certificates, independent of the pipeline.
cd('C:\users\ylucet\ai\CCA2');
syms x y
V = [x, y];

% --- maxLinear -------------------------------------------------------------------------
[v,st] = region.maxLinear([1 0; 0 1; -1 -1], [1;1;0], [1 1]);
fprintf('maxLinear bounded   : val=%g st=%d (want 2, 0)\n', v, st);
[v,st] = region.maxLinear([-1 0], [0], [1 0]);
fprintf('maxLinear unbounded : val=%g st=%d (want inf, 1)\n', v, st);
[v,st] = region.maxLinear([1 0; -1 0], [-1; -1], [1 0]);
fprintf('maxLinear infeasible: val=%g st=%d (want -inf, -1)\n', v, st);

% --- linearForm on a mixed region -------------------------------------------------------
r = region([-x, -y, x+y-2, x^2+y^2-9], V);
[A,b,lin] = r.linearForm;
fprintf('linearForm lin flags: %s (want 1 1 1 0)\n', mat2str(lin));
fprintf('linearForm rows     : %s | %s\n', mat2str(A(lin,:)), mat2str(b(lin)'));

% --- redundancy: x<=3 is implied by x<=1 -------------------------------------------------
r = region([-x, -y, x-1, x-3, y-1], V);
fprintf('redundant of [4]    : %s (want 4)\n', mat2str(r.redundantSubset(1:5)));
% nothing is redundant in a plain triangle
r = region([-x, -y, x+y-2], V);
fprintf('redundant of [1 2 3]: %s (want empty)\n', mat2str(r.redundantSubset(1:3)));
% two identical constraints: exactly one may go, never both
r = region([-x, -y, x+y-2, x+y-2], V);
fprintf('dup ineqs n=%d, redundant: %s\n', size(r.ineqs,2), mat2str(r.redundantSubset(1:size(r.ineqs,2))));

% --- unionIsExact -----------------------------------------------------------------------
% Two unit squares sharing the edge x=1: union IS the rectangle [0,2]x[0,1].
A1 = region([-x, -y, x-1, y-1], V);          % [0,1]x[0,1], facet x-1<=0 is ineq 3
A2 = region([1-x, -y, x-2, y-1], V);         % [1,2]x[0,1], facet 1-x<=0
i1 = 0; i2 = 0;
for i=1:size(A1.ineqs,2), for j=1:size(A2.ineqs,2)
    if A1.ineqs(i) == -A2.ineqs(j), i1=i; i2=j; end
end, end
fprintf('squares share facet at (%d,%d); unionIsExact=%d (want 1)\n', i1, i2, A1.unionIsExact(A2,i1,i2));
[l, m] = A1.merge(A2);
fprintf('merge squares: l=%d nIneqs=%d\n', l, size(m.ineqs,2));
if l, fprintf('  contains (1.5,0.5)=%d (want 1), (1.5,1.5)=%d (want 0)\n', ...
        m.ptFeasible(V,[1.5 0.5]), m.ptFeasible(V,[1.5 1.5])); end

% Same shared edge, but the second box is TALLER: union is an L, not convex. Must refuse,
% and the old recipe would have returned [0,2]x[0,2] -- claiming (0.5,1.5), in neither.
B2 = region([1-x, -y, x-2, y-2], V);
i1 = 0; i2 = 0;
for i=1:size(A1.ineqs,2), for j=1:size(B2.ineqs,2)
    if A1.ineqs(i) == -B2.ineqs(j), i1=i; i2=j; end
end, end
fprintf('L-shape share facet at (%d,%d); unionIsExact=%d (want 0)\n', i1, i2, A1.unionIsExact(B2,i1,i2));
[l, m] = A1.merge(B2);
fprintf('merge L-shape: l=%d (want 0)\n', l);
