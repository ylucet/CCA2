function edgeConvexIndexProbe()
% Is `k = sum(obj.F(j,:)>0)` really used as a FACE INDEX in isEdgeConvex, and does it give a wrong
% answer? Rung 1 of the ladder: observe, do not edit.
%
% The mesh is the closed upper half-plane split by the ray x=0, y>=0. Edge 2 is the negative
% x-axis, whose only finite face is face 1; edge 3 is the positive x-axis, whose only finite face
% is face 2. Both are boundary edges, so both take the k==1 branch, and both then ask
% isFaceConvex(1) -- accidentally right for edge 2 and wrong for edge 3.
%
% Face 1 is x^2+y^2 (convex). Face 2 is x^2-y^2 (indefinite). So the correct answers are
% isEdgeConvex(2) = true and isEdgeConvex(3) = FALSE.
    V = [0 0; 0 1; -1 0; 1 0];
    E = [1 2 0; 1 3 0; 1 4 0];
    F = [1 2; 0 1; 2 0];
    f = [2 0 2 0 0 0;      % face 1: x^2 + y^2, Hessian [2 0;0 2], PSD
         2 0 -2 0 0 0];    % face 2: x^2 - y^2, Hessian [2 0;0 -2], indefinite
    p = QuaPol(V, E, f, F);

    fprintf('isFaceConvex(1) = %d  (expected 1)\n', p.isFaceConvex(1));
    fprintf('isFaceConvex(2) = %d  (expected 0)\n', p.isFaceConvex(2));
    for j = 2:3
        [c, cont] = p.isEdgeConvex(j);
        fprintf('isEdgeConvex(%d) = %d, continuous=%d   F(%d,:) = [%d %d]\n', ...
                j, c, cont, j, p.F(j,1), p.F(j,2));
    end
    fprintf('\nedge 3 is the one to watch: its finite face is 2, which is NOT convex.\n');
end
