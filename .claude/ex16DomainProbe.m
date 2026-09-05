function ex16DomainProbe()
% Is the domain of QuaPol.examples(){16} convex? isConvexQ says no and the legacy dom.isConvex
% says yes, so one of them is wrong. Decided WITHOUT either of them: pick two points the mesh
% itself says are in the domain and ask whether their midpoint is.
%
% The domain is a union of cones at the origin, so membership is tested through eval, which
% returns +Inf outside the domain.
    P = QuaPol.examples();  f = P{16};
    fprintf('V =\n'); disp(f.V);
    fprintf('E =\n'); disp(f.E);
    fprintf('F =\n'); disp(f.F);
    fprintf('dom.isConvex (legacy, doubles) = %d\n', f.dom.isConvex);

    X = [-1 -1; 1 0; 0 -0.5; -0.5 -0.5; 0 1; -1 0; 1 1; 0.5 -0.5];
    v = f.eval(X);
    fprintf('\n%-14s %-10s %s\n','point','f value','in domain');
    for i = 1:size(X,1)
        fprintf('(%5.2f,%5.2f)  %-10.4g %d\n', X(i,1), X(i,2), v(i), isfinite(v(i)));
    end
    fprintf(['\nIf (-1,-1) and (1,0) are in the domain and their midpoint (0,-0.5) is not, the ' ...
             'domain is NOT convex.\n']);
end
