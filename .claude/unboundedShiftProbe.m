function unboundedShiftProbe()
% The three expectations that changed when the unbounded envelope landed. For each: what does
% biconjQ do now, and is that the right thing? An expectation is only edited once this says so.
    fprintf('=== 1. the "fin" fixture: cone((1,0),(1,1)), H = diag(1,-1) ===\n');
    fin = QuaPol([0 0; 1 0; 1 1], [1 2 0; 1 3 0], ...
                 [0 0 0 0, 1, 0, -1, 0, 0, 0], [1 0; 0 1]);
    report(fin);
    fprintf(['  the ray (1,0) has d''Hd = %d. Positive means q grows QUADRATICALLY along it, so\n' ...
             '  the piece is not edge-concave there and the lifted hull is the wrong formula --\n' ...
             '  which is an edge-convexity obstruction, not an unboundedness one.\n'], 1);

    fprintf('\n=== 2. q = x*y on the first quadrant ===\n');
    q1 = QuaPol([0 0; 1 0; 0 1], [1 2 0; 1 3 0], ...
                [0 0 0 0, 0, 1, 0, 0, 0, 0], [1 0; 0 1]);
    report(q1);
    fprintf(['  H*(1,0)'' = (0,1) is nonzero while (1,0)''H(1,0) = 0, so the rate along that ray\n' ...
             '  varies across the piece: exactly the sheared case, and co f = 0 here is real but\n' ...
             '  not reachable by a fixed lifted ray.\n']);

    fprintf('\n=== 3. examples{16}, the half-star (five rays, f = 0) ===\n');
    P = QuaPol.examples();  hs = P{16};
    report(hs);
    try
        h = biconjQ(hs);
        X = [-1 -1; 1 0; 0 -0.5; 2 -3; -4 1; 0 0; 5 5];
        fprintf('  x            f          co f\n');
        for i = 1:size(X,1)
            fprintf('  (%5.2g,%5.2g) %8.4g %10.4g\n', X(i,1), X(i,2), hs.eval(X(i,:)), h.eval(X(i,:)));
        end
        fprintf(['\n  f is 0 on its domain and +inf off it. The five rays span from 0 to 225\n' ...
                 '  degrees -- MORE than a half turn -- so conv(dom f) is the whole plane and\n' ...
                 '  co f = 0 everywhere is the correct envelope, including at (0,-0.5), the\n' ...
                 '  midpoint witness that is outside dom f.\n']);
        bad = max(abs(h.eval(X) - 0));
        fprintf('  max |co f - 0| over those points: %.3e\n', bad);
    catch e
        fprintf('  biconjQ raises %s\n', e.identifier);
    end
end

function report(f)
    try
        h = biconjQ(f);
        fprintf('  biconjQ ANSWERS, nf = %d\n', h.nf);
    catch e
        fprintf('  biconjQ raises %s\n', e.identifier);
    end
end
