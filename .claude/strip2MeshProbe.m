function strip2MeshProbe()
% Is the sweep's `strip2` fixture the function it was meant to be, or is its split edge oriented
% backwards? Asked BEFORE looking at biconjQ, because a malformed mesh would make the routine's
% answer meaningless rather than wrong.
%
% The intended function: q0 on the unit square {0<=x<=1, 0<=y<=1} and q0 + c*(y-1) on the strip
% above it {0<=x<=1, y>=1}. Face 2 is the upper one, so on edge 4 -- which runs (0,1) -> (1,1),
% direction (1,0), whose LEFT is +y -- face 2 must sit in column 1.
    c  = -2;
    q0 = [0 0 0 0, 2, 0, 0, 0, 0, 0];                 % H = diag(2,0): q0 = x^2
    q1 = q0;  q1(9) = q1(9) + c;  q1(10) = q1(10) - c; % + c*(y-1)

    V = [0 0; 1 0; 0 1; 1 1; 0 2; 1 2];
    E = [1 2 1; 1 3 1; 2 4 1; 3 4 1; 3 5 0; 4 6 0];
    Fbad  = [1 0; 0 1; 1 0; 1 2; 0 2; 2 0];           % what the sweep had
    Fgood = [1 0; 0 1; 1 0; 2 1; 0 2; 2 0];           % edge 4 flipped

    X = [0.5 0.5; 0.5 3; 0.25 1; 0.5 1];              % lower, upper, and two ON the split
    want = zeros(size(X,1),1);
    for i = 1:size(X,1)
        x = X(i,1);  y = X(i,2);
        if y <= 1, want(i) = x^2; else, want(i) = x^2 + c*(y-1); end
    end

    tags = {'as written (F(4,:) = [1 2])', 'flipped    (F(4,:) = [2 1])'};
    Fs   = {Fbad, Fgood};
    for t = 1:2
        f = QuaPol(V, E, [q0; q1], Fs{t});
        got = f.eval(X);
        fprintf('\n%s\n', tags{t});
        fprintf('   %-14s %10s %10s\n', 'point', 'intended', 'f.eval');
        for i = 1:size(X,1)
            fprintf('   (%4.2f,%5.2f) %10.4g %10.4g\n', X(i,1), X(i,2), want(i), got(i));
        end
        fprintf('   max |f.eval - intended| = %.3g\n', max(abs(got - want)));
    end
    fprintf(['\nThe orientation rule: edge 4 goes (0,1) -> (1,1), so its direction is (1,0) and its\n' ...
             'LEFT is +y, which is the UPPER face. F(j,1) is the left face, so it must be 2.\n']);
end
