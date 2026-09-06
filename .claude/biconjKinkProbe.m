function biconjKinkProbe()
% WHICH refusal does a BOUNDED, continuous, two-piece f with CONVEX pieces and a non-convex kink
% get? biconjAllHSweep shows poly2/PD answering only 2 of 6 and square2/PD 4 of 6, and the
% distinguishing knob is the sign of c -- so this is the shape of gap the sweep found. Naming the
% identifier says whether it is a known gap or a new one.
    for shape = ["square2", "poly2"]
        for c = [-2 0 2]
            f = fixture(shape, c);
            conv = 'ERR';
            try
                [t, ~] = isConvexQ(f);
                if t, conv = 'convex'; else, conv = 'NOT convex'; end
            catch
            end
            try
                h = biconjQ(f);
                res = sprintf('OK nf=%d', h.nf);
            catch e
                res = e.identifier;
            end
            fprintf('%-9s c=%+d  isConvexQ: %-11s  biconjQ: %s\n', shape, c, conv, res);
        end
    end
end

function f = fixture(shape, c)
    H = [2 0; 0 2];  L = [0; 0];  k0 = 0;
    q0 = [0 0 0 0, H(1,1), H(1,2), H(2,2), L(1), L(2), k0];
    if shape == "square2"
        V = [0 0; 1 0; 1 1; 0 1];
        E = [1 2 1; 2 3 1; 3 4 1; 4 1 1; 1 3 1];
        F = [1 0; 1 0; 2 0; 2 0; 2 1];
        lin = [-1 1 0];
    else
        V = [-2 0; 2 0; 1 2; -1 2];
        E = [1 2 1; 2 3 1; 3 4 1; 4 1 1; 1 3 1];
        F = [1 0; 1 0; 2 0; 2 0; 2 1];
        lin = [0 1 0];
    end
    q1 = q0;  q1(8) = q1(8) + c*lin(1);  q1(9) = q1(9) + c*lin(2);  q1(10) = q1(10) + c*lin(3);
    f = QuaPol(V, E, [q0; q1], F);
end
