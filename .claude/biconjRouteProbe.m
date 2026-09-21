function biconjRouteProbe()
% Does `biconj` now ANSWER a two-piece input whose pieces meet along a curved tie, and is the answer
% right? The oracle is the assembled envelope's own evaluator, which was verified against the
% definition; `biconj` is supposed to be the same function wearing a type.
    for c = [0 -2]
        f = fixture(c);
        try
            h = biconjQ(f);
            kind = class(h);
        catch ME
            fprintf('c = %2d: REFUSED %s\n', c, ME.identifier);
            continue
        end
        env = envelopePrune(envelopeAssemble(conjQ(tri([0 0; 1 0; 1 1], 0)), ...
                                            conjQ(tri([0 0; 1 1; 0 1], c))));
        rng(37);
        X = [rand(300,1), rand(300,1); rand(100,1)*3-1, rand(100,1)*3-1];
        ref = envelopeEval(env, X);
        got = h.eval(X);
        both = isfinite(ref) & isfinite(got);
        fprintf('c = %2d: %-8s %d faces | domain mismatches %d | max |diff| %.3e over %d points\n', ...
                c, kind, h.nf, sum(isfinite(ref) ~= isfinite(got)), ...
                max(abs(got(both) - ref(both))), sum(both));
    end
end

function f = fixture(c)
% biconjQTest's own fixture: the unit square split by its diagonal, both pieces x^2 + y^2 up to a
% multiple of the diagonal's linear form, so the sign of c alone decides convexity.
    q0 = [0 0 0 0, 2, 0, 2, 0, 0, 0];
    q1 = q0;  q1(8) = q1(8) - c;  q1(9) = q1(9) + c;
    f = QuaPol([0 0; 1 0; 1 1; 0 1], ...
               [1 2 1; 2 3 1; 3 4 1; 4 1 1; 1 3 1], [q0; q1], ...
               [1 0; 1 0; 2 0; 2 0; 2 1]);
end

function p = piece(f, k)
% The k-th piece alone, as the triangle it lives on.
    if k == 1
        p = tri([0 0; 1 0; 1 1], 0);
    else
        p = tri([0 0; 1 1; 0 1], evalin('caller', 'c'));
    end
end

function f = tri(V, c)
    q = [0 0 0 0, 2, 0, 2, 0, 0, 0];  q(8) = q(8) - c;  q(9) = q(9) + c;
    f = QuaPol(V, [1 2 1; 2 3 1; 3 1 1], q, [1 0; 1 0; 1 0]);
end
