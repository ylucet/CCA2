function big_fair()
% My first BigInteger benchmark built each operand with sprintf -- string parsing per number, which
% is not how it would be used. valueOf(long) is the fair comparison.
    rng(7); N = 500;
    A = randi([-10000 10000], N, 6);  B = randi([-10000 10000], N, 6);

    t = tic;
    for i = 1:N, d = viaString(A(i,:), B(i,:)); end %#ok<NASGU>
    fprintf('BigInteger via sprintf   %8.1f us/op   (what I benchmarked first)\n', 1e6*toc(t)/N);

    t = tic;
    for i = 1:N, d = viaValueOf(A(i,:), B(i,:)); end %#ok<NASGU>
    fprintf('BigInteger via valueOf   %8.1f us/op   (the fair comparison)\n', 1e6*toc(t)/N);

    % and confirm it is exact
    want = detSym(sym(A(1,:)), sym(B(1,:)));
    got  = viaValueOf(A(1,:), B(1,:));
    fprintf('exact against sym: %d   (value %s)\n', sym(char(got.toString())) == want, char(want));
end

function d = viaString(a, b)
    B = @(x) java.math.BigInteger(sprintf('%d', x));
    d = combine(B, a, b);
end
function d = viaValueOf(a, b)
    B = @(x) java.math.BigInteger.valueOf(int64(x));
    d = combine(B, a, b);
end
function d = combine(B, a, b)
    t1 = B(a(1)).multiply(B(b(3))).multiply(B(a(5))).multiply(B(b(6)));
    t2 = B(a(2)).multiply(B(b(2))).multiply(B(a(4))).multiply(B(b(4)));
    t3 = B(a(3)).multiply(B(b(1))).multiply(B(a(6))).multiply(B(b(5)));
    t4 = B(a(4)).multiply(B(b(6))).multiply(B(a(1))).multiply(B(b(3)));
    s  = t1.subtract(t2).add(t3).subtract(t4);
    d  = s.multiply(B(a(1)).multiply(B(b(1))).subtract(B(a(2)).multiply(B(b(2)))));
end
function d = detSym(a,b)
    d = a(1)*b(3)*a(5)*b(6) - a(2)*b(2)*a(4)*b(4) + a(3)*b(1)*a(6)*b(5) - a(4)*b(6)*a(1)*b(3);
    d = d * (a(1)*b(1) - a(2)*b(2));
end
