function backend_bench()
% How fast are the candidate exact-integer backends, on the operation that actually overflows:
% a 4x4 determinant of quadratic polynomials (conicMeet's resultant).
    N = 200;
    rng(7);
    A = randi([-10000 10000], N, 6);
    B = randi([-10000 10000], N, 6);

    % --- baseline: doubles (WRONG past 2^53, but the speed floor) ------------------------------
    t = tic;
    for i = 1:N, dummy = detLike(A(i,:), B(i,:), @(x)x); end %#ok<NASGU>
    tD = toc(t);
    fprintf('doubles            %8.1f us/op   (the floor; silently wrong past 2^53)\n', 1e6*tD/N);

    % --- java.math.BigInteger: unbounded, NO COMPILER NEEDED ------------------------------------
    t = tic;
    for i = 1:N, dummy = detLikeBig(A(i,:), B(i,:)); end %#ok<NASGU>
    tJ = toc(t);
    fprintf('java BigInteger    %8.1f us/op   (%.0fx doubles) -- ships with MATLAB\n', ...
        1e6*tJ/N, tJ/tD);

    % --- sym: unbounded, needs the Symbolic Toolbox + VPN ---------------------------------------
    n = 20;
    t = tic;
    for i = 1:n, dummy = detLike(sym(A(i,:)), sym(B(i,:)), @(x)x); end %#ok<NASGU>
    tS = toc(t);
    fprintf('sym                %8.1f us/op   (%.0fx doubles) -- licence + VPN on the compute path\n', ...
        1e6*tS/n, (tS/n)/(tD/N));
end

function d = detLike(a, b, ~)
% a stand-in with the same shape as the resultant: products of degree 8 in the inputs
    d = a(1)*b(3)*a(5)*b(6) - a(2)*b(2)*a(4)*b(4) + a(3)*b(1)*a(6)*b(5) - a(4)*b(6)*a(1)*b(3);
    d = d * (a(1)*b(1) - a(2)*b(2));
end

function d = detLikeBig(a, b)
    B = @(x) java.math.BigInteger(sprintf('%d', x));
    t1 = B(a(1)).multiply(B(b(3))).multiply(B(a(5))).multiply(B(b(6)));
    t2 = B(a(2)).multiply(B(b(2))).multiply(B(a(4))).multiply(B(b(4)));
    t3 = B(a(3)).multiply(B(b(1))).multiply(B(a(6))).multiply(B(b(5)));
    t4 = B(a(4)).multiply(B(b(6))).multiply(B(a(1))).multiply(B(b(3)));
    s  = t1.subtract(t2).add(t3).subtract(t4);
    d  = s.multiply(B(a(1)).multiply(B(b(1))).subtract(B(a(2)).multiply(B(b(2)))));
end
