function crt_exact()
% Is the CRT-reconstructed value EXACT, or only exact-looking? Compare against sym as an INTEGER,
% not as a double -- the previous check compared two doubles that were both already rounded.
    P = [67108859 67108837 67108819 67108777];
    rng(7);
    A = randi([-10000 10000], 1, 6);  B = randi([-10000 10000], 1, 6);
    want = detSym(sym(A), sym(B));                 % the true integer, exactly

    r = zeros(1,numel(P));
    for k = 1:numel(P)
        r(k) = double(mod(want, sym(P(k))));       % the residues, however they were obtained
    end

    % accumulate the reconstruction in a DOUBLE, as a pure-MATLAB route must
    d = garner(r, P);
    x = 0; mul = 1;
    for k = 1:numel(d), x = x + d(k)*mul; mul = mul*P(k); end

    fprintf('true integer   : %s\n', char(want));
    fprintf('double recon   : %.20g\n', x);
    fprintf('EXACT as integers? %d\n', sym(x) == want);
    fprintf('absolute error : %s\n', char(abs(sym(x) - want)));
    fprintf('2^53 = %.0f;  |true| = %s  (%.3g x bigger)\n', 2^53, char(abs(want)), ...
        double(abs(want))/2^53);
    fprintf('\nSO: the residues are exact and cheap; the RECONSTRUCTED VALUE cannot live in a double.\n');
end

function d = garner(r, P)
    n = numel(P);  d = zeros(1,n);  d(1) = r(1);
    for k = 2:n
        p = P(k);  acc = d(1);  mul = 1;
        for j = 2:k-1
            mul = mod(mul*P(j-1), p);
            acc = mod(acc + d(j)*mul, p);
        end
        mul = mod(mul*P(k-1), p);
        if k == 2, mul = mod(P(1), p); end
        d(k) = mod((r(k) - acc) * modinv(mul, p), p);
    end
end
function y = modinv(a,p), [~,x,~] = egcd(mod(a,p),p); y = mod(x,p); end
function [g,x,y] = egcd(a,b)
    if b == 0, g=a; x=1; y=0; return, end
    [g,x1,y1] = egcd(b, mod(a,b));  x = y1;  y = x1 - floor(a/b)*y1;
end
function d = detSym(a,b)
    d = a(1)*b(3)*a(5)*b(6) - a(2)*b(2)*a(4)*b(4) + a(3)*b(1)*a(6)*b(5) - a(4)*b(6)*a(1)*b(3);
    d = d * (a(1)*b(1) - a(2)*b(2));
end
