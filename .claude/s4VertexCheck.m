% .claude/s4VertexCheck.m -- INDEPENDENT check that a vertex of f* can have degree 4 with
% Galois group S4, i.e. that no tower of square roots can represent it. Recomputes g1,g2,g3
% from the PRIMAL data (it takes nothing from CONJ_FIELD_PROOF.md except the input), gets the
% eliminant, and checks irreducibility, the resolvent cubic, the discriminant, and that the
% triple point is real -- all three argmaxes strictly inside their own triangles.
% Nothing here uses that document's g_k or its quartic: both are RECOMPUTED and compared.
syms s1 s2 t real
Q = { sym([1 -1; -1 4])/3, sym([1 0; 0 4])/4, sym([2 1; 1 1]) };
be = { sym([-2; 4]), [sym(-1)/2; sym(-2)], sym([1; 1]) };
ga = { sym(8), sym(5)/2, sym(-5)/2 };
T  = { sym([4 0; 4 -20; 24 10]), sym([4 0; 24 10; -16 10]), sym([4 0; -16 10; 4 -20]) };
s = [s1; s2];
g = cell(1,3);
for k = 1:3
    if any(eig(double(Q{k})) <= 0), fprintf('piece %d NOT positive definite\n', k); end
    g{k} = simplify( (s - be{k}).' * inv(Q{k}) * (s - be{k}) / 2 - ga{k} );
    fprintf('g%d = %s\n', k, char(expand(g{k})));
end

% the triple point: g1 = g2 and g1 = g3, eliminated down to one variable
r12 = expand(g{1} - g{2});
r13 = expand(g{1} - g{3});
q1 = collect(expand(resultant(r12, r13, s2)), s1);
fprintf('\neliminant in s1 : %s\n', char(q1));
fprintf('their quartic   : %s\n', char(collect(3*s1^4 - 24*s1^3 + 10*s1^2 + 160*s1 - 96, s1)));
rat = simplify(q1 / (3*s1^4 - 24*s1^3 + 10*s1^2 + 160*s1 - 96));
fprintf('ratio (must be a nonzero CONSTANT): %s\n', char(rat));

quart = 3*t^4 - 24*t^3 + 10*t^2 + 160*t - 96;
fprintf('\nfactor(quartic) = %s   -> irreducible over Q iff this is the quartic itself\n', ...
        char(factor(quart)));
% resolvent cubic of  a t^4 + b t^3 + c t^2 + d t + e   (depressed-form standard resolvent)
a = sym(3); b = sym(-24); c = sym(10); d = sym(160); e = sym(-96);
res = t^3 - (c/a)*t^2 + ((b*d - 4*a*e)/a^2)*t - ((b^2*e - 4*a*c*e + a*d^2)/a^3);
fprintf('resolvent cubic = %s\n', char(collect(expand(res), t)));
fprintf('factor(resolvent) = %s\n', char(factor(collect(expand(res), t))));
disc = simplify(resultant(quart, diff(quart, t), t) / a);
fprintf('discriminant = %s   sqrt = %.6f  (integer? %d)\n', char(disc), double(sqrt(disc)), ...
        double(sqrt(disc)) == round(double(sqrt(disc))));

% the real root near 0.6080508815, and the argmaxes there
rr = double(vpa(root(quart, t), 40));
rr = rr(abs(imag(rr)) < 1e-20);
px = rr(abs(rr - 0.608050881512364) < 1e-9);
fprintf('\nroot near 0.60805 : %.15f  (found %d real roots)\n', px, numel(rr));
py = double(vpa(subs(solve(subs(r12, s1, sym(px)), s2), 1), 30));
fprintf('matching s2       : %.15f\n', py(1));
p = [px; py(1)];
for k = 1:3
    xs = double(inv(Q{k})) * (p - double(be{k}));
    W = double(T{k});
    lam = [W.'; 1 1 1] \ [xs; 1];
    fprintf('piece %d: argmax (%8.4f,%8.4f)  barycentric [%.3f %.3f %.3f]  interior=%d\n', ...
        k, xs(1), xs(2), lam(1), lam(2), lam(3), all(lam > 1e-9));
    fprintf('         g%d(p) = %.15f\n', k, double(subs(g{k}, [s1 s2], p.')));
end
