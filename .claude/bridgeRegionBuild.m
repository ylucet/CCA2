function bridgeRegionBuild()
% ITEM 2: THE REGION A BRIDGE FACE OWNS, constructed EXACTLY rather than traced.
%
% `.claude/bridgeRegionProbe.m` established WHAT the region is, by classification and curve fitting:
% its boundaries are conics, line fits missing by 0.30 to 0.71 where conic fits sit at 1e-14. This
% builds them, in integers, from the cells' own data.
%
% THE CONSTRUCTION, in one line. At a tie point the subdifferential of a max is the segment between
% the two gradients, so
%
%       bridge region  =  union over s on the tie arc of  conv{ grad g1(s), grad g2(s) } ,
%
% and each `grad g_k(s) = A_k s + b_k` is AFFINE. So each boundary curve is the IMAGE of the tie
% conic under an affine map -- and the image of a conic under an invertible affine map is a conic
% whose coefficients are an exact substitution away:
%
%       image of { q = 0 } under x = A s + b   is   { q(A^-1 (x - b)) = 0 } ,
%
% which is `ratQ.substAffine` with `M = adj(A)`, `t = -adj(A) b` and map denominator `det(A)`,
% both numerator and denominator integral because `A` is. No fitting, no sampling, no tolerance.
%
% WHY A CONIC GOES THROUGH A ROUTINE THAT TAKES A FACE. `substAffine` substitutes into RatCon's
% 10-wide weighted cubic basis, where the value is `c5 x^2/2 + c6 xy + c7 y^2/2 + c8 x + c9 y +
% c10`. A conic `[a b c d e f]` is the same FUNCTION spelled `c5 = 2a, c6 = b, c7 = 2c` and the
% rest unchanged, and a function's zero set is its zero set, so the conversion is exact and the
% substitution is the one already tested.
%
% THE DEGENERATE CASE IS THE EASY ONE, and it is not rare. Over a VERTEX cell a conjugate is
% AFFINE, so `A_k = 0` and `grad g_k` is the CONSTANT `b_k`: the whole boundary curve collapses to
% a single POINT and the region is a fan from it over the other boundary. `substAffine` refuses a
% singular map by name -- correctly, since it is not a change of variables -- so that case is
% detected before it is called rather than by catching its error.
%
% THE TEST IS DIFFERENTIAL AND EXACT IN INTENT. Take points on the tie arc, map them through each
% gradient, and evaluate the constructed conic there: it must vanish. That checks the construction
% against the definition it came from and shares no arithmetic with it.

    fprintf('\n=== isotropic: g1 = |s|^2/2, g2 = <s,a> + |s|^2/4 - c ===\n');
    a = [2; 1];  c = -1;
    check(struct('A', 2*eye(2),   'b', [0;0], 'k',  0), ...
          struct('A', 0.5*eye(2), 'b', a,     'k', -c));

    fprintf('\n=== anisotropic: g2 with a genuinely different Hessian ===\n');
    check(struct('A', 2*eye(2), 'b', [0;0], 'k', 0), ...
          struct('A', [1 0.25; 0.25 1.5], 'b', a, 'k', -c));

    fprintf('\n=== VERTEX cell: g2 AFFINE, so its boundary collapses to a point ===\n');
    check(struct('A', 2*eye(2), 'b', [0;0], 'k', 0), ...
          struct('A', zeros(2), 'b', a, 'k', -c));
end

% ------------------------------------------------------------------------------------------------
function check(g1, g2)
% Build both boundaries, then verify each against tie points mapped through its own gradient.
    q = diffConic(g1, g2);
    fprintf('  tie conic q = g1 - g2 : [%s]   kind %s\n', num2str(q), ratQ.conicKind(ratQ.conic(q)));

    S = tiePoints(q, 400);
    if size(S,1) < 20
        fprintf('  the tie has too few real points here; nothing to verify.\n');  return
    end
    fprintf('  %d tie points found\n', size(S,1));

    for k = 1:2
        if k == 1, gk = g1; else, gk = g2; end
        [cimg, kind] = imageConic(q, gk);
        X = S * gk.A.' + gk.b.';                        % grad g_k(s) = A_k s + b_k, row-wise
        if isempty(cimg)
            % A_k singular: the image is not a conic. Report what it actually is.
            spread = max(max(X,[],1) - min(X,[],1));
            fprintf('  boundary %d: %s (image spread %.3e)\n', k, kind, spread);
            continue
        end
        v = evalConic(cimg, X);
        scale = max(1, max(abs(evalConic(abs(cimg), abs(X)))));
        fprintf('  boundary %d: conic [%s] kind %s | max |q_img| on the mapped tie %.3e (rel %.3e)\n', ...
                k, num2str(cimg), ratQ.conicKind(ratQ.conic(cimg)), max(abs(v)), max(abs(v))/scale);
    end
end

function [cimg, kind] = imageConic(q, gk)
% The conic through grad g_k(tie). Empty when A_k is singular, with `kind` saying what it is.
    A = gk.A;  b = gk.b;
    D = det(A);
    if abs(D) < 1e-12
        if all(abs(A(:)) < 1e-12), kind = 'a single POINT (the cell is a VERTEX cell, g_k affine)';
        else,                      kind = 'a LINE (A_k has rank 1: an EDGE cell)';
        end
        cimg = [];  return
    end
    % Clear A and b to integers so substAffine's exact path applies. The map is x = A s + b, so the
    % substitution needed is its INVERSE, s = adj(A)(x - b)/det(A).
    [An, ad] = intify(A);  [bn, bd] = intify(b);
    L = lcm(ad, bd);  An = An * (L/ad);  bn = bn * (L/bd);   % A = An/L, b = bn/L
    adjA = [An(2,2), -An(1,2); -An(2,1), An(1,1)];           % adj(An), so adj(A) = adjA/L
    dAn  = ratQ.detExact(An);                                % det(A) = dAn/L^2
    % s = adj(A)(x-b)/det(A) = (adjA/L)(x - bn/L) * L^2/dAn = (L*adjA*x - adjA*bn)/dAn
    Mn = ratQ.chk(L * adjA,      'image map');
    tn = ratQ.chk(-adjA * bn,    'image translation');
    [n2, d2] = ratQ.substAffine(faceOf(q), 1, Mn, tn, dAn);
    cimg = ratQ.conic(conicOf(n2));                          % d2 scales the whole form; the zero
    kind = '';                                               % set is unchanged, so it is dropped
end

% ---- conic <-> the 10-wide weighted face basis -------------------------------------------------
function n = faceOf(q),  n = [0 0 0 0, 2*q(1), q(2), 2*q(3), q(4), q(5), q(6)]; end
function q = conicOf(n)
% TWICE the face, so the result is integral: the weighted basis carries halves on the squared
% terms, and a conic has NO SCALE -- { v = 0 } and { 2v = 0 } are one curve -- so doubling is free
% here and dividing by two is not. Halving instead produced a non-integer conic and ratQ.conic
% refused it by name, which is the guard working.
    q = [n(5), 2*n(6), n(7), 2*n(8), 2*n(9), 2*n(10)];
end

function q = diffConic(g1, g2)
% g_k(s) = s'A_k s/2 + b_k's + k_k, so g1 - g2 is the conic below, cleared to integers.
    A = g1.A - g2.A;  b = g1.b - g2.b;  k = g1.k - g2.k;
    q = [A(1,1)/2, A(1,2), A(2,2)/2, b(1), b(2), k];
    [q, ~] = intify(q);
end

function [N, d] = intify(X)
% Clear a small rational array to integers by the lcm of the denominators its entries need.
    d = 1;
    for v = X(:).'
        [~, dv] = rat(v, 1e-12);
        d = lcm(d, abs(dv));
    end
    N = round(X * d);
end

function S = tiePoints(q, n)
% Real points of { q = 0 }, found by solving the conic in y for each x on a grid. Solving beats
% bisecting: a ray-based search misses a conic whose branch the ray does not cross, which is the
% oracle defect DECISIONS 2026-09-06 records costing two wrong conclusions.
    S = zeros(0,2);
    for x = linspace(-12, 12, n)
        A = q(3);  B = q(2)*x + q(5);  C = q(1)*x^2 + q(4)*x + q(6);
        if abs(A) < 1e-14
            if abs(B) > 1e-14, S(end+1,:) = [x, -C/B]; end %#ok<AGROW>
            continue
        end
        disc = B^2 - 4*A*C;
        if disc < 0, continue, end
        r = sqrt(disc);
        S(end+1,:) = [x, (-B + r)/(2*A)]; %#ok<AGROW>
        if r > 0, S(end+1,:) = [x, (-B - r)/(2*A)]; end %#ok<AGROW>
    end
    S = S(all(abs(S) < 1e6, 2), :);
end

function v = evalConic(c, X)
    v = c(1)*X(:,1).^2 + c(2)*X(:,1).*X(:,2) + c(3)*X(:,2).^2 + ...
        c(4)*X(:,1)      + c(5)*X(:,2)        + c(6);
end
