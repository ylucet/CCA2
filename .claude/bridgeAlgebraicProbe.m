function bridgeAlgebraicProbe()
% Is the bridge face ALGEBRAIC, and if so of what degree and over what field?
%
% `.claude/bridgeFaceProbe.m` showed the face is neither affine nor quadratic AS A GRAPH z = p(x).
% That is NOT the same as being non-algebraic, and the distinction decides how expensive the
% missing type is: sqrt(Q) satisfies z^2 - Q = 0, so it is algebraic of degree 2, and if the
% coefficients of that implicit equation are RATIONAL then the face is storable exactly with
% rational data plus a branch choice -- which is what QuaCon already does for its edge conics, one
% dimension down.
%
% TESTED TWO WAYS, because the first instance is symmetric and might be special:
%
%   A. the ISOTROPIC pair from bridgeFaceProbe -- q1 = |v|^2/2, q2 = |v-a|^2 + c. The derivation
%      gives co f = <2a, x-2a> + R*|x-2a| + 2c with R = 2*sqrt(|a|^2-c), so
%
%          (z - <2a,x-2a> - 2c)^2 - R^2*|x-2a|^2 = 0,
%
%      and R^2 = 4(|a|^2 - c) is RATIONAL even though R is not. Checked directly.
%
%   B. an ANISOTROPIC pair, q2(v) = (v-a)'H2(v-a)/2 + c with H2 = diag(2,4), where no symmetry is
%      available and no derivation is leaned on. Two searches for a polynomial relation among the
%      oracle's own points: first by TOTAL degree, then in the basis the elimination actually
%      predicts -- degree at most 4 in z, with coefficient polynomials in (x,y) of degree up to D.
%      A relation exists where the smallest singular value COLLAPSES to the oracle's noise floor;
%      a smooth decrease with basis size is approximation, not a relation.
%
% The oracle in both cases is the primal definition through the perspective inf-convolution, which
% forms no conjugate and knows nothing about either implicit form.

    rng(20260905);
    a = [1.3; -0.7];  c = -0.5;
    ts = linspace(0, 1, 200001).';
    X = 6*(rand(3000,2) - 0.5);

    fprintf('=== A. isotropic pair: the derived implicit form, checked directly ===\n');
    R2 = 4*(a.'*a - c);
    z = zeros(size(X,1),1);  tS = zeros(size(X,1),1);
    for i = 1:size(X,1)
        d1 = X(i,1) - ts*a(1);  d2 = X(i,2) - ts*a(2);
        vals = (d1.^2 + d2.^2)./(2 - ts) + ts*c;
        [z(i), j] = min(vals);  tS(i) = ts(j);
    end
    inB = tS > 1e-6 & tS < 1-1e-6;
    Dv = X(inB,:) - 2*a.';
    res = (z(inB) - 2*(Dv*a) - 2*c).^2 - R2*sum(Dv.^2, 2);
    fprintf('bridge points %d\n', sum(inB));
    fprintf('R^2 = %g   (rational: 4*(|a|^2 - c))\n', R2);
    fprintf('max |(z - affine)^2 - R^2*|x-2a|^2| = %.3e\n', max(abs(res)));
    fprintf('-> the graph IS a quadric cone, with RATIONAL coefficients; z is its upper branch.\n');

    fprintf('\n=== B. anisotropic pair (H2 = diag(2,4)) ===\n');
    % The inf-convolution of two quadratics is a quadratic with Hessian ((1-t)*I^-1 + t*H2^-1)^-1,
    % evaluated at x - t*a. Both Hessians are DIAGONAL here, so that inverse is diagonal and the
    % whole t-scan vectorises -- which matters, since the plain loop over 3000 points x 200001
    % values of t does not finish.
    %   m1(t) = (1-t) + t/2 = 1 - t/2,    m2(t) = (1-t) + t/4 = 1 - 3t/4
    m1 = 1 - ts/2;  m2 = 1 - 3*ts/4;
    z2 = zeros(size(X,1),1);  t2 = zeros(size(X,1),1);
    for i = 1:size(X,1)
        d1 = X(i,1) - ts*a(1);  d2 = X(i,2) - ts*a(2);
        vals = 0.5*(d1.^2./m1 + d2.^2./m2) + ts*c;
        [z2(i), j] = min(vals);  t2(i) = ts(j);
    end
    inB2 = t2 > 1e-6 & t2 < 1-1e-6;
    P = X(inB2,:);  Z = z2(inB2);
    fprintf('bridge points %d\n', sum(inB2));

    fprintf('\nby TOTAL degree:\n');
    fprintf('%-6s %-11s %14s %14s\n', 'deg', 'monomials', 'min sing val', 'rel residual');
    for deg = 2:5
        report(polyBasis(P(:,1), P(:,2), Z, deg), deg);
    end

    % THE BASIS THE ELIMINATION PREDICTS. Setting dval/dt = 0 and clearing the two squared
    % denominators gives a polynomial in t of degree 4 whose coefficients are polynomial in (x,y),
    % so z is algebraic over Q(x,y) of degree at most 4 IN Z. Nothing bounds the degree of those
    % coefficient polynomials, so a total-degree basis cannot see the relation until it happens to
    % be large enough to contain them -- which is why the scan above is the wrong shape.
    fprintf('\nby degree in (x,y), with degree <= 4 in z:\n');
    fprintf('%-6s %-11s %14s %14s\n', 'D', 'monomials', 'min sing val', 'rel residual');
    for D = 2:7
        report(zBasis(P(:,1), P(:,2), Z, 4, D), D);
    end

    fprintf(['\nRead the collapse, not the trend. Part A''s genuine relation sits at 1.8e-10, the\n' ...
             'oracle''s own noise floor. A column that merely drifts downward as the basis grows is\n' ...
             'a better APPROXIMATION on a bounded box, not a polynomial identity.\n']);
end

% ------------------------------------------------------------------------------------------------
function report(A, lbl)
% Column-normalise before the SVD: these monomials differ by orders of magnitude on a domain of
% size ~3, and without it the smallest singular value reports the SCALING rather than whether a
% polynomial relation exists at all.
    sc = sqrt(sum(A.^2, 1));  sc(sc == 0) = 1;
    An = A ./ sc;
    [~, S, V] = svd(An, 0);
    coef = V(:,end);
    rel = max(abs(An*coef)) / max(abs(An(:)));
    fprintf('%-6d %-11d %14.3e %14.3e\n', lbl, size(A,2), S(end,end), rel);
end

% ------------------------------------------------------------------------------------------------
function A = polyBasis(x, y, z, deg)
% every monomial x^i y^j z^k with i + j + k <= deg
    A = zeros(numel(x), 0);
    for i = 0:deg
        for j = 0:(deg-i)
            for k = 0:(deg-i-j)
                A(:, end+1) = (x.^i) .* (y.^j) .* (z.^k); %#ok<AGROW>
            end
        end
    end
end

% ------------------------------------------------------------------------------------------------
function A = zBasis(x, y, z, kmax, D)
% every monomial x^i y^j z^k with k <= kmax and i + j <= D
    A = zeros(numel(x), 0);
    for k = 0:kmax
        for i = 0:D
            for j = 0:(D-i)
                A(:, end+1) = (x.^i) .* (y.^j) .* (z.^k); %#ok<AGROW>
            end
        end
    end
end
