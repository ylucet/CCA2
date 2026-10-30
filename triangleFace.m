function tf = triangleFace(tri)
% triangleFace  A triangle face of the biconjugate, held EXACTLY by the name of its dual point --
%   step (b) of Example 2's subdivision (TODO item 1, DECISIONS.md 2026-10-30 (fifth)).
%
% objective: from one entry of biconjAtoms(f, 'triangles'), the face z(x) = <s0,x> - g(s0) on the
%   contact triangle conv{m_1(s0), m_2(s0), m_3(s0)}, with s0 NAMED (never written down) and every
%   membership question decided exactly by `triangleFaceSide`.
%
% [input]  tri : struct with .Gn (3x6 integer conics) and .Gd (3x1 positive): the three forms g_i
%                exactly (biconjAtoms' formQ), and .s (the double s0 the forms tie at)
% [output] tf  : struct
%                  .kind    'rational' -- s0 is rational (both ties are LINES: three vertex forms,
%                           or any pair of equal Hessians); 'algebraic' -- s0 needs a name
%                  .s0n .s0d  s0 = s0n/s0d exactly, for 'rational'
%                  .Tn .Td    the corners Tn(i,:)/Td(i) exactly, for 'rational'
%                  .dp      dualPoint's name of s0, for 'algebraic'
%                  .Gn .Gd  the forms, as given
%                  .s .c .T the doubles: s0, g(s0) and the 3x2 contact triangle -- the CACHE
%                  .orient  the sign of the triangle's area, +1 or -1, decided exactly
%
% THE CORNERS ARE AFFINE IN s0. Each form's contact point is its gradient, m_i(s) = grad g_i(s), a
% linear map plus a constant (biconjAtoms' formQ); so each side of the triangle, at a rational x,
% is a polynomial of degree 2 in s0 with rational coefficients -- dualPointSign's question. The
% value is the same: affine in x over Q(s0), read in doubles from the cache.
%
% WHICH TWO TIES NAME s0. g1 = g2 and g1 = g3; the third tie g2 = g3 is their difference and holds
% at s0 automatically -- triangleFaceTest checks it EXACTLY, which is the test that the name is
% the name of THIS point. When g1 = g2 is a line and g1 = g3 is too, s0 is their crossing, rational.
    if ~isstruct(tri) || ~all(isfield(tri, {'Gn', 'Gd', 's'}))
        error('triangleFace:shape', 'expects a biconjAtoms(f, ''triangles'') entry with .Gn .Gd .s.');
    end
    if isfield(tri, 'G') && ~isempty(tri.G)
        error('triangleFace:unbounded', ['the face holds on T + cone(G) (erratum Case E); only ' ...
              'the bounded triangle is carried so far.']);
    end
    Gn = tri.Gn;  Gd = tri.Gd(:);
    if ~isequal(size(Gn), [3 6]) || numel(Gd) ~= 3
        error('triangleFace:noExactForms', ['the forms have no exact conic (formQ overflowed ' ...
              '2^53), so s0 cannot be named.']);
    end
    t12 = ratQ.conic(Gn(1,:)*Gd(2) - Gn(2,:)*Gd(1));
    t13 = ratQ.conic(Gn(1,:)*Gd(3) - Gn(3,:)*Gd(1));
    tf = struct('kind', '', 's0n', [], 's0d', [], 'Tn', [], 'Td', [], 'dp', [], 'Gn', Gn, 'Gd', Gd, ...
                's', [], 'c', [], 'T', [], 'orient', 0);
    if all(t12(1:3) == 0) && all(t13(1:3) == 0)
        [xn, xd] = ratQ.solve2([t12(4:5); t13(4:5)], [-t12(6); -t13(6)]);
        if xd < 0, xn = -xn;  xd = -xd; end
        tf.kind = 'rational';  tf.s0n = xn(:);  tf.s0d = xd;
        tf.s = xn(:) / xd;
        % the corners EXACTLY: m_i = (H_i s0n + [d e] s0d) / (s0d Gd_i), H_i = [2a b; b 2c]
        tf.Tn = zeros(3,2);  tf.Td = zeros(3,1);
        for i = 1:3
            g = Gn(i,:);
            H = [2*g(1) g(2); g(2) 2*g(3)];
            tn = ratQ.chkB((H * xn(:) + g(4:5).' * xd).', (abs(H) * abs(xn(:)) + abs(g(4:5).') * xd).', ...
                           'corner');
            [tn, td] = ratQ.canon(tn, ratQ.chk(xd * Gd(i), 'corner denominator'));
            tf.Tn(i,:) = tn;  tf.Td(i) = td;
        end
    else
        tf.kind = 'algebraic';
        tf.dp = dualPoint(t12, t13, tri.s);
        tf.s = tf.dp.s;
    end
    G = Gn(1,:) / Gd(1);  s = tf.s;
    tf.c = G(1)*s(1)^2 + G(2)*s(1)*s(2) + G(3)*s(2)^2 + G(4)*s(1) + G(5)*s(2) + G(6);
    tf.T = zeros(3,2);
    for i = 1:3
        G = Gn(i,:) / Gd(i);
        tf.T(i,:) = ([2*G(1) G(2); G(2) 2*G(3)] * s + G(4:5).').';
    end
    tf.orient = triangleFaceSign(tf, areaPoly(Gn, Gd));
    if tf.orient == 0
        error('triangleFace:flat', 'the contact triangle has zero area: its corners are collinear.');
    end
end

function A = areaPoly(Gn, Gd)
% twice the signed area of the contact triangle, times Gd1 Gd2 Gd3 > 0, as a polynomial in s:
% cross(m2 - m1, m3 - m1) = c(m1,m2) + c(m2,m3) + c(m3,m1), each c(mi,mj) over Gd_i Gd_j.
    M = arrayfun(@(i) contactPoly(Gn(i,:)), 1:3, 'UniformOutput', false);
    A = 0;
    for p = [1 2 3; 2 3 1]
        i = p(1);  j = p(2);  o = setdiff(1:3, [i j]);
        A = addArr(A, Gd(o) * crossArr(M{i}, M{j}));
    end
end

function M = contactPoly(g)
% grad of the conic [a b c d e f] as two affine polynomials in s, polyN arrays (s1 down, s2 across)
    M = {[g(4) g(2); 2*g(1) 0], [g(5) 2*g(3); g(2) 0]};
end

function C = crossArr(a, b)
    P = conv2(a{1}, b{2});  Q = conv2(a{2}, b{1});
    PB = conv2(abs(a{1}), abs(b{2}));  QB = conv2(abs(a{2}), abs(b{1}));
    sz = max(size(P), size(Q));
    C = zeros(sz);  B = zeros(sz);
    C(1:size(P,1), 1:size(P,2)) = P;   B(1:size(P,1), 1:size(P,2)) = PB;
    C(1:size(Q,1), 1:size(Q,2)) = C(1:size(Q,1), 1:size(Q,2)) - Q;
    B(1:size(Q,1), 1:size(Q,2)) = B(1:size(Q,1), 1:size(Q,2)) + QB;
    C = ratQ.chkB(C, B, 'cross product');
end

function C = addArr(A, B)
    sz = max(size(A), size(B));
    C = zeros(sz);  M = zeros(sz);
    C(1:size(A,1), 1:size(A,2)) = A;  M(1:size(A,1), 1:size(A,2)) = abs(A);
    C(1:size(B,1), 1:size(B,2)) = C(1:size(B,1), 1:size(B,2)) + B;
    M(1:size(B,1), 1:size(B,2)) = M(1:size(B,1), 1:size(B,2)) + abs(B);
    C = ratQ.chkB(C, M, 'triangle polynomial');
end
