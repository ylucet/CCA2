function sd = triangleFaceSide(tf, xn, xd)
% triangleFaceSide  On which side of each edge of a triangle face's contact triangle a RATIONAL
%   point lies, decided exactly although the corners are irrational.
%
% [input]  tf     : as triangleFace returns it
%          xn, xd : the point xn/xd, xn 2 integers, xd a positive integer
% [output] sd     : 1 x 3 in {-1, 0, +1}, for the edges (m1,m2), (m2,m3), (m3,m1): +1 on the side
%                   of the opposite corner, 0 ON the edge's line, -1 beyond it. The point is in the
%                   closed triangle exactly when all three are >= 0.
%
% THE SIDE IS A POLYNOMIAL IN s0. With m_i = M_i(s)/Gd_i (M_i affine, the gradient of form i),
%
%   cross(m_j - m_i, x - m_i) * Gd_i Gd_j xd = Gd_i c(M_j, xn) - Gd_j c(M_i, xn) - xd c(M_j, M_i),
%
% c(a,b) = a1 b2 - a2 b1: of degree 2 in s, integer coefficients; its sign at s0, times the
% triangle's orientation, is the side. triangleFaceSign decides it.
    xn = ratQ.chk(double(xn(:)), 'point numerator');
    xd = ratQ.chk(double(xd), 'point denominator');
    if numel(xn) ~= 2 || ~isscalar(xd) || xd <= 0
        error('triangleFaceSide:point', 'expects xn (2 integers) and xd > 0.');
    end
    M = arrayfun(@(i) contactPoly(tf.Gn(i,:)), 1:3, 'UniformOutput', false);
    sd = zeros(1,3);
    for e = 1:3
        i = e;  j = mod(e, 3) + 1;
        O = addArr(addArr(tf.Gd(i) * crossX(M{j}, xn), -tf.Gd(j) * crossX(M{i}, xn)), ...
                   -xd * crossArr(M{j}, M{i}));
        sd(e) = tf.orient * triangleFaceSign(tf, O);
    end
end

function M = contactPoly(g)
% grad of the conic [a b c d e f] as two affine polynomials in s (polyN arrays, s1 down, s2 across)
    M = {[g(4) g(2); 2*g(1) 0], [g(5) 2*g(3); g(2) 0]};
end

function C = crossX(M, xn)
    C = addArr(ratQ.chk(xn(2) * M{1}, 'side'), ratQ.chk(-xn(1) * M{2}, 'side'));
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
