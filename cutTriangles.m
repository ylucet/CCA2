function [subdivision, cuts] = cutTriangles(subdivision, tfs, mode)
% cutTriangles  Put the biconjugate's TRIANGLE faces into a subdivision of the pairs' minimum -- step
%   (c) of Example 2's subdivision (TODO item 1, DECISIONS.md 2026-10-30 (fifth)).
%
% objective: co f = min( pairs' envelopes, triangle planes on their triangles ) (biconjAtoms). Inside
%   a triangle T its plane is co f and every pair's value is >= it; outside T the pairs alone are
%   co f. So the subdivision is the pairs' subdivision with T cut OUT of every piece, plus T itself.
%   Cutting is mandatory, not tidying: every evaluator takes the MAX over the faces owning a point,
%   and inside T a pair piece is above the plane.
%
% Two calls, around envelopeSubdivisionAlg:
%   cuts = cutTriangles([], tfs, 'lines')    the RATIONAL lines the arrangement must contain
%   sub  = cutTriangles(sub, tfs)            the subdivision with the triangles cut in
%
% [input]  subdivision : envelopeSubdivisionAlg's output, built WITH `cuts` among its curves
%          tfs  : struct array from triangleFace
% [output] subdivision : the same faces, each with two more fields --
%                  .tri  k x 3 rows [t e side]: side * s_e(x) >= 0, s_e the e-th side of triangle t
%                        as triangleFaceSide orients it (+1 toward the opposite corner)
%                  .tf   the triangleFace struct, for kind 'triangle'
%                plus one face per triangle: kind 'quad' for a RATIONAL s0 (an ordinary affine face
%                on three rational lines), kind 'triangle' for a named one. Every face also carries
%                .tfs, the triangle list its .tri rows index.
%          cuts : n x 3 integer lines [a b c] (a x + b y + c), for the 'lines' call
%
% A RATIONAL TRIANGLE IS THREE CURVES. Its corners are rational, so its edges are rational lines;
% in the arrangement every piece lies on one side of each, and a piece inside all three is simply
% dropped. A NAMED TRIANGLE has irrational edges, which no arrangement curve can be; it is enclosed
% in a rational BOX instead -- the corners' doubles widened, and the containment DECIDED corner by
% corner (triangleFaceSign) -- so only the pieces inside the box are cut, each into the three
% pieces { s1 < 0 }, { s1 >= 0, s2 < 0 }, { s1 >= 0, s2 >= 0, s3 < 0 } that cover its part outside
% T. Some of those are empty; that changes no value (the emptiness certificate is TODO item 3's).
    if nargin > 2 && strcmp(mode, 'lines')
        cuts = zeros(0,3);
        for t = 1:numel(tfs)
            cuts = [cuts; triangleLines(tfs(t))]; %#ok<AGROW>
        end
        subdivision = [];
        return
    end
    cuts = [];
    for k = 1:numel(subdivision)
        subdivision(k).tri = zeros(0,3);  subdivision(k).tf = [];  subdivision(k).tfs = tfs;
    end
    keep = true(1, numel(subdivision));
    out = subdivision([]);
    for t = 1:numel(tfs)
        L = triangleLines(tfs(t));
        keys = arrayfun(@(i) ratQ.conic([0 0 0 L(i,:)]), 1:size(L,1), 'UniformOutput', false);
        flip = cellfun(@(c, i) sign(c(find(c, 1))) * sign(L(i, find(L(i,:), 1))), keys, ...
                       num2cell(1:size(L,1)));       % canonical conic vs the line's own orientation
        for k = find(keep)
            sd = zeros(1, size(L,1));
            for i = 1:size(L,1)
                r = find(all(subdivision(k).con(:,1:6) == keys{i}, 2), 1);
                if isempty(r)
                    error('cutTriangles:missingCut', ['face %d has no row for cut line %s: the ' ...
                          'subdivision was not built with the cut lines among its curves.'], ...
                          k, mat2str(L(i,:)));
                end
                sd(i) = subdivision(k).con(r,7) * flip(i);
            end
            if strcmp(tfs(t).kind, 'rational')
                if all(sd > 0), keep(k) = false; end              % inside T: the plane owns it
                continue
            end
            if any(sd < 0), continue, end                         % outside the box: T misses it
            keep(k) = false;
            base = subdivision(k);
            for e = 1:3
                p = base;
                p.tri = [p.tri; [t*ones(e-1,1), (1:e-1).', ones(e-1,1)]; t, e, -1];
                out(end+1) = p; %#ok<AGROW>
            end
        end
        out(end+1) = triangleItself(tfs(t), t, L, keys, flip, tfs, subdivision); %#ok<AGROW>
    end
    subdivision = [subdivision(keep), out];
end

function L = triangleLines(tf)
% The rational lines oriented >= 0 toward T: its three edges when s0 is rational, else a box.
    if strcmp(tf.kind, 'rational')
        L = zeros(3,3);
        for e = 1:3
            i = e;  j = mod(e,3) + 1;  o = 6 - i - j;
            d = ratQ.chk(lcm(tf.Td(i), tf.Td(j)), 'corner denominator');
            a = tf.Tn(i,:) * (d / tf.Td(i));  b = tf.Tn(j,:) * (d / tf.Td(j));   % corners times d
            n = [-(b(2) - a(2)), b(1) - a(1)];                       % normal of the edge, times d
            r = ratQ.chkB([n * d, -n * a.'], [abs(n) * d, abs(n) * abs(a.')], 'edge line');   % n.(x - a/d) d
            r = r / gcd(gcd(r(1), r(2)), r(3));
            if ratQ.chkB(r(1:2) * tf.Tn(o,:).' + r(3) * tf.Td(o), ...
                         abs(r(1:2)) * abs(tf.Tn(o,:).') + abs(r(3)) * tf.Td(o), 'edge side') < 0
                r = -r;
            end
            L(e,:) = r;
        end
        return
    end
    L = boxAround(tf);
end

function L = boxAround(tf)
% A dyadic box containing the three corners, DECIDED: m_ik(s0) <= hi_k and >= lo_k for each corner
% i and coordinate k, as the sign of an affine polynomial in s0. Widened until it holds.
    T = tf.T;  w = max(1, max(max(T) - min(T)));
    for e = -8:2:0
        del = 2^(floor(log2(w)) + e);            % a power of two: the ends stay SHORT dyadics, so the
                                                % decided polynomials stay inside 2^53 (2^-30 did not)
        lo = floor((min(T) - del) / del) * del;  hi = ceil((max(T) + del) / del) * del;
        if allInside(tf, lo, hi)
            L = zeros(4,3);
            for k = 1:2
                [ln, ld] = ratQ.toRational(lo(k));  [hn, hd] = ratQ.toRational(hi(k));
                L(2*k-1, [k 3]) = [ld, -ln];        % ld x_k - ln >= 0
                L(2*k,   [k 3]) = [-hd, hn];        % hn - hd x_k >= 0
            end
            return
        end
    end
    error('cutTriangles:box', 'no dyadic box up to the triangle''s own size provably contains it.');
end

function tf = allInside(tf0, lo, hi)
    tf = true;
    for i = 1:3
        g = tf0.Gn(i,:);
        M = {[g(4) g(2); 2*g(1) 0], [g(5) 2*g(3); g(2) 0]};    % grad, over Gd(i)
        for k = 1:2
            [ln, ld] = ratQ.toRational(lo(k));  [hn, hd] = ratQ.toRational(hi(k));
            % m_ik >= lo_k  <=>  ld*M_ik - ln*Gd >= 0 ;  m_ik <= hi_k  <=>  hn*Gd - hd*M_ik >= 0
            A = ld * M{k};  A(1,1) = A(1,1) - ln * tf0.Gd(i);
            B = -hd * M{k};  B(1,1) = B(1,1) + hn * tf0.Gd(i);
            if triangleFaceSign(tf0, A) < 0 || triangleFaceSign(tf0, B) < 0
                tf = false;  return
            end
        end
    end
end

function f = triangleItself(tf, t, L, keys, flip, tfs, like)
% T as one face: rational -> an affine 'quad' on its three lines; named -> kind 'triangle'.
    f = like(1);
    f.FP = [];  f.FB = 0;  f.alg = like(1).alg([]);  f.tri = zeros(0,3);  f.tf = [];  f.tfs = tfs;
    if strcmp(tf.kind, 'rational')
        [num, den] = planeQ(tf);
        f.kind = 'quad';  f.num = num;  f.den = den;
        f.con = [cell2mat(keys(:)), flip(:), zeros(numel(keys),1)];
        f.con = f.con(:, 1:size(like(1).con, 2));
        return
    end
    f.kind = 'triangle';  f.num = zeros(1,10);  f.den = 1;  f.tf = tf;
    f.con = [cell2mat(keys(:)), flip(:), zeros(numel(keys),1)];      % the box, then the sides
    f.con = f.con(:, 1:size(like(1).con, 2));
    f.tri = [t 1 1; t 2 1; t 3 1];
    assert(size(L,1) == 4);
end

function [num, den] = planeQ(tf)
% z = <s0,x> - g1(s0) exactly, s0 = s0n/s0d, in the weighted basis (n8, n9 linear; n10 constant).
    sn = tf.s0n;  sd = tf.s0d;  g = tf.Gn(1,:);  gd = tf.Gd(1);
    c = ratQ.chkB(g(1)*sn(1)^2 + g(2)*sn(1)*sn(2) + g(3)*sn(2)^2 + (g(4)*sn(1) + g(5)*sn(2))*sd + g(6)*sd^2, ...
                  abs(g(1))*sn(1)^2 + abs(g(2)*sn(1)*sn(2)) + abs(g(3))*sn(2)^2 ...
                  + (abs(g(4)*sn(1)) + abs(g(5)*sn(2)))*sd + abs(g(6))*sd^2, 'plane constant');
    num = zeros(1,10);
    num(8) = sn(1) * sd * gd;  num(9) = sn(2) * sd * gd;  num(10) = -c;
    [num, den] = ratQ.canon(ratQ.chk(num, 'plane'), ratQ.chk(sd^2 * gd, 'plane denominator'));
end
