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
%   cuts = cutTriangles([], tfs, 'lines')    the RATIONAL lines the arrangement may contain
%   sub  = cutTriangles(sub, tfs)            the subdivision with the triangles cut in
%
% AFTER THE ARRANGEMENT (2026-11-01). The lines need not be among the arrangement's curves: a face
% without a row for one has its side DECIDED -- its straight rows' polygon in doubles when clear of
% the line by a margin, else cellInterior on its conic rows with the line added on each side -- and
% a face the line provably crosses is SPLIT along it, into the pieces covering its part outside T.
% Erratum Example 2's pairs are subdivided apart and overlaid (envelopeSubdivisionAlg's minOverlay),
% and 7 global lines in each pair's arrangement were cells no face needed (TODO item 1 step d).
%
% [input]  subdivision : envelopeSubdivisionAlg's output, with or without `cuts` among its curves
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
    % EACH TRIANGLE SEES THE PIECES THE EARLIER ONES LEFT (2026-11-01): until then a face split for
    % T1 was never looked at for T2, so a piece of it inside T2 kept a pair value above T2's plane.
    triF = subdivision([]);
    for t = 1:numel(tfs)
        keep = true(1, numel(subdivision));
        out = subdivision([]);
        L = triangleLines(tfs(t));
        keys = arrayfun(@(i) ratQ.conic([0 0 0 L(i,:)]), 1:size(L,1), 'UniformOutput', false);
        flip = cellfun(@(c, i) sign(c(find(c, 1))) * sign(L(i, find(L(i,:), 1))), keys, ...
                       num2cell(1:size(L,1)));       % canonical conic vs the line's own orientation
        for k = find(keep)
            sd = zeros(1, size(L,1));
            for i = 1:size(L,1)
                sd(i) = sideOf(subdivision(k).con, keys{i}, flip(i), L(i,:));
                if sd(i) < 0, break, end                          % outside one line: T misses it
            end
            if any(sd < 0), continue, end                         % outside T, or outside its box
            if strcmp(tfs(t).kind, 'rational')
                keep(k) = false;
                if all(sd > 0), continue, end                     % inside T: the plane owns it
                % CROSSED by some edge: the pieces { L_1..L_{e-1} >= 0, L_e <= 0 } cover the face
                % outside T; a piece whose L_e the face is known to be inside of is empty.
                base = subdivision(k);
                for e = 1:3
                    if sd(e) > 0, continue, end
                    p = base;
                    add = [cell2mat(keys(1:e).'), [flip(1:e-1).'; -flip(e)], zeros(e,1)];
                    p.con = [p.con; add(:, 1:size(p.con,2))];
                    out(end+1) = p; %#ok<AGROW>
                end
                continue
            end
            keep(k) = false;
            base = subdivision(k);
            for e = 1:3
                p = base;
                p.tri = [p.tri; [t*ones(e-1,1), (1:e-1).', ones(e-1,1)]; t, e, -1];
                out(end+1) = p; %#ok<AGROW>
            end
        end
        triF(end+1) = triangleItself(tfs(t), t, L, keys, flip, tfs, subdivision); %#ok<AGROW>
        subdivision = [subdivision(keep), out];
    end
    subdivision = [subdivision, triF];
end

function sd = sideOf(con, key, flip, line)
% The face's side of a cut line: +1 inside (toward T), -1 outside, 0 when the line provably or
% possibly crosses it -- a 0 only costs a split. Read off the face's own row when it has one; else
% its straight rows' polygon (a superset of the face) clear of the line by a margin; else exactly.
    r = find(all(con(:,1:6) == key, 2) & (size(con,2) < 8 | con(:, min(8, end)) == 0), 1);
    if ~isempty(r), sd = con(r,7) * flip; return, end
    sd = 0;
    V = polygonOf(con);
    if ~isempty(V)
        v = V * line(1:2).' + line(3);
        m = 1e-9 * max(1, norm(line)) * max(1, max(abs(V(:))));
        if all(v > m), sd = 1; return, end
        if all(v < -m), sd = -1; return, end
    end
    % EXACTLY: the face with the line added on one side has no interior -> it lies on the other
    keep = ~all(con(:,1:5) == 0, 2);
    C = [con(keep, 1:6), branchCol(con(keep,:))];
    rows = [(1:nnz(keep)).', con(keep, 7)];
    C(end+1,:) = [key, 0];
    j = size(C,1);
    if ~cellInterior(C, [rows; j, -flip]), sd = 1; return, end
    if ~cellInterior(C, [rows; j, flip]), sd = -1; return, end
end

function b = branchCol(con)
    if size(con,2) >= 8, b = con(:,8); else, b = zeros(size(con,1), 1); end
end

function V = polygonOf(con)
% The vertices of the face's straight, non-branch rows' polygon, within |x|, |y| <= 1e6; [] if empty
% to rounding (then the exact test decides).
    br = branchCol(con);
    st = all(con(:,1:3) == 0, 2) & br == 0 & any(con(:,4:5) ~= 0, 2);
    H = con(st, 7) .* con(st, 4:6);
    H = [H; 1 0 1e6; -1 0 1e6; 0 1 1e6; 0 -1 1e6];
    V = zeros(0,2);
    for i = 1:size(H,1)
        for j = i+1:size(H,1)
            M = H([i j], 1:2);
            if abs(det(M)) < 1e-12 * max(1, norm(M(:)))^2, continue, end
            x = (M \ -H([i j], 3)).';
            if all(H(:,1:2) * x.' + H(:,3) >= -1e-9 * max(1, norm(x)) * max(1, max(abs(H(:))))), V(end+1,:) = x; end %#ok<AGROW>
        end
    end
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
    if isfield(f, 'ii'), f.ii = like(1).ii([]); end      % T is one face: no tie row of the overlay
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
