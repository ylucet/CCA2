% Repeat of the same-hyperplane measurement, now that Step 2 is EXACT, and with the question
% the first version conflated split in two:
%   symbolic  : does merge's own test, ineqs(i) == -ineqs(j), find a shared facet?
%   hyperplane: do the two carry the same hyperplane with OPPOSITE orientation (numeric rows)?
%   touching  : is the intersection of the two linear relaxations a SEGMENT, not a point or empty?
% Only "touching" pairs are ones a merge could ever be right about.
cd(getenv('CCA2DIR'));
warning('off','symbolic:sym:isAlways:TruthUnknown');
setenv('CCA2_A45_SPLIT','');
V = [0 0; 3 3; 1 2];
E = [1 2 1; 2 3 1; 3 1 1]; F = [1 0; 1 0; 1 0];
p = ratPolToPlq(convEnvCPLQ(QuaPol(V, E, [0 1 0 0 0 0], F)));
for i = 1:p.nPieces
    p.pieces(i) = p.pieces(i).convexEnvelope.conjugate.maximumConjugate;
end
acc = p.pieces(1).maxConjugate * p.pieces(2).maxConjugate;
acc = acc.maximumP(true);
n = numel(acc);
grp = zeros(1,n); g = 0;
for i = 1:n
    if grp(i), continue, end
    g = g + 1; grp(i) = g;
    for j = i+1:n
        if grp(j), continue, end
        if isAlways(simplifyFraction(acc(i).f(1).f - acc(j).f(1).f) == 0, 'Unknown', 'false')
            grp(j) = g;
        end
    end
end
fprintf('cells = %d, distinct functions = %d\n', n, g);

cnt = containers.Map();
for i = 1:n
    for j = i+1:n
        if grp(i) ~= grp(j), continue, end
        A = acc(i).d; B = acc(j).d;
        if isempty(A) || isempty(B), continue, end
        symYes = false;
        for a = 1:size(A.ineqs,2)
            for b = 1:size(B.ineqs,2)
                if A.ineqs(a) == -B.ineqs(b), symYes = true; break, end
            end
            if symYes, break, end
        end
        [hypYes, dirn] = sharedHyperplane(A, B);
        touch = touching(A, B, dirn);
        key = sprintf('sym=%d hyp=%d touch=%s', symYes, hypYes, touch);
        if isKey(cnt,key), cnt(key) = cnt(key)+1; else, cnt(key) = 1; end
    end
end
k = keys(cnt);
for t = 1:numel(k)
    fprintf('  %-28s : %d\n', k{t}, cnt(k{t}));
end

function [tf, dirn] = sharedHyperplane(A, B)
    tf = false; dirn = [];
    [Aa, ba, la] = A.linearForm; [Ab, bb, lb] = B.linearForm;
    for i = 1:numel(ba)
        if ~la(i), continue, end
        u = [Aa(i,:), ba(i)]; if norm(u) < 1e-12, continue, end; u = u/norm(u);
        for j = 1:numel(bb)
            if ~lb(j), continue, end
            v = [Ab(j,:), bb(j)]; if norm(v) < 1e-12, continue, end; v = v/norm(v);
            if norm(u + v) < 1e-7
                tf = true; dirn = [-Aa(i,2), Aa(i,1)]; dirn = dirn/norm(dirn); return
            end
        end
    end
end

function s = touching(A, B, dirn)
% 'no' if the two linear relaxations do not meet, 'point' if they meet in a single point,
% 'seg' if the meeting has positive extent -- the only case a merge could be right about.
    [Aa, ba, la] = A.linearForm; [Ab, bb, lb] = B.linearForm;
    M = [Aa(la,:); Ab(lb,:)]; q = [ba(la); bb(lb)];
    if isempty(dirn), dirn = [1 0]; end
    [hi, sh] = region.maxLinear(M, q,  dirn);
    [lo, sl] = region.maxLinear(M, q, -dirn);
    if sh == -1 || sl == -1, s = 'no'; return, end
    if sh ~= 0 || sl ~= 0, s = 'seg'; return, end     % unbounded along the line
    if hi + lo > 1e-7, s = 'seg'; else, s = 'point'; end
end
