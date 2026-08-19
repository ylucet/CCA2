% 5a differential: does the rewritten isconvex agree with the committed one, and where it
% differs, which is right?  The captured merge pairs supply real operands.
%
% GROUND TRUTH for this routine is not "the union is convex" -- isconvex is only a LOCAL
% necessary probe. What it is trying to answer is: stepping along each region's OTHER edge at the
% shared vertex, does the midpoint of the two steps lie in A or in B? So the truth is computed
% directly: take the two edges, walk a small step along each INSIDE its own region, and test the
% midpoint. Done at several step sizes, since a local probe at one fixed 0.1 can be fooled.
cd(getenv('CCA2DIR')); addpath(getenv('SPDIR'));
warning('off','symbolic:sym:isAlways:TruthUnknown');
D = dir(fullfile(getenv('CCA2_DUMP_MERGE'), 'mg_*.mat'));
nPair = 0; nAgree = 0; nDiff = 0;
for i = 1:numel(D)
    S = load(fullfile(D(i).folder, D(i).name));
    A = S.A; B = S.B;
    for a = 1:size(A.ineqs,2)
        for b = 1:size(B.ineqs,2)
            if ~(A.ineqs(a) == -B.ineqs(b)), continue, end
            [nv, vxs, vys] = A.vertexOfEdge(a);
            for k = 1:min(nv,2)
                vx = vxs(k); vy = vys(k);
                got = A.isconvex(B, a, b, vx, vy);
                tru = truth(A, B, a, b, vx, vy);
                nPair = nPair + 1;
                if isempty(tru)
                    fprintf('  %s v(%s,%s): code=%d truth=UNDECIDED\n', D(i).name, char(vx), char(vy), got);
                elseif got == tru
                    nAgree = nAgree + 1;
                else
                    nDiff = nDiff + 1;
                    fprintf('  %s v(%s,%s): code=%d TRUTH=%d\n', D(i).name, char(vx), char(vy), got, tru);
                end
            end
        end
    end
end
fprintf('\npairs %d  agree %d  differ %d\n', nPair, nAgree, nDiff);

function t = truth(A, B, a, b, vx, vy)
% The midpoint test, computed WITHOUT choosing a root: walk along each region's other edge at the
% vertex by taking points of the region on that edge at several distances, and ask whether the
% midpoint lands in A or B. Undecided ([]) if either edge yields no point of its region.
    t = [];
    ea = A.getOtherEdgeAtVertex(a, [vx,vy]);
    eb = B.getOtherEdgeAtVertex(b, [vx,vy]);
    if ea == 0 || eb == 0, return, end
    if isAlways(A.ineqs(ea) == B.ineqs(eb)), t = true; return, end
    votes = [];
    for h = [0.1 0.05 0.02]
        pa = edgePt(A, ea, vx, vy, h);
        pb = edgePt(B, eb, vx, vy, h);
        if isempty(pa) || isempty(pb), continue, end
        m = (pa + pb)/2;
        votes(end+1) = A.ptFeasible(A.vars, m) || B.ptFeasible(B.vars, m); %#ok<AGROW>
    end
    if isempty(votes), return, end
    t = all(votes);            % a local probe that flips with step size is not a clean 'true'
end

function p = edgePt(R, e, vx, vy, h)
% A point of R on constraint e, a step h from (vx,vy) -- searched over BOTH offsets and, for a
% vertical edge, in y. Independent of the implementation under test.
    p = [];
    v = R.vars;
    if R.slopeIneq(e, [vx,vy]) == inf
        for dy = [h -h]
            q = [double(vx), double(vy)+dy];
            if R.ptFeasible(v, q), p = q; return, end
        end
        return
    end
    for dx = [h -h]
        xq = vx + dx;
        ey = subs(R.ineqs(e).f, v(1), xq);
        rts = region.rootsIn(ey, v(2));
        for r = 1:numel(rts)
            try, yq = double(rts(r)); catch, continue, end
            if ~isreal(yq) || ~isfinite(yq), continue, end
            if R.ptFeasible(v, [double(xq), yq]), p = [double(xq), yq]; return, end
        end
    end
end
