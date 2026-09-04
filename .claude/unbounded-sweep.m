function unb_sweep()
% conjQ on UNBOUNDED pieces, against a sampled sup.
%
% THE ORACLE IS ONE-SIDED ON PURPOSE, and plqCheck.m already uses this discipline: a sup taken over
% sampled points of the domain is a LOWER bound on the true sup, so `f* < sampled` is a definite
% defect while the other direction is expected. On top of that, where conjQ reports +infinity the
% sampled values must actually GROW along the ray -- otherwise the routine has thrown away a finite
% answer, which the lower bound alone cannot catch.
    cases = {};
    cases{end+1} = {'oneNorm (|x|_1)', QuaPol.oneNorm()};
    cases{end+1} = {'energy on a wedge', wedge([0 0], [1 0], [0 1], [0 0 0 0 0 0 0 0 0 0])};
    cases{end+1} = {'affine on a wedge', wedge([0 0], [1 0], [0 1], [0 0 0 0 0 0 0 -1 -2 3])};
    cases{end+1} = {'concave on a wedge', wedge([0 0], [1 0], [0 1], [0 0 0 0 -2 0 -2 1 1 0])};
    cases{end+1} = {'affine on a half-strip', halfStrip()};

    for c = 1:numel(cases)
        nm = cases{c}{1};  f = cases{c}{2};
        if ~f.isExact(), fprintf('%-24s NOT EXACT\n', nm); continue, end
        try
            g = conjQ(f);
        catch e
            fprintf('%-24s REFUSED %s\n', nm, e.identifier);
            continue
        end
        rng(20260904);
        S = [randn(400,2)*2; 0 0; 1 1; -1 -1];
        val = g.eval(S);

        % sample the domain: vertices, convex combinations, and long rides down every ray
        X = sampleDomain(f);
        low = -inf(size(S,1),1);
        for i = 1:size(S,1)
            low(i) = max(S(i,:)*X.' - evalQ(f, X).');
        end
        under = val < low - 1e-7*max(1,abs(low));
        % where conjQ says +inf, the sampled sup must be growing
        % Where conjQ reports +inf the sup must really diverge, so compare two reaches: the
        % sampled value must keep RISING, and by a factor, not by a fixed amount -- an absolute
        % threshold is meaningless when |s| is small.
        f4 = reachMax(f, S, 1e4);
        f6 = reachMax(f, S, 1e6);
        infs = isinf(val);
        notGrowing = infs & ~(f6 > f4 + 1e-6*max(1,abs(f4)));
        fprintf('%-24s nf=%-3d ne=%-3d  BELOW-SAMPLED %d   INF-not-growing %d\n', ...
            nm, g.nf, g.ne, sum(under), sum(notGrowing));
    end
end

function v = evalQ(f, X)
% The FUNCTION, not face 1's formula. Using f.f(1,:) for every point was wrong on any multi-piece
% input -- it evaluated |x|_1 as the single affine piece of the first quadrant, so the "sampled
% sup" exceeded the true one and 74 correct values looked like defects.
    v = f.eval(X);
end

function m = reachMax(f, S, reach)
    X = sampleDomain(f, reach);
    m = zeros(size(S,1),1);
    for i = 1:size(S,1)
        m(i) = max(S(i,:)*X.' - evalQ(f, X).');
    end
end

function X = sampleDomain(f, reach)
    if nargin < 2, reach = 8; end
    V = f.V;
    X = V;
    for a = 1:size(V,1)
        for b = a+1:size(V,1)
            for t = [0.25 0.5 0.75]
                X(end+1,:) = (1-t)*V(a,:) + t*V(b,:); %#ok<AGROW>
            end
        end
    end
    for j = 1:size(f.E,1)
        if f.E(j,3) ~= 0, continue, end
        base = V(f.E(j,1),:);  d = V(f.E(j,2),:) - V(f.E(j,1),:);
        d = d / norm(d);
        for t = [1 10 100 reach]
            X(end+1,:) = base + t*d; %#ok<AGROW>
        end
    end
end

function f = wedge(apex, d1, d2, fc)
% a two-ray cone: vertices are the apex and one point along each ray
    % The two rays bound the wedge from OPPOSITE sides, so their F columns are opposite -- F(j,:)
    % is [left, right] of edge j and 0 means +inf. Writing [1 0; 1 0] puts the face on the left of
    % both, which no convex wedge satisfies, and QuaPol.eval then cannot locate its own interior
    % (measured: q(1,1) = Inf inside the first quadrant).
    V = [apex; apex + d1; apex + d2];
    E = [1 2 0; 1 3 0];                       % two rays leaving the apex
    F = [1 0; 0 1];
    f = QuaPol(V, E, fc, F);
end

function f = halfStrip()
% {0 <= x <= 1, y >= 0}: two parallel rays and the segment joining their bases
    V = [0 0; 1 0; 0 1; 1 1];
    E = [1 2 1; 1 3 0; 2 4 0];
    F = [1 0; 0 1; 1 0];
    f = QuaPol(V, E, [0 0 0 0 0 0 0 1 1 0], F);
end
