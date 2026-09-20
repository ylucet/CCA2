function cornerWins2()
% Does a CORNER face ever carry the envelope's value? Asked pointwise, so no cell can be missed:
% at each point take the faces whose own region contains it and see which is largest. If a corner
% never wins, the corner/quadratic boundary never has to be drawn and no new kind of edge is needed.
    L = load(fullfile(tempdir, 'branchRowPruning_env.mat'));  env = L.env;
    rng(3);
    X = [rand(4000,1)*8-4, rand(4000,1)*8-4];
    nCov = 0;  nCornerTop = 0;  worst = 0;  nCornerOnly = 0;
    for i = 1:size(X,1)
        x = X(i,:);
        best = -inf;  bestCorner = false;  bestQuad = -inf;  bestCornerVal = -inf;
        for k = 1:numel(env.quad)
            if ~inRows(env.quad(k).con, x), continue, end
            v = quadVal(env.quad(k), x);
            bestQuad = max(bestQuad, v);
            if v > best, best = v;  bestCorner = false; end
        end
        for k = 1:numel(env.bridge)
            f = env.bridge{k};
            if ~strcmp(f.kind, 'corner') || isempty(f.con), continue, end
            if ~inRows(f.con, x), continue, end
            v = cornerVal(f, x);
            bestCornerVal = max(bestCornerVal, v);
            if v > best, best = v;  bestCorner = true; end
        end
        if isinf(best), continue, end
        nCov = nCov + 1;
        if bestCorner
            nCornerTop = nCornerTop + 1;
            if isinf(bestQuad), nCornerOnly = nCornerOnly + 1; end
        end
        if ~isinf(bestCornerVal)
            % THE DECISIVE NUMBER: how far a corner's value ever rises above the best quadratic
            % face owning the same point. At or below zero everywhere means a corner never carries
            % the envelope's VALUE -- it only witnesses that its pair is covered.
            worst = max(worst, bestCornerVal - bestQuad);
        end
    end
    fprintf(['%d of %d points covered | a CORNER is the largest at %d of them ' ...
             '(no quadratic face there at all: %d) | largest excess over a quadratic face %.3e\n'], ...
            nCov, size(X,1), nCornerTop, nCornerOnly, worst);
end

function tf = inRows(con, x)
    tf = true;
    for t = 1:size(con,1)
        c = con(t,1:6);
        if all(c(1:5) == 0), if con(t,7)*c(6) < 0, tf = false; return, end, continue, end
        br = 0;
        if size(con,2) >= 8, br = con(t,8); end
        if br ~= 0
            v = ratQ.pairLineEval(c, br, x);
        else
            v = c(1)*x(1)^2 + c(2)*x(1)*x(2) + c(3)*x(2)^2 + c(4)*x(1) + c(5)*x(2) + c(6);
        end
        if con(t,7)*v < 0, tf = false; return, end
    end
end

function v = quadVal(f, x)
    n = f.num;
    v = (n(5)*x(1)^2 + n(6)*x(1)*x(2) + n(7)*x(2)^2 + n(8)*x(1) + n(9)*x(2) + n(10)) / f.den;
end

function v = cornerVal(f, x)
    n = f.src.n;  d = f.src.d;  s = f.s;
    N = [n(5) n(6); n(6) n(7)];  m = [n(8); n(9)];  q = n(10);
    v = x*s - (s.'*N*s/2 + m.'*s + q)/d;
end
