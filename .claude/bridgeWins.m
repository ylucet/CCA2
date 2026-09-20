function bridgeWins()
% The same question for a BRIDGE face as was asked for a corner: does it ever carry the envelope's
% VALUE, or is it there only so that its pair counts as covered? If a bridge only ever ties the best
% quadratic face, then neither the bridge's value nor -- for the value -- its predicate region is
% needed, and what remains is expressing COVERAGE.
    L = load(fullfile(tempdir, 'branchRowPruning_env.mat'));  env = L.env;
    rng(5);
    X = [rand(1500,1)*5-2.5, rand(1500,1)*5-2.5];
    nB = 0;  nTop = 0;  worst = -inf;  nOnly = 0;
    for k = 1:numel(env.bridge)
        f = env.bridge{k};
        if strcmp(f.kind, 'corner'), continue, end
        in = bridgeMember(f.pair, X);
        for i = find(in(:).')
            v = algOrRat(f, X(i,:));
            if ~isfinite(v), continue, end
            nB = nB + 1;
            bq = -inf;
            for j = 1:numel(env.quad)
                if ~inRows(env.quad(j).con, X(i,:)), continue, end
                n = env.quad(j).num;  x = X(i,:);
                bq = max(bq, (n(5)*x(1)^2 + n(6)*x(1)*x(2) + n(7)*x(2)^2 + ...
                              n(8)*x(1) + n(9)*x(2) + n(10)) / env.quad(j).den);
            end
            if isinf(bq), nOnly = nOnly + 1; continue, end
            worst = max(worst, v - bq);
            if v > bq + 1e-9, nTop = nTop + 1; end
        end
    end
    fprintf(['%d (point, bridge face) pairs | bridge STRICTLY above the best quadratic at %d | ' ...
             'no quadratic face there at all %d | largest excess %.3e\n'], nB, nTop, nOnly, worst);
end

function v = algOrRat(f, x)
    switch f.kind
        case 'algebraic'
            v = algBranch(f.FP, x, f.FB);
        case 'rational'
            q = f.fN;
            num = q(5)/2*x(1)^2 + q(6)*x(1)*x(2) + q(7)/2*x(2)^2 + q(8)*x(1) + q(9)*x(2) + q(10);
            den = f.den(1)*x(1) + f.den(2)*x(2) + f.den(3);
            if den == 0, v = NaN; return, end
            v = num/den / f.fD;
        otherwise
            v = NaN;
    end
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
