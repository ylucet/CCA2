function cutProbe()
% Does `bridgeQuadCut`'s curve actually separate a bridge face from a quadratic one? The test is the
% definition: at points where both own x, sign(R) must agree with sign(bridge value - quad value),
% up to ONE global orientation per (bridge, quad) pair -- R contains the boundary, and which side is
% which is what a consumer has to record.
    L = load(fullfile(tempdir, 'branchRowPruning_env.mat'));  env = L.env;
    rng(77);
    X = [rand(6000,1)*5-2.5, rand(6000,1)*5-2.5];
    nPair = 0;
    for k = 1:numel(env.bridge)
        f = env.bridge{k};
        if strcmp(f.kind, 'corner') || ~strcmp(f.kind, 'algebraic'), continue, end
        in = bridgeMember(f.pair, X);
        if ~any(in), continue, end
        P = X(in,:);
        vb = arrayfun(@(i) algBranch(f.FP, P(i,:), f.FB), 1:size(P,1)).';
        for j = 1:numel(env.quad)
            q = env.quad(j);
            own = ownsRows(q.con, P);
            if sum(own) < 20, continue, end
            vq = ratQ.evalFace(q.num, q.den, P(own,:));
            d = vb(own) - vq;
            [R, info] = bridgeQuadCut(f.FP, q.num, q.den);
            r = evalPoly2(R, P(own,:));
            good = abs(d) > 1e-9 & abs(r) > 1e-12;
            if sum(good) < 20, continue, end
            agree = max(mean(sign(r(good)) == sign(d(good))), ...
                        mean(sign(r(good)) == -sign(d(good))));
            % THE READING UNDER TEST: R vanishes where q meets EITHER root of FP, so its sign should
            % be the PRODUCT over both branches times the leading coefficient. If that holds, the
            % product agrees with R everywhere while the single branch does not -- which says the
            % edge needs a BRANCH INDEX, exactly as a degenerate conic's line does one degree down.
            Pg = P(own,:);  Pg = Pg(good,:);
            prod = zeros(size(Pg,1),1);
            for i = 1:size(Pg,1)
                [~, i2] = algBranch(f.FP, Pg(i,:), f.FB);
                rr = i2.allReal{1};
                qq = ratQ.evalFace(q.num, q.den, Pg(i,:));
                lead = polyAt(f.FP(:,:,3), Pg(i,:));
                pr = lead;
                for t = 1:numel(rr), pr = pr * (qq - rr(t)); end
                prod(i) = pr;
            end
            rg = r(good);
            ok2 = abs(prod) > 0;
            both = max(mean(sign(rg(ok2)) == sign(prod(ok2))), ...
                       mean(sign(rg(ok2)) == -sign(prod(ok2))));
            nPair = nPair + 1;
            if nPair <= 8
                fprintf(['bridge %d vs quad %2d: deg %d, %4d pts | ONE branch %.3f | ' ...
                         'BOTH branches %.3f\n'], ...
                        k, j, info.degree, sum(good), agree, both);
            end
        end
    end
    fprintf('%d (bridge, quadratic) pairs compared\n', nPair);
end

function v = evalPoly2(R, X)
    v = zeros(size(X,1), 1);
    for i = 1:size(R,1)
        for j = 1:size(R,2)
            if R(i,j) == 0, continue, end
            v = v + R(i,j) * X(:,1).^(i-1) .* X(:,2).^(j-1);
        end
    end
end

function own = ownsRows(con, X)
    own = true(size(X,1),1);
    for t = 1:size(con,1)
        c = con(t,1:6);
        if all(c(1:5) == 0)
            if con(t,7)*c(6) < 0, own(:) = false; end
            continue
        end
        br = 0;
        if size(con,2) >= 8, br = con(t,8); end
        if br ~= 0
            v = ratQ.pairLineEval(c, br, X);
        else
            v = c(1)*X(:,1).^2 + c(2)*X(:,1).*X(:,2) + c(3)*X(:,2).^2 + ...
                c(4)*X(:,1) + c(5)*X(:,2) + c(6);
        end
        own = own & (con(t,7)*v >= -1e-9);
    end
end

function v = polyAt(P, x)
    v = 0;
    for i = 1:size(P,1)
        for j = 1:size(P,2)
            if P(i,j) ~= 0, v = v + P(i,j) * x(1)^(i-1) * x(2)^(j-1); end
        end
    end
end
