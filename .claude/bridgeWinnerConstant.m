function bridgeWinnerConstant()
% Does the WINNER change inside a cell, once bridge faces are in the max?
%
% A corner never carries the value, so it needed no boundary against a quadratic face. A bridge DOES
% (22 of 60 points, up to 1.2e-1), so if the winner changes inside a cell the mesh needs the curve
% {bridge value = quadratic value} -- which is FP(x, y, q(x,y)) = 0, a rational curve of degree well
% above two, and not an edge any conic subdivision can hold. If instead each cell has ONE winner
% throughout, the arrangement already separates them and nothing new is needed.
    L = load(fullfile(tempdir, 'branchRowPruning_env.mat'));  env = L.env;
    C = curveTable(env);
    rng(41);
    X = [rand(30000,1)*6-3, rand(30000,1)*6-3];
    for k = 1:numel(env.bridge)
        f = env.bridge{k};
        if strcmp(f.kind, 'corner') || numel(f.wit) < 3 || f.wit(3) == 0, continue, end
        w = f.wit(1:2) / f.wit(3);
        for sc = [1 0.1 0.01]
            X = [X; w + sc*(rand(4000,2)*2-1)]; %#ok<AGROW>
        end
    end
    S = signVectors(C, X);
    keep = all(S ~= 0, 2);
    S = S(keep,:);  X = X(keep,:);
    [~, ~, cellOf] = unique(S, 'rows');

    val = zeros(size(X,1), 0);  isBridge = [];  pairOf = [];
    for i = 1:numel(env.quad)
        val(:,end+1) = faceVal(env.quad(i), X); %#ok<AGROW>
        isBridge(end+1) = false; %#ok<AGROW>
        pairOf(end+1) = env.quad(i).pair(1); %#ok<AGROW>
    end
    for k = 1:numel(env.bridge)
        f = env.bridge{k};
        pairOf(end+1) = f.pairIndex; %#ok<AGROW>
        v = -inf(size(X,1),1);
        in = ownsIt(f, X);
        for i = find(in(:).'), v(i) = faceValAlg(f, X(i,:)); end
        val(:,end+1) = v; %#ok<AGROW>
        isBridge(end+1) = ~strcmp(f.kind, 'corner'); %#ok<AGROW>
    end
    % a face that does not own the point cannot win there
    for i = 1:numel(env.quad)
        own = ownsRows(env.quad(i).con, X);
        val(~own, i) = -inf;
    end

    [best, who] = max(val, [], 2);
    nCell = 0;  mixed = 0;  mixedWithBridge = 0;  bridgeWins = 0;
    for c = 1:max(cellOf)
        idx = find(cellOf == c);
        if numel(idx) < 2, continue, end
        if all(isinf(best(idx))), continue, end
        nCell = nCell + 1;
        u = unique(who(idx));
        if any(isBridge(u)), bridgeWins = bridgeWins + 1; end
        if numel(u) > 1
            mixed = mixed + 1;
            if any(isBridge(u))
                mixedWithBridge = mixedWithBridge + 1;
                % WHICH competition is it? Two faces of ONE pair would mean the region rows are
                % what is missing -- a pair's candidates partition its own region -- while faces of
                % DIFFERENT pairs mean a genuine max between pairs, whose boundary is the curve
                % {bridge value = quadratic value} and is not a conic.
                fprintf(['  cell %d: winners %s pairs %s bridge %s' newline], c, ...
                        mat2str(u(:).'), mat2str(pairOf(u)), mat2str(isBridge(u)));
            end
        end
    end
    fprintf(['%d cells with a finite value | more than ONE winner inside %d ' ...
             '(a bridge among them %d) | a bridge wins somewhere in %d\n'], ...
            nCell, mixed, mixedWithBridge, bridgeWins);
end

function v = faceVal(f, X)
    v = ratQ.evalFace(f.num, f.den, X);        % the WEIGHTED basis, not a hand-rolled one
end

function v = faceValAlg(f, x)
    switch f.kind
        case {'algebraic', 'corner'}
            v = algBranch(f.FP, x, f.FB);
        case 'rational'
            q = f.fN;
            num = q(5)/2*x(1)^2 + q(6)*x(1)*x(2) + q(7)/2*x(2)^2 + q(8)*x(1) + q(9)*x(2) + q(10);
            den = f.den(1)*x(1) + f.den(2)*x(2) + f.den(3);
            if den == 0, v = -inf; return, end
            v = num/den / f.fD;
        otherwise
            v = -inf;
    end
    if ~isfinite(v), v = -inf; end
end

function in = ownsIt(f, X)
    if strcmp(f.kind, 'corner')
        in = false(size(X,1),1);
        for i = 1:size(X,1)
            mu = f.M \ (f.p - X(i,:).');
            in(i) = all(mu >= -1e-9);
        end
    else
        in = bridgeMember(f.pair, X);
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

function S = signVectors(C, X)
    S = zeros(size(X,1), size(C,1));
    for j = 1:size(C,1)
        c = C(j,1:6);  br = C(j,7);
        if br ~= 0
            v = ratQ.pairLineEval(c, br, X);
        else
            v = c(1)*X(:,1).^2 + c(2)*X(:,1).*X(:,2) + c(3)*X(:,2).^2 + ...
                c(4)*X(:,1) + c(5)*X(:,2) + c(6);
        end
        v(abs(v) < 1e-9) = 0;
        S(:,j) = sign(v);
    end
end

function C = curveTable(env)
    C = zeros(0,7);
    for i = 1:numel(env.quad), C = addRows(C, env.quad(i).con); end
    for k = 1:numel(env.bridge)
        f = env.bridge{k};
        if strcmp(f.kind, 'corner')
            C = addRows(C, f.con);
        elseif isfield(f, 'cells') && ~isempty(f.cells)
            b = f.cells;
            if ~iscell(b), b = {b}; end
            for t = 1:numel(b), C = addRows(C, b{t}); end
        end
    end
end

function C = addRows(C, con)
    for t = 1:size(con,1)
        c = con(t,1:6);
        if all(c(1:5) == 0), continue, end
        br = 0;
        if size(con,2) >= 8, br = con(t,8); end
        key = [c, br];
        if isempty(find(all(C == key, 2), 1)), C(end+1,:) = key; end %#ok<AGROW>
    end
end
