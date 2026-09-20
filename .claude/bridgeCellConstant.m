function bridgeCellConstant()
% Is a BRIDGE face's membership CONSTANT on the cells of the arrangement?
%
% Two constructions that tried to write the bridge's region as sign conditions are refuted
% (DECISIONS.md 2026-09-17), and both failed the same way: the cut boundary is a SEGMENT, and they
% described it with the unbounded curve it lies on. An ARRANGEMENT does not have that problem -- a
% segment is a piece of a curve bounded by the OTHER curves -- but it only helps if the region is a
% UNION OF CELLS, i.e. if the predicate never changes inside one cell.
%
% So that is the measurement, and it needs no cell enumeration: take many points, give each its SIGN
% VECTOR against every curve (the cell it lies in), and ask whether `bridgeMember` ever disagrees
% between two points of the same cell. A mixed cell refutes the approach and names the curve that is
% missing; none over many points is the evidence to build on.
    L = load(fullfile(tempdir, 'branchRowPruning_env.mat'));  env = L.env;
    C = curveTable(env);
    fprintf('%d curves (quadratic faces, corner wedges, and the bridges'' own outer bounds)\n', ...
            size(C,1));

    rng(31);
    % A uniform box samples the BIG cells; the ones that decide this question are the small ones
    % around the cut, so each bridge face's own witness seeds a dense cloud at several scales.
    X = [rand(20000,1)*6-3, rand(20000,1)*6-3];
    for k = 1:numel(env.bridge)
        f = env.bridge{k};
        if strcmp(f.kind, 'corner') || ~isfield(f, 'wit') || isempty(f.wit), continue, end
        if numel(f.wit) < 3 || f.wit(3) == 0, continue, end
        w = f.wit(1:2) / f.wit(3);              % the witness is stored as a RATIONAL point
        for sc = [1 0.1 0.01 1e-3 1e-4]
            X = [X; w + sc*(rand(3000,2)*2-1)]; %#ok<AGROW>
        end
    end
    S = signVectors(C, X);
    keep = all(S ~= 0, 2);                       % points ON a curve belong to no cell
    S = S(keep,:);  X = X(keep,:);
    [~, ~, cellOf] = unique(S, 'rows');
    fprintf('%d points in %d distinct cells\n', size(X,1), max(cellOf));

    for k = 1:numel(env.bridge)
        f = env.bridge{k};
        if strcmp(f.kind, 'corner'), continue, end
        in = bridgeMember(f.pair, X);
        mixed = 0;  nCell = 0;  worst = 0;  nIn = 0;  nOut = 0;
        for c = 1:max(cellOf)
            idx = find(cellOf == c);
            if numel(idx) < 2, continue, end
            nCell = nCell + 1;
            v = in(idx);
            if all(v)
                nIn = nIn + 1;
            elseif ~any(v)
                nOut = nOut + 1;
            else
                mixed = mixed + 1;
                worst = max(worst, min(sum(v), sum(~v)));
            end
        end
        % A count of constant cells alone can be VACUOUS -- a predicate that is false everywhere
        % is constant on every cell and proves nothing -- so both sides are reported.
        fprintf(['bridge %d (%s): %d cells | all IN %d | all OUT %d | MIXED %d | ' ...
                 'smallest minority %d\n'], k, f.kind, nCell, nIn, nOut, mixed, worst);
    end
end

% --------------------------------------------------------------------------------------------------
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

% --------------------------------------------------------------------------------------------------
function C = curveTable(env)
    C = zeros(0,7);
    for i = 1:numel(env.quad), C = addRows(C, env.quad(i).con); end
    for k = 1:numel(env.bridge)
        f = env.bridge{k};
        if strcmp(f.kind, 'corner')
            C = addRows(C, f.con);
        elseif isfield(f, 'cells') && ~isempty(f.cells)
            % The outer bound is a UNION of row blocks, so every block's curves go in.
            blocks = f.cells;
            if ~iscell(blocks), blocks = {blocks}; end
            for b = 1:numel(blocks), C = addRows(C, blocks{b}); end
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
