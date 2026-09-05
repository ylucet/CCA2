function [ok, report] = checkQuaConConsistent(g, nSample)
% checkQuaConConsistent  Does this QuaCon's subdivision carry a single well-defined function?
%
% objective: detect cells that OVERLAP while carrying DIFFERENT functions. Such a mesh does not
%   define a function at all -- `eval` resolves it by first match, so the answer depends on cell
%   ORDER -- and the failure is silent, which is why this exists.
%
% [input]  g       : QuaCon;  nSample : (optional) how many points to probe, default 4000
% [output] ok      : true when no inconsistent overlap was found
%          report  : struct with .pairs (k x 2 offending cell indices), .worst (largest value gap),
%                    .point (a witness), .exactPairs (pairs proved to overlap, not merely sampled)
%
% TWO TESTS, and the split is the usual one in this codebase.
%
%   EXACT, for a pair whose constraints are ALL LINEAR: their intersection is a polyhedron, and
%   ratQ.feasible2 decides whether it has interior exactly. If it does and the two face functions
%   differ as exact rationals, the mesh is inconsistent -- decided, not estimated.
%
%   SAMPLED, as a backstop for pairs involving a CURVED constraint, which the exact test cannot
%   reach until Phase 2c's kernel lands. Sampling can only CONFIRM an overlap, never rule one out,
%   so a clean sampled result is evidence and not proof -- `ok` says what was found, not what is
%   absent.
%
% This is a TEST-side check, not a compute-path one: it is O(nf^2) and its sampled half is
% probabilistic. Use it on fixtures, not inside conj.

    if nargin < 2, nSample = 4000; end
    report = struct('pairs', zeros(0,2), 'worst', 0, 'point', [], 'exactPairs', zeros(0,2));
    ok = true;
    if g.nf < 2, return, end

    % ---- exact half: pairs whose constraints are all straight --------------------------------
    for a = 1:g.nf
        for b = a+1:g.nf
            if ratQ.eqRat(g.fN(a,:), g.fD(a), g.fN(b,:), g.fD(b)), continue, end
            ra = g.FC{a};  rb = g.FC{b};
            rows = [ra; rb];
            if any(any(g.EcQ(rows(:,1), 1:3) ~= 0)), continue, end   % curved: not decidable here
            P = rows(:,2) .* g.EcQ(rows(:,1), 4:6);
            if ratQ.feasible2(P, true)
                ok = false;
                report.exactPairs(end+1,:) = [a b]; %#ok<AGROW>
                report.pairs(end+1,:) = [a b];      %#ok<AGROW>
            end
        end
    end

    % ---- sampled half: everything else --------------------------------------------------------
    V = g.vertexCoords();
    if isempty(V), lo = [-4 -4];  hi = [4 4];
    else,          c = mean(V,1);  r = max(2, 2*max(max(abs(V - c))));  lo = c - r;  hi = c + r;
    end
    rng(20260904);
    X = lo + rand(nSample, 2) .* (hi - lo);
    for i = 1:size(X,1)
        hit = [];
        for k = 1:g.nf
            if cellHolds(g, k, X(i,:)), hit(end+1) = k; end %#ok<AGROW>
        end
        if numel(hit) < 2, continue, end
        vals = arrayfun(@(k) ratQ.evalFace(g.fN(k,:), g.fD(k), X(i,:)), hit);
        gap = max(vals) - min(vals);
        if gap > 1e-9 * max(1, max(abs(vals)))
            ok = false;
            if gap > report.worst
                report.worst = gap;  report.point = X(i,:);
            end
            [~, iw] = max(vals);  [~, ib] = min(vals);
            pr = sort([hit(iw) hit(ib)]);
            if ~any(all(report.pairs == pr, 2)), report.pairs(end+1,:) = pr; end %#ok<AGROW>
        end
    end
end

function tf = cellHolds(g, k, x)
    fc = g.FC{k};
    tf = true;
    for r = 1:size(fc,1)
        c = g.EcQ(fc(r,1), :);
        v = c(1)*x(1)^2 + c(2)*x(1)*x(2) + c(3)*x(2)^2 + c(4)*x(1) + c(5)*x(2) + c(6);
        scale = max(1, max(abs(c)) * max(1, max(abs(x))^2));
        if fc(r,2) == 0
            if abs(v) > sqrt(eps)*scale, tf = false; return, end
        elseif fc(r,2) * v < -sqrt(eps)*scale
            tf = false;  return
        end
    end
end
