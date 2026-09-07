function quadSignOracle(varargin)
% THE DIFFERENTIAL ORACLE for ratQ.quadSignOnPoly, and the only reason to trust it.
%
% The routine's contract has a SAFE direction and an UNSAFE one, and they are not symmetric:
%
%   * returning 0 ("takes both signs, or undecided") when the truth is one sign costs ONE
%     elimination that did not have to run. Harmless.
%   * returning +1 or -1 when the quadratic actually changes sign DISCARDS A BRIDGE FACE, and the
%     resulting envelope is wrong with no symptom at all.
%
% So this checks the unsafe direction hard and reports the safe one as a rate. The oracle is
% sampling: random points inside the polyhedron, plus the polyhedron's own vertices and points
% pushed far along its recession directions, since a sign change that only happens far out is
% exactly what a bounded sample misses and exactly what the recession analysis is for.
%
%   quadSignOracle              % 20000 random (conic, polyhedron) pairs
%   quadSignOracle(200000)      % longer

    N = 20000;
    if nargin >= 1, N = varargin{1}; end
    rng(2026);

    nDecided = 0;  nBoth = 0;  nDecl = 0;  nViol = 0;  nEmpty = 0;
    worst = {};
    for it = 1:N
        [c, P] = randomCase();
        if ~ratQ.feasible2(P, true), nEmpty = nEmpty + 1; continue, end
        s = ratQ.quadSignOnPoly(c, P);

        X = samplePoly(P, 400);
        if isempty(X), nEmpty = nEmpty + 1; continue, end
        q = evalConic(c, X);
        tol = 1e-9 * max(1, max(abs(q)));
        hasPos = any(q >  tol);
        hasNeg = any(q < -tol);

        if s == 0
            if hasPos && hasNeg, nBoth = nBoth + 1; else, nDecl = nDecl + 1; end
        else
            nDecided = nDecided + 1;
            bad = (s > 0 && hasNeg) || (s < 0 && hasPos);
            if bad
                nViol = nViol + 1;
                if numel(worst) < 5, worst{end+1} = {c, P, s, min(q), max(q)}; end %#ok<AGROW>
            end
        end
    end

    fprintf('\n%d cases, %d with empty or thin interior (skipped)\n', N, nEmpty);
    fprintf('  decided one-signed  : %6d\n', nDecided);
    fprintf('  0, and it DOES change sign : %6d   (the honest 0)\n', nBoth);
    fprintf('  0, but one-signed   : %6d   (%.1f%% of the one-signed cases were declined)\n', ...
            nDecl, 100 * nDecl / max(1, nDecl + nDecided));
    fprintf('  UNSOUND (claimed one sign, sample says otherwise) : %d\n', nViol);
    if nViol == 0
        fprintf('\nNo violation. The unsafe direction is clean on this many cases; the declines\n');
        fprintf('are the price of the named conservative guards, and the rate above is that price.\n');
    else
        fprintf('\nVIOLATIONS -- the first few:\n');
        for i = 1:numel(worst)
            w = worst{i};
            fprintf('  c = [%s], s = %d, sampled q in [%.3g, %.3g]\n', ...
                    strjoin(string(w{1}), ' '), w{3}, w{4}, w{5});
            disp(w{2});
        end
    end
end

% ------------------------------------------------------------------------------------------------
function [c, P] = randomCase()
% Small integer coefficients on purpose: the interesting behaviour is combinatorial (which face
% holds the minimum, whether the cone is pointed, whether a direction is exactly flat) and small
% coefficients make the exactly-flat and exactly-parallel cases actually OCCUR rather than being
% measure zero. A range that never produces a degenerate case cannot test the guards for them.
    c = randi([-4 4], 1, 6);
    k = randi([2 5]);
    P = randi([-3 3], k, 3);
    P(all(P(:,1:2) == 0, 2), 1) = 1;      % a row with no normal is not a half-plane
end

function X = samplePoly(P, n)
% Points of P: a rejection sample from a box, its vertices, and points pushed FAR along each
% recession direction. The far points matter most -- a quadratic that dips negative only at
% distance is precisely the case a box sample calls one-signed and the routine must not.
    X = zeros(0,2);
    B = 6;
    Y = (rand(20*n, 2) * 2 - 1) * B;
    Y = Y(all(Y * P(:,1:2).' + P(:,3).' >= 0, 2), :);
    X = [X; Y(1:min(end,n), :)];

    for i = 1:size(P,1)-1
        for j = i+1:size(P,1)
            A = [P(i,1:2); P(j,1:2)];
            if abs(det(A)) < 1e-12, continue, end
            v = (A \ [-P(i,3); -P(j,3)]).';
            if all(v * P(:,1:2).' + P(:,3).' >= -1e-12), X = [X; v]; end %#ok<AGROW>
        end
    end
    if isempty(X), return, end

    D = zeros(0,2);
    for i = 1:size(P,1)
        for u = [[P(i,2), -P(i,1)]; [-P(i,2), P(i,1)]].'
            if all(P(:,1:2) * u >= 0), D = [D; u.']; end %#ok<AGROW>
        end
    end
    base = X(1:min(end,20), :);
    for t = [1e2 1e4 1e6]
        for i = 1:size(D,1)
            X = [X; base + t * D(i,:)]; %#ok<AGROW>
        end
    end
end

function v = evalConic(c, X)
    v = c(1)*X(:,1).^2 + c(2)*X(:,1).*X(:,2) + c(3)*X(:,2).^2 + ...
        c(4)*X(:,1)      + c(5)*X(:,2)        + c(6);
end
