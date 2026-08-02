function checkBoxEnvelopeForSCIP()
% checkBoxEnvelopeForSCIP  What would a SCIP + QPLIB run need from CCA2?  Run it; read the table.
% Standing summary and the conclusions drawn from it: SUPPORT_MATRIX.md section 0.0.1.
%
% What would SCIP+QPLIB need from CCA2?  QPLIB's viable family (Sahinidis 1913/1922/1931/1940) has
% UNIT-BOX variable domains and objectives that are sums of x_i*x_j terms, so per term SCIP needs
% the convex envelope of a quadratic over an AXIS-ALIGNED BOX (and, after branching, a sub-box).
%
% Three questions, each measured:
%   1. Does biconj get the box envelope RIGHT, for the function shapes QPLIB actually produces?
%   2. What does it RETURN -- can the existing bridge (which reads V/E/F/f/den off a RatPol)
%      consume it?
%   3. How long does one term take?

    box = @(a,b,c,d) deal([a c; b c; b d; a d], [1 2 1; 2 3 1; 3 4 1; 4 1 1], [1 0; 1 0; 1 0; 1 0]);

    tests = {
      'unit box, f = x*y            (QPLIB off-diagonal term)', [0 1 0 1], [0 1 0 0 0 0], @(x,y) x.*y
      'unit box, f = 3xy+7x-2y+5    (term + linear part)     ', [0 1 0 1], [0 3 0 7 -2 5], @(x,y) 3*x.*y+7*x-2*y+5
      'unit box, f = x^2 - y^2      (indefinite, diagonal)   ', [0 1 0 1], [2 0 -2 0 0 0], @(x,y) x.^2-y.^2
      'unit box, f = (x^2+y^2)/2    (convex diagonal)        ', [0 1 0 1], [1 0 1 0 0 0], @(x,y) 0.5*(x.^2+y.^2)
      'SUB-box [.25,.75]x[0,1], x*y (after one branching)    ', [0.25 0.75 0 1], [0 1 0 0 0 0], @(x,y) x.*y
      'wide box [-2,3]x[-1,4], x*y  (general variable bounds)', [-2 3 -1 4], [0 1 0 0 0 0], @(x,y) x.*y
    };

    fprintf('\n%-56s %-9s %-11s %-16s %s\n', 'case', 'verdict', 'worst err', 'returns', 'time');
    fprintf('%s\n', repmat('-', 1, 115));
    for k = 1:size(tests,1)
        nm = tests{k,1}; bb = tests{k,2}; coef = tests{k,3}; f = tests{k,4};
        [V, E, F] = box(bb(1), bb(2), bb(3), bb(4));
        p = QuaPol(V, E, coef, F);
        t0 = tic;
        try
            h = p.biconj('cplq');
        catch ME
            fprintf('%-56s %-9s %-11s %-16s %.0fs\n', nm, 'ERROR', '-', ME.identifier, toc(t0));
            continue
        end
        el = toc(t0);
        [worst, at] = checkAgainstLowerHull(V, f, h);
        shape = describeShape(h);
        v = 'OK'; if worst > 2e-3, v = 'WRONG'; end
        fprintf('%-56s %-9s %-11.3g %-16s %.0fs\n', nm, v, worst, shape, el);
        if worst > 2e-3
            fprintf('%-56s worst at (%g,%g)\n', '', at(1), at(2));
        end
    end
    fprintf('\n');
end

% ------------------------------------------------------------------------------------------------
function s = describeShape(h)
% Can the SCIP bridge, which reads V/E/F/f/den off a RatPol, consume this?
    s = class(h);
    try
        hasMesh = ~isempty(h.V) && ~isempty(h.E) && ~isempty(h.F);
    catch
        hasMesh = false;
    end
    if ~hasMesh
        s = [s ' NO-MESH'];
    end
end

% ------------------------------------------------------------------------------------------------
function [worst, at] = checkAgainstLowerHull(V, f, h)
    n = 40;
    [uu,vv] = meshgrid(linspace(min(V(:,1)), max(V(:,1)), n), linspace(min(V(:,2)), max(V(:,2)), n));
    S = [uu(:), vv(:)];
    F = f(S(:,1), S(:,2));
    A = [S, ones(size(S,1),1)];
    cf = A \ F;
    if max(abs(A*cf - F)) < 1e-12 * max(1, max(abs(F)))
        env = @(x) cf(1)*x(1) + cf(2)*x(2) + cf(3);
    else
        P = [S, F]; K = convhulln(P); planes = zeros(0,3);
        for t = 1:size(K,1)
            a = P(K(t,1),:); b = P(K(t,2),:); c = P(K(t,3),:);
            nn = cross(b-a, c-a);
            if abs(nn(3)) < 1e-10, continue; end
            aa = -nn(1)/nn(3); bb = -nn(2)/nn(3); cc = a(3) - aa*a(1) - bb*a(2);
            if all(F >= aa*S(:,1) + bb*S(:,2) + cc - 1e-9*max(1,max(abs(F))))
                planes(end+1,:) = [aa bb cc]; %#ok<AGROW>
            end
        end
        env = @(x) max(planes(:,1)*x(1) + planes(:,2)*x(2) + planes(:,3));
    end
    lo = min(V); hi = max(V);
    X = lo + (hi-lo) .* [0.25 0.25; 0.5 0.5; 0.75 0.25; 0.3 0.8; 0.9 0.6; 0.1 0.1];
    worst = 0; at = [0 0];
    for i = 1:size(X,1)
        if isa(h, 'QuaParCPLQ'), got = evalFunctionNDomain(h.fnd, X(i,:)); else, got = h.eval(X(i,:)); end
        if isnan(got), got = inf; end
        e = abs(got - env(X(i,:)));
        if ~isfinite(e), e = inf; end
        if e > worst, worst = e; at = X(i,:); end
    end
end
