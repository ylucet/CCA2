% End-to-end round-by-round check of Step 3's assembly for the 4-face reference case
% (f = xy over T = conv{(0,0),(3,3),(1,2)}), against the pointwise max of the per-piece
% Step 2 conjugates -- the same f*, computed the other way. Repo under test: $CCA2DIR.
DIR = getenv('CCA2DIR');
cd(DIR);
fprintf('=== repo: %s ===\n', DIR);

V = [0 0; 3 3; 1 2]; E = [1 2 1; 2 3 1; 3 1 1]; F = [1 0; 1 0; 1 0];
q = QuaPol(V, E, [0 1 0 0 0 0], F);
tAll = tic;
env = convEnvCPLQ(q);
fprintf('envelope nf=%d\n', env.nf);

p = ratPolToPlq(env);
for i = 1:p.nPieces
    tt = tic;
    p.pieces(i) = p.pieces(i).convexEnvelope;
    p.pieces(i) = p.pieces(i).conjugate;
    p.pieces(i) = p.pieces(i).maximumConjugate;
    fprintf('piece %d: %d regions (%.0f s)\n', i, size(p.pieces(i).maxConjugate,2), toc(tt));
end

R = 2*max(abs(env.V), [], 'all') + 1;
t = linspace(-R, R, 17);
[G1, G2] = meshgrid(t, t);
S = [G1(:), G2(:)];
N = size(S,1);
PV = nan(N, p.nPieces);
for k = 1:p.nPieces
    for i = 1:N
        PV(i,k) = evalFunctionNDomain(p.pieces(k).maxConjugate, S(i,:));
    end
end
FSTAR = arrayfun(@(i) convEnvCPLQTest.supBilinearOverPoly(S(i,:), V), (1:N)');
refAll = max(fillmissing(PV, 'constant', -inf), [], 2);
fprintf('CHECK per-piece max vs true f*: %d of %d disagree\n', ...
    sum(abs(refAll - FSTAR) > 1e-7 .* max(1, abs(FSTAR))), N);

acc = p.pieces(1).maxConjugate;
report('seed', acc, 1, S, PV, FSTAR);
for j = 2:p.nPieces
    tt = tic; acc = acc * p.pieces(j).maxConjugate;
    report(sprintf('* piece %d (%.0fs)', j, toc(tt)), acc, j, S, PV, FSTAR);
    tt = tic; acc = acc.maximumP(true);
    report(sprintf('maximumP %d (%.0fs)', j, toc(tt)), acc, j, S, PV, FSTAR);
end
fprintf('TOTAL %.0f s\n', toc(tAll));

function report(tag, acc, upto, S, PV, FSTAR)
    N = size(S,1);
    ref = max(fillmissing(PV(:,1:upto), 'constant', -inf), [], 2);
    nd = 0; ngap = 0; nextra = 0; first = '';
    for ii = 1:N
        got = evalMulti(acc, S(ii,:));
        haveRef = isfinite(ref(ii)); haveGot = isfinite(got);
        if haveRef && ~haveGot
            ngap = ngap + 1;
        elseif ~haveRef && haveGot
            nextra = nextra + 1;
        elseif haveRef && haveGot
            if abs(got-ref(ii)) > 1e-7*max(1,abs(ref(ii)))
                nd = nd + 1;
                if isempty(first)
                    first = sprintf('  first s=(%.4f,%.4f) got %.6f ref %.6f f*=%.6f', ...
                        S(ii,1), S(ii,2), got, ref(ii), FSTAR(ii));
                end
            end
        end
    end
    fprintf('CHECK %-18s nRegions=%3d disagree=%3d gap=%3d extra=%3d%s\n', ...
        tag, size(acc,2), nd, ngap, nextra, first);
end

% evalFunctionNDomain's exact semantics (FIRST matching region wins -- what the production
% evaluator and conjCPLQ's gate both do), extended to the state between mtimes and maximumP
% where a domain still carries two candidate functions.
function val = evalMulti(fnd, s)
    tol = 1e-6; val = NaN;
    for i = 1:numel(fnd)
        r = fnd(i).d;
        if isempty(r) || isempty(r.ineqs), continue; end
        inside = true;
        for j = 1:numel(r.ineqs)
            if double(subs(r.ineqs(j).f, r.vars, s)) > tol, inside = false; break; end
        end
        if ~inside, continue; end
        val = -inf;
        for k = 1:numel(fnd(i).f)
            val = max(val, double(subs(fnd(i).f(k).f, r.vars, s)));
        end
        return
    end
end
