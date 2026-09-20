function sideProbe()
% Does `algCurveSide` answer the definition? The definition is the sign of q(x) - value_FB(x), and
% `algBranch` computes that value in doubles -- so the two must agree wherever the value exists, at
% points simple enough to be exact rationals.
    L = load(fullfile(tempdir, 'branchRowPruning_env.mat'));  env = L.env;
    rng(101);
    % SIMPLE rationals. MEASURED (`ceilingProbe`): FP is degree 4 in x and 4 in y with coefficients
    % of size 4 to 16, so clearing a point with denominator 2^e costs about |coef| (3*2^e)^4 2^(4e),
    % which stays inside 2^53 up to e = 4 and overflows at e = 6. k/8 is comfortably inside; an
    % arbitrary double is a rational too, but one with a 2^52 denominator.
    X = round((rand(4000,2)*5-2.5) * 8) / 8;
    total = 0;  agree = 0;  onCurve = 0;  refused = 0;
    for k = 1:numel(env.bridge)
        f = env.bridge{k};
        if ~strcmp(f.kind, 'algebraic'), continue, end
        in = bridgeMember(f.pair, X);
        if ~any(in), continue, end
        P = X(in,:);
        for j = 1:min(6, numel(env.quad))
            q = env.quad(j);
            try
                [s, ~] = algCurveSide(f.FP, f.FB, q.num, q.den, P);
            catch ME
                if strcmp(ME.identifier, 'ratQ:overflow'), refused = refused + 1; continue, end
                rethrow(ME)
            end
            vb = arrayfun(@(i) algBranch(f.FP, P(i,:), f.FB), 1:size(P,1)).';
            vq = ratQ.evalFace(q.num, q.den, P);
            d = vq - vb;
            ok = isfinite(d) & ~isnan(s) & abs(d) > 1e-7;
            total = total + sum(ok);
            agree = agree + sum(sign(d(ok)) == s(ok));
            onCurve = onCurve + sum(s == 0);
        end
    end
    fprintf('exact side vs the definition: %d of %d agree | on the curve %d | refused %d\n', ...
            agree, total, onCurve, refused);
end
