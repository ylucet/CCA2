function ceilingProbe()
% What does the EXACT side test actually cost? `algCurveSide` clears a rational point's denominators
% through FP, and FP is a discriminant -- so its degree and coefficient size decide which points are
% affordable under ratQ's 2^53 ceiling. Measured, not guessed.
    L = load(fullfile(tempdir, 'branchRowPruning_env.mat'));  env = L.env;
    for k = 1:numel(env.bridge)
        f = env.bridge{k};
        if ~strcmp(f.kind, 'algebraic'), continue, end
        FP = f.FP;
        I = size(FP,1) - 1;  J = size(FP,2) - 1;
        mx = max(abs(FP(:)));
        fprintf('bridge %d: FP degree %d in x, %d in y, %d in z | largest |coefficient| %.3e\n', ...
                k, I, J, size(FP,3)-1, mx);
        % what a point with denominator 2^e costs: |coef| * (|n|)^deg * d^(deg) at |x| <= 3
        for e = [0 1 2 4 6]
            d = 2^e;
            worst = mx * (3*d)^I * d^J;
            fprintf('   denominator 2^%d: worst term ~ %.3e  (%s 2^53)\n', ...
                    e, worst, ternary(worst < 2^53, 'inside', 'OVER'));
        end
    end
end

function s = ternary(c, a, b)
    if c, s = a; else, s = b; end
end
