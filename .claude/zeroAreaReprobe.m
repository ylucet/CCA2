% zeroAreaReprobe -- re-measure the c = -2 diagonal-kink zero-area faces (TODO "Open now").
%
% Same procedure as the 2026-11-05 measurement: biconjQ of the kinked side of
% biconjQTest.diagonalKink (unit square split by its diagonal, face 1 x^2 + y^2,
% face 2 x^2 + y^2 - 2 (y - x)), then 100000 uniform interior points plus 8000
% points on the four boundary edges; count how many points each face OWNS (eval
% keeps the first face on a tie) and whether any point is +inf.
function zeroAreaReprobe()
    root = fileparts(fileparts(mfilename('fullpath')));      % CCA2
    addpath(fullfile(root, 'MATLAB'));

    % diagonalKink(-2), exactly as MATLAB/biconjQTest.m builds it.
    q0 = [0 0 0 0, 2, 0, 2, 0, 0, 0];
    q1 = [0 0 0 0, 2, 0, 2, 2, -2, 0];
    f = QuaPol([0 0; 1 0; 1 1; 0 1], [1 2 1; 2 3 1; 3 4 1; 4 1 1; 1 3 1], [q0; q1], ...
               [1 0; 1 0; 2 0; 2 0; 2 1]);

    t0 = tic;
    h = biconjQ(f);
    fprintf('class=%s  nf=%d  (build %.1f s)\n', class(h), h.nf, toc(t0));
    for k = 1:h.nf
        fprintf('  face %d: FP=%d  fN=[%s]  fD=%g\n', k, ~isempty(h.FP{k}), ...
            mat2str(h.fN(k,:), 2), h.fD(k));
    end

    % 100000 interior + 8000 edge (2000 per boundary edge), as on 2026-11-05.
    rng(20261116);
    X = rand(100000, 2);
    t = (1:2000)'/2001;
    X = [X; t, zeros(2000,1);  t, ones(2000,1);  zeros(2000,1), t;  ones(2000,1), t];

    [val, idx] = h.eval(X);
    fprintf('+inf: %d of %d\n', sum(isinf(val)), numel(val));
    for k = 1:h.nf
        raw = ratQ.evalFace(h.fN(k,:), h.fD(k), X);
        ties = sum(isfinite(raw) & raw - val > -1e-9);
        excess = max(raw - val);
        fprintf('face %d: owned=%d  ties=%d  maxExcess=%.3e\n', k, sum(idx == k), ties, excess);
    end
end
