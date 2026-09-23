function notDegenerateProbe(shape, c, op)
% IS maxQuaPar:notDegenerate REACHABLE FROM AN IN-CONTRACT INPUT? One fixture, one operator.
%
% The guard's own message says a pair built to trip it is out of contract, so only the wired
% symbolic pipeline on a legal QuaPol can settle it (DECISIONS 2026-09-06). conjCPLQ catches
% maxQuaPar: identifiers and falls back, so the top-level outcome cannot show a firing; the guard
% appends to $CCA2_GUARD_LOG instead. Driven one process per (fixture, op) by
% .claude/notDegenerateRun.sh, because solve/simplify cannot be interrupted from M-code.
%
%   notDegenerateProbe('square2', -2, 'conjCPLQ')
    f = kink(shape, c);
    t0 = tic;
    try
        g = feval(op, f);
        fprintf('OUTCOME ok %s  (%.0f s)\n', class(g), toc(t0));
    catch e
        fprintf('OUTCOME %s  (%.0f s)\n  %s\n', e.identifier, toc(t0), e.message);
    end
end

function f = kink(shape, c)
% Verbatim from .claude/reachabilityProbe.m: two PD pieces that agree on the shared edge for
% every c, the kink swept convex / flat / non-convex by c. DECISIONS 2026-09-05.
    q0 = [0 0 0 0, 2, 0, 2, 0, 0, 0];
    if shape == "square2"
        V = [0 0; 1 0; 1 1; 0 1];  lin = [-1 1 0];
    else
        V = [-2 0; 2 0; 1 2; -1 2];  lin = [0 1 0];
    end
    E = [1 2 1; 2 3 1; 3 4 1; 4 1 1; 1 3 1];
    F = [1 0; 1 0; 2 0; 2 0; 2 1];
    q1 = q0;  q1(8) = q1(8) + c*lin(1);  q1(9) = q1(9) + c*lin(2);  q1(10) = q1(10) + c*lin(3);
    f = QuaPol(V, E, [q0; q1], F);
end
