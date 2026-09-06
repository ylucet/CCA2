function kinkDualProbe()
% RUNG 1: OBSERVE. What does `conjQ` actually produce for the fixtures `biconjQ` refuses with
% notImplemented? If f* comes back with only LINE edges, it is a polyhedral subdivision wearing
% QuaCon's clothing, and the double conjugation f** = (f*)* is reachable with machinery that
% already exists -- conjQ itself -- rather than needing the general "conjugate of a QuaCon".
%
% No editing, no fixing. Just the numbers.
    fixtures = {};
    fixtures{end+1} = {'square2 c=-2 (PD kink)',   kink('square2', -2)};
    fixtures{end+1} = {'poly2   c=+2 (PD kink)',   kink('poly2',   +2)};
    fixtures{end+1} = {'poly2   c=-2 (PD kink)',   kink('poly2',   -2)};
    fixtures{end+1} = {'square1 indefinite',       single([2 0; 0 -2])};
    fixtures{end+1} = {'square1 PD (convex)',      single([2 0; 0  2])};

    for i = 1:numel(fixtures)
        name = fixtures{i}{1};  f = fixtures{i}{2};
        fprintf('\n=== %s ===\n', name);
        try
            g = conjQ(f);
        catch e
            fprintf('  conjQ REFUSED: %s\n', e.identifier);
            continue
        end
        kinds = cell(g.ne,1);
        for j = 1:g.ne, kinds{j} = g.edgeKind(j); end
        u = unique(kinds);
        fprintf('  f* : nv=%d ne=%d nf=%d   edge kinds:', g.nv, g.ne, g.nf);
        for k = 1:numel(u), fprintf(' %s x%d', u{k}, sum(strcmp(kinds,u{k}))); end
        fprintf('\n');
        allLines = all(strcmp(kinds, 'line'));
        fprintf('  ALL EDGES ARE LINES: %d\n', allLines);
        if allLines
            fprintf('  -> f* is polyhedral, so (f*)* is reachable by conjQ if f* can be\n');
            fprintf('     re-expressed as a QuaPol. Is dom f* the whole plane? ');
            fprintf('nv=%d (0 would mean unbounded/full)\n', g.nv);
        end
    end
end

function f = single(H)
    V = [0 0; 1 0; 1 1; 0 1];
    f = QuaPol(V, [1 2 1; 2 3 1; 3 4 1; 4 1 1], ...
               [0 0 0 0, H(1,1), H(1,2), H(2,2), 0, 0, 0], [1 0; 1 0; 1 0; 1 0]);
end

function f = kink(shape, c)
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
