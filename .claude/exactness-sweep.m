function exactness_sweep()
% Does the CURRENT double pipeline produce face coefficients / edge conics that are exactly
% rational with a small denominator, and does it produce vertices that are not?
%
% This decides the scope of Phase 1: making `f`/`den`/`Ec` exact is safe only if every value the
% existing pipeline stores in them round-trips through ratQ.fromDouble. CONJ_FIELD_PROOF.md
% Theorem 1 predicts exactly this split (coefficients rational, vertices not), so the sweep is
% also a check on the theorem against the code.

    cases = {};
    sq = [0 0; 1 0; 1 1; 0 1];  sqE = [1 2 1; 2 3 1; 3 4 1; 4 1 1];  sqF = [1 0;1 0;1 0;1 0];
    cases{end+1} = {'energy on the unit square',  sq, sqE, [0 0 0 0 1 0 1 0 0 0], sqF};
    cases{end+1} = {'xy on the unit square',      sq, sqE, [0 0 0 0 0 1 0 0 0 0], sqF};
    cases{end+1} = {'x^2-y^2 on the unit square', sq, sqE, [0 0 0 0 1 0 -1 0 0 0], sqF};
    cases{end+1} = {'3xy+7x-2y+5 on the square',  sq, sqE, [0 0 0 0 0 3 0 7 -2 5], sqF};
    tri = [0 0; 2 0; 0 3];  triE = [1 2 1; 2 3 1; 3 1 1];  triF = [1 0;1 0;1 0];
    cases{end+1} = {'energy on a triangle',       tri, triE, [0 0 0 0 1 0 1 0 0 0], triF};
    cases{end+1} = {'xy on a triangle',           tri, triE, [0 0 0 0 0 1 0 0 0 0], triF};
    cases{end+1} = {'concave on a triangle',      tri, triE, [0 0 0 0 -1 0 -1 0 0 0], triF};
    wide = [0 0; 5 0; 5 2; 0 2];
    cases{end+1} = {'xy on a wide box',           wide, sqE, [0 0 0 0 0 1 0 0 0 0], sqF};

    fprintf('%-30s %-10s %8s %8s %8s\n', 'case', 'kind', 'f bad', 'Ec bad', 'V bad');
    fprintf('%s\n', repmat('-', 1, 70));
    tot = zeros(1,3);
    for i = 1:numel(cases)
        c = cases{i};
        try
            f = QuaPol(c{2}, c{3}, c{4}, c{5});
            g = f.conj();
        catch e
            fprintf('%-30s %-10s   conj failed: %s\n', c{1}, '-', e.identifier);
            continue
        end
        nf = countBad(g.f);
        ne = 0;
        if ~isempty(g.Ec), ne = countBad(g.Ec); end
        nv = 0;
        if ~isempty(g.V),  nv = countBad(g.V);  end
        tot = tot + [nf ne nv];
        fprintf('%-30s %-10s %4d/%-3d %4d/%-3d %4d/%-3d\n', c{1}, g.kind(), ...
            nf, size(g.f,1), ne, size(g.Ec,1), nv, size(g.V,1));
    end
    fprintf('%s\nTOTAL BAD  f=%d  Ec=%d  V=%d\n', repmat('-', 1, 70), tot(1), tot(2), tot(3));
end

function n = countBad(M)
% how many ROWS of M fail to convert to an exact rational vector
    n = 0;
    for k = 1:size(M,1)
        try
            ratQ.fromDouble(M(k,:));
        catch
            n = n + 1;
        end
    end
end
