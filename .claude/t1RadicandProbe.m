% .claude/t1RadicandProbe.m -- which quadratic extensions ONE polygon needs (DECISIONS 2026-08-20).
% Which quadratic extensions does ONE polygon's Step 1 actually need?  exactQ allows ONE
% squarefree d and RAISES when two are mixed (DECISIONS 2026-08-19 night, T1/T2).
tris = { [0 0; 2 0; 2.5 1.5], [2.5 1.5; 0 0; 0.5 1], [0 0; 1 0; 2 1], [0 0; 1 1; 3 2] };
names = {'quad tri A', 'quad tri B', 'A.4 fixture', 'DESIGN 3-convex'};
allrad = [];
for t = 1:numel(tris)
    sub = splitTightTriangleSym(tris{t});
    rad = [];
    for k = 1:numel(sub)
        V = sym(sub{k});
        c = children(V(:).');
        for i = 1:numel(V)
            e = V(i);
            [n, d] = numden(e);
            ch = symvar(e);            % no variables; use string scan instead
            s = char(simplify(e));
            % pull every  ^(1/2)  radicand out of the printed form
            tok = regexp(s, '(\d+)\^\(1/2\)', 'tokens');
            for q = 1:numel(tok), rad(end+1) = str2double(tok{q}{1}); end %#ok<AGROW>
        end
    end
    rad = unique(rad);
    allrad = unique([allrad, rad]);
    fprintf('%-16s  %d sub-triangles  radicands: %s\n', names{t}, numel(sub), mat2str(rad));
end
fprintf('DISTINCT RADICANDS ACROSS THESE FIXTURES: %s\n', mat2str(allrad));
