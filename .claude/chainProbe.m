function chainProbe()
% RUNG 1 for the CHAIN item: which fixtures are chains, what shape are they, and which of the two
% chain refusals does each actually reach? `ambiguousChain` fires before `thinEnvelope`, so it is
% worth knowing whether `chainNotImplemented` is reachable at all before writing anything for it.
    P = QuaPol.examples2();
    fprintf('%-14s %3s %3s %3s %6s %5s  %-28s %s\n', ...
            'fixture', 'nv', 'ne', 'nf', 'rowsfN', 'dim', 'biconjQ', 'collinear?');
    fprintf('%s\n', repmat('-', 1, 95));
    for k = 1:numel(P)
        f = P{k};
        col = '-';
        if f.dom.dim < 2 && f.nv >= 2
            if QuaPol.isCollinear(f.V, sqrt(eps)), col = 'collinear'; else, col = 'NOT collinear'; end
        end
        try
            h = biconjQ(f);  res = sprintf('OK nf=%d', h.nf);
        catch e
            id = e.identifier;  res = id(max(1,find(id==':',1,'last')+1):end);
        end
        fprintf('examples2{%2d}  %3d %3d %3d %6d %5d  %-28s %s\n', ...
                k, f.nv, f.ne, f.nf, size(f.fN,1), f.dom.dim, res, col);
    end

    fprintf(['\nA chain is dim < 2 with more than one edge. `ambiguousChain` fires when the mesh\n' ...
             'carries more than one function row, and eval NEEDS one row per edge for a chain, so\n' ...
             'the question is whether any chain can reach thinEnvelope at all.\n']);
end
