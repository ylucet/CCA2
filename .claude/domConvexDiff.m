function domConvexDiff()
% Where does the new doubles domain test disagree with the one createDom records today
% (orderEdges(0))? Measured over the whole fixture corpus BEFORE wiring it in, so that every
% changed expectation is a decision rather than a surprise.
%
% The exact isConvexQ is shown alongside where it can run, since it settles any disagreement.
    groups = {'examples', QuaPol.examples(); 'examples2', QuaPol.examples2(); ...
              'examplesNonconvex', QuaPol.examplesNonconvex(); ...
              'examplesDiscontinuous', QuaPol.examplesDiscontinuous()};
    nDiff = 0;
    fprintf('%-26s %-8s %-8s %-8s\n', 'fixture', 'old', 'new', 'exact');
    fprintf('%s\n', repmat('-', 1, 54));
    for g = 1:size(groups,1)
        L = groups{g,2};
        for k = 1:numel(L)
            f = L{k};
            old = f.dom.isConvex;
            new = domainIsConvexD(f.V, f.E, f.F);
            ex = '-';
            try
                [t, ~] = isConvexQ(f);
                if t, ex = 'true'; else, ex = 'false'; end   % NOTE: this is f, not just dom
            catch
            end
            mark = ' ';
            if logical(old) ~= logical(new), mark = '!'; nDiff = nDiff + 1; end
            fprintf('%s %-24s %-8d %-8d %-8s\n', mark, ...
                    sprintf('%s{%d}', groups{g,1}, k), old, new, ex);
        end
    end
    fprintf('\ndisagreements: %d\n', nDiff);
    fprintf(['(the "exact" column is isConvexQ on the WHOLE function, so a false there does not ' ...
             'by itself\n mean the domain is bad -- read it only for the rows marked !)\n']);
end
