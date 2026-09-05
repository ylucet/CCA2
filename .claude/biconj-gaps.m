function biconj_gaps()
% For each fixture the legacy envelope answers and the exact one does not: WHICH refusal, and would
% asserting convexity (fIsConvex) unblock it? That separates "needs an algorithm" from "needs a
% convexity test", which are very different amounts of work.
    names = {}; objs = {};
    for n = {'oneNorm'}
        names{end+1} = n{1};  objs{end+1} = QuaPol.(n{1})(); %#ok<AGROW>
    end
    P = QuaPol.examples();
    for k = [1 5 8 16 18 19]
        names{end+1} = sprintf('examples(%d)',k);  objs{end+1} = P{k}; %#ok<AGROW>
    end
    P2 = QuaPol.examples2();
    names{end+1} = 'examples2(9)';  objs{end+1} = P2{9};

    fprintf('%-16s %-8s %-7s %-22s %s\n','fixture','pieces','rays','refusal','with fIsConvex=true');
    fprintf('%s\n', repmat('-',1,78));
    for i = 1:numel(names)
        f = objs{i};
        nray = sum(f.E(:,3)==0);
        e1 = idOf(@() biconjQ(f));
        g = f;  g.fIsConvex = true;
        e2 = idOf(@() biconjQ(g));
        fprintf('%-16s %-8d %-7d %-22s %s\n', names{i}, size(f.fN,1), nray, e1, e2);
    end
end

function s = idOf(fn)
    try, h = fn(); s = sprintf('OK nf=%d', h.nf);
    catch e, id = e.identifier; s = id(max(1,find(id==':',1,'last')+1):end); end
end
