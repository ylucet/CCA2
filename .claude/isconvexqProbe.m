function isconvexqProbe()
% What does the new exact isConvexQ answer on every fixture in the corpus, and does it agree with
% the legacy doubles QuaPol.isConvex where the legacy one runs at all?
%
% The legacy routine is the DIFFERENTIAL ORACLE here (umbrella CLAUDE.md 6), not the reference: two
% latent smells in it are recorded in TODO.md 2026-09-05, so a disagreement is a question, not a
% verdict. The verdict comes from the corpus names -- examplesNonconvex must be false,
% examplesDiscontinuous must be false, and the 8 fixtures biconjQ refuses must be true.
    groups = {};
    P = QuaPol.examples();            groups{end+1} = {'examples', P, []};
    P2 = QuaPol.examples2();          groups{end+1} = {'examples2', P2, []};
    PN = QuaPol.examplesNonconvex();  groups{end+1} = {'examplesNonconvex', PN, false};
    PD = QuaPol.examplesDiscontinuous(); groups{end+1} = {'examplesDiscontinuous', PD, false};

    nAgree = 0; nDisagree = 0; nErr = 0;
    for g = 1:numel(groups)
        name = groups{g}{1};  L = groups{g}{2};  expect = groups{g}{3};
        fprintf('\n== %s (%d) ==\n', name, numel(L));
        for k = 1:numel(L)
            f = L{k};
            [nw, w] = tryNew(f);
            old = tryOld(f);
            flag = ' ';
            if ischar(nw)
                nErr = nErr + 1; flag = 'E';
            elseif ischar(old)
                flag = '-';
            elseif nw == old
                nAgree = nAgree + 1;
            else
                nDisagree = nDisagree + 1; flag = '!';
            end
            bad = '';
            if ~isempty(expect) && islogical(nw) && nw ~= expect, bad = '   <== EXPECTED FALSE'; end
            fprintf('%s %-18s new=%-6s old=%-6s %s%s\n', flag, sprintf('%s{%d}',name,k), ...
                    str(nw), str(old), w, bad);
        end
    end

    fprintf('\n== the 8 fixtures biconjQ refuses ==\n');
    names = {'oneNorm'}; objs = {QuaPol.oneNorm()};
    for k = [1 5 8 16 18 19]
        names{end+1} = sprintf('examples(%d)',k);  objs{end+1} = P{k}; %#ok<AGROW>
    end
    names{end+1} = 'examples2(9)';  objs{end+1} = P2{9};
    for i = 1:numel(objs)
        [nw, w] = tryNew(objs{i});
        bad = ''; if ~isequal(nw,true), bad = '   <== EXPECTED TRUE'; end
        fprintf('  %-16s new=%-6s %s%s\n', names{i}, str(nw), w, bad);
    end

    fprintf('\nagree=%d disagree=%d newErrors=%d\n', nAgree, nDisagree, nErr);
end

function [tf, why] = tryNew(f)
    try, [tf, why] = isConvexQ(f);
    catch e, tf = shortId(e); why = e.message(1:min(60,end));
    end
end

function tf = tryOld(f)
    try, tf = logical(f.isConvex());
    catch e, tf = shortId(e);
    end
end

function s = shortId(e)
    id = e.identifier;
    if isempty(id), s = 'ERR'; else, s = id(max(1,find(id==':',1,'last')+1):end); end
end

function s = str(v)
    if ischar(v), s = ['<' v '>']; elseif v, s = 'true'; else, s = 'false'; end
end
