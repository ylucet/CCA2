function coverage_probe()
% What does conjQ actually accept, and what does it refuse? Enumerates the input space of a QuaPol
% along the axes that the dispatch branches on -- DOMAIN shape x HESSIAN class, plus the degenerate
% inputs -- and reports OK or the refusal identifier for each cell. No claim of completeness is
% made here; the table is the claim.
    Hs = {'PD', 'PSD-sing', 'indefinite', 'ND', 'NSD-sing', 'affine'};
    Hm = { [2 0; 0 2], [1 0; 0 0], [1 0; 0 -1], [-2 0; 0 -2], [-1 0; 0 0], [0 0; 0 0] };

    fprintf('%-22s', 'domain \ Hessian');
    for i = 1:numel(Hs), fprintf('%-16s', Hs{i}); end
    fprintf('\n%s\n', repmat('-', 1, 22 + 16*numel(Hs)));

    doms = {'full plane', 'bounded triangle', 'bounded square', 'unbounded wedge', 'half-strip'};
    for d = 1:numel(doms)
        fprintf('%-22s', doms{d});
        for i = 1:numel(Hm)
            H = Hm{i};
            fc = [0 0 0 0, H(1,1), H(1,2), H(2,2), 1, -1, 2];
            f = buildDomain(doms{d}, fc);
            fprintf('%-16s', tryConj(f));
        end
        fprintf('\n');
    end

    fprintf('\n');
    % multi-piece and degenerate inputs
    extra = {};
    V = [0 0; 1 0; 1 1; 0 1];
    extra{end+1} = {'multi-face bounded', QuaPol(V, [1 2 1; 2 3 1; 3 4 1; 4 1 1; 1 3 1], ...
        [0 0 0 0 1 0 1 0 0 0; 0 0 0 0 4 1 3 -2 1 0], [1 0; 1 0; 2 0; 2 0; 1 2])};
    extra{end+1} = {'multi-face unbounded', QuaPol.oneNorm()};
    extra{end+1} = {'needle (dim 0)',   QuaPol([0 0], zeros(0,3), [0 0 0 0 0 0 0 0 0 1], zeros(0,2))};
    extra{end+1} = {'segment (dim 1)',  QuaPol([0 0; 1 0], [1 2 1], [0 0 0 0 0 0 0 0 0 1], [0 0])};
    extra{end+1} = {'cubic numerator',  QuaPol(V, [1 2 1; 2 3 1; 3 4 1; 4 1 1], ...
        [1 0 0 0 1 0 1 0 0 0], [1 0; 1 0; 1 0; 1 0])};
    extra{end+1} = {'inexact input',    QuaPol([sqrt(2) 0 1 0 0 0])};
    for e = 1:numel(extra)
        fprintf('%-22s%s\n', extra{e}{1}, tryConj(extra{e}{2}));
    end
end

function s = tryConj(f)
    if isempty(f), s = '(no fixture)'; return, end
    try
        g = conjQ(f);
        s = sprintf('OK nf=%d', g.nf);
    catch e
        id = e.identifier;
        s = id(max(1, find(id == ':', 1, 'last')+1):end);
    end
end

function f = buildDomain(name, fc)
    switch name
        case 'full plane'
            f = QuaPol(fc(5:10));
        case 'bounded triangle'
            f = QuaPol([0 0; 2 0; 0 3], [1 2 1; 2 3 1; 3 1 1], fc, [1 0; 1 0; 1 0]);
        case 'bounded square'
            f = QuaPol([0 0; 1 0; 1 1; 0 1], [1 2 1; 2 3 1; 3 4 1; 4 1 1], fc, ...
                       [1 0; 1 0; 1 0; 1 0]);
        case 'unbounded wedge'
            f = QuaPol([0 0; 1 0; 0 1], [1 2 0; 1 3 0], fc, [1 0; 0 1]);
        case 'half-strip'
            f = QuaPol([0 0; 1 0; 0 1; 1 1], [1 2 1; 1 3 0; 2 4 0], fc, [1 0; 0 1; 1 0]);
        otherwise
            f = [];
    end
end
