function biconj_consist()
% Does the ENVELOPE ever produce overlapping cells carrying different functions? conjQ had that
% defect on two fixtures; biconjQ shares assembleQuaConCells, so the invariant is worth asking of
% it too rather than assumed.
    fx = {};
    sq = [0 0; 1 0; 1 1; 0 1];  Esq = [1 2 1; 2 3 1; 3 4 1; 4 1 1];  Fsq = [1 0;1 0;1 0;1 0];
    Hs = {[2 0;0 2],[1 0;0 0],[-2 0;0 -2],[-1 0;0 0],[0 0;0 0],[0 1;1 0]};
    nm = {'PD','PSDsing','ND','NSDsing','affine','xy-McCormick'};
    for i = 1:numel(Hs)
        H = Hs{i};
        fx{end+1} = {['square-' nm{i}], QuaPol(sq, Esq, ...
            [0 0 0 0, H(1,1), H(1,2), H(2,2), 1, -1, 2], Fsq)}; %#ok<AGROW>
    end
    fx{end+1} = {'pentagon-concave', poly([0 0; 4 0; 5 3; 2 5; -1 3], [0 0 0 0 -2 0 -2 0 0 0])};
    fx{end+1} = {'L-shape-concave', poly([0 0; 2 0; 2 1; 1 1; 1 2; 0 2], [0 0 0 0 -1 0 -1 0 0 0])};
    fx{end+1} = {'tri-indef-edgeconc', poly([0 0; 2 0; 0 2], [0 0 0 0 -1 2 -1 1 -1 2])};
    V = [0 0; 1 0; 1 1; 0 1];
    fx{end+1} = {'two-triangles', QuaPol(V, [1 2 1; 2 3 1; 3 4 1; 4 1 1; 1 3 1], ...
        [0 0 0 0 -2 0 0 0 0 0; 0 0 0 0 0 0 -2 0 0 0], [1 0;1 0;2 0;2 0;2 1])};
    fx{end+1} = {'needle', QuaPol([1 2], zeros(0,3), [0 0 0 0 0 0 0 0 0 5], zeros(0,2))};
    fx{end+1} = {'segment', QuaPol([0 0; 2 0], [1 2 1], [0 0 0 0 -1 0 0 0 0 0], [0 0])};

    nbad = 0; nrun = 0;
    for i = 1:numel(fx)
        try, h = biconjQ(fx{i}{2}); catch, continue, end
        nrun = nrun + 1;
        [ok, rep] = checkQuaConConsistent(h, 1500);
        if ok
            fprintf('%-22s nf=%-4d CONSISTENT\n', fx{i}{1}, h.nf);
        else
            nbad = nbad + 1;
            fprintf('%-22s nf=%-4d OVERLAP pairs %s worst %.4g exact %d\n', fx{i}{1}, h.nf, ...
                mat2str(rep.pairs(1:min(2,size(rep.pairs,1)),:)), rep.worst, size(rep.exactPairs,1));
        end
    end
    fprintf('BICONJ-CONSISTENCY: %d of %d envelopes inconsistent\n', nbad, nrun);
end

function f = poly(V, fc)
    m = size(V,1);
    f = QuaPol(V, [(1:m).', [2:m,1].', ones(m,1)], fc, [ones(m,1), zeros(m,1)]);
end
