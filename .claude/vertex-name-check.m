function vname_check()
% Every named vertex must LIE ON both curves it names, and be realisable. That is the contract of
% the name: [edgeA edgeB rootIdx] means "the rootIdx-th intersection of those two curves", so if
% the realised point is not on both, the name is wrong.
    E = [1 2 1; 2 3 1; 3 1 1];  F = [1 0; 1 0; 1 0];
    cases = {};
    cases{end+1} = {'fold of two triangles', maxQ( ...
        conjQ(QuaPol([0 0; 1 0; 0 1], E, [0 0 0 0 1 0 1 0 0 0], F)), ...
        conjQ(QuaPol([1 0; 1 1; 0 1], E, [0 0 0 0 4 1 3 -2 1 0], F)))};
    sq = [0 0; 1 0; 1 1; 0 1];  Esq = [1 2 1; 2 3 1; 3 4 1; 4 1 1];  Fsq = [1 0;1 0;1 0;1 0];
    cases{end+1} = {'PSD-sing square', conjQ(QuaPol(sq, Esq, [0 0 0 0 1 0 0 1 -1 2], Fsq))};
    cases{end+1} = {'indefinite square', conjQ(QuaPol(sq, Esq, [0 0 0 0 1 0 -1 1 -1 2], Fsq))};
    u = [2;1]; H = u*u.';  V5 = [0 0; 3 0; 3 2; 1 3; 0 2];  m = size(V5,1);
    cases{end+1} = {'PSD-sing pentagon', conjQ(QuaPol(V5, [(1:m).', [2:m,1].', ones(m,1)], ...
        [0 0 0 0, H(1,1), H(1,2), H(2,2), 1, -1, 2], [ones(m,1), zeros(m,1)]))};

    for c = 1:numel(cases)
        nm = cases{c}{1};  g = cases{c}{2};
        V = g.vertexCoords();
        worst = 0;  bad = 0;  curved = 0;
        for i = 1:g.nv
            a = g.Vname(i,1);  b = g.Vname(i,2);
            if any(g.EcQ(a,1:3) ~= 0) || any(g.EcQ(b,1:3) ~= 0), curved = curved + 1; end
            for e = [a b]
                cc = g.EcQ(e,:);
                v = cc(1)*V(i,1)^2 + cc(2)*V(i,1)*V(i,2) + cc(3)*V(i,2)^2 ...
                  + cc(4)*V(i,1) + cc(5)*V(i,2) + cc(6);
                s = max(1, max(abs(cc)) * max(1, max(abs(V(i,:)))^2));
                r = abs(v)/s;
                worst = max(worst, r);
                if r > 1e-9, bad = bad + 1; end
            end
        end
        fprintf('%-24s nv=%-4d curved %-4d  OFF-CURVE %d  worst residual %.2e\n', ...
            nm, g.nv, curved, bad, worst);
    end
end
