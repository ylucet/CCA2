function biconj_probe()
% Where does biconjQ actually stand? Same discipline as conjQ's coverage probe: enumerate the axes
% the dispatch branches on and report OK or the refusal identifier, with the answer checked against
% the envelope's DEFINITION wherever an oracle exists.
    fprintf('%-40s %s\n', 'case', 'result');
    fprintf('%s\n', repmat('-', 1, 78));

    sq = [0 0; 1 0; 1 1; 0 1];
    Esq = [1 2 1; 2 3 1; 3 4 1; 4 1 1];  Fsq = [1 0;1 0;1 0;1 0];
    tri = [0 0; 2 0; 0 2];
    Etr = [1 2 1; 2 3 1; 3 1 1];  Ftr = [1 0;1 0;1 0];

    Hs = {'PD',[2 0;0 2]; 'PSD-sing',[1 0;0 0]; 'ND',[-2 0;0 -2]; 'NSD-sing',[-1 0;0 0]; ...
          'affine',[0 0;0 0]; 'indefinite edge-concave',[0 1;1 0]; ...
          'indefinite edge-CONVEX',[1 0;0 -1]};
    for i = 1:size(Hs,1)
        H = Hs{i,2};
        fc = [0 0 0 0, H(1,1), H(1,2), H(2,2), 1, -1, 2];
        run(sprintf('square, %s', Hs{i,1}), QuaPol(sq, Esq, fc, Fsq), sq, H, [1;-1], 2);
    end
    run('triangle, indefinite edge-concave', ...
        QuaPol(tri, Etr, [0 0 0 0, -1, 2, -1, 0, 0, 0], Ftr), tri, [-1 2; 2 -1], [0;0], 0);

    % ---- the axes the table holds fixed -------------------------------------------------------
    L = [0 0; 2 0; 2 1; 1 1; 1 2; 0 2];
    m = size(L,1);
    run('NON-CONVEX L-shaped face, concave', ...
        QuaPol(L, [(1:m).', [2:m,1].', ones(m,1)], [0 0 0 0 -1 0 -1 0 0 0], ...
              [ones(m,1), zeros(m,1)]), [], [], [], []);

    run('UNBOUNDED wedge, concave', ...
        QuaPol([0 0; 1 0; 0 1], [1 2 0; 1 3 0], [0 0 0 0 -1 0 -1 0 0 0], [1 0; 0 1]), ...
        [], [], [], []);

    E2 = [1 2 1; 2 3 1; 3 4 1; 4 1 1; 1 3 1];
    F2 = [1 0; 1 0; 2 0; 2 0; 2 1];
    run('MULTI-piece, both convex', ...
        QuaPol(sq, E2, [0 0 0 0 1 0 1 0 0 0; 0 0 0 0 2 0 2 0 0 0], F2), [], [], [], []);
    run('MULTI-piece, one concave', ...
        QuaPol(sq, E2, [0 0 0 0 1 0 1 0 0 0; 0 0 0 0 -2 0 -2 0 0 0], F2), [], [], [], []);

    run('needle (dim 0)', QuaPol([1 2], zeros(0,3), [0 0 0 0 0 0 0 0 0 5], zeros(0,2)), ...
        [], [], [], []);
    run('segment (dim 1)', QuaPol([0 0; 2 0], [1 2 1], [0 0 0 0 1 0 0 0 0 0], [0 0]), ...
        [], [], [], []);
end

function run(name, f, V, H, L, k0)
    if ~f.isExact()
        fprintf('%-40s INPUT NOT EXACT\n', name);  return
    end
    try
        h = biconjQ(f);
    catch e
        id = e.identifier;
        fprintf('%-40s %s\n', name, id(max(1, find(id == ':', 1, 'last')+1):end));
        return
    end
    if isempty(V)
        fprintf('%-40s OK nf=%d (no oracle here)\n', name, h.nf);  return
    end
    rng(20260904);
    X = inHull(V, 400);
    qx = 0.5*sum((X*H).*X, 2) + X*L + k0;
    hx = h.eval(X);
    above = sum(hx > qx + 1e-9*max(1,abs(qx)));
    A = inHull(V, 200);  B = inHull(V, 200);
    nonconv = sum(h.eval((A+B)/2) > (h.eval(A) + h.eval(B))/2 + 1e-9);
    qv = 0.5*sum((V*H).*V, 2) + V*L + k0;
    vert = sum(abs(h.eval(V) - qv) > 1e-9*max(1,abs(qv)));
    fprintf('%-40s OK nf=%-3d above-f %d  nonconvex %d  vertex-gap %d\n', ...
        name, h.nf, above, nonconv, vert);
end

function X = inHull(V, n)
    m = size(V,1);  X = zeros(n,2);
    for i = 1:n
        w = rand(m,1);  w = w/sum(w);
        X(i,:) = w.' * V;
    end
end
