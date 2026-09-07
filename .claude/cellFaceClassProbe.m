function cellFaceClassProbe()
% ITEM 3's RESIDUE, MEASURED BEFORE ANYTHING IS BUILT: what do the NON-BRIDGE faces actually need?
%
% WHERE THIS SITS. The bridge machinery is done and verified. What stands between it and a complete
% `f**` is the other kind of face: over a single dual cell, with no tie, the envelope is
%
%       h_k(x) = sup over s in P_k of  <s,x> - g_k(s) ,
%
% i.e. THE CONJUGATE OF A QUADRATIC OVER A POLYHEDRON -- which is exactly what `conjQ`'s Case B
% already computes. That is the useful observation: the cells of a convex quadratic's conjugate are
% POLYHEDRAL (every cell boundary is a straight line -- it is what `cellRows` relies on throughout
% this line of work), so no new geometry is involved. Only two hypotheses of Case B fail:
%
%   * Case B takes a BOUNDED convex polygon; most dual cells are UNBOUNDED.
%   * Case B takes a STRICTLY CONVEX quadratic; over a VERTEX cell the conjugate is AFFINE
%     (Hessian 0) and over an EDGE cell it is rank 1.
%
% So the residue is not "a new algorithm" but "Case B, minus two hypotheses", and the size of the
% job is the count of cells that violate them. This measures that count rather than guessing it,
% and it also says which sub-case to build first.
%
% THE VERTEX CELL IS THE EASY ONE AND IT IS WORTH SAYING WHY. With Hessian 0 the face is affine,
% `g_k(s) = <b,s> + c`, so
%       h_k(x) = sup over s in P_k of <s, x - b>  -  c ,
% the SUPPORT FUNCTION of the cell at `x - b`. For a polyhedral cone that is piecewise LINEAR, with
% the pieces given by the cone's own normal fan, and `ratQ.recCone2` already returns the extreme
% rays that fan is built from. No conjugation machinery needed at all.

    fixtures = { ...
        'triangle (kink piece A)',  tri([0 0; 1 0; 1 1], 0), ...
        'triangle (kink piece B)',  tri([0 0; 1 1; 0 1], -2), ...
        'unit square',              poly([0 0; 1 0; 1 1; 0 1]), ...
        'hexagon',                  poly([0 0; 2 0; 3 2; 2 4; 0 4; -1 2]) };

    fprintf('\n%-26s %6s %8s %8s %8s %8s %8s %8s\n', 'piece', 'cells', 'bounded', 'unbnd', ...
            'rank2', 'rank1', 'rank0', 'CaseB');
    fprintf('%s\n', repmat('-', 1, 88));
    tot = zeros(1,6);
    for i = 1:2:numel(fixtures)
        nm = fixtures{i};  g = conjQ(fixtures{i+1});
        r = classify(g);
        fprintf('%-26s %6d %8d %8d %8d %8d %8d %8d\n', nm, g.nf, r(1), r(2), r(3), r(4), r(5), r(6));
        tot = tot + r;
    end
    fprintf('%s\n', repmat('-', 1, 88));
    fprintf('%-26s %6s %8d %8d %8d %8d %8d %8d\n', 'TOTAL', '', tot(1), tot(2), tot(3), tot(4), ...
            tot(5), tot(6));

    fprintf(['\nCaseB counts the cells that are BOUNDED and rank 2 -- the ones conjQ could take\n' ...
             'as they stand. Everything else needs one of the two hypotheses lifted, and the\n' ...
             'column that dominates says which to build first.\n']);
end

% ------------------------------------------------------------------------------------------------
function r = classify(g)
% r = [bounded, unbounded, rank2, rank1, rank0, ready-for-CaseB]
    r = zeros(1,6);
    for k = 1:g.nf
        P = cellRows(g, k);
        if isempty(P), continue, end

        % BOUNDED means the recession cone is {0}. recCone2 returns the extreme rays of that cone
        % and declines a cone that is not pointed -- a decline is counted as unbounded, which is
        % the safe direction here: it is certainly not bounded.
        [U, ok] = ratQ.recCone2(P);
        bounded = ok && isempty(U);
        if bounded, r(1) = r(1) + 1; else, r(2) = r(2) + 1; end

        % The rank of the face's Hessian, decided exactly on integers. A vertex cell gives 0 (the
        % conjugate is affine there), an edge cell 1, an interior cell 2.
        n = g.fN(k,:);
        H = [n(5), n(6); n(6), n(7)];
        if any(H(:) ~= 0)
            if ratQ.detExact(H) ~= 0, rk = 2; else, rk = 1; end
        else
            rk = 0;
        end
        r(3 + (2 - rk)) = r(3 + (2 - rk)) + 1;
        if bounded && rk == 2, r(6) = r(6) + 1; end
    end
end

function rows = cellRows(g, k)
    rows = zeros(0,3);
    FC = g.FC{k};
    for t = 1:size(FC,1)
        cc = g.EcQ(FC(t,1),:);
        if any(cc(1:3) ~= 0), rows = zeros(0,3); return, end
        if FC(t,2) == 0, continue, end
        rows(end+1,:) = FC(t,2) * cc(4:6); %#ok<AGROW>
    end
end

function f = tri(V, c)
    q = [0 0 0 0, 2, 0, 2, 0, 0, 0];
    q(8) = q(8) - c;   q(9) = q(9) + c;
    f = QuaPol(V, [1 2 1; 2 3 1; 3 1 1], q, [1 0; 1 0; 1 0]);
end

function f = poly(V)
    n = size(V,1);
    E = [(1:n).', [2:n, 1].', ones(n,1)];
    f = QuaPol(V, E, [0 0 0 0, 2, 0, 2, 0, 0, 0], repmat([1 0], n, 1));
end
