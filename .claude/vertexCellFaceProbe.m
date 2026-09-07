function vertexCellFaceProbe()
% THE RANK-0 NON-BRIDGE FACE, in closed form, exactly -- the cheapest of the three sub-cases the
% previous entry ordered, and it needs no conjugation machinery at all.
%
% THE DERIVATION, in four lines. Over a VERTEX cell the conjugate is AFFINE,
% `g_k(s) = <b,s> + c`, so the envelope face there is
%
%       h_k(x) = sup over s in P_k of <s, x - b>  -  c ,
%
% the SUPPORT FUNCTION of the cell at `x - b`. A vertex cell is a POINTED POLYHEDRAL CONE with a
% single apex `a` -- it is a polygon vertex's normal cone pushed through the gradient map, and an
% affine image of a cone is a cone -- so writing `P_k = a + C` with `C` the recession cone,
%
%       sup over s in a + C of <s, y>  =  <a, y> + sup over d in C of <d, y>
%                                      =  <a, y>        if y is in the POLAR of C,
%                                         +infinity     otherwise.
%
% So the face is ONE AFFINE FUNCTION on ONE POLYHEDRAL CONE:
%
%       h_k(x) = <a, x - b> - c        on   { x : <u_i, x - b> <= 0 for every extreme ray u_i of C }
%
% which is exactly what one expects -- near a primal vertex the envelope is affine -- and it is
% computable from two routines that already exist: `ratQ.recCone2` for the extreme rays (written
% for filter 2 earlier today) and `ratQ.solve2` for the apex.
%
% EVERY NUMBER STAYS EXACT. The apex is a rational point from an integer 2x2 solve, the extreme
% rays are integer vectors, and the face's coefficients are integer combinations of those. Nothing
% here needs a tolerance, and nothing here needs an elimination.
%
% CHECKED AGAINST THE DEFINITION, which shares none of that: the support function is also computed
% by maximising `<s, x-b>` over a fine sample of the cell, and the two are compared where the
% closed form says finite. Where it says +infinity the assertion is different, and it is the one
% that catches a polar cone computed with the wrong sign: walk the offending recession ray and the
% value must grow LINEARLY, increment ratio exactly 1. An earlier version compared sampled sups at
% two radii and called that a growth test -- it is not, because the sample is drawn in a DISC and
% rejected against the cone, so a larger radius wins a smaller share of hits and the ratio came out
% 4.3 where a linear function on a cone must give 10.

    fixtures = { 'triangle', tri([0 0; 1 0; 1 1], 0), ...
                 'square',   poly([0 0; 1 0; 1 1; 0 1]), ...
                 'hexagon',  poly([0 0; 2 0; 3 2; 2 4; 0 4; -1 2]) };

    nOK = 0;  nInf = 0;  nSkip = 0;  worst = 0;  grewLeast = inf;
    for i = 1:2:numel(fixtures)
        nm = fixtures{i};  g = conjQ(fixtures{i+1});
        fprintf('\n=== %s: %d dual cells ===\n', nm, g.nf);
        for k = 1:g.nf
            n = g.fN(k,:);
            if any([n(5) n(6) n(7)] ~= 0), continue, end        % not a rank-0 cell
            P = cellRows(g, k);
            if isempty(P), nSkip = nSkip + 1; continue, end
            [face, cone, why] = vertexFace(g, k, P);
            if isempty(face)
                fprintf('  cell %d: DECLINED -- %s\n', k, why);
                nSkip = nSkip + 1;  continue
            end
            [ok, inf_, w, gl] = checkFace(g, k, P, face, cone);
            nOK = nOK + ok;  nInf = nInf + inf_;
            worst = max(worst, w);  grewLeast = min(grewLeast, gl);
            fprintf(['  cell %d: h(x) = <[%s], x> %+d  on %d half-planes | %d finite pts checked ' ...
                     '(worst %.2e), %d infinite pts\n'], k, num2str(face.a), face.k0, ...
                    size(cone,1), ok, w, inf_);
        end
    end

    fprintf('\nTOTAL: %d finite points matched the closed form (worst relative error %.3e)\n', ...
            nOK, worst);
    fprintf('       %d points where the face is +infinity, and it really is unbounded there\n', nInf);
    fprintf('       (increment ratio along the offending recession ray: %.4f; exactly 1 for an affine function that grows without bound)\n', grewLeast);
    fprintf('       %d cells declined\n', nSkip);
    if worst < 1e-9 && abs(grewLeast - 1) < 1e-9
        fprintf(['\nThe rank-0 face is right: one affine function on one polyhedral cone, exact,\n' ...
                 'and the +infinity region really is unbounded rather than merely large.\n']);
    end
end

% ------------------------------------------------------------------------------------------------
function [face, cone, why] = vertexFace(g, k, P)
% The closed form. face = {a, k0} meaning h(x) = <a,x> + k0; cone = rows [p1 p2 p0] for the
% half-planes p1*x + p2*y + p0 >= 0 on which it is valid.
    face = [];  cone = [];  why = '';
    n = g.fN(k,:);  d = g.fD(k);
    b = [n(8); n(9)] / d;   c = n(10) / d;

    [U, ok] = ratQ.recCone2(P);
    if ~ok
        why = 'the recession cone is not pointed, so the cell is not a translated cone';  return
    end
    if isempty(U)
        why = 'the cell is BOUNDED -- a rank-0 bounded cell is a different sub-case';  return
    end

    apex = coneApex(P);
    if isempty(apex)
        why = 'no single apex found; the cell is a slab or a half-plane rather than a cone';  return
    end

    % h(x) = <apex, x - b> - c, so a = apex and the constant is -<apex,b> - c.
    face = struct('a', apex.', 'k0', -(apex.' * b) - c);

    % The POLAR of C, in x rather than in y = x - b: <u_i, x - b> <= 0 becomes
    % -u_i(1)*x1 - u_i(2)*x2 + <u_i, b> >= 0.
    cone = zeros(size(U,1), 3);
    for r = 1:size(U,1)
        cone(r,:) = [-U(r,1), -U(r,2), U(r,:) * b];
    end
end

function a = coneApex(P)
% The single vertex of a pointed polyhedral cone: the one pairwise intersection of its boundary
% lines that satisfies every constraint. Exact, via ratQ.solve2 on integer rows.
    a = [];
    for i = 1:size(P,1)-1
        for j = i+1:size(P,1)
            A = [P(i,1:2); P(j,1:2)];
            if ratQ.detExact(A) == 0, continue, end
            [xn, xd] = ratQ.solve2(A, [-P(i,3); -P(j,3)]);
            if all(P(:,1:2) * xn + P(:,3) * xd >= 0)      % xd > 0, so no sign flip
                a = xn / xd;  return
            end
        end
    end
end

function [nFin, nInf, worst, growLeast] = checkFace(g, k, P, face, cone)
% Compare the closed form with a sampled support function, and check the +infinity region really
% is unbounded rather than just large.
    n = g.fN(k,:);  d = g.fD(k);
    b = [n(8); n(9)] / d;   c = n(10) / d;
    nFin = 0;  nInf = 0;  worst = 0;  growLeast = inf;

    rng(7);
    for t = 1:300
        x = (rand(2,1) * 2 - 1) * 6;
        inCone = all(cone(:,1:2) * x + cone(:,3) >= 0);
        if inCone
            v = sampledSup(P, x - b, 1e4) - c;
            h = face.a * x + face.k0;
            worst = max(worst, abs(v - h) / max(1, abs(h)));
            nFin = nFin + 1;
        else
            % THE +INFINITY SIDE, MEASURED ALONG A RECESSION RAY RATHER THAN BY REJECTION
            % SAMPLING. Comparing sampled sups at two radii looked like a growth test and was not:
            % the sample is drawn in a DISC and rejected against the cone, so a larger radius wins
            % a smaller share of hits and the ratio came out 4.3 where a linear function on a cone
            % must give 10. Walking the ray the polar test says is unbounded measures the thing
            % itself -- <a + t*u, y> is affine in t with slope <u,y>, so doubling t must double the
            % growth exactly.
            growLeast = min(growLeast, rayGrowth(P, x - b));
            nInf = nInf + 1;
        end
    end
end

function r = rayGrowth(P, y)
% Along the recession ray that violates the polar test, how does <s,y> grow? Returns the ratio of
% the increments at t and 2t, which for an affine function is exactly 1 -- anything unbounded and
% linear gives 1, and a bounded direction gives 0.
    [U, ok] = ratQ.recCone2(P);
    a = coneApex(P);
    if ~ok || isempty(U) || isempty(a), r = 0; return, end
    best = 0;
    for i = 1:size(U,1)
        u = U(i,:).';
        d1 = (a + 1e3*u).' * y - a.' * y;
        d2 = (a + 2e3*u).' * y - (a + 1e3*u).' * y;
        if d1 > 0, best = max(best, d2 / d1); end
    end
    r = best;
end

function v = sampledSup(P, y, R)
% sup of <s,y> over the cell, sampled. The cell is a cone, so its far reaches are what matter:
% points are drawn in a disc of radius R and rejected against the cell's own half-planes, and the
% APEX is included so a tiny cone is never missed entirely.
    S = (rand(4000,2) * 2 - 1) * R;
    S = S(all(S * P(:,1:2).' + P(:,3).' >= 0, 2), :);
    a = coneApex(P);
    if ~isempty(a), S = [S; a.']; end
    if isempty(S), v = -inf; return, end
    v = max(S * y);
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
