function [Vname, info] = algCrvVertices(g)
% algCrvVertices  Name the VERTICES of an AlgCrv: where two of a face's own edges meet.
%
% objective: an AlgCrv evaluates without vertices, but everything that wants its TOPOLOGY -- an edge
%   list, a figure, the paper's classification tables -- needs them. `assembleAlgCrv` leaves the
%   table empty; this fills it.
%
% [input]  g     : AlgCrv
% [output] Vname : nv x 5, each row [edgeA tableA edgeB tableB rootIdx] -- see below
%          info  : .candidates, .named, .refused, .seconds
%
% ==================================================================================================
% A VERTEX IS A NAME, NOT A COORDINATE, and that is the same convention one degree up
%
% `QuaCon.Vname` names a vertex as [edgeA edgeB rootIdx]: the rootIdx-th intersection of two edge
% CONICS, in `conicMeet`'s canonical order, with no coordinate stored because the coordinate is
% irrational. Here an edge may come from either of two tables, so the name carries which:
%
%       [edgeA tableA edgeB tableB rootIdx] ,   table 0 = EcQ (a conic), 1 = EcA (a curve).
%
% and the intersections are ordered lexicographically by (x, y), which is `conicMeet`'s order too.
%
% WHY THE EXACT ROUTE IS AFFORDABLE, measured before it was written (DECISIONS.md 2026-09-21). The
% algebraic edges of the real fixture have degree 4, so the three kinds of pair cost:
%
%       conic x conic   `conicMeet`, the existing degree-4 kernel
%       conic x curve   resultant of degree 4, coefficients of size 1 to 4, instant
%       curve x curve   resultant of degree 12 to 16, coefficients about 1e5, 1.5 to 2.6 s
%
% -- all of it inside ratQ's 2^53 ceiling, which was the open question. `ratQ.resultantVarBig` does
% the elimination in BigInteger and the answer comes back as ordinary integers.
%
% MEMBERSHIP IS DECIDED IN DOUBLES, and that is deliberate rather than a gap. `assembleQuaConCells`
% already does the same for a curved corner, with the same argument: a vertex that is wrongly kept
% or wrongly dropped changes the mesh's TOPOLOGY, never any face's VALUE, because `eval` reads the
% face functions and their sign conditions and neither mentions a vertex. The exact test would be
% the sign of a rational polynomial at a degree-16 algebraic number, which `ratQ.signAt` can do for
% a conic row but not for a BRANCH row -- that one is a comparison against a root, not a sign.

    if nargin ~= 1
        error('algCrvVertices:nargin', 'expects (g).');
    end
    if ~isa(g, 'AlgCrv')
        error('algCrvVertices:type', 'expects an AlgCrv; got a %s.', class(g));
    end
    t0 = tic;
    Vname = zeros(0,5);
    info = struct('candidates', 0, 'named', 0, 'refused', 0, 'seconds', 0, 'why', {{}});

    curves = curveTable(g);                      % one polyN per (table, index) actually used
    seen = containers.Map('KeyType', 'char', 'ValueType', 'logical');

    for k = 1:g.nf
        rows = g.FC{k};
        if size(rows,1) < 2, continue, end
        for a = 1:size(rows,1)-1
            for b = a+1:size(rows,1)
                [ia, ta] = rowEdge(rows, a);
                [ib, tb] = rowEdge(rows, b);
                if ia == ib && ta == tb, continue, end
                key = sprintf('%d,%d,%d,%d', ta, ia, tb, ib);
                if seen.isKey(key), continue, end
                seen(key) = true;

                P = curves{ta+1}{ia};  Q = curves{tb+1}{ib};
                if isempty(P) || isempty(Q), continue, end
                try
                    X = meetCurves(P, Q);
                catch ME
                    if ~startsWith(ME.identifier, 'ratQ:'), rethrow(ME), end
                    info.refused = info.refused + 1;
                    % The identifier AND the shape, because a refusal on a degree-1 curve and one
                    % on a degree-4 curve are different problems wearing the same name.
                    if info.refused <= 2
                        for si = 1:min(5, numel(ME.stack))
                            fprintf('   STACK %s line %d\n', ME.stack(si).name, ME.stack(si).line);
                        end
                    end
                    if info.refused <= 6
                        fprintf('REFUSAL %s | deg %dx%d and %dx%d | tables %d,%d\n', ...
                            ME.identifier, size(P,1)-1, size(P,2)-1, size(Q,1)-1, size(Q,2)-1, ...
                            ta, tb);
                    end
                    info.why{end+1} = sprintf('%s | deg %dx%d and %dx%d | tables %d,%d', ...
                        ME.identifier, size(P,1)-1, size(P,2)-1, size(Q,1)-1, size(Q,2)-1, ...
                        ta, tb); %#ok<AGROW>
                    continue
                end
                info.candidates = info.candidates + size(X,1);
                for r = 1:size(X,1)
                    if ~onFace(g, k, X(r,:)), continue, end
                    Vname(end+1,:) = [ia, ta, ib, tb, r]; %#ok<AGROW>
                    info.named = info.named + 1;
                end
            end
        end
    end
    info.seconds = toc(t0);
end

% ==================================================================================================
function X = meetCurves(P, Q)
% The real intersections of two integer curves, in `conicMeet`'s canonical order -- lexicographic by
% (x, y).
%
% WHICH VARIABLE IS ELIMINATED IS NOT FREE. A resultant needs BOTH arguments to have positive degree
% in the variable it removes, and a VERTICAL LINE has no y at all -- so eliminating y refuses it
% (`ratQ:resultantDegree`, measured on 66 of 176 conic-by-curve pairs before this chose). Such a
% pair is perfectly ordinary and has a vertex; it just has to be taken the other way round. So the
% variable is chosen, and when x is the one eliminated the roots are y-coordinates and the fibre is
% read in the other direction.
    X = zeros(0,2);
    dyP = size(P,2) - 1;  dyQ = size(Q,2) - 1;
    dxP = size(P,1) - 1;  dxQ = size(Q,1) - 1;
    if dyP >= 1 && dyQ >= 1
        k = 2;
    elseif dxP >= 1 && dxQ >= 1
        k = 1;
    else
        return                                   % one of them is a constant: no curve, no vertex
    end

    Rb = ratQ.resultantVarBig(P, Q, k, 2);
    R = bigToDoubles(Rb);
    R = R(:).';
    if isempty(R) || all(R == 0)
        return                                   % a shared component: no isolated vertex
    end
    ts = realRoots(R);
    for i = 1:numel(ts)
        if k == 2
            others = fibreRoots(P, ts(i));       % eliminated y: ts are x, fibre gives y
            if isempty(others), others = fibreRoots(Q, ts(i)); end
            for o = others, X = keepIfOn(X, [ts(i), o], P, Q); end
        else
            others = fibreRootsX(P, ts(i));      % eliminated x: ts are y, fibre gives x
            if isempty(others), others = fibreRootsX(Q, ts(i)); end
            for o = others, X = keepIfOn(X, [o, ts(i)], P, Q); end
        end
    end
    if ~isempty(X), X = unique(round(sortrows(X), 10), 'rows', 'stable'); end
end

function X = keepIfOn(X, pt, P, Q)
% A candidate is kept only if BOTH curves vanish there -- the resultant's roots include the
% projections of complex intersections, which this is what removes.
    if abs(evalPoly2(P, pt)) < 1e-7 * scaleOf(P) && abs(evalPoly2(Q, pt)) < 1e-7 * scaleOf(Q)
        X(end+1,:) = pt; %#ok<AGROW>
    end
end

function xs = fibreRootsX(P, y)
% The real x with P(x, y) = 0 -- the fibre over one y, the counterpart of fibreRoots.
    c = zeros(1, size(P,1));
    for i = 1:size(P,1)
        c(i) = polyval(flip(P(i,:)), y);
    end
    c = flip(c);
    nz = find(abs(c) > 1e-12, 1);
    if isempty(nz), xs = []; return, end
    rr = roots(c(nz:end));
    xs = real(rr(abs(imag(rr)) < 1e-8)).';
end

function R = bigToDoubles(Rb)
% The resultant's coefficients as doubles, EXACTLY when they fit and rounded when they do not.
%
% Two curves of degree 4 meet in a resultant of degree up to 16, and measured on the real object its
% coefficients pass 2^53 for 49 of the pairs -- `ratQ.fromBig` raises `ratQ:overflow` there, which is
% the honest answer to "is this an exact integer array" and the wrong question here.
%
% WHAT IS AND IS NOT LOST. The resultant itself is computed in BigInteger and is exact; what this
% rounds is the copy used to LOCATE its roots. The vertex's NAME is [edges..., rootIdx] and its
% meaning is "the rootIdx-th real intersection in the canonical order" -- an exact statement about
% an exact polynomial, with no coordinate stored. `conicMeet` computes its own roots numerically for
% the same reason one level down, so this is the established standard and not a new compromise. What
% could go wrong is the ORDER of two roots closer together than the rounding, which would rename a
% vertex rather than invent one.
    R = zeros(Rb.sz);
    over = false;
    for i = 1:numel(Rb.c)
        v = Rb.c{i};
        if isempty(v) || v.signum() == 0, continue, end
        if v.abs().bitLength() > 52, over = true; end
        R(i) = v.doubleValue();
    end
    if ~over
        R = ratQ.trimN(R);                       % it fits: the array IS exact
    else
        R = trimApprox(R);
    end
end

function A = trimApprox(A)
% `ratQ.trimN` insists on exact integers, so a rounded array is trimmed here instead.
    while numel(A) > 1 && A(end) == 0, A(end) = []; end
end

function r = realRoots(P)
% The real roots of a univariate polynomial, ascending.
%
% EXACT ISOLATION WHEN THE COEFFICIENTS FIT, and numeric roots when they do not. `ratQ.isolateRoots`
% is the better route -- Sturm gives intervals that provably each hold one root -- but it goes
% through `ratQ.chk`, so a coefficient past 2^53 raises `ratQ:overflow`. That is not a rare corner
% here: a resultant of two degree-4 curves has degree up to 16, and 49 of the real object's pairs
% land past the ceiling. Refusing them would leave the vertex table silently short, so those go to
% `roots`, which is what `conicMeet` uses for its own quartic anyway.
    P = P(:).';
    if isempty(P) || all(P == 0), r = []; return, end
    if max(abs(P)) < 2^53 && all(P == round(P))
        % THE POLYNOMIAL FITTING IS NOT ENOUGH, and that was measured rather than assumed. Sturm's
        % chain is built from PSEUDO-REMAINDERS, whose coefficients grow far past the polynomial's
        % own: on the degree-16 resultant of two quartic edges the chain overflows 2^53 although
        % every coefficient of the resultant is inside it. The stack said so --
        % isolateRoots -> sturmChain -> polyPseudoRem -> ratQ:overflow -- after two wrong guesses
        % about where the ceiling was being hit.
        try
            I = ratQ.isolateRoots(ratQ.polyTrim(P));
            r = zeros(1, size(I,1));
            for i = 1:size(I,1)
                lo = I(i,1);  hi = I(i,2);
                for it = 1:80
                    mid = (lo + hi) / 2;
                    if polyval(P, lo) * polyval(P, mid) <= 0, hi = mid; else, lo = mid; end
                end
                r(i) = (lo + hi) / 2;
            end
            return
        catch ME
            % `isolateNoConverge` is the same message in the other voice: the exact route did not
            % finish on this polynomial, whether because the chain overflowed or because the
            % bisection bound ran out. Both mean "use the numeric route", which is what `conicMeet`
            % uses for its own quartic.
            if ~any(strcmp(ME.identifier, {'ratQ:overflow', 'ratQ:isolateNoConverge'}))
                rethrow(ME)
            end
        end
    end
    nz = find(abs(P) > 0, 1);
    rr = roots(P(nz:end));
    sc = max(1, max(abs(P)))^(1/max(1, numel(P)-1));
    r = sort(real(rr(abs(imag(rr)) < 1e-8 * max(1, sc))).');
end

function ys = fibreRoots(P, x)
% The real y with P(x, y) = 0, in doubles -- the fibre over one x.
    c = zeros(1, size(P,2));
    for j = 1:size(P,2)
        c(j) = polyval(flip(P(:,j)).', x);
    end
    c = flip(c);
    nz = find(abs(c) > 1e-12, 1);
    if isempty(nz), ys = []; return, end
    rr = roots(c(nz:end));
    ys = real(rr(abs(imag(rr)) < 1e-8)).';
end

function v = evalPoly2(P, x)
    v = 0;
    for i = 1:size(P,1)
        for j = 1:size(P,2)
            if P(i,j) ~= 0, v = v + P(i,j) * x(1)^(i-1) * x(2)^(j-1); end
        end
    end
end

function s = scaleOf(P)
    s = max(1, max(abs(P(:))));
end

% --------------------------------------------------------------------------------------------------
function tf = onFace(g, k, x)
% Does this point satisfy face k's own rows? In doubles, for the reason the header gives.
    tf = true;
    rows = g.FC{k};
    for t = 1:size(rows,1)
        [idx, tbl] = rowEdge(rows, t);
        if tbl == 1
            a = g.EcA(idx);
            v = algBranch(a.FP, x, a.FB);
            if ~isfinite(v), tf = false; return, end
            q = ratQ.evalFace(a.num, a.den, x);
            if rows(t,2) * (q - v) < -1e-7, tf = false; return, end
        else
            c = g.EcQ(idx, :);
            br = 0;
            if size(rows,2) >= 3, br = rows(t,3); end
            if br ~= 0
                v = ratQ.pairLineEval(c, br, x);
            else
                v = c(1)*x(1)^2 + c(2)*x(1)*x(2) + c(3)*x(2)^2 + c(4)*x(1) + c(5)*x(2) + c(6);
            end
            if rows(t,2) ~= 0 && rows(t,2) * v < -1e-7, tf = false; return, end
        end
    end
end

function [idx, tbl] = rowEdge(rows, t)
    idx = rows(t,1);
    tbl = 0;
    if size(rows,2) >= 4, tbl = rows(t,4); end
end

% --------------------------------------------------------------------------------------------------
function C = curveTable(g)
% Every edge as a polyN in (x, y): the conics directly, the algebraic edges through
% `bridgeQuadCut`, which is the rational curve containing the boundary they name.
    C = {cell(1, size(g.EcQ,1)), cell(1, numel(g.EcA))};
    for i = 1:size(g.EcQ,1)
        c = g.EcQ(i,:);
        M = zeros(3,3);
        M(3,1) = c(1);  M(2,2) = c(2);  M(1,3) = c(3);
        M(2,1) = c(4);  M(1,2) = c(5);  M(1,1) = c(6);
        C{1}{i} = ratQ.trimN(M);
    end
    for j = 1:numel(g.EcA)
        a = g.EcA(j);
        C{2}{j} = bridgeQuadCut(a.FP, a.num, a.den);
    end
end
