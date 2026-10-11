function [rows, info] = cornerWedgeRows(f, n, d, tie, side)
% cornerWedgeRows  A corner face's WEDGE as [conic, side, branch] rows, so that a face whose region
%   is a predicate can enter an arrangement after all.
%
% objective: `envelopeMesh` reads a cell's applicable faces off its signs, which needs every face's
%   region spelled as rows over the INTEGERS. A corner's region is `M \ (p - x) >= 0`, two half-planes
%   whose coefficients are quadratic irrationals -- so no row of it exists over Q on its own. Its
%   CONJUGATE PAIR does, and `ratQ.pairLineData` names one line of such a pair exactly.
%
% [input]  f    : one face of `cornerFaceAlg`, with .sigma, .M, .p, .s
%          n, d : the cell's face, as `cornerFaceAlg` took it -- 1x10 numerator over d > 0
%          tie  : the tie conic, integer [a b c d e f]
%          side : +1 or -1, as passed to `cornerFaceAlg`
% [output] rows : k x 8, each [conic(1:6), side, branch] with the same meaning as `QuaCon.FC` --
%                 x is in the wedge iff row(7)*pairLineEval(row(1:6), row(8), x) >= 0 for every row.
%                 A branch of 0 means the row is an ordinary rational line, carried in the same
%                 columns with no pair around it.
%          info : .witness, .kind per row, .refused (empty when rows were produced)
%
% ==================================================================================================
% THE ELIMINATION, AND WHY IT IS DONE THIS WAY
%
% An earlier attempt reduced the multipliers modulo the facet's parameter quadratic BY HAND, and was
% measured wrong: 450 of 3200 points disagreed with the predicate and half the corners named no
% branch at all (DECISIONS.md 2026-09-18). The defect was arithmetic, not geometry, so the arithmetic
% is now the thing that is written down ONCE:
%
%       every quantity here lives in Q(sigma) = { (u1 + u2 sigma)/D },
%
% with `sigma` a root of `al sigma^2 + be sigma + ga`, and the local `qmul` is the only place the
% relation `al sigma^2 = -be sigma - ga` is ever used. The corner, the tie gradient there, `M`, its
% adjugate and `p` are then built by COMPOSITION, and no index can be dropped in a step that no
% longer exists.
%
% THE ROW. Membership is `mu = M \ (p - x) >= 0`, i.e. `sgn(det M) * (adj(M)(p - x))_i >= 0` -- so
% each edge is the affine form
%
%       L_i(x) = A_i(x) + sigma B_i(x) ,   A_i, B_i affine over Z ,
%
% and the two conjugate corners give the two lines `L_i^{(sigma)}` and `L_i^{(sigma')}`. Their
% product is SYMMETRIC in the two roots, so it is rational, and clearing the denominator with `al`:
%
%       al * L^{(sigma)} L^{(sigma')}  =  al A^2 - be A B + ga B^2 ,
%
% since `sigma + sigma' = -be/al` and `sigma sigma' = ga/al`. That is the pair, over the integers, in
% one line and with nothing to reduce.
%
% B == 0 IS NOT A FAILURE. It means the edge is already RATIONAL -- which the FACET edge can be, its
% normal being rational -- and then the pair above is a perfect square, i.e. one line counted twice.
% It is emitted as itself, with branch 0, rather than handed to `pairLineData`, which rightly refuses
% a doubled line. The other half of the same fact is that the facet edge's pair, when it is not
% rational, is a pair of PARALLEL lines: `ratQ.parallelPairData` covers exactly that case.
%
% THE SIDE AND THE BRANCH ARE READ AT A WITNESS, and the witness needs no search: `mu = [1; 1]` is
% strictly interior, so `x_w = p - M [1; 1]` is an interior point by construction.

    if nargin ~= 5
        error('cornerWedgeRows:nargin', 'expects (f, n, d, tie, side).');
    end
    if ~isstruct(f) || ~all(isfield(f, {'sigma', 'M', 'p', 's'}))
        error('cornerWedgeRows:shape', 'f must be a face of cornerFaceAlg.');
    end
    tie = ratQ.chk(tie(:).', 'tie conic');
    n = ratQ.chk(double(n(:).'), 'face');
    N = [n(5) n(6); n(6) n(7)];  m = [n(8); n(9)];
    al = f.sigma.quad(1);  be = f.sigma.quad(2);  ga = f.sigma.quad(3);
    Q  = [al be ga];

    % ---- the corner, exactly: s = (s0n + sigma w)/na ---------------------------------------------
    na = f.sigma.na;
    s1 = qred([f.sigma.s0n(1), f.sigma.w(1), na]);
    s2 = qred([f.sigma.s0n(2), f.sigma.w(2), na]);

    % ---- the tie gradient there, and M = [a, side grad q] ----------------------------------------
    g1 = qadd(qadd(qmul([2*tie(1) 0 1], s1, Q), qmul([tie(2) 0 1], s2, Q)), [tie(4) 0 1]);
    g2 = qadd(qadd(qmul([tie(2) 0 1], s1, Q), qmul([2*tie(3) 0 1], s2, Q)), [tie(5) 0 1]);
    a  = f.M(:,1);                              % the facet normal, RATIONAL and already integral
    M  = {[a(1) 0 1], qmul([side 0 1], g1, Q); ...
          [a(2) 0 1], qmul([side 0 1], g2, Q)};

    % ---- p = grad g at the corner, EXACTLY ---------------------------------------------------------
    % `cornerFaceAlg` holds p only in doubles, so it is rebuilt here in Q(sigma) from the same data.
    p1 = qred(qadd(qadd(qmul([N(1,1) 0 1], s1, Q), qmul([N(1,2) 0 1], s2, Q)), [m(1) 0 1]));
    p2 = qred(qadd(qadd(qmul([N(2,1) 0 1], s1, Q), qmul([N(2,2) 0 1], s2, Q)), [m(2) 0 1]));
    p1(3) = p1(3) * d;   p2(3) = p2(3) * d;   p1 = qred(p1);  p2 = qred(p2);

    % ---- adj(M) (p - x), one affine form in Q(sigma) per edge -------------------------------------
    adjM = {M{2,2}, qneg(M{1,2}); qneg(M{2,1}), M{1,1}};
    rows = zeros(0,8);
    info = struct('witness', [], 'kind', {{}}, 'refused', '');

    sg = (-be + f.sigma.root*sqrt(be^2 - 4*al*ga)) / (2*al);   % OUR root, in doubles
    xw = f.p - f.M * [1; 1];                    % strictly interior: mu = [1; 1] > 0 there
    info.witness = xw(:).';
    info.pExact = [qval(p1, sg), qval(p2, sg)];        % the corner's gradient, rebuilt
    info.sExact = [qval(s1, sg), qval(s2, sg)];

    for i = 1:2
        [A, B] = edgeForm(adjM(i,:), {p1, p2}, Q);
        [r, kind] = rowFor(A, B, Q, xw, sg);
        if isempty(r)
            rows = zeros(0,8);
            info.refused = kind;
            return
        end
        rows(end+1,:) = r; %#ok<AGROW>
        info.kind{end+1} = kind;
    end
end

% ==================================================================================================
function [A, B] = edgeForm(adjRow, p, Q)
% L(x) = adjRow(1)*(p1 - x) + adjRow(2)*(p2 - y), as A(x) + sigma B(x) over Z.
% The x-coefficients are MINUS the adjugate entries, and the constant is the adjugate applied to p.
    c1 = adjRow{1};  c2 = adjRow{2};
    k  = qred(qadd(qmul(c1, p{1}, Q), qmul(c2, p{2}, Q)));       % the constant term
    D  = lcm2(lcm2(c1(3), c2(3)), k(3));
    A  = [-c1(1)*(D/c1(3)), -c2(1)*(D/c2(3)), k(1)*(D/k(3))];
    B  = [-c1(2)*(D/c1(3)), -c2(2)*(D/c2(3)), k(2)*(D/k(3))];
    g  = gcdv([A B]);
    if g > 0, A = A/g;  B = B/g; end
    A = ratQ.chk(A, 'wedge edge form');  B = ratQ.chk(B, 'wedge edge form');
end

% --------------------------------------------------------------------------------------------------
function [row, kind] = rowFor(A, B, Q, xw, sg)
% The edge as a row: its conjugate pair, which branch of it, and which side the witness is on.
    row = [];
    if all(B == 0)
        % RATIONAL edge -- the pair would be this line squared, which names no branch.
        g = gcdv(A);
        if g == 0, kind = 'degenerate';  return, end
        L = A / g;
        v = L(1)*xw(1) + L(2)*xw(2) + L(3);
        if abs(v) <= 1e-9 * max(1, max(abs(L)))
            kind = 'witnessOnEdge';  return
        end
        row = [0 0 0 L, sign(v), 0];
        kind = 'rational';
        return
    end
    if all(ratQ.chkB([A(1)*B(2) - A(2)*B(1), A(1)*B(3) - A(3)*B(1), A(2)*B(3) - A(3)*B(2)], ...
                     [abs(A(1)*B(2)) + abs(A(2)*B(1)), abs(A(1)*B(3)) + abs(A(3)*B(1)), ...
                      abs(A(2)*B(3)) + abs(A(3)*B(2))], 'wedge edge parallel test') == 0)
        % A PARALLEL to B: L = (lambda + sigma) B, the RATIONAL line B scaled by an irrational
        % number. Its conjugate product is a DOUBLE line, which no branch of a pair can name
        % (ratQ:notALinePair, 2026-10-11) -- but it needs none: the edge IS the line B, and the
        % witness, strictly inside the wedge, says which side.
        g = gcdv(B);
        L = B / g;
        v = L(1)*xw(1) + L(2)*xw(2) + L(3);
        if abs(v) <= 1e-9 * max(1, max(abs(L)))
            kind = 'witnessOnEdge';  return
        end
        row = [0 0 0 L, sign(v), 0];
        kind = 'rationalScaled';
        return
    end
    al = Q(1);  be = Q(2);  ga = Q(3);
    C = ratQ.chkB(al*linMul(A,A) - be*linMul(A,B) + ga*linMul(B,B), ...
                  abs(al)*linMul(abs(A),abs(A)) + abs(be)*linMul(abs(A),abs(B)) + abs(ga)*linMul(abs(B),abs(B)), ...
                  'wedge pair conic');
    if all(C == 0), kind = 'degenerate';  return, end
    C = ratQ.conic(C);
    try
        pd = ratQ.pairLineData(C);
    catch err
        kind = ['pair:' err.identifier];  return
    end
    [b, ok] = branchOf(C, A, B, sg);
    if ~ok, kind = 'branch';  return, end
    v = ratQ.pairLineEval(C, b, xw(:).');
    if abs(v) <= 1e-9 * max(1, max(abs(C)))
        kind = 'witnessOnEdge';  return
    end
    row = [C, sign(v), b];
    kind = pd.kind;
end

% --------------------------------------------------------------------------------------------------
function [b, ok] = branchOf(C, A, B, sg)
% WHICH line of the pair is ours: the one that VANISHES where ours does. Two points ON our own line
% decide it -- the two lines of the pair are distinct, the discriminant of the corner's quadratic
% being positive, so the other branch is bounded away from zero at those points.
    b = 0;  ok = false;
    L = A + sg*B;                               % our own edge, in doubles
    nrm = hypot(L(1), L(2));
    if nrm == 0, return, end
    x0 = -L(3) * [L(1), L(2)] / nrm^2;          % the foot of the perpendicular from the origin
    t  = [-L(2), L(1)] / nrm;                   % along the edge
    X  = [x0 + t; x0 - t; x0 + 3*t];
    sc = max(abs(C)) * max(1, max(abs(X(:))))^2;
    r  = [max(abs(ratQ.pairLineEval(C, 1, X))), max(abs(ratQ.pairLineEval(C, 2, X)))];
    [lo, b] = min(r);
    if lo > 1e-7 * sc || r(3-b) < 1e-4 * sc
        b = 0;  return                          % neither vanishes, or both do: not decidable here
    end
    ok = true;
end

function v = qval(u, sg)
    v = (u(1) + u(2)*sg) / u(3);
end

% ==================================================================================================
% Q(sigma) arithmetic: [u1 u2 D] means (u1 + u2 sigma)/D. It lives in `ratQ` now -- `qsRed`, `qsAdd`,
% `qsMul` -- because `cellInterior` needs the same field and one sign convention wrong in two places
% is the defect this file was rewritten for (DECISIONS.md 2026-09-18).
function u = qadd(u, v)
    u = ratQ.qsAdd(u, v);
end

function u = qneg(u)
    u = [-u(1), -u(2), u(3)];
end

function w = qmul(u, v, Q)
    w = ratQ.qsMul(u, v, Q);
end

function u = qred(u)
    u = ratQ.qsRed(u);
end

% --------------------------------------------------------------------------------------------------
function c = linMul(A, B)
% (A.x)(B.x) for two affine forms, as a conic [a b c d e f].
    c = [A(1)*B(1), A(1)*B(2) + A(2)*B(1), A(2)*B(2), ...
         A(1)*B(3) + A(3)*B(1), A(2)*B(3) + A(3)*B(2), A(3)*B(3)];
end

function g = gcdv(v)
    v = abs(v(v ~= 0));
    if isempty(v), g = 0; return, end
    g = v(1);
    for k = 2:numel(v), g = gcd(g, v(k)); end
end

function l = lcm2(a, b)
    l = abs(a*b) / gcd(abs(a), abs(b));
end
