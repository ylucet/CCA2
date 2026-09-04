classdef QuaCon < RatCon & Qua
% QuaCon  A quadratic function on a CONIC subdivision, stored EXACTLY: rational faces, primitive
% integer edge conics, and vertices that are NAMES rather than coordinates.
%
% WHAT IT IS FOR. This is the type `conj` produces once it computes exactly. `CONJ_FIELD_PROOF.md`
% Theorem 1 says precisely where the irrationality of `f*` is and is not: **the face functions and
% the edge conics are always rational, and only the VERTICES leave Q** (degree <= 4, S4 in the
% worst case). So the representation splits along that line -- everything rational is stored
% exactly, and the one thing that is not is never stored at all. A vertex is the label "the
% intersection of edge conic i and edge conic j, root r"; its coordinates are computed on demand
% by `conicMeet` and thrown away again.
%
% WHY IT IS NOT A QuaPar. A conjugate's edge between faces g_i and g_j lies on {g_i = g_j}, whose
% discriminant is -det(H_i - H_j). [COAP]/[JOGO] Theorem 6 makes that vanish for two ADJACENT
% pieces of f, but Step 3 maximises over ALL pairs, so from three pieces up some compared pair is
% non-adjacent and the edge is a genuine ELLIPSE. `CONJ_FIELD_PROOF.md` 7.4b traces such an arc of
% positive length against the exact definition of the conjugate, every traced point on the conic to
% <= 2e-15; `doc/QuaConExample.md` 2 reduces it to three pieces. `QuaPar.assertParabolicEdges`
% correctly refuses to store that as a parabola, and this class is where it goes instead.
%
% ------------------------------------------------------------------------------------------------
% THE FOUR THINGS STORED, AND WHY EACH IS THE SHAPE IT IS
%
%   fN, fD    the face functions, exactly. Row k is a numerator in RatCon's 10-wide weighted cubic
%             basis over the positive integer fD(k), in `ratQ.canon` form. A face HAS a scale
%             (doubling it is a different function), so this is the canonical form with the gcd
%             removed and the denominator positive.
%
%   EcQ       the edge conics, exactly. Row j is `ratQ.conic` form -- PRIMITIVE, first nonzero
%             entry positive. A conic has NO scale: {c'm = 0} and {1000 c'm = 0} are one curve, so
%             two edges lie on the same curve exactly when these two integer vectors are BITWISE
%             EQUAL. That is the single property this whole design is built to obtain, and it is
%             the one neither doubles nor `sym` can provide: `DECISIONS.md` 2026-08-17 records an
%             ULP hiding a shared facet from `merge`, after which Step 3's cell count grew without
%             bound, and `sym` has no canonical form at all (equal quantities can compare Unknown).
%             EVERY edge carries its own equation, including a straight one, which is [0 0 0 d e f].
%
%   Vname     the vertices, as NAMES. Row i is [edgeA edgeB rootIdx]: the point where the conics of
%             those two edges meet, rootIdx selecting among the (at most four) real intersections in
%             `conicMeet`'s canonical order. No coordinate is stored, ever -- see `vertexCoords`.
%
%   FC        the faces, in H-FORM. FC{k} is an m x 2 list of [edgeIndex, side], meaning
%                     face k = { x : side_i * Q_{e_i}(x) >= 0 for every row i },
%             and a side of 0 means ON the curve, i.e. Q = 0 -- an EQUALITY rather than an
%             inequality. That third value is what lets a face be THIN: dom f* is a line when the
%             primal is a positive semidefinite singular quadratic on the whole plane, and a single
%             point when the primal is affine there. Both are correct conjugates, and without the
%             equality side there is nowhere to put either. It costs one value in a column that
%             already existed rather than a new class -- see TODO.md 2026-09-04.
%             This is NOT derivable from the combinatorics. For a straight edge the side can be
%             recovered from the face's other vertices; for a conic it cannot -- the two sides of an
%             ellipse are not distinguished by the edge's endpoints -- so the sign is stored.
%             `CONJ_FIELD_PROOF.md` 8.6 calls exactly this the H-representation.
%
% WHAT IS INHERITED AND USED. The COMBINATORIAL half of RatCon's mesh carries over untouched,
% because none of it mentions a coordinate: `nv`/`ne`/`nf`, `E` (edge endpoints as vertex INDICES,
% plus the segment/ray flag), `F` (edge -> the face on each side), `P` (the face adjacency lists),
% `dom` and `fIsConvex`.
%
% WHAT IS INHERITED AND DELIBERATELY LEFT EMPTY: `V`, `f` and `Ec`. Those are RatCon's DOUBLE
% coordinate/coefficient fields and a QuaCon has no doubles to put in them -- `V` cannot be filled
% without storing an irrational, and filling `f`/`Ec` would be a second, rounded copy of data that
% is already here exactly, i.e. two sources of truth for one value. That is the same objection
% RatCon.m records against a stored type flag: a copy that can disagree with the object is worse
% than no copy. Reading `.f`, `.V` or `.Ec` off a QuaCon therefore yields [] BY DESIGN, and the
% legacy double routines are not meant to consume one; use `fN`/`fD`, `vertexCoords()` and `EcQ`.
% `den` IS filled, with the honest value: a QuaCon is a Qua, so its denominator is identically 1.
%
% NOTHING HERE IS APPROXIMATE, WITH ONE NAMED EXCEPTION. Every stored number is an exact integer and
% every mesh-deciding comparison is exact. `vertexCoords` and `eval` return doubles on purpose:
% `CONJ_FIELD_PROOF.md` 8.0 lists the operations that must leave Q -- building a vertex, the sign of
% a rational polynomial at one, and comparing two of them -- and EVALUATING the finished function at
% a point is not among them. `ratQ.evalFace` says the same thing. So point location inside `eval`
% screens with a tolerance; the exact kernel behind that screen is Phase 2c and is NOT yet built,
% which is why `eval` documents its own limit rather than pretending otherwise.

   properties
        % See the header for why each is shaped as it is. All four are EXACT; none is a cache.
        fN (:,10){mustBeNumeric}   % nf x 10 integer numerators, ratQ.canon form against fD
        fD (:,1) {mustBeNumeric}   % nf x 1  positive integer denominators
        EcQ (:,6){mustBeNumeric}   % ne x 6  primitive integer conics, ratQ.conic form
        Vname (:,3){mustBeInteger,mustBeNonnegative} % nv x 3 [edgeA edgeB rootIdx]
        FC                         % cell(nf,1); FC{k} is m x 2 of [edgeIndex, sign in {-1,+1}]
   end

   methods
       function obj = QuaCon(Vname, E, EcQ, fN, fD, F, FC)
       % objective: build a QuaCon from exact data, validating every invariant the rest of the
       %            class then relies on.
       % [input]  Vname : nv x 3 [edgeA edgeB rootIdx]   E   : ne x 3 [v1 v2 isSegment]
       %          EcQ   : ne x 6 integer conics          fN  : nf x 10 integer numerators
       %          fD    : nf x 1 positive integers       F   : ne x 2 faces either side of an edge
       %          FC    : cell(nf,1) of m x 2 [edgeIndex, sign]
       %
       % No-argument path writes NOTHING and returns -- RatCon.m's CONSTRUCTOR PROTOCOL. MATLAB
       % re-runs a shared base constructor once per inheritance path and the second run overwrites
       % the first, so all state-writing happens here, in the leaf, and nowhere above.
       %
       % THE VALIDATION IS THE POINT. Every downstream routine assumes canonical forms and
       % in-range names; a QuaCon that lies about either is a silent wrong answer rather than a
       % crash. Canonicalisation is APPLIED rather than merely checked (that is what a canonical
       % form is for), but a name that cannot be realised RAISES -- there is no canonical repair
       % for "the third intersection of two conics that meet twice".
            if nargin == 0, return; end
            if nargin ~= 7
                error('QuaCon:nargin', ...
                    'QuaCon needs 7 arguments (Vname, E, EcQ, fN, fD, F, FC); got %d.', nargin);
            end

            nv = size(Vname,1);  ne = size(E,1);  nf = size(fN,1);
            if size(E,2) ~= 3 || size(EcQ,1) ~= ne || size(F,1) ~= ne || size(F,2) ~= 2
                error('QuaCon:size', ...
                    ['E must be ne x 3 and EcQ / F must have one row per edge; got E %dx%d, ' ...
                     'EcQ %dx%d, F %dx%d.'], size(E,1), size(E,2), size(EcQ,1), size(EcQ,2), ...
                     size(F,1), size(F,2));
            end
            if numel(fD) ~= nf || ~iscell(FC) || numel(FC) ~= nf
                error('QuaCon:size', ...
                    'fD and FC need one entry per face: %d faces, %d denominators, %d constraint lists.', ...
                    nf, numel(fD), numel(FC));
            end
            if nv > 0 && max(E(:,1:2), [], 'all') > nv
                error('QuaCon:badEdge', 'E names vertex %d but there are only %d vertices.', ...
                    max(E(:,1:2), [], 'all'), nv);
            end

            % ---- the exact face functions, canonically -------------------------------------
            for k = 1:nf
                [fN(k,:), fD(k)] = ratQ.canon(fN(k,:), fD(k));
                if any(fN(k,1:4) ~= 0)
                    error('QuaCon:cubicFace', ...
                        ['face %d has a cubic numerator. A QuaCon holds the QUADRATIC output of ' ...
                         'conj; a cubic has no conic level set, so its edges could not be ' ...
                         'stored here at all (see ratQ.diffConic).'], k);
                end
            end

            % ---- the exact edge conics, canonically ----------------------------------------
            for j = 1:ne
                EcQ(j,:) = ratQ.conic(EcQ(j,:));   % raises ratQ:zeroConic on an all-zero row
            end

            % ---- the vertex NAMES, and that each one can actually be realised ---------------
            for i = 1:nv
                a = Vname(i,1);  b = Vname(i,2);  r = Vname(i,3);
                if a < 1 || b < 1 || a > ne || b > ne
                    error('QuaCon:badVertexName', ...
                        'vertex %d names edges %d and %d, but there are %d edges.', i, a, b, ne);
                end
                if a == b
                    error('QuaCon:badVertexName', ...
                        ['vertex %d names edge %d twice. A vertex is where TWO curves meet; one ' ...
                         'curve does not name a point.'], i, a);
                end
                [P, info] = conicMeet(EcQ(a,:), EcQ(b,:));
                if info.degenerate
                    error('QuaCon:degenerateVertex', ...
                        ['vertex %d names edges %d and %d, whose conics share a component, so ' ...
                         'they meet in infinitely many points and name no vertex.'], i, a, b);
                end
                if r < 1 || r > size(P,1)
                    error('QuaCon:badRootIndex', ...
                        ['vertex %d asks for real intersection %d of edges %d and %d, which meet ' ...
                         'in %d real points.'], i, r, a, b, size(P,1));
                end
            end

            % ---- the H-form face constraints ------------------------------------------------
            for k = 1:nf
                Ck = FC{k};
                if isempty(Ck), Ck = zeros(0,2); end
                if size(Ck,2) ~= 2
                    error('QuaCon:badConstraint', ...
                        'FC{%d} must be m x 2 of [edgeIndex, sign]; got %d columns.', k, size(Ck,2));
                end
                if ~isempty(Ck)
                    if any(Ck(:,1) < 1) || any(Ck(:,1) > ne)
                        error('QuaCon:badConstraint', ...
                            'FC{%d} names an edge outside 1..%d.', k, ne);
                    end
                    if ~all(ismember(Ck(:,2), [-1 0 1]))
                        error('QuaCon:badConstraint', ...
                            ['FC{%d} carries a side that is not -1, 0 or +1. The side says WHICH ' ...
                             'side of the conic the face lies on and cannot be omitted; 0 means ' ...
                             'ON the curve.'], k);
                    end
                end
                FC{k} = Ck;
            end

            obj.nv = nv;  obj.ne = ne;  obj.nf = nf;
            obj.E = E;    obj.F = F;
            obj.fN = fN;  obj.fD = fD(:);  obj.EcQ = EcQ;  obj.Vname = Vname;  obj.FC = FC;
            % A QuaCon IS a Qua, so the denominator is identically 1 -- the honest value, not a
            % placeholder. `f`, `V` and `Ec` stay EMPTY; see the header for why that is by design.
            obj.den = repmat([0 0 1], nf, 1);
       end

       function tf = isMeshed(obj) %#ok<MANU>
       % objective: whether this object carries an explicit mesh. It does.
       %
       % OVERRIDDEN because RatCon's version tests `~isempty(obj.f)`, and a QuaCon's `f` is empty
       % by design (the faces live in fN/fD). The inherited test would report a fully meshed object
       % as unmeshed, which is the opposite of the truth and would send callers down
       % QuaParCPLQ's "mesh not yet reconstructed" path.
            tf = true;
       end

       function V = vertexCoords(obj)
       % objective: realise the named vertices as double coordinates, ON DEMAND.
       % [output] V : nv x 2 double
       %
       % NOT STORED, NOT CACHED, AND THAT IS THE DESIGN. A vertex of f* has degree up to 4 over Q
       % and is generically irrational (of twelve feasible continuous three-piece configurations
       % the vertex quartic is irreducible over Q in TEN -- CONJ_FIELD_PROOF.md 8.0), so any stored
       % coordinate would be a rounded copy of the truth, and the mesh would then be deciding its
       % own combinatorics from rounded data. That is the failure this class exists to remove.
       % Coordinates are an OUTPUT convenience -- plotting, evaluation, a user's eye -- and the
       % caller who wants them repeatedly should hold the result, not ask the object to.
            V = zeros(obj.nv, 2);
            for i = 1:obj.nv
                P = conicMeet(obj.EcQ(obj.Vname(i,1),:), obj.EcQ(obj.Vname(i,2),:));
                V(i,:) = P(obj.Vname(i,3), :);
            end
       end

       function tf = sameEdgeCurve(obj, j1, j2)
       % objective: do edges j1 and j2 lie on the SAME curve. Exact; no tolerance anywhere.
       %
       % One line, and it is the payoff of the whole representation: because EcQ rows are primitive
       % and sign-normalised, "same curve" is bitwise equality of two integer vectors. Two arcs of
       % one ellipse are recognised as such however they were built, which is what `merge` needs
       % and what an ULP or a symbolic expression cannot deliver.
            tf = isequal(obj.EcQ(j1,:), obj.EcQ(j2,:));
       end

       function k = edgeKind(obj, j)
       % objective: name edge j's conic type from its EXACT discriminant.
       % [output] k : 'line' | 'parabola' | 'hyperbolic' | 'elliptic'
            k = ratQ.conicKind(obj.EcQ(j,:));
       end

       function [val, idx] = eval(obj, X)
       % objective: the value of the function at points X, and which face supplied it.
       % [input]  X   : k x 2 double
       % [output] val : k x 1 double, +Inf off the domain; idx : k x 1 face index, 0 where +Inf
       %
       % EVALUATION IS DELIBERATELY NOT EXACT, and this is the one place in the class where a
       % tolerance appears. CONJ_FIELD_PROOF.md 8.0 lists the operations that must leave Q --
       % building a vertex, the sign of a rational polynomial at one, comparing two of them -- and
       % evaluating the finished function at a point is not among them; ratQ.evalFace records the
       % same split. The mesh is decided exactly, then evaluated in doubles.
       %
       % WHAT IS NOT YET BUILT, STATED PLAINLY. The sign test below is the FILTER of
       % CONJ_FIELD_PROOF.md 8.7 without its exact fallback: a point within `tol` of a face
       % boundary is accepted, and a point that straddles two faces' tolerances takes the first.
       % On the boundary itself both faces give the same value (f* is continuous), so the answer is
       % right there; the case this cannot yet decide is a point whose membership is genuinely
       % undecided at double precision, which is what Phase 2c's exact kernel is for.
            if isempty(X), val = zeros(0,1); idx = zeros(0,1); return, end
            n = size(X,1);
            val = inf(n,1);  idx = zeros(n,1);
            for k = 1:obj.nf
                Ck = obj.FC{k};
                inFace = true(n,1);
                for r = 1:size(Ck,1)
                    q = QuaPar.evalConic(obj.EcQ(Ck(r,1),:), X);
                    scale = max(1, max(abs(obj.EcQ(Ck(r,1),:))) * max(1, max(abs(X(:)))^2));
                    if Ck(r,2) == 0
                        % ON the curve. A THIN face, and the honest caveat is that a point picked
                        % at random in floating point is never exactly on a line, so such a face
                        % answers +inf almost surely -- the mesh is still exactly right, and the
                        % exact predicates that build it are unaffected. See this class's header.
                        inFace = inFace & (abs(q) <= sqrt(eps) * scale);
                    else
                        inFace = inFace & (Ck(r,2) * q >= -sqrt(eps) * scale);
                    end
                end
                take = inFace & (idx == 0);
                if any(take)
                    val(take) = ratQ.evalFace(obj.fN(k,:), obj.fD(k), X(take,:));
                    idx(take) = k;
                end
            end
       end

       function disp(obj)
       % objective: a readable summary that never pretends a vertex is a number.
            if isempty(obj.fN)
                fprintf('  QuaCon (empty)\n');  return
            end
            fprintf('  QuaCon: %d vertices (named), %d edges, %d faces\n', obj.nv, obj.ne, obj.nf);
            kinds = cell(obj.ne,1);
            for j = 1:obj.ne, kinds{j} = obj.edgeKind(j); end
            u = unique(kinds);
            for i = 1:numel(u)
                fprintf('    %-10s %d edge(s)\n', u{i}, sum(strcmp(kinds, u{i})));
            end
            for k = 1:min(obj.nf, 8)
                fprintf('    face %d: [%s] / %d\n', k, ...
                    strtrim(sprintf('%d ', obj.fN(k,5:10))), obj.fD(k));
            end
            if obj.nf > 8, fprintf('    ... %d more faces\n', obj.nf - 8); end
       end
   end
end
