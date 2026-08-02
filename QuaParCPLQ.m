classdef QuaParCPLQ < RatPar & Qua
    % QuaParCPLQ  Wraps a cPLQ `functionNDomain` array (a list of (region, symbolicFunction)
    % pairs) so it can be composed by the same generic convex-analysis operators (conj/add/
    % scalarMul/addQuadratic/addScaledEnergy/eval) that QuaPol/QuaPar already support -- the
    % "QuaPar-like return type" conjCPLQ's Case C needs (DESIGN.md II.5.1,
    % .claude/SESSION_HANDOFF.md "Next steps"). infConv/moreau/proxAverage/QuaPar.biconj call
    % these operators by plain function name (conj(f,engine), add(f,g), scalarMul(f,c), ...), so
    % MATLAB dispatches to whichever class the object at hand actually is -- providing the same
    % method names here is enough for a Case C result to compose, with NO changes needed to
    % infConv.m/moreau.m/proxAverage.m themselves.
    %
    % WHY NOT A REAL QuaPar: Case C's conjugate generally lives on an UNBOUNDED (all of R^2)
    % parabolic subdivision produced directly by cPLQ's own exact symbolic pipeline, not by
    % CCA2's V/E/Ec/F/P geometry -- reconstructing that geometry (vertex dedup, ray-vs-segment
    % edges, conic-edge orientation, face adjacency P) from the exact symbolic regions is a
    % separate, larger task, left open. QuaParCPLQ instead keeps the cPLQ representation as-is and
    % re-exposes only the operator surface conj/add/scalarMul/addQuadratic/eval need, driving the
    % SAME per-piece machinery `plq.biconjugateF` already uses (functionNDomain's own
    % conjugateOfPiecePoly/mergeL/addEq) -- generalized here into a standalone, recursively-
    % composable operator instead of being fused into the `plq` class.
    %
    % STATUS: `conj` reuses the exact conjugateOfPiecePoly/mergeL/addEq recipe biconjugateF calls;
    % that recipe has a known open bug for some domains (region.getNormalConeVertexQ -- see
    % conjCPLQ.m's own header) not fixed here. `add` overlays two full-plane piecewise functions
    % pairwise (same domain-intersection idiom as functionNDomain.times/mtimes) and sums their
    % values on each overlap cell -- correct whenever both operands' pieces jointly tile the SAME
    % region (true for two conjugates of compact-domain primal functions, where dom f*=R^2, the
    % expected use case for infConv/proxAverage). Mixed QuaParCPLQ/QuaPol/QuaPar composition
    % (e.g. infConv of a Case A/B function with a Case C one) is NOT supported -- add errors with
    % an unrecognized-class MATLAB error rather than silently doing the wrong thing.

    properties
        fnd = functionNDomain.empty();  % array of (region, symbolicFunction) pairs
    end

    methods
        function obj = QuaParCPLQ(fnd)
            if nargin == 0, return; end
            obj.fnd = fnd;
        end

        function assertOperable(obj) %#ok<MANU>
            % objective: interface parity with QuaPol/QuaPar.assertOperable -- a QuaParCPLQ is
            % always the quadratic-on-parabolic result of a 'cplq' conjugate, never cubic, so
            % there is nothing to reject here.
        end

        function [val, idx] = eval(obj, x)
            % objective: numerically evaluate at 1+ points x (kx2), mirroring QuaPar.eval's
            % [val,region] signature -- idx plays region's role (index into obj.fnd, 0 if none).
            val = zeros(size(x,1),1); idx = zeros(size(x,1),1);
            for i = 1:size(x,1)
                [val(i), idx(i)] = evalFunctionNDomain(obj.fnd, x(i,:));
            end
        end

        function h = conj(obj, engine)
            % objective: Legendre-Fenchel conjugate of this (already quadratic-on-parabolic)
            % piecewise function -- reuses plq.biconjugateF's own recipe verbatim.
            if nargin < 2, engine = 'cplq'; end
            if ~strcmpi(engine, 'cplq')
                error('QuaParCPLQ:conj:engine', ...
                    'QuaParCPLQ conjugate only supports the ''cplq'' engine (got ''%s'').', engine);
            end
            % Same shape as plq.biconjugateF: conjugate each piece, then take the MAX of
            % those conjugates -- the last step of the algorithm. This used to be
            % conjugateOfPiecePoly -> mergeL -> addEq, which performs no max at all, and addEq
            % then grouped by equal function ACROSS pieces and intersected their regions, which
            % collapses. The max is functionNDomain.maxOfList, shared with plq.maximumConjugate
            % and plq.biconjugateF so there is exactly one implementation of it.
            [bcAll, iaB] = obj.fnd.conjugateOfPiecePoly;
            groups = cell(1, numel(iaB)-1);
            for k = 1:numel(iaB)-1
                groups{k} = bcAll(iaB(k):iaB(k+1)-1);
            end
            bc = functionNDomain.maxOfList(groups);
            % An EMPTY piece list is not a function -- it evaluates to NaN everywhere, i.e. it
            % claims f* is +inf on all of R^2, which no conjugate of a proper function is. Refuse
            % it rather than hand it back.
            %
            % This became reachable when conjCPLQ's isDomBounded gate came out: biconj is
            % conj-of-conj, so a bounded triangle's f** now gets as far as this second
            % conjugation, where it used to stop at the first. The first conjugation is right --
            % for f = x*y on conv{(0,0),(1,0),(0,1)}, f* = max(0,s1,s2) and (f*)* comes out as
            % the indicator of the simplex {s>=0, s1+s2<=1}, verified exact at 9 probes. It is
            % conjugateOfPiecePoly that then returns nothing for that indicator, whose conjugate
            % is its support function max(0,x,y) -- 3 affine vertex pieces. That is a Step 2 gap,
            % unrelated to boundedness, and it is now visible instead of masked.
            if isempty(bc)
                error('QuaParCPLQ:conj:emptyResult', ...
                    ['conjugateOfPiecePoly returned no pieces, which would represent f* as ' ...
                     '+inf everywhere. Refusing: this is the Step 2 gap on an indicator-like ' ...
                     'piece, not a valid conjugate.']);
            end
            h = QuaParCPLQ(bc);
        end

        function h = scalarMul(obj, c)
            % objective: (c*f)(x) = c*f(x) -- scale every piece's value, domain unchanged.
            fnd2 = obj.fnd;
            for i = 1:numel(fnd2)
                fnd2(i).f = symbolicFunction(c * fnd2(i).f.f);
            end
            h = QuaParCPLQ(fnd2);
        end

        function h = negate(obj)
            h = scalarMul(obj, -1);
        end

        function h = addQuadratic(obj, A, b, c)
            % objective: h(x) = f(x) + (1/2 x'Ax + b'x + c) -- a full-domain quadratic is finite
            % everywhere, so it is simply added to every piece's own value symbolically.
            b = b(:);
            fnd2 = obj.fnd;
            for i = 1:numel(fnd2)
                v = fnd2(i).d.vars;
                q = 0.5*(A(1,1)*v(1)^2 + 2*A(1,2)*v(1)*v(2) + A(2,2)*v(2)^2) + b(1)*v(1) + b(2)*v(2) + c;
                fnd2(i).f = symbolicFunction(fnd2(i).f.f + q);
            end
            h = QuaParCPLQ(fnd2);
        end

        function h = addScaledEnergy(obj, alpha)
            h = addQuadratic(obj, 2*alpha*eye(2), [0;0], 0);
        end

        function h = add(obj, obj2)
            % objective: (f+g)(x) = f(x)+g(x) -- overlay every pair of pieces (same domain-
            % intersection idiom as functionNDomain.times/mtimes), summing the two values on each
            % nonempty overlap cell.
            if ~isa(obj2, 'QuaParCPLQ')
                error('QuaParCPLQ:add:unsupportedType', ...
                    'QuaParCPLQ.add: expected a QuaParCPLQ operand, got %s.', class(obj2));
            end
            f1 = obj.fnd; f2 = obj2.fnd;
            out = functionNDomain.empty(); n = 0;
            for i = 1:numel(f1)
                for j = 1:numel(f2)
                    rf = f1(i).d + f2(j).d;
                    if isempty(rf), continue; end
                    rf = rf.simplifyUnboundedRegion;
                    if isempty(rf), continue; end
                    n = n + 1;
                    out(n) = functionNDomain([symbolicFunction(f1(i).f.f + f2(j).f.f)], rf);
                end
            end
            if n == 0
                error('QuaParCPLQ:add:noOverlap', ...
                    'QuaParCPLQ.add: the two operands'' domains do not overlap anywhere.');
            end
            out = out.mergeL;
            h = QuaParCPLQ(out);
        end
    end
end
