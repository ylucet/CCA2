function [val, who] = envelopeEval(env, X)
% envelopeEval  The value of an assembled convex envelope: the MAX over the faces whose region
%   contains the point.
%
% objective: `envelopeAssemble` returns the envelope as its FACES rather than as one mesh, because
%   folding them pairwise explodes (its header has the measurement). This is how such an envelope is
%   evaluated, and the rule is the one the assembly rests on:
%
%       f**(x) = MAX over the faces that own x of that face's value.
%
%   No gluing argument is needed. Every face's value is `<s,x> - G(s)` at a point s that is FEASIBLE
%   for its own cell pair, so no face can exceed f**; and the pairs cover the dual plane, so some
%   face attains it. `.claude/assemblyProbe.m` measured both halves directly.
%
% [input]  env : as `envelopeAssemble` returns it -- `.quad` (candidates with conic regions) and
%                `.bridge` (curved-tie faces, as `bridgeFaceAlg` returns them)
%          X   : n x 2 double
% [output] val : n x 1, +Inf where no face owns the point -- which is the answer, not a gap: the
%                sup diverges there, so the envelope IS +infinity
%          who : n x 1 struct with `kind` ('quad' | 'bridge' | 'none') and `idx`, for a caller that
%                needs to know WHICH face answered -- the paper's classification table wants it.
%
% EVALUATION IS IN DOUBLES, and deliberately: CONJ_FIELD_PROOF.md 8.0 lists the operations that must
% leave Q -- building a vertex, the sign of a rational polynomial at one, comparing two of them --
% and evaluating a finished function at a point is not among them. Every REGION here was decided
% exactly; a bridge face's membership is decided by `bridgeMember`, which is the double-precision
% counterpart of an exact Sturm test.

    if nargin ~= 2
        error('envelopeEval:nargin', 'expects (env, X).');
    end
    if ~isstruct(env) || ~all(isfield(env, {'quad', 'bridge'}))
        error('envelopeEval:shape', ...
            'env must be what envelopeAssemble returns, with .quad and .bridge.');
    end
    if ~isnumeric(X) || size(X,2) ~= 2
        error('envelopeEval:points', 'X must be n x 2; got %s.', mat2str(size(X)));
    end

    n = size(X,1);
    val = inf(n,1);
    who = repmat(struct('kind', 'none', 'idx', 0), n, 1);
    if isempty(env.quad) && isempty(env.bridge), return, end

    % ---- THE MAX IS OVER PAIRS, AND A PAIR THAT DIVERGES CONTRIBUTES +INFINITY -----------------
    % f**(x) = max over the cell pairs of (sup over that pair), and a pair whose own sup DIVERGES at
    % x makes the whole max +infinity. That is how the envelope's DOMAIN is expressed: outside
    % dom f** some pair always diverges. Taking the max of the finite faces alone would answer a
    % finite value everywhere any single pair happened to be bounded -- measured, that reported
    % finite values outside the unit square for a fixture whose envelope is +infinity there.
    %
    % A pair diverges at x exactly where NONE of its candidates owns x: a candidate exists precisely
    % where the sup is attained (envelopeCellRegions' own header states that, and its test asserts
    % the complement).
    % ONLY A PAIR WHOSE CANDIDATE SET IS COMPLETE MAY DECLARE DIVERGENCE. For a poly, tangential or
    % straight-tie pair the candidates are exactly `envelopeCellRegions`', whose union IS the set
    % where that pair's sup is finite -- so "nothing covers x" means "+infinity" there. A CURVED-tie
    % pair is not yet complete: the restricted problem can also be maximised at a CORNER, where the
    % tie meets a facet of the pair, and that candidate is an affine face with quadratic-irrational
    % coefficients which nothing here builds yet. Measured: the two curved pairs of the c = -2
    % fixture cover 320 and 280 of 600 interior points, so letting them vote made the envelope
    % +infinity everywhere. They still contribute VALUES -- every one of which is a lower bound --
    % they simply do not get to say where the sup diverges.
    pairs = unique([[env.quad.pair], cellfun(@(f) f.pairIndex, env.bridge)]);
    curved = unique(cellfun(@(f) f.pairIndex, env.bridge));
    best = -inf(n,1);
    bestKind = repmat({'none'}, n, 1);
    bestIdx = zeros(n,1);
    diverges = false(n,1);

    for p = pairs
        cover = false(n,1);
        pv = -inf(n,1);  pk = repmat({'none'}, n, 1);  pi = zeros(n,1);
        for k = find([env.quad.pair] == p)
            f = env.quad(k);
            in = inRegion(f.con, X);
            if ~any(in), continue, end
            v = ratQ.evalFace(f.num, f.den, X(in,:));
            idx = find(in);
            upd = v > pv(idx);
            pv(idx(upd)) = v(upd);
            pk(idx(upd)) = {'quad'};
            pi(idx(upd)) = k;
            cover(idx) = true;
        end
        for k = 1:numel(env.bridge)
            face = env.bridge{k};
            if face.pairIndex ~= p, continue, end
            in = bridgeMember(face.pair, X);
            idx = find(in);
            for j = idx(:).'
                v = bridgeValue(face, X(j,:));
                if ~isfinite(v), continue, end
                cover(j) = true;
                if v > pv(j), pv(j) = v;  pk{j} = 'bridge';  pi(j) = k; end
            end
        end
        if ~ismember(p, curved)
            diverges = diverges | ~cover;
        end
        upd = cover & (pv > best);
        best(upd) = pv(upd);
        bestKind(upd) = pk(upd);
        bestIdx(upd) = pi(upd);
    end

    fin = ~diverges & isfinite(best);
    val(fin) = best(fin);
    for i = find(fin).'
        who(i).kind = bestKind{i};
        who(i).idx  = bestIdx(i);
    end
end

% ==================================================================================================
function v = bridgeValue(face, x)
% A bridge face carries its value either implicitly (a polynomial and a branch) or as a rational
% function -- `face.kind` says which, and both are read here so a caller never has to ask.
    switch face.kind
        case 'algebraic'
            v = algBranch(face.FP, x, face.FB);
        case 'rational'
            q = face.fN;
            num = q(5)/2*x(1)^2 + q(6)*x(1)*x(2) + q(7)/2*x(2)^2 + q(8)*x(1) + q(9)*x(2) + q(10);
            den = face.den(1)*x(1) + face.den(2)*x(2) + face.den(3);
            if den == 0, v = NaN; return, end
            v = num/den / face.fD;
        otherwise
            error('envelopeEval:faceKind', ...
                'a bridge face must be ''algebraic'' or ''rational''; got ''%s''.', face.kind);
    end
end

% --------------------------------------------------------------------------------------------------
function in = inRegion(con, X)
% Does each point satisfy every row of this face's region? The tolerance is the one QuaCon.eval
% already documents: a point within it of a boundary is accepted, and on the boundary two faces
% agree anyway, because the envelope is continuous there.
    n = size(X,1);
    in = true(n,1);
    for r = 1:size(con,1)
        c = con(r,1:6);
        if size(con,2) >= 8 && con(r,8) ~= 0
            q = ratQ.pairLineEval(c, con(r,8), X);
            scale = max(1, max(abs(X(:))));
        else
            q = QuaPar.evalConic(c, X);
            scale = max(1, max(abs(c)) * max(1, max(abs(X(:)))^2));
        end
        in = in & (con(r,7) * q >= -sqrt(eps) * scale);
        if ~any(in), return, end
    end
end
