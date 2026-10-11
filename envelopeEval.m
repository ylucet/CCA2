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
    % EVERY PAIR MAY DECLARE DIVERGENCE, because every pair's candidate set is now complete: the
    % straight candidates from `envelopeCellRegions`, the curved facet's interior from
    % `bridgeFaceAlg`, and the CORNERS where the two meet from `cornerFaceAlg`. Their union is the
    % set where that pair's sup is finite, so "nothing covers x" means "+infinity" -- which is how
    % dom f** is expressed.
    %
    % It was not always so, and the history is the reason the check exists: before the corners were
    % built, the two curved pairs of the c = -2 fixture covered 320 and 280 of 600 interior points,
    % and letting them vote made the envelope +infinity everywhere. Measured again after: 600 of
    % 600, for every pair.
    % A face may belong to SEVERAL pairs once the prune has run: a dominating face inherits the
    % pairs of everything it replaced, so that the coverage those faces provided is not lost.
    allPairs = [];
    for k = 1:numel(env.quad), allPairs = [allPairs, env.quad(k).pair(:).']; end %#ok<AGROW>
    pairs = unique([allPairs, cellfun(@(f) f.pairIndex, env.bridge)]);
    best = -inf(n,1);
    bestKind = repmat({'none'}, n, 1);
    bestIdx = zeros(n,1);
    diverges = false(n,1);

    for p = pairs
        cover = false(n,1);
        pv = -inf(n,1);  pk = repmat({'none'}, n, 1);  pi = zeros(n,1);
        for k = find(arrayfun(@(f) any(f.pair == p), env.quad))
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
            in = faceOwns(face, X);
            idx = find(in);
            for j = idx(:).'
                v = bridgeValue(face, X(j,:));
                if ~isfinite(v), continue, end
                cover(j) = true;
                if v > pv(j), pv(j) = v;  pk{j} = 'bridge';  pi(j) = k; end
            end
        end
        diverges = diverges | ~cover;
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
function in = faceOwns(face, X)
% Which of the curved-tie faces owns each point. A BRIDGE face is swept, so its membership is the
% existence of a sweeping segment (`bridgeMember`); a CORNER face is a WEDGE whose apex is a fixed
% irrational point, so its membership is the KKT multipliers of the two constraints active there
% being nonnegative -- two linear tests against a fixed matrix.
    switch face.kind
        case 'corner'
            in = false(size(X,1),1);
            for i = 1:size(X,1)
                mu = face.M \ (face.p - X(i,:).');
                in(i) = all(mu >= -1e-9 * max(1, max(abs(mu))));
            end
        otherwise
            in = bridgeMember(face.pair, X);
    end
end

% --------------------------------------------------------------------------------------------------
function v = bridgeValue(face, x)
% A bridge face carries its value either implicitly (a polynomial and a branch) or as a rational
% function -- `face.kind` says which, and both are read here so a caller never has to ask.
    switch face.kind
        case 'corner'
            % A corner face is AFFINE, <s,x> - g(s), and is read as such. NOT through FP/FB: FP's two
            % roots are the two conjugate corners' PLANES, which CROSS, so "the k-th root from the
            % top" names this corner on one side of their crossing line and the other corner on the
            % other side -- measured 2026-10-11 on the poly2 kink, -0.382 returned for -2.642,
            % above f** (envelopeAssembleTest/anIRRATIONALStraightTieGetsItsCorners).
            n = face.src.n;  s = face.s(:);
            g = ((s.'*[n(5) n(6); n(6) n(7)]*s)/2 + [n(8) n(9)]*s + n(10)) / face.src.d;
            v = x(:).'*s - g;
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
                ['a curved-tie face must be ''algebraic'', ''rational'' or ''corner''; got ' ...
                 '''%s''.'], face.kind);
    end
end

% --------------------------------------------------------------------------------------------------
function in = inRegion(con, X)
% Does each point satisfy every row of this face's region? The tolerance is the one QuaCon.eval
% already documents: a point within it of a boundary is accepted, and on the boundary two faces
% agree anyway, because the envelope is continuous there.
%
% The scale is PER POINT (2026-10-10, envelopeEvalTest): scaled by the whole call's largest |x|, a
% far point in the same call widened every other point's band, and a point outside a face's
% region took that face's value -- QuaCon.eval's defect, fixed there the same day.
    n = size(X,1);
    in = true(n,1);
    xmax = max(1, max(abs(X), [], 2));
    for r = 1:size(con,1)
        c = con(r,1:6);
        if size(con,2) >= 8 && con(r,8) ~= 0
            q = ratQ.pairLineEval(c, con(r,8), X);
            scale = xmax;
        else
            q = QuaPar.evalConic(c, X);
            scale = max(1, max(abs(c)) * xmax.^2);
        end
        in = in & (con(r,7) * q(:) >= -sqrt(eps) * scale);
        if ~any(in), return, end
    end
end
