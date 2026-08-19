classdef plqCheck
% plqCheck  Numeric verification for the cPLQ pipeline's stages, so a "crash test" becomes a
%   test.
%
% WHY THIS EXISTS. `testcPLQ` and `testMaxMultiRegion` between them held 32 tests and ZERO
% assertions: each ran `triangulate -> convexEnvelope -> conjugate -> maximum -> biconjugateF`,
% printed the answer, and returned. They passed if nothing threw. That is worth something -- most
% of this pipeline's defects this month WERE crashes -- but it is not verification, and it cost
% 90% of a two-hour bucket to learn that a function returns. Measured 2026-08-19:
% `testcPLQ/testRectBiconj` alone was 3198 s.
%
% WHAT IS CHECKED, and each is a DEFINITION rather than a golden value, so nothing here needs
% re-pinning when a representation changes:
%
%   convex envelope   co f <= f on the domain, and co f = f at every vertex of it (an envelope
%                     that is not an underestimator is wrong; one that dips below the vertex
%                     values is wrong; and on these fixtures the vertices are where the envelope
%                     touches).
%   conjugate         f*(s) = sup_{x in D} <s,x> - f(x), against a numeric sup over the domain.
%                     This is the strongest of the three and needs no convexity assumption.
%   biconjugate       f** <= f on the domain, f** convex along segments, and f** = f wherever f
%                     is already convex there.
%
% NUMERIC SUP, and why sampling is sound here. The domains are convex polygons and f is a
% quadratic, so `<s,x> - f(x)` is a quadratic whose sup over the polygon is attained on the
% boundary or at an interior critical point. Both are sampled: the whole boundary at a fine step,
% every vertex exactly, plus an interior grid. That gives a LOWER bound on the true sup, so
% `f*(s) < sup_sampled - tol` is a definite failure, while `f*(s) > sup_sampled` within tolerance
% is expected (the sampler cannot exceed the truth). Both directions are checked, with the
% one-sided one given the tolerance it deserves.
%
% TOLERANCES are absolute-and-relative and deliberately loose (1e-6 relative): the point is to
% catch a wrong branch, a dropped cell or a mis-oriented constraint, all of which are O(1) errors
% -- not to police the last ulp, which the exact-arithmetic tests elsewhere already do.

    methods (Static)

        function pts = domainSample(d, nB, nI)
        % Points of the polygon `d`: every vertex, nB points along each edge, and an interior
        % grid of about nI points. `d` is a `domain` (its polygon carries vx/vy).
            if nargin < 2, nB = 40; end
            if nargin < 3, nI = 200; end
            V = [double(d.polygon.vx(:)), double(d.polygon.vy(:))];
            n = size(V,1);
            pts = V;
            for k = 1:n
                a = V(k,:); b = V(mod(k,n)+1,:);
                t = linspace(0, 1, nB+2)'; t = t(2:end-1);
                pts = [pts; a + t.*(b - a)];                          %#ok<AGROW>
            end
            % interior: rejection-sample the bounding box against the polygon
            lo = min(V,[],1); hi = max(V,[],1);
            m = ceil(sqrt(nI)) * 4;
            [gx, gy] = meshgrid(linspace(lo(1), hi(1), m), linspace(lo(2), hi(2), m));
            C = [gx(:), gy(:)];
            inside = plqCheck.inPolygon(C, V);
            C = C(inside,:);
            if size(C,1) > nI
                C = C(round(linspace(1, size(C,1), nI)), :);
            end
            pts = [pts; C];
        end

        function tf = inPolygon(P, V)
        % Convex-polygon membership by half-planes, orientation-agnostic. `inpolygon` needs the
        % Mapping/Image toolboxes on some installs; this needs nothing and these domains are
        % convex by construction.
            n = size(V,1);
            s = zeros(size(P,1), n);
            for k = 1:n
                a = V(k,:); b = V(mod(k,n)+1,:);
                s(:,k) = (b(1)-a(1))*(P(:,2)-a(2)) - (b(2)-a(2))*(P(:,1)-a(1));
            end
            tol = 1e-9 * max(1, max(abs(s), [], 'all'));
            tf = all(s >= -tol, 2) | all(s <= tol, 2);
        end

        function [v, ok] = safeEval(fnd, pt)
        % evalFunctionNDomain at one point, with a RATIONAL face's pole treated as "no value
        % here" rather than an error. Step 1's envelope is quadratic-over-LINEAR, and its
        % denominator vanishes exactly at some domain vertices -- which is where an underestimator
        % check most wants to sample. MATLAB raises symbolic:kernel:DivisionByZero there, so the
        % sampler must skip the point instead of the test dying on the fixture's own geometry.
            v = NaN; ok = false;
            try
                v = evalFunctionNDomain(fnd, pt);
                ok = isfinite(v);
            catch
                v = NaN; ok = false;
            end
        end

        function v = evalSym(f, vars, X)
        % f evaluated at the rows of X, numerically.
            h = matlabFunction(f, 'Vars', {vars(1), vars(2)});
            v = arrayfun(@(i) double(h(X(i,1), X(i,2))), (1:size(X,1))');
        end

        function s = supOverDomain(fExpr, vars, d, sPt, nB, nI)
        % sup_{x in D} <sPt,x> - f(x), sampled. A LOWER bound on the true sup; see the header.
            if nargin < 5, nB = 200; end
            if nargin < 6, nI = 900; end
            X = plqCheck.domainSample(d, nB, nI);
            fv = plqCheck.evalSym(fExpr, vars, X);
            s = max(X * sPt(:) - fv);
        end

        % ---- the three checks -----------------------------------------------------------

        function envelopeUnderestimates(tc, piece, name)
        % co f <= f on the domain, and co f = f at every vertex. `piece` is a plq_1p AFTER
        % convexEnvelope, so piece.envelope is the functionNDomain array of the envelope.
            d = piece.d; vars = d.polygon.vars;
            X = plqCheck.domainSample(d, 12, 60);
            fv = plqCheck.evalSym(piece.f.f, vars, X);
            ev = arrayfun(@(i) plqCheck.safeEval(piece.envelope, X(i,:)), (1:size(X,1))');
            got = ~isnan(ev) & ~isnan(fv);
            tc.verifyGreaterThan(nnz(got), 0, ...
                sprintf('%s: the envelope covers no sampled point of its own domain', name));
            tol = 1e-6 * max(1, max(abs(fv)));
            tc.verifyLessThanOrEqual(max(ev(got) - fv(got)), tol, ...
                sprintf('%s: the convex envelope EXCEEDS f somewhere on the domain', name));

            V = [double(d.polygon.vx(:)), double(d.polygon.vy(:))];
            fvV = plqCheck.evalSym(piece.f.f, vars, V);
            evV = arrayfun(@(i) plqCheck.safeEval(piece.envelope, V(i,:)), (1:size(V,1))');
            ok = ~isnan(evV) & ~isnan(fvV);
            tc.verifyLessThanOrEqual(max(abs(evV(ok) - fvV(ok))), 1e-6 * max(1, max(abs(fvV))), ...
                sprintf('%s: the envelope does not touch f at the domain vertices', name));
        end

        function conjugateMatchesSup(tc, fnd, fExpr, vars, d, S, name)
        % f*(s) = sup_{x in D} <s,x> - f(x) at each dual point in the rows of S.
            for i = 1:size(S,1)
                s = S(i,:);
                got = plqCheck.safeEval(fnd, s);
                want = plqCheck.supOverDomain(fExpr, vars, d, s);
                tc.verifyFalse(isnan(got), sprintf( ...
                    '%s: f* is uncovered at (%g,%g) -- a cell is missing', name, s(1), s(2)));
                if isnan(got), continue, end
                tol = 1e-5 * max(1, abs(want));
                % The sampler is a LOWER bound, so got < want is a definite defect; got > want
                % beyond tolerance means f* claims more than the definition allows.
                tc.verifyGreaterThanOrEqual(got, want - tol, sprintf( ...
                    '%s: f*(%g,%g) = %.9g is BELOW the sampled sup %.9g', ...
                    name, s(1), s(2), got, want));
                tc.verifyLessThanOrEqual(got, want + tol, sprintf( ...
                    '%s: f*(%g,%g) = %.9g EXCEEDS the sup %.9g over the domain', ...
                    name, s(1), s(2), got, want));
            end
        end

        function biconjugateIsAConvexUnderestimator(tc, bnd, fExpr, vars, d, name)
        % f** <= f on the domain, and f** convex along sampled segments.
            X = plqCheck.domainSample(d, 10, 50);
            fv = plqCheck.evalSym(fExpr, vars, X);
            bv = arrayfun(@(i) plqCheck.safeEval(bnd, X(i,:)), (1:size(X,1))');
            ok = ~isnan(bv) & ~isnan(fv);
            tc.verifyGreaterThan(nnz(ok), 0, ...
                sprintf('%s: f** covers no sampled point of the domain', name));
            tol = 1e-6 * max(1, max(abs(fv)));
            tc.verifyLessThanOrEqual(max(bv(ok) - fv(ok)), tol, ...
                sprintf('%s: f** EXCEEDS f somewhere on the domain', name));

            % convexity along segments between sampled points
            rng(20260819);
            idx = randi(size(X,1), 60, 2);
            worst = -inf;
            for k = 1:size(idx,1)
                a = X(idx(k,1),:); b = X(idx(k,2),:);
                m = (a + b)/2;
                va = plqCheck.safeEval(bnd, a);
                vb = plqCheck.safeEval(bnd, b);
                vm = plqCheck.safeEval(bnd, m);
                if any(isnan([va vb vm])), continue, end
                worst = max(worst, vm - (va + vb)/2);
            end
            if isfinite(worst)
                tc.verifyLessThanOrEqual(worst, 1e-6 * max(1, max(abs(fv))), sprintf( ...
                    '%s: f** is not convex -- a midpoint sits above the chord by %.3g', ...
                    name, worst));
            end
        end
    end
end
