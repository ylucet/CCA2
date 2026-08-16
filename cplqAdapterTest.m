classdef cplqAdapterTest < matlab.unittest.TestCase
% cplqAdapterTest  End-to-end tests for the Phase 1 cPLQ integration adapter (quaPolToPlq.m /
%   evalFunctionNDomain.m): CCA2 QuaPol -> cPLQ plq -> triangulate/conjugate/maximum ->
%   evalFunctionNDomain, validated against numeric sup-sampling ground truth (this codebase's
%   standard convention, e.g. conjPieceCPLQTest.m). See DESIGN.md II.5.1 and
%   .claude/SESSION_HANDOFF.md.

    methods (Test)

        function singleTriangleMatchesConjPieceCPLQ(testCase)
            % Sanity check on the INPUT conversion only: f=xy over a single one-convex-edge
            % triangle, run through the cPLQ pipeline (.conjugate, no .maximum needed for one
            % piece), must agree with CCA2's own existing numeric conjPieceCPLQ implementation.
            V = [0 0; 2 0; 1 1]; E = [1 2 1; 2 3 1; 3 1 1]; F = [1 0; 1 0; 1 0];
            q = QuaPol(V, E, [0 1 0 0 0 0], F);
            g = conjPieceCPLQ(q);   % existing numeric implementation (ground truth-checked elsewhere)

            p = quaPolToPlq(q);
            testCase.verifyEqual(p.nPieces, 1);
            p = p.triangulate;
            p = p.conjugate;        % Step 1 (envelope) + Step 2 (conjugate) via cPLQ, symbolic

            S = [0.5 0.5; 3 -1; -2 3; 1 1; 0 -3; 4 4; -3 -3; 0.5 1.5];
            for i = 1:size(S,1)
                gv = evalFunctionNDomain(p.pieces(1).conjugates, S(i,:));
                testCase.verifyEqual(gv, g.eval(S(i,:)), 'AbsTol', 1e-6, sprintf('s=%d', i));
            end
        end

        function twoTriangleSquareMaxMatchesNumericSup(testCase)
            % The case CCA2's own conjCPLQ.m cannot handle yet (nf>1, general Step 3): f=xy over
            % the square [0,1]^2 split into two triangles along the diagonal. Each triangle's own
            % conjugate is individually correct (already exercised by conjPieceCPLQ), but
            % maxQuaPar refuses to combine them (curved edges) -- this is exactly the gap the cPLQ
            % integration closes. Validate against the true numeric sup over the WHOLE square.
            V = [0 0; 1 0; 1 1; 0 1];
            E = [1 2 1; 2 3 1; 3 1 1; 3 4 1; 4 1 1];
            F = [1 0; 1 0; 1 2; 2 0; 2 0];
            f = [0 1 0 0 0 0; 0 1 0 0 0 0];   % xy on both faces
            q = QuaPol(V, E, f, F);

            p = quaPolToPlq(q);
            testCase.verifyEqual(p.nPieces, 2);
            p = p.triangulate;
            p = p.maximum;           % Steps 1+2+3 via cPLQ: envelope, conjugate, max

            nt = 220; [uu,vv] = meshgrid(linspace(0,1,nt));
            Xg = uu(:); Yg = vv(:); xyg = Xg.*Yg;
            % s=(0.5,0.5) is an exact symmetric tie point where both triangles' own
            % off-diagonal-vertex cones meet (true sup 0.5, attained at BOTH (1,0) [triangle 1]
            % and (0,1) [triangle 2] simultaneously). Previously excluded here: region.maxArray's
            % overlap-domain vertices are themselves exactly on the f1==f2 tie line at this
            % junction, so its interior-probe heuristic gave up (lSing) and functionNDomain.
            % maximumP discarded the whole overlap region instead of falling back to splitmax3 --
            % fixed in maximumP (see its own HISTORY comment / .claude/SESSION_HANDOFF.md).
            S = [3 -1; -2 3; 1 1; 0 -3; 4 4; -3 -3; 6 2; -1 6; 2 2; 0.5 0.5];
            for i = 1:size(S,1)
                sup = max(S(i,1)*Xg + S(i,2)*Yg - xyg);
                gv = evalFunctionNDomain(p.maxConjugate, S(i,:));
                testCase.verifyEqual(gv, sup, 'AbsTol', 2e-3, sprintf('s=%d', i));
            end
        end

        function generalQuadrilateralStep1IsTheEnvelopeNotAMinorant(testCase)
        % f = x*y over conv{(0,0),(2,0),(2.5,1.5),(0.5,1)}. Step 1 must produce, for EVERY piece,
        % the CONVEX ENVELOPE of x*y on that piece -- not merely a convex minorant, and not
        % nothing at all. Today it does neither, in two different ways, and that is the whole of
        % the general-quadrilateral failure (SUPPORT_MATRIX.md 7.1).
        %
        % triangulate splits this quadrilateral into conv{(2.5,1.5),(2,0),(0,0)}, which has two
        % convex edges, and conv{(2.5,1.5),(0,0),(0.5,1)}, which has three.
        %
        %   * THREE convex edges: [COAP] has no single-quadratic form for that case at all, and
        %     plq_1p.convexEnvelope1 branches on nCE == 0, 1, 2 and then falls off the end -- so
        %     the piece keeps an EMPTY envelope, and plq_1p.conjugate's
        %     `for i = 1:max(1,size(envelope,2))` indexes envelope(1) and raises
        %     MATLAB:badsubscript.
        %   * TWO convex edges: the single-quadratic form touches x*y along both convex edges and
        %     is a valid convex MINORANT, but Appendix A.4 shows it is tight only over a
        %     sub-region -- and this branch never tests. Both coordinates are >= 0 on that
        %     triangle, so x*y >= 0 there, the affine minorant 0 is admissible, and the true
        %     envelope is >= 0 EVERYWHERE on it. The form returned reaches -0.2835 at (1,0).
        %
        % Both checks below are exact statements about the envelope, not comparisons against a
        % sampled reference: an envelope must exist, must be <= x*y (a minorant), and on a
        % triangle where x*y >= 0 must be >= 0.
            V = [0 0; 2 0; 2.5 1.5; 0.5 1];
            E = [1 2 1; 2 3 1; 3 4 1; 4 1 1]; F = [1 0; 1 0; 1 0; 1 0];
            q = QuaPol(V, E, [0 1 0 0 0 0], F);

            p = quaPolToPlq(q);
            p = p.triangulate;
            p = p.convexEnvelope;

            x = sym('x'); y = sym('y');
            for i = 1:p.nPieces
                testCase.verifyGreaterThanOrEqual(size(p.pieces(i).envelope,2), 1, sprintf( ...
                    'piece %d has no convex envelope at all', i));
                for k = 1:size(p.pieces(i).envelope,2)
                    W = double([p.pieces(i).envelope(k).d.vx(:), ...
                                p.pieces(i).envelope(k).d.vy(:)]);
                    e = matlabFunction(p.pieces(i).envelope(k).f.f, 'Vars', [x y]);
                    G = cplqAdapterTest.samplePolygon(W, 60);
                    v = e(G(:,1), G(:,2));
                    % A one-convex-edge envelope is RATIONAL ([COAP] A.3 eq.16) and its
                    % denominator vanishes at one vertex of the face, where the numerator vanishes
                    % too -- a REMOVABLE singularity, and the value there is the limit. Sampling
                    % lands on it exactly, so drop the non-finite samples rather than pretend the
                    % envelope is infinite there. Everything else is checked.
                    fin = isfinite(v);
                    testCase.verifyGreaterThan(nnz(fin), 0.9*numel(v), sprintf( ...
                        'piece %d face %d: the envelope is non-finite on more than a curve', i, k));
                    G = G(fin,:); v = v(fin);
                    testCase.verifyLessThanOrEqual(max(v - G(:,1).*G(:,2)), 1e-7, sprintf( ...
                        'piece %d face %d: the envelope must be a MINORANT of x*y', i, k));
                    if all(W(:,1) >= -1e-12) && all(W(:,2) >= -1e-12)
                        % x*y >= 0 on this face, so the affine minorant 0 is admissible and the
                        % true convex envelope is >= 0 everywhere on it.
                        testCase.verifyGreaterThanOrEqual(min(v), -1e-7, sprintf( ...
                            ['piece %d face %d: x*y >= 0 here, so the ENVELOPE is >= 0 -- ' ...
                             'anything below it is a minorant that is not tight'], i, k));
                    end
                end
            end
        end

        function generalQuadrilateralConjugateMatchesTheSup(testCase)
        % The end of the same story: with Step 1 right, f* of x*y over that quadrilateral must
        % match the sup over the domain. Each value below is attained at a VERTEX -- x*y is
        % bilinear, so its sup against a linear form over a polygon is -- which makes the
        % reference exact rather than sampled: the max over the four vertices of <s,v> − v1·v2.
        %
        % WHY THE PER-PIECE MAX AND NOT `q.conj('cplq')`. They are the same function:
        % f* = (conv f)* = max_k (env_k + I_{face k})*, which is what Step 3 assembles. But
        % assembling it is what costs, and the cost is Step 3's, not this fix's. MEASURED on this
        % input, with the split in: Step 1 and Step 2 take about 25 s for all six pieces, and the
        % cross-piece fold then takes **73 minutes** -- 93 s, 294 s, 647 s, 1273 s, 2087 s per
        % fold, with the cell count running 5, 14, 29, 45, 70, 86. The assembled answer is exact
        % (8 of 8 probe points), so this test loses no coverage of the fix by taking the max
        % itself; what it avoids is paying for the known Step 3 blow-up, which is the same one
        % SUPPORT_MATRIX.md records for a pentagon (885 s, 41 regions) and is tracked separately.
        %
        % The uncovered count matters as much as the values: f* of a BOUNDED domain is finite
        % everywhere, so every piece must answer at every point, and a hole would mean a piece's
        % conjugate domain is too small.
            V = [0 0; 2 0; 2.5 1.5; 0.5 1];
            E = [1 2 1; 2 3 1; 3 4 1; 4 1 1]; F = [1 0; 1 0; 1 0; 1 0];
            q = QuaPol(V, E, [0 1 0 0 0 0], F);

            p = quaPolToPlq(q);
            p = p.triangulate;
            for i = 1:p.nPieces
                p.pieces(i) = p.pieces(i).convexEnvelope;
                p.pieces(i) = p.pieces(i).conjugate;
                p.pieces(i) = p.pieces(i).maximumConjugate;
            end

            S = [0 0; 1 0.5; 0.5 1; 2 1; -1 0.5; 1.5 -0.5; 3 2; -2 -2; 0.25 0.25; 4 -1];
            fv = V(:,1).*V(:,2);
            for i = 1:size(S,1)
                want = max(S(i,1)*V(:,1) + S(i,2)*V(:,2) - fv);
                got = -inf;
                for k = 1:p.nPieces
                    v = evalFunctionNDomain(p.pieces(k).maxConjugate, S(i,:));
                    testCase.verifyFalse(isnan(v), sprintf( ...
                        'piece %d does not cover (%g,%g), but a bounded piece''s conjugate is finite everywhere', ...
                        k, S(i,1), S(i,2)));
                    got = max(got, v);
                end
                testCase.verifyEqual(got, want, 'AbsTol', 2e-3, sprintf( ...
                    'f* at (%g,%g)', S(i,1), S(i,2)));
            end
        end

    end

    methods (Static)
        function G = samplePolygon(W, n)
        % Points of the convex polygon W (rows = vertices), as barycentric combinations of its
        % vertices -- which needs no inpolygon and no bounding box, and lands exactly on the
        % vertices and edges as well as the interior.
            G = zeros(0,2);
            m = size(W,1);
            if m == 3
                for a = 0:n
                    for b = 0:(n-a)
                        w = [a, b, n-a-b]/n;
                        G(end+1,:) = w*W;              %#ok<AGROW>
                    end
                end
            else
                for a = 0:n
                    for b = 0:n
                        w = [a*b, a*(n-b), (n-a)*(n-b), (n-a)*b]/n^2;
                        G(end+1,:) = (w/sum(w))*W;     %#ok<AGROW>
                    end
                end
            end
        end
    end
end
