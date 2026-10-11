classdef cornerFaceAlgTest < matlab.unittest.TestCase
% [fast] cornerFaceAlg, checked against the DEFINITION of a corner face rather than stored output:
% at the corner s -- a point of the tie on one of the cell's facets -- the face's value is the
% AFFINE function <s,x> - g(s), so the implicit polynomial must vanish on it and the chosen branch
% must equal it.
%
% The inputs are the corner calls biconjQ makes on the poly2 kink at c = -2 (captured 2026-09-25).
% Five of them used to raise ratQ:overflow: the square in the implicit polynomial passed 2^53 while
% the polynomial, content removed, does not.

    methods (Test)
        function aFACETParallelToTheParabolasAXISHasOneRATIONALCorner(t)
        % Found 2026-10-10 (biconjAtoms on an indefinite triangle): a vertex's affine conjugate
        % against an edge's rank-1 one ties on a PARABOLA whose axis is the null direction of the
        % rank-1 Hessian -- PARALLEL to the strip's own facets. Restricted to such a facet the tie
        % is LINEAR in the parameter (al = 0), so the facet meets the parabola exactly ONCE, at a
        % rational point. The routine skipped every such facet as "along the asymptote", the pair
        % lost its corners, and envelopeEval answered +inf on the whole envelope.
        %
        % g(s) = s1^2/2 on the strip -1 <= s1 <= 1, tie g - s2 = 0, i.e. 2 s2 = s1^2: the corners
        % are (+-1, 1/2), and each face's value is <s,x> - g(s), affine.
            n = [0 0 0 0 1 0 0 0 0 0];  P = [1 0 1; -1 0 1];  tie = [1 0 0 0 -2 0];
            faces = cornerFaceAlg(n, 1, P, tie, 1);
            t.verifyEqual(numel(faces), 2, 'one corner on each facet');
            S = sortrows([faces.s].');
            t.verifyEqual(S, [-1 0.5; 1 0.5], 'AbsTol', 1e-14);
            rng(11);  X = randn(20,2)*3;
            for k = 1:numel(faces)
                s = faces(k).s;
                want = X*s - s(1)^2/2;
                got = arrayfun(@(i) algBranch(faces(k).FP, X(i,:), faces(k).FB), (1:20).');
                t.verifyEqual(got, want, 'AbsTol', 1e-10, 'the value is the affine <s,x> - g(s)');
                t.verifyEqual(size(faces(k).con, 1), 2, 'and the wedge is spelled as two rows');
                % the rows must say what the KKT predicate says: M \ (p - x) >= 0
                Y = randn(200,2)*3;
                byRows = true(200,1);
                for r = 1:2
                    byRows = byRows & faces(k).con(r,7) * QuaPar.evalConic(faces(k).con(r,1:6), Y) >= 0;
                end
                byKkt = arrayfun(@(i) all(faces(k).M \ (faces(k).p - Y(i,:).') >= 0), (1:200).');
                t.verifyEqual(byRows, byKkt, 'the wedge rows and the KKT test agree');
            end
        end

        function aWEDGEEdgeThatIsARationalLineTimesAnIrrationalScalarGetsItsRow(t)
        % Found 2026-10-11, the poly2 kink at c = -2 (pair 15: an affine face against a rank-1 one,
        % tied on irrational parallel lines). A wedge edge A + sigma*B with A PARALLEL to B is the
        % RATIONAL line B scaled by (lambda + sigma); its conjugate product is a DOUBLE line, which
        % pairLineData rejects (ratQ:notALinePair), so the corner got no rows and the mesh could
        % not place it. Every corner of the pair must come with two rows that agree with its KKT
        % predicate.
            P = [-1 2 4; 1 -2 6; 2 1 -8; 3 2 -10; 1 0 -2];
            q = [1 -4 4 12 -24 -44];
            calls = {{[0 0 0 0 1 -2 4 16 8 -32], 10, P, q,  1}, ...
                     {[0 0 0 0 0 0 0 1 2 -1],      1, P, q, -1}};
            rng(13);  Y = randn(300,2)*4;
            n = 0;
            for c = 1:numel(calls)
                faces = cornerFaceAlg(calls{c}{:});
                for k = 1:numel(faces)
                    n = n + 1;
                    t.verifyEqual(size(faces(k).con, 1), 2, sprintf('call %d corner %d has its two rows', c, k));
                    if size(faces(k).con, 1) ~= 2, continue, end
                    byRows = true(size(Y,1),1);
                    for r = 1:2
                        cr = faces(k).con(r,:);
                        if cr(8) ~= 0, qv = ratQ.pairLineEval(cr(1:6), cr(8), Y);
                        else, qv = QuaPar.evalConic(cr(1:6), Y); end
                        byRows = byRows & cr(7) * qv(:) >= 0;
                    end
                    byKkt = arrayfun(@(i) all(faces(k).M \ (faces(k).p - Y(i,:).') >= 0), (1:size(Y,1)).');
                    t.verifyEqual(byRows, byKkt, 'the wedge rows and the KKT test agree');
                end
            end
            t.verifyGreaterThan(n, 0, 'the fixture has corners');
        end

        function theCORNERSOfThePOLY2KinkBuildAndCarryTheirValue(t)
            calls = {{[0 0 0 0 1 0 0 0 0 0],     2, [1 0 4;-1 0 4;0 -1 0;3 2 16;-3 -2 10;2 -3 2], [8 -24 -8 16 -128 -96],  1}, ...
                     {[0 0 0 0 9 6 4 -4 32 24], 26, [1 0 4;-1 0 4;0 -1 0;3 2 16;-3 -2 10;2 -3 2], [8 -24 -8 16 -128 -96], -1}, ...
                     {[0 0 0 0 1 0 1 0 0 0],     2, [0 1 0;-2 -1 8;2 -3 8;1 0 -2;3 2 -10],         [1 0 1 -4 -8 4],          1}, ...
                     {[0 0 0 0 0 0 0 1 2 -1],    1, [0 1 0;-2 -1 8;2 -3 8;1 0 -2;3 2 -10],         [1 0 1 -4 -8 4],         -1}, ...
                     {[0 0 0 0 1 0 1 0 0 0],     2, [0 1 0;-2 -1 8;2 -3 8;3 2 16;-3 -2 10;2 -3 2], [8 -24 18 16 -128 -96],  1}, ...
                     {[0 0 0 0 9 6 4 -4 32 24], 26, [0 1 0;-2 -1 8;2 -3 8;3 2 16;-3 -2 10;2 -3 2], [8 -24 18 16 -128 -96], -1}};
            nFaces = 0;
            for k = 1:numel(calls)
                a = calls{k};
                faces = cornerFaceAlg(a{:});
                n = a{1};  d = a{2};
                N = [n(5) n(6); n(6) n(7)];  m = [n(8); n(9)];  q = n(10);
                g = @(s) ((s.'*N*s)/2 + m.'*s + q) / d;
                for f = faces
                    nFaces = nFaces + 1;
                    s = f.s;
                    for off = [0 0; 0.1 0; 0 0.1; -0.07 0.05].'
                        x = f.p + off;                               % around grad g at the corner
                        z = s.'*x - g(s);
                        t.verifyLessThan(relResidual(f.FP, x, z), 1e-10, ...
                            sprintf('call %d: the polynomial vanishes on the corner''s value', k));
                    end
                    % THE BRANCH only at the corner's own point: the index names this corner inside
                    % its region, and an offset can leave it -- for an affine g the region at p is a
                    % cone, and beside it the top root is the CONJUGATE corner's value. The
                    % polynomial vanishes on both, which is why the residual above is asked everywhere.
                    t.verifyEqual(algBranch(f.FP, f.p.', f.FB), s.'*f.p - g(s), 'AbsTol', 1e-7, ...
                        sprintf('call %d: the chosen branch IS the corner''s value at grad g(s)', k));
                end
            end
            t.verifyGreaterThan(nFaces, 0, 'the fixture must produce corners, or this proves nothing');
        end
    end
end

function r = relResidual(P, x, z)
% |P(x,y,z)| relative to the sum of its terms' magnitudes.
    v = 0;  a = 0;
    for i = 1:size(P,1)
        for j = 1:size(P,2)
            for k = 1:size(P,3)
                term = P(i,j,k) * x(1)^(i-1) * x(2)^(j-1) * z^(k-1);
                v = v + term;  a = a + abs(term);
            end
        end
    end
    r = abs(v) / max(a, realmin);
end
