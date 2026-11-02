classdef triangleFaceTest < matlab.unittest.TestCase
% triangleFace / triangleFaceSide / triangleFaceSign: Example 2's triangle faces held exactly by the
% name of their dual point -- step (b) of TODO item 1 (DECISIONS 2026-10-30 (fifth)).
%
% BUCKET: fast. biconjAtoms(f, 'triangles') skips the pair assembly (~1 s); the rest is integer
% arithmetic and Sturm counts.

    properties
        f
        tris
        faces
    end

    methods (TestClassSetup)
        function build(t)
            t.f = erratumExamplesTest.example2();
            t.tris = biconjAtoms(t.f, 'triangles');
            t.faces = arrayfun(@triangleFace, t.tris);
        end
    end

    methods (Test)

        function theEXACTFormsReproduceTheDoublePipeline(t)
        % formQ's conics against the forms the double pipeline tied at s0: each g_i(s0) is the
        % triangle's c, and each gradient is its contact corner.
            for k = 1:numel(t.tris)
                tr = t.tris(k);  s = tr.s(:);
                for i = 1:3
                    G = tr.Gn(i,:) / tr.Gd(i);
                    g = G(1)*s(1)^2 + G(2)*s(1)*s(2) + G(3)*s(2)^2 + G(4)*s(1) + G(5)*s(2) + G(6);
                    t.verifyEqual(g, tr.c, 'RelTol', 1e-12, 'AbsTol', 1e-9);
                    t.verifyEqual([2*G(1) G(2); G(2) 2*G(3)] * s + G(4:5).', tr.T(i,:).', ...
                                  'AbsTol', 1e-9);
                end
            end
        end

        function exampleTwoHasOneDEGREEFourAndOneRATIONALTriangle(t)
        % conv{G1,G2,G3} ties three EDGE forms (two conics, eliminant of degree 4); the other ties
        % three VERTEX forms, which are affine, so its s0 is rational.
            t.verifyEqual(sort({t.faces.kind}), {'algebraic', 'rational'});
            a = t.faces(strcmp({t.faces.kind}, 'algebraic'));
            t.verifyEqual(numel(a.dp.R), 5, 'the eliminant is a quartic');
            for k = 1:numel(t.faces)
                t.verifyEqual(t.faces(k).s, t.tris(k).s(:), 'AbsTol', 1e-9);
            end
        end

        function theTHIRDTieHoldsExactlyAtTheNamedPoint(t)
        % s0 is named by g1 = g2 and g1 = g3; g2 = g3 is not in the name, so its exact vanishing
        % says the name is of THIS point. One unit off, it must not vanish.
            for k = 1:numel(t.faces)
                tf = t.faces(k);
                c23 = tf.Gn(2,:) * tf.Gd(3) - tf.Gn(3,:) * tf.Gd(2);
                t.verifyEqual(triangleFaceSign(tf, conicArr(c23)), 0, tf.kind);
                t.verifyNotEqual(triangleFaceSign(tf, conicArr(c23 + [0 0 0 0 0 1])), 0, tf.kind);
            end
        end

        function membershipAgreesWithBarycentricOffTheEdges(t)
        % Property: at dyadic points clear of the edges by a margin, the exact sides say what the
        % double barycentric coordinates say; both sides of every edge are reached.
            rng(51);
            for k = 1:numel(t.faces)
                tf = t.faces(k);  T = tf.T;
                lo = min(T) - 0.3 * (max(T) - min(T));  hi = max(T) + 0.3 * (max(T) - min(T));
                nIn = 0;  nOut = 0;
                for r = 1:150
                    x = round((lo + rand(1,2) .* (hi - lo)) * 64) / 64;
                    lam = [T.' ; 1 1 1] \ [x(:); 1];
                    if min(abs(lam)) < 1e-6, continue, end
                    sd = triangleFaceSide(tf, x * 64, 64);
                    t.verifyEqual(all(sd >= 0), all(lam > 0), sprintf('%s at %s', tf.kind, mat2str(x)));
                    t.verifyFalse(any(sd == 0), 'off the edges no side is zero');
                    nIn = nIn + all(lam > 0);  nOut = nOut + ~all(lam > 0);
                end
                t.verifyGreaterThan([nIn nOut], [5 5]);
            end
        end

        function aRATIONALCornerIsOnTwoEdgesExactly(t)
            tf = t.faces(strcmp({t.faces.kind}, 'rational'));
            for i = 1:3
                t.verifyEqual(tf.Tn(i,:) / tf.Td(i), tf.T(i,:), 'AbsTol', 1e-12);
                sd = triangleFaceSide(tf, tf.Tn(i,:), tf.Td(i));
                t.verifyEqual(sort(sd), [0 0 1], sprintf('corner %d', i));
            end
        end

        function theValueAtEachCornerIsF(t)
        % Fenchel equality at a contact point: <s0, m_i> - g_i(s0) = f(m_i). The face is the plane
        % through three points OF THE GRAPH of f, which is what makes it a chord.
            for k = 1:numel(t.faces)
                tf = t.faces(k);
                z = tf.T * tf.s - tf.c;
                t.verifyEqual(z, t.f.eval(tf.T), 'RelTol', 1e-10, 'AbsTol', 1e-8, tf.kind);
            end
        end

        function anUNBOUNDEDOrUnnamedTriangleIsRefusedByName(t)
            tr = t.tris(1);  tr.G = [1; 0];
            t.verifyError(@() triangleFace(tr), 'triangleFace:unbounded');
            tr = t.tris(1);  tr.Gn = tr.Gn(1:2,:);  tr.Gd = tr.Gd(1:2);
            t.verifyError(@() triangleFace(tr), 'triangleFace:noExactForms');
            t.verifyError(@() triangleFace(42), 'triangleFace:shape');
            t.verifyError(@() triangleFaceSide(t.faces(1), [1 2 3], 1), 'triangleFaceSide:point');
            t.verifyError(@() biconjAtoms(t.f, 'pairs'), 'biconjAtoms:mode');
        end

        function cuttingPutsTheTRIANGLEPlaneExactlyOnTAndNowhereElse(t)
        % A stand-in for the pairs' subdivision: one constant face per sign pattern of the cut lines
        % (every face then lies on one side of each, as envelopeSubdivisionAlg guarantees), valued
        % far ABOVE the plane, as a pair's value is inside T. After cutting, the evaluator must read
        % the plane inside T and the constant outside -- the max over owners would read the
        % constant everywhere if any piece inside T survived.
            rng(52);
            for k = 1:numel(t.faces)
                tf = t.faces(k);
                [~, L] = cutTriangles([], tf, 'lines');
                keys = arrayfun(@(i) ratQ.conic([0 0 0 L(i,:)]), 1:size(L,1), 'UniformOutput', false);
                flip = arrayfun(@(i) sign(keys{i}(find(keys{i}, 1))) * sign(L(i, find(L(i,:), 1))), ...
                                1:size(L,1));
                sub = struct('kind', {}, 'num', {}, 'den', {}, 'FP', {}, 'FB', {}, 'con', {}, 'alg', {});
                for pat = 0:2^size(L,1) - 1
                    sg = 2 * bitget(pat, 1:size(L,1)) - 1;          % the face's side of each line
                    con = [cell2mat(keys(:)), (sg .* flip).', zeros(size(L,1), 1)];
                    sub(end+1) = struct('kind', 'quad', 'num', [zeros(1,9) 1e6], 'den', 1, ...
                                        'FP', [], 'FB', 0, 'con', con, 'alg', []); %#ok<AGROW>
                end
                cut = cutTriangles(sub, tf);
                T = tf.T;  lo = min(T) - 0.5 * (max(T) - min(T));  hi = max(T) + 0.5 * (max(T) - min(T));
                nIn = 0;
                for r = 1:600
                    x = lo + rand(1,2) .* (hi - lo);
                    lam = [T.' ; 1 1 1] \ [x(:); 1];
                    if min(abs(lam)) < 1e-6, continue, end
                    want = 1e6;
                    if all(lam > 0), want = x * tf.s - tf.c;  nIn = nIn + 1; end
                    t.verifyEqual(subdivisionAlgEval(cut, x), want, 'RelTol', 1e-12, ...
                                  sprintf('%s at %s', tf.kind, mat2str(x)));
                end
                t.verifyGreaterThan(nIn, 20);
            end
        end

        function aSUBDIVISIONBuiltWithoutTheCutLinesIsCutAfterwards(t)
        % 2026-11-01: the cut lines need not be among the subdivision's curves. A face without a row
        % for one has its side decided (a polygon clear of the line, else exactly), and a face the
        % line crosses is split. Faces here: the two sides of a vertical line through BOTH
        % triangles, and one far away -- valued far above the planes. Both triangles are cut in
        % one call, so the pieces T1 leaves must still be cut by T2.
            T = [t.faces.T];  T = reshape(T, 3, 2, []);
            allT = [t.faces(1).T; t.faces(2).T];
            c = round(mean(allT(:,1)));
            mk = @(con) struct('kind', 'quad', 'num', [zeros(1,9) 1e6], 'den', 1, 'FP', [], 'FB', 0, ...
                               'con', con, 'alg', []);
            line = [0 0 0 1 0 -c];                           % x - c
            sub = [mk([line, 1, 0]), mk([line, -1, 0]), mk([0 0 0 1 0 1e4, -1, 0])];
            cut = cutTriangles(sub, t.faces);
            rng(53);
            lo = min(allT) - 2;  hi = max(allT) + 2;  nIn = [0 0];
            for r = 1:800
                x = lo + rand(1,2) .* (hi - lo);
                want = 1e6;
                for k = 1:2
                    Tk = T(:,:,k);  lam = [Tk.' ; 1 1 1] \ [x(:); 1];
                    if min(abs(lam)) < 1e-6, want = NaN; break, end
                    if all(lam > 0), want = x * t.faces(k).s - t.faces(k).c;  nIn(k) = nIn(k) + 1; end
                end
                if isnan(want) || abs(x(1) - c) < 1e-6, continue, end
                t.verifyEqual(subdivisionAlgEval(cut, x), want, 'RelTol', 1e-12, sprintf('at %s', mat2str(x)));
            end
            t.verifyGreaterThan(nIn, [20 20]);
            t.verifyEqual(subdivisionAlgEval(cut, [-2e4 0]), 1e6, 'the far face is kept whole');
            % the typed object reads the same: table 3 rows and the 'triangle' face (AlgCrv.FT)
            g = assembleAlgCrv(cut);
            X = lo + rand(400, 2) .* (hi - lo);
            t.verifyEqual(g.eval(X), subdivisionAlgEval(cut, X), 'RelTol', 1e-12);
            t.verifyEqual(nnz(g.FT), 1, 'the named triangle is a face of its own; the rational one is affine');
        end
    end
end

function A = conicArr(c)
% [a b c d e f] as the polyN array A(i+1,j+1) = coefficient of s1^i s2^j
    A = [c(6) c(5) c(3); c(4) c(2) 0; c(1) 0 0];
end
