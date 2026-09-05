classdef biconjQTest < matlab.unittest.TestCase
% Unit tests for biconjQ, the EXACT biconjugate f** = cl co f.
%
% BUCKET: fast (integer arithmetic only; no symbolic call anywhere).
%
% THE ORACLE IS THE DEFINITION OF A CONVEX ENVELOPE, in the three parts plqCheck.m also uses:
%   co f <= f with EQUALITY at the extreme points of the domain; co f is CONVEX; and co f is the
%   LARGEST convex minorant, so no affine function lying below f may exceed it anywhere. None of
%   the three mentions a lower hull, which is how biconjQ computes -- so they check the answer
%   rather than restate the method.

    methods (Test)

        function aNEEDLEIsItsOwnEnvelope(testCase)
        % A single point is convex, so co f = f: the value q(p) at p and +infinity elsewhere. The
        % domain is THIN and is carried by the H-form's equality side, the same machinery conj uses
        % for a point-supported conjugate.
            h = biconjQ(QuaPol([1 2], zeros(0,3), [0 0 0 0 0 0 0 0 0 5], zeros(0,2)));
            testCase.verifyEqual(h.nf, 1);
            testCase.verifyEqual(h.eval([1 2]), 5, 'AbsTol', 1e-12);
            testCase.verifyTrue(all(isinf(h.eval([1 2.5; 1.5 2; 0 0]))));
        end

        function aSEGMENTEnvelopesInOneDimension(testCase)
        % The problem is one-dimensional in the segment parameter, so the answer turns on the
        % curvature ALONG the segment and not on H in the plane. Both directions are pinned.
            X = [0 0; 0.5 0; 1 0; 2 0];

            % q = x^2/2 curves UP along the segment, so it is already convex there: co f = f.
            up = biconjQ(QuaPol([0 0; 2 0], [1 2 1], [0 0 0 0 1 0 0 0 0 0], [0 0]));
            testCase.verifyEqual(up.eval(X), X(:,1).^2/2, 'AbsTol', 1e-12);

            % q = -x^2/2 curves DOWN, so co f is the CHORD through (0,0) and (2,-2), i.e. -x.
            dn = biconjQ(QuaPol([0 0; 2 0], [1 2 1], [0 0 0 0 -1 0 0 0 0 0], [0 0]));
            testCase.verifyEqual(dn.eval(X), -X(:,1), 'AbsTol', 1e-12);
            testCase.verifyEqual(dn.eval([0 0; 2 0]), [0; -2], 'AbsTol', 1e-12, ...
                'the chord interpolates the endpoint values exactly');

            % and the domain really is the segment, not its line and not the plane
            testCase.verifyTrue(all(isinf(dn.eval([1 0.3; 3 0; -1 0]))), ...
                'off the line, and beyond either end, the envelope is +infinity');
        end

        function theDimensionCheckRunsBEFORETheConvexShortCircuit(testCase)
        % ORDERING IS LOAD-BEARING HERE, and the wrong order is silent about which case it is in:
        % an affine or PSD q on a low-dimensional domain IS convex, so the convex branch would
        % claim it and then try to read a FACE -- and such a mesh has nf = 0 and an empty F.
        % Measured while building this: the needle and the convex segment both raised noFace, while
        % the CONCAVE segment, which falls past the short-circuit, worked.
            needle = QuaPol([1 2], zeros(0,3), [0 0 0 0 0 0 0 0 0 5], zeros(0,2));
            convexSeg = QuaPol([0 0; 2 0], [1 2 1], [0 0 0 0 1 0 0 0 0 0], [0 0]);
            testCase.verifyEqual(biconjQ(needle).eval([1 2]), 5, 'AbsTol', 1e-12);
            testCase.verifyEqual(biconjQ(convexSeg).eval([1 0]), 0.5, 'AbsTol', 1e-12);
        end

        function aRAYDomainIsRefusedByName(testCase)
        % Where q curves down along a ray the envelope runs to -infinity, a correct answer with
        % nowhere to be stored -- the gap conjCPLQ records as convEnvUnbounded:minusInfinity.
            r = QuaPol([0 0; 1 0], [1 2 0], [0 0 0 0 -1 0 0 0 0 0], [0 0]);
            testCase.verifyError(@() biconjQ(r), 'PLQ:biconjQ:unbounded');
        end

        function aConvexFunctionIsItsOwnEnvelope(testCase)
        % co f = f when f is convex, so there is nothing to compute -- and getting this wrong is
        % the silent failure ALGORITHM.md warns about: biconj returning a non-convex f unchanged.
            f = biconjQTest.square([0 0 0 0, 2, 0, 2, -1, 3, 1]);   % strictly convex
            h = biconjQ(f);
            X = biconjQTest.insideSquare(200);
            testCase.verifyEqual(h.eval(X), f.eval(X), 'AbsTol', 1e-12);
        end

        function anAffineFunctionIsItsOwnEnvelopeToo(testCase)
        % The H = 0 edge of both branches: affine is convex AND concave, so whichever route runs
        % must return f unchanged.
            f = biconjQTest.square([0 0 0 0, 0, 0, 0, -1, 2, 3]);
            h = biconjQ(f);
            X = biconjQTest.insideSquare(200);
            testCase.verifyEqual(h.eval(X), f.eval(X), 'AbsTol', 1e-12);
            testCase.verifyEqual(h.nf, 1, 'an affine function needs exactly one piece');
        end

        function theEnvelopeOfAConcaveQuadraticInterpolatesItsVertexValues(testCase)
        % The defining property of this case: a concave function on a polytope has its whole
        % envelope determined by the VERTEX values, so co f must match f at every corner and lie
        % strictly below it in between.
            f = biconjQTest.square([0 0 0 0, -1, 0, -1, 0, 0, 0]);   % -(x^2+y^2)/2
            h = biconjQ(f);
            V = [0 0; 1 0; 1 1; 0 1];
            testCase.verifyEqual(h.eval(V), f.eval(V), 'AbsTol', 1e-12, ...
                'equality at every extreme point');

            X = biconjQTest.insideSquare(300);
            testCase.verifyLessThanOrEqual(h.eval(X), f.eval(X) + 1e-12, 'co f <= f');
            testCase.verifyLessThan(h.eval([0.5 0.5]), f.eval([0.5 0.5]) - 1e-9, ...
                'and strictly below in the interior, since q is strictly concave there');
        end

        function theEnvelopeIsConvexAndIsTheLargestConvexMinorant(testCase)
        % Two definition-level properties that no cell-by-cell check implies. The second is what
        % separates the true envelope from any convex function that merely sits below f.
            f = biconjQTest.polygon([0 0; 3 0; 4 3; 1 4], [0 0 0 0, -1, 0, -1, 0, 0, 0]);
            h = biconjQ(f);
            V = [0 0; 3 0; 4 3; 1 4];

            rng(20260904);
            A = biconjQTest.inHull(V, 200);  B = biconjQTest.inHull(V, 200);
            testCase.verifyLessThanOrEqual(h.eval((A+B)/2), (h.eval(A) + h.eval(B))/2 + 1e-9, ...
                'the envelope must be convex');

            qv = f.eval(V);
            X  = biconjQTest.inHull(V, 300);
            hx = h.eval(X);
            for t = 1:25
                g = randn(1,2);
                b = min(qv - V*g.');                 % the highest affine minorant with slope g
                testCase.verifyLessThanOrEqual(X*g.' + b, hx + 1e-9, ...
                    'an affine function below f may never exceed co f');
            end
        end

        function aCollinearVertexIsDominatedAndTheEnvelopeSaysSo(testCase)
        % NOT a defect, and it cost a diagnosis to establish: a vertex lying ON the segment between
        % its neighbours is not an extreme point, so for a strictly concave q its lifted point sits
        % strictly ABOVE the chord and is not on the lower hull. co f is then strictly BELOW q
        % there. The dual of conjQ's aDominatedVertexContributesNoCellAtAll.
            V = [0 0; 1 0; 2 0; 2 2; 0 2];           % (1,0) is collinear with (0,0) and (2,0)
            f = biconjQTest.polygon(V, [0 0 0 0, -2, 0, -2, 0, 0, 0]);
            h = biconjQ(f);
            ext = [1 3 4 5];
            testCase.verifyEqual(h.eval(V(ext,:)), f.eval(V(ext,:)), 'AbsTol', 1e-12, ...
                'equality at the four EXTREME points');
            testCase.verifyLessThan(h.eval(V(2,:)), f.eval(V(2,:)) - 1e-9, ...
                'and strictly below at the collinear one, which the hull does not use');
        end

        function theEnvelopeCoefficientsAreExactRationals(testCase)
        % The point of the exact route: the stored coefficients are integers over integers, not
        % doubles that happen to be close.
            f = biconjQTest.square([0 0 0 0, -1, 0, -1, 0, 0, 0]);
            h = biconjQ(f);
            testCase.verifyEqual(h.fN, round(h.fN), 'AbsTol', 0);
            testCase.verifyTrue(all(h.fD > 0));
            testCase.verifyTrue(all(h.fN(:,5) == 0 & h.fN(:,6) == 0 & h.fN(:,7) == 0), ...
                'the envelope of a concave q is piecewise AFFINE -- no curvature may appear');
        end

        % ---- what is refused, by name ------------------------------------------------------

        function aNONCONVEXFaceEnvelopesOverItsCONVEXHULL(testCase)
        % THE ENVELOPE'S DOMAIN IS conv(P), NOT P -- co f must be convex, so its domain is convex.
        % For a non-convex face those differ, and using the face's own edges produced a region that
        % was not even the face: biconjQ returned +Inf at (2,0), (2,1), (1,2) and (1.5,0.5), every
        % one of which is IN the L and where the envelope is finite.
        %
        % The FACE FUNCTIONS needed no change -- the lower hull of the lifted vertices was already
        % right, reflex vertices included -- so only the domain was wrong. Unlike the conjugate this
        % could NOT have been fixed by splitting: the convex hull of a union is not the union of
        % hulls, so the envelope does not decompose over a subdivision the way conj does.
            V = [0 0; 2 0; 2 1; 1 1; 1 2; 0 2];         % an L, reflex at (1,1)
            m = size(V,1);
            f = biconjQTest.polygon(V, [0 0 0 0, -1, 0, -1, 0, 0, 0]);
            [~, isConv] = f.orderEdges(1);
            testCase.verifyFalse(logical(isConv), 'the fixture really is non-convex');

            h = biconjQ(f);

            % Every vertex of the L, reflex included, is a point of the domain with a finite value.
            testCase.verifyTrue(all(isfinite(h.eval(V))), ...
                'no vertex of the face may fall outside the envelope''s domain');
            % Points of conv(L) that are NOT in L are still in the envelope's domain.
            testCase.verifyTrue(all(isfinite(h.eval([1.2 1.2; 1.4 1.4]))), ...
                'the domain is the convex HULL, so it covers the reflex notch');
            % And the envelope lies below f wherever f is finite, with equality at the hull corners.
            qv = -0.5*sum(V.^2, 2);
            testCase.verifyLessThanOrEqual(h.eval(V), qv + 1e-9, 'co f <= f at every vertex');
            ext = [1 2 3 5 6];                          % (1,1) is the reflex one, not a hull corner
            testCase.verifyEqual(h.eval(V(ext,:)), qv(ext), 'AbsTol', 1e-9, ...
                'and equality at the corners of the hull');
        end

        function aNonConvexFaceIsRefusedWhenTheAnswerWouldHaveToCarryIt(testCase)
        % The other side of the same coin. When f is CONVEX the envelope is f itself, so the answer
        % has to carry f's own non-convex region -- which an intersection of half-planes cannot
        % express. Refused by name rather than silently replaced by the hull.
            V = [0 0; 2 0; 2 1; 1 1; 1 2; 0 2];
            f = biconjQTest.polygon(V, [0 0 0 0, 1, 0, 1, 0, 0, 0]);   % convex on a non-convex face
            testCase.verifyError(@() biconjQ(f), 'PLQ:biconjQ:nonConvexPiece');
        end

        function theBILINEAREnvelopeOverABoxIsExactlyMcCORMICK(testCase)
        % THE case ALGORITHM.md singles out, and it now falls out of the concave branch rather than
        % needing one of its own. q = xy on the unit square has H = [0 1; 1 0], which is genuinely
        % INDEFINITE -- so the old H-negative-semidefinite guard refused it -- but d'Hd = 0 along
        % BOTH edge directions, so q is affine (hence concave) along every edge and its envelope is
        % the lower hull of the four corner values. That hull is max(0, x+y-1): the McCormick /
        % Al-Khayyal-Falk envelope, written here independently of anything biconjQ computes.
            f = biconjQTest.square([0 0 0 0, 0, 1, 0, 0, 0, 0]);
            h = biconjQ(f);
            testCase.verifyEqual(h.nf, 2, 'McCormick has exactly two affine pieces');
            rng(20260904);
            X = rand(400,2);
            testCase.verifyEqual(h.eval(X), max(0, X(:,1) + X(:,2) - 1), 'AbsTol', 1e-12);
            V = [0 0; 1 0; 1 1; 0 1];
            testCase.verifyEqual(h.eval(V), [0; 0; 1; 0], 'AbsTol', 1e-12, ...
                'and it interpolates xy at all four corners');
        end

        function anINDEFINITEButEdgeConcavePieceStillHasAVertexEnvelope(testCase)
        % The criterion in general, not just for xy. What the envelope needs is concavity along the
        % polygon's EDGE DIRECTIONS, which is weaker than concavity of H: on a triangle, q - L
        % vanishes at the corners and is concave along each edge, so q >= L on the boundary, and
        % with H not positive semidefinite q - L has no interior minimum -- so the minimum over the
        % closed triangle is on that boundary. Checked here against the DEFINITION rather than
        % against another implementation.
            V = [0 0; 2 0; 0 2];
            H = [-1 2; 2 -1];                       % det = -3: indefinite
            for j = 1:3
                d = V(mod(j,3)+1,:) - V(j,:);
                testCase.verifyLessThanOrEqual(d*H*d.', 0, 'the fixture must be edge-concave');
            end
            testCase.verifyLessThan(det(H), 0, 'and genuinely indefinite');

            f = biconjQTest.polygon(V, [0 0 0 0, H(1,1), H(1,2), H(2,2), 1, -1, 2]);
            h = biconjQ(f);

            qv = 0.5*sum((V*H).*V,2) + V*[1;-1] + 2;
            testCase.verifyEqual(h.eval(V), qv, 'AbsTol', 1e-9, 'equality at the corners');
            X = biconjQTest.inHull(V, 400);
            qx = 0.5*sum((X*H).*X,2) + X*[1;-1] + 2;
            testCase.verifyLessThanOrEqual(h.eval(X), qx + 1e-9, 'co f <= f');
            A = biconjQTest.inHull(V, 200);  B = biconjQTest.inHull(V, 200);
            testCase.verifyLessThanOrEqual(h.eval((A+B)/2), (h.eval(A) + h.eval(B))/2 + 1e-9);
        end

        function aPieceThatCurvesUPWARDAlongAnEdgeIsRefusedByName(testCase)
        % The condition is EDGE-concavity, so what is refused is an edge along which q curves
        % upward -- there the vertex values no longer determine the envelope. This is the trigger
        % for AlgAlg: the first face that cannot be written rationally appears here, since an
        % affine cell of a general f** names a dual vertex of degree up to 4.
            f = biconjQTest.square([0 0 0 0, 1, 0, -1, 0, 0, 0]);   % +1 along (1,0)
            testCase.verifyError(@() biconjQ(f), 'PLQ:biconjQ:notImplemented');
        end

        function aMultiPieceNonConvexInputIsRefusedByName(testCase)
        % The envelope COUPLES pieces -- the convex hull of a union is not determined piecewise --
        % so unlike the conjugate it is not a fold, and refusing is the honest answer.
            V = [0 0; 1 0; 1 1; 0 1];
            E = [1 2 1; 2 3 1; 3 4 1; 4 1 1; 1 3 1];
            F = [1 0; 1 0; 2 0; 2 0; 1 2];
            f = QuaPol(V, E, [0 0 0 0, -1, 0, -1, 0, 0, 0; 0 0 0 0, -2, 0, -2, 0, 0, 0], F);
            testCase.verifyError(@() biconjQ(f), 'PLQ:biconjQ:notImplemented');
        end

        function anInexactInputIsRefusedBeforeAnythingElse(testCase)
            testCase.verifyError(@() biconjQ(QuaPol([sqrt(2) 0 1 0 0 0])), 'PLQ:QuaPol:notExact');
        end

        function aFalseConvexityAssertionIsRefusedLoudly(testCase)
        % ALGORITHM.md's rule: the free necessary condition is still checked, because the failure
        % it prevents is silent -- biconj returning a non-convex f as its own convex envelope.
            f = biconjQTest.square([0 0 0 0, 1, 0, -1, 0, 0, 0]);
            f.fIsConvex = true;
            testCase.verifyError(@() biconjQ(f), 'PLQ:biconjQ:notConvex');
        end
    end

    methods (Static)
        function f = square(fc)
            f = QuaPol([0 0; 1 0; 1 1; 0 1], [1 2 1; 2 3 1; 3 4 1; 4 1 1], fc, ...
                       [1 0; 1 0; 1 0; 1 0]);
        end

        function f = polygon(V, fc)
            m = size(V,1);
            f = QuaPol(V, [(1:m).', [2:m,1].', ones(m,1)], fc, [ones(m,1), zeros(m,1)]);
        end

        function X = insideSquare(n)
            rng(20260904);
            X = rand(n,2);
        end

        function X = inHull(V, n)
            m = size(V,1);  X = zeros(n,2);
            for i = 1:n
                w = rand(m,1);  w = w/sum(w);
                X(i,:) = w.' * V;
            end
        end
    end
end
