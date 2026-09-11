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

        function anUnboundedRefusalSaysWHICHKindItIs(testCase)
        % Two different situations hid behind one identifier. Along a recession direction d, q
        % behaves as t*<grad q, d> + (t^2/2) d'Hd, so the envelope runs to -infinity exactly when
        % some recession direction has d'Hd < 0, or d'Hd = 0 with a negative linear rate. Otherwise
        % q is bounded below on its own recession cone and the envelope is FINITE -- a real answer,
        % simply not computed yet.
        %
        % The distinction is not cosmetic: it says which refusals are REPRESENTATIONAL (a correct
        % -infinity with nowhere to go) and which are a genuine gap. Measured over the fixture
        % corpus: all 8 fixtures the legacy envelope answers and this one refuses are the FINITE
        % kind, so all 8 are reachable.
            % strictly concave on a cone: q runs DOWN along both rays, so co f is -infinity
            down = QuaPol([0 0; 1 0; 0 1], [1 2 0; 1 3 0], ...
                          [0 0 0 0, -2, 0, -2, 0, 0, 0], [1 0; 0 1]);
            testCase.verifyError(@() biconjQ(down), 'PLQ:biconjQ:minusInfinity');

            % An AFFINE q on the same cone is CONVEX, so it is answered rather than refused --
            % co f = f. Worth pinning, because it was the first fixture tried here and the refusal
            % expected of it was simply wrong.
            aff = QuaPol([0 0; 1 0; 0 1], [1 2 0; 1 3 0], ...
                         [0 0 0 0, 0, 0, 0, 1, 1, 0], [1 0; 0 1]);
            testCase.verifyEqual(biconjQ(aff).nf, 1, 'an affine q is its own envelope');

            % Non-convex, but BOUNDED BELOW on its recession cone: H = diag(1,-1) is indefinite,
            % yet on cone((1,0),(1,1)) the curvature is 1 - t^2 >= 0 throughout. So the envelope is
            % FINITE and the refusal must be an honest "not computed yet", not minusInfinity.
            %
            % CHANGED TWICE on 2026-09-05, and it is worth saying why rather than only that. It
            % first asserted PLQ:biconjQ:unbounded, meaning "finite but not computed"; the
            % unbounded envelope removed that. It then asserted shearedRecession; the minimum-rate
            % lifted ray removed that too. What is left is the REAL obstruction and it was there
            % all along: the ray (1,0) has d'Hd = +1, so q grows QUADRATICALLY along it, the piece
            % is not edge-concave, and its own envelope is the unimplemented single-piece case.
            % The distinction this test exists for -- a correct -infinity against an honest gap --
            % is unchanged through all three.
            fin = QuaPol([0 0; 1 0; 1 1], [1 2 0; 1 3 0], ...
                         [0 0 0 0, 1, 0, -1, 0, 0, 0], [1 0; 0 1]);
            testCase.verifyError(@() biconjQ(fin), 'PLQ:biconjQ:notImplemented');
        end

        function anUNBOUNDEDPieceGetsTheLIFTEDHULLWithItsRECESSIONRaysInIt(testCase)
        % The unbounded envelope, on a fixture whose answer is known in closed form rather than
        % only checkable: q = -x^2 on the half-strip {0 <= x <= 1, y >= 0}.
        %
        % q is concave in x and FLAT in y, and the strip recedes along (0,1), which is in the
        % kernel of H = diag(-2,0) -- so the lifted set recedes in the fixed direction (0,1,0) and
        % the hull is conv(lifted vertices) + cone(lifted rays). Along x the envelope is the chord
        % from q(0)=0 to q(1)=-1 and it cannot depend on y, so
        %
        %       co f(x,y) = -x   on the whole half-strip.
        %
        % Getting the recession rays WRONG shows up here twice over: drop them from the hull and
        % the domain collapses to the bounded segment y = 0, so co f returns +infinity at every
        % point with y > 0; keep them but lift them to the wrong height and the face tilts in y.
            f = QuaPol([0 0; 1 0; 0 1; 1 1], [1 2 1; 1 3 0; 2 4 0], ...
                       [0 0 0 0, -2, 0, 0, 0, 0, 0], [1 0; 0 1; 1 0]);
            h = biconjQ(f);

            X = [0 0; 1 0; 0.5 0; 0.25 3; 0.5 10; 1 7; 0 100];
            testCase.verifyEqual(h.eval(X), -X(:,1), 'AbsTol', 1e-12, ...
                'co f = -x, at y = 0 and far out along the recession direction alike');

            % and it really is the ENVELOPE of this f, not some other convex function
            fx = f.eval(X);
            testCase.verifyLessThanOrEqual(h.eval(X), fx + 1e-12, 'co f <= f');
            testCase.verifyEqual(h.eval([0 0; 1 0]), f.eval([0 0; 1 0]), 'AbsTol', 1e-12, ...
                'with equality at the two extreme points of the strip');

            % off the domain it is still +infinity: the hull grew a recession cone, not a plane
            testCase.verifyTrue(all(isinf(h.eval([-1 5; 2 5; 0.5 -3]))), ...
                'conv(dom f) is the half-strip itself, so x < 0, x > 1 and y < 0 stay +infinity');
        end

        function theUNBOUNDEDHullIsCHECKEDAgainstTheDefinitionNotOnlyAKnownValue(testCase)
        % The same machinery on a fixture with no closed form to lean on: the third defining
        % property, co f(sum lam_i x_i) <= sum lam_i f(x_i), which holds for every convex minorant
        % and is the one that catches an envelope sitting too HIGH. It needs no vertex enumeration,
        % so unlike a vertex-equality check it is valid on an unbounded domain.
            rng(20260905);
            f = QuaPol([0 0; 2 0; 0 1; 2 1], [1 2 1; 1 3 0; 2 4 0], ...
                       [0 0 0 0, -2, 0, 0, 1, 3, -1], [1 0; 0 1; 1 0]);
            h = biconjQ(f);

            X = [2*rand(300,1), 12*rand(300,1)];          % inside the strip by construction
            fx = f.eval(X);  hx = h.eval(X);
            testCase.verifyTrue(all(isfinite(fx)) && all(isfinite(hx)), ...
                'the sample and the envelope are both finite on the strip');
            testCase.verifyLessThanOrEqual(hx, fx + 1e-9, 'co f <= f');

            for t = 1:200
                k = randi([2 4]);
                idx = randi(numel(fx), k, 1);
                lam = rand(k,1);  lam = lam/sum(lam);
                z = lam.' * X(idx,:);
                testCase.verifyLessThanOrEqual(h.eval(z), lam.'*fx(idx) + 1e-9, ...
                    'co f is a minorant of every convex combination of values of f');
            end
        end

        function theDIVERGENCETestReadsTheRateONThePieceNotAtTheOrigin(testCase)
        % Along a recession direction d with d'H d = 0, q(x + t d) is AFFINE in t with rate
        % (H x + L)'d -- which depends on x whenever H*d is nonzero, i.e. whenever H is indefinite.
        % The test used to read L'd, the rate at the ORIGIN, and so missed every piece whose apex
        % is not there.
        %
        % q = x*y on the TRANSLATED quadrant {x >= 0, y >= -5}. Both recession directions have
        % d'H d = 0 and L = 0, so the curvature and the origin-rate tests both pass it; the cone
        % test has a = c = 0 and b = +1, so that passes it too. But f(t,-5) = -5t, so f is
        % unbounded BELOW on its own domain and co f is -infinity -- which needs no envelope
        % machinery to see, and is what makes this a WITNESS rather than an opinion.
            f = QuaPol([0 -5; 1 -5; 0 -4], [1 2 0; 1 3 0], ...
                       [0 0 0 0, 0, 1, 0, 0, 0, 0], [1 0; 0 1]);
            testCase.verifyEqual(f.eval([1000 -5]), -5000, 'AbsTol', 1e-9, ...
                'the witness: f marches to -infinity along its own recession ray');
            testCase.verifyError(@() biconjQ(f), 'PLQ:biconjQ:minusInfinity', ...
                'so the refusal must be minusInfinity, not the finite-but-uncomputed one');
        end

        function theSAMECaseAtTheORIGINIsANSWEREDBecauseItIsBoundedBelow(testCase)
        % The control, so the test above is about the TRANSLATION and not about x*y or about cones.
        % On the first quadrant x*y >= 0, and co f is not merely finite -- it is ZERO, which is the
        % smallest instance of a lifted ray whose height is a MINIMUM rather than a constant.
        %
        % r(x,(1,0)) = y and r(x,(0,1)) = x, both minimised at the apex (0,0) where they vanish. So
        % the lifted rays are (1,0,0) and (0,1,0), the hull is the quadrant at height 0, and
        % co f = 0. Reading the rate at any other point would tilt the answer upward, and reading
        % L'd -- which is also 0 here -- would be right by accident, which is what the shifted
        % fixture below is for.
            f = QuaPol([0 0; 1 0; 0 1], [1 2 0; 1 3 0], ...
                       [0 0 0 0, 0, 1, 0, 0, 0, 0], [1 0; 0 1]);
            testCase.verifyGreaterThanOrEqual(f.eval([1000 3; 7 2; 0 9]), 0, ...
                'x*y is nonnegative on the first quadrant');
            h = biconjQ(f);
            X = [0 0; 1 0; 0 1; 3 4; 100 0.5; 0.25 80];
            testCase.verifyEqual(h.eval(X), zeros(size(X,1),1), 'AbsTol', 1e-12, ...
                'co f = 0 on the whole quadrant');
            testCase.verifyLessThanOrEqual(h.eval(X), f.eval(X) + 1e-12, 'co f <= f');
        end

        function aSHIFTEDApexMakesTheLiftedRayHeightNONZERO(testCase)
        % The fixture that separates the minimum-rate rule from the old constant-L'd one, since
        % L = 0 in both and only the apex moves. q = x*y on {x >= 0, y >= 2}: r(x,(1,0)) = y is
        % minimised at the apex (0,2) where it is 2, and r(x,(0,1)) = x is minimised there at 0.
        % So the lifted rays are (1,0,2) and (0,1,0), and
        %
        %       co f(x,y) = 2x   on {x >= 0, y >= 2}.
        %
        % L'd would have given height 0 for BOTH rays and hence co f = 0 -- a strict underestimate,
        % and one that is still a convex minorant, so only a tightness check catches it. That is
        % why equality along y = 2 is asserted and not just co f <= f.
            f = QuaPol([0 2; 1 2; 0 3], [1 2 0; 1 3 0], ...
                       [0 0 0 0, 0, 1, 0, 0, 0, 0], [1 0; 0 1]);
            h = biconjQ(f);
            X = [0 2; 1 2; 5 2; 0.5 2; 3 9; 0 40; 7 100];
            testCase.verifyEqual(h.eval(X), 2*X(:,1), 'AbsTol', 1e-12, 'co f = 2x');
            testCase.verifyLessThanOrEqual(h.eval(X), f.eval(X) + 1e-12, 'co f <= f');
            testCase.verifyEqual(h.eval([0 2; 1 2; 5 2]), f.eval([0 2; 1 2; 5 2]), ...
                'AbsTol', 1e-12, 'and it is TIGHT along y = 2, where f is already affine');
        end

        function aCHAINOfEdgesIsEnvelopedOverIts2DHull(testCase)
        % A chain of two non-collinear segments is NOT thin: its convex hull is a triangle, so co f
        % is an ordinary lower-hull answer over a 2-D domain. This used to be refused as
        % `ambiguousChain`, on the stated grounds that several function rows with no F left it
        % undecided which edge took which -- but for a chain ROW j IS EDGE j's function, which is
        % what QuaPol.eval already reads, so there was never an ambiguity.
        %
        % The fixture has a closed form. Edge 1 is (0,0)-(2,0) carrying -x^2, which is CONCAVE
        % along itself, so its lifted arc contributes only its endpoints: (0,0,0) and (2,0,-4).
        % Edge 2 is (0,0)-(0,2) carrying 0, contributing (0,2,0). Three lifted points span one
        % plane, so
        %
        %       co f(x,y) = -2x   on conv{(0,0), (2,0), (0,2)}.
        %
        % The interior point (0.5,0.5) is the load-bearing one: f is +infinity there, because it is
        % off the chain, while co f is finite -- which is what makes this a 2-D answer and not a
        % thin one.
            f = biconjQTest.chainL();
            testCase.verifyTrue(isinf(f.eval([0.5 0.5])), ...
                'f is +infinity off the chain');

            h = biconjQ(f);
            X = [0 0; 2 0; 0 2; 1 0; 0.5 0.5; 0.25 1.5; 1 0.5];
            testCase.verifyEqual(h.eval(X), -2*X(:,1), 'AbsTol', 1e-12, ...
                'co f = -2x on the triangle');
            testCase.verifyEqual(h.eval([0 0; 2 0; 0 2]), [0; -4; 0], 'AbsTol', 1e-12, ...
                'and it interpolates the three lifted vertices exactly');

            % co f <= f wherever f is finite, i.e. on the chain itself
            C = [0 0; 0.5 0; 1 0; 2 0; 0 0.5; 0 2];
            testCase.verifyLessThanOrEqual(h.eval(C), f.eval(C) + 1e-12, 'co f <= f on the chain');

            % and outside the hull it is still +infinity
            testCase.verifyTrue(all(isinf(h.eval([-1 0; 0 -1; 2 2]))), ...
                'the envelope lives on conv(dom f), not on the whole plane');
        end

        function aCHAINThatCurvesUPWARDAlongAnEdgeIsRefusedByName(testCase)
        % The companion, and the reason the previous fixture uses -x^2. A chain edge carrying +x^2
        % has a lifted arc that dips BELOW its own chord, so the vertex values no longer determine
        % the envelope -- the same obstruction as a face edge, and the same refusal.
        %
        % Every chain in examples2 is of this kind (they carry x^2 or y^2 along their edges), which
        % is why the corpus alone cannot exercise the case above.
            f = QuaPol([0 0; 2 0; 0 2], [1 2 1; 1 3 1], ...
                       [0 0 0 0, 2, 0, 0, 0, 0, 0; 0 0 0 0, 0, 0, 0, 0, 0, 0], [0 0; 0 0]);
            testCase.verifyError(@() biconjQ(f), 'PLQ:biconjQ:notImplemented');

            P = QuaPol.examples2();
            for k = [5 6 7 10]
                testCase.verifyError(@() biconjQ(P{k}), 'PLQ:biconjQ:notImplemented', ...
                    sprintf('examples2{%d} carries x^2 or y^2 along its edges', k));
            end
        end

        function aCOLLINEARChainIsTheONEDIMENSIONALHullAlongItsOwnSpan(testCase)
        % A collinear chain's hull is still a SEGMENT, so the answer is thin: the 1-D lower convex
        % hull of the lifted vertices, piecewise affine, every cell carrying the line as an H-form
        % EQUALITY. thinEnvelope's single segment is the one-piece case of exactly this.
        %
        % Fixture, whose hull can be read off by hand: (-1,0)-(0,0) carries -x^2 and (0,0)-(2,0)
        % carries 0, so the lifted points are (-1,-1), (0,0) and (2,0) in (t,z). The middle one
        % lies ABOVE the chord from (-1,-1) to (2,0) -- that chord is at -2/3 when t = 0 -- so it
        % is dominated and the hull is the single chord
        %
        %       co f(x,0) = -1 + (x+1)/3   on  -1 <= x <= 2.
            f = QuaPol([-1 0; 0 0; 2 0], [1 2 1; 2 3 1], ...
                       [0 0 0 0, -2, 0, 0, 0, 0, 0; 0 0 0 0, 0, 0, 0, 0, 0, 0], [0 0; 0 0]);
            h = biconjQ(f);
            testCase.verifyEqual(h.nf, 1, 'the middle vertex is dominated, so one cell');

            X = [-1 0; -0.5 0; 0 0; 1 0; 2 0];
            testCase.verifyEqual(h.eval(X), -1 + (X(:,1)+1)/3, 'AbsTol', 1e-12, ...
                'co f is the chord');
            testCase.verifyLessThanOrEqual(h.eval(X), f.eval(X) + 1e-12, 'co f <= f');
            testCase.verifyEqual(h.eval([-1 0; 2 0]), f.eval([-1 0; 2 0]), 'AbsTol', 1e-12, ...
                'with equality at the two extreme points of the span');

            testCase.verifyTrue(all(isinf(h.eval([0 1; -0.5 0.2; 3 0]))), ...
                'the domain stays THIN: off the line, and beyond either end, +infinity');
        end

        function aCOLLINEARChainWhoseHullHasTWOPiecesGetsTWOCells(testCase)
        % The multi-cell path. Two affine edges meeting in a V: -2x on [0,1] and 2x-4 on [1,2],
        % so the lifted points are (0,0), (1,-2), (2,0) and ALL THREE are on the lower hull. f is
        % already convex, so co f = f, and the answer must be two cells rather than one chord.
            f = QuaPol([0 0; 1 0; 2 0], [1 2 1; 2 3 1], ...
                       [0 0 0 0, 0, 0, 0, -2, 0,  0; ...
                        0 0 0 0, 0, 0, 0,  2, 0, -4], [0 0; 0 0]);
            h = biconjQ(f);
            testCase.verifyEqual(h.nf, 2, 'both hull pieces survive');
            X = [0 0; 0.5 0; 1 0; 1.5 0; 2 0];
            testCase.verifyEqual(h.eval(X), f.eval(X), 'AbsTol', 1e-12, ...
                'f is already convex along the span, so co f = f');
        end

        function aCOLLINEARChainCurvingUPWARDIsStillRefused(testCase)
        % examples2{4}, which carries x^2 along both its edges: the lifted arcs dip below their own
        % chords, so the vertex values do not determine the envelope, exactly as in the 2-D case.
            P = QuaPol.examples2();
            testCase.verifyError(@() biconjQ(P{4}), 'PLQ:biconjQ:notImplemented');
        end

        function biconjQREFUSESANonQuaPolByName(testCase)
            testCase.verifyError(@() biconjQ(42), 'PLQ:biconjQ:input');
            testCase.verifyError(@() biconjQ('x'), 'PLQ:biconjQ:input');
        end

        function aDirectionINSIDETheConeCanDivergeWhileBothGeneratorsDoNot(testCase)
        % The cone test, which exists because checking the two GENERATORS is not enough. On
        % cone(d1,d2) the curvature is a*l1^2 + 2b*l1*l2 + c*l2^2, nonnegative for every l >= 0
        % exactly when a >= 0, c >= 0 and (b >= 0 or b^2 <= a*c).
        %
        % q = x*y on the FOURTH quadrant, cone((1,0),(0,-1)): both generators have d'Hd = 0, so
        % neither is refused on its own, and both rates vanish at the apex. But b = -1 < 0 with
        % b^2 - a*c = 1 > 0, and indeed x*y -> -infinity along (1,-1). The witness needs no
        % machinery: f(t,-t) = -t^2.
            f = QuaPol([0 0; 1 0; 0 -1], [1 2 0; 1 3 0], ...
                       [0 0 0 0, 0, 1, 0, 0, 0, 0], [0 1; 1 0]);
            testCase.verifyEqual(f.eval([10 -10]), -100, 'AbsTol', 1e-9, ...
                'x*y marches to -infinity along the cone''s interior direction (1,-1)');
            testCase.verifyError(@() biconjQ(f), 'PLQ:biconjQ:minusInfinity');
        end

        function aCHAINWithARAYEdgeIsEnvelopedOverItsUnboundedHull(testCase)
        % The chain gather's ray branch. A chain of a SEGMENT (0,0)-(2,0) and a RAY (0,0)->(0,2),
        % both carrying 0: not collinear, so the hull is 2-D, and unbounded, so it needs the lifted
        % recession ray as well as the lifted vertices. co f = 0 on that hull.
            f = QuaPol([0 0; 2 0; 0 2], [1 2 1; 1 3 0], ...
                       [0 0 0 0, 0, 0, 0, 0, 0, 0; 0 0 0 0, 0, 0, 0, 0, 0, 0], [0 0; 0 0]);
            h = biconjQ(f);
            X = [0 0; 1 0; 2 0; 0 5; 1 40; 0.5 1];
            testCase.verifyEqual(h.eval(X), zeros(size(X,1),1), 'AbsTol', 1e-12, ...
                'co f = 0 on the hull, far out along the ray included');
            testCase.verifyTrue(all(isinf(h.eval([-1 1; 3 -1]))), ...
                'and +infinity outside it');
        end

        function aRAYDomainIsRefusedByName(testCase)
        % Where q curves down along a ray the envelope runs to -infinity, a correct answer with
        % nowhere to be stored -- the gap conjCPLQ records as convEnvUnbounded:minusInfinity.
            r = QuaPol([0 0; 1 0], [1 2 0], [0 0 0 0 -1 0 0 0 0 0], [0 0]);
            testCase.verifyError(@() biconjQ(r), 'PLQ:biconjQ:unbounded');
        end

        function noEnvelopeInTheFixtureSetHasOverlappingCells(testCase)
        % The invariant conjQ's half-plane defects broke, asked of the ENVELOPE too rather than
        % assumed: no two cells may overlap while carrying different functions. Such a mesh does
        % not define a function at all -- eval resolves it by first match, so the answer depends on
        % cell ORDER, and the failure is silent. biconjQ shares assembleQuaConCells with conjQ, so
        % it is exposed to the same class. Measured across the shapes below: 0 of 12 inconsistent
        % (.claude/biconj-consistency-sweep.m runs the wider set).
            Hs = {[2 0;0 2], [1 0;0 0], [-2 0;0 -2], [-1 0;0 0], [0 0;0 0], [0 1;1 0]};
            for i = 1:numel(Hs)
                H = Hs{i};
                h = biconjQ(biconjQTest.square([0 0 0 0, H(1,1), H(1,2), H(2,2), 1, -1, 2]));
                testCase.verifyTrue(checkQuaConConsistent(h, 400), ...
                    sprintf('square with H = %s', mat2str(H)));
            end
            % a non-convex face, a multi-piece input, and both thin domains
            h = biconjQ(biconjQTest.polygon([0 0; 2 0; 2 1; 1 1; 1 2; 0 2], ...
                                            [0 0 0 0, -1, 0, -1, 0, 0, 0]));
            testCase.verifyTrue(checkQuaConConsistent(h, 400), 'non-convex L');

            V = [0 0; 1 0; 1 1; 0 1];
            h = biconjQ(QuaPol(V, [1 2 1; 2 3 1; 3 4 1; 4 1 1; 1 3 1], ...
                [0 0 0 0 -2 0 0 0 0 0; 0 0 0 0 0 0 -2 0 0 0], [1 0;1 0;2 0;2 0;2 1]));
            testCase.verifyTrue(checkQuaConConsistent(h, 400), 'two triangles');

            testCase.verifyTrue(checkQuaConConsistent( ...
                biconjQ(QuaPol([1 2], zeros(0,3), [0 0 0 0 0 0 0 0 0 5], zeros(0,2))), 200), ...
                'needle');
            testCase.verifyTrue(checkQuaConConsistent( ...
                biconjQ(QuaPol([0 0; 2 0], [1 2 1], [0 0 0 0 -1 0 0 0 0 0], [0 0])), 200), ...
                'segment');
        end

        function aFULLPLANEConvexQuadraticIsItsOwnEnvelope(testCase)
        % QuaPol.energy -- (x^2+y^2)/2 on all of R^2 -- is the most elementary convex input there
        % is, and it used to RAISE noFace: quaPolAsQuaCon looks for a FACE, and a full-domain mesh
        % has none. co f = f, with no constraints at all.
            f = QuaPol.energy();
            h = biconjQ(f);
            testCase.verifyEqual(h.nf, 1);
            testCase.verifyEqual(h.ne, 0, 'the whole plane needs no bounding curve');
            X = [0 0; 1 1; -1 2; 0.5 -0.5; 7 -3];
            testCase.verifyEqual(h.eval(X), f.eval(X), 'AbsTol', 1e-12);
        end

        function anUNBOUNDEDConvexInputIsItsOwnEnvelope(testCase)
        % QuaPol.oneNorm is |x|_1: four AFFINE pieces on four unbounded cones, and convex. co f = f,
        % but re-expressing it in H-form was refused, because an unbounded piece's edge SIDE cannot
        % be read off its vertices -- a recession direction settles it, and for a half-plane piece
        % only F does. 11 of the 13 fixtures the LEGACY envelope answered and this one did not were
        % refused for that reason.
        %
        % fIsConvex is set by the caller here to exercise the SHORT-CIRCUIT specifically: it is
        % checked before isConvexQ and returns the same answer without the edge loop. The flag is
        % no longer the only route -- the next test is the same fixture without it -- and the free
        % necessary condition is still checked, see aFalseConvexityAssertionIsRefusedLoudly.
            f = QuaPol.oneNorm();
            f.fIsConvex = true;
            h = biconjQ(f);
            X = [0 0; 1 1; -1 2; 0.5 -0.5; 3 -4];
            testCase.verifyEqual(h.eval(X), f.eval(X), 'AbsTol', 1e-12, ...
                'the envelope of a convex function is the function');
            testCase.verifyTrue(checkQuaConConsistent(h, 400), 'and its cells do not overlap');
        end

        function theConvexBranchIsREACHEDWithNoFlagAtAll(testCase)
        % The capability isConvexQ adds, stated as the property that used to fail: an honest convex
        % multi-piece input is answered without the caller asserting anything. Before it, every
        % fixture here raised PLQ:biconjQ:unbounded unless fIsConvex was set by hand.
        %
        % examples{16} is deliberately NOT in this list: its domain is not convex, so it must stay
        % refused. See DECISIONS.md 2026-09-05.
            P  = QuaPol.examples();
            P2 = QuaPol.examples2();
            L = {'oneNorm', QuaPol.oneNorm(); 'examples{5}', P{5}; 'examples{8}', P{8}; ...
                 'examples{18}', P{18}; 'examples{19}', P{19}; 'examples2{9}', P2{9}};
            X = [0 0; 1 1; -1 2; 0.5 -0.5; 3 -4; 2 2; -3 1];
            for i = 1:size(L,1)
                f = L{i,2};
                testCase.verifyEmpty(f.fIsConvex, ...
                    sprintf('%s must carry no assertion for this test to mean anything', L{i,1}));
                h = biconjQ(f);
                testCase.verifyEqual(h.eval(X), f.eval(X), 'AbsTol', 1e-12, ...
                    sprintf('%s: co f = f, so the envelope must agree with f', L{i,1}));
            end
        end

        function aNonConvexDOMAINIsEnvelopedOverItsHULLNotReturnedUnchanged(testCase)
        % The half-star, examples{16}: f = 0 on five cones meeting at the origin, whose domain is
        % NOT convex. Every piece is affine, so per-piece PSD passes, and a convexity test reading
        % the legacy dom.isConvex would call f convex and hand it back unchanged -- wrong, because
        % co f must be defined on conv(dom f), which is strictly larger.
        %
        % CHANGED 2026-09-05, when the unbounded envelope landed. This used to assert
        % PLQ:biconjQ:unbounded, which recorded only that biconjQ could not do unbounded domains
        % yet -- an incidental fact about coverage, not the property the input is here to pin. Now
        % the case is computed, so the test asserts the PROPERTY: the envelope extends past dom f.
        %
        % The five rays sweep from 0 to 225 degrees, MORE than a half turn, so conv(dom f) is the
        % whole plane and co f = 0 everywhere. (0,-0.5) is the midpoint witness from
        % isConvexQTest -- outside dom f, where f is +infinity -- and it is exactly where an
        % envelope that had been handed back unchanged would still say +infinity.
            P = QuaPol.examples();
            h = biconjQ(P{16});
            testCase.verifyTrue(isinf(P{16}.eval([0 -0.5])), ...
                'f itself is +infinity at the witness, since the domain is not convex there');
            testCase.verifyEqual(h.eval([0 -0.5]), 0, 'AbsTol', 1e-12, ...
                'but co f is finite there, because its domain is conv(dom f)');
            testCase.verifyEqual(h.eval([-1 -1; 1 0; 2 -3; -4 1; 5 5; 0 0]), zeros(6,1), ...
                'AbsTol', 1e-12, 'and co f = 0 on the whole plane');
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
        % upward -- there the vertex values no longer determine the envelope. This piece is
        % INDEFINITE, so its own single-piece envelope is the unimplemented one, and that case is
        % [COAP] A.2-A.5.
            f = biconjQTest.square([0 0 0 0, 1, 0, -1, 0, 0, 0]);   % +1 along (1,0)
            testCase.verifyError(@() biconjQ(f), 'PLQ:biconjQ:notImplemented');
        end

        function CONVEXPiecesWithANonConvexKinkGetTheirOWNRefusal(testCase)
        % The same edge-curvature test trips here, and until 2026-09-05 it gave the same message --
        % which was wrong about this case. Every piece is convex, so NO piece needs [COAP]: each
        % one's own envelope is itself. What is missing is the hull ACROSS pieces.
        %
        % The fixture is the unit square split by its diagonal, both pieces x^2+y^2 up to a
        % multiple of the diagonal's own linear form l = y-x. l vanishes on the diagonal, so f is
        % continuous there for EVERY c, and the sign of c alone decides convexity -- which is what
        % makes the pair below differ in exactly one thing.
            convex = biconjQTest.diagonalKink(2);
            testCase.verifyTrue(isConvexQ(convex), 'c > 0 is the convex side of the knob');
            h = biconjQ(convex);
            X = [0.2 0.1; 0.9 0.4; 0.3 0.8; 0.5 0.5];
            testCase.verifyEqual(h.eval(X), convex.eval(X), 'AbsTol', 1e-12, ...
                'and there co f = f, so it is answered rather than refused');

            kinked = biconjQTest.diagonalKink(-2);
            testCase.verifyFalse(isConvexQ(kinked), 'c < 0 is the non-convex side');
            testCase.verifyError(@() biconjQ(kinked), 'PLQ:biconjQ:convexPieceCoupling', ...
                'convex pieces with a non-convex coupling get their own identifier');
        end

        function aMULTIPIECEEnvelopeIsOneHullOverAllTheLiftedVertices(testCase)
        % WAS a refusal. The envelope COUPLES its pieces -- the convex hull of a union is not the
        % union of hulls -- so unlike the conjugate it is not a fold. But that coupling IS just one
        % hull: each face's graph has its lower hull spanned by that face's lifted vertices, and
        % the hull of a union is the hull of the union of the spanning sets, so the whole envelope
        % is the lower hull of ALL the lifted vertices at once. One step, not a fold.
        %
        % f = -x^2 below the diagonal and -y^2 above it, on the unit square: continuous, and each
        % piece edge-concave. The oracle is the envelope's DEFINITION by linear programming --
        % co f(x) as the least value of sum w_i f(x_i) over convex combinations reaching x -- which
        % shares no machinery with the lower hull.
            V = [0 0; 1 0; 1 1; 0 1];
            E = [1 2 1; 2 3 1; 3 4 1; 4 1 1; 1 3 1];
            F = [1 0; 1 0; 2 0; 2 0; 2 1];        % see conjQTest for why F(5,:) is [2 1]
            f = QuaPol(V, E, [0 0 0 0 -2 0 0 0 0 0; 0 0 0 0 0 0 -2 0 0 0], F);
            testCase.verifyEqual(f.eval([0.8 0.2]), -0.64, 'AbsTol', 1e-12, ...
                'the fixture must be well formed: -x^2 applies BELOW the diagonal');
            testCase.verifyEqual(f.eval([0.2 0.8]), -0.64, 'AbsTol', 1e-12);

            h = biconjQ(f);

            [X, Y] = meshgrid(linspace(0, 1, 121));
            P = [X(:), Y(:)];
            q = f.eval(P);  keep = isfinite(q);  P = P(keep,:);  q = q(keep);
            opts = optimoptions('linprog', 'Display', 'none');
            probes = [0.5 0.5; 0.25 0.75; 0.9 0.1; 0 0; 1 1; 0.5 0.2; 1 0; 0 1];
            for i = 1:size(probes,1)
                [~, fv, flag] = linprog(q, [], [], [P.'; ones(1,size(P,1))], ...
                                        [probes(i,:).'; 1], zeros(size(P,1),1), [], opts);
                testCase.assumeEqual(flag, 1);
                testCase.verifyEqual(h.eval(probes(i,:)), fv, 'RelTol', 1e-6, 'AbsTol', 1e-6, ...
                    sprintf('at (%g,%g)', probes(i,1), probes(i,2)));
            end
        end


        function anInexactInputIsRefusedBeforeAnythingElse(testCase)
            testCase.verifyError(@() biconjQ(QuaPol([sqrt(2) 0 1 0 0 0])), 'PLQ:QuaPol:notExact');
        end

        function underBICONJ_ASSERTAFalseFlagIsREFUTEDRatherThanTrusted(testCase)
        % examples{16} is the case that shows what trusting the flag costs: every piece is affine,
        % so the free necessary condition passes, but the DOMAIN is not convex, and with the flag
        % set biconjQ returned that f as its own convex envelope. isConvexQ can refute it, and
        % under the opt-in global it does.
        %
        % Both halves are asserted, because the point is the DIFFERENCE: off, the flag is still
        % trusted and nothing is raised (that is what makes the flag a short-circuit worth having);
        % on, the same call is refused.
            global BICONJ_ASSERT %#ok<GVMIS>
            prev = BICONJ_ASSERT;
            testCase.addTeardown(@() biconjQTest.setBiconjAssert(prev));

            f = QuaPol.examples();  f = f{16};
            f.fIsConvex = true;

            BICONJ_ASSERT = [];
            err = '';  h = [];
            try, h = biconjQ(f); catch e, err = e.identifier; end
            testCase.verifyEmpty(err, ...
                'with the assert off the flag is trusted, so nothing is raised');
            testCase.verifyNotEmpty(h, 'and an answer comes back -- the wrong one, which is why');

            BICONJ_ASSERT = 1;
            testCase.verifyError(@() biconjQ(f), 'PLQ:biconjQ:notConvex', ...
                'with the assert on the exact test refutes the flag');
        end

        function aFalseConvexityAssertionIsRefusedLoudly(testCase)
        % ALGORITHM.md's rule: the free necessary condition is still checked, because the failure
        % it prevents is silent -- biconj returning a non-convex f as its own convex envelope.
            f = biconjQTest.square([0 0 0 0, 1, 0, -1, 0, 0, 0]);
            f.fIsConvex = true;
            testCase.verifyError(@() biconjQ(f), 'PLQ:biconjQ:notConvex');
        end

        function theExactEngineIsReachableFromTheClassAPI(testCase)
        % biconjQ was callable only as a free function: `QuaPol.biconj` dispatched to biconjCPLQ
        % for 'cplq' and to conj(conj(.)) otherwise, with no name for the exact route.
        % `biconj(f,'exact')` is that route. Same three guarantees as conj's: it reaches biconjQ,
        % it leaves the other engine names alone, and the refusals stay named.
            f = biconjQTest.square([0 0 0 0, 0, 1, 0, 0, 0, 0]);   % q = x*y on the unit square

            hExact = f.biconj('exact');
            testCase.verifyEqual(hExact.kind(), 'QuaCon', '''exact'' must reach biconjQ');
            hDirect = biconjQ(f);
            testCase.verifyEqual(hExact.nf, hDirect.nf);
            X = [0 0; 1 1; 0.5 0.5; 0.25 0.75; 1 0];
            testCase.verifyEqual(hExact.eval(X), hDirect.eval(X), 'AbsTol', 0);
            testCase.verifyEqual(f.biconj('EXACT').kind(), 'QuaCon');

            % 'symbolic' is the legacy engine; 'cplq' is unchanged and still the default
            testCase.verifyEqual(class(f.biconj('symbolic')), class(f.biconj('cplq')));
            testCase.verifyEqual(class(f.biconj()),           class(f.biconj('cplq')));
        end

        function theExactEnginesREFUSALSSurviveTheClassAPI(testCase)
        % A refusal routed through `biconj` must arrive with its own identifier -- the whole point
        % of an engine that never falls back. `x^2 - y^2` on the square is the edge-CONVEX piece,
        % which is the one genuine gap left in biconjQ.
            f = biconjQTest.square([0 0 0 0, 2, 0, -2, 0, 0, 0]);
            testCase.verifyError(@() f.biconj('exact'), 'PLQ:biconjQ:notImplemented');
        end
    end

    methods (Static)
        function f = chainL(c) %#ok<INUSD>
        % Two non-collinear segments meeting at the origin, with NO faces (F all zero, which is
        % what makes the constructor call it a chain). Edge 1 carries -x^2 and edge 2 carries 0.
        % A chain needs one function row PER EDGE -- QuaPol.eval indexes f(i,:) by edge -- and that
        % convention is the whole content of the case.
            f = QuaPol([0 0; 2 0; 0 2], [1 2 1; 1 3 1], ...
                       [0 0 0 0, -2, 0, 0, 0, 0, 0; 0 0 0 0, 0, 0, 0, 0, 0, 0], [0 0; 0 0]);
        end

        function f = diagonalKink(c)
        % The unit square split by the diagonal (0,0)-(1,1). Face 1 (below it) carries x^2+y^2 and
        % face 2 (above) carries x^2+y^2 + c*(y-x). Edge 5 is the diagonal; its direction is (1,1)
        % so its LEFT is the upper triangle, which is why F(5,:) = [2 1].
            q0 = [0 0 0 0, 2, 0, 2, 0, 0, 0];
            q1 = q0;  q1(8) = q1(8) - c;  q1(9) = q1(9) + c;      % + c*(y - x)
            f = QuaPol([0 0; 1 0; 1 1; 0 1], ...
                       [1 2 1; 2 3 1; 3 4 1; 4 1 1; 1 3 1], [q0; q1], ...
                       [1 0; 1 0; 2 0; 2 0; 2 1]);
        end

        function setBiconjAssert(v)
        % Restore the global to whatever it was, including empty. A plain assignment inside the
        % teardown closure would not reach the global without this declaration.
            global BICONJ_ASSERT %#ok<GVMIS>
            BICONJ_ASSERT = v;
        end

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
