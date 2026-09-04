classdef QuaConTest < matlab.unittest.TestCase
% Unit tests for QuaCon, the exact quadratic-on-conic storage: rational faces, primitive integer
% edge conics, and vertices stored as NAMES rather than coordinates.
%
% BUCKET: fast (closed-form integer arithmetic and conicMeet; no symbolic call anywhere).
%
% The oracle for the polyhedral cases is the LEGACY DOUBLE PIPELINE -- the same function built as a
% QuaPol -- which is what CLAUDE.md section 6 asks for: the thing being replaced is the thing that
% judges the replacement. The curved cases have no legacy oracle, because a QuaPar cannot hold them,
% so they are pinned against the definitions instead (the conic's own equation, and the exact
% discriminant that classifies it).

    methods (Test)

        % ---- against the legacy pipeline, which is the oracle ---------------------------------

        function aPolyhedralQuaConAgreesWithTheSameFunctionAsAQuaPol(testCase)
        % THE differential test. The unit square carrying (x^2+y^2)/2, built twice: once the legacy
        % way with double coordinates, once exactly with named vertices and H-form faces. The two
        % must agree at the vertices, inside, and outside.
            [qc, qp] = QuaConTest.unitSquareBothWays();

            testCase.verifyEqual(qc.vertexCoords(), qp.V, 'AbsTol', 0, ...
                'a vertex realised from its NAME must be the same point the legacy mesh stores');

            rng(20260903);
            X = [ 0.25 0.25;  0.5 0.5;  0.9 0.1;  0 0;  1 1;  0.5 0; ...
                  rand(40,2) ];
            [vc, ic] = qc.eval(X);
            vp = qp.eval(X);
            testCase.verifyEqual(vc, vp, 'AbsTol', 1e-12, ...
                'the exact face value must agree with the legacy one inside the domain');
            testCase.verifyTrue(all(ic > 0), 'every point of the square lands in a face');

            outside = [ -1 0.5;  2 0.5;  0.5 -1;  0.5 2;  -3 -3 ];
            [vo, io] = qc.eval(outside);
            testCase.verifyTrue(all(isinf(vo)), 'off the domain the value is +Inf');
            testCase.verifyEqual(io, zeros(size(outside,1),1), 'and no face is reported');
        end

        function everyEdgeOfAPolyhedralQuaConIsExactlyALine(testCase)
        % The classification is DECIDED by the exact discriminant, not estimated.
            qc = QuaConTest.unitSquareBothWays();
            for j = 1:qc.ne
                testCase.verifyEqual(qc.edgeKind(j), 'line', ...
                    sprintf('edge %d of a square is a line', j));
            end
        end

        % ---- the property the whole representation exists to obtain ---------------------------

        function twoSpellingsOfOneCurveBecomeTheSameIntegerRowAndCompareBitwise(testCase)
        % DECISIONS.md 2026-08-17 measured what the alternative costs: two doubles of one exact
        % number, one ULP apart, made a shared facet invisible to merge and Step 3's cell count
        % then grew without bound. Here the two spellings are not merely close, they are ONE ROW.
            base   = [0 0 0 0 1 0];        % y = 0
            scaled = [0 0 0 0 1000 0];     % the same line, written 1000x larger
            negated = -[0 0 0 0 7 0];      % and again, negative and at another scale

            qc = QuaConTest.squareWithEdgeConics([base; [0 0 0 1 0 -1]; scaled - [0 0 0 0 0 1000]; ...
                                                 [0 0 0 1 0 0]]);
            testCase.verifyEqual(ratQ.conic(scaled), ratQ.conic(base), ...
                'a conic has no scale: the two spellings canonicalise to one row');
            testCase.verifyEqual(ratQ.conic(negated), ratQ.conic(base), ...
                'and negating it does not change the curve either');
            testCase.verifyTrue(qc.sameEdgeCurve(1, 1));
            testCase.verifyFalse(qc.sameEdgeCurve(1, 2), ...
                'y = 0 and x = 1 are different curves');
        end

        % ---- what a QuaPar cannot hold, which is why this class exists -------------------------

        function anEllipticEdgeIsStoredHereAndRefusedByQuaPar(testCase)
        % The reason for the whole type. A conjugate's edge between two NON-adjacent pieces lies on
        % a genuine ellipse (CONJ_FIELD_PROOF.md 7.4b), QuaPar.assertParabolicEdges correctly
        % refuses to store that as a parabola, and this is where it goes instead.
            circle = [1 0 1 0 0 -1];                       % x^2 + y^2 - 1 = 0
            testCase.verifyError(@() QuaPar.assertParabolicEdges(circle), 'QuaPar:notParabola', ...
                'a QuaPar must still refuse an ellipse -- this test also pins that it does');

            qc = QuaConTest.halfDisk();
            testCase.verifyEqual(qc.edgeKind(1), 'elliptic');
            testCase.verifyEqual(qc.edgeKind(2), 'line');
            testCase.verifyEqual(ratQ.discriminant(qc.EcQ(1,:)), -4, ...
                'b^2 - 4ac of the unit circle is exactly -4, decided not estimated');
        end

        function aVertexOnACurvedEdgeIsRealisedFromItsNameAlone(testCase)
        % The half disk's two corners are where the circle meets y = 0, i.e. (-1,0) and (1,0), in
        % conicMeet's canonical (lexicographic) order. Nothing about those coordinates is stored.
            qc = QuaConTest.halfDisk();
            testCase.verifyEqual(qc.vertexCoords(), [-1 0; 1 0], 'AbsTol', 1e-12);
            testCase.verifyTrue(isempty(qc.V), ...
                'the inherited double vertex array stays EMPTY by design');
        end

        function aParabolicEdgeIsClassifiedAndRealisedExactly(testCase)
        % y = x^2 capped by y = 1: the two corners are (-1,1) and (1,1).
            qc = QuaConTest.parabolicCap();
            testCase.verifyEqual(qc.edgeKind(1), 'parabola');
            testCase.verifyEqual(ratQ.discriminant(qc.EcQ(1,:)), 0, ...
                'a parabola has discriminant exactly zero');
            testCase.verifyEqual(qc.vertexCoords(), [-1 1; 1 1], 'AbsTol', 1e-12);

            % inside the cap, and outside it on both counts
            [vin, iin] = qc.eval([0 0.5; 0 0.9; 0.5 0.5]);
            testCase.verifyTrue(all(iin == 1) && all(isfinite(vin)));
            [vout, iout] = qc.eval([0 -0.5; 0 2; 2 2]);
            testCase.verifyTrue(all(isinf(vout)) && all(iout == 0));
        end

        % ---- the object cannot be made to lie about itself -------------------------------------

        function aQuaConReportsItselfAsMeshedAlthoughItsDoubleFieldsAreEmpty(testCase)
        % RatCon.isMeshed tests ~isempty(obj.f), and a QuaCon's f is empty BY DESIGN, so the
        % inherited answer would be exactly backwards and would send callers down QuaParCPLQ's
        % "mesh not yet reconstructed" path.
            qc = QuaConTest.unitSquareBothWays();
            testCase.verifyTrue(qc.isMeshed());
            testCase.verifyEqual(qc.kind(), 'QuaCon');
            testCase.verifyTrue(isa(qc, 'RatCon') && isa(qc, 'Qua') && isa(qc, 'Conic'));
            testCase.verifyFalse(isa(qc, 'RatPar'), ...
                'a QuaCon is deliberately NOT a RatPar -- its edges need not be parabolas');
            testCase.verifyTrue(isempty(qc.f) && isempty(qc.Ec) && isempty(qc.V), ...
                'the inherited DOUBLE fields stay empty; the exact ones are the truth');
            testCase.verifyEqual(qc.den, repmat([0 0 1], qc.nf, 1), ...
                'a QuaCon is a Qua, so its denominator is identically 1');
        end

        function aQuaConCannotBeBuiltWithADenominatorOtherThanOne(testCase)
        % The Qua axis is enforced from RatCon's single set.den validator, and adding a new type to
        % the lattice must not have created a hole in it.
            qc = QuaConTest.unitSquareBothWays();
            testCase.verifyError(@() setfield(qc, 'den', [1 0 1]), ...
                'RatPar:denMustBeTrivial'); %#ok<SFLD>
        end

        % ---- the constructor refuses what it cannot represent, by name --------------------------

        function anUnrealisableVertexNameIsRefusedRatherThanStored(testCase)
        % Two lines meet in ONE point, so asking for the second is not a rounding question.
            [V, E, Ec, fN, fD, F, FC] = QuaConTest.squareParts();
            V(1,3) = 2;
            testCase.verifyError(@() QuaCon(V, E, Ec, fN, fD, F, FC), 'QuaCon:badRootIndex');
        end

        function aVertexNamingOneEdgeTwiceIsRefused(testCase)
            [V, E, Ec, fN, fD, F, FC] = QuaConTest.squareParts();
            V(1,1:2) = [1 1];
            testCase.verifyError(@() QuaCon(V, E, Ec, fN, fD, F, FC), 'QuaCon:badVertexName');
        end

        function twoEdgesOnOneCurveNameNoVertexAndAreRefused(testCase)
        % Coincident conics meet in infinitely many points; that is a degenerate configuration the
        % caller must handle, not a vertex.
            [V, E, Ec, fN, fD, F, FC] = QuaConTest.squareParts();
            Ec(4,:) = Ec(1,:);                       % edge 4 becomes the same line as edge 1
            V(1,1:2) = [1 4];
            testCase.verifyError(@() QuaCon(V, E, Ec, fN, fD, F, FC), 'QuaCon:degenerateVertex');
        end

        function anAllZeroEdgeConicIsRefusedBecauseItNamesNoCurve(testCase)
        % QuaPar spells a straight edge as an all-zero Ec row and recovers the line from stored
        % endpoint COORDINATES. That is exactly the dependence on coordinates the H-form removes,
        % so here every edge must carry its own equation.
            [V, E, Ec, fN, fD, F, FC] = QuaConTest.squareParts();
            Ec(2,:) = zeros(1,6);
            testCase.verifyError(@() QuaCon(V, E, Ec, fN, fD, F, FC), 'ratQ:zeroConic');
        end

        function aCubicFaceIsRefusedBecauseItHasNoConicLevelSet(testCase)
            [V, E, Ec, fN, fD, F, FC] = QuaConTest.squareParts();
            fN(1,1) = 1;
            testCase.verifyError(@() QuaCon(V, E, Ec, fN, fD, F, FC), 'QuaCon:cubicFace');
        end

        function aFaceConstraintCarriesOneOfTheThreeSidesAndNothingElse(testCase)
        % CHANGED 2026-09-04, and the change is a widening rather than a weakening. This used to
        % assert that a side of 0 was REFUSED, on the reading that 0 meant "omitted". It now means
        % ON the curve -- an EQUALITY -- which is how a THIN face is expressed: dom f* is a line
        % when the primal is a positive semidefinite singular quadratic on the plane, and a point
        % when the primal is affine there. Both are correct conjugates that had nowhere to go.
        %
        % What the test exists for is unchanged: a side must be one of the three meanings and
        % nothing else, because for an ellipse the side cannot be recovered from the endpoints.
            [V, E, Ec, fN, fD, F, FC] = QuaConTest.squareParts();

            FC{1}(1,2) = 0;                       % an equality is now legitimate
            qc = QuaCon(V, E, Ec, fN, fD, F, FC);
            testCase.verifyEqual(qc.FC{1}(1,2), 0);

            for badSide = [2, -3, 0.5]
                FC{1}(1,2) = badSide;
                testCase.verifyError(@() QuaCon(V, E, Ec, fN, fD, F, FC), ...
                    'QuaCon:badConstraint', sprintf('side %g must be refused', badSide));
            end
        end

        function anEqualitySideMakesTheFaceTHINAndEvalRespectsIt(testCase)
        % The behaviour the new side buys, asserted directly: a face carrying an equality is only
        % entered ON the curve. The unit square's face is cut down to the segment y = 0 within it.
            [V, E, Ec, fN, fD, F, FC] = QuaConTest.squareParts();
            FC{1} = [FC{1}; 1 0];                 % edge 1 is y = 0; demand ON it
            qc = QuaCon(V, E, Ec, fN, fD, F, FC);
            testCase.verifyTrue(isfinite(qc.eval([0.5 0])), 'on the line, inside the square');
            testCase.verifyTrue(isinf(qc.eval([0.5 0.5])), 'off the line: not in the face');
            testCase.verifyTrue(isinf(qc.eval([2 0])), 'on the line but outside the square');
        end

        function theFaceCoefficientsAreStoredCanonicallyWhateverSpellingIsPassed(testCase)
        % A face HAS a scale, so canon must preserve the VALUE while fixing the spelling.
            [V, E, Ec, fN, fD, F, FC] = QuaConTest.squareParts();
            qc1 = QuaCon(V, E, Ec, fN, fD, F, FC);
            qc2 = QuaCon(V, E, Ec, 6*fN, 6*fD, F, FC);
            testCase.verifyEqual(qc2.fN, qc1.fN, 'the same function, canonically, is one row');
            testCase.verifyEqual(qc2.fD, qc1.fD);
        end
    end

    methods (Static)

        function [qc, qp] = unitSquareBothWays()
        % (x^2+y^2)/2 on the unit square, built exactly and built the legacy way.
            [V, E, Ec, fN, fD, F, FC] = QuaConTest.squareParts();
            qc = QuaCon(V, E, Ec, fN, fD, F, FC);
            qp = QuaPol([0 0; 1 0; 1 1; 0 1], [1 2 1; 2 3 1; 3 4 1; 4 1 1], ...
                        [0 0 0 0 1 0 1 0 0 0], [1 0; 1 0; 1 0; 1 0]);
        end

        function [V, E, Ec, fN, fD, F, FC] = squareParts()
        % The unit square in H-form. Edges, in order: y=0, x=1, y=1, x=0.
        % Vertices are NAMES: (0,0) is edges 1 and 4, (1,0) is 1 and 2, and so on.
            Ec = [0 0 0 0 1  0;      % y = 0
                  0 0 0 1 0 -1;      % x = 1
                  0 0 0 0 1 -1;      % y = 1
                  0 0 0 1 0  0];     % x = 0
            V  = [1 4 1;             % (0,0)
                  1 2 1;             % (1,0)
                  2 3 1;             % (1,1)
                  3 4 1];            % (0,1)
            E  = [1 2 1; 2 3 1; 3 4 1; 4 1 1];
            fN = [0 0 0 0 1 0 1 0 0 0];   fD = 1;
            F  = [1 0; 1 0; 1 0; 1 0];
            % the face is  y >= 0,  x <= 1,  y <= 1,  x >= 0
            FC = { [1 1; 2 -1; 3 -1; 4 1] };
        end

        function qc = squareWithEdgeConics(Ec)
        % the square again, with the caller's own spelling of the four edge conics
            [V, E, ~, fN, fD, F, FC] = QuaConTest.squareParts();
            qc = QuaCon(V, E, Ec, fN, fD, F, FC);
        end

        function qc = halfDisk()
        % The upper half of the unit disk: x^2 + y^2 <= 1 and y >= 0. Edge 1 is the CIRCLE, which
        % is what a QuaPar cannot hold; edge 2 is the diameter.
            Ec = [1 0 1 0 0 -1;      % x^2 + y^2 - 1 = 0   (elliptic)
                  0 0 0 0 1  0];     % y = 0
            V  = [1 2 1;             % (-1,0), the first of the two crossings lexicographically
                  1 2 2];            % ( 1,0)
            E  = [1 2 1; 2 1 1];
            fN = [0 0 0 0 1 0 1 0 0 0];   fD = 1;
            F  = [1 0; 1 0];
            FC = { [1 -1; 2 1] };    % inside the circle, above the diameter
            qc = QuaCon(V, E, Ec, fN, fD, F, FC);
        end

        function qc = parabolicCap()
        % { (x,y) : y >= x^2, y <= 1 }. Edge 1 is the parabola, edge 2 the cap.
            Ec = [1 0 0 0 -1  0;     % x^2 - y = 0        (parabola)
                  0 0 0 0  1 -1];    % y = 1
            V  = [1 2 1;             % (-1,1)
                  1 2 2];            % ( 1,1)
            E  = [1 2 1; 2 1 1];
            fN = [0 0 0 0 1 0 1 0 0 0];   fD = 1;
            F  = [1 0; 1 0];
            FC = { [1 -1; 2 -1] };   % x^2 - y <= 0  and  y - 1 <= 0
            qc = QuaCon(V, E, Ec, fN, fD, F, FC);
        end
    end
end
