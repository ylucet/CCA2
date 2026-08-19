classdef testcPLQ < matlab.unittest.TestCase
% The cPLQ pipeline on two bilinear fixtures, STAGE BY STAGE, each stage verified numerically.
%
% WHAT THIS USED TO BE, and why it changed (2026-08-19). Eight tests, ZERO assertions: each ran
% some prefix of `triangulate -> convexEnvelope -> conjugate -> maximum -> biconjugateF`, printed
% the result and returned. They passed if nothing threw. That cost 3562 s -- half the slow bucket
% -- to establish that the functions return, and when `testRectBiconj` began failing there was
% nothing to read but an exception after 53 minutes.
%
% Now every stage:
%   * asserts its own output against a DEFINITION (plqCheck: the envelope underestimates f and
%     touches it at the vertices; f* equals the sup over the domain; f** is a convex
%     underestimator) -- not against a golden value, so nothing needs re-pinning;
%   * caches that output (plqStage) so the NEXT stage starts where this one finished. A cold run
%     costs what it always did; a re-run after a failure re-does only the broken stage onward.
%
% THE FIXTURES are unchanged, deliberately -- this is a re-verification of the same inputs, not a
% new test. PRect is x*y over two adjacent polygons, PRect3 is x*y over one quadrilateral.

    properties
        PRect
        PRect3
        Poly
    end

    methods (TestMethodSetup)
       function setUpTestData(testCase)
            % Initialize common data for all test methods
            x=sym('x');
            y=sym('y');
            f=symbolicFunction(x*y);

            d(1)=domain([-5,-4;0,-4;1,3;-5,5],x,y);
            d(2)=domain([0,-4;2,0;2,1;1,3],x,y);
            d(3)=domain([-1,1;-3,-3;-4,-3],x,y);
            d(4)=domain([1,0;3,1;2,2;0,1],x,y);
            d(5)=domain([-5,-4;0,-4;2,0;2,1;1,3;-5,5],x,y);

            p(1) = plq_1p(d(1),f);
            p(2) = plq_1p(d(2),f);
            testCase.PRect = plq(p);

            p(3) = plq_1p(d(5),f);
            testCase.Poly = plq(p(3));

            d(6)=domain([0,0;1,1;2,1;2,0],x,y);
            testCase.PRect3 = plq([plq_1p(d(6),symbolicFunction(x*y))]);
       end
    end

    methods (Static)
        function S = dualPoints()
        % Dual points for the conjugate check: a spread that puts the sup at different faces of
        % the domain, including the axes and both diagonals.
            S = [0 0; 1 0; 0 1; 1 1; -1 1; 1 -1; -1 -1; 2 0.5; -0.5 2; 3 -2];
        end
    end

    methods (Test)

        % ================= PRect: x*y over two adjacent polygons =========================

        function rectTriangulates(testCase)
            p = plqStage.get('PRect', 'tri', @() testCase.PRect.triangulate);
            testCase.verifyGreaterThanOrEqual(p.nPieces, testCase.PRect.nPieces, ...
                'triangulation must not lose pieces');
            for i = 1:p.nPieces
                testCase.verifyEqual(p.pieces(i).d.polygon.nv, 3, ...
                    sprintf('piece %d is not a triangle after triangulate', i));
            end
        end

        function rectConvexEnvelopeUnderestimates(testCase)
            p = plqStage.get('PRect', 'ce', @() ...
                plqStage.get('PRect', 'tri', @() testCase.PRect.triangulate).convexEnvelope);
            for i = 1:p.nPieces
                plqCheck.envelopeUnderestimates(testCase, p.pieces(i), ...
                    sprintf('PRect piece %d', i));
            end
        end

        function rectConjugateMatchesTheSup(testCase)
        % Per PIECE: the conjugate of q + I_T is the sup over that triangle alone.
        %
        % WHICH PROPERTY HOLDS THE ANSWER, and it is not the obvious one. `plq.conjugate` leaves
        % each piece's result in `.conjugates` -- one cell per ENVELOPE FACE, which still have to
        % be maxed against each other. `.maxConjugate` is what `maximumConjugate` writes, and it
        % is the object that equals the sup. Comparing `.conjugates` to the sup reports every
        % point "uncovered", which is what the first version of this test did.
            p = plqStage.get('PRect', 'conj', @() ...
                plqStage.get('PRect', 'tri', @() testCase.PRect.triangulate).conjugate);
            for i = 1:p.nPieces
                q = p.pieces(i).maximumConjugate;
                plqCheck.conjugateMatchesSup(testCase, q.maxConjugate, ...
                    q.f.f, q.d.polygon.vars, q.d, ...
                    testcPLQ.dualPoints(), sprintf('PRect piece %d f*', i));
            end
        end

        function rectMaximumIsTheConjugateOfTheWholeDomain(testCase)
        % Step 3. The max ACROSS pieces is f* of the union, so it must equal the sup over the
        % ORIGINAL (pre-triangulation) domain -- which is the property that makes Step 3 correct
        % and the one a crash test could never see.
            p = plqStage.get('PRect', 'max', @() ...
                plqStage.get('PRect', 'tri', @() testCase.PRect.triangulate).maximum);
            testCase.verifyNotEmpty(p.maxConjugate, 'Step 3 produced no cells');
            S = testcPLQ.dualPoints();
            for i = 1:size(S,1)
                s = S(i,:);
                got = evalFunctionNDomain(p.maxConjugate, s);
                want = -inf;
                for k = 1:testCase.PRect.nPieces
                    q = testCase.PRect.pieces(k);
                    want = max(want, plqCheck.supOverDomain(q.f.f, q.d.polygon.vars, q.d, s));
                end
                testCase.verifyFalse(isnan(got), sprintf( ...
                    'PRect f* uncovered at (%g,%g)', s(1), s(2)));
                if isnan(got), continue, end
                tol = 1e-5 * max(1, abs(want));
                testCase.verifyEqual(got, want, 'AbsTol', tol, sprintf( ...
                    'PRect f*(%g,%g): got %.9g, sup over the domain is %.9g', ...
                    s(1), s(2), got, want));
            end
        end

        function rectBiconjugateIsAConvexUnderestimator(testCase)
        % The stage that used to be the 3198 s crash test. It now starts from the cached Step 3
        % result, so a failure here is attributable to biconjugateF alone.
            p = plqStage.get('PRect', 'biconj', @() ...
                plqStage.get('PRect', 'max', @() ...
                    plqStage.get('PRect', 'tri', @() testCase.PRect.triangulate).maximum ...
                ).biconjugateF);
            testCase.verifyNotEmpty(p.biconjugate, 'biconjugateF produced no cells');
            for k = 1:testCase.PRect.nPieces
                q = testCase.PRect.pieces(k);
                plqCheck.biconjugateIsAConvexUnderestimator(testCase, p.biconjugate, ...
                    q.f.f, q.d.polygon.vars, q.d, sprintf('PRect f** over piece %d', k));
            end
        end

        % ================= PRect3: x*y over one quadrilateral ============================

        function rect3TriangulatesAndEnvelopeUnderestimates(testCase)
            p = plqStage.get('PRect3', 'ce', @() ...
                plqStage.get('PRect3', 'tri', @() testCase.PRect3.triangulate).convexEnvelope);
            for i = 1:p.nPieces
                plqCheck.envelopeUnderestimates(testCase, p.pieces(i), ...
                    sprintf('PRect3 piece %d', i));
            end
        end

        function rect3ConjugateMatchesTheSup(testCase)
            p = plqStage.get('PRect3', 'conj', @() ...
                plqStage.get('PRect3', 'tri', @() testCase.PRect3.triangulate).conjugate);
            for i = 1:p.nPieces
                % .maxConjugate, not .conjugates -- see rectConjugateMatchesTheSup.
                q = p.pieces(i).maximumConjugate;
                plqCheck.conjugateMatchesSup(testCase, q.maxConjugate, ...
                    q.f.f, q.d.polygon.vars, q.d, ...
                    testcPLQ.dualPoints(), sprintf('PRect3 piece %d f*', i));
            end
        end

        function rect3MaximumIsTheConjugateOfTheWholeDomain(testCase)
            p = plqStage.get('PRect3', 'max', @() ...
                plqStage.get('PRect3', 'tri', @() testCase.PRect3.triangulate).maximum);
            testCase.verifyNotEmpty(p.maxConjugate, 'Step 3 produced no cells');
            q = testCase.PRect3.pieces(1);
            S = testcPLQ.dualPoints();
            for i = 1:size(S,1)
                s = S(i,:);
                got = evalFunctionNDomain(p.maxConjugate, s);
                want = plqCheck.supOverDomain(q.f.f, q.d.polygon.vars, q.d, s);
                testCase.verifyFalse(isnan(got), sprintf( ...
                    'PRect3 f* uncovered at (%g,%g)', s(1), s(2)));
                if isnan(got), continue, end
                testCase.verifyEqual(got, want, 'AbsTol', 1e-5 * max(1, abs(want)), sprintf( ...
                    'PRect3 f*(%g,%g): got %.9g, sup over the domain is %.9g', ...
                    s(1), s(2), got, want));
            end
        end

        function rect3BiconjugateIsAConvexUnderestimator(testCase)
            p = plqStage.get('PRect3', 'biconj', @() ...
                plqStage.get('PRect3', 'max', @() ...
                    plqStage.get('PRect3', 'tri', @() testCase.PRect3.triangulate).maximum ...
                ).biconjugateF);
            testCase.verifyNotEmpty(p.biconjugate, 'biconjugateF produced no cells');
            q = testCase.PRect3.pieces(1);
            plqCheck.biconjugateIsAConvexUnderestimator(testCase, p.biconjugate, ...
                q.f.f, q.d.polygon.vars, q.d, 'PRect3 f**');
        end
    end
end
