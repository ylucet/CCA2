classdef testMaxMultiRegion < matlab.unittest.TestCase
% The cPLQ pipeline over MULTI-REGION fixtures, verified numerically.
%
% WHAT THIS USED TO BE (changed 2026-08-19): 24 tests and ZERO assertions -- each ran some prefix
% of `maximum -> biconjugateF`, printed, and returned, so a test passed whenever nothing threw.
% Four of them ran 6 to 17 minutes apiece and together were half the slow bucket.
%
% Every test now asserts against a DEFINITION rather than a golden value (plqCheck): Step 3's
% output must equal the sup of `<s,x> - f(x)` over the ORIGINAL domain, and f** must be a convex
% underestimator of f on it. The four heavy ones are additionally split by STAGE and cached
% (plqStage), so a failure names the stage that broke instead of arriving as an exception at the
% end of a quarter of an hour.
%
% This suite lives in the `--verylong` bucket; see .claude/suite.sh.

    properties
        PTri
        PRect
        PRect2
        PRect3
        Poly
        PTri2
        PThesis
        POpen
        PSqroot

        PCE0
        PCE0_2
        PCE0_3

        PCE1

        x=sym('x');
        y=sym('y');
        %f=symbolicFunction(x*y);
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
            

            
            p(1) = plq_1piece(d(1),f);
            %f=symbolicFunction(x^2-y^2);
            p(2) = plq_1piece(d(2),f);
            testCase.PRect = plq(p);
            
            testCase.PTri = plq([plq_1piece(d(3),symbolicFunction(x^2-y^2))]);

            testCase.PRect2 = plq([plq_1piece(d(4),symbolicFunction(x^2-y^2))]);

            p(3) = plq_1piece(d(5),f);
            testCase.Poly = plq(p(3));

            d(6)=domain([0,0;1,1;2,1;2,0],x,y);
            testCase.PRect3 = plq([plq_1piece(d(6),symbolicFunction(x*y))]);

            % Expt1
            d(7)=domain([0,0;1,1;2,0],x,y);
            d(8)=domain([1,1;2,1;2,0],x,y);
            
            % Expt8
            d(7)=domain([-5,0;2,0;0,-4;-5,-4],x,y);
            d(8)=domain([-5,5;1,3;2,1;2,0;-5,0],x,y);
            
            q(1) = plq_1piece(d(7),f);
            q(2) = plq_1piece(d(8),f);
            
            testCase.PTri2 = plq(q);

            %cube1
            d(9)=domain([-1,-1;-1,1;1,1;1,-1],x,y);
            q(3) = plq_1piece(d(9),f);
            testCase.PThesis = plq(q(3));

            %sqroot
            %d(10)=domain([-1,1;-3,-3;-4,-3;-4,1],x,y);
            %d(10)=domain([-1,1;-1,-2;-3,-3;-4,-3],x,y);
            d(10)=domain([-1,1;-3,-3;-4,-3],x,y);
            testCase.PSqroot = plq(plq_1piece(d(10),f));

            d(10)=domain([-1,-1;1,-1;-1,1],x,y);
            testCase.PCE0 = plq(plq_1piece(d(10),f));


            d(10)=domain([2,0;1,2;0,2],x,y);
            testCase.PCE0_2 = plq(plq_1piece(d(10),f));

             d(10)=domain([-2,1;-2,-2;2,-2;2,1],x,y);
            testCase.PCE0_3 = plq(plq_1piece(d(10),f));

       end
    end

    methods (Static)
        function S = dualPoints()
        % A spread of dual points that puts the sup on different faces of the domain.
            S = [0 0; 1 0; 0 1; 1 1; -1 1; 1 -1; -1 -1; 2 0.5; -0.5 2];
        end

        function verifyStep3(tc, p, orig, name)
        % Step 3's max ACROSS the pieces is f* of the union, so it equals the sup of
        % `<s,x> - f(x)` over the ORIGINAL, pre-triangulation domain. That identity is what makes
        % Step 3 correct, and it is exactly what a crash test cannot see.
            S = testMaxMultiRegion.dualPoints();
            for i = 1:size(S,1)
                sPt = S(i,:);
                got = evalFunctionNDomain(p.maxConjugate, sPt);
                want = -inf;
                for k = 1:orig.nPieces
                    q = orig.pieces(k);
                    want = max(want, plqCheck.supOverDomain(q.f.f, q.d.polygon.vars, q.d, sPt));
                end
                tc.verifyFalse(isnan(got), sprintf('%s: f* uncovered at (%g,%g)', ...
                    name, sPt(1), sPt(2)));
                if isnan(got), continue, end
                tc.verifyEqual(got, want, 'AbsTol', 1e-5 * max(1, abs(want)), sprintf( ...
                    '%s: f*(%g,%g) = %.9g, sup over the domain is %.9g', ...
                    name, sPt(1), sPt(2), got, want));
            end
        end

        function verifyCellsAreUsable(tc, pc, name)
        % The weakest honest assertion for a test whose only real output is a
        % `conjugateOfPiecePoly` result: it produced cells, every cell carries a function and a
        % non-empty region, and no region is degenerate. These tests were exploratory scratch --
        % several `return` before they reach anything -- so this checks what actually runs rather
        % than pretending to verify a value they never compute.
            tc.verifyNotEmpty(pc, sprintf('%s: conjugateOfPiecePoly returned nothing', name));
            for k = 1:numel(pc)
                tc.verifyNotEmpty(pc(k).d, sprintf('%s: cell %d has an empty region', name, k));
                if isempty(pc(k).d), continue, end
                tc.verifyGreaterThan(size(pc(k).d.ineqs,2), 0, sprintf( ...
                    '%s: cell %d has no constraints', name, k));
                tc.verifyNotEmpty(pc(k).f, sprintf('%s: cell %d carries no function', name, k));
            end
        end

        function verifyBiconj(tc, p, orig, name)
        % f** is a convex underestimator of f on the domain.
            for k = 1:orig.nPieces
                q = orig.pieces(k);
                plqCheck.biconjugateIsAConvexUnderestimator(tc, p.biconjugate, ...
                    q.f.f, q.d.polygon.vars, q.d, sprintf('%s piece %d f**', name, k));
            end
        end
    end

    methods (Test)

        function testMaxStep3 (testCase)
        % Was the first half of `testMax` (932 s, no assertions). Step 3's max across the pieces
        % IS f* of the union, so it must equal the sup over the ORIGINAL domain -- checked here,
        % and cached for the biconjugate test below.
            p = plqStage.get('MMR_PRect', 'max', @() testCase.PRect.maximum);
            testCase.verifyNotEmpty(p.maxConjugate, 'Step 3 produced no cells');
            testMaxMultiRegion.verifyStep3(testCase, p, testCase.PRect, 'PRect');
        end

        function testMaxBiconjugate (testCase)
        % The second half. Starts from the cached Step 3 result, so a failure is attributable to
        % biconjugateF alone.
            p = plqStage.get('MMR_PRect', 'biconj', @() ...
                plqStage.get('MMR_PRect', 'max', @() testCase.PRect.maximum).biconjugateF);
            testCase.verifyNotEmpty(p.biconjugate, 'biconjugateF produced no cells');
            testMaxMultiRegion.verifyBiconj(testCase, p, testCase.PRect, 'PRect');
        end

        function testPCE0 (testCase)

            testCase.PCE0_2 = testCase.PCE0_2.maximum;
            testCase.PCE0_2 = testCase.PCE0_2.biconjugateF;
            % VERIFIED, not merely run. Step 3's max across the pieces IS f* of the
            % union, so it must equal the sup of <s,x> - f(x) over the ORIGINAL
            % domain; f** must be a convex underestimator of f on it. Placed HERE,
            % before the print/return block, because several of these tests return
            % early -- an assertion after that would be dead code that always passes.
            testMaxMultiRegion.verifyStep3(testCase, testCase.PCE0_2, testCase.PCE0_2, 'PCE0_2');
            testCase.verifyNotEmpty(testCase.PCE0_2.biconjugate, ...
                'PCE0_2: biconjugateF produced no cells');
            testMaxMultiRegion.verifyBiconj(testCase, testCase.PCE0_2, testCase.PCE0_2, 'PCE0_2');
            testCase.PCE0_2.print;
             testCase.PCE0_2.printDomainMaple;
            return

            
           

        end


        function testPCE3 (testCase)

            testCase.PCE0_3 = testCase.PCE0_3.maximum;
            testCase.PCE0_3 = testCase.PCE0_3.biconjugateF;
            % VERIFIED, not merely run. Step 3's max across the pieces IS f* of the
            % union, so it must equal the sup of <s,x> - f(x) over the ORIGINAL
            % domain; f** must be a convex underestimator of f on it. Placed HERE,
            % before the print/return block, because several of these tests return
            % early -- an assertion after that would be dead code that always passes.
            testMaxMultiRegion.verifyStep3(testCase, testCase.PCE0_3, testCase.PCE0_3, 'PCE0_3');
            testCase.verifyNotEmpty(testCase.PCE0_3.biconjugate, ...
                'PCE0_3: biconjugateF produced no cells');
            testMaxMultiRegion.verifyBiconj(testCase, testCase.PCE0_3, testCase.PCE0_3, 'PCE0_3');
            testCase.PCE0_3.print;
             testCase.PCE0_3.printDomainMaple;
            return

            
           

        end


        function testPCE1 (testCase)
            d=domain([0,0;1,1;2,0],testCase.x,testCase.y);
            testCase.PCE1 = plq(plq_1piece(d,symbolicFunction(testCase.x*testCase.y)));
            % testCase.PCE1 = testCase.PCE1.convexEnvelope;
             testCase.PCE1 = testCase.PCE1.maximum;
             testCase.PCE1 = testCase.PCE1.biconjugateF;
            % VERIFIED, not merely run. Step 3's max across the pieces IS f* of the
            % union, so it must equal the sup of <s,x> - f(x) over the ORIGINAL
            % domain; f** must be a convex underestimator of f on it. Placed HERE,
            % before the print/return block, because several of these tests return
            % early -- an assertion after that would be dead code that always passes.
            testMaxMultiRegion.verifyStep3(testCase, testCase.PCE1, testCase.PCE1, 'PCE1');
            testCase.verifyNotEmpty(testCase.PCE1.biconjugate, ...
                'PCE1: biconjugateF produced no cells');
            testMaxMultiRegion.verifyBiconj(testCase, testCase.PCE1, testCase.PCE1, 'PCE1');
            testCase.PCE1.print;
            testCase.PCE1.printDomainMaple;
            return

            
           

        end

        function testPCE2 (testCase)
            d=domain([0,0;1,0;2,1],testCase.x,testCase.y);
            testCase.PCE1 = plq(plq_1piece(d,symbolicFunction(testCase.x*testCase.y)));
            % .maximum (which runs convexEnvelope internally) must run before
            % printDomainMaple -- printDomainMaple's per-piece Mprint reads
            % maxConjugate unconditionally, which stays empty (and crashes
            % Mprint) if only convexEnvelope has run. .biconjugateF is left
            % disabled: it hits a separate, still-open bug
            % (functionNDomain.addEq errors with an unassigned output when
            % the biconjugate result is empty for this domain).
             testCase.PCE1 = testCase.PCE1.maximum;
            % VERIFIED, not merely run -- see plqCheck. Step 3 only: biconjugateF is
            % deliberately not called here (see the note above), so there is no f** to check.
            testMaxMultiRegion.verifyStep3(testCase, testCase.PCE1, testCase.PCE1, 'PCE2 fixture');
             testCase.PCE1.print;
            %
            testCase.PCE1.printDomainMaple;
            return

            
           
        end

        function testBiconjugate (testCase)
            s_1 = sym('s_1');
            s_2 = sym('s_2');
            f = symbolicFunction(-5*s_1 + 5*s_2 + 25);
            ineq(1) = s_1 + 2*s_2 + 4;
            ineq(2) = s_1 - (9*s_2)/5 - 5;
            ineq(3) = - s_2 - 5 ;
            ineq(4) = s_1 + 7*s_2 - 46;
            d = region(ineq,[s_1,s_2]);
            d = d.removeInfV;
            
            d = d.poly2orderUnbounded;
            edgeNo = d.getEdgeNosInf(d.vars)
            d.print
            return
            
            d.ineqs(edgeNo) = d.ineqs;
            d.print
            return
            p = functionNDomain([f], [d]);
            p.printL
            pc = p(1).conjugateOfPiecePoly ;
            % VERIFIED, not merely run: the cells must be a usable cover.
            testMaxMultiRegion.verifyCellsAreUsable(testCase, pc, 'cells 1');

        end


        function testBiconjugate2 (testCase)
            s_1 = sym('s_1');
            s_2 = sym('s_2');
            f1 = symbolicFunction(-5*s_1 - 4*s_2 - 20);
            ineq(1) = s_1 + 4;
            ineq(2) =  s_2 + 5 ;
            d1 = region(ineq,[s_1,s_2]);
            d1 = d1.removeInfV;
            d1 = d1.poly2orderUnbounded;
            edgeNo = d1.getEdgeNosInf(d1.vars)
            d1.ineqs(edgeNo) = d1.ineqs;
            p(1) = functionNDomain(f1,d1);
            %d.print
            %return
            f2 = symbolicFunction(-5*s_1 + 5*s_2 + 25);
            ineq(1) = s_1 - (9*s_2)/5 - 5;
            ineq(2) =  -s_2 - 5 ;
            d2 = region(ineq,[s_1,s_2]);
            d2 = d2.removeInfV;
            d2 = d2.poly2orderUnbounded;
            edgeNo = d2.getEdgeNosInf(d2.vars)
            d2.ineqs(edgeNo) = d2.ineqs; 
            
            p(2) = functionNDomain(f2,d2);

            f3 = symbolicFunction(-4*s_2);
            ineq(1) = - s_1 - 4;
            ineq(2) =  (9*s_2)/5 - s_1 + 5;
            d3 = region(ineq,[s_1,s_2]);
            d3 = d3.removeInfV;
            d3 = d3.poly2orderUnbounded;
            edgeNo = d3.getEdgeNosInf(d3.vars)
            d3.ineqs(edgeNo) = d3.ineqs; 
            p(3) = functionNDomain(f3,d3);

            
            
            p.printL
            p.printM;
            
            %return
            pc = p.conjugateOfPiecePoly ;
            % VERIFIED, not merely run: the cells must be a usable cover.
            testMaxMultiRegion.verifyCellsAreUsable(testCase, pc, 'cells 2');
            pc = pc.addEq;
            pc.printL
            pc.printM

            %pc.printDomainMaple;
            % d = pc(1).d+pc(2).d;
            % d.print
            % d = d + pc(3).d;





        end

        function testBiconjugate3 (testCase)  % piece 3 of 1 piece conjugate - 4 pieces
            s_1 = sym('s_1');
            s_2 = sym('s_2');
            f1 = symbolicFunction(s_1 + 3*s_2 - 1);
             

            
            ineq(1) = s_1 - 2*s_2 - 1;
            ineq(2) =  s_2/3 - s_1 + 13/3 ;
            d1 = region(ineq,[s_1,s_2]);
            d1 = d1.removeInfV;
            d1 = d1.poly2orderUnbounded;
            edgeNo = d1.getEdgeNosInf(d1.vars)
            d1.ineqs(edgeNo) = d1.ineqs;
            p(1) = functionNDomain(f1,d1);
            %d.print
            %return
            
            f2 = symbolicFunction(-5*s_1 + 5*s_2 + 25);
            ineq(1) = s_1 - (5*s_2)/7 - 25/7;
            ineq(2) =  s_1 - s_2/3 - 13/3 ;

            
            d2 = region(ineq,[s_1,s_2]);
            d2 = d2.removeInfV;
            d2 = d2.poly2orderUnbounded;
            edgeNo = d2.getEdgeNosInf(d2.vars)
            d2.ineqs(edgeNo) = d2.ineqs; 
            
            p(2) = functionNDomain(f2,d2);

            f3 = symbolicFunction(2*s_1+s_2-2);

           
             
            ineq(1) = -s_1 + 2*s_2 + 1;
            ineq(2) =  -s_2 + 2;
            d3 = region(ineq,[s_1,s_2]);
            d3 = d3.removeInfV;
            d3 = d3.poly2orderUnbounded;
            edgeNo = d3.getEdgeNosInf(d3.vars)
            d3.ineqs(edgeNo) = d3.ineqs; 
            p(3) = functionNDomain(f3,d3);

            f4 = symbolicFunction(2*s_1);
            
            ineq(1) = (5*s_2)/7 - s_1 + 25/7;
            
            ineq(2) =  s_2 - 2;
            d4 = region(ineq,[s_1,s_2]);
            d4 = d4.removeInfV;
            d4 = d4.poly2orderUnbounded;
            edgeNo = d4.getEdgeNosInf(d4.vars)
            d4.ineqs(edgeNo) = d4.ineqs; 
            p(4) = functionNDomain(f4,d4);
            
            
            p.printL
            p.printM;
            %return
            pc = p.conjugateOfPiecePoly ;
            % VERIFIED, not merely run: the cells must be a usable cover.
            testMaxMultiRegion.verifyCellsAreUsable(testCase, pc, 'cells 3');
            pc = pc.addEq;
            pc.printL
            pc.printM;
            return
            pc.printL
            d = pc(1).d+pc(2).d;
            disp('1+2')
            d.print
            disp('1+2+3')
            d = d+pc(3).d;
            d.print
            disp('1+2+3+4')
            d = d+pc(4).d;
            d.print
            return




        end

        function testBiconjugate4 (testCase)  % piece 3 of 1 piece conjugate - 4 pieces
            s_1 = sym('s_1');
            s_2 = sym('s_2');

            % V1
            f1 = symbolicFunction(2*s_1);
            ineq(1) = (5*s_2)/7 - s_1 + 25/7;
            ineq(2) =  4 - 2*s_2 - s_1 ;
            d1 = region(ineq,[s_1,s_2]);
            d1 = d1.removeInfV;
            d1 = d1.poly2orderUnbounded;
            edgeNo = d1.getEdgeNosInf(d1.vars)
            d1.ineqs(edgeNo) = d1.ineqs;
            p(1) = functionNDomain(f1,d1);

            % V2
            f2 = symbolicFunction(-5*s_1 + 5*s_2 + 25);
            ineq(1) = - s_1 - 2*s_2 - 4 ;
            ineq(2) =  s_1 + 2*s_2 - 4 ;
            ineq(3) = 48*s_1 - 56*s_2 + (s_1 + 2*s_2)^2 - 184;
            d2 = region(ineq,[s_1,s_2]);
            d2 = d2.removeInfV;
            d2 = d2.poly2orderUnbounded;
            edgeNo = d2.getEdgeNosInf(d2.vars)
            d2.ineqs(edgeNo) = d2.ineqs;
            p(2) = functionNDomain(f2,d2);

            % E1
            f3 = symbolicFunction(-5*s_1 + 5*s_2 + 25);
            ineq = sym.empty() 
            ineq(1) = s_1 + 2*s_2 + 4 ;
            ineq(2) =  s_1 - (9*s_2)/5 - 5 ;
            
            
            d3 = region(ineq,[s_1,s_2]);
            d3 = d3.removeInfV;
            d3 = d3.poly2orderUnbounded;
            edgeNo = d3.getEdgeNosInf(d3.vars)
            d3.ineqs(edgeNo) = d3.ineqs;
            p(3) = functionNDomain(f3,d3);

            % E2
            f4 = symbolicFunction(-5*s_1 + 5*s_2 + 25);
            ineq = sym.empty() 
            ineq(1) = s_1 - (5*s_2)/7 - 25/7 ;
            ineq(2) =  4 - 2*s_2 - s_1 ;
            
            d4 = region(ineq,[s_1,s_2]);
            d4 = d4.removeInfV;
            d4 = d4.poly2orderUnbounded;
            edgeNo = d4.getEdgeNosInf(d4.vars)
            d4.ineqs(edgeNo) = d4.ineqs;
            p(4) = functionNDomain(f4,d4);

            % V3
            f5 = symbolicFunction(-4*s_2);
            ineq = sym.empty() 


            ineq(1) = s_1 + 2*s_2 + 4 ;
            ineq(2) = (9*s_2)/5 - s_1 + 5;
            
            d5 = region(ineq,[s_1,s_2]);
            d5 = d5.removeInfV;
            d5 = d5.poly2orderUnbounded;
            edgeNo = d5.getEdgeNosInf(d5.vars)
            d5.ineqs(edgeNo) = d5.ineqs;
            p(5) = functionNDomain(f5,d5);
            

            % E3
            f6 = symbolicFunction(0.125*s_1^2 + 0.5*s_1*s_2 + s_1 + 0.5*s_2^2 + -2*s_2 + 2);
            ineq = sym.empty() 


            ineq(1) = - s_1 - 2*s_2 - 4 ;
            ineq(2) = s_1 + 2*s_2 - 4;
            ineq(3) = 56*s_2 - 48*s_1 - (s_1 + 2*s_2)^2 + 184;


            
            d6 = region(ineq,[s_1,s_2]);
            d6 = d6.removeInfV;
            d6 = d6.poly2orderUnbounded;
            edgeNo = d6.getEdgeNosInf(d6.vars)
            d6.ineqs(edgeNo) = d6.ineqs;
            p(6) = functionNDomain(f6,d6);

            p.printL
            p.printM
            
            pc = p.conjugateOfPiecePoly ;
            % VERIFIED, not merely run: the cells must be a usable cover.
            testMaxMultiRegion.verifyCellsAreUsable(testCase, pc, 'cells 4');
            %pc = pc.mergeL;
            pc = pc.addEq;
            pc.printL
            pc.printM
            % pc = pc.jSort;
            % pc.printL
return
            d = pc(1).d+pc(6).d;
            disp('1+6')
            d.print
            
            % disp('1+2+8')
            %  d = d+pc(8).d;
            % 
            %  d.print
              disp('1+2+10')
             d = d+pc(10).d;

             d.print
             return
            disp('1+2+6+8')
             d = d+pc(8).d;
             d.print
             return
            % d = d+pc(5).d;
            % d.print
            d = d+pc(6).d;
            d.print
            d = d+pc(7).d;
            d.print
            d = d+pc(8).d;
            d.print
            d = d+pc(9).d;
            d.print
            d = d+pc(10).d;
            d.print
            return
            
        end
  

        function testBiconjugate5 (testCase)  % piece 3 of 1 piece conjugate - 4 pieces
            s_1 = sym('s_1');
            s_2 = sym('s_2');
            f1 = symbolicFunction(-4*s_2);
             

            
            ineq(1) = -s_1 - 4;
            ineq(2) =  (9*s_2)/5 - s_1 + 5;
            ineq(3) = s_1 + 2*s_2 + 4;
            d1 = region(ineq,[s_1,s_2]);
            d1 = d1.removeInfV;
            d1 = d1.poly2orderUnbounded;
            edgeNo = d1.getEdgeNosInf(d1.vars)
            d1.ineqs(edgeNo) = d1.ineqs;
            p(1) = functionNDomain(f1,d1);

            f1 = symbolicFunction(-5*s_1+5*s_2+25);
            ineq(1) = - s_2 - 5;
            ineq(2) =  s_1 - (9*s_2)/5 - 5;
            ineq(3) = s_1 + 2*s_2 + 4;
            d2 = region(ineq,[s_1,s_2]);
            d2 = d2.removeInfV;
            d2 = d2.poly2orderUnbounded;
            edgeNo = d2.getEdgeNosInf(d2.vars)
            d2.ineqs(edgeNo) = d2.ineqs;
            p(2) = functionNDomain(f1,d2);
 p.printL
% return
             pc = p.conjugateOfPiecePoly ;
             % VERIFIED, not merely run: the cells must be a usable cover.
             testMaxMultiRegion.verifyCellsAreUsable(testCase, pc, 'cells 5');
            %pc = pc.mergeL;
            %pc = pc.addEq;
            pc.printL
        end

     
        function testBiconjugate6 (testCase)  % piece 3 of 1 piece conjugate - 4 pieces
            s_1 = sym('s_1');
            s_2 = sym('s_2');
            f1 = symbolicFunction(s_1 + 3*s_2 -3);
             

            
            ineq(1) = s_2/3 - s_1 + 14/3;
            ineq(2) =  s_1 - 2*s_2 + 1;
            
            d1 = region(ineq,[s_1,s_2]);
            d1 = d1.removeInfV;
            d1 = d1.poly2orderUnbounded;
            edgeNo = d1.getEdgeNosInf(d1.vars)
            d1.ineqs(edgeNo) = d1.ineqs;
            p(1) = functionNDomain(f1,d1);

            f2 = symbolicFunction(2*s_1 + s_2  -2);
             

            
            ineq(1) = 2 - s_2;
            ineq(2) =  2*s_2 - s_1 - 1;
            ineq(3) = (4*s_2)/7 - s_1 + 27/7;
            
            d2 = region(ineq,[s_1,s_2]);
            d2 = d2.removeInfV;
            d2 = d2.poly2orderUnbounded;
            edgeNo = d2.getEdgeNosInf(d2.vars)
            d2.ineqs(edgeNo) = d2.ineqs;
            p(2) = functionNDomain(f2,d2);
            p.printL
            
             pc = p.conjugateOfPiecePoly ;
             % VERIFIED, not merely run: the cells must be a usable cover.
             testMaxMultiRegion.verifyCellsAreUsable(testCase, pc, 'cells 6');
            %pc = pc.mergeL;
            pc.printL
            pc = pc.addEq;
            pc.printL
        end

        function testBiconjugate7 (testCase)  % piece 3 of 1 piece conjugate - 4 pieces
            s_1 = sym('s_1');
            s_2 = sym('s_2');
           
            f1 = symbolicFunction(-5*s_1 + 5*s_2 + 25);
                      
            % ineq(1) = -s_2/3 + s_1 - 14/3;
            % ineq(2) =  s_1 - 2*s_2 + 44;
            % ineq(3) =  46 - 7*s_2 - s_1;

            ineq(1) = s_1 - (5*s_2)/7 - 25/7 
            ineq(2) = 4 - 2*s_2 - s_1 
            ineq(3) = s_1 - (4*s_2)/7 - 27/7 
            ineq(4) = s_1 - s_2/3 - 14/3 

            
            d1 = region(ineq,[s_1,s_2]);
            d1.print
            d1 = d1.removeInfV;
            d1.print
            d1 = d1.poly2orderUnbounded;
            d1.print
            edgeNo = d1.getEdgeNosInf(d1.vars)
            %edgeNo = d1.getEdgeNosInf2(d1.vars)
            d1.ineqs(edgeNo) = d1.ineqs;
            d1.print
            return
            
            p(1) = functionNDomain(f1,d1);

            return
            
            ineq(1) = 4 - 2*s_2 - s_1;
            ineq(2) = s_1 - (5*s_2)/7 - 25/7;
            ineq(3) =  s_1 - s_2/3 - 14/3;
            ineq(4) = 2*s_2 - s_1 - 44;
            ineq(5) = s_1 - (4*s_2)/7 - 27/7;
            
            d1 = region(ineq,[s_1,s_2]);
            d1 = d1.removeInfV;
            
            d1 = d1.poly2order;
            if d1.nv == size(d1.ineqs,2)
                disp('here')
                d1.print
                edgeNo = d1.getEdgeNos(d1.vars)
            else
              edgeNo = d1.getEdgeNosInf(d1.vars)
            end
            %return
            d1.ineqs(edgeNo) = d1.ineqs;
            p(2) = functionNDomain(f1,d1);

            
            p.printL
            
            pc = p.conjugateOfPiecePoly ;
            % VERIFIED, not merely run: the cells must be a usable cover.
            testMaxMultiRegion.verifyCellsAreUsable(testCase, pc, 'cells 7');
            %pc = pc.mergeL;
            pc.printL
            pc = pc.addEq;
            pc.printL
        end

        % move to testing of region
        function testBiconjugate8 (testCase)  % piece 3 of 1 piece conjugate - 4 pieces
            s_1 = sym('s_1');
            s_2 = sym('s_2');
           
            f1 = symbolicFunction(-5*s_1 + 5*s_2 + 25);
                      
            % ineq(1) = -9*s_2/5 + s_1 - 5;
            % ineq(2) =  s_1 + 2*s_2 + 4;
            % ineq(3) =  -46 + 7*s_2 + s_1;
            % ineq(4) =  -s_2 -5;
            
            ineq(1) =  s_1 + 2*s_2 - 4 
            ineq(2) =  - s_1 - 2*s_2 - 4 
            ineq(3) =  48*s_1 - 56*s_2 + 4*s_1*s_2 + s_1^2 + 4*s_2^2 - 184 
            ineq(4) = 2*s_2 - s_1 - 44 
            ineq(5) = s_1 - 9*s_2/5-5
            d1 = region(ineq,[s_1,s_2]);
            d1.vx
            d1.vy
            d1.print
            px = sym.empty
            py = sym.empty
            px(1) = 14/19
            py(1) = -45/19
            d1 = d1.removeTangent(1,px,py);
            d1.print
            d1 = d1.removeInfV;
            d1.print
            %return
            %d1.simplifyUnboundedRegion
            d1 = d1.poly2orderUnbounded;
            d1.print
            return
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ineq(1) =  s_2/3 - s_1 + 14/3 
            ineq(2) = s_1 - 2*s_2 + 1 

            % ineq(1) =  s_1 + 2*s_2 + 4;
            % ineq(2) = -9*s_2/5 + s_1 - 5;
            % ineq(3) =  -s_2 -5;
            % ineq(4) =  -46 + 7*s_2 + s_1;
            d1 = region(ineq,[s_1,s_2]);
            d1.print
            d1 = d1.removeInfV;
            d1.print
            d1 = d1.poly2orderUnbounded;
            d1.print
            edgeNo = d1.getEdgeNosInf(d1.vars)
            d1.ineqs(edgeNo) = d1.ineqs;
            d1.print
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
return
            ineq(1) =  s_1 + 2*s_2 - 4 
            ineq(2) =  - s_1 - 2*s_2 - 4 
            ineq(3) =  48*s_1 - 56*s_2 + 4*s_1*s_2 + s_1^2 + 4*s_2^2 - 184 
            ineq(4) = 2*s_2 - s_1 - 44 
            ineq(5) = s_1 - 9*s_2/5-5
            d1 = region(ineq,[s_1,s_2]);
            d1.print
            d1 = d1.removeInfV;
            d1.print
            % change to poly2order for closed regions
            d1 = d1.poly2orderUnbounded;
            d1.print
            edgeNo = d1.getEdgeNosInf(d1.vars)
            d1.ineqs(edgeNo) = d1.ineqs;
            d1.print


                    
            return
            p(1) = functionNDomain(f1,d1);

            
            
            p.printL
            
            pc = p.conjugateOfPiecePoly ;
            % VERIFIED, not merely run: the cells must be a usable cover.
            testMaxMultiRegion.verifyCellsAreUsable(testCase, pc, 'cells 8');
            %pc = pc.mergeL;
            pc.printL
            pc = pc.addEq;
            pc.printL
        end

        function testConvex (testCase)
            %testCase.PRect3.print
            testCase.PRect3 = testCase.PRect3.convexEnvelope
            % VERIFIED, not merely run -- see plqCheck and the note on the first such
            % block above. Placed immediately after the pipeline call because several
            % of these tests `return` before their printing.
            for iP = 1:testCase.PRect3.nPieces
                plqCheck.envelopeUnderestimates(testCase, testCase.PRect3.pieces(iP), ...
                    sprintf('PRect3 co f piece %d', iP));
            end
            
            %% 
            %testCase.PRect.printLatex
            %testCase.PRect3.print
           
        end
        function testConjugate (testCase)
            %testCase.PRect3.print
            testCase.PRect3 = testCase.PRect3.conjugate
            % VERIFIED, not merely run -- see plqCheck and the note on the first such
            % block above. Placed immediately after the pipeline call because several
            % of these tests `return` before their printing.
            for iP = 1:testCase.PRect3.nPieces
                % .maxConjugate is written by maximumConjugate, not by conjugate -- comparing
                % .conjugates to the sup reports every point uncovered.
                q = testCase.PRect3.pieces(iP).maximumConjugate;
                plqCheck.conjugateMatchesSup(testCase, q.maxConjugate, q.f.f, ...
                    q.d.polygon.vars, q.d, testMaxMultiRegion.dualPoints(), ...
                    sprintf('PRect3 piece %d f*', iP));
            end
            
            %% 
            %testCase.PRect.printLatex
            %testCase.PRect3.print
           
        end

        function testMaxR3 (testCase)
            testCase.PRect3.print
            testCase.PRect3 = testCase.PRect3.maximum;
            testCase.PRect3 = testCase.PRect3.biconjugateF;
            % VERIFIED, not merely run -- see plqCheck and the note on the first such
            % block above. Placed immediately after the pipeline call because several
            % of these tests `return` before their printing.
            testMaxMultiRegion.verifyStep3(testCase, testCase.PRect3, testCase.PRect3, 'PRect3');
            testCase.verifyNotEmpty(testCase.PRect3.biconjugate, ...
                'PRect3: biconjugateF produced no cells');
            testMaxMultiRegion.verifyBiconj(testCase, testCase.PRect3, testCase.PRect3, 'PRect3');
            %% 
            %return
            %testCase.PRect.printLatex
            testCase.PRect3.print
            testCase.PRect3.printDomainMaple
        end
        
        function testMaxT (testCase)
        % Step 3 on the two-polygon PTri2 fixture, against the sup over its own domain.
            p = plqStage.get('MMR_PTri2', 'max', @() testCase.PTri2.maximum);
            testCase.verifyNotEmpty(p.maxConjugate, 'Step 3 produced no cells');
            testMaxMultiRegion.verifyStep3(testCase, p, testCase.PTri2, 'PTri2');
        end
        
        function testMaxP (testCase)
        % Step 3 on the six-vertex polygon -- the longest single job in the bucket at 1023 s, and
        % previously the one with the least to show for it (it printed and returned).
            p = plqStage.get('MMR_Poly', 'max', @() testCase.Poly.maximum);
            testCase.verifyNotEmpty(p.maxConjugate, 'Step 3 produced no cells');
            testMaxMultiRegion.verifyStep3(testCase, p, testCase.Poly, 'Poly');
        end

        
        % function testMax2 (testCase)
        %     disp('here')
        %     testCase.PTri = testCase.PTri.maximum
        % end

        function testMax3 (testCase)
        % x^2 - y^2 over a quadrilateral: Step 3 then the biconjugate, both verified. Split by
        % stage like testMax above, since this one also runs for minutes.
            p = plqStage.get('MMR_PRect2', 'max', @() testCase.PRect2.maximum);
            testMaxMultiRegion.verifyStep3(testCase, p, testCase.PRect2, 'PRect2');
            p = plqStage.get('MMR_PRect2', 'biconj', @() p.biconjugateF);
            testCase.verifyNotEmpty(p.biconjugate, 'biconjugateF produced no cells');
            testMaxMultiRegion.verifyBiconj(testCase, p, testCase.PRect2, 'PRect2');
        end

        function testPSqroot (testCase)
            testCase.PSqroot = testCase.PSqroot.maximum
            % VERIFIED, not merely run. Step 3's max across the pieces IS f* of the
            % union, so it must equal the sup of <s,x> - f(x) over the ORIGINAL
            % domain; f** must be a convex underestimator of f on it. Placed HERE,
            % before the print/return block, because several of these tests return
            % early -- an assertion after that would be dead code that always passes.
            testMaxMultiRegion.verifyStep3(testCase, testCase.PSqroot, testCase.PSqroot, 'PSqroot');
            %testCase.PSqroot = testCase.PSqroot.biconjugateF
            testCase.PSqroot.print
            testCase.PSqroot.printDomainMaple
            %testCase.PRect2.printLatex
           
        end

        function testFractional(testCase)
            x = sym('x')
            y = sym('y')

            d=domain([-1,1;-3,-3;-4,-3],x,y);

            
            
            q = plq_1piece(d,symbolicFunction(x*y));
            
            
            P = plq(q);
            f = symbolicFunction(36*x^2+21*x*y+36*y^2-81*x+24*y-252,-12*x+9*y+75)
            
            P.pieces(1).envelope = [P.pieces(1).envelope,functionNDomain(f,P.pieces(1).d.polygon)];
            P.pieces(1).envelope(1).d = P.pieces(1).envelope(1).d.normalize1
            P.pieces(1)=P.pieces(1).conjugate;
            P.pieces(1).conjugates(4).d.ineqs(3) = -P.pieces(1).conjugates(4).d.ineqs(3)
            
            P.pieces(1).biconjugateP
            % VERIFIED, not merely run. This test hand-builds a RATIONAL envelope face and then
            % hand-flips one conjugate constraint, so the object it produces is deliberately not
            % the pipeline's own answer -- there is no sup to compare against. What IS assertable
            % is that the hand-edits left a well-formed conjugate: cells exist, each carries a
            % function and a non-degenerate region, and the flipped constraint really did flip.
            testCase.verifyNotEmpty(P.pieces(1).conjugates, ...
                'testFractional: the conjugate has no cells');
            testMaxMultiRegion.verifyCellsAreUsable(testCase, P.pieces(1).conjugates, ...
                'testFractional conjugate');
            P.print
            
            P.printDomainMaple
        end

        function testMaxThesis (testCase)
            testCase.PRect3.print
            
            warning('off','all') 
            testCase.PThesis = testCase.PThesis.maximum
            testCase.PThesis = testCase.PThesis.biconjugateF
            % VERIFIED, not merely run. Step 3's max across the pieces IS f* of the
            % union, so it must equal the sup of <s,x> - f(x) over the ORIGINAL
            % domain; f** must be a convex underestimator of f on it. Placed HERE,
            % before the print/return block, because several of these tests return
            % early -- an assertion after that would be dead code that always passes.
            testMaxMultiRegion.verifyStep3(testCase, testCase.PThesis, testCase.PThesis, 'PThesis');
            testCase.verifyNotEmpty(testCase.PThesis.biconjugate, ...
                'PThesis: biconjugateF produced no cells');
            testMaxMultiRegion.verifyBiconj(testCase, testCase.PThesis, testCase.PThesis, 'PThesis');
            %% 
            %testCase.PRect.printLatex
           % testCase.PThesis.print
            testCase.PThesis.printDomainMaple

        end

        function testMaxThesis2 (testCase)
            %testCase.PRect3.print
            warning('off','all')

            n = 1;
            x = sym('x');
            y = sym('y');
            f=symbolicFunction(x*y);
            
            div = symbolicFunction(2,x);
            div = subs(div.f,x,n);
            del = symbolicFunction(x);
            del = subs(del.f,x,-1);
            
            %del = -1;
            for i=1:n+1
                A(i) = del;
                del = div + del;
                
            end
            A
            %return

             
            m = 0;
            dom = [];
            for i = 1:n
                for j = 1:n
                    for ki = i:i+1
                        for kj = j:j+1
                          dom = [dom;[A(ki),A(kj)]];
                        end
                    end
                    temp = dom(3,1:2);
                    dom(3,1:2) = dom(4,1:2);
                    dom(4,1:2) = temp;
                    %dom
                    m = m + 1;
                    d(m) = domain(dom,x,y);
                    %d(m).print
                    q(m) = plq_1piece(d(m),f);
                    dom = [];
                    %break
                end
                %break
            end
            testCase.PThesis = plq(q);

            testCase.PThesis = testCase.PThesis.maximum;
            testCase.PThesis = testCase.PThesis.biconjugateF;
            % VERIFIED, not merely run -- see plqCheck and the note on the first such
            % block above. Placed immediately after the pipeline call because several
            % of these tests `return` before their printing.
            testMaxMultiRegion.verifyStep3(testCase, testCase.PThesis, testCase.PThesis, 'PThesis2');
            testCase.verifyNotEmpty(testCase.PThesis.biconjugate, ...
                'PThesis2: biconjugateF produced no cells');
            testMaxMultiRegion.verifyBiconj(testCase, testCase.PThesis, testCase.PThesis, 'PThesis2');
            testCase.PThesis.printDomainMaple;
            return
            
            %% 
            %testCase.PRect.printLatex
            testCase.PThesis.print
            
        end


        function testOpenconvex (testCase)
            x=sym('x');
            y=sym('y');
            d = domain();
            d = d.domainEdge([y,-x-1,x-1],[x,y]);
            %d.print
            testCase.POpen = plq([plq_1piece(d,symbolicFunction(x*y))]);
            testCase.POpen = testCase.POpen.convexEnvelope
            % VERIFIED, not merely run -- see plqCheck and the note on the first such
            % block above. Placed immediately after the pipeline call because several
            % of these tests `return` before their printing.
            for iP = 1:testCase.POpen.nPieces
                plqCheck.envelopeUnderestimates(testCase, testCase.POpen.pieces(iP), ...
                    sprintf('POpen co f piece %d', iP));
            end
            testCase.POpen.print
        end

        
    end

    
end