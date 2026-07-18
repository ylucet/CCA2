classdef testcPLQ < matlab.unittest.TestCase

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

    methods (Test)

        function testRect (testCase)
            testCase.PRect = testCase.PRect.triangulate;

            % add to split 3 convex edge case
            testCase.PRect =  testCase.PRect.convexEnvelope;
             testCase.PRect.print;
             testCase.PRect.printDomainMaple;
            
           
        end

        function testRectConj (testCase)
            testCase.PRect = testCase.PRect.triangulate;
            testCase.PRect =  testCase.PRect.conjugate;
             testCase.PRect.print;
             testCase.PRect.printDomainMaple;
            
           
        end

        function testRectMax (testCase)
            testCase.PRect = testCase.PRect.triangulate;
            testCase.PRect =  testCase.PRect.maximum;
             testCase.PRect.print;
             testCase.PRect.printDomainMaple;
            
           
        end

        function testRectBiconj (testCase)


            
            testCase.PRect = testCase.PRect.triangulate;
            testCase.PRect =  testCase.PRect.maximum;
            testCase.PRect =  testCase.PRect.biconjugateF;
             % testCase.PRect.print;
             % testCase.PRect.printDomainMaple;
            
           
        end
        

        function testRect3CE (testCase)
            testCase.PRect3 = testCase.PRect3.triangulate;
            
             testCase.PRect3 = testCase.PRect3.convexEnvelope;
            
             testCase.PRect3.print;
             testCase.PRect3.printDomainMaple;
            
           
        end

        function testRect3Conj (testCase)
            testCase.PRect3 = testCase.PRect3.triangulate;
            
            testCase.PRect3 = testCase.PRect3.conjugate;
            %return
             testCase.PRect3.print;
             testCase.PRect3.printDomainMaple;
            
           
        end

         function testRect3Max (testCase)
            testCase.PRect3 = testCase.PRect3.triangulate;
            
            testCase.PRect3 = testCase.PRect3.maximum;
            %return
             testCase.PRect3.print;
             testCase.PRect3.printDomainMaple;
            
           
         end

         function testRect3Biconj (testCase)
            testCase.PRect3 = testCase.PRect3.triangulate;
            testCase.PRect3 =  testCase.PRect3.maximum;
            testCase.PRect3 =  testCase.PRect3.biconjugateF;
             testCase.PRect3.print;
             testCase.PRect3.printDomainMaple;
            
           
        end
    end
end