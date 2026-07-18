classdef testSymbolicFunction < matlab.unittest.TestCase

    methods (Test)
        function noParameter (testCase)
          f1 = symbolicFunction();
          
          testCase.verifyEqual(f1.getNum(),0);
          testCase.verifyEqual(f1.getDen(),1);
          testCase.verifyEqual(f1.f,0);
        end
        function oneParameter (testCase)
          x = sym('x'); 
          f1 = symbolicFunction(x^2);
          testCase.verifyEqual(f1.getNum(),x^2);
          testCase.verifyEqual(f1.getDen(),sym(1));
          testCase.verifyEqual(f1.f,x^2);
        end
        function twoParameter (testCase)
          x = sym('x'); 
          f1 = symbolicFunction(x^2,x+1);
          testCase.verifyEqual(f1.getNum(),x^2);
          testCase.verifyEqual(f1.getDen(),x+1);
          testCase.verifyEqual(f1.f,x^2/(x+1));
        end

        function testGetf (testCase)
            x = sym('x');
            y = sym('y');
            f = 148*x - 196*y + (x + 7*y)^2 - 684;
            fs = symbolicFunction(f);
            testCase.verifyEqual(f,fs.getF);
        
        end

       

         function testgetNum (testCase)
            x = sym('x');
            y = sym('y');
            f = 148*x - 196*y + x + 7*y^2 - 684;
            fs = symbolicFunction(f);
            testCase.verifyEqual(f,fs.getNum);
         end
        
         function testgetNum2 (testCase)
            x = sym('x');
            y = sym('y');
            f = 148*x - 196*y + x + 7*y^2 - 684;
            fs = symbolicFunction(f,x^2);
            testCase.verifyEqual(f,fs.getNum);
         end

         function testgetNum3 (testCase)
            x = sym('x');
            y = sym('y');
            f = 148*x - 196*y + x + 7*y^2 - 684;
            fs = symbolicFunction(0,x^2);
           
            testCase.verifyEqual(sym(0),fs.getNum);
         end

          function testgetNum4 (testCase)
            % x = sym('x');
            % y = sym('y');
            % f = 148*x - 196*y + x + 7*y^2 - 684;
            fs = symbolicFunction(0);
            testCase.verifyEqual(0.0,fs.getNum);
         end

         function testgetDen (testCase)
            x = sym('x');
            y = sym('y');
            f = 148*x - 196*y + x + 7*y^2 - 684;
            fs = symbolicFunction(f);
            
            testCase.verifyEqual(sym(1),fs.getDen);
         end
        
         function testgetDen2 (testCase)
            x = sym('x');
            y = sym('y');
            f = 148*x - 196*y + x + 7*y^2 - 684;
            fs = symbolicFunction(f,x^2);
            testCase.verifyEqual(x^2,fs.getDen);
         end

         function testgetDen3 (testCase)
            x = sym('x');
            y = sym('y');
            f = 148*x - 196*y + x + 7*y^2 - 684;
            fs = symbolicFunction(0,x^2);
           
            testCase.verifyEqual(sym(1),fs.getDen);
         end

         function testgetDen4 (testCase)
            % x = sym('x');
            % y = sym('y');
            % f = 148*x - 196*y + x + 7*y^2 - 684;
            fs = symbolicFunction(0);
            testCase.verifyEqual(1,fs.getDen);
         end


        
        function testPlus (testCase)
            x = sym('x');
            y = sym('y');
            f1 = symbolicFunction(x + y);
            f2 = symbolicFunction(x - y);
            testCase.verifyEqual(f1+f2,symbolicFunction(2*x));

        end

        function testMinus (testCase)
            x = sym('x');
            y = sym('y');
            f1 = symbolicFunction(x + y);
            f2 = symbolicFunction(x - y);
            testCase.verifyEqual(f1-f2,symbolicFunction(2*y));

        end

        function testuMinus (testCase)
            x = sym('x');
            y = sym('y');
            f1 = symbolicFunction(x + y);
            f2 = symbolicFunction(x - y);
            testCase.verifyEqual(-f2,symbolicFunction(-x+y));

        end

        function testmTimes (testCase)
            x = sym('x');
            y = sym('y');
            f1 = symbolicFunction(x + y);
            f2 = symbolicFunction(x - y);
            g = f1*f2;
            
            testCase.verifyEqual(f1*f2,symbolicFunction((x+y)*(x-y)));

        end

        function testdfdx (testCase)
            x = sym('x');
            y = sym('y');
            f1 = symbolicFunction(x^2 + 2*y);
            f1.dfdx(x);
            testCase.verifyEqual(f1.dfdx(x),symbolicFunction(2*x));
           
            
            % really stupid 
            testCase.verifyEqual(f1.dfdx(y) ==            symbolicFunction(2), true);

        end
 
        
        function testTangent(testCase)
            x = sym('x');
            y = sym('y');
            f2 = symbolicFunction(148*x - 196*y + (x + 7*y)^2 - 684)
            %f2.tangent(3.159091,-1.022727)
            f3 = symbolicFunction(x - (9*y)/5 - 5)
            [tx,ty] = solve(f3.f==0,f2.f==0)
            
            t = f2.tangent(tx,ty)
            t = t.normalize1
            f2 = -f2  
            t = f2.tangent(tx,ty)
            t = t.normalize1
        end
        
        function testGradient(testCase)
           x = sym('x');
           y = sym('y');
           f = symbolicFunction(x*y)
           g = f.gradient(f.getVars)
           testCase.verifyEqual(g(1).f,y);
           testCase.verifyEqual(g(2).f,x);
        end


         % function testGetf (testCase)
        % end

    end

    
end