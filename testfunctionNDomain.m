classdef testfunctionNDomain < matlab.unittest.TestCase

    properties
        r
    end
    
    methods(Test)

        function testMerge(testCase)
          x = sym('x');
          y = sym('y');
          l1 = [-x - 7*y - 4, x + 7*y-10, 196*y-148*x-(x+7*y)^2 + 684, 4-2*y-x ] 
          r1 = region(l1,[x,y]);
                  
          l2 = [(y/3)-x+ 14/3, 10-7*y-x, y-2]
          r2 = region(l2,[x,y]);
         
          l3 = [x+7*y+4,4-2*y-x]
          r3 = region(l3,[x,y]);

          l4 = [-y/3+x-14/3, 10-7*y-x, y-2, 5*y/7-x+25/7]
          r4 = region(l4,[x,y]);
         
         
          l5 = [x + 7*y-10, -196*y+148*x+(x+7*y)^2 - 684, 4-2*y-x, x - y/3-14/3,5*y/7-x+25/7 ] 
          r5 = region(l5,[x,y]);
          fL = functionNDomain.empty();
          fL = [fL,functionNDomain(symbolicFunction (x),r1)];
          fL = [fL,functionNDomain(symbolicFunction (x),r2)];
          fL = [fL,functionNDomain(symbolicFunction (x),r3)];
          fL = [fL,functionNDomain(symbolicFunction (x),r4)];
          fL = [fL,functionNDomain(symbolicFunction (x),r5)];
         
          
          fL.printL
          fL.printM
          %return
          fL2 = fL.mergeL;
          fL2.printL
         fL2.printM
        end

        function testUnique(testCase)
            x = sym('x');
          y = sym('y');
          l1 = [-x - 7*y - 4, x + 7*y-10, 196*y-148*x-(x+7*y)^2 + 684, 4-2*y-x ] 
          r1 = region(l1,[x,y]);
                  
          l2 = [(y/3)-x+ 14/3, 10-7*y-x, y-2]
          r2 = region(l2,[x,y]);
         
          l3 = [x+7*y+4,4-2*y-x]
          r3 = region(l3,[x,y]);

          l4 = [-y/3+x-14/3, 10-7*y-x, y-2, 5*y/7-x+25/7]
          r4 = region(l4,[x,y]);
         
         
          l5 = [x + 7*y-10, -196*y+148*x+(x+7*y)^2 - 684, 4-2*y-x, x - y/3-14/3,5*y/7-x+25/7 ] 
          r5 = region(l5,[x,y]);
          fL = functionNDomain.empty();
          fL = [fL,functionNDomain(symbolicFunction (x),r1)];
          fL = [fL,functionNDomain(symbolicFunction (x),r2)];
          fL = [fL,functionNDomain(symbolicFunction (x),r3)];
          fL = [fL,functionNDomain(symbolicFunction (x),r4)];
          fL = [fL,functionNDomain(symbolicFunction (x),r5)];
          fL = [fL,functionNDomain(symbolicFunction (x),r1)];
          fL = [fL,functionNDomain(symbolicFunction (x),r1)];
         
          
          fL.printL
            fL =  fL.unique  
             fL.printL
             fL = [fL,functionNDomain(symbolicFunction (x^2),r1)];
             fL =  fL.unique  
             fL.printL
            

    end
    end
end