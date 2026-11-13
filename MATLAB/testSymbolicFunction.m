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
        % `tangent(x0,y0)` is the tangent LINE to the curve {f = 0} at (x0,y0), so it is pinned by
        % two properties and nothing else:
        %     t(x0,y0) = 0                      the line passes through the point
        %     grad t  parallel to  grad f       it has the curve's own direction there
        % Both are checked at every root of the conic-and-line system, for the conic AND for its
        % NEGATION -- the same curve, so the normalised tangent must come out identical. That last
        % one is the real content: `tangent` divides by partial derivatives, and negating f flips
        % every sign in the quotient, which is exactly where a sign slip hides.
        %
        % Previously this test computed all of the above and displayed it.
            x = sym('x');
            y = sym('y');
            f2 = symbolicFunction(148*x - 196*y + (x + 7*y)^2 - 684);
            f3 = symbolicFunction(x - (9*y)/5 - 5);
            [tx, ty] = solve(f3.f == 0, f2.f == 0);
            testCase.verifyGreaterThanOrEqual(numel(tx), 1, ...
                'the line must meet the conic; with no root nothing below is tested');

            for k = 1:numel(tx)
                p = [tx(k), ty(k)];
                for s = [1 -1]
                    fk = f2;
                    if s < 0, fk = -f2; end
                    t = fk.tangent(tx(k), ty(k));
                    t = t.normalize1;

                    testCase.verifyEqual(double(subs(t.f, [x y], p)), 0, 'AbsTol', 1e-9, ...
                        sprintf('root %d, sign %+d: the tangent does not pass through the point', k, s));

                    gT = double([subs(diff(t.f, x), [x y], p),  subs(diff(t.f, y), [x y], p)]);
                    gF = double([subs(diff(fk.f, x), [x y], p), subs(diff(fk.f, y), [x y], p)]);
                    testCase.verifyGreaterThan(norm(gT), 1e-9, ...
                        sprintf('root %d, sign %+d: the tangent is degenerate (zero gradient)', k, s));
                    cross = gT(1)*gF(2) - gT(2)*gF(1);
                    testCase.verifyEqual(cross / (norm(gT)*norm(gF)), 0, 'AbsTol', 1e-9, ...
                        sprintf(['root %d, sign %+d: the tangent is not parallel to the curve ' ...
                                 'there (grad t = %s, grad f = %s)'], k, s, mat2str(gT,4), mat2str(gF,4)));

                    % The tangent is AFFINE: a "tangent line" that still carries a quadratic term
                    % is the failure mode a value check at the touch point cannot see.
                    testCase.verifyEqual(double(diff(t.f, x, 2)), 0, 'AbsTol', 1e-12, ...
                        sprintf('root %d, sign %+d: the tangent is not affine in x', k, s));
                    testCase.verifyEqual(double(diff(t.f, y, 2)), 0, 'AbsTol', 1e-12, ...
                        sprintf('root %d, sign %+d: the tangent is not affine in y', k, s));
                end

                % Same curve, opposite sign of f: after normalize1 the two tangents must be the
                % same line.
                fNeg = -f2;
                tPos = f2.tangent(tx(k), ty(k));    tPos = tPos.normalize1;
                tNeg = fNeg.tangent(tx(k), ty(k));  tNeg = tNeg.normalize1;
                testCase.verifyTrue(isAlways(simplify(tPos.f - tNeg.f) == 0), sprintf( ...
                    'root %d: tangent(f) = %s but tangent(-f) = %s -- the same curve must give the same line', ...
                    k, char(tPos.f), char(tNeg.f)));
            end
        end

        function tangentOfACurveMissingOneAmbientVariable(testCase)
        % Regression (2026-08-30): obj.vars is populated from symvar(obj.f) at construction, the
        % expression's OWN free variables -- not the ambient 2-variable dual space. A degenerate
        % conic that genuinely does not depend on one ambient variable (a repeated/parallel-line
        % pair like s1^2+101*s1+239=0, independent of s2) then has obj.vars of length 1, and the
        % old tangent(obj,x,y) crashed on obj.vars(2) with MATLAB:badsubscript. Reached in
        % production via region.removeTangent on doc/QuaConExample.md's 3-piece witness.
            s1 = sym('s_1'); s2 = sym('s_2');
            f = symbolicFunction(101*s1 + s1^2 + 239);
            t = f.tangent(0, 0, [s1 s2]);
            testCase.verifyTrue(isAlways(t.f == s1), ...
                'tangent to a curve independent of s2 is the vertical line s1=x');
            % Backward compatible: omitting vars still works when the function genuinely uses both.
            g = symbolicFunction(s1^2 + s2^2 - 4);
            t2 = g.tangent(2, 0);
            testCase.verifyTrue(isAlways(t2.f == s1 - 2));
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