classdef PLQVCTest < matlab.unittest.TestCase
    % HausdorffTest tests computation of Hausdorff distance between
    % 2 convex polytopes
    
    methods (Test)                
        function testNorm1(testCase)
            p = PLQVC.oneNorm();
            actP = p.P;
            expP = {[-1 4]; [-2 1]; [-3 2]; [-4 3]};
            testCase.verifyEqual(actP,expP);
        end      
        function testIsDomBounded1(testCase)
            p = PLQVC.oneNorm();
            testCase.verifyFalse(p.isDomBounded,false);
        end  
        function testIsDomBounded2(testCase)
            p = PLQVC.oneNormConjugate();%indicator of unit square
            testCase.verifyTrue(p.isDomBounded,true);
        end        
        function testIsDomBounded3(testCase)
            p = PLQVC.energy();%indicator of unit square
            testCase.verifyFalse(p.isDomBounded,false);
        end
        function testIsDomBounded4(testCase)
            V = [0 0; 0 1; 1 0];
            E = [1 2 0; 1 3 0];
            f = [0.5 0 0.5 0 0 0];
            F = [0 1;1 0];
            p = PLQVC(V,E,f,F); 
            testCase.verifyFalse(p.isDomBounded,false);
        end
        function testDisp0(testCase)
            p = PLQVC.energy();
            s = evalc('disp(p)');
            testCase.verifyEqual(s(1:end-1),' 0.5 x^2 + 0.5 y^2  if (x,y) in P( 1 )');
        end         
        function testDisp1(testCase)
            p = PLQVC.cubic1();
            s = evalc('disp(p)');
            expS=sprintf('%s\n',...
                '"- x^3 - y^3 +  x^2 +  y^2  if (x,y) in P( 1 )"',...
                '"- x^3 + 2 y^3 +  x^2 + 2 y^2  if (x,y) in P( 2 )"',...
                '" 3 x^3 + 2 y^3 + 2 x^2 + 2 y^2  if (x,y) in P( 3 )"',...
                '" 3 x^3 - y^3 + 2 x^2 +  y^2  if (x,y) in P( 4 )"');
            testCase.verifyEqual(strrep(s(1:end-1)," ",""),strrep(expS," ",""));
        end   
        function testDisp2(testCase)
            p = PLQVC.oneNormConjugate();
            s = evalc('disp(p)');
            testCase.verifyEqual(s(1:end-1),'0 if (x,y) in P( 1 )');
        end           
        function testCreateP1(testCase)
            p = PLQVC.cubic1();
            P = {[-1 4];[-2 1];[-3 2];[-4 3]};
            testCase.verifyEqual(p.P,P);
        end      
        function testCreateP2(testCase)
            V = [1 0;0 0;0 1;1 1;2 0;2 1];
            E = [1 2 1;2 3 1;3 4 1;4 1 1;1 5 0;4 6 0];
            f = [1 0 0 1 0 0 1 0 0 0;
                 1 0 1 1 0 0 0 0 0 0];
            P = {[1 2 3 4]; [-5 -4 6]};
            F = [0 1; 0 1; 0 1; 2 1; 2 0; 0 2];
            p = PLQVC(V,E,f,F);
            testCase.verifyEqual(p.P,P);
        end
        function testCreateP3(testCase)
            V = [0 0; 0 1; -1 0; 1 0];
            E = [1 2 0; 1 3 0; 1 4 0];
            f = [-1 0 0 1 0 0 0 0 0 0; 
                  1 0 0 1 0 0 0 0 0 0];
            P = {[-1 2];[-3 1]};
            F = [1 2;0 1;2 0];
            p = PLQVC(V,E,f,F);
            testCase.verifyEqual(p.P,P);
        end
        function testCreateP4(testCase)
            V = [0 0; 0 1; 1 0; 1 1; 2 0; 2 1];
            E = [1 2 1; 3 4 1; 5 6 1; 1 3 1; 2 4 1; 3 5 1; 4 6 1];
            f = [1 0 1 1 0 0 0 0 1 0; 
                 1 0 0 1 0 1 1 0 0 0];
            P = {[2 7 -3 -6 ];[1 5 -2 -4]};
            F = [0 2;2 1;1 0;2 0;0 2;1 0;0 1];
            p = PLQVC(V,E,f,F);
            testCase.verifyEqual(p.P,P);
        end
        function testCreateP5(testCase)
            V = [1 0; 1 2; 2 0; 0 0];
            E = [1 2 0; 1 3 0; 1 4 0];
            f = [3 1 1 1 0 0;
                 4 1 1 0 0 0];
            F = [1 2;2 0;0 1];
            P = {[-1 3];[-2 1]};
            p = PLQVC(V,E,f,F);
            testCase.verifyEqual(p.P,P);
        end        
        function testCreateP6(testCase)
            V = [0 0; 0 1; -1 0; 1 0];
            E = [1 2 0; 1 3 0; 1 4 0];
            f = [1 0 1 0 0 0;
                 2 1 1 0 0 0];
            F = [1 2;0 1;2 0];
            P = {[-1 2];[-3 1]};
            p = PLQVC(V,E,f,F);
            testCase.verifyEqual(p.P,P);
        end       
        function testCreateP7(testCase)
            V = [-1 -1;1 -1; 1  1; -1 1];
            E = [1 2 1; 2 3 1; 3 4 1; 4 1 1];
            F = [1 0;1 0;1 0;1 0];
            f = [1 2 1 0 0 0];
            P = {[-1 -4 -3 -2]};
            p = PLQVC(V,E,f,F);
            testCase.verifyEqual(p.P,P);
        end               
        function testCreateP9(testCase)
            V = [0 0;0 1;0 -1];
            E = [1 2 0;1 3 0];
            f = [1 1 1 0 0 0;
                 1 0 1 0 0 0];
            F = [1 2;2 1];
            P = {[-1 2];[-2 1]};
            p = PLQVC(V,E,f,F);
            testCase.verifyEqual(p.P,P);
        end    
        function testCreateP10(testCase)
            V = [0 0;1 0;0.5 1;];
            E = [1 2 1;2 3 1;1 3 1];
            F = [1 0; 1 0; 0 1];
            f = [1 2 1 0 0 0];
            P = {[-1 3 -2]};
            p = PLQVC(V,E,f,F);
            testCase.verifyEqual(p.P,P);
        end    
        function testCreateP11(testCase)
            V = [0 0;1 0.1;.5 .6];
            E = [1 2 0;1 3 0];
            F = [1 0; 0 1];
            f = [1 2 1 0 0 0];
            P = {[-1 2]};
            p = PLQVC(V,E,f,F);
            testCase.verifyEqual(p.P,P);        
        end
        function testCreateP12(testCase)
            V = [0 0; 1 0; 1 1; 0 1; 2 0; 2 1];
            E = [1 2 1; 3 4 1; 4 1 1; 2 3 1; 2 5 0; 3 6 0];
            f = [1 0 0 1 0 0 1 0 0 0;
                 1 0 1 1 0 0 0 0 0 0];
            P = {[-1 -3 -2 -4];[-5 4  6]};
            F = [1 0;1 0;1 0;1 2;2 0;0 2];
            p = PLQVC(V,E,f,F);
            testCase.verifyEqual(p.P,P);        
        end
        function testEvalPoly(testCase)
           %use symbolic computation to verify in all cases both evaluation methods give
           %the same results
           syms x y
           assume(x,'real');assume(y,'real');
           X=[x;y];
           c = sym('c',[1,10]); assume(c,'real');
           c1 = simplify(PLQVC.evalMatrixForm(c,X'));
           c2 = PLQVC.evalPoly(c,X');
           res = double(simplify(c1-c2));
           testCase.verifyEqual(res,0);
        end
        function testEval(testCase)
            p = PLQVC.oneNorm();
            X=[0,0;1,1;-1,-1;-1,1;1,-1;1,2];
            expZ=[0;2;2;2;2;3];
            expRegion=[0;3;1;2;4;3];
            [fVal,region] = p.eval(X);
            testCase.verifyEqual(fVal,expZ);
            testCase.verifyEqual(region,expRegion);
            p = PLQVC.energy();
            X=[0,0; 1,1; -1,-1; -1,1; 1,-1; 1,2];
            expZ=[0;1;1;1;1;2.5];
            expRegion=[1;1;1;1;1;1];
            [fVal,region] = p.eval(X);
            testCase.verifyEqual(fVal,expZ);
            testCase.verifyEqual(region,expRegion);
        end
        function testEval1(testCase)
            Ex = PLQVC.examples();
            p = Ex{5};
            P = {[-1 4];[-2 1];[-3 2];[-4 3]};
            testCase.verifyEqual(p.P,P);
            X = [1  1;1  0;0  1;-1 -1;4  5;3 -1;-5  2;-9 -1];
            fVal = eval(p,X);
            testCase.verifyEqual(fVal,[4 2 2 2 82 19 33 82]');
            p = Ex{6};
            P = {[-1 4];[-2 1];[-3 2];[-4 3]};
            testCase.verifyEqual(p.P,P);
            X = [1  1;1  0;0  1;-1 -1];
            fVal = eval(p,X);
            testCase.verifyEqual(fVal,-[4 2 2 2]');
            p = Ex{7};
            P = {[-1 -4 -3 -2]};
            testCase.verifyEqual(p.P,P);
            X = [1  1;2  2;0  0;-0.5  0];
            fVal = eval(p,X);
            testCase.verifyEqual(fVal,[4 Inf 0 .25]');
            p = Ex{8};
            P = {[-1 3];[-2 1]};
            testCase.verifyEqual(p.P,P);
            X = [1  1;-1 -1;-1  1];
            fVal = eval(p,X);
            testCase.verifyEqual(fVal,[3 Inf 2]');
            p = Ex{9};
            P = {[-1 3 -2]};
            testCase.verifyEqual(p.P,P);
            X = [.5 .5;1  0;-1  1];
            fVal = eval(p,X);
            testCase.verifyEqual(fVal,[1 1 Inf]');
            p = Ex{10};
            P = {[-1 2]};
            testCase.verifyEqual(p.P,P);
            X = [.5  .5;1  .2;10  10;1   0];
            fVal = eval(p,X);
            testCase.verifyEqual(fVal,[1 1.44 400 Inf]');
            p = Ex{11};
            P = {[-2 1];[-1 2]};
            testCase.verifyEqual(p.P,P);
            X = [1  1;1  0;0  1;-1 -1;4  5;3 -1;-5  2;-9 -1];
            fVal = eval(p,X);
            testCase.verifyEqual(fVal,[0.25 0.25 0.375 0.25 6.375 4.25 7.125 8.25]');
            p = Ex{12};
            P = {[-1 2]; [-2 1]};
            testCase.verifyEqual(p.P,P);
            X = [1  1;1  0;0  1;-1 -1;4  5;3 -1;-5  2;-9 -1];
            fVal = eval(p,X);
            testCase.verifyEqual(fVal,[4 1 1 5 81 4 69 253]');
            p = Ex{13};
            P = {[-1 4];[-2 1];[-3 2];[-4 3]};
            testCase.verifyEqual(p.P,P);
            X = [1  1;1  0;0  1;-1 -1;4  5;3 -1;-5  2;-9 -1];
            fVal = eval(p,X);
            testCase.verifyEqual(fVal,[6 3 3 4 271 47 166 812]');
            p = Ex{14};
            P = {[-1 -3 -2 -4];[-5  4  6]};
            testCase.verifyEqual(p.P,P);
            X = [1/2 1/2;1     0;0     1;2   0.5;-1    -1;4     5;3    -1;-5     2;-9    -1];
            fVal = eval(p,X);
            testCase.verifyEqual(fVal,[0.5 1 2 8.625 Inf Inf Inf Inf Inf]');
            p = Ex{15};
            P = {[1 2 3 4];[-5 -4 6]};
            testCase.verifyEqual(p.P,P);
            X = [1/2 1/2;1  0;0  1;2 0.5;-1 -1;4  5;3 -1;-5  2;-9 -1];
            fVal = eval(p,X);
            testCase.verifyEqual(fVal,[0.5 1 2 8.625 Inf Inf Inf Inf Inf]');
        end  
        function testEval2(testCase)
            p = PLQVC([1,1,1,1, 1,1,1, 1,1, 1]);
            X=[0,0;1,1;-1,-1;-1,1;1,-1;1,2];
            z1 = PLQVC.evalPoly(p.f,X);
            z2 = p.eval(X);
            testCase.verifyEqual(z1,z2);
            p = PLQVC([0,0,0,0, 1,0,1, 0,0, 42]);
            X=[0,0;10,0;0,10;-11,4;-5,-7;2.5,-4.56];
            z1 = PLQVC.evalPoly(p.f,X);
            z2 = p.eval(X);
            testCase.verifyEqual(z1,z2);
            p = PLQVC([-3,4.5,-64.8,42.89, 3.56,-2.55,9.45, 42.42,-56.346, 123.321]);
            X=[0,0;10,0;0,10;-11,4;-5,-7;2.5,-4.56];
            z1 = PLQVC.evalPoly(p.f,X);
            z2 = p.eval(X);
            testCase.verifyEqual(z1,z2);
        end
        function testEvalMatrixForm(testCase)
            p = PLQVC([1,1,1,1, 1,1,1, 1,1, 1]);
            X=[0,0; 1,1; -1,-1; -1,1; 1,-1; 1,2];
            z1 = PLQVC.evalPoly(p.f,X);
            z2 = PLQVC.evalMatrixForm(p.f,X);
            testCase.verifyEqual(z1,z2);
            p = PLQVC([0,0,0,0, 1,0,1, 0,0, 42]);
            X=[0,0;10,0;0,10;-11,4;-5,-7;2.5,-4.56];
            z1 = PLQVC.evalPoly(p.f,X);
            z2 = PLQVC.evalMatrixForm(p.f,X);
            testCase.verifyEqual(z1,z2);
            p = PLQVC([-3,4.5,-64.8,42.89, 3.56,-2.55,9.45, 42.42,-56.346, 123.321]);
            X=[0,0;10,0;0,10;-11,4;-5,-7;2.5,-4.56];
            z1 = PLQVC.evalPoly(p.f,X);
            z2 = PLQVC.evalMatrixForm(p.f,X);
            testCase.verifyEqual(z1,z2,'RelTol',1e-15);
        end        
        function testCreateDom_isDomBounded(testCase)
            p = PLQVC.energy;%finite everywhere             
            testCase.verifyEqual(p.dom.dim,2);
            testCase.verifyTrue(p.dom.isConvex,true);
            P2 = {PLQVC.examples(), PLQVC.examples2()};
            D2 = {2*ones(1,19),[0, 1, 1, 1, 1, 1, 1, 2,2,1]};
            isConvex2 = {true(1,19),[true,true,true, true, false, false, false, true, true, false]};
            isDB2 = {false(1,19),[true, true, false, false, true, false, false, false, false,true]};
            isDB2{1}(2)=true;isDB2{1}(7)=true;isDB2{1}(9)=true;isDB2{1}(17)=true;
            for i=1:2                
                Pall=P2{i};D=D2{i};isConvex=isConvex2{i};isDB=isDB2{i};
                k1=1;k2=length(Pall);%k1=k2;
                for j=k1:k2
                    %disp([i, j]);
                    testCase.verifyEqual(Pall{j}.dom.dim,D(j));
                    testCase.verifyEqual(Pall{j}.dom.isConvex,isConvex(j));
                    testCase.verifyEqual(Pall{j}.isDomBounded,isDB(j));
                end            
            end
        end
        function testDomainEntityCount(testCase)
            Pall=PLQVC.examples2();
            D=[1,0,0; 2,1,0; 2,1,0; 3,2,0; 3,2,0; 3,2,0; 3,2,0;3,2,1;8,7,2; 6,5,0];
            for i=1:length(Pall)
                %disp(i);
                p=Pall{i};
                testCase.verifyEqual([p.nv, p.ne, p.nf],D(i,:));            
            end
        end
        function testCreateDom1(testCase)
            Pall = PLQVC.examples2();
            p=Pall{9};               
            P = {[-7, -1, -6, 4];[-3, 5, 2, 7]};
            testCase.verifyEqual(p.dom.dim,2);
            testCase.verifyEqual(p.P,P);
            testCase.verifyEqual(p.dom.isConvex,true);
            testCase.verifyEqual(p.isDomBounded,false);           
        end
        function testOrderEdges(testCase)
            Pall = PLQVC.examples2();
            p=Pall{9};
            Pf = {[-7, -1, -6, 4],[-3, 5, 2, 7]};
            PDomain =[-4, 6, 1, -2, -5, 3];
            [P, isConvex] = orderEdges(p,1);            
            testCase.verifyEqual(P,Pf{1});
            testCase.verifyEqual(isConvex,true);
            [P, isConvex] = orderEdges(p,2);
            testCase.verifyEqual(P,Pf{2});
            testCase.verifyEqual(isConvex,true);
            [P, isConvex] = orderEdges(p,0);
            testCase.verifyEqual(P,PDomain);
            testCase.verifyEqual(isConvex,true);   
            p=Pall{1};%needle function
            [~, isConvex] = orderEdges(p,0);
            testCase.verifyEqual(isConvex,true);            
        end
        function testExample16(testCase)
           Pall=PLQVC.examples();
           nv=10;
           Dim =[2,2,2,2];
           isConvex =[true,true,true,true];
           isDomBounded = [false, true,false,false];
           entityCount = [6,5,4; nv,nv,1; 9,8,4; 12,12, 9];
           for i=1:length(Pall)-15
                p = Pall{i+15}; 
                testCase.verifyEqual(p.dom.dim,Dim(i));
                testCase.verifyEqual(p.dom.isConvex,isConvex(i));
                testCase.verifyEqual(p.isDomBounded,isDomBounded(i));
                testCase.verifyEqual([p.nv, p.ne, p.nf],entityCount(i,:));  
           end
        end
        function testNonconvex1(testCase)            
            V = [0,0;1,1;2,1;0,2;2,2];
            E = [1,2,1; 2,3,1; 1,4,1; 3,5,1; 4,5,1];
            F = [1,0  ; 1,0  ; 0,1  ; 1,0  ; 0,1  ];
            f = 0;
            p = PLQVC(V,E,f,F);                       
            testCase.verifyEqual(p.dom.dim,2);
            testCase.verifyEqual(p.dom.isConvex,false);
            %testCase.verifyEqual(p.dom.P,[-1,3,5,-4,-2]);
            testCase.verifyEqual(p.dom.P,[1,2, 4, -5, -3]);%remember domain edge order is inverted
            testCase.verifyEqual(p.isDomBounded,true);
            testCase.verifyEqual([p.nv, p.ne, p.nf],[5,5,1]);
            [P,isConvex] = p.orderEdges(1);
            testCase.verifyEqual(P,[-1,3,5,-4,-2]);
            testCase.verifyEqual(isConvex,false);
        end
        function testNonconvex2(testCase)            
            V = [0,0;1,1;2,1;0,2;2,2];
            E = [2,1,1; 2,3,1; 1,4,1; 3,5,1; 4,5,1];
            F = [0,1  ; 1,0  ; 0,1  ; 1,0  ; 0,1  ];
            f = 0;
            p = PLQVC(V,E,f,F);                       
            testCase.verifyEqual(p.dom.dim,2);
            testCase.verifyEqual(p.dom.isConvex,false);
            %testCase.verifyEqual(p.dom.P,[-1,3,5,-4,-2]);
            testCase.verifyEqual(p.dom.P,[-1,2, 4, -5, -3]);%remember domain edge order is inverted
            testCase.verifyEqual(p.isDomBounded,true);
            testCase.verifyEqual([p.nv, p.ne, p.nf],[5,5,1]);
            [P,isConvex] = p.orderEdges(1);
            testCase.verifyEqual(P,[1,3,5,-4,-2]);
            testCase.verifyEqual(isConvex,false);
        end
        function testNonconvex3(testCase)            
            V = [0,0;1,1;2,1;0,2;2,2];
            E = [2,1,1; 3,2,1; 1,4,1; 3,5,1; 4,5,1];
            F = [0,1  ; 0,1  ; 0,1  ; 1,0  ; 0,1  ];
            f = 0;
            p = PLQVC(V,E,f,F);                       
            testCase.verifyEqual(p.dom.dim,2);
            testCase.verifyEqual(p.dom.isConvex,false);
            %testCase.verifyEqual(p.dom.P,[-1,3,5,-4,-2]);
            testCase.verifyEqual(p.dom.P,[-1,-2, 4, -5, -3]);%remember domain edge order is inverted
            testCase.verifyEqual(p.isDomBounded,true);
            testCase.verifyEqual([p.nv, p.ne, p.nf],[5,5,1]);
            [P,isConvex] = p.orderEdges(1);
            testCase.verifyEqual(P,[1,3,5,-4,2]);
            testCase.verifyEqual(isConvex,false);
        end
        function testNonconvex4(testCase)            
            V = [0,0;1,1;2,1;0,2;2,2];
            E = [1,2,1; 3,2,1; 1,4,1; 3,5,1; 4,5,1];
            F = [1,0  ; 0,1  ; 0,1  ; 1,0  ; 0,1  ];
            f = 0;
            p = PLQVC(V,E,f,F);                      
            testCase.verifyEqual(p.dom.dim,2);
            testCase.verifyEqual(p.dom.isConvex,false);
            %testCase.verifyEqual(p.dom.P,[-1,3,5,-4,-2]);
            testCase.verifyEqual(p.dom.P,[1,-2, 4, -5, -3]);%remember domain edge order is inverted
            testCase.verifyEqual(p.isDomBounded,true);
            testCase.verifyEqual([p.nv, p.ne, p.nf],[5,5,1]);
            [P,isConvex] = p.orderEdges(1);
            testCase.verifyEqual(P,[-1,3,5,-4,2]);
            testCase.verifyEqual(isConvex,false);
        end        
        function testNonconvex5(testCase)            
            V = [-2,1;-1,0;0,1;1,0;2,1];
            E = [2,1,0;2,3,1;3,4,1;4,5,0];
            F = [0,1  ; 1,0 ; 1,0 ; 1,0 ];
            f = 0;
            p = PLQVC(V,E,f,F);                 
            testCase.verifyEqual(p.dom.dim,2);
            testCase.verifyEqual(p.dom.isConvex,false);
            %testCase.verifyEqual(p.dom.P,[-1,3,5,-4,-2]);
            testCase.verifyEqual(p.dom.P,[-1,2, 3, 4]);%remember domain edge order is inverted
            testCase.verifyEqual(p.isDomBounded,false);
            testCase.verifyEqual([p.nv, p.ne, p.nf],[5,4,1]);
            [P,isConvex] = p.orderEdges(1);
            testCase.verifyEqual(P,[-4, -3, -2, 1]);
            testCase.verifyEqual(isConvex,false);
        end   
        function testLineEquation_isCubicConvexOnEdge(testCase)
            v1=[0, 0];v2=[1, 0];
            c = PLQVC.lineEquation(v1,v2);%y=0
            testCase.verifyEqual(rank([c;[0,1,0]]),1);%use rank since coefficients are not normalized
            isConvex = PLQVC.isCubicConvexOnEdge([6,0,0,0,2,0,0,1,0,1],v1,v2,true);%x^3+x^2+x+1 on segment
            testCase.verifyEqual(isConvex,true);
            isConvex = PLQVC.isCubicConvexOnEdge([6,0,0,0,2,0,0,1,0,1],v1,v2,false);%x^3+x^2+x+1 on ray
            testCase.verifyEqual(isConvex,true);
            isConvex = PLQVC.isCubicConvexOnEdge(-[6,0,0,0,2,0,0,1,0,1],v1,v2,true);%-(x^3+x^2+x+1) on segment
            testCase.verifyEqual(isConvex,false);       
            
            v1=[0, 0];v2=[0, 1];
            c = PLQVC.lineEquation(v1,v2);%x=0
            testCase.verifyEqual(rank([c;[1,0,0]]),1);%use rank since coefficients are not normalized
            isConvex = PLQVC.isCubicConvexOnEdge([0,0,0,6,0,0,2,1,0,1],v1,v2,true);%y^3+y^2+y+1 on segment
            testCase.verifyEqual(isConvex,true);
            isConvex = PLQVC.isCubicConvexOnEdge([0,0,0,6,0,0,2,0,1,1],v1,v2,false);%y^3+y^2+y+1 on ray
            testCase.verifyEqual(isConvex,true);
            isConvex = PLQVC.isCubicConvexOnEdge(-[0,0,0,6,0,0,2,0,1,1],v1,v2,true);%-(y^3+y^2+y+1) on segment
            testCase.verifyEqual(isConvex,false); 
            
            v1=[0, 0];v2=[1, 1];
            c = PLQVC.lineEquation(v1,v2);%y=x
            testCase.verifyEqual(rank([c;[1,-1,0]]),1);
            isConvex = PLQVC.isCubicConvexOnEdge([1,0,0,1,0,0,0,0,0,0],v1,v2,true);%x^3+y^3
            testCase.verifyEqual(isConvex,true);
            isConvex = PLQVC.isCubicConvexOnEdge([1,0,0,-1,0,0,0,0,0,0],v1,v2,true);%x^3-y^3
            testCase.verifyEqual(isConvex,true);            
            v1=v1-[1;1];v2=v2+0.5;
            c = PLQVC.lineEquation(v1,v2);
            testCase.verifyEqual(rank([c;[1,-1,0]]),1);           
            v1=[0, 1];v2=[1, 0];
            c = PLQVC.lineEquation(v1,v2);
            testCase.verifyEqual(rank([c;[-1,-1,1]]),1);
            v1=[0, 0];v2=[0, 1];
            c = PLQVC.lineEquation(v1,v2);
            testCase.verifyEqual(rank([c;[1,0,0]]),1);  
            v1=[-1, 0];v2=[0, 1];
            isConvex = PLQVC.isCubicConvexOnEdge([1,0,0,-1,0,0,0,0,0,0],v1,v2,true);%x^3-y^3
            testCase.verifyEqual(isConvex,false);    
        end
        function testIsFaceConvex(testCase)
            p=PLQVC.energy;
            isConvex = p.isFaceConvex(1);            
            testCase.verifyEqual(isConvex,true);  
            Pall=PLQVC.examplesNonconvex;
            p=Pall{1};
            isConvex = p.isFaceConvex(1);            
            testCase.verifyEqual(isConvex,false);  
            p=Pall{2};
            isConvex = p.isFaceConvex(1);            
            testCase.verifyEqual(isConvex,false);  
            p = PLQVC.cubic1;
            isConvex = p.isFaceConvex(1);            
            testCase.verifyEqual(isConvex,true);                          
        end
        function testIsEdgeConvex(testCase)
            Pall = PLQVC.examples2();
            %domain is dimension 1
            p = Pall{2};% 2 segment
            [isConvex, isContinuous] = isEdgeConvex(p,1);
            testCase.verifyEqual([isConvex, isContinuous],[true,true]);
            p = Pall{3};%3 ray
            [isConvex, isContinuous] = isEdgeConvex(p,1);
            testCase.verifyEqual([isConvex, isContinuous],[true,true]);   
            p = Pall{4};%4 line
            [isConvex, isContinuous] = isEdgeConvex(p,1);
            testCase.verifyEqual([isConvex, isContinuous],[true,true]);   
            [isConvex, isContinuous] = isEdgeConvex(p,2);
            testCase.verifyEqual([isConvex, isContinuous],[true,true]);   
            p = Pall{5};%5 slice boundary
            [isConvex, isContinuous] = isEdgeConvex(p,1);
            testCase.verifyEqual([isConvex, isContinuous],[true,true]);   
            [isConvex, isContinuous] = isEdgeConvex(p,2);
            testCase.verifyEqual([isConvex, isContinuous],[true,true]);   
            %domain is dimension 2
            p = Pall{8};%8 half-space
            [isConvex, isContinuous] = isEdgeConvex(p,1);
            testCase.verifyEqual([isConvex, isContinuous],[true,true]);   
            [isConvex, isContinuous] = isEdgeConvex(p,2);
            testCase.verifyEqual([isConvex, isContinuous],[true,true]);   
            p = PLQVC.oneNorm();
            [isConvex, isContinuous] = isEdgeConvex(p,1);
            testCase.verifyEqual([isConvex, isContinuous],[true,true]);  
            p = PLQVC.oneNormConjugate();
            [isConvex, isContinuous] = isEdgeConvex(p,1);
            testCase.verifyEqual([isConvex, isContinuous],[true,true]);
            p = PLQVC.energy();
            verifyError(testCase,@() p.isEdgeConvex(1),'isEdgeConvex:NoEdge');
            Pall=PLQVC.examples();
            p = Pall{6};
            [isConvex, isContinuous] = isEdgeConvex(p,1);
            testCase.verifyEqual([isConvex, isContinuous],[false,true]);
            p = Pall{8};
            [isConvex, isContinuous] = isEdgeConvex(p,1);
            testCase.verifyEqual([isConvex, isContinuous],[true,true]);
            p = Pall{11};
            [isConvex, isContinuous] = isEdgeConvex(p,1);
            testCase.verifyEqual([isConvex, isContinuous],[true,true]);
            p = Pall{12};
            [isConvex, isContinuous] = isEdgeConvex(p,1);
            testCase.verifyEqual([isConvex, isContinuous],[true,true]);
            p = PLQVC.cubic1();
            [isConvex, isContinuous] = isEdgeConvex(p,1);
            testCase.verifyEqual([isConvex, isContinuous],[true,true]);
        end
        function testIsEdgeContinuous(testCase)
            Pall = PLQVC.examples2();
            %domain is dimension 1
            p = Pall{2};% 2 segment
            b = isEdgeContinuous(p,1);
            testCase.verifyEqual(b,true);
            p = Pall{3};%3 ray
            b = isEdgeContinuous(p,1);
            testCase.verifyEqual(b,true);  
            p = Pall{4};%4 line
            b = isEdgeContinuous(p,1);
            testCase.verifyEqual(b,true); 
            b = isEdgeContinuous(p,2);
            testCase.verifyEqual(b,true);  
            p = Pall{5};%5 slice boundary
            b = isEdgeContinuous(p,1);
            testCase.verifyEqual(b,true);   
            b = isEdgeContinuous(p,2);
            testCase.verifyEqual(b,true);
            %domain is dimension 2
            p = Pall{8};%8 half-space
            b = isEdgeContinuous(p,1);
            testCase.verifyEqual(b,true);   
            b = isEdgeContinuous(p,2);
            testCase.verifyEqual(b,true);   
            p = PLQVC.oneNorm();
            b = isEdgeContinuous(p,1);
            testCase.verifyEqual(b,true);   
            Pall = PLQVC.examplesDiscontinuous();
            p=Pall{1};
            b = isEdgeContinuous(p,1);
            testCase.verifyEqual(b,false);   
            p=Pall{2};
            b = isEdgeContinuous(p,1);
            testCase.verifyEqual(b,false);
            p=Pall{3};
            b = isEdgeContinuous(p,1);
            testCase.verifyEqual(b,false);
            p=Pall{3};
            b = isEdgeContinuous(p,1);
            testCase.verifyEqual(b,false);
            b = isEdgeContinuous(p,2);
            testCase.verifyEqual(b,false);
            
        end
        function testIsQuadratic(testCase)
            p=PLQVC.oneNorm();
            testCase.verifyEqual(p.isQuadratic,true);
            p=PLQVC.oneNormConjugate();
            testCase.verifyEqual(p.isQuadratic,true);
            p=PLQVC.energy();
            testCase.verifyEqual(p.isQuadratic,true);
            p=PLQVC.cubic1();
            testCase.verifyEqual(p.isQuadratic,false);
        end
        function testIsPositiveSemidefinite(testCase)            
            testCase.verifyEqual(PLQVC.isPositiveSemidefinite(eye(2,2)),true);
            testCase.verifyEqual(PLQVC.isPositiveSemidefinite([-1,0;0,1]),false);
            testCase.verifyEqual(PLQVC.isPositiveSemidefinite([0 1; 1 0]),false);
            testCase.verifyEqual(PLQVC.isPositiveSemidefinite([1 0; 0 0]),true);
            H=cat(3,[0,0;0,0],[10,1;1,20],[-2,0;0,-4]);
            testCase.verifyEqual(PLQVC.isPositiveSemidefinite(H),false);
        end        
        function testEvalHessian(testCase)
            %single Hessian
            syms x y
            C = [x^3, x^2*y, x*y^2, y^3, x^2, x*y, y^2, x, y, 1];
            C = C .* [1/6, 1/2, 1/2, 1/6, 1/2, 1, 1/2, 1, 1, 1];
            c = [ 0 ,   0  ,   0  ,  0 ,  1 ,  0 ,  1 , 0, 0, 0];
            p = c*C.';
            H = double(hessian(p,[x,y]));
            Hc = PLQVC.evalHessian(c);
            testCase.verifyEqual(Hc,H);
            c = [ 1 ,   zeros(1,9)];            
            verifyError(testCase,@()PLQVC.evalHessian(c),'evalHessian:xNotProvidedForCubic');
            %multiple Hessians
            v = [0,0;1,0;0,1;1,1;-1,-1;-2,3;3,5];
            Hc = PLQVC.evalHessian(c,v);
            p = c*C.';
            H = hessian(p,[x,y]);
            f = @(v)double(subs(H,[x,y],[v(1),v(2)]));   
            Hs = zeros(2,2,size(v,1));
            for i=1:size(v,1)
                Hs(:,:,i) = f(v(i,:));
            end
            testCase.verifyEqual(Hc,Hs);     
        end
        function testEvalHessian2(testCase)
            %single Hessian
            syms x y
            C = [x^3, x^2*y, x*y^2, y^3, x^2, x*y, y^2, x, y, 1];
            C = C .* [1/6, 1/2, 1/2, 1/6, 1/2, 1, 1/2, 1, 1, 1];
            c = [ 0 ,   0  ,   0  ,  0 ,  1 ,  0 ,  1 , 0, 0, 0];
            p = c*C.';
            H = double(hessian(p,[x,y]));
            g = gradient(p,[x,y]);
            v = [0,0;1,0;0,1;1,1;-1,-1;-2,3;3,5];%need points to evaluate gradient
            [Hc,gc] = PLQVC.evalHessian(c,v);
            testCase.verifyEqual(Hc,repmat(H,1,1,size(v,1)));
            fg = @(v)double(subs(g,[x,y],[v(1),v(2)]));
            gs = zeros(2,size(v,1));
            for i=1:size(v,1)
                gs(:,i) = fg(v(i,:));
            end
            testCase.verifyEqual(gc,gs);
            c = [ 1 ,   zeros(1,9)];            
            verifyError(testCase,@()PLQVC.evalHessian(c),'evalHessian:xNotProvidedForCubic');
            %multiple Hessians
            v = [0,0;1,0;0,1;1,1;-1,-1;-2,3;3,5];
            [Hc,gc] = PLQVC.evalHessian(c,v);
            p = c*C.';
            H = hessian(p,[x,y]);
            g = gradient(p,[x,y]);
            fH = @(v)double(subs(H,[x,y],[v(1),v(2)]));
            fg = @(v)double(subs(g,[x,y],[v(1),v(2)]));
            Hs = zeros(2,2,size(v,1));gs = zeros(2,size(v,1));
            for i=1:size(v,1)
                Hs(:,:,i) = fH(v(i,:));
                gs(:,i) = fg(v(i,:));
            end
            testCase.verifyEqual(Hc,Hs);
            testCase.verifyEqual(gc,gs);
            c = (1:10);
            v = [0,0;1,0;0,1;1,1;-1,-1;-2,3;3,5];
            [Hc,gc] = PLQVC.evalHessian(c,v);
            p = c*C.';
            H = hessian(p,[x,y]);
            g = gradient(p,[x,y]);
            fH = @(v)double(subs(H,[x,y],[v(1),v(2)]));
            fg = @(v)double(subs(g,[x,y],[v(1),v(2)]));
            Hs = zeros(2,2,size(v,1));gs = zeros(2,size(v,1));
            for i=1:size(v,1)
                Hs(:,:,i) = fH(v(i,:));
                gs(:,i) = fg(v(i,:));
            end
            testCase.verifyEqual(Hc,Hs);
            testCase.verifyEqual(gc,gs);
        end
        
        function testisCollinear(testCase)
           P=[1, 1; 3.5, 3.5; -7.2, -7.2];
           b = PLQVC.isCollinear(P,1E-10);
           testCase.verifyTrue(b,true);
           P=[0, 1;0, 3.5;0 -7.2];
           b = PLQVC.isCollinear(P,sqrt(eps));
           testCase.verifyTrue(b,true);
           P=[0 1; 14 5; 8 4];
           b = PLQVC.isCollinear(P,1E-10);
           testCase.verifyFalse(b,false);
           V1=[0,0];V2=[1,0];
           X = [0,0;0.5,0;1,0;-1E-16,0;0.5,1E-12;1,1;1+1E-4,0];
           %check whether each point in X(:,1:2) belong to [V1,V2]
           
        end
        function testBelongToEdge(testCase)
            verifyError(testCase,@()PLQVC.belongToEdge([0,0],[1, 2;3 4],[1,1]),'belongToEdge:dimensionMismatch');
            V1=[0,0];V2=[1,0];
            X = [0,0;0.5,0;1,0;-1E-16,0;0.5,1E-12;1,1;1+1E-4,0];
            b = PLQVC.belongToEdge(V1,V2,X);
            testCase.verifyEqual(b,[true;true;true;true;true;false;false]);
            V1=[0,0];V2=[0,1];
            X = X(:,[2,1]);%swap columns
            b = PLQVC.belongToEdge(V1,V2,X);
            testCase.verifyEqual(b,[true;true;true;true;true;false;false]);
            V1=[0,0];V2=[1,1];
            X = [0,0;0.5,0.5;1,1;-1E-16,0;0.5,0.5+1E-12;1,2;1+1E-4,1;-1,0];
            b = PLQVC.belongToEdge(V1,V2,X);
            testCase.verifyEqual(b,[true;true;true;true;true;false;false;false]);
            
        end
        function testEvalPoly_Vectorized(testCase)
           c = 0;X=[1 2; 3 4];
           z = PLQVC.evalPoly(c,X);
           testCase.verifyEqual(z,[0;0]);
           c = [0;0];
           verifyError(testCase,@()PLQVC.evalPoly(c,X),'evalPoly:notVectorizedInC');
        end
        function testEdgeChainFunction(testCase)
           Pall=PLQVC.examples2();
           p=Pall{10};
           testCase.verifyEqual(p.eval([0.5,0]),0.5^2);
           X=[0,0;0.5,0;1,0;1,1;-1,-1];
           fexp=[0;0.5^2;1;nan;inf];
           testCase.verifyEqual(p.eval(X),fexp);
           
        end
%       examples of plotting; not automated testing so commented out        
%         function mytestPlotDomainGRaph(testCase)
%             V = [0 2; -0.5 2.5; 3 1.25; 2 1; 5 1.75; 5.5 3.25; 6 2.75; 2.5 1.75];
%             E = [4 3 1; 5 3 1; 7 6 0; 1 2 0; 7 5 1; 1 4 1; 3 8 0];
%             F = [1 0; 2 0; 2 0; 0 1; 0 2; 1 0; 1 2];
%             f = [0;0];            
%             p = PLQVC(V,E,f,F);%3 rays
%             p.plotDomainGraph            
%             p = PLQVC.oneNormConjugate;%no ray
%             p.plotDomainGraph;
%         end
    end
end 