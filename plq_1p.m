classdef plq_1p
    properties
        f;
        d;
        envelope = functionNDomain.empty();
        conjfia = [];
        conjugates = functionNDomain.empty();
        maxConjugate = functionNDomain.empty();
        % For printouts 
        lCE = false
        lConj = false
        lMConj = false
        lPrintEta = false %true
    end




    methods % creation & print
         function obj = plq_1p(d,f)
            % put checks for type of f and d
            if nargin > 0
              obj.f = f;
              obj.d = d;
            end 
         end

         function print(obj)
         
           disp("Domain")
           obj.d.print
           fprintf("\n")
           disp("Function")
          
           obj.f.f = simplifyFraction(obj.f.f);
           obj.f.print
           fprintf("\n\n\n")
           
           disp("Convex Envelope")
           for j=1: size(obj.envelope,2)
             if obj.lCE  
                 disp('Function')  
                
                 obj.envelope(j).f.f = simplifyFraction(obj.envelope(j).f.f);
                 obj.envelope(j).f.print
                 disp('Domain')
                 obj.envelope(j).d.print
             end
             if obj.lConj
                 if (size(obj.conjfia,1) > 0)
                     obj.conjugates(obj.conjfia(j):obj.conjfia(j+1)-1).printL
                 end
             end
           end
           fprintf("\n\n\n\n\n")
           if ~obj.lMConj
               return
           end
           disp("Maximum conjugate")
           obj.maxConjugate.printL
         end

         function printLatex(obj)
         
           disp("Domain")
           %obj.d.printLatex
           obj.d.polygon.printLatex
           fprintf("\n")
           disp("Function")
           
           obj.f.printLatex
           fprintf("\n\n\n")
              
           disp("Convex Envelope")
           disp(" ")
           for j=1:size(obj.envelope,2) 
             disp('Function')  
             obj.envelope(j).f.f = simplifyFraction(obj.envelope(j).f.f);
             obj.envelope(j).f.printLatex
             disp('Domain')
             disp(" ")
             obj.envelope(j).d.printLatex
             if (size(obj.conjfia,1) > 0)
                 obj.conjugates(obj.conjfia(j):obj.conjfia(j+1)-1).printLLatex
             end
           end
           fprintf("\n\n\n\n\n")
           disp("Maximum conjugate")
           obj.maxConjugate.printLLatex
         end


         
         function Mprint(obj)
           fprintf("display(inequal(("); 
           obj.d.polygon.printMaple
           
           fprintf("),x=-5..5,y=-5..5,color=[red,blue,yellow,green],nolines)) \n")
           if obj.lCE  
             
           disp("Convex Envelope")
           fprintf("display(inequal((");
           for j=1:size(obj.envelope,2)-1 
             obj.envelope(j).d.printMaple
             fprintf(",");
           end
           obj.envelope(size(obj.envelope,2)).d.printMaple
           fprintf("),x=-5..5,y=-5..5,color=[red,blue,yellow,green],nolines)) \n")
           end
           if obj.lConj
             
           

           for j=1:size(obj.envelope,2) 
           if (size(obj.conjfia,1) > 0)
                  obj.conjugates(obj.conjfia(j):obj.conjfia(j+1)-1).printM
           end
           end
           end
           if ~ obj.lMConj
             return
           end
           disp('max in piece')
           obj.maxConjugate.printM;
         end

         function plotMaxConjugateDomain(obj)
             figure;
             colors = ['b', 'r', 'g', 'm', 'c', 'y'];
             n = 0
             f = obj.maxConjugate (1).f
             c = colors(mod(n,6)+1)
             for i =1:size(obj.maxConjugate,2)
                if (f.f ~= obj.maxConjugate (i).f.f)
                  n = n + 1
                  c = colors(mod(n,6)+1)
                  f = obj.maxConjugate (i).f
                end
                obj.maxConjugate (i).d.plot;
                textR = "R"+num2str(i);
                textR="";
                obj.maxConjugate (i).d.plotRegionC(textR,c);
             end
         end

         function plotDomain(obj)
             figure;
             obj.d.plot;
         end
    end

    methods
        function ps = triangulate (obj, ps)
            
            d = obj.d;
            vars = d.polygon.vars;
            if d.polygon.nv == 3
                ps = [ps,obj];
                return
            end
            vx = d.polygon.vx;
            vy = d.polygon.vy;

            % (mvx,mvy) : vertex in up left corner
            mvy = max(vy);
            ind = [];
            for i = 1:d.polygon.nv
                if vy(i) ~= mvy
                    ind = [ind,i];
                end
            end
            vx(ind) = [];
            mvx = min(vx);
            vx = d.polygon.vx;
            for i = 1:d.polygon.nv
                if vx(i) == mvx & vy(i) == mvy
                    break;
                end
            end

            start = i;
            

            for i = start+1:d.polygon.nv-1
                triangle = [start,i,i+1];
                 t(:,1) = vx(triangle);
                 t(:,2) = vy(triangle);
                 ps = [ps,plq_1p(domain(t,vars(1),vars(2)),obj.f)];
                 
                 
            end
            if start ~= 1 & start ~= d.polygon.nv
                triangle = [start,d.polygon.nv,1];
                 t(:,1) = vx(triangle);
                 t(:,2) = vy(triangle);
                 ps = [ps,plq_1p(domain(t,vars(1),vars(2)),obj.f)];
            end
            for i = 1: start-2
                triangle = [start,i,i+1];
                 t(:,1) = vx(triangle);
                 t(:,2) = vy(triangle);
                 ps = [ps,plq_1p(domain(t,vars(1),vars(2)),obj.f)];
            end
           
        end
    end
    
    methods % convex envelope
        function obj = convexEnvelope(obj)
          vars = obj.f.getVars;
          x = vars(1);
          y = vars(2);
          obj = convexEnvelope1 (obj,x,y);
        end

        function obj = convexEnvelope1 (obj,x,y)
            % a=sym('a');
            % b=sym('b');
              
            % etaV : eta functions corresponding to set obj.d.V
            % etaE : eta functions corresponding to set obj.d.E
            % etaR : domain of etaE - stored as etaR(i,1:3) : [function,lb,ub] => lb <= function <= ub 
            
            % [etaV, etaE, etaR] =  getEtaFunctions (obj,x,y);
            % 
            % if obj.lPrintEta
            %   obj.d.V
            %   obj.d.E
            %   disp('etaV')
            %   etaV.printL
            %   disp('etaE')
            %   etaE.printL
            %   disp('etaR')
            %   etaR.printL
            % end

            nCE = obj.d.nE;  %size(etaE,1)
     
            
            if nCE == 0
              x1 = obj.d.polygon.vx(1);  
              x2 = obj.d.polygon.vx(2);  
              x3 = obj.d.polygon.vx(3);  
              y1 = obj.d.polygon.vy(1);  
              y2 = obj.d.polygon.vy(2);  
              y3 = obj.d.polygon.vy(3);  
              a = ((x1*y1*y2 - x1*y1*y3 - x2*y1*y2 + x2*y2*y3 + x3*y1*y3 - x3*y2*y3))/((x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2));
              b = ((x1*x2*y2 - x1*x2*y1 + x1*x3*y1 - x1*x3*y3 - x2*x3*y2 + x2*x3*y3))/((x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2));
              c = ((x1*x2*y1*y3 - x1*x3*y1*y2 - x1*x2*y2*y3 + x2*x3*y1*y2 + x1*x3*y2*y3 - x2*x3*y1*y3))/((x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2));
              f = symbolicFunction(a*x+b*y+c);
              obj.envelope = [obj.envelope,functionNDomain(f, obj.d.polygon)];

              obj.lCE = true;
              return
            end
            if nCE == 1
              
              x2 = obj.d.polygon.vx(2);  
              x3 = obj.d.polygon.vx(3);  
              y1 = obj.d.polygon.vy(1);  
              y2 = obj.d.polygon.vy(2);  
              y3 = obj.d.polygon.vy(3);  

              m =  obj.d.mE(1);
              q =  obj.d.cE(1);
              x1 = obj.d.polygon.vx(obj.d.V(1))  ;
              y1 = obj.d.polygon.vy(obj.d.V(1))  ;

              a = -m*y1;
              b = q;
              c = x1;
              d = -q*y1 + m*x1*y1;
              e = -q*x1 - x1*y1;
              f0 = q*x1*y1;
              g = -m;
              h = 1;
              k = -y1 + m*x1;
              f = symbolicFunction(a*x^2+b*x*y+c*y^2+d*x+e*y+f0,g*x+h*y+k);
              obj.envelope = [obj.envelope,functionNDomain(f, obj.d.polygon)];

              obj.lCE = true;
              return
            end

            if nCE == 2
                m_h = sym('m_h');
                 m_w = sym('m_w');
                 q_h = sym('q_h');
                 q_w = sym('q_w');
                % m_h =  obj.d.mE(1);
                % q_h =  obj.d.cE(1);
                % m_w =  obj.d.mE(2);
                % q_w =  obj.d.cE(2);
                
                a =  (m_h*m_w)/(m_h + m_w + 2*sqrt(m_h*m_w));
                b =  (2*sqrt(m_h*m_w))/(m_h + m_w + 2*sqrt(m_h*m_w));
                c =  (1)/(m_h + m_w + 2*sqrt(m_h*m_w));
                d =  ((m_h*q_w + m_w*q_h))/(m_h + m_w + 2*sqrt(m_h*m_w));
                e =  ((- (q_h + q_w)))/(m_h + m_w + 2*sqrt(m_h*m_w));
                f0 =  (q_h*q_w)/(m_h + m_w + 2*sqrt(m_h*m_w));

                
                f = symbolicFunction(a*x^2+b*x*y+c*y^2+d*x+e*y+f0);
                f = f.subsF([m_h,m_w],obj.d.mE(1:2));
                f = f.subsF([q_h,q_w],obj.d.cE(1:2));
                
                obj.envelope = [obj.envelope,functionNDomain(f, obj.d.polygon)];
                
                obj.lCE = true;
                return
            end
           
             
        end
           % disp('out')

         

        function [etaV, etaE, etaR] = getEtaFunctions (obj,x,y)

            a=sym('a');
            b=sym('b');
            eta = obj.f - symbolicFunction(a*x+b*y);
            
            % Eta for Edges
            for i = 1:obj.d.nE
                xv1 = obj.d.polygon.vx(obj.d.E(i,1));
                yv1 = obj.d.polygon.vy(obj.d.E(i,1));
                
                xv2 = obj.d.polygon.vx(obj.d.E(i,2));
                yv2 = obj.d.polygon.vy(obj.d.E(i,2));
                
                edgey = obj.d.mE(i) * x + obj.d.cE(i);
                etaT = eta.subsVarsPartial([y],[edgey]);
                df = etaT.dfdx(x);

                f0 = obj.f.subsVarsPartial([y],[edgey]);
                df0 = f0.dfdx(x);
                df1 = df0.subsF([x],[xv1]);
                df2 = df0.subsF([x],[xv2]);
                if df1 < df2
                    etaE(i,1) = eta.subsVarsPartial([x,y],[xv1,yv1]);
                    etaE(i,3) = eta.subsVarsPartial([x,y],[xv2,yv2]);
                    etaR(i,2) = df1; 
                    etaR(i,3) = df2; 
                else
                    etaE(i,3) = eta.subsVarsPartial([x,y],[xv1,yv1]);
                    etaE(i,1) = eta.subsVarsPartial([x,y],[xv2,yv2]);
                    etaR(i,3) = df1; 
                    etaR(i,2) = df2; 
                end
                etaR(i,1) = symbolicFunction(a+obj.d.mE(i)*b);
                etaE(i,2) =  symbolicFunction((-(a+obj.d.mE(i)*b-obj.d.cE(i))^2/(4*obj.d.mE(i)))-b*obj.d.cE(i))
               
                
                xp = df.solve(x);
                etaE(i,2) = etaT.subsVarsPartial([x],[xp]);
                etaR(i,1) = symbolicFunction(a+obj.d.mE(i)*b);
                
            end

            % Eta for Vertices
            for i = 1:obj.d.nV
                xv = obj.d.polygon.vx(obj.d.V(i));
                yv = obj.d.polygon.vy(obj.d.V(i));
                etaV(i) = eta.subsVarsPartial([x,y],[xv,yv]);
            end
            if obj.d.nV == 0
                etaV=symbolicFunction.empty();
            end
            if obj.d.nE == 0
                etaE=symbolicFunction.empty();
                etaR=symbolicFunction.empty();
            end
            
        end

    end

    methods % conjugate
        function obj = conjugate (obj)
            obj.conjfia(1) = 1;  
            for i=1:max(1 , size(obj.envelope,2)) %For triangles where convex envelope is not computed
              obj = obj.conjugateFunction(i);
              obj.conjfia(i+1) = size(obj.conjugates,2)+1;
            end
        end

        function obj = conjugateFunction (obj,i)
            
            nCE = obj.d.nE;
            vars = obj.f.getVars;
            s1 = sym('s_1');
            s2 = sym('s_2');
            dualVars = [s1,s2];
            if nCE == 0

                
                NCV = obj.envelope(i).d.getNormalConeVertex(s1, s2);
                [subdV,undV] =  obj.envelope(i).getSubdiffVertexT1 (NCV, dualVars);
                
                x1 = obj.d.polygon.vx(1);  
                x2 = obj.d.polygon.vx(2);  
                x3 = obj.d.polygon.vx(3);  
                y1 = obj.d.polygon.vy(1);  
                y2 = obj.d.polygon.vy(2);  
                y3 = obj.d.polygon.vy(3);  
                a = ((x1*y1*y2 - x1*y1*y3 - x2*y1*y2 + x2*y2*y3 + x3*y1*y3 - x3*y2*y3))/((x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2));
                b = ((x1*x2*y2 - x1*x2*y1 + x1*x3*y1 - x1*x3*y3 - x2*x3*y2 + x2*x3*y3))/((x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2));
                c = ((x1*x2*y1*y3 - x1*x3*y1*y2 - x1*x2*y2*y3 + x2*x3*y1*y2 + x1*x3*y2*y3 - x2*x3*y1*y3))/((x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2));
                for j = 1:obj.envelope(i).d.nv
                  x1 = obj.d.polygon.vx(j);   % same as obj.envelope(i).d.vx(1) for triangles
                  y1 = obj.d.polygon.vy(j);  
                    
                  conjf = symbolicFunction(s1 * x1 + s2 * y1  - (a*x1+b*y1+c));
                  conjd = region(subdV(j,:), dualVars);
                  obj.conjugates = [obj.conjugates,functionNDomain([conjf],conjd)];
                end
                    
              obj.lConj = true;
              return
            end
            
            if nCE == 1
                
              m =  obj.d.mE(1);
              q =  obj.d.cE(1);
              x1 = obj.d.polygon.vx(obj.d.V(1))  ;
              y1 = obj.d.polygon.vy(obj.d.V(1))  ;

              a = -1;
              b = -2*m;
              d = 2*q+4*m*x1;
              c = -m^2;
              e = -(2*m*q - 4*m*y1);
              f = -(q^2 + 4*m*x1*y1);
             

             
              crs = a*s1^2+b*s1*s2 + c* s2^2 + d*s1 + e*s2 + f;   % nonconvex part


              NCV = obj.envelope(i).d.getNormalConeVertex(s1, s2);
              [subdV,undV] =  obj.envelope(i).getSubdiffVertexT1 (NCV, dualVars);
              subdV = obj.envelope(i).getSubDiffVertexSpT1(subdV, undV, -crs);


              NCE = obj.envelope(i).d.getNormalConeEdge(s1, s2);
              [subdE,unR] = obj.envelope(i).getSubdiffVertexT2 (NCE, dualVars);
              
              %obj.envelope(i).d.print

              edgeNo = obj.envelope(i).d.getEdgeNos(vars);
             
              [subdE, unR, crs] = obj.envelope(i).getSubDiffEdgeT1(subdE, edgeNo, undV, crs, dualVars);
              
              
              %expr = obj.envelope(i).conjugateExprVerticesT1 (dualVars, undV )

                x1 = obj.d.polygon.vx(1);  
                x2 = obj.d.polygon.vx(2);  
                x3 = obj.d.polygon.vx(3);  
                y1 = obj.d.polygon.vy(1);  
                y2 = obj.d.polygon.vy(2);  
                y3 = obj.d.polygon.vy(3);  
                a = ((x1*y1*y2 - x1*y1*y3 - x2*y1*y2 + x2*y2*y3 + x3*y1*y3 - x3*y2*y3))/((x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2));
                b = ((x1*x2*y2 - x1*x2*y1 + x1*x3*y1 - x1*x3*y3 - x2*x3*y2 + x2*x3*y3))/((x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2));
                c = ((x1*x2*y1*y3 - x1*x3*y1*y2 - x1*x2*y2*y3 + x2*x3*y1*y2 + x1*x3*y2*y3 - x2*x3*y1*y3))/((x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2));
                for j = 1:obj.envelope(i).d.nv
                  if undV(j)
                      if j == obj.envelope(i).d.nv
                        e0 = 1;
                      else    
                        e0 = j+1;
                      end  
                      x1 = obj.d.polygon.vx(j);   % same as obj.envelope(i).d.vx(1) for triangles
                      y1 = obj.d.polygon.vy(j);  
                  
                      conjf = symbolicFunction(s1 * x1 + s2 * y1  - (a*x1+b*y1+c));
                      conjd = region([subdE(e0,1:2),-subdE(e0,3)], dualVars);
                      obj.conjugates = [obj.conjugates,functionNDomain([conjf],conjd)];
                    %%%%%%%%%%%%%
                      r = subdV(e0,:);
                      r(2) = -r(2);
                       conjd = region(r, dualVars);
                       obj.conjugates = [obj.conjugates,functionNDomain([conjf],conjd)];
                       %%%%%%%%%%%%%
                       if j == 1
                        e0 = obj.envelope(i).d.nv;
                       else    
                        e0 = j-1;
                       end  
                       r = subdV(e0,:);
                       r(1) = -r(1);
                       conjd = region(r, dualVars);
                       obj.conjugates = [obj.conjugates,functionNDomain([conjf],conjd)];
                  else



                      x1 = obj.d.polygon.vx(j);   % same as obj.envelope(i).d.vx(1) for triangles
                      y1 = obj.d.polygon.vy(j);  
                        
                      conjf = symbolicFunction(s1 * x1 + s2 * y1  - (a*x1+b*y1+c));
                      conjd = region(subdV(j,:), dualVars);
                      obj.conjugates = [obj.conjugates,functionNDomain([conjf],conjd)];
                  end
                end
                


               for j = 1:obj.envelope(i).d.nv

                  if unR(j)
                       a = 1/(4*m);
                       b = 1/2;
                       c = m/4;
                       d = -q/(2*m);
                       e =  q/2;
                       f = q^2/(4*m);
                       conjf = symbolicFunction(a*s1^2+b*s1*s2+c*s2^2+d*s1+e*s2+f);
                       % obj.envelope(1).f
                       % edge = vars(2) - m*vars(1) -q
                       % edge
                       % conjugateExpr (edge, obj.envelope(1).f.f,vars(1),vars(2))
                       % pause
                       conjd = region(subdE(j,:), dualVars);
                       obj.conjugates = [obj.conjugates,functionNDomain([conjf],conjd)];
                  
                  end
               end
               obj.lConj = true;
               return
            end
            if nCE == 2
                
              NCV = obj.envelope(i).d.getNormalConeVertex(s1, s2);
              [subdV,undV] =  obj.envelope(i).getSubdiffVertexT1 (NCV, dualVars);
              

               [coef,terms] = coeffs(obj.envelope(i).f.f,obj.envelope(i).d.vars);
               x = obj.d.polygon.vars(1);
               y = obj.d.polygon.vars(2);
               for j = 1:obj.envelope(i).d.nv
                 x1 = obj.d.polygon.vx(j); 
                 y1 = obj.d.polygon.vy(j);  
                 % change this for direct computation
                 conjf = symbolicFunction(s1 * x1 + s2 * y1)  - obj.envelope(i).f.subsF([x,y],[x1,y1]);
                 %subdV(j,:)
                 conjd = region(subdV(j,:), dualVars);
                 %conjd.print
                 if isempty(conjd)
                     continue
                 end
                 obj.conjugates = [obj.conjugates,functionNDomain([conjf],conjd)];
                end
                %pause
              for j = 1:size(terms,2)
                  if isAlways (terms(j) == x^2)
                      a = coef(j);
                  end
                  if isAlways (terms(j) == x*y)
                      b = coef(j);
                  end
                  if isAlways (terms(j) == y^2)
                      c = coef(j);
                  end
                  if isAlways (terms(j) == x)
                      d = coef(j);
                  end
                  if isAlways (terms(j) == y)
                      e = coef(j);
                  end
                  if isAlways (terms(j) == 1)
                      f = coef(j);
                  end
              end
              %obj.envelope(i).d.ineqs(1)
              grad = b*s1 - 2*a*s2 - b*d + 2*a*e;

              NCE = obj.envelope(i).d.getNormalConeEdge(s1, s2);
              [subdE,unR] = obj.envelope(i).getSubdiffVertexT2 (NCE, dualVars); 
              
              % check signs and then put gradient ineq properly
              % simplifyFraction(grad)
              % simplifyFraction(grad)
              % simplifyFraction(subdE(3,1))
              % simplifyFraction(subdE(3,2))
              % pause
              % obj.envelope(i).d.ineqs(1)
              % obj.envelope(i).d.ineqs(2)
              % obj.envelope(i).d.ineqs(3)
              [c0,t0] = (coeffs(simplifyFraction(grad),[s1,s2]));
              %double(c0)
              %pause
             subdE(1,3) = grad;
             subdE(2,3) = -grad;
             subdE(3,3) = grad;
 %            for j = 1:obj.envelope(i).d.nv
 % %                 j
 %  %                obj.envelope(i).d.ineqs(j)
 %                  d1(j) = 2*a*obj.envelope(i).d.vx(j) + b* obj.envelope(i).d.vy(j) +d;
 %                  d2(j) = b*obj.envelope(i).d.vx(j) + 2*c* obj.envelope(i).d.vy(j) +e;
 % 
 % 
 % 
 %            end
 %            double(d1)
 %            double(d2)
              for j = 1:obj.envelope(i).d.nv
 %                 j
                    coeffs0 = obj.envelope(i).d.ineqs(j).getLinearCoeffs(obj.envelope(i).d.vars);

                    % if coeffs0(2)/coeffs0(1) > 0
                    %     subdE(j,3) = -grad;
                    % else
                    %     subdE(j,3) = grad;
                    % end

                    % if j < obj.envelope(i).d.nv
                    %    double(subs(grad,[s1,s2],[(d1(j)+d1(j+1))/2,(d2(j)+d2(j+1))/2]))
                    % else
                    %    double(subs(grad,[s1,s2],[(d1(j)+d1(1))/2,(d2(j)+d2(1))/2]))
                    % end
                    % pause
  %                obj.envelope(i).d.ineqs(j)
                  
                  c0 = obj.envelope(i).d.ineqs(j).getLinearCoeffs ([x,y]);

                  m = -c0(1)/c0(2);   % put check for zero d
                  q = -c0(3)/c0(2);
% conj = (mh^2*y^2 + 2*mh*qh*y + 2*mh*x*y + qh^2 - 2*qh*x + x^2)/(4*mh)


                  av = 1/(4*m);
                   bv = (2*m)/(4*m);
                   cv = (m^2)/(4*m);
                   dv = (-(2*q))/(4*m);
                   ev = ((2*m*q))/(4*m);
                   fv = (q^2)/(4*m);

                   if m > 0
                    conjf = symbolicFunction(av*s1^2+bv*s1*s2+cv*s2^2+dv*s1+ev*s2+fv);
                   else
                    conjf =  symbolicFunction( conjugateExpr(obj.envelope(i).d.ineqs(j).f,obj.envelope(i).f.f,x,y));
                    conjf = conjf.subsF([x,y],[s1,s2]);
                   end
                   %obj.envelope(i).d.ineqs(j).f
                   %obj.envelope(i).f.f
                   %conjf =  symbolicFunction( conjugateExpr(obj.envelope(i).d.ineqs(j).f,obj.envelope(i).f.f,x,y))
                   %pause
                   %conjf = conjf.subsF([x,y],[s1,s2]);
                   %pause
                   conjd = region(subdE(j,:), dualVars);
                   
                   conjd = conjd.simplifyUnboundedRegion;
                  % conjd.print
                  % pause
                    if isempty(conjd)
                     continue
                 end
                   obj.conjugates = [obj.conjugates,functionNDomain([conjf],conjd)];
              end
              %pause
            end  

               obj.lConj = true;
              return
              
            end
        end

    

    methods % max

        function obj = maximumConjugate(obj)
          for k = obj.conjfia(1):obj.conjfia(2)-1 
             obj.maxConjugate(k) = obj.conjugates(k);
          end
          for i = 2:size(obj.envelope,2) %size(obj.conjfia,2)-1
              obj.maxConjugate = obj.maxConjugate * obj.conjugates(obj.conjfia(i):obj.conjfia(i+1)-1);
              %obj.maxConjugate.printM2
              %obj.maxConjugate.printL
             obj.maxConjugate = obj.maxConjugate.maximumP(true);
          end
          obj.lMConj = true;
        end
    end


end