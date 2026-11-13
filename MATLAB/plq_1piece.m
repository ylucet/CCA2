classdef plq_1piece
    properties
        f;
        d;
        envelope = functionNDomain.empty();
        conjfia = [];
        conjugates = functionNDomain.empty();
        maxConjugate = functionNDomain.empty();
    end

% 36 methods
    methods % creation & print
         function obj = plq_1piece(d,f)
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
          % size(obj.envf)
           for j=1:size(obj.envelope,2) 
             disp('Function')  
             obj.envelope(j).f.f = simplifyFraction(obj.envelope(j).f.f);
             obj.envelope(j).f.print
             % disp("Expr")
             % obj.envelopeExpr(j).print
             disp('Domain')
             obj.envelope(j).d.print
           
             %disp("Conjugate Expr")
             %obj.conjf.printL(obj.conjfia(j),obj.conjfia(j+1)-1)
             if (size(obj.conjfia,1) > 0)
                 obj.conjugates(obj.conjfia(j):obj.conjfia(j+1)-1).printL
             for k = obj.conjfia(j):obj.conjfia(j+1)-1
               %disp("Conjugate Expr")
               %obj.conjugates(k).print
               %obj.conjf(k).print
               %disp('Conjugate Domain')
               %obj.conjd(k).print
             end
             end
           end
           fprintf("\n\n\n\n\n")
           disp("Maximum conjugate")
           obj.maxConjugate.printL

           % size(obj.maxf,2)
           % for i = 1:size(obj.maxf,1)
           %   disp(i)
           %   obj.maxf(i).print
           %   obj.maxd(i).print
           % end
           % %disp("Conjugate Expr")
            % obj.conjf.printL
            % disp('Conjugate Domain')
             %obj.conjd(j).print

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
             %disp("Expr")
             %obj.envExpr(j).printLatex
             disp('Domain')
             disp(" ")
             obj.envelope(j).d.printLatex
           
             %disp("Conjugate Expr")
             %obj.conjf.printL(obj.conjfia(j),obj.conjfia(j+1)-1)
             if (size(obj.conjfia,1) > 0)
                 obj.conjugates(obj.conjfia(j):obj.conjfia(j+1)-1).printLLatex
             for k = obj.conjfia(j):obj.conjfia(j+1)-1
               %disp("Conjugate Expr")
               %obj.conjugates(k).print
               %obj.conjf(k).print
               %disp('Conjugate Domain')
               %obj.conjd(k).print
             end
             end
           end
           fprintf("\n\n\n\n\n")
           disp("Maximum conjugate")
           obj.maxConjugate.printLLatex

           % size(obj.maxf,2)
           % for i = 1:size(obj.maxf,1)
           %   disp(i)
           %   obj.maxf(i).print
           %   obj.maxd(i).print
           % end
           % %disp("Conjugate Expr")
            % obj.conjf.printL
            % disp('Conjugate Domain')
             %obj.conjd(j).print

         end


         
         function Mprint(obj)
           fprintf("display(inequal({"); 
           obj.d.polygon.printMaple
           fprintf("},x=-5..5,y=-5..5,color=[red,blue,yellow,green],nolines)) \n")
           %return
           disp("Convex Envelope")
%           size(obj.envf)
           fprintf("display(inequal({");
             
             
           for j=1:size(obj.envelope,2)-1 
             
             obj.envelope(j).d.printMaple
           
             
             fprintf(",");
           end
           obj.envelope(size(obj.envelope,2)).d.printMaple

           fprintf("},x=-5..5,y=-5..5,color=[red,blue,yellow,green],nolines)) \n")
           
            
           for j=1:size(obj.envelope,2)

           if (size(obj.conjfia,1) > 0)
                  obj.conjugates(obj.conjfia(j):obj.conjfia(j+1)-1).printM

           end
           end
           % HISTORY: maxConjugate is only populated by .maximum
           % (maximumConjugate); a piece that only went through
           % .convexEnvelope/.conjugate(/.biconjugateP), never .maximum, has
           % an empty maxConjugate here, and functionNDomain.printM
           % unconditionally indexes objL(1) -- errors on an empty array
           % rather than just having nothing further to print.
           if ~isempty(obj.maxConjugate)
             obj.maxConjugate.printM;
           end
         end

         function plotMaxConjugateDomain(obj)
             
             figure;
         
         
             colors = ['b', 'r', 'g', 'm', 'c', 'y'];
             n = 0
             f = obj.maxConjugate (1).f
             c = colors(mod(n,6)+1)

             for i =1:size(obj.maxConjugate,2)
                i
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

%          function plot(obj)
%              %figure;
%              %obj.d.plot;
% 
%            for j=1:size(obj.envf,2) 
% 
%              limits = [min(obj.envd(j).vx),max(obj.envd(j).vx),min(obj.envd(j).vy),max(obj.envd(j).vy)];
%              for i = 1:4
%                  if limits(i) > 10
%                      limits(i) = 10;
%                  end
%                  if limits(i) < -10
%                      limits(i) = -10;
%                  end
%              end
% 
%              obj.envf(j).f = simplify(obj.envf(j).f);
% %              figure;  
% %              obj.envf(j).plot3d (limits);
% %              figure;
% %              obj.d.plot;
% %              obj.envd(j).plot;
% %              obj.envd(j).plotRegion ("Env");
%              if (size(obj.conjfia,1) > 0)
%              for k = obj.conjfia(j):obj.conjfia(j+1)-1
%                %disp("Conjugate Expr")
%                limits = [min(obj.conjd(k).vx),max(obj.conjd(k).vx),min(obj.conjd(k).vy),max(obj.conjd(k).vy)];
%                for i = 1:4
%                  if limits(i) > 10
%                      limits(i) = 10;
%                  end
%                  if limits(i) < -10
%                      limits(i) = -10;
%                  end
%                end
%                %obj.conjd(k).vx
%                %obj.conjd(k).vy
%                %limits
% %                figure;
% %                obj.conjf(k).plot3d (limits);
%              end
%              disp("plot conj domain")
%              figure;
% 
%              for k = obj.conjfia(j):obj.conjfia(j+1)-1
%                limits = [min(obj.conjd(k).vx),max(obj.conjd(k).vx),min(obj.conjd(k).vy),max(obj.conjd(k).vy)];
% %                for i = 1:4
% %                  if limits(i) > 6
% %                      limits(i) = 6;
% %                  end
% %                  if limits(i) < -6
% %                      limits(i) = -6;
% %                  end
% %                end
% %                %disp('Conjugate Domain')
%                %obj.conjd(k).plotByVertex;
% 
%                obj.conjd(k).plot;
%                obj.conjd(k).plotRegion ("Conj");
% 
%              end
% 
% 
%              end
%            end
% %            disp("Maximum conjugate")
%             figure;
%             for i = 1:size(obj.maxf,1)
%               disp(i)
%               %obj.maxd(i).plot;
%               %obj.maxd(i).plotRegion("M"+num2str(i));
%             end
% %            
%          end

         function plotDomain(obj)
             figure;
             obj.d.plot;
             
           
         end
    end

  
    methods % convex
  
        % function obj = setC
        %     obj.envelope(1) = functionNDomain()
        % end

        function obj = setConv(obj)
            vars = obj.f.getVars;
            x = vars(1);
            y = vars(2);
          
            f = symbolicFunction(30*x - 5*y + 4*x*y + 10*x^2 + 5*y^2 - 100,2*x - y + 15);
            l = [x - y/2 - 2, - x - (5*y)/9 - 20/9, x + (7*y)/5 - 2  ];

            r = region(l,vars);
            %obj.envelope(1) = functionNDomain(f,r);
            s0 = (46*x)/19 - (23*y)/19 - 130/19;
            s1 = (7*x)/38 + (3*y)/19 + 5/38;
            s2 = x/76 - y/152 + 15/152;
            dl = x;
            %obj.envelopeExpr(1) = convexExpr(1, s0, s1, s2, dl);
            obj.envelope = [obj.envelope,functionNDomain(f,r)];
            %obj.type = [obj.type,1];
            %obj.envelopeExpr = [obj.envelopeExpr,convexExpr(1, s0, s1, s2, dl)];
            


            f = symbolicFunction(5*x + 2*y  -10);
            l = [x - 2, x + y/2 - 5/2,x + 3*y - 10, 2 - (7*y)/5 - x ];



            r = region(l,vars);
            %obj.envelope(2) = functionNDomain(f,r);
            s0 = 2*y;
            s1 = x - 2;
            s2 = x;
            dl = 0*x+5;
            
            %obj.envelopeExpr(2) = convexExpr(3, s0, s1, s2, dl);
            
           obj.envelope = [obj.envelope,functionNDomain(f,r)];
            %obj.envelopeExpr = [obj.envelopeExpr,convexExpr(3, s0, s1, s2, dl)];
            %obj.type = [obj.type,3];

            f = symbolicFunction(-4*x  -5*y  -20);
            l = [- y - 4 ,- x - 5 ,x + (5*y)/9 + 20/9 ];

            r = region(l,vars);
            %obj.envelope(3) = functionNDomain(f,r);
            s0 = -4*x;
            s1 = y + 4;
            s2 = x;
            dl = 0*x-5;
           % obj.envelopeExpr(3) = convexExpr(3, s0, s1, s2, dl);
            obj.envelope = [obj.envelope,functionNDomain(f,r)];
           % obj.envelopeExpr = [obj.envelopeExpr,convexExpr(3, s0, s1, s2, dl)];
           % obj.type = [obj.type,3];
           
        end

        function obj = convexEnvelope(obj)
          %  obj = obj.setConv;
          %  return
          vars = obj.f.getVars;
          if (size(vars,2)==2)
            x = vars(1);
            y = vars(2);
          else
            disp("not bivariate in 'plq.m'")
            return
          end
          %disp('in convexE 1')
          obj = convexEnvelope1 (obj,x,y);
          % %obj.print
          % return
          %disp('hh')
          obj.envelope =  obj.envelope.unique  
          %obj.print
          [obj.envelope,index,lCh] = obj.envelope.maxEqFun;
          %expr = obj.envelopeExpr;
          %obj.envelopeExpr = expr(index); 
          
          obj.print
          %return
          % fix this loop
          nx = size(obj.envelope,2)
          for ik = 1:5
               ik
               
            lCh = true;
            while lCh
              [obj.envelope,index,lCh] = obj.envelope.maxEqDom;
              obj.envelope.printL
              %pause
              %expr = obj.envelopeExpr;
              %obj.envelopeExpr = expr(index); 
              
            end  
             disp('eq domain')
             size(obj.envelope,2)
            % obj.print
            % return
            % 
            [objL2,index2] = obj.envelope .* obj.envelope;
           
              size(objL2,2)
            [obj.envelope,index] = objL2.maximumPC(index2) ;
            %expr = obj.envelopeExpr;
            %obj.envelopeExpr = expr(index); 
           
              disp('max')
              size(obj.envelope,2)
         %     obj.print
            lCh = true;
            while lCh
              [obj.envelope,index,lCh] = obj.envelope.maxEqDom;
              %expr = obj.envelopeExpr;
              %obj.envelopeExpr = expr(index);
               
          
            end  
            disp('eqd2')
            size(obj.envelope,2)
            if nx == size(obj.envelope,2)
                break;
            end
            nx = size(obj.envelope,2)
          %return 
          end
          %disp('check this')
          %obj.envelope.printL
          %return
          %obj.envelope.printL
          [obj.envelope, index] = obj.envelope.mergeL;
          %expr = obj.envelopeExpr;
          %obj.envelopeExpr = expr(index); 
          
        end

        function obj = convexEnvelope1 (obj,x,y)
            a=sym('a');
            b=sym('b');
              
            % etaV : eta functions corresponding to set obj.d.V
            % etaE : eta functions corresponding to set obj.d.E
            % etaR : domain of etaE - stored as etaR(i,1:3) : [function,lb,ub] => lb <= function <= ub 
            %disp("getEtaFunctions")
            [etaV, etaE, etaR] =  getEtaFunctions (obj,x,y,a,b);
            %disp('eta done')
              obj.d.V
              obj.d.E
              disp('etaV')
              etaV.printL
              disp('etaE')
              etaE.printL
              disp('etaR')
              etaR.printL
             % pause
            % return
            [ix,jx,vix, vjx, ixd, jxd] = feasiblePairs (obj,etaR, a,b);
          

            [envfs,  envds] = solveC (obj, ix,jx,vix, vjx,ixd, jxd, etaV, etaE, etaR,a, b, x, y);
            %disp('solved')
             %  return
            for i = 1:size(envfs,2)
             %   i,size(envfs,2)

              %r = envds(i).simplify;
              r = envds(i).simplifyUnboundedRegion;
              if isempty(r)
                  continue
              end
              %if (r.isFeasible & r.nv > 2)  % added on 29 oct % removed for unbounded region
                  obj.envelope = [obj.envelope,functionNDomain(envfs(i), r)];
                  %obj.envelopeExpr = [obj.envelopeExpr,envxs(i)];
              %end
             
            end
           % disp('out')

        end 

         % returns etah, slope and y intercept of edge
        function [etah, m, yintercept] = edgeInfoInSolve(obj, etaE, etaR, ix, ixd) %, nc, c)
           etah = etaE(ix,ixd);
           if size(obj.d.mE,2) == 1
             m = obj.d.mE;
             yintercept = obj.d.cE;
           else  
             m = obj.d.mE(ix);
             yintercept = obj.d.cE(ix);
           end
           
        end

        %change name
        % gives bound of linear inequality c <= 0
        function [nb,lb, ub, linfeasible] = getBound1 (obj, c,b,nb,lb,ub)
            %disp('getBound1')
          c1 =  c.solve(b);
          
          linfeasible = isempty(c1);
          if (linfeasible)
              disp('here')
                      return;
          end

          s = coeffs(c.f,b);
          %size(c1)
          if (size(c1,1) > 1)
              %isreal(c1(1))
            if isreal(c1(1))
              if(s(end) > 0)
                nb = nb+1;
                lb(nb)=c1(1);
                ub(nb)=c1(2);
              end
              else
                  c.f
                  s
                  s(end)
                  % changed on 18-Aug-2024
                  %disp('here2')
                  if s(end) < 0
                    linfeasible = true;
                    
              end
            end
          else
            if (isreal(c1(1)))
              nb = nb+1;
              if s(end) < 0
                ub(nb)= inf;
                lb(nb) = c1(1);
              else
                lb(nb)= -inf;
                ub(nb) = c1(1);
              end
              
              %if s(end)/s(1) < 0
              %  ub(nb)= inf;
              %  lb(nb) = c1(1);
              %else
              %  lb(nb)= -inf;
              %  ub(nb) = c1(1);
              %end
              else
                  disp('here3')
                  linfeasible = true;
              end
          end
          
        end

        % get bounds depending on region
        function [nb, lb, ub, linfeasible] = getBoundsLinear (obj, etaR,ixd, a,b,av, nb, lb, ub)
            if ixd == 1
              nc = 1  ;
              c(1) = etaR(1)-etaR(2);
            elseif ixd == 2

                nc = 2;
              c(1) = etaR(2)-etaR(1);
              c(2) = etaR(1)-etaR(3);
            elseif ixd == 3
                nc = 1;
              c(1) = etaR(3)-etaR(1);
              
            end
            for i = 1:nc
               % c(i)
              c(i) = c(i).subsVarsPartial([a],[av]);
              [nb,lb, ub, linfeasible] = getBound1 (obj, c(i),b,nb,lb,ub);
              if (linfeasible) 
                return;
              end
            end
            
        end
    
        

        function [envfs, envds, lSol] = solveLinearLinear1(obj, obj0, etah, etaw, etaV, lV, etaE, lE, etaRix, ixd, etaRjx, jxd, x, y, a, b, envfs,  envds) 
            %if (etah == etaw) wont work
            f0 = etah-etaw;
            % if a gets eliminated we exit this routine and try again
            % exchanging a and b
            linfeasible = false;
            av = f0.solve(a);
            % Get a in terms of b 
            if (isempty(av))
                lSol = false;
                %disp('returning')
                return;
            end
            lSol = true;
            %disp('here')
            %bcoeffs = coeffs(av,b); 
            %alpha0 = bcoeffs(1);
            %alpha1 = bcoeffs(2);

            %substitute a in objective and get coefficients of b
            %obj0
            obj0 = obj0.subsVarsPartial([a],[av]);
            objfacts = coeffs(obj0.f,b);
           if size(objfacts,2) == 1
              eta0 = 0;
              eta1 = objfacts(1);  
           else
              eta0 = objfacts(1);
              eta1 = objfacts(2);  
            end
            % obj = b * eta1 + eta0
            %return
            nb = 0 ;
            lb = [];
            ub = [];
            %ixd
            if (ixd > 0)
                %etaRix.printL
               [nb, lb, ub, linfeasible] = getBoundsLinear (obj, etaRix,ixd, a,b,av, nb, lb, ub)
               if (linfeasible) 
                      return;
                  end
            end
            %jxd
            if (jxd > 0)
                %etaRjx.printL
              [nb, lb, ub, linfeasible] = getBoundsLinear (obj, etaRjx,jxd, a,b,av, nb, lb, ub)
              if (linfeasible)  
                      return;
                  end
            end
            % disp("bounds")
            % nb
            % lb
            % ub

%             nb = 1;
%             lb(1) = max(lb);
%             ub(1) = min(ub);
%             if lb(1) >= ub(1)
%                 return
%             end
%             nb0 = nb;
%disp('b4 etak')
            for j=1:size(etaV,2)
               if (lV(j)) 
                 continue;
               end
               etak = etaV(j);
               c = etah - etak; %; % <= 0   easier for substitution
               c = c.subsVarsPartial([a],[av]);
                 [nb,lb, ub, linfeasible] = getBound1 (obj,c,b,nb,lb,ub);
               if (linfeasible) 
                   disp('not feasible1')
                 return;
               end
            end
            
            for j=1:size(etaE,1)
               if (lE(j))
                 continue;
               end
               for k = 1:size(etaE,2) 
                  % disp('in E loop')
                  etak = etaE(j,k)
                  c = etah - etak  % <= 0   easier for substitution
                  c = c.subsVarsPartial([a],[av])
                  if (isZero(c)) 
                      continue
                  end
                  [nb,lb, ub, linfeasible] = getBound1 (obj, c,b,nb,lb,ub);
                 
               
                  if (linfeasible) 
                      disp('not feasible2')
                      return;
                  end
                end
            end
            
            

          %   for i = 1:nb
          %    mlb = lb(i);
          %    mub = ub(i);
            mlb = max(lb)
           mub = min(ub)
           %eta0
           %eta1
           if size(mlb,1) == 0
               envfs = [symbolicFunction(eta1)];
               %envxs = [convexExpr(3,eta0,eta1,intmax)];
               
               envds = [envds, obj.d.polygon ];
               return
           end

            %size(envfs)

            % soln
            % max (lb * eta1 + eta0, ub * eta1 + eta0)
            % if eta1 >= 0   :  ub * eta1 + eta0
            % if eta1 <= 0   :  lb * eta1 + eta0   

            % skipping inf and -inf cases

            if (mlb > mub)
                %disp('infeasible')
            else
               
              if (mub == inf)

                r00  = obj.d.polygon + region(eta1, [x,y]);
                if ~isempty(r00) 
                  
              
                  envfs = [envfs, obj0.subsVarsPartial([b],[mlb])];
                %envxs = [envxs, convexExpr(3,eta0,eta1,mlb)];
                %envds = [envds, region(eta1, [x,y]) ];
                envds = [envds, r00 ];
              end
  %              disp("inf")
  %              envfs(end).print 
  %              envfs(end).print
              elseif (mlb == -inf)
                  %eta1
                  r00 = region(-eta1, [x,y]);
                  %r00.print
                  r00  = obj.d.polygon + region(-eta1, [x,y]);
              if ~isempty(r00) 
              
                envfs = [envfs, obj0.subsVarsPartial([b],[mub])];
               % envxs = [envxs, convexExpr(3,eta0,eta1,mub)];
                envds = [envds, r00 ];
                %envds = [envds, region(-eta1, [x,y]) ];
              end
   %             disp("-inf")
   %             mub
   %             obj0.print
   %             envfs(end).print
              elseif (mlb ~= mub)
                   r00  = obj.d.polygon + region(eta1, [x,y]);
              if ~isempty(r00) 
              
                envfs = [envfs, obj0.subsVarsPartial([b],[mlb])];
                %envxs = [envxs, convexExpr(3,eta0,eta1,mlb)];
                %envds = [envds, region(eta1, [x,y]) ];
                envds = [envds, r00 ];
              end
    %          disp("notinf1")
     %           envfs(end).print
              r00  = obj.d.polygon + region(-eta1, [x,y]);
              if ~isempty(r00) 
              
                envfs = [envfs, obj0.subsVarsPartial([b],[mub])];
                %envxs = [envxs, convexExpr(3,eta0,eta1,mub)];
                %envds = [envds, region(-eta1, [x,y]) ];
                envds = [envds, r00 ];
              end
      %          disp("notinf2")
      %          envfs(end).print
              
              end
              %disp("in linear")
              %envds(end).print
            end
           %  end
           %envfs
            return
            
        end

        
        
        function [envfs, envds] = solveQuadLinear1 (obj, m, a, b, x, y, etah, etaw, etaRix, ixd, etaRjx, jxd, etaV, lV, etaE, lE, ix, envfs,  envds)
          %disp("one quad")
          z = sym('z');
          %av = z - m*b;
          %z = a+m*b;
          av = solve(z-a-m*b,a)
          eq = etah-etaw
          eq = eq.subsVarsPartial([a],[av])
          bv = eq.solve(b)
          linfeasible = false;

          if (isempty(bv))
              return
              % opp case gives same ans so skip
          end
          nb = 0 ;
          lb = [];
          ub = [];
          %%%%%%%%%%%%%%%%%%%%%%55
          
          %%%%%%%%%%%%%%%%%%%%%%
            if ixd == 1  
              nb = nb+1;
              lb(nb) = -inf;
              ub(nb) = etaRix(2).f;
            elseif ixd == 2
              nb = nb+1;
              lb(nb) = etaRix(2).f;
              ub(nb) = etaRix(3).f;
            elseif ixd == 3
              nb = nb+1;
              lb(nb) = etaRix(3).f;
              ub(nb) = inf;

            end

            % jxd will be 0 - vertex, 1,3 vertex of edge
            % if jxd == 1
            %     nb = nb+1;
            %     lb(nb) = -inf;
            %     ubf = etaRjx(1).subsVarsPartial([a],[av])
            %     ubf = ubf.subsVarsPartial([b],[bv])
            %     ubf = ubf.f - etaRjx(2).f
            %     zubf = solve(ubf,z)
            %     %lb(nb) = zubf(1)
            %     ub(nb) = zubf(2)
            %     %pause
            % elseif jxd == 3
            %     nb = nb+1;
            %     ub(nb) = inf;
            %     ubf = etaRjx(1).subsVarsPartial([a],[av])
            %     ubf = ubf.subsVarsPartial([b],[bv])
            %     ubf = etaRjx(2).f - ubf.f
            %     zubf = solve(ubf,z)
            %     lb(nb) = zubf(1)
            %     %pause
            % end

          % if jxd == 1  
          %     nb = nb+1;
          %     lb(nb) = -inf;
          %     ub(nb) = etaRjx(2).f;
          %   elseif jxd == 2
          %     nb = nb+1;
          %     lb(nb) = etaRjx(2).f;
          %     ub(nb) = etaRjx(3).f;
          %   elseif jxd == 3
          %     nb = nb+1;
          %     lb(nb) = etaRjx(3).f;
          %     ub(nb) = inf;
          % 
          %   end  
          nb0 = nb;
          lb
          ub
          

          % solve these to get bounds
          for j=1:size(etaV,2)
            if (lV(j))
              continue;
            end
            etak = etaV(j);
            c = etah - etak;  % <= 0   easier for substitution
            c = c.subsVarsPartial([a],[av]);
            c = c.subsVarsPartial([b],[bv]);
            if isNegativeSqr(c,z)
              %disp ('negative sqr')
              continue
            end
            z0 = c.solve(z);
            % check this part
            % as z is always satisfied
            if isreal(z0(1))
              nb = nb+1;
              lb(nb) = min(z0);
              if (lb(nb)  < lb(1))
                  return
              end
              ub(nb) = max(z0);
              if (ub(nb)  > ub(1))
                  return
              end
              
            end
          end
          
          for j=1:size(etaE,1)
            if (lE(j))
              continue;
            end
            for k = 1:size(etaE,2)
              etak = etaE(j,k);
              %nc = nc+1;
              c = etah - etak;  % <= 0   easier for substitution""
              c = c.subsVarsPartial([a],[av]);
              c = c.subsVarsPartial([b],[bv]);
              %c.print

              % -x^2 always <= 0 hence skipping
              if isNegativeSqr(c,z)
                  %disp ('negative sqr')
                  continue
                  
              end 
              cz = coeffs(c.f,z);
              z1 = roots(cz);
              %z0 = c.solve(z)
              %disp('size')
              %size(z0)
              nz = 0;
              for iz = 1:size(z1)
                  if (~isreal(z1(iz)))
                      continue
                  end
                  nz = nz + 1;
                  z0(nz) = z1(iz);
              end
              if nz > 0
                nb = nb+1;
                lb(nb) = min(z0);
                ub(nb) = max(z0);
              end
              
            end
          end
          %lb
          %ub
          
          nb1 = nb0; 
          
            for i = nb0+1:nb
                wB = false;
                for j = 1:nb0
                    if (lb(i) < lb(j))
                        wB = true;
                        break;
                    end
                    if (ub(i) > ub(j))
                        wB = true;
                        break;
                    end
                end
                if wB
                    continue;
                end
                nb1 = nb1+1;
                lb(nb1) = lb(i);
                ub(nb1) = ub(i);
            end
            nb = nb1;
          %lb
          %ub
           
          %lb
          %ub
          av
          bv
          
          
          
          obj0 = etah + symbolicFunction(a*x+b*y)
          obj0 = obj0.subsVarsPartial([a],[av])
          obj0 = obj0.subsVarsPartial([b],[bv])
          simplifyFraction(obj0.f)
          
          
          deg = polynomialDegree(obj0.f,z);
          if (deg ~= 2) 
              disp("check degree of polynomial in z in quad-lin")
              return
          end
          objfacts = coeffs(obj0.f,z)
          
          if (deg+1 == size(objfacts,2))
            psi0 = objfacts(1);
            psi1 = objfacts(2)/2;
            psi2 = -objfacts(3);
          elseif (deg == size(objfacts,2))
            z0 = simplify(obj0.f - objfacts(end)*z^2);
            %polynomialDegree(z0,z)
            if polynomialDegree(z0,z) == 1
              psi0=0;
              psi1 = objfacts(1)/2;
            else
              psi0=objfacts(1);
              psi1 =0;

            end

            psi2 = -objfacts(2);
          else
            psi0=0;
            psi1=0;
            psi2 = -objfacts(1);
          end

          psi0
          psi1
          psi2
          %pause
          lb
          ub

          mlb = max(lb);
          mub = min(ub);
          %disp("in here")  
            %for i = 1:size(lb,2)
           
          if positivePsi (obj,-psi2,x,y)==1
            disp("-psi2 +ve check this")
            mlb,mub,psi2,psi0
           % psi1
           % psi2
           if mlb < mub
            f0 = -mub^2*psi2 +2*mub*psi1+psi0 ;
              r0 = simplify(psi2);

              r00  = obj.d.polygon + r0;
              
              if ~isempty(r00) 
              
              envfs = [envfs, symbolicFunction(f0)];
              %envxs = [envxs, convexExpr(2,psi0,psi1,psi2)];
              %envxs = [envxs, convexExpr(3,psi0,-mub^2*psi2 +2*mub*psi1,1)];
              
              %envds = [envds, region(r0, [x,y])];
              envds = [envds, region(r00, [x,y])];
              end

              f0 = -mlb^2*psi2 +2*mlb*psi1+psi0;
              r0 = simplify(psi2);
              r00  = obj.d.polygon + r0;
              if ~isempty(r00) 
              
              envfs = [envfs, symbolicFunction(f0)];
              %envxs = [envxs, convexExpr(2,psi0,psi1,psi2)];
              %envxs = [envxs, convexExpr(3,psi0,-mlb^2*psi2 +2*mlb*psi1,1)];
              
              %envds = [envds, region(r0, [x,y])];
              envds = [envds, region(r00, [x,y])];
              end
           end
          else
          
%           mlb = min(lb);
%           mub = max(ub);

          % Why returning here - check
          % fishy but needed
          %for i = 1:size(lb,2)
          %    mlb = lb(i);
          %     mub = ub(i);

          
%           if mlb == -inf
%               %continue
%               return
%               disp('returning')
%               return
%               
%           elseif mub == inf
%               %continue
%               return
%               
%           end
%             %for i = 1:size(lb,2)
             %   mlb = lb(i);
             %   mub = ub(i);
          %   mlb = max(lb)
          %   mub = min(ub)
            if (mlb >= mub)
%             %   disp('infeasible')
            else
          

            %for i = 1:nb
              %f0 = simplify(psi1^2/psi2+psi0);
              f0 = symbolicFunction(psi1^2,psi2) + symbolicFunction(psi0)
              r0 = -simplify(mub*psi2-psi1);
              r1 = simplify(mlb*psi2-psi1);
              r00  = obj.d.polygon + region([r0,r1], [x,y]);
              if ~isempty(r00) 
                envfs = [envfs, f0]
                envds = [envds, r00];
              end
              f0 = -mub^2*psi2 +2*mub*psi1+psi0 ;
              r0 = simplify(mub*psi2-psi1);
              r00  = obj.d.polygon + region(r0, [x,y]);
              disp('lq')
              r00.print
              r00 = r00.simplify
              if ~isempty(r00) 
                envfs = [envfs, symbolicFunction(f0)];
                envds = [envds, r00];
              end 
              f0 = -mlb^2*psi2 +2*mlb*psi1+psi0;
              r0 = simplify(-mlb*psi2+psi1);
              r00  = obj.d.polygon + region(r0, [x,y]);
              disp('lq2')
              r00.print
              r00 = r00.simplify
            %  pause
               envfs
              if ~isempty(r00) 
              
                envfs = [envfs, symbolicFunction(f0)];
                envds = [envds, r00];
              end
              %size(envds)
            end
           % end
            
            %
            %envfs = [envfs, f0];
            %envds = [envds, r0];
            %r0 = simplify(etaR(ix,3).f*psi2-psi1);
            %envfs = [envfs, f1];
            %envds = [envds, r0];
          end
         
        end
           
        function [envfs, envds] = solveQuadQuad1(obj, etah, x, y, a, b, alpha0, alpha1, mh, qh, mw, ix, jx, etaR, etaV, lV, etaE, lE, envfs, envds)
     %       disp("quadquad")
      return
          av = alpha1*b + alpha0
          % alpha1 = sqrt (mh * mw) so > 0, thus coef of b is + - no need
          lb(1) = (etaR(ix,2).f - alpha0)/(alpha1+mh)
          ub(1) = (etaR(ix,3).f - alpha0)/(alpha1+mh)
          nb = 2;
          lb(2) = (etaR(jx,2).f - alpha0)/(alpha1+mw)
          ub(2) = (etaR(jx,3).f - alpha0)/(alpha1+mw)
          double(lb)
          double(ub)
          for j=1:size(etaV,2)
            if (lV(j))
              continue;
            end
            etak = etaV(j);
           % nc = nc+1;
            % make this a function
            c = etah - etak;  % <= 0   easier for substitution
            c = c.subsVarsPartial([a],[av]);
           [l,qc] = c.quadterm (b);

            if ~l
                disp('no quad term in quadquad fishy')
            end
            z0 = c.solve(b);
            if isreal(z0(1))
              if qc > 0
                nb = nb+1;
                if z0(1) < z0(2)
                  lb(nb) = z0(1);
                  ub(nb) = z0(2);
                else
                  ub(nb) = z0(1);
                  lb(nb) = z0(2);
                end
              

            else
              
                if z0(1) < z0(2)
                  nb = nb+1;  
                  ub(nb) = z0(1);
                  lb(nb) = -inf;
                  nb = nb+1;  
                  ub(nb) = inf;
                  lb(nb) = z0(2);
                else
                   nb = nb+1;  
                  ub(nb) = z0(2);
                  lb(nb) = -inf;
                  nb = nb+1;  
                  ub(nb) = inf;
                  lb(nb) = z0(1);
                end
              end
            end
          end
          for j=1:size(etaE,1)
            if (lE(j))
              continue;
            end
            for k = 1:size(etaE,2)
              etak = etaE(j,k);
              c = etah - etak;  % <= 0   easier for substitution
              c = c.subsVarsPartial([a],[av]);
              [l,qc] = c.quadterm (b);

            if ~l
                disp('no quad term in quadquad fishy')
            end
            z0 = c.solve(b);
            if isreal(z0(1))
              if qc > 0
                nb = nb+1;
                if z0(1) < z0(2)
                  lb(nb) = z0(1);
                  ub(nb) = z0(2);
                else
                  ub(nb) = z0(1);
                  lb(nb) = z0(2);
                end
              

            else
              
                if z0(1) < z0(2)
                  nb = nb+1;  
                  ub(nb) = z0(1);
                  lb(nb) = -inf;
                  nb = nb+1;  
                  ub(nb) = inf;
                  lb(nb) = z0(2);
                else
                   nb = nb+1;  
                  ub(nb) = z0(2);
                  lb(nb) = -inf;
                  nb = nb+1;  
                  ub(nb) = inf;
                  lb(nb) = z0(1);
                end
              end
            end
            end
          end
          lb
          ub
          mlb = max(lb)
          mub = min(ub)
          
          %if mlb < mub
          %checked on 29/8
          mh,qh,alpha0,alpha1
          psi0 = simplifyFraction (-(alpha0-qh)^2 /(4*mh) + alpha0*x)
          psi1 = simplifyFraction (-(alpha1+mh)*(alpha0-qh)/(4*mh) + (y - qh + alpha1*x))/2
          psi2 = simplifyFraction ((alpha1+mh)^2/(4*mh))
          pause
          %return
          %for i = 1: nb
          %    mlb = lb(i);
           %   mub = ub(i);
           %mlb
           %mub
             if (lb(2) > ub(1))
                 return
             end
              if double(mlb) >= double(mub)
            %      disp('infeasible')
                  return
              end
              %else
              f0 = symbolicFunction(psi1^2,psi2) + symbolicFunction(psi0)
              r0 = -simplifyFraction(mub*psi2-psi1)
              r1 = simplifyFraction(mlb*psi2-psi1)
              r00  = obj.d.polygon + region([r0], [x,y]) ;
%              r0.print
%              r1.print
              r00.print
              r00  = r00 + region([r1], [x,y]);
              r00.print
              if ~isempty(r00) 
                  envfs = [envfs, f0];
                  envds = [envds, r00];
              end
              f0 = simplify(-mlb^2*psi2 + 2*mlb*psi1 + psi0)
              r0 = simplify(psi1-mlb*psi2);
              r00  = obj.d.polygon + region([r0], [x,y]);
              r00.print
              if ~isempty(r00) 
                  envfs = [envfs, symbolicFunction(f0)];
                  envds = [envds, r00];
              end

              f0 = simplify(-mub^2*psi2 + 2*mub*psi1 + psi0)
              r0 = simplify(mub*psi2-psi1);
              r00  = obj.d.polygon + region([r0], [x,y]);
              r00.print
              if ~isempty(r00) 
                  envfs = [envfs, symbolicFunction(f0)];
                  envds = [envds, r00];
              end
              
              pause
          
        end


        function [envfs,  envds] = solveConstLinear(obj, obj0, etah, etaw, x, y, a, b, ixd, jxd, etaRix, etaRjx, etaV, lV, etaE, lE, envfs, envds)
            %return
          vars = etaw.getVars;
          linfeasible = false;
          nb = 0 ;
          lb = [];
          ub = [];
          if size(vars,2) == 1
              
              if vars(1) == a
               %   etaw
               %   etah
                  t0 = etaw-etah;
                  av = solve(t0.f,a);
                av = etah.f;
                obj0 = obj0.subsVarsPartial([a],[av]);
                objfacts = coeffs(obj0.f,b);
                if size(objfacts,2) == 1
                eta0 = 0;
                eta1 = objfacts(1);  
                
                else
                eta0 = objfacts(1);
                eta1 = objfacts(2);  
                end
                if (ixd > 0)
                  [nb, lb, ub, linfeasible] = getBoundsLinear (obj, etaRix,ixd, a,b,av, nb, lb, ub);
                  %linfeasible
                  if (linfeasible)  
                    return;
                  end
                end
                if (jxd > 0)
                  [nb, lb, ub, linfeasible] = getBoundsLinear (obj, etaRjx,jxd, a,b,av, nb, lb, ub);
                  %linfeasible
                  if (linfeasible)  
                    return;
                  end
                end
                for j=1:size(etaV,2)
                  if (lV(j)) 
                    continue;
                  end
                  etak = etaV(j);
                  c = etah - etak ; % <= 0   easier for substitution
                  c = c.subsVarsPartial([a],[av]);
                  [nb,lb, ub, linfeasible] = getBound1 (obj, c,b,nb,lb,ub);
                  if (linfeasible) 
                    return;
                  end
                end
                for j=1:size(etaE,1)
                  if (lE(j))
                    continue;
                  end
                  for k = 1:size(etaE,2) 
                    etak = etaE(j,k);
                    c = etah - etak;  % <= 0   easier for substitution
                    c = c.subsVarsPartial([a],[av]);
                    if (isZero(c)) 
                      continue
                    end
                    [nb,lb, ub, linfeasible] = getBound1 (obj,c,b,nb,lb,ub);
                 
               
                    if (linfeasible) 
                      return;
                    end
                  end
                end
                mlb = max(lb);
                mub = min(ub);
          
                % soln
                % max (lb * eta1 + eta0, ub * eta1 + eta0)
                % if eta1 >= 0   :  ub * eta1 + eta0
                % if eta1 <= 0   :  lb * eta1 + eta0   

                %objfacts
                if (double(mlb) > double(mub))
          %      disp('infeasible')
                else
                if size(objfacts,2) == 2
                if (mub == inf)
                    r00  = obj.d.polygon + region(objfacts(2), [x,y]);
              if ~isempty(r00)
                  envfs = [envfs, obj0.subsVarsPartial([b],[mlb])];
                  %envxs = [envxs, convexExpr(3,eta0,eta1,mlb)];
                  %envds = [envds, region(objfacts(2), [x,y]) ];
                  envds = [envds, r00 ];
              end
                elseif (mlb == -inf)
                     r00  = obj.d.polygon + region(-objfacts(2), [x,y]);
              if ~isempty(r00)
                  envfs = [envfs, obj0.subsVarsPartial([b],[mub])];
                  %envxs = [envxs, convexExpr(3,eta0,eta1,mub)];
                  %envds = [envds, region(-objfacts(2), [x,y]) ];
                  envds = [envds, r00 ];
              end
                elseif (mlb ~= mub)
                     r00  = obj.d.polygon + region(objfacts(2), [x,y]);
              if ~isempty(r00)
                  envfs = [envfs, obj0.subsVarsPartial([b],[mlb])];
                 % envxs = [envxs, convexExpr(3,eta0,eta1,mlb)];
                 % envds = [envds, region(objfacts(2), [x,y]) ];
                  envds = [envds, r00 ];
              end 
                 r00  = obj.d.polygon + region(-objfacts(2), [x,y]);
              if ~isempty(r00)
                  envfs = [envfs, obj0.subsVarsPartial([b],[mub])];
                  %envxs = [envxs, convexExpr(3,eta0,eta1,mub)];
                 % envds = [envds, region(-objfacts(2), [x,y]) ];
                   envds = [envds, r00 ];
              end
                end
                else
                if (mub == inf)
                    r00  = obj.d.polygon + region(objfacts(1), [x,y]);
              if ~isempty(r00)
                  envfs = [envfs, obj0.subsVarsPartial([b],[mlb])];
                 % envxs = [envxs, convexExpr(3,eta0,eta1,mlb)];
                  %envds = [envds, region(objfacts(1), [x,y]) ];
                   envds = [envds, r00 ];
              end
                elseif (mlb == -inf)
                    r00  = obj.d.polygon + region(-objfacts(1), [x,y]);
              if ~isempty(r00)
                  envfs = [envfs, obj0.subsVarsPartial([b],[mub])];
                  %envxs = [envxs, convexExpr(3,eta0,eta1,mub)];
                  %envds = [envds, region(-objfacts(1), [x,y]) ];
                   envds = [envds, r00 ];
              end
                elseif (mlb ~= mub)
                    r00  = obj.d.polygon + region(objfacts(1), [x,y]);
              if ~isempty(r00)
                  envfs = [envfs, obj0.subsVarsPartial([b],[mlb])];
                  %envxs = [envxs, convexExpr(3,eta0,eta1,mlb)];
                  %envds = [envds, region(objfacts(1), [x,y]) ];
                   envds = [envds, r00 ];
              end
                  
              r00  = obj.d.polygon + region(-objfacts(1), [x,y]);
              if ~isempty(r00)
                  envfs = [envfs, obj0.subsVarsPartial([b],[mub])];
                  %envxs = [envxs, convexExpr(3,eta0,eta1,mub)];
                %  envds = [envds, region(-objfacts(1), [x,y]) ];
                   envds = [envds, r00 ];
              end
                end
                
                end
              %disp("in linear")
              %envds(end).print
            end
            
            
              else
                 %etaw
                 % etah
                  t0 = etaw-etah;
                  bv = solve(t0.f,b);
                  %return
                obj0 = obj0.subsVarsPartial([b],[bv]);
                objfacts = coeffs(obj0.f,a);
                if size(objfacts,2) == 1
                eta0 = 0;
                eta1 = objfacts(1);  
                
                else
                eta0 = objfacts(1);
                eta1 = objfacts(2);  
                end
                if (ixd > 0)
                  [nb, lb, ub, linfeasible] = getBoundsLinear (obj, etaRix,ixd, b,a,bv, nb, lb, ub);
                  %linfeasible
                  if (linfeasible)  
                    return;
                  end
                end
                if (jxd > 0)
                  [nb, lb, ub, linfeasible] = getBoundsLinear (obj, etaRjx,jxd, b,a,bv, nb, lb, ub);
                  %linfeasible
                  if (linfeasible)  
                    return;
                  end
                end
                for j=1:size(etaV,2)
                  if (lV(j)) 
                    continue;
                  end
                  etak = etaV(j);
                  c = etah - etak ; % <= 0   easier for substitution
                  c = c.subsVarsPartial([b],[bv]);
                  [nb,lb, ub, linfeasible] = getBound1 (obj, c,a,nb,lb,ub);
                  if (linfeasible) 
                    return;
                  end
                end
                for j=1:size(etaE,1)
                  if (lE(j))
                    continue;
                  end
                  for k = 1:size(etaE,2) 
                    etak = etaE(j,k);
                    c = etah - etak;  % <= 0   easier for substitution
                    c = c.subsVarsPartial([b],[bv]);
                    if (isZero(c)) 
                      continue
                    end
                    [nb,lb, ub, linfeasible] = getBound1 (obj,c,a,nb,lb,ub);
                 
               
                    if (linfeasible) 
                      return;
                    end
                  end
                end
                mlb = max(lb);
                mub = min(ub);
          
                % soln
                % max (lb * eta1 + eta0, ub * eta1 + eta0)
                % if eta1 >= 0   :  ub * eta1 + eta0
                % if eta1 <= 0   :  lb * eta1 + eta0   

                %objfacts
                if (double(mlb) > double(mub))
          %      disp('infeasible')
                else
                if size(objfacts,2) == 2
                if (mub == inf)
                    r00  = obj.d.polygon + region(objfacts(2), [x,y]);
              if ~isempty(r00)
                  envfs = [envfs, obj0.subsVarsPartial([a],[mlb])];
                %  envxs = [envxs, convexExpr(3,eta0,eta1,mlb)];
                 % envds = [envds, region(objfacts(2), [x,y]) ];
                  envds = [envds, r00 ];
              end
                elseif (mlb == -inf)
                    r00  = obj.d.polygon + region(-objfacts(2), [x,y]);
              if ~isempty(r00)
                  envfs = [envfs, obj0.subsVarsPartial([a],[mub])];
                 % envxs = [envxs, convexExpr(3,eta0,eta1,mub)];
                  %envds = [envds, region(-objfacts(2), [x,y]) ];
                   envds = [envds, r00 ];
              end
                elseif (mlb ~= mub)
                    r00  = obj.d.polygon + region(objfacts(2), [x,y]);
              if ~isempty(r00)
                  envfs = [envfs, obj0.subsVarsPartial([a],[mlb])];
                  %envxs = [envxs, convexExpr(3,eta0,eta1,mlb)];
                  %envds = [envds, region(objfacts(2), [x,y]) ];
                  envds = [envds, r00 ];
              end 
              r00  = obj.d.polygon + region(-objfacts(2), [x,y]);
              if ~isempty(r00)
                  envfs = [envfs, obj0.subsVarsPartial([a],[mub])];
                  %envxs = [envxs, convexExpr(3,eta0,eta1,mub)];
                  %envds = [envds, region(-objfacts(2), [x,y]) ];
                   envds = [envds, r00 ];
              end
                end
                else
                if (mub == inf)
                    r00  = obj.d.polygon + region(objfacts(1), [x,y]);
              if ~isempty(r00)
              
                  envfs = [envfs, obj0.subsVarsPartial([a],[mlb])];
                  %envxs = [envxs, convexExpr(3,eta0,eta1,mlb)];
                  %envds = [envds, region(objfacts(1), [x,y]) ];
                   envds = [envds, r00 ];
              end
              
              elseif (mlb == -inf)
                    r00  = obj.d.polygon + region(-objfacts(1), [x,y]);
              if ~isempty(r00)
              
                  envfs = [envfs, obj0.subsVarsPartial([a],[mub])];
                  %envxs = [envxs, convexExpr(3,eta0,eta1,mub)];
                  %envds = [envds, region(-objfacts(1), [x,y]) ];
                   envds = [envds, r00];
              end
              
              elseif (mlb ~= mub)
                    r00  = obj.d.polygon + region(objfacts(1), [x,y]);
              if ~isempty(r00)
              
                  envfs = [envfs, obj0.subsVarsPartial([a],[mlb])];
                  %envxs = [envxs, convexExpr(3,eta0,eta1,mlb)];
                  %envds = [envds, region(objfacts(1), [x,y]) ];
                   envds = [envds, r00 ];
              end
                  
              r00  = obj.d.polygon + region(-objfacts(1), [x,y]);
              if ~isempty(r00)
              
                  envfs = [envfs, obj0.subsVarsPartial([a],[mub])];
                  %envxs = [envxs, convexExpr(3,eta0,eta1,mub)];
                  %envds = [envds, region(-objfacts(1), [x,y]) ];
                   envds = [envds, r00 ];
              end
              
              end
                
                end
              %disp("in linear")
              %envds(end).print
            end
              end
                
          else
          
              %disp('cl2')
              linfeasible = false;
            f0 = etah-etaw;
            av = f0.solve(a); 
            obj0 = obj0.subsVarsPartial([a],[av]);
            objfacts = coeffs(obj0.f,b);
                if size(objfacts,2) == 1
                eta0 = 0;
                eta1 = objfacts(1);  
                
                else
                eta0 = objfacts(1);
                eta1 = objfacts(2);  
                end
            if (ixd > 0)
              [nb, lb, ub, linfeasible] = getBoundsLinear (obj, etaRix,ixd, a,b,av, nb, lb, ub);
              %linfeasible
              if (linfeasible)  
                return;
              end
            end
            if (jxd > 0)
              [nb, lb, ub, linfeasible] = getBoundsLinear (obj, etaRjx,jxd, a,b,av, nb, lb, ub);
              %linfeasible
              if (linfeasible)  
                return;
              end
            end
            for j=1:size(etaV,2)
              if (lV(j)) 
                continue;
              end
              etak = etaV(j);
              c = etah - etak ; % <= 0   easier for substitution
              c = c.subsVarsPartial([a],[av]);
              [nb,lb, ub, linfeasible] = getBound1 (obj,c,b,nb,lb,ub);
              if (linfeasible) 
                return;
              end
            end
            for j=1:size(etaE,1)
              if (lE(j))
                continue;
              end
              for k = 1:size(etaE,2) 
                etak = etaE(j,k);
                c = etah - etak;  % <= 0   easier for substitution
                c = c.subsVarsPartial([a],[av]);
                if (isZero(c)) 
                  continue
                end
                [nb,lb, ub, linfeasible] = getBound1 (obj,c,b,nb,lb,ub);
                
               
                if (linfeasible) 
                  return;
                end
              end
            end
            mlb = max(lb);
            mub = min(ub);
          
                % soln
                % max (lb * eta1 + eta0, ub * eta1 + eta0)
                % if eta1 >= 0   :  ub * eta1 + eta0
                % if eta1 <= 0   :  lb * eta1 + eta0   

            
            if (mlb > mub)
          %      disp('infeasible')
            else
               
                if (mub == inf)
               %     disp('ub inf')
               r00  = obj.d.polygon + region(eta1, [x,y]);
              if ~isempty(r00)
              
                  envfs = [envfs, obj0.subsVarsPartial([b],[mlb])];
                %  envxs = [envxs, convexExpr(3,eta0,eta1,mlb)];
                  %envds = [envds, region(eta1, [x,y]) ];
                  envds = [envds, r00 ];
              end
                %  envfs(end).print
                %  envds(end).print
                elseif (mlb == -inf)
                    r00  = obj.d.polygon + region(-objfacts(2), [x,y]);
              if ~isempty(r00)
              
                  envfs = [envfs, obj0.subsVarsPartial([b],[mub])];
              %    envxs = [envxs, convexExpr(3,eta0,eta1,mub)];
              %    envds = [envds, region(-objfacts(2), [x,y]) ];
                  envds = [envds, r00 ];
              end
              
              elseif (mlb ~= mub)
                  r00  = obj.d.polygon + region(objfacts(2), [x,y]);
              if ~isempty(r00)
              
                  envfs = [envfs, obj0.subsVarsPartial([b],[mlb])];
              %    envxs = [envxs, convexExpr(3,eta0,eta1,mlb)];
              %    envds = [envds, region(objfacts(2), [x,y]) ];
                  envds = [envds, r00 ];
              end
                  r00  = obj.d.polygon + region(-objfacts(2), [x,y]);
              if ~isempty(r00)
              
              
                  envfs = [envfs, obj0.subsVarsPartial([b],[mub])];
             %     envxs = [envxs, convexExpr(3,eta0,eta1,mub)];
              %    envds = [envds, region(-objfacts(2), [x,y]) ];
                  envds = [envds, r00 ];
              end
                  
              
                end
              end
              return
          end
       
            
            
        end

        function [envfs, envds] = solveConstQuad  (obj, etah, etaw, x, y, a, b, alpha0, alpha1, mh, qh, ix, jx, etaR, etaV, lV, etaE, lE, envfs, envds)
          linfeasible = false;
          f0 = etah-etaw;
          av = f0.solve(a); 
          obj0 = obj0.subsVarsPartial([a],[av]);
          objfacts = coeffs(obj0.f,b);
          if size(objfacts,2) == 1
            eta0 = 0;
            eta1 = objfacts(1);  
                
          else
            eta0 = objfacts(1);
            eta1 = objfacts(2);  
          end
          if (ixd > 0)
            [nb, lb, ub, linfeasible] = getBoundsLinear (obj, etaRix,ixd, a,b,av, nb, lb, ub);
            %linfeasible
            if (linfeasible)  
               return;
            end
          end
          if (jxd > 0)
            [nb, lb, ub, linfeasible] = getBoundsLinear (obj, etaRjx,jxd, a,b,av, nb, lb, ub);
            %linfeasible
            if (linfeasible)  
              return;
            end
          end
            for j=1:size(etaV,2)
              if (lV(j)) 
                continue;
              end
              etak = etaV(j);
              c = etah - etak ; % <= 0   easier for substitution
              c = c.subsVarsPartial([a],[av]);
              [nb,lb, ub, linfeasible] = getBound1 (obj,c,b,nb,lb,ub);
              if (linfeasible) 
                return;
              end
            end
            for j=1:size(etaE,1)
              if (lE(j))
                continue;
              end
              for k = 1:size(etaE,2) 
                etak = etaE(j,k);
                c = etah - etak;  % <= 0   easier for substitution
                c = c.subsVarsPartial([a],[av]);
                if (isZero(c)) 
                  continue
                end
                [nb,lb, ub, linfeasible] = getBound1 (obj,c,b,nb,lb,ub);
                
               
                if (linfeasible) 
                  return;
                end
              end
            end
            mlb = max(lb);
            mub = min(ub);
          
         end

        % solves all subproblems
        function [envfs, envds] = solveC (obj, ix,jx,vix, vjx, ixd, jxd, etaV, etaE, etaR, a, b, x, y)
          envfs = [];
          %envxs = [];
          envds = [];
          ix
          jx
%           disp("ix")
%           size(ix,2)
          for i=1:size(ix,2)
              i
             % if i ~= 18
             %     continue
             % end
             i00 = size(envfs,2);
             lV = []; 
             for j = 1:size(etaV,2)
               lV(j) = false;
             end
             lE=[];
             for j = 1:size(etaE,1)
               lE(j) = false;
             end
            
              if (vix(i)==1)
                  [etah, mh, qh] = edgeInfoInSolve(obj, etaE, etaR, ix(i), ixd(i)) 

                  lE(ix(i)) = true;
              else
                  etah = etaV(ix(i));
                  lV(ix(i)) = true;
              end
              if (vjx(i)==1)
                  
                  [etaw, mw, qw] = edgeInfoInSolve(obj, etaE, etaR, jx(i), jxd(i)) 
                  
                  lE(jx(i)) = true;
              else
                  etaw = etaV(jx(i));
                  lV(jx(i)) = true;
              end
              
              % etah
              % etaw
             % ix(i)
             % jx(i)
             % ixd(i)
             % jxd(i)
              if (etah == etaw)
                  continue;
              end
              degreeh = polynomialDegree(etah.f);
              degreew = polynomialDegree(etaw.f);
            %  etah
            %  etaw
              if (degreeh==0 & degreew==1)
             %       disp("const-lin")
            %        continue
                    obj0 = etah + symbolicFunction(a*x+b*y);
                    if ixd(i) == 0
                        etaRi=symbolicFunction;
                    else
                        etaRi = etaR(ix(i),:);
                    end
                    if jxd(i) == 0
                        etaRj=symbolicFunction;
                    else
                        etaRj = etaR(jx(i),:);
                    end
                   
                    [envfs, envds]  = solveConstLinear(obj, obj0, etah, etaw, x, y, a, b, ixd(i), jxd(i), etaRi, etaRj, etaV, lV, etaE, lE, envfs,  envds);
                    disp('const-lin')
%                    envfs.printL
             
                  
              end 

              if (degreeh==1 & degreew==0)
              %      disp("const-lin")
            %        continue
                    obj0 = etah + symbolicFunction(a*x+b*y);
                    if ixd(i) == 0
                        etaRi=symbolicFunction;
                    else
                        etaRi = etaR(ix(i),:);
                    end
                    if jxd(i) == 0
                        etaRj=symbolicFunction;
                    else
                        etaRj = etaR(jx(i),:);
                    end
                   
                    [envfs, envds]  = solveConstLinear(obj, obj0, etaw, etah, x, y, a, b, jxd(i), ixd(i), etaRj, etaRi, etaV, lV, etaE, lE, envfs,  envds);
                    disp('lin-const')
             %       envfs.printL
             
                  
              end 

              if (degreeh==1 & degreew==1)
                    % disp("lin-lin")
                   % continue
                    %objective function set here as we can exchange a and b
                    %if required

                    obj0 = etah + symbolicFunction(a*x+b*y);
                    if ixd(i) == 0
                        etaRi=symbolicFunction;
                    else
                        etaRi = etaR(ix(i),:);
                    end
                    if jxd(i) == 0
                        etaRj=symbolicFunction;
                    else
                        etaRj = etaR(jx(i),:);
                    end
                   
                    [envfs,  envds, lSol] = solveLinearLinear1(obj, obj0, etah, etaw, etaV, lV, etaE, lE, etaRi, ixd(i), etaRj, jxd(i), x, y, a, b, envfs, envds) ;
                    
                    if (~lSol)
                        [envfs,  envds, lSol] = solveLinearLinear1(obj, obj0, etah, etaw, etaV, lV, etaE, lE, etaRi, ixd(i), etaRj, jxd(i), x, y, b, a, envfs,  envds) ;
                    end

                  
                end 
                
                if ixd(i) == 0
                  etaRi=symbolicFunction;
                else
                  etaRi = etaR(ix(i),:);
                end
                if jxd(i) == 0
                  etaRj=symbolicFunction;
                else
                  etaRj = etaR(jx(i),:);
                end
                if (degreeh==2 & degreew==1)
                    %disp("quad-lin")
                    %continue;
                    
                    [envfs,  envds] = solveQuadLinear1 (obj, mh, a, b, x, y, etah, etaw, etaRi, ixd(i), etaRj, jxd(i), etaV, lV, etaE, lE, ix(i), envfs, envds);
            % flipping a,b gives same answer
                end
                if (degreeh==2 & degreew==0)
                %    disp("const-quad")
                    
                    [envfs,  envds] = solveQuadLinear1 (obj, mh, a, b, x, y, etah, etaw, etaRi, ixd(i), etaRj, jxd(i), etaV, lV, etaE, lE, ix(i), envfs, envds);
                    
                end
                
                if (degreeh==1 & degreew==2)
                    disp("lin-quad")
                    %continue;
                

                    
                    [envfs, envds] = solveQuadLinear1 (obj, mw, a, b, x, y, etaw, etah, etaRj, jxd(i), etaRi, ixd(i), etaV, lV, etaE, lE, jx(i), envfs, envds);
                    
                end

                if (degreeh==0 & degreew==2)
                  %  disp("const-quad")
                    
                    [envfs, envds] = solveQuadLinear1 (obj, mw, a, b, x, y, etaw, etah, etaRj, jxd(i), etaRi, ixd(i), etaV, lV, etaE, lE, jx(i), envfs, envds);
                    
                end
                            
                if (degreeh==2 & degreew==2)
                   % disp("quad-quad")
        %            etah.print
        %            etaw.print
        %            mh
        %            mw

                    %continue;
                    %obj0 = etah + symbolicFunction(a*x+b*y);

                    %checking av = alpha1*b+alpha0
                    %f = etah - etaw
                    %av = f.solve(a)

                    %checking etah
                    %simplify(-((a+mh*b-qh)^2)/(4*mh) -b*qh)
                    if (mh == mw)
                     %  disp('equal')
                       %av = (2*mh*b+qh+qw)/2;
                       alpha1 = mh;
                       alpha0 = (qh+qw)/2;
                       [envfs, envds] = solveQuadQuad1(obj, etah, x, y, a, b, alpha0, alpha1, mh, qh, mw, ix(i), jx(i), etaR, etaV, lV, etaE, lE, envfs, envds);
                       
                    else
                    %   disp("mh /= mw")
                     %  disp('slopes')
                     %  qh
                     %  qw
                     %  mh
                     %  mw


                       av = (qw*mh-qh*mw+sqrt(mh*mw)*((mh-mw)*b+qh-qw))/(mh-mw)
                       alpha = coeffs(av,b)


                       if (size(alpha,2) == 2)
                           alpha1 = alpha(2);
                           alpha0 = alpha(1);
                       else
                           alpha1 = alpha(1);
                           alpha0 = 0;
                       end

                       alpha0
                       alpha1
                       %alpha0 = (qw-qh)*(1-sqrt(mh*mw)/(mh-mw))
                       %alpha1 = sqrt(mh*mw)
                       %disp('checking av')
                       % should evaluate to 0
                       %f = etah - etaw
                       %f = f.subsF([a],[alpha1*b + alpha0])
                       %simplify(f.f)
                %       size(envfs)
                                        
                       [envfs, envds] = solveQuadQuad1(obj, etah, x, y, a, b, alpha0,  alpha1, mh, qh, mw, ix(i), jx(i), etaR, etaV, lV, etaE, lE, envfs, envds);
                %       size(envfs)
                       disp("second")
                       % dont recalculate alphas
                       av = (qw*mh-qh*mw-sqrt(mh*mw)*((mh-mw)*b+qh-qw))/(mh-mw);
                       alpha = coeffs(av,b);
                       if (size(alpha,2) == 2)
                           alpha1 = alpha(2);
                           alpha0 = alpha(1);
                       else
                           alpha1 = alpha(1);
                           alpha0 = 0;
                       end
%                        disp('checking av2')
%                        %should evaluate to 0
%                        f = etah - etaw
%                        f = f.subsF([a],[alpha1*b + alpha0])
%                        simplify(f.f)
%                        
                 %      size(envfs)
                                         
                       [envfs,  envds] = solveQuadQuad1(obj, etah,  x, y, a, b, alpha0, alpha1, mh, qh, mw, ix(i), jx(i), etaR, etaV, lV, etaE, lE, envfs,  envds);
                 %      size(envfs)
                 %return
                    end
                end
         %       for i0 = i00+1:size(envfs,2)
         %     envfs(i0)     
         %       end
         %       for ienv = 1:size(envfs,2)
         %         simplifyFraction(envfs(ienv).f)
         %       end
         envfs
%         envfs.printL
         for i = 1:size(envds,2)
             envds(i).print
         end
            end 
            envfs.printL
        end

        

        function l=positivePsi (obj,p,x,y)
            l = 0;
            for i = 1: obj.d.polygon.nv
               if (subs(p,[x,y],[obj.d.polygon.vx(i),obj.d.polygon.vy(i)]) < 0)
                   return
               end
            end
            l = 1;
        end


                %currently edited to take all pairs - to be fixed
        %feasible pairs
        % (ix,jx) feasible pair
        % (vix,vjx) 0 : from V, 1 : from E
        % ixd,jxd : if from E - gives region no - 1,2,3
        function [ix,jx,vix, vjx, ixd, jxd] = feasiblePairs (obj, etaR,a,b)
            n = 0;
            x1 = sym('x1');
            y1 = sym('y1');
            for i = 1:size(etaR,1)
                vj1 = obj.d.E(i,1);  
                vertex = [obj.d.getVertex(vj1)];
                vjx1 = vertex(1);
                vj2 = obj.d.E(i,2);  
                vertex = [obj.d.getVertex(vj2)];
                vjx2 = vertex(1);
                if (vertex(1) < vjx1)
                    vjx2 = vjx1;
                    vjx1 = vertex(1);
                end
                eq1 = y1 - obj.d.mE(i)*x1 - obj.d.cE(i);
                for j = 1:size(obj.d.V,2)
                    vertex = [obj.d.getVertex(obj.d.V(j))];

                    c = vertex(2) + vertex(1) * (1/obj.d.mE(i));
                    eq2 = y1 + (1/obj.d.mE(i))*x1 - c;
                    % eq1
                    % eq2
                    proj = solve(eq1,eq2);
                    s = proj.x1;
                    subs(eq1,proj);
                    subs(eq2,proj);
                    %if ((vjx1<=s)&(s<=vjx2))
                      n = n + 1;
                      ix(n) = i;
                      jx(n) = j;
                      vix(n) = 1;
                      vjx(n) = 0;
                      ixd(n) = 2;
                      jxd(n) = 0;
                    %end  
                    %if (s <= vjx2)
                      n = n + 1;
                      ix(n) = i;
                      jx(n) = j;
                      vix(n) = 1;
                      vjx(n) = 0;
                      ixd(n) = 1;
                      jxd(n) = 0;
                    %end  
                    %if (s >= vjx1)
                      n = n + 1;
                      ix(n) = i;
                      jx(n) = j;
                      vix(n) = 1;
                      vjx(n) = 0;
                      ixd(n) = 3;
                      jxd(n) = 0;
                    %end
                    
                end
            end
            
            for i = 1:size(etaR,1)
                for j = i+1:size(etaR,1)
                      for i1 = 1:3
                          for j1=1:3
                              %if (i1==1 & j1==1)
                              %    continue
                              %end
                              %if (i1==3 & j1==3)
                              %    continue
                              %end
                            n = n + 1;
                            ix(n) = i;
                            jx(n) = j;
                            vix(n) = 1;
                            vjx(n) = 1;
                            ixd(n) = i1;
                            jxd(n) = j1;
                          end
                      end
                end
            end
            
            %return
            % vertex vertex cases
            %size(obj.d.V,2)
            
            for i = 1:size(obj.d.V,2)
              for j = i+1:size(obj.d.V,2)
                 n = n + 1;
                 ix(n) = i;
                 jx(n) = j;
                 vix(n) = 0;
                 vjx(n) = 0;
                 ixd(n) = 0;
                 jxd(n) = 0;
                    
              end
            end
            

        end 
        
        % combine conditions for etaR and etaE
        function [etaV, etaE, etaR] = getEtaFunctions (obj,x,y,a,b)
            eta = obj.f - symbolicFunction(a*x+b*y);
            
            % Eta for Edges
            for i = 1:obj.d.nE
                xv1 = obj.d.polygon.vx(obj.d.E(i,1));
                yv1 = obj.d.polygon.vy(obj.d.E(i,1));
                %etaE(i,1) = eta.subsVarsPartial([x,y],[xv1,yv1]);
            
                xv2 = obj.d.polygon.vx(obj.d.E(i,2));
                yv2 = obj.d.polygon.vy(obj.d.E(i,2));
                %etaE(i,3) = eta.subsVarsPartial([x,y],[xv2,yv2]);
                % check f'(xv1,yv1) < f1(xv2,yv2)
                if xv1 < xv2
                    etaE(i,1) = eta.subsVarsPartial([x,y],[xv1,yv1]);
                    etaE(i,3) = eta.subsVarsPartial([x,y],[xv2,yv2]);
                else
                    etaE(i,3) = eta.subsVarsPartial([x,y],[xv1,yv1]);
                    etaE(i,1) = eta.subsVarsPartial([x,y],[xv2,yv2]);
                end
               % etaE(i,1)
               % etaE(i,3)

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
               
                %disp("1")
                %expand(etaE(i,2).f)

                xp = df.solve(x);
                %(a+obj.d.mE(i)*b - obj.d.cE(i))/(2*obj.d.mE(i))
                %etaT.print
                %etaT.f
                %subs(etaT.f,[x],[xp])
                etaE(i,2) = etaT.subsVarsPartial([x],[xp]);
               
                %disp("2")
                %etaE(i,2).f
               % etaE(i,2)
               % etaE(i,2).print
                
                %obj.f
                %f1 = obj.f.subsVarsPartial([y],[edgey]);
                %f1.f
                %df1 = f1.dfdx(x);
                %df1
                
                etaR(i,1) = symbolicFunction(a+obj.d.mE(i)*b);
                %df1 = symbolicFunction(a+obj.d.mE(i)*b)

                %b1 = df1.subsVarsPartial ([x,y],[xv1,yv1]);
                %b2 = df1.subsVarsPartial ([x,y],[xv2,yv2]);
                %if (double(b1.f) < double(b2.f))
                if xv1 < xv2
                %  etaR(i,2) = b1;  %subs(df1,[x,y],[xv1,yv1]);
                %  etaR(i,3) = b2; %subs(df1,[x,y],[xv2,yv2]);
                else
                %  etaR(i,3) = b1;  %subs(df1,[x,y],[xv1,yv1]);
                %  etaR(i,2) = b2; %subs(df1,[x,y],[xv2,yv2]);
                  %t = etaE(i,1);
                  %etaE(i,1) = etaE(i,3);
                  %etaE(i,3) = t;
                end
                
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
            %disp("in conjugate")
            %size(obj.envf)
            %obj.conjf = sym(zeros(size(obj.envf,2),2*obj.envd(1).nv))
            obj.conjfia(1) = 1;  
            %for i=1:size(obj.envf,2)
            for i=1:size(obj.envelope,2)
              %if i > 1
              %    obj.conjfia(i+1) = size(obj.conjf,2)+1
              %    continue
              %end
            %  i
            %  obj.print
              obj = obj.conjugateFunction(i);
              obj.conjfia(i+1) = size(obj.conjugates,2)+1;
             % return
              %conjd = obj.envd(i).conjugate;
          end
          %obj.conjfia
        end

        function obj = conjugateFunction (obj,i)
            %disp("in conjugateF")
            vars = obj.f.getVars;
            s1 = sym('s_1');
            s2 = sym('s_2');
            dualVars = [s1,s2];
            obj.envelope(i).d = obj.envelope(i).d.simplify;
            %disp('check order')
            %obj.envelope(i).d = obj.envelope(i).d.poly2order;
            %obj.envelope(i).d.print
            if obj.envelope(i).d.nv == size(obj.envelope(i).d.ineqs,2)
                obj.envelope(i).d = obj.envelope(i).d.poly2order;
            else
                obj.envelope(i).d = obj.envelope(i).d.poly2orderUnbounded;
            end
            % disp('inf edge')
            % obj.envelope(i).d.print
            edgeNo = obj.envelope(i).d.getEdgeNosInf(obj.envelope(i).d.vars);
            %d.print
            obj.envelope(i).d.ineqs(edgeNo) = obj.envelope(i).d.ineqs;
            %d.print
%             disp('ord')
             %obj.envelope(i).d.print
% return
            %if obj.envExpr(i).type == 1
            %if obj.envelopeExpr(i).type == 1    
            %if obj.type(i) == 1    
            if ~obj.envelope(i).f.isLinear     
              %disp('type 1')
              t = sym('t');
              % % [x, y, const]
              % %psi2 = psi2(1)*x1+psi22*x2+psi20
              % % hence use index 3 for psi20 terms
              % cpsi2 = obj.envelopeExpr(i).vpsi2.getLinearCoeffs (vars);
              % cpsi1 = obj.envelopeExpr(i).vpsi1.getLinearCoeffs (vars);
              % cpsi0 = obj.envelopeExpr(i).vpsi0.getLinearCoeffs (vars);
              % 
              % vs1 = s1 - (2*cpsi1(1)*t - cpsi2(1)*t^2 + cpsi0(1));
              % vs2 = s2 - (2*cpsi1(2)*t - cpsi2(2)*t^2 + cpsi0(2));
              % 
              % vt = solve (cpsi2(2)*vs1 - cpsi2(1)*vs2, t );    %  cpsi2(2)*vs1 - cpsi2(1)*vs2  cancels t^2 term - then solve for t 
              % crs = subs(vs1,t, vt)    % crs = 0  equation of parabola 
              
              
              num = obj.envelope(i).f.getNum
              den = obj.envelope(i).f.getDen
              [r,q] = polynomialReduce(num,den)
              size(r)
              
              r1 = factor(r)
              size(r1)
              
              r2 = r1(2)
              
              
              if size(r1,2) == 3
                den = symbolicFunction((1/r1(1)) * den)
              else
                den = symbolicFunction(den)  
              end
              q = symbolicFunction(q)
              r = symbolicFunction(r2)
              cpsi2 = den.getLinearCoeffs (vars)
              cpsi0 = q.getLinearCoeffs (vars)
              cpsi1 = r.getLinearCoeffs (vars)
              
              vs1 = s1 - (2*cpsi1(1)*t - cpsi2(1)*t^2 + cpsi0(1));
              vs2 = s2 - (2*cpsi1(2)*t - cpsi2(2)*t^2 + cpsi0(2));

              vt = solve (cpsi2(2)*vs1 - cpsi2(1)*vs2, t );    %  cpsi2(2)*vs1 - cpsi2(1)*vs2  cancels t^2 term - then solve for t 
              crs = simplifyFraction(subs(vs1,t, vt))    % crs = 0  equation of parabola 

              %simplify(crs - crs2)
              
             % pause;
%%%%%%%%%%%%%
% Checking parabola
%              crs2 = cpsi2(1)^2 * cpsi2(2)^2 * s1^2 -2 * cpsi2(1)^3*cpsi2(2)*s1*s2 + cpsi2(1)^4*s2^2
%%%%%%%%%%%%%
              
              NCV = obj.envelope(i).d.getNormalConeVertex(s1, s2);
              %disp('ncv')
              edgeNo = obj.envelope(i).d.getEdgeNos(vars);
              % for i = 1:obj.envelope(i).d.nv
              %     edgeNo(i) = i;
              % end
              NCE = obj.envelope(i).d.getNormalConeEdge(s1, s2);

  
            % check eta1(1)=eta2(1)=0  page 68/136
              [subdV,undV] =  obj.envelope(i).getSubdiffVertexT1 (NCV, dualVars);
             % disp('subd')
              [subdE,unR] = obj.envelope(i).getSubdiffVertexT2 (NCE, dualVars);
            %  disp('subdt2')
              %[subdE, unR, crs] = obj.getSubDiffEdgeT1(i, subdE, edgeNo, undV, crs, dualVars);
              [subdE, unR, crs] = obj.envelope(i).getSubDiffEdgeT1(subdE, edgeNo, undV, crs, dualVars)
              
              subdV = obj.envelope(i).getSubDiffVertexSpT1(subdV, undV, crs)
              %disp('type1');
              expr = obj.envelope(i).conjugateExprVerticesT1 (dualVars, undV )
              %expr = obj.envelope(i).conjugateExprEdgesT1Poly (dualVars, edgeNo, cpsi0, cpsi1, cpsi2, expr );
              %disp('Checking conjugate expr')
              unR
              %expr
              disp('temp code')
              %
              % pause
              %unR(1) = 0
              
              %pause
              for j = 1:obj.envelope(i).d.nv

                  if unR(j)
                      % rt = region(subdE(j,:), dualVars)
                      % rt.print
                      % pause
                      ex1 = conjugateExpr(obj.envelope(i).d.ineqs(j).f,obj.envelope(i).f.f,obj.envelope(i).d.vars(1),obj.envelope(i).d.vars(2));
                      expr = [expr, subs(ex1,obj.envelope(i).d.vars,dualVars)];
                  else
                      expr = [expr, inf];
                  end
              end
              %obj.envelope(i).d.conjugateExprEdgesT1Poly2 (obj.envelope(i).f, obj.envelope(i).d.vars)   % test this and replace
            %elseif obj.envelopeExpr(i).type == 3   
            else  %if obj.type(i) == 3   
              % cpsi0 = obj.envelopeExpr(i).vpsi0.getLinearCoeffs (vars)  ;
              % vs1 = cpsi0(1);
              % vs2 = cpsi0(2);
              NCV = obj.envelope(i).d.getNormalConeVertex(s1, s2);
              [subdV,undV] =  obj.envelope(i).getSubdiffVertexT1 (NCV, dualVars);
              %disp('type3')
              expr = obj.envelope(i).conjugateExprVerticesT1 (dualVars, undV );
              for j = 1:obj.envelope(i).d.nv
                unR(j) = false;
              end
              % corollary 4.15
              % corollary 4.20
              % corollary 4.26 - vertex
            end
           conjf=symbolicFunction.empty;
           conjd=region.empty;
           
            for j = 1:obj.envelope(i).d.nv
                
                if undV(j)
                  if j == obj.envelope(i).d.nv
                    e0 = 1;
                  else    
                    e0 = j+1;
                  end  
                  conjf = [conjf,expr(j)];
                  conjd = [conjd, region([subdE(e0,1:2),-subdE(e0,3)], dualVars)];
                   %%%%%%%%%%%%%
                  r = subdV(e0,:);
                  r(2) = -r(2);
                   conjf = [conjf,expr(j)];
                   conjd = [conjd, region(r, dualVars)];
                   if j == 1
                    e0 = obj.envelope(i).d.nv;
                   else    
                    e0 = j-1;
                   end  
                   r = subdV(e0,:);
                   r(1) = -r(1);
                   conjf = [conjf,expr(j)];
                   conjd = [conjd, region(r, dualVars)];
                else
                  conjf = [conjf,expr(j)];  
                  conjd = [conjd, region(subdV(j,:), dualVars)];
                end
                %conjd(j).print
                %expr(j)
                
            end
           
            for j = 1:obj.envelope(i).d.nv
                if (unR(j))
                    r = region(subdE(j,:), dualVars);
                    r = r.simplifyUnboundedRegion;
                    if isempty(r)
                        continue
                    end
                conjf = [conjf,expr(obj.envelope(i).d.nv+j)];
                conjd = [conjd, region(subdE(j,:), dualVars)];
                end
            end
            conjugates = functionNDomain.empty;
            for i = 1:size(conjf,2)
                %disp('in conjugate')
                %conjd(i).print
                %conjf(i).print
                conjugates = [conjugates,functionNDomain([conjf(i)],conjd(i))];
            end
            %disp('b4 merge')
            % HISTORY: this leftover debug print crashed on a rational
            % (non-polynomial-denominator) conjugate expression --
            % symbolicFunction.print's isPolynomial/degreeDen chain calls
            % polynomialDegree on the denominator unconditionally, which
            % errors ("Polynomial expression expected") for a genuinely
            % rational expression rather than just reporting isPolynomial
            % false. Purely diagnostic output with no effect on conjugates
            % itself (testFractional, whose conjugate is exactly such a
            % rational expression) -- commented out like its sibling debug
            % lines above.
            %conjugates.printL
            conjugates = conjugates.mergeL;

            for i = 1:size(conjugates,2)
                obj.conjugates = [obj.conjugates,conjugates(i)];
            end
%           
            
            return 

        end

       

        


        % function [subdE, unR, crs] = getSubDiffEdgeT1(obj, i, subdE, edgeNo, unDV, crs, dualvars)
        %     %subdE = sym(zeros(obj.envd(i).nv,4));
        %         vars =  obj.envelope(i).d.vars;
        %         drx1 = obj.envelope(i).f.dfdx(vars(1));
        %         drx2 = obj.envelope(i).f.dfdx(vars(2));
        % 
        %     unR = zeros(obj.envelope(i).d.nv,1);
        %     for j = 1:obj.envelope(i).d.nv-1
        %         if unDV(j)
        %             unR(j) = false;
        %             continue
        %         end
        %         if unDV(j+1)
        %             continue
        %         end
        % 
        %         unR(j) = true;
        %        %obj.envelope(i).d.vx(j),obj.envelope(i).d.vx(j+1)
        %        % Put checks for div by 0
        %         dv11 = subs(drx1.f,vars,[obj.envelope(i).d.vx(j),obj.envelope(i).d.vx(j+1)]);
        %         dv12 = subs(drx2.f,vars,[obj.envelope(i).d.vx(j),obj.envelope(i).d.vx(j+1)]);
        %         dv21 = subs(drx1.f,vars,[obj.envelope(i).d.vy(j),obj.envelope(i).d.vy(j+1)]);
        %         dv22 = subs(drx2.f,vars,[obj.envelope(i).d.vy(j),obj.envelope(i).d.vy(j+1)]);
        % 
        %         if subs(crs,dualvars,[(dv11+dv21)/2,(dv12+dv22)/2]) < 0
        %           crs = -crs;
        %         end
        %         subdE(j,3) = crs;
        % 
        %     end    
        %     j = obj.envelope(i).d.nv;
        %     if unDV(obj.envelope(i).d.nv)
        %         unR(j) = false;
        %       return
        %     end
        %     if unDV(1) 
        %       return
        %     end
        %     unR(j) = true;
        % 
        % 
        % 
        %         dv11 = subs(drx1.f,vars,[obj.envelope(i).d.vx(j),obj.envelope(i).d.vx(1)]);
        %         dv12 = subs(drx2.f,vars,[obj.envelope(i).d.vx(j),obj.envelope(i).d.vx(1)]);
        %         dv21 = subs(drx1.f,vars,[obj.envelope(i).d.vy(j),obj.envelope(i).d.vy(1)]);
        %         dv22 = subs(drx2.f,vars,[obj.envelope(i).d.vy(j),obj.envelope(i).d.vy(1)]);
        % 
        %         if subs(crs,dualvars,[(dv11+dv21)/2,(dv12+dv22)/2]) < 0
        %           crs = -crs;
        %         end
        %         subdE(j,3) = crs;
        % end

        % function [subdV,undV] = getSubdiffVertexT1 (obj, i, NCV, dualVars)
        % 
        %     [subdV,undV] =  obj.envelope(i).getSubdiffVertexT1 (NCV, dualVars);
        %     return
        % 
        % end
        % 
        % function [subdV,undV] = getSubdiffVertexT2 (obj, i, NCV, dualVars)
        %     [subdV,undV] = obj.envelope(i).getSubdiffVertexT2 (NCV, dualVars);
        %     return
        % 
        % end
        % 
        % 
        % function [edgeNo] = getEdgeNos(obj, i)
        %   vars = obj.f.getVars;  
        %   edgeNo = obj.envelope(i).d.getEdgeNos(vars)  
        % end

    

    
     function obj = maximumConjugate(obj)
       %obj.conjfia(1):obj.conjfia(2)-1 
       for k = obj.conjfia(1):obj.conjfia(2)-1 
         obj.maxConjugate(k) = obj.conjugates(k);
        end
          % obj.conjugates(obj.conjfia(1):obj.conjfia(2)-1 ).printM2
          % obj.conjugates(obj.conjfia(2):obj.conjfia(3)-1 ).printM2
          % disp("mConjM")
          %     obj.maxConjugate.printM
          for i = 2:size(obj.envelope,2) %size(obj.conjfia,2)-1
             % obj.conjfia(i):obj.conjfia(i+1)-1
              %obj.maxConjugate.printL
              obj.maxConjugate = obj.maxConjugate * obj.conjugates(obj.conjfia(i):obj.conjfia(i+1)-1);
              % to be put back
              %disp("mConj")
              %obj.maxConjugate.printM2
              %obj.maxConjugate.printL
              %return
              %obj.maxConjugate.printM2
              obj.maxConjugate = obj.maxConjugate.maximumP(true);
              %obj.maxConjugate.printM
              %disp("mConjM")
              %obj.maxConjugate.printL
              %obj.maxConjugate.printM
          end
     end

     function biconjugateP(obj)
         
              
           for j=1:size(obj.envelope,2) 
             if (size(obj.conjfia,1) > 0)
                 bic = functionNDomain.empty ;
                bic = obj.conjugates(obj.conjfia(j):obj.conjfia(j+1)-1).conjugateOfPiecePoly;
                bic = bic.addEq;
                bic.printL
             end
           end
           
           bc = obj.maxConjugate.conjugateOfPiecePoly
           bc.printL
           bc = bc.addEq;
           bc.printL

       

         end

    end

    
    

 
end