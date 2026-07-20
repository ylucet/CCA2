classdef symbolicFunction
    properties (Access=private)
        % to be re written restricting functions to be in terms of defined
        % variables only
        nv; 
        vars; 
    end    
    properties  
        f = sym('f') ;
    end
% 57 methods
   
    methods  % init 
        function obj = symbolicFunction(num0, den0)
 
            
            if nargin == 0
              num=0;
              den=1;
              obj.f= num / den;
            elseif nargin == 1
               num=num0;  
               den=1;  
            elseif nargin == 2
              num=num0;
              den=den0;
            end
            if nargin ~= 0
            obj.f= num / den;

            % temp change 19/4/24
            % class(num0)
            % isreal(num0)
            if num0 == 0 | isreal(num0) 
                obj.vars = [];
                obj.nv = 0;
                return
            end

            obj.vars = symvar(obj.f);
            obj.nv = size(obj.vars,2);
            end
        end
        
        function f = getF(obj)
            f = obj.f;
        end  
        
        function num = getNum(obj)
            if isequal(class(obj.f),'double')
              num = obj.f;  
            else
                %obj.f
              [num,den] = numden(obj.f);
            end
        end   
        function den = getDen(obj)
            if isequal(class(obj.f),'double')
              den = 1;  
            else
              [num,den] = numden(obj.f);
            end
        end   
        
    end
    methods % display
        function print(obj)
            
          if obj.isPolynomial
          [coef,terms] = coeffs(obj.f);
         
          for i=1:length(terms)
                fprintf(num2str(double(coef(i))));
              if terms(i) ~= 1
                fprintf(char(terms(i)));
              end
              if i == length(terms)
                  break;
              end
              fprintf(" + ");
          end
          fprintf("\n")
          else
              fprintf(char(simplify(obj.f))); 
          fprintf("\n")
          
          end
        end

        function printLatexWB(obj)
          %  disp('obj.isPol')
           % obj.isPolynomial
          if obj.isPolynomial
          [coef,terms] = coeffs(obj.f);
          
          for i=1:length(terms)
              if abs(double(coef(i))-1) > 1.0d-8
                if i == 1 & double(coef(i)) < 0
                      fprintf(" - ")
                  end
                    
                [n,d] = numden(coef(i));
                %isAlways(d == 1)
                if  isAlways(d == 1)
                  fprintf(num2str(abs(double(coef(i)))));
                  
                else
                  fprintf("\\frac{" + num2str(abs(double(n)))+"}{"+ num2str(double(d))+"}");  
                end
              end
              if terms(i) ~= 1
                fprintf(char(terms(i)));
                
              end
              if i == length(terms)
                  break;
              end
              if double(coef(i+1)) > 0
                fprintf(" + ");
              else  
                fprintf(" - ");  
              end
          end
          else
              [n,d] = numden(obj.f);
              fprintf("\\frac{" + char(n)+"}{"+ char(d)+"}");  
          
          end
        end

        function printLatex(obj)
     
          fprintf("\\[");
          obj.printLatexWB;
          fprintf("\\]\n")
          
        end

        function fprint(obj, uNo)
          fprintf(uNo, char(simplify(obj.f))); 
          fprintf(uNo, "\n")
        end

        
        function plot3d(obj, limits)
            if limits(1) == limits(2)
                limits(1) = limits(1)-30;
                limits(2) = limits(2)+30;
            end
            if limits(3) == limits(4)
                limits(3) = limits(3)-30;
                limits(4) = limits(4)+30;
            end
            ezsurf(obj.f, limits);
        end

        function plot(obj, vars, limits)
            colors = ['b', 'r', 'g', 'm', 'c', 'y', 'k'];
            polyvars = obj.vars;
            if  size(vars,2) == size(polyvars,2)
                vx = linspace(limits(1), limits(2), 100);
                fy = solve(obj.f==0,vars(2));
                vy = subs(fy,polyvars(1),vx);
                %[vx,vy] = getPoints (obj, polyvars, limits)
                plot(vx,vy);
                %fill(vx,vy,colors(1+mod(obj.getGlobalParameter,7)),'FaceAlpha',0.5);
            elseif ismember (vars(1), polyvars)
                fx = solve(obj.f==0,vars(1));
                vy = zeros(1,100);
                vx = double(fx) + vy;
                vy = linspace(limits(1), limits(2), 100);
                plot(vx,vy);
            elseif ismember (vars(2), polyvars)
                fy = solve(obj.f==0,vars(2));
                vx = zeros(1,100);
                vy = double(fy) + vx;
                vx = linspace(limits(1), limits(2), 100);
                plot(vx,vy);
            else
            end
            %obj.setGlobalParameter; 
        end

        function plotIneq (obj, limits)
            ezsurf(obj.f <=0, limits)
        end 

        
        function printIneq(obj)
            fprintf(char(obj.f));
            fprintf(" <= 0 \n")
        end

        function printIneqLatex(obj)
            fprintf("\\[")
            obj.printLatexWB
            fprintf(" \\le 0")
            fprintf("\\] ")
            fprintf("\n")
        end

        function printIneqM(obj)
            fprintf(char(obj.f));
            fprintf(" <= 0 ")
        end

        function printL (l, first, last)
            %l
            %nargin
            %size(l,1)
            if nargin == 1
            
            for i = 1: size(l,1)
                for j = 1: size(l,2)
                    
                    % if isinf(l(i,j)) 
                    %     disp('Inf')
                    % else
                    %class (l(i,j).f)
                    %if isSymType(l(i,j).f, 'sym')
                    l(i,j).f = simplifyFraction(l(i,j).f);
                    
                    l(i,j).print;
                    %end
                end
            end
            else
                for i = 1: size(l,1)
                for j = first: last
                    l(i,j).f = simplifyFraction(l(i,j).f);
                    l(i,j).print;
                end
            end
            
            end
        end

        function printLLatex0 (l, first, last)

            if nargin == 1
            
            for i = 1: size(l,1)
                for j = 1: size(l,2)
                    l(i,j).printLatex;
                end
            end
            else
                for i = 1: size(l,1)
                for j = first: last
                    l(i,j).printLatex;
                end
            end
            
            end
        end

        function printLIneq (l)
            
            for i = 1: size(l,1)
                for j = 1: size(l,2)
                    l(i,j).printIneq;
                end
            end
        end

        function printLIneqLatex (l)
            
            for i = 1: size(l,1)
                for j = 1: size(l,2)
                    l(i,j).printIneqLatex;
                end
            end
        end

         function printLIneqM (l)
            
             fprintf("{")
                for j = 1: size(l,2)-1
                    l(j).printIneqM;
                    fprintf(",");
                end
                l(size(l,2)).printIneqM;
                    fprintf("}");
            
        end
        function fprintLIneq (l,uNo)
            
            for i = 1: size(l,1)
                for j = 1: size(l,2)
                    l(i,j).fprint(uNo);
                end
            end
        end

        function plotLIneq (l, vars,limit)
            % change this
            %l2 = [];
            for i = 1: size(l,1)
                for j = 1: size(l,2)
                    %l(i,j)
                     
                    l(i,j).plot(vars,limit);
                    hold on;
             %       l2 = [l2,l(i,j).f<=0]
                end
            end
            %plot(intersect(l2))
        end

        function plot2 (l)
            figure;
            plot(intersect(l))
        end
    
    end
    
    methods % operations
        function f = plus(obj1,obj2)
            f = symbolicFunction(obj1.f+obj2.f);
        end
        function f = minus(obj1,obj2)
            f = symbolicFunction(obj1.f-obj2.f);
        end

        % to be disabled
        function f = unaryminus(obj)
            f = symbolicFunction(-obj.f);
        end

        function f = uminus(obj)
            f = symbolicFunction(-obj.f);
        end

        
        function f = mtimes(f1,f2)
            f = symbolicFunction(f1.f*f2.f);
        end

    end

    methods % derivatives pass variables here

         function f = dfdx (obj,x)
             
            f = symbolicFunction(simplify(diff(obj.f,x)));
         end 

         function f = limit (obj, vars, pt)
             f0 = obj.f;
             for i=1:size(vars,2)
               f0 = limit(f0,vars(i),pt(i));
             end
             f = symbolicFunction(f0);
         end

         function g = gradient(obj, vars)
           
           for i = 1:size(vars,2)
             g(i) = obj.dfdx(vars(i));
            % g(i).print
           end 
         end

         % this works only for bivariate - put checks in place and change
         % name
         
         function f = tangentOfSlope (obj, m)
             dx = obj.dfdx(obj.vars(1))
             dy = obj.dfdx(obj.vars(2))
             
             f = symbolicFunction(m * dy.f - dx.f);
             
             
         end


         function f = tangent (obj, x, y)
             
         
             dx = obj.dfdx(obj.vars(1));
             dy = obj.dfdx(obj.vars(2));
             if isAlways(dy.subsF(obj.vars,[x,y]).f==0)
                 f = symbolicFunction(obj.vars(1)-x);
             else
               m =  - dx.subsF(obj.vars,[x,y]).f / dy.subsF(obj.vars,[x,y]).f;
               c = y - m*x;
               f = symbolicFunction(obj.vars(2) - m * obj.vars(1) -c);
             end
             
         end

    end

    methods % inquiry

        function l = isPolynomial(obj)
            % symType(obj.f)
            % if symType(obj.f) == 'expression'
            %     l = false;
            % else
            %     disp('hhh')
            %     obj.degreeDen
              l = obj.degreeDen == 0;
            % end
        end
        
        function l = isQuad(obj)
        
          if (obj.degreeDen ~= 0)
              l = false;
              return;
          end
        
          if (obj.degreeNum == 2)
              l = true;
          else
              l = false;
          end
        
        end

        function l = isParabolic(obj)
            l = false;
            if ~ obj.isQuad
                return;
            end
            c = coeffs(obj.f, obj.vars)
        end
        
        function l = isLinear(obj)
         
          if (obj.degreeDen ~= 0)
              l = false;
              return;
          end
          if (obj.degreeNum == 1)
              l = true;
          else
              l = false;
          end
        end

        % fix this
        function l = isConst(obj)

            if size(obj.vars,2) > 0
                l = false;
                return
            end
            cn = coeffs(obj.f,obj.vars);
            l = all(cn(2:end)==0);
        end
    end
    
    methods
        function vars = getVars(obj)
            vars = obj.vars;
        end

        % [x, y, const]
        function c = getLinearCoeffs (obj,vars)
           %obj.f
           %vars
           if obj.isZero
               c(1) = 0;
               c(2) = 0;
               c(3) = 0;
               return
           end
           cvars = obj.getVars;
           if (isempty(cvars))
               
               c(1) = 0;
               c(2) = 0;
               c(3) = coeffs(obj.f,vars(1));
               return;
           end
           if size(cvars,2) == size(vars,2)
              ct = coeffs(obj.f,vars(1));
              c(1) = ct(2) ;
              ct = coeffs(ct(1),vars(2));

              if (size(ct,2) == 2)
                c(3) = ct(1);
              else
                  c(3) = 0;
              end
              ct = coeffs(obj.f,vars(2));
              c(2) = ct(2);
           elseif size(cvars,2) == 1
              if (cvars(1) == vars(1))
                ct = coeffs(obj.f,vars(1));
                if (size(ct,2) == 2)
                  c(3) = ct(1);
                  c(1) = ct(2);
                  c(2) = 0;
                else
                  c(2) = 0;
                  c(3) = 0;
                  c(1) = ct(1);
                end
               else
                ct = coeffs(obj.f,vars(2));
                if (size(ct,2) == 2)
                  c(3) = ct(1);
                  c(1) = 0;
                  c(2) = ct(2);
                else
                  c(2) = ct(1);
                  c(3) = 0;
                  c(1) = 0;
                end
               
              end
           end

        end

        function [l,c] = quadterm (obj, x)
            
            [qc,qt] = coeffs(obj.f,x);
            l = false;
            
            for i = 1:size(qt,2)
                if isAlways (qt(i) == x^2)
                    l = true;
                    c = qc(i);
                    return
                end
            end
        end
        

        % not for rational functions
        function obj = normalize1 (obj)
            if obj.getDen ~= 1
                %disp('Rational in normalize1')
            end
            %obj.f
            %obj.vars
            
          c = coeffs(obj.f,obj.vars);
          obj = symbolicFunction(simplify((1/abs(c(end)))*obj.f));
          
        end

        function obj = normalize (obj,vars)
          c = obj.getLinearCoeffs (vars);
             
          if (c(2) == 0)
            obj.f = (1/c(1)) * obj.f;
          else
             obj.f = (1/c(2)) * obj.f;
          end
        end

        function f = subsF (obj,vars,vals)
            
            den = simplifyFraction(subs(obj.getDen, vars, vals));
            num = simplifyFraction(subs(obj.getNum, vars, vals));

            if (isAlways(den == 0))
                if (num == 0)
                  f = symbolicFunction(sym(nan),1);  
                elseif (isAlways(num > 0) )   
                  f = symbolicFunction(sym(intmax),1);
                else
                  f = symbolicFunction(sym(-intmax),1);
                end  
                return;
            end
            % vars
            % vals
            %subs(obj.f, vars, vals)
            f = symbolicFunction(num,den);
            % if (isAlways(subs(obj.getDen, vars, vals) == 0))
            %     if (subs(obj.getNum, vars, vals) == 0)
            %       f = symbolicFunction(sym(nan),1);  
            %     elseif (isAlways(subs(obj.getNum, vars, vals) > 0) )   
            %       f = symbolicFunction(sym(intmax),1);
            %     else
            %       f = symbolicFunction(sym(-intmax),1);
            %     end  
            %     return;
            % end
            % % vars
            % % vals
            % %subs(obj.f, vars, vals)
            % f = symbolicFunction(subs(obj.f, vars, vals));
        end    

           
        function l = isZero(obj)
            l = false;
            if obj.getDen ~= 1
                return
            end
            
            if obj.nv == 0
                % change to isAlways
                if (abs(double(obj.f)) < 1.0e-6)
                    l = true;
                end
                return
            end
            c = coeffs(obj.f,obj.vars);

            n = size(c,2);
            for i = 1:n
                if isAlways(abs(c(i)) > 0)
                %if (abs(double(c(i))) > 1.0e-6)
                    return
                end
            end
            l = true;
            
        end
        
        function f = subsVarsPartial (obj,vars,varVals)
            f0 = simplify(subs(obj.f, vars, varVals));
            f = symbolicFunction(f0);
            
        end    

        function d = double (obj)
            c = obj.f.coeffs();
            if (size(c) == 0)
                d = 0;
                return
            end
            if (size(c) ~= 1)
                disp("Error in double in symbolicFunction");
                return
            end
            d = c(1);
        end

       
        function f = solve (obj,x)
            f = solve(obj.f,x);
        end    
        
   
        % change to isAlways
        function res = eq(obj1,obj2)
            res = false;
            if (obj1.f==obj2.f)
                res = true;
            end 
        end
        
         function res = le(obj1,obj2)
            res = false;
            if (obj1.f<=obj2.f)
                res = true;
            end 
        end
        function res = ge(obj1,obj2)
            res = false;
            if (obj1.f>=obj2.f)
                res = true;
            end 
        end
        
        function res = lt(obj1,obj2)
            
            if isequal(class(obj2),'double')
                
                res = ltd(obj1,obj2);
            else
            res = false;
            if (obj1.f<obj2.f)
                res = true;
            end 
            end
        end
        
        function res = gt(obj1,obj2)
            if isequal(class(obj2),'double')
                res = gtd(obj1,obj2);
            else
            res = false;
            
            if (obj1.f>obj2.f)
                res = true;
            end 
            end
        end
        
        
        
        %obj2 double
        function res = ltd(obj1,obj2)
            res = false;
            if (obj1.f<obj2)
                res = true;
            end 
        end
        
        function res = gtd(obj1,obj2)
            res = false;
            if (obj1.f>obj2)
                res = true;
            end 
        end
        

        function [vx,vy] = solveF (f2)
            f1x = subs(obj.f, obj.vars,[x,y]);
            f2x = subs(f2.f, obj.vars,[x,y]);
            s = solve ([f1==0,f2==0],[x,y]);
            vx = s.x;
            vy = s.y;
        end
            
        
        function f = removeDenominator2 (obj)
            %f = obj.num;
            %num = f.num
           
            cx = coeffs(obj.getNum,obj.vars);
            cz=[];
            for i = 1:size(cx,2)
              cz = [cz,1/cx(i)];
            end
            mult = 1;
            if (size(cz)>0)
                mult= lcm(cz);
            end
            f = obj.f* abs(mult);
            return

            x = obj.vars(1);
            y = obj.vars(2);
            cy = [];
            cx = coeffs(obj.getNum,x);
            for i = 1:size(cx,2)
                g  = cx(i);
                cy = [cy,coeffs(cx(i),y)];
                
            end
            cz=[];
            for i = 1:size(cy,2)
              cz = [cz,1/cy(i)];
            end
            if (size(cz)>0)
                mult= lcm(cz);
            end
            f = obj.f* abs(mult);
        end

        function l = isNegativeSqr(obj,z)
            l = false;
            z0 = simplify(obj.f);
            c = coeffs(obj.f);
            if c(end) > 0
                return
            end
            z1 = obj.solve(z);
            s = z1(1);
            for i = 2:size(z1,2)
                if (s ~= z1(i)) 
                    return;
                end
               
            end
             l = true;
        end
        
        

            function d = degreeNum(obj)
                if isequal(class(obj.getNum),'double')
                    d = 0;
                else
                  % HISTORY: polynomialDegree errors ("Polynomial expression
                  % expected") on a genuinely non-polynomial numerator (e.g.
                  % one involving a fractional power/sqrt term, as in
                  % testFractional's conjugate expressions) rather than
                  % returning a value -- but degreeNum/degreeDen are boolean
                  % predicates' building blocks (isPolynomial, isQuad), which
                  % already handle "the degree isn't 0/2" correctly; they
                  % just need a degree value to compare against, not a crash.
                  % Inf is unambiguously "not a polynomial of any finite
                  % degree", matching every caller's existing comparison.
                  try
                    d = polynomialDegree(obj.getNum);
                  catch
                    d = Inf;
                  end
                end
            end

            function d = degreeDen(obj)
                if isequal(class(obj.getDen),'double')
                    d = 0;
                else
                  try
                    d = polynomialDegree(obj.getDen);
                  catch
                    d = Inf;
                  end
                end
            end


            
    end 

    %     function [l,obj] = removeSum(obj, lprint)
    %         rm = [];
    %         l = false;
    %         for i = 1:size(obj,2)
    %           o1 = obj(i);
    %           for j = i+1:size(obj,2)
    %             o2 = obj(j) + o1;
    %             % removing = 0 also although only >0 are invalid
    %             if lprint
    %                 disp("o2 b4")
    %                 o2
    %             end
    % 
    %             o2 = simplify(o2.f < 0);
    %             if o2 == symfalse
    %                l = true;
    %               return
    %             end
    %             if lprint
    %                 disp("o2")
    %                 o2
    %             end
    %             for k = 1:size(obj,2)
    %               o = simplify(obj(k).f<0);
    %               if (o2 == o)
    %                   rm = [rm,k];
    %               end
    %             end
    % 
    %           end
    % 
    %         end
    %         if lprint
    %         rm 
    %         end
    %         obj(rm) = [];
    %     end
    % 
    % end

end