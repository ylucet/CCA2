
classdef domain
  properties
    % vertex information  
    %nVertices {mustBeInteger}
    %vx;
    %vy;
    
    % ineq representation
    polygon=region;
    % edge information
    nE=0; % number of convex edges
    E;    % index set 
    % SYMBOLIC, and the declaration is the whole of it. These were plain `[]`, i.e. DOUBLE
    % arrays, so `obj.mE(k) = m` silently converted an exact slope to a double -- and
    % plq_1p.conjugateFunction's nCE == 1 branch builds that piece's entire conjugate out of
    % mE(1) and cE(1).
    %
    % MEASURED 2026-08-17 on the A.4 sub-triangle conv{(1/2,1), (sqrt(30)/6, sqrt(30)/10),
    % (sqrt(30)/12 - sqrt(15)/6 + 5/4, ...)}: the exact slope
    % (sqrt(30)/10 - 1)/(sqrt(30)/6 - 1/2) arrived as the double 0.6, whose `sym` is
    % 5404319552844595/9007199254740992 -- the 2^53 denominator that then multiplied out to
    % 145-digit coefficients in the conjugate cells. The y-intercept, exactly 0, arrived as
    % -9.0557e-72. That is a WRONG VALUE, not only an expensive one: an intercept that should
    % be zero is not.
    %
    % Downstream, region.merge finds a shared facet with `ineqs(i) == -ineqs(j)` and compares
    % quadratics with `~=`; neither can decide when one side is exact and the other is its own
    % double, which is why Step 3's cross-piece maximum stopped merging anything at all.
    % DECISIONS.md 2026-08-17 has the measurements.
    mE = sym.empty(); % slope
    cE = sym.empty(); % y intercept
    ind;   % misplaced parameter - needed for convex envelope - remove later
    % remaining vertices
    nV=0;
    V;    % index set

  end 
% 9 methods
  methods  % testing




    end
 
    methods

        function obj = domainEdge (obj,l,vars)
            obj.polygon = region(l,vars);
            obj = getEdges (obj);;
        end

      function obj = domain(v, x,y)
          
          if nargin > 0
            %obj.nVertices = size(v,1) ;
            %[obj.vx,obj.vy] = poly2cw(v(:,1),v(:,2));

            obj.polygon.nv=size(v,1) ;
            % Fix order to clockwise
            %[obj.polygon.vx,obj.polygon.vy] = v(:,1),v(:,2));
            obj.polygon.vx = v(:,1);
            obj.polygon.vy = v(:,2);
            
            obj.polygon.vars = [x,y];
          else
              return
          end
          %obj.polygon.vx
          %obj.nVertices
          %obj.vx
          %obj.vy
           %obj.polygon.print;
          obj = getEdges (obj);
          % obj.polygon.print;
          %obj.E
          %obj.mE
          %obj.cE
          %obj.V
          obj = getAllEdges (obj, x, y);
          %obj.polygon.print;
      end
      
      function vertex = getVertex(obj,i)
          vertex = [obj.polygon.vx(i), obj.polygon.vy(i)];
      end

      function obj = getEdges (obj)
          for i = 1:obj.polygon.nv
              l(i) = 0;
          end
          for i = 1:obj.polygon.nv-1
            m = obj.slope (i,i+1);
             if (m > 0 & m < inf)
                obj.nE = obj.nE+1;
                obj.E(obj.nE,1) =  i;
                obj.E(obj.nE,2) =  i+1;
                obj.mE(obj.nE) = m;
                obj.cE(obj.nE) = yIntercept (obj,i,m);
                l(i)=1;
                l(i+1)=1;
            end
          end 
          
          m = obj.slope (obj.polygon.nv,1);
          if (m > 0 & m < inf)
                
            obj.nE = obj.nE+1;
            obj.E(obj.nE,1) =  obj.polygon.nv ;
            obj.E(obj.nE,2) =  1;
            obj.mE(obj.nE) = m;
            obj.cE(obj.nE) = yIntercept (obj,1,m);
            l(1)=1;
            l(obj.polygon.nv) = 1;
          end

          
          for i = 1:obj.polygon.nv
              if (l(i) == 1) 
                  continue
              end    
              obj.nV = obj.nV+1;    
              obj.V(obj.nV) = i;
          end    
          
      end

      function obj = getAllEdges (obj, x, y)
        cx = mean(obj.polygon.vx);
        cy =  mean(obj.polygon.vy);
        for i = 1:obj.polygon.nv
          if (i == obj.polygon.nv) 
            m = obj.slope (i,1);
          else
            m = obj.slope (i,i+1);
          end
          if (m == inf | m == -inf)
            obj.polygon.ineqs(i) = symbolicFunction(x  - obj.polygon.vx(i));
          else
            obj.polygon.ineqs(i) = symbolicFunction(y - m*x - yIntercept (obj,i,m));
          end
          if obj.polygon.ineqs(i).subsVarsPartial([x,y],[cx,cy]) > 0
              %obj.ineqs(i) = -obj.ineqs(i)
              obj.polygon.ineqs(i) = obj.polygon.ineqs(i).unaryminus;
          end
        end
           
         
      end
      
      function m = slope (obj,i,j)
          m = (obj.polygon.vy(i)-obj.polygon.vy(j))/(obj.polygon.vx(i)-obj.polygon.vx(j));
          % change to
          %m = obj.polygon.slope(i,j);
      end
      
      function c = yIntercept (obj,i,m)
          c = obj.polygon.vy(i)-m*obj.polygon.vx(i);   
          % change to c = obj.yIntercept (i,m)
      end

      function print(obj)
        %disp(["nVertices = ", num2str(obj.polygon.nv)]);
        %fprintf("vx =  ")
        %fprintf("%d  ", obj.polygon.vx);
        %fprintf("\n")
        %fprintf("vy =  ")
        %fprintf("%d  ", obj.polygon.vy);
        %fprintf("\n")
        fprintf("Polygon Ineqs <= 0 \n")
        obj.polygon.print
        disp(["Number of Convex edges = ", num2str(obj.nE)])
        disp("Edges joining vertex numbers")
        disp(obj.E)
        fprintf("Slopes =  ")
        for iP = 1:numel(obj.mE), fprintf("%s  ", string(obj.mE(iP))); end
        fprintf("\n")
        fprintf("y-intercepts =  ")
        for iP = 1:numel(obj.cE), fprintf("%s  ", string(obj.cE(iP))); end
        fprintf("\n")
        disp(["Remaining vertices = ", num2str(obj.nV)])
        disp("Vertex number")
        disp(obj.V)
      end

      function plot(obj)
          obj.polygon.plot
      end
      
      
  end
end