
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
    mE; % slope
    cE; % y intercept
    ind;   % misplaced parameter - needed for convex envelope - remove later
    % remaining vertices
    nV=0;
    V;    % index set

  end 
% 9 methods
  methods  % testing
        function l = checkPiece1 (obj)
           %obj.print
           
           
           l = false;
           if ~ obj.polygon.checkPiece1
               return;
           end
           if obj.nE ~= 1
               return;
           end
           if obj.E ~= [4,1]
               return;
           end
           if obj.mE ~= 1
               return;
           end
           if obj.cE ~= 0
               return;
           end
           if obj.nV ~= 2
               return;
           end
           
           if obj.V ~= [2,3]
               return;
           end
              
           
           l = true;
           return
          
        end

        function l = checkPiece2 (obj)
           %obj.print
           
           
           l = false;
           if ~ obj.polygon.checkPiece2
               return;
           end
           
           if obj.nE ~= 2
               return;
           end
           if obj.E(1,:) ~= [2     3]
               return;
           end
           if obj.E(2,:) ~= [ 3     4]
               return;
           end
           
           if obj.mE ~= [1.500000e+00  1]
               return;
           end
           if obj.cE ~= [1.500000e+00  1]
               return;
           end
           if obj.nV ~= 1
               return;
           end
           
           if obj.V ~= [1]
               return;
           end
              
           
           l = true;
           return
          
        end


        function l = checkPiece3 (obj)
           %obj.print
           
           
           l = false;
           if ~ obj.polygon.checkPiece3
               return;
           end
           if obj.nE ~= 2
               return;
           end
           if obj.E(1,:) ~= [2,3]
               return;
           end
           if obj.E(2,:) ~= [4,1]
               return;
           end
           
           if obj.mE ~= [1, 1.5]
               return;
           end
           if obj.cE ~= [0, 1.5]
               return;
           end
           if obj.nV ~= 0
               return;
           end
           
              
           
           l = true;
           return
          
        end

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
        fprintf("%d  ", obj.mE);
        fprintf("\n")
        fprintf("y-intercepts =  ")
        fprintf("%d  ", obj.cE);
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