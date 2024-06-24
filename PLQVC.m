classdef PLQVC
   % A value class that implements a data type for piecewise linear-quadratic functions
   % based on 
   %    - Shambhavi Singh, Yves Lucet, Linear-time convexity test for low-order piecewise polynomials, 2019
   %    - Convexity and structure of piecewise polynomial functions on polyhedral domains, Singh, Shambhavi, 
   %       MSc Thesis, UBC, 2019. http://hdl.handle.net/2429/70675 
   % the name PLC stands for piecewise linear-cubic (domain is build from linear functions-polyhedral sets; 
   % the function is up to cubic) while V indicates a vertex representation of the convex polyhedral sets defining the
   % domain, and C indicates a coefficient representation of the cubic functions. 
   % Other data structures may use H for half-space representation of polyhedral sets and P for pointwise representation 
   % of quadratic functions. The class is an extension of piecewise linear-quadratic functions to bivariate CUBIC 
   % polynomials with CONNECTED domain forming a polyhedral partition. Disconnected domains are not supported.
   
   properties
        nv {mustBeInteger,mustBeNonnegative}%number of vertices
        ne {mustBeInteger,mustBeNonnegative}%number of edges
        nf {mustBeInteger,mustBeNonnegative}%number of faces
        V (:,2){mustBeNumeric} % nv x 2 matrix storing unique vertices
        E (:,3){mustBeInteger,mustBeNonnegative} % ne x 3 matrix storing edge indices where Edge j is [V(E(j,1)), V(E(j,2))] 
            %and E(j,3) = 1 for a segment and 0 for a ray.
        f (:,10){mustBeNumeric} % nf x 10 matrix storing the coefficients of the cubic
            %f(k,:)=c means at the point (x,y) in the face k  the function value 
            % is C*([x^3; x^2 y; x y^2; y^3; x^2; x y; y^2; x; y; 1] .* [1/6; 1/2; 1/2; 1/6; 1/2; 1; 1/2; 1; 1; 1]
            % the coefficients vector [1/6; 1/2; 1/2; 1/6; 1/2; 1; 1/2; 1; 1; 1] is required to easily manipulate
            % Hessians; it is the same reason we work with 0.5 x^2 instead of x^2, i.e. so that diff(0.5 x^2,1) = 1
        F (:,2){mustBeInteger,mustBeNonnegative} % ne x 2 matrix storing indicating function indices. F(j)=[k1, k2] means
            %the quadratic on the left (resp. right) of edge j has coefficients f(k1) (resp. f(k2)).
            %an index value of 0 indicates the function has value +infinity on that side
        % Special cases: functions with domain dimension less than 2
            %the needle function is stored as V = [0 0];E = [];f = [1 2 3];F = [];%V has 1 row, f has 1 row (and should
            %       be a constant but any cubic is accepted
            %the edge chain function
            %
        P %{cell} %cell array of size nf representing an adjacency list. Each
            %element P{i} is an array of indices k with 1<= abs(k) <= ne(i) 
            %where ne(i) is the number of edges in face i. 
            %If k>0, then the face lies on the right of the edge; otherwise k<0. 
            %The edge indices are ordered to obtain a unique
            %representation. If the face is bounded, the smallest index
            %is on the left, then go clockwise. If the face is
            %unbounded, the first index is the unbounded edge with the
            %smallest index, then go clockwise (the last index is the
            %remaining unbounded edge).
        dom {struct}% struct that stores the domain information (set on which the function is finite); PLC assumes the domain is CONNECTED
            %   dim: nonnegative integer storing the dimension of the domain
            %       dim=0 means the domain is a single vertex (only row in V); PLC is a needle function
            %       dim=1 has 3 cases: segment (2 vertexes, 1 edge), ray (2 vertexes, 1 edge), and boundary of a slice
            %       (3 vertexes, 2 edges). COULD ALSO BE A CHAIN OF EDGES OR ANY GRAPH WITH FUNCTION ONLY DEFINED
            %       ON EDGES; NOT TESTED SINCE THE CLASS FOCUSES MOSTLY ON CONVEX FUNCTIONS
            %   P: if nonempty, the domain has dimension two; in that case P stores the indexes of the face 
            %       representing the domain with the INVERSE convention as P above to obtain a unique representation.
            %   isConvex: boolean. true if the domain is convex
   end
   
   % Class methods
   methods
       function obj = PLQVC(varargin) % constructor    
            ind=3; if nargin==1, ind=1; end%handle Property f in both cases: quadratic and cubic
            [n1,n2] = size(varargin{ind});
            obj.f = [zeros(n1,10-n2), varargin{ind}];%pad with zeros for noncubic functions
            if nargin ==1%quadratic function entered as f, which is the first and only argument            
                obj.nv=0;obj.ne=0;obj.nf=1;
                obj.V=[];obj.E=[];obj.F=[];
            elseif nargin == 4
                obj.V = varargin{1}; obj.E = varargin{2}; obj.F = varargin{4}; 
                [nv,n2] = size(obj.V);[ne,n3] = size(obj.E);[nf, n6] = size(obj.f);[ne1, n22] = size(obj.F); 
                i = [(n2==2), (n3==3) && (ne==ne1), ((n6==6)||(n6==10)) && (nf>0), n22==2  && (ne==ne1)];
                if ~all(i), error(sprintf("Invalid matrix size in PLQVC constructor. Validation vector: [V, E, f, F]=[%i, %i, %i, %i]",i));end %#ok<SPERR>
                obj.nv=nv; obj.ne=ne; obj.nf=nf;
                obj.P = obj.createP();
                %fix nf for degenerate cases
                if ismember([obj.nv,obj.ne],[1,0;2,1],'rows') 
                    %function       nv  ne	nf  dim
                    %needle         1   0   0   0
                    %segment/ray    2   1   0   1
                    obj.nf=0;
                end
                if max(obj.F,[],'all')==0 %domain has dimension < 2, e.g. slice boundary or chain of edges
                    %function       nv  ne	nf  dim
                    %slice boundary 3   2   0   1
                    %chain of edges ?   ?   0   1
                    obj.nf=0;
                end
                if obj.nf>0%for nondegenerate cases, additional checks
                    if (obj.ne>1) && (max(max(obj.F)) ~= size(obj.P,1))
                        error('max index in F must be equal to size(P,1), i.e., number of faces');
                    end
                    if obj.nv>1 && max(max(obj.E)) ~= size(obj.V,1)
                        error('max index in E must be equal to size(V,1)');
                    end
                    if (max(cellfun(@max,cellfun(@abs,obj.P,'UniformOutput',false))) ~= size(obj.F,1))
                        error('max index in P must be equal to size(F,1)');
                    end
                end
             else
                error("Invalid input arguments in PLQVC constructor; need 4 arguments but got %i", nargin);
            end        
            obj.dom = obj.createDom;
       end % constructor
        
       function P = createP(obj)%create object P
            P = cell(obj.nf,1);
            for j = 1:obj.ne
                if obj.F(j,1) > 0 % since zero index is out of domain
                    P{obj.F(j,1)} = [P{obj.F(j,1)},-j];
                end
                if obj.F(j,2) > 0
                    P{obj.F(j,2)} = [P{obj.F(j,2)},j];
                end
            end

            for i = 1:size(P)
                %P{i} = obj.sortRow(P{i},obj.E);
                if ~isempty(P{i})%only sort nonempty arrays
                    P{i} = obj.orderEdges(i);
                end
            end
            if size(P,1)~=obj.nf,error("P must have nf=%i rows, but instead found %i rows",obj.nf,size(P,1));end
       end   
       function [fVal,region] = eval(obj,x)
        % objective:  Evaluates PC function on x
        % [input]  obj of type PLQVC to access V, E, f, P
        %          x : kx2 matrix containing points to evaluate at
        % [output]  fval : kx2 values where fval(i) = PC(x(i)) where fval(i)=nan if PC is discontinuous at x(i)
        %           region : returns region in which the point was found; used for
        %               plot color. Default value is infinity. If it doesn't satisfy
        %               the constraints for any function, the function must be in a
        %               polytope where the value is infinity.
           
            if obj.nv==0%cubic function defined everywhere
                region = ones(size(x,1),1);
                fVal = PLQVC.evalPoly(obj.f,x);
                return;
            elseif max(obj.F,[],'all')==0 %edge chain   
                fVal=Inf * ones(size(x,1),1);
                for i=1:obj.ne%loop on each edge
                    b = PLQVC.belongToEdge(obj.V(obj.E(i,1),:),obj.V(obj.E(i,2),:),x);%true if x belongs to edge i                    
                    if any(b)
                        R = PLQVC.evalPoly(obj.f(i,:),x(b,:)); 
                        ib=find(b);
                        for j=1:length(ib)
                           if isinf(fVal(ib(j)))
                               fVal(ib(j)) = R(j);
                           elseif abs(R(j)-fVal(ib(j))) > sqrt(eps)
                               fVal(ib(j)) = nan;
                           end
                        end
                    end
                end
                region = ones(size(x,1),1);
                return;
            end
            fVal = inf*ones(size(x,1),1);
            region = -ones(size(x,1),1);
            % stores coefficients of the type a^Tx<b as [a1,a2,b] for 2D
            EL = obj.E(:,1:2);  % only the indices of the vertex
            edges = obj.V(EL(:),:);
            m = size(edges,1);
            %need equation of the line going thru (E1,E3) and (E2,E4)
            E1 = edges(1:m/2,1);
            E2 = edges(1:m/2,2);
            E3 = edges(m/2+1:m,1);
            E4 = edges(m/2+1:m,2);
            coefficients = PLQVC.lineEquation([E1 E2],[E3 E4]);%[-(E4-E2), E3-E1, E1.*E4-E3.*E2];

            % to figure out the polytope that contains the point:  
            %   since the number of functions represent the number of polytopes, the
            %   function may be stored twice if two pieces have the same function
            %   (this would only arise in 2 non-adjacent pieces, which has a test case)

            for i = 1:size(obj.P,1)
                % constraint violation for any x will result in the expression to
                % evaluate as 0 for the point
                flags = all(([x,ones(size(x,1),1)]*coefficients(abs(obj.P{i}),:)').*sign(obj.P{i})<=0,2);

                if any(flags==1)
                    % Checks if any x lies in the current polytope, evaluates them and
                    % marks border points as a different colour
                    %fVal(flags==1) = obj.plqvc_funcval2d(f(i,:),x(flags==1,:)');
                    fVal(flags==1) = PLQVC.evalPoly(obj.f(i,:),x(flags==1,:));
                    region(region~=-1 & flags==1) = 0;
                    region(region==-1 & flags==1) = i;
                end
            end
       end
       function dom = createDom(obj)
           % return the domain of the PLQVC function; in addition store it in obj.Dom
           % output: Dom: struct as explained in the class properties
           % dim=0
           dom.dim=0;dom.P=[];dom.isConvex=false;
           if (obj.nv==0)%cubic function with full domain
               if obj.nf~=1, error("Error in Dom: a function with no vertex must have exactly 1 face");end
               dom.dim=2;dom.isConvex=true;return;
           end
           if obj.nv==1,dom.isConvex=true;return;end%dim(domain)=0: needle function
           if obj.nf==0%dim(domain)=1: segment, ray, boundary of slice
               dom.dim=1;
               if (obj.nv==2) && (obj.ne==1) %segment or ray
                   dom.isConvex=true;
               else %chain of edges
                   if PLQVC.isCollinear(obj.V,sqrt(eps))
                       dom.isConvex=true;
                   end
                   %Error("Found 0 face but not full domain, not needle, not segment/ray, or not boundary of slice. What is it?");
               end
               return;
           end
           jMin = find(obj.F==0, 1);%warning return indexes from 1 to 2*ne
           if isempty(jMin) %function is finite everywhere
               dom.dim=2;dom.isConvex=true;return;
           end
           %build boundary of domain, which is a face (dim=2)
           dom.dim=2;
           [dom.P, dom.isConvex] = obj.orderEdges(0);
           %NEED TO REORDER THE EDGE INDEXES IF WE WANT THE SAME CONVENTION AS P
           %for now the result is a list of edges counter-clockwise with negative edge indicating the face is on the
           %right, i.e. it is the inverse of P
       end       
       function [P, isConvex] = orderEdges(obj,k)
        % Return P edge indexes in clockwise order of face k, and determine whether face k is convex
        %
        % input:    k: face index whose boundary you wish to compute and convexity to determine
        %
        % output:   P: 1xnep array of nep edge indexes ordered clockwise and starting with ray of lower index
        %               value or, when no ray are present, from edge index of lower value
        %           isConvex: boolean; true when the polyhedral set defined by the edges list is convex
            if obj.nv==1,isConvex=true;P={};return;end%needle function; cannot describe the domain as a list of edges (no edge)
            edges = find(any(obj.F==k,2));%indexes of all the edges on the boundary of Face k
            rays = edges(obj.E(edges,3)==0);
            Angles = zeros(1,length(edges)-1);%integer vector storing whether angles between edges(j),edges(j+1) is acute (1), collinear (0), or obtuse (-1)
            if ~isempty(rays)
                if length(rays)~=2
                    error("orderEdges: Face %i has % rays, but must be zero or 2",k,length(rays));
                end
                %there are 2 rays on the boundary; one has the Face k on its left and the other on its right
                %because we describe the edges of the boundary clockwise, we start with the
                %ray having the Face k on its right
                if obj.F(rays(1),1)==k%domain on the right
                    j = rays(1);
                else
                    j = rays(2);
                end
           else%start with edges of minimal index in absolute value
                j = find(any(obj.F==k,2), 1 );%find first edge on boundary of Part/Face k
            end
           %traverse the chain from j until jEnd 
           nedges = length(edges);%number of edges on Face/Part k
           P=zeros(1,nedges);%vector storing indices of edges on the boundary in clockwise order
           A=obj.E(edges,1:2);%only consider edge indexes on Part/Face k
           jj=1;
           if obj.F(j,1)==k%face on the left of edge j so store -j
                P(jj)=-j;
           else%face on right of edge j so store +j
               P(jj)=j;
           end

           for jj=2:nedges 
               %consider edge j find i the vertex index of its endpoint
               if obj.E(j,3)==0 %ray so take base point
                   i = obj.E(j,1);%vertex index
               else
                   if obj.F(j,1)==k %face on the left, take previous edge to move clockwise so take base point of edge
                        i = obj.E(j,1);
                   else %face on the right, take next edge to move clockwise so take end point of edge
                        i = obj.E(j,2);
                   end
               end               
               %looking at edge j between vertex index i and vertex index iEndpoint
               iEndpoint=sum(obj.E(j,1:2),2) - i;
               if jj==2, iFirst=iEndpoint;iEndpointFirst=i;end %needed for testing last angle in bounded polyhedral set below
               
               Ind=find(A==i);%should only return 2 indexes
               if length(Ind)~=2,error("The boundary must be a chain of edges so only 2 edges share the same vertex on the boundary; expected 2 but got %i",length(Ind));end
               [r,c]=ind2sub(size(A),Ind);%r and c have both 2 elements
               if A(r(1),mod(c(1),2)+1)==iEndpoint %the modulo + 1 function maps 1->2, and 2->1
                   iNext = A(r(2),mod(c(2),2)+1);
               elseif A(r(2),mod(c(2),2)+1)==iEndpoint
                   iNext = A(r(1),mod(c(1),2)+1);
               else
                   error("Unable to compute iNext");
               end
                   
               jNext = find(ismember(obj.E(:,1:2),[i, iNext],'rows'));
               if isempty(jNext)
                   jNext = find(ismember(obj.E(:,1:2),[iNext, i],'rows'));
               end
                              
               %compute sign of angle between edges j and jNext; 1 means acute, -1 means obtuse, 0 means collinear
               %https://matlabgeeks.com/tips-tutorials/computational-geometry/check-convexity-of-polygon/
               v1 = obj.V(iEndpoint,1:2) - obj.V(i,1:2);
               v2 = obj.V(iNext,1:2) - obj.V(i,1:2);  
               Angles(jj-1)= sign(det([v1;v2]));
            
               if obj.F(jNext,1)==k%face on the left of edge j so store -j
                    P(jj)=-jNext;
               else%face on right of edge j so store +j
                   P(jj)=jNext;
               end               
               j=jNext;
           end
           %for bounded polyhedral set, need to compute the angle between the first edge and the last edge
           if iNext==iFirst
               v1 = obj.V(i,1:2) - obj.V(iFirst,1:2);
               v2 = obj.V(iEndpointFirst,1:2) - obj.V(iFirst,1:2);  
               Angles(nedges+1)= sign(det([v1;v2]));
           end
           
           if k==0
               isConvex= ~any(Angles==1);%boundary of domain is complement; it is convex if all edge angles are obtuse
           else
               isConvex= ~any(Angles==-1);%all edge angles are acute or edges are collinear 
           end                  
       end
       function b = isDomBounded(obj)
           if obj.ne==0 % full domain or needed function
               if obj.nv==1 %needle function
                   b=true;
               else%full domain
                   b=false;
               end
               return
           end
           b = all(obj.E(:,3));%all segments so no ray hence bounded
       end       
       function isConvex = isFaceConvex(obj,k)
       %determines whether the function restricted to Face k is convex
       %input: k: face index
       %output: isConvex: true if the face is convex and the function is convex on that face
            c = obj.f(k,:);
            [~,Q,C] = PLQVC.matrixForm(c);
            isConvex = obj.dom.isConvex;if ~isConvex, return;end%check that domain is convex first
            
            if sum(abs(C),'all') < sqrt(eps)  %quadratic function; treat as special case for clarity AND performance
                                    % e.g. face with a large number of vertexes
                isConvex = PLQVC.isPositiveSemidefinite(Q);
                return;
            else%cubic function; need Hessian to be positive semi-definite on all vertexes & all ray directions
                J = abs(obj.P{k})'; %edge indexes of Face k                
                M = obj.E(J,:);
                Iv = M(:,1:2);%all vertex indexes
                I = unique ([Iv(:,1);Iv(:,2)]);
                Vect = obj.V(I,:);
                Ir = M(M(:,3)==0,1:2);%vertex index of rays; Ir is empty or a 2x2 matrix
                if ~isempty(Ir)
                    if ~all(size(Ir)==[2, 2]), error("Each face can only have 0 or 2 rays but found %i",size(Ir,1));end
                    Dir = obj.V(Ir(:,2),:) - obj.V(Ir(:,1),:);
                    Vect = [Vect; Dir];
                end
                H = PLQVC.evalHessian(obj.f(k,:),Vect);%returns hypermatrix of size 2x2xsize(Vect,1)
                isConvex = PLQVC.isPositiveSemidefinite(H);
            end
        
       end
       function [isConvex, isContinuous] = isEdgeConvex(obj,j)
       %determines if edge j breaks the convexity or continuity of the piecewise function. Checks continuity and the
       %convexity of the 2 faces (in the domain) bordering the edge
       %input: j: edge index
       %output: isConvex: true if the function is continuous accross the edge and its convex subdifferential is nonempty
       %            at every point on the edge (the edge can be a segment or a ray)
       %        isContinuous: true if the function is continuous accross the edge
            if isempty(obj.F), error("isEdgeConvex:NoEdge","The current function has no edge");end
            k = sum(obj.F(j,:)>0);%number of nonzero face indexes for edge j
            if k==0%case function is only defined on the edge
                isContinuous=true;%continuity is always true for cubic polynomial
                isConvex = PLQVC.isCubicConvexOnEdge(obj.f(j,:),obj.V(obj.E(j,1),:),obj.V(obj.E(j,2),:),obj.E(j,3));
                return;
            end
            
            if k==1%case bordering the domain
                isContinuous = true;
                isConvex = obj.isFaceConvex(k(k~=0));
                return;
            end
            if k~=2,error("Number of nonzero indexes must be 0, 1, or 2 but found %i.",k);end
            %now the function is finite on both side of the edge
            isContinuous = obj.isEdgeContinuous(j);
            if ~isContinuous, isConvex=false; return;end
            %check face convexity
            isConvex = obj.isFaceConvex(obj.F(j,1)) && obj.isFaceConvex(obj.F(j,2));
            if ~isConvex, return;end

            % We use notations from
            %cubic function with edge between 2 faces on which the function is finite
            %Singh, S. (2019). Convexity and structure of piecewise polynomial functions on polyhedral domains (T).
            %Algorithm 3, p. 40. University of British Columbia. Retrieved from http://hdl.handle.net/2429/70675
            % The common boundary is a^T x = beta; it goes through V(1,:) and V(2,:)            
            V = [obj.V(obj.E(j,1),:);obj.V(obj.E(j,2),:)];%V1=V(1,:); V2=V(2,:);
            c = PLQVC.lineEquation(V(1,:),V(2,:));a=c(1,1:2)'; beta=c(1,3);
            
            [~,g1] = PLQVC.evalHessian(obj.f(obj.F(j,1),:),V);
            [~,g2] = PLQVC.evalHessian(obj.f(obj.F(j,2),:),V);
            delta1 = (g1(:,1)-g2(:,1)) ./ a;
            if norm(delta1(1)-delta1(2)) > sqrt(eps), error("Unable to compute delta1");end
            delta2 = (g1(:,2)-g2(:,2)) ./ a;
            if norm(delta2(1)-delta2(2)) > sqrt(eps), error("Unable to compute delta2");end
            delta1=rmmissing(delta1);delta2=rmmissing(delta2);%division by 0 results in nan; those values are filtered out here
            delta=[delta1(1);delta2(1)];

            isConvex = delta(1)>-sqrt(eps) && delta(2) > -sqrt(eps);
            if obj.E(j,3)==0 %ray
               isConvex = isConvex && (delta(1) < delta(2)+sqrt(eps));
            end
            
            if norm(obj.f(obj.F(j,1),:))<sqrt(eps) && norm(obj.f(obj.F(j,2),:))<sqrt(eps) %function is quadratic on adjacent faces
%                 %Theorem 3.8, p. 28.                
%                 [L1,Q1] = PLQVC.matrixForm(obj.f(obj.F(j,1),:));
%                 [L2,Q2] = PLQVC.matrixForm(obj.f(obj.F(j,2),:));
% 
%                 delta1 = ((Q1-Q2)*V(1,:)'+L1-L2) ./ a;
%                 if norm(delta1(1)-delta1(2)) > sqrt(eps), error("Unable to compute delta1");end
%                 delta2 = ((Q1-Q2)*V(2,:)'+L1-L2) ./ a;
%                 if norm(delta2(1)-delta2(2)) > sqrt(eps), error("Unable to compute delta2");end
%                 delta1=rmmissing(delta1);delta2=rmmissing(delta2);%division by 0 results in nan; those values are filtered out here
%                 delta=[delta1(1);delta2(1)];
%                 
%                 isConvex = delta(1)>-sqrt(eps) && delta(2) > -sqrt(eps);
%                 if obj.E(j,3)==0 %ray
%                    isConvex = isConvex && (delta(1) < delta(2)+sqrt(eps));
%                 end
                return                
            end
            %now truly cubic function
            %Given coefficients c, we can compute the gradient with
%             syms x y
%             assume(x,'real');assume(y,'real');
%             X=[x, y];
%             c = sym('c',[1,10]); assume(c,'real');
%             [H,g] = PLQVC.evalHessian(c,X);
%             expand(g)
            %and we obtain
            %g = [(c1*x^2)/2 + c2*x*y + c5*x + (c3*y^2)/2 + c6*y + c8; (c2*x^2)/2 + c3*x*y + c6*x + (c4*y^2)/2 + c7*y + c9];
            %so the gradient can be written as
            %g = [c1;c2] .* x^2/2 + [c2;c3] .* x*y + [c3;c4] .* y^2/2 + [c5;c6] .* x + [c6;c7] .* y + [c8;c9]
            %now the gradient restricted to a^T [x;y] = beta is
            %
            
       end
       function b = isQuadratic(obj,tol)
       %true if the function is quadratic up to floating point errors
       %note: the function is not cleaned up so it may still have small cubic coefficients
       %use eps instead of sqrt(eps) to minimize floating point issues 
          if nargin < 2, tol = eps;end
          b = norm(obj.f(:,1:4)) < tol;
       end
       function b = isEdgeContinuous(obj,j)
       %check whether the PC function is continuous across edge j
       %input: j: edge index
       %output: b: boolean; true if the PC function is continuous accross edge j
            %if edge is on the boundary of the domain, the function is continuous; nothing to check
            k = sum(obj.F(j,:)>0);
            if k==0 || k==1, b=true; return;end
            
            %Since PC is a cubic polynomial, we only need to check 4 points on the edge to ensure continuity
            V1=obj.V(obj.E(j,1),:);V2=obj.V(obj.E(j,2),:);t=linspace(0,1,4);
            X=(V1' + t.*(V2-V1)')';%4 equi-spaced points between V1 and V2 included
            c1 = obj.f(obj.F(j,1),:);c2=obj.f(obj.F(j,2),:);
            fVal1 = PLQVC.evalPoly(c1,X);fVal2=PLQVC.evalPoly(c2,X);
            b = norm(fVal1-fVal2)<sqrt(eps);
       end
       function b = isConvex(obj)%TO BE TESTED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       %returns true is the function is convex
            if ~obj.dom.isConvex, b=false;end
            if obj.nf==0%special cases
                error("Not implemented yet");
            end
            b=false;
            for k=1:obj.nf%interior convexity, i.e. the restriction of the function to each must be convex                
                if ~obj.isFaceConvex(k), return;end
            end
            %edge convexity, i.e. check that the function is convex when accross each edge
            for j=1:obj.ne
                if ~obj.isEdgeConvex(j), return;end
            end
            b = true;
       end
       function disp(obj)
           strPow=["x^3","x^2y","xy^2","y^3","x^2","xy","y^2","x","y",""];%unweighted basis
           S=strings(obj.nf,1);
           for i=1:obj.nf
               c=obj.f(i,:);
               Cw= [1/6, 1/2, 1/2, 1/6, 1/2, 1, 1/2, 1, 1, 1];%weights to avoid constant coefficients in computation of Hessian
               c = c.*Cw;
               %if size(f,2)<10, f=[zeros(1,10-size(f,2)), f];end%pad with zeros for noncubic functions
               str="";
               for j=1:size(c,2)
                   v=c(1,j);
                   if v~=0
                       s=string(v);
                       if s == "-1", s="-";elseif s == "1", s="";end
                       if v>0, s="+ "+s;end
                       str=str + s + " " + strPow(j) + " ";
                   end
               end
               if startsWith(str,"+")
                   str=extractBetween(str,2,strlength(str));
               end
               if str=="", str="0";end
               S(i)=str;
           end
           disp(join([S,repmat("if (x,y) in P(",obj.nf,1),(1:obj.nf)',repmat(")",obj.nf,1)]));
       end             
       function plotDomain(obj)
            if ~obj.isDomBounded, error("plotDomain only implemented for bounded domain for now");end       
            hold on;
            warning('off','MATLAB:polyshape:repairedBySimplify');%may produce warning when an edge is split in 2 along boundary; in that case, polyshape simplifies the polyhedral set; ignore warning
            for iFace = 1:obj.nf % loop on each face and plot the associated polygon
                % build a list of vertices clockwise that make up the face
                face = obj.P{iFace};
                iVs = zeros(size(face));
               for i = 1:length(face)
                   j = face(i);
                   if j>0
                       iVs(i) = obj.E(j,1);
                   else
                       iVs(i) = obj.E(-j,2);
                   end
               end
               %disp(obj.V(iVs,:));%for debugging: display vertices coordinates in order
               pgon = polyshape(obj.V(iVs,:));
               plot(pgon);
            end 
            hold off;
            warning('on','MATLAB:polyshape:repairedBySimplify');           
       end
       function plotVertexes(obj)
           plot(obj.V(:,1),obj.V(:,2),'+');
       end
       function plotDomainGraph(obj)
            G = digraph(obj.E(:,1),obj.E(:,2),1:obj.ne);
            H=plot(G,'XData',obj.V(:,1),'YData',obj.V(:,2),'EdgeLabel',G.Edges.Weight,'EdgeLabelColor','b');
            %highlight rays
            j = find(obj.E(:,3)==0);
            highlight(H,obj.E(j,1),obj.E(j,2));
            %number faces
            x=zeros(obj.nf,1);y=x;
            for i=1:obj.nf
                J1 = obj.E(abs(obj.P{i}),1:2);
                J2 = unique ([J1(:,1);J1(:,2)]);
                n=length(J2);
                x(i)=sum(obj.V(J2,1))/n;y(i)=sum(obj.V(J2,2))/n;
            end
            text(x,y,string(1:obj.nf),'Color','m');
       end
       function plot(obj,xAxis,yAxis)
            % Plots a plq function
            n = 201;edgecolor = 'none';
            if nargin==1, xlim = [-1,1];ylim = [-1,1];
            else, xlim = xAxis; ylim = yAxis;
            end
            colormap winter;
            [X,Y]=meshgrid(linspace(xlim(1),xlim(2),n),linspace(ylim(1),ylim(2),n));
            [fVal,region] = obj.eval([X(:),Y(:)]);
            %plq.x = [X(:),Y(:)];[e,c]=plqvc_evaluate2d(plq);
            Z=reshape(fVal,[],size(X,2));
            C=reshape(region,[],size(X,2));
            surf(X,Y,Z,C,'EdgeColor',edgecolor);%'[0.5 0.5 0.5 ]');
            xlabel('x');ylabel('y');
       end
       
    end % methods
   
   methods(Static)
        function [H,g] = evalHessian(c,x)
        %returns Hessians and Gradients of the cubic function given by coefficients c at the points x
        %input  c: 1x10 coefficient matrix defining a cubic function (leading zeros can be ignored)        
        %       x: Nx2 array storing points x(i,1:2). If x is ommitted c must be a quadratic function
        %output H: 2x2xN hypermatrix such that H(:,:,i) is the Hessian at point x(i,:)
        %       g: 2xN matrix such that g(:,i) is the gradient at point x(i,:)
            [L,Q,C] = PLQVC.matrixForm(c);
            if (nargin < 2) && (nargout<2)
                if (sum(abs(C),'all') < sqrt(eps))%quadratic function does not need x when no gradient is requested
                    H = Q;g=[];return;
                else
                    error('evalHessian:xNotProvidedForCubic','A nonquadratic function requires input x to evaluate its Hessian at x, but x was not provided; or gradient was requested but x was not provided');
                end
            end
            %H = zeros(2,2,size(x,1));g = zeros(2,size(x,1));%preallocate; Preallocating prevent running symbolically
            for i=1:size(x,1)                
                if isempty(C)%quadratic function but Hessian is required
                    H(:,:,i) = Q;
                    g(:,i) = L + Q*x(i,:)';
                else
                    H(:,:,i) = Q + (C(:,:,1)*x(i,1)+C(:,:,2)*x(i,2));
                    g(:,i) = L + Q*x(i,:)' + 0.5*(C(:,:,1)*x(i,1)+C(:,:,2)*x(i,2))*x(i,:)';
                end                
            end
        end
        function b = isPositiveSemidefinite(H)
        %returns true if all H(:,:,i) matrices are positive semidefinite
        %input H: hypermatrix n x n x N 
        %output b: boolean; true if for all i=1:N, H(:,:,i) is positive semidefinite     
            b = false;
            for i=1:size(H,3)
                if ~all(eig(H(:,:,i)) > -sqrt(eps)), return;end %account for floating point errors
            end
            b = true;
        end
        function isConvex = isCubicConvexOnEdge(c,v1,v2,isSegment)
        %Given a cubic polynomial determine whether it is convex on the edge between v1 and v2
        %input  c: 1x10 array of coefficients of the cubic polynomial; if size(c,1)<10 then c is padded with zeros
        %       v1,v2: 1x2 vectors; coordinates of the boundaries of the edge
        %       isSegment: true: segment between v1 and v2, false: ray from v1 passing thru v2 till infinity
            syms x y;
            C = [x.^3, x.^2.*y, x.*y.^2, y.^3, x.^2, x.*y, y.^2, x, y, 1];
            Cw= [1/6, 1/2, 1/2, 1/6, 1/2, 1, 1/2, 1, 1, 1];%weights to avoid constant coefficients in computation of Hessian
            C = C.*Cw;
            c = [zeros(1,10-size(c,2)),c];
            P = C*c';%cubic polynomial
            coeff = PLQVC.lineEquation(v1,v2);
            if abs(coeff(2))<sqrt(eps) %line x = constant                
                f = subs(P,x,-coeff(3)/coeff(1));%function on y only
                ddf = diff(f,2);%this is a linear function
                ddf1 = double(subs(ddf,y,v1(2))); ddf2 = double(subs(ddf,y,v2(2)));%value of 2nd derivative at vertexes
            else
                f = subs(P,y,-coeff(1)/coeff(2)*x-coeff(3)/coeff(2));%function of x only
                ddf = diff(f,2);%this is a linear function
                ddf1 = double(subs(ddf,x,v1(1))); ddf2 = double(subs(ddf,x,v2(1)));%value of 2nd derivative at vertexes
            end
           if isSegment
                isConvex =  (ddf1 > -sqrt(eps)) && (ddf2 > -sqrt(eps));%2nd derivative is nonnegative at both endpoints
            else
                isConvex = (ddf1 < ddf2 + sqrt(eps));%line going downward will be come negative eventually
            end             
        end
        function c = lineEquation(v1,v2)
        %Coompute the equation of the line going thru v1 and v2
        %input  v1,v2: 1x2 vectors; coordinates of the points thru which the line goes thru. 
        %               If v1,v2 are nx2 arrays then c returns a matrix nx3 of the coefficients of each line.
        %output c: 1x3 coefficients of the equation of the line as c(1) x + c(2) y + c(3) = 0
        %           if v1,v2 are nx2 arrays then c is a nx3 matrix
        %
        %   need equation of the line going thru v1 and v2 (E1,E3) and (E2,E4)
        %   normal vector n=(-(v2(2)-v1(2), v2(1)-v1(1)) (-(E4-E2),(E3-E1)) going thru point v1, v2 (E1,E2). So equation is dot(n,((x,y)-v1))=0, 
        %   which gives n1 x+ n2 y -n1 v1(1)E1 - n2 v1(2) E2=0. Stored as a x + b y  + c =0, 
        %   we get a = -(v2(2)-v1(2)E4-E2), b=(v2(1)-v1(1)E3-E1), and c=v1(:,1).*v2(:,2)-v2(:,1).*v1(:,2)E1.*E4-E3.*E2
            c = [-(v2(:,2)-v1(:,2)), v2(:,1)-v1(:,1), v1(:,1).*v2(:,2)-v2(:,1).*v1(:,2)];
        end       
        function fVal = evalMatrixForm(c,x)
        % objective: evaluates cubic polynomial given by coefficients by first converting to matrix form and then 
        %           with matrix multiplications.
        %           WARNING: for testing purposes mostly; use eval for evaluating the PLC function directly from its
        %           coefficient vector c
        %
        % [input]  c:   1x6 coefficient matrix fVal = c*[x^2 xy y^2 x y const] for quadratic or 
        %               lx10 matrix associated f(x) = c*[x^3; x^2 y; x y^2; y^3; x^2; x y; y^2; x; y; const] for cubic
        %           x : kx2 vector of points 
        %
        % [output]  fVal : kx1 vector of evaluated points at x such that f(x(i,:)) = c' * x
        %               
            if size(c,2)>10
                error("Not implemented for coefficient vector of size greater than 10, i.e. for polynomial of degree higher than cubic");
            end 
            [L,Q,C] = PLQVC.matrixForm(c);
            fVal = 0.5 * dot((x*Q)',x')' + x * L + c(end);%vectorized version
            if ~isempty(C)%cubic function
                for i=1:size(x,1)
                    fVal(i) = fVal(i) + 1/6 * x(i,:) * (C(:,:,1)*x(i,1)+C(:,:,2)*x(i,2))*x(i,:)';
                end
                %the above formula is equivalent to the following formula, which is more cryptic but generalizes to 
                %higher dimensions if needed. Couldn't find any vectorization in pure linear algebra like for quadratic.
                %fVal = fVal + 1/6 * x' * (sum(C.*permute(x,[3,2,1]),3))*x;
            end                     
        end
        function [L,Q,C] = matrixForm(c)
        % computes matrix form of cubic/quadratic function given by coefficient vector
        % f(x) = c'*x = 1/2 x'*Q*x + L'*x + c(10) for quadratic; similar for cubic with 2x3x3 tensor matrix C
        % input c : 1x6 coefficient matrix in the order of [x^2 xy y^2 x y const] for quadratic or
        %             1x10 matrix for [x^3, x^2y, xy^2, y^3, x^2, xy, y^2, x, y, const]
        % output L : L'*[x;y] = c(8:9)*[x;y]
        %        Q : [x;y]' * Q' * [x;y] = [x;y]' * c(5:7) * [x;y]
        %        C : C(1,:,:)*[x;y]+C(2,:,:)*[x;y] =  x ([x;y]' * c(1:4) * [x,y]) + y ([x;y]' * c(1:4) * [x,y])
            if size(c,2)<10, c=[zeros(1,10-size(c,2)), c];end%padd with zeros for noncubic functions
            Cst = [1 1;1 1];
            C=cat(3,[c(1) c(2);c(2) c(3)].*Cst,[c(2) c(3);c(3) c(4)].*Cst);%tensor matrix
            Q=[c(5) c(6);c(6) c(7)];%quadratic part
            L=[c(8) c(9)]';%linear part
            if sum(abs(C),'all')==0,C=[];end%remove padding for noncubic functions
        end       
        function z = evalPoly(c,X)
        %evaluate cubic polynomial p (given by its coefficients c) at points X
        %input  c 1xn vector where n<=10; n=6 for quadratic, and n=10 for cubic
        %       X kx2 matrix storing points
        %
        %output z kx1 vector such that y(i)=p(X(i,:))
            if size(c,1)>1,error("evalPoly:notVectorizedInC","evalPoly is vectorized for multiple points, but not for multiple polynomials");end
            x=X(:,1);y=X(:,2);
            C = [x.^3, x.^2.*y, x.*y.^2, y.^3, x.^2, x.*y, y.^2, x, y, ones(size(x))];   
            Cw= [1/6, 1/2, 1/2, 1/6, 1/2, 1, 1/2, 1, 1, 1];%weights to avoid constant coefficients in computation of Hessian
            C = C.*Cw;
            z = C*[zeros(1,10-size(c,2)), c]';    
        end
        function b = belongToEdge(V1,V2,X,isSegment,tol)
        % returns true if X belongs to segment [V1,V2] when isSegment=true; to ray [V1,V2) when isSegment=false
        % input V1, V2: 1x2 row vectors defining the segment [V1,V2]
        %       X: n x d matrix of points so that the ith point is X(i,:)
        %       isSegment: boolean; true if considering segment [V1,V2], false if considering ray [V1,V2)
        %       tol: floating point number indicating tolerance. Default to sqrt(eps)=1.4901e-08
        % output b: n x 1 boolean vector with b(i) = true if X(i,:) belongs to segment [V1,V2]
            if nargin < 4, isSegment=true;end
            if nargin < 5, tol=sqrt(eps);end
            if any(size(V1)~=[1,2]) || any(size(V1)~=size(V2)),error("belongToEdge:dimensionMismatch","V1 and V2 must be the same dimension");end
            b = false(size(X,1),1);
            for i=1:size(X,1)%isCollinear is vectorized on X, but not on V1 and V2 so need to loop
                b(i) = PLQVC.isCollinear([V1;V2;X(i,:)],tol);
            end
            Xb = X(b,:);%only check withing range for collinear points
            if isSegment
                b(b) =  all(min(V1,V2) <= Xb + tol,2) & all( Xb - tol <= max(V1,V2),2);
            else
                b(b) = all(V1 <= Xb + tol,2);
            end
        end
        function b = isCollinear(P, tol)
        % returns true if points P are collinear
        %   b = 1 if the points P(i,:),  i=1,...,size(P,1) are collinear up to a given
        %   tolerance TOL. Works for any dimension (number of columns in P) for any 
        %   number of points (rows of P). If not given, TOL = 0.
        %
        %   Example:
        %       % Two points are always collinear
        %       collinear([1 0; 1 5]); % returns true
        %       % Three points in 3D which are supposed to be collinear
        %       collinear([0 0 0; 1 1 1; 5 5 5]); % returns false due to numerical error
        %       % The previous example with looser tolerance
        %       collinear([0 0 0; 1 1 1; 5 5 5], 1e-14); % returns true
        %
        %   See also:  RANK
        %
        %   Source: https://www.mathworks.com/matlabcentral/fileexchange/62737-collinear
        %
        %   The algorithm for three points is from Tim Davis:
        %       http://blogs.mathworks.com/loren/2008/06/06/collinearity/#comment-29479
        %
        %   Cite as
        %   Zoltán Csáti (2020). collinear (https://www.mathworks.com/matlabcentral/fileexchange/62737-collinear), 
        %   MATLAB Central File Exchange. Retrieved July 7, 2020.
            if nargin < 2
                tol = 0;
            end
            b = rank(bsxfun(@minus, P, P(1,:)), tol) < 2;
        end
        function p = oneNorm()
            V = [0 0; -1 0; 0 1; 1 0; 0 -1];
            E = [1 2 0; 1 3 0; 1 4 0; 1 5 0];
            f = [0 0 0 -1 -1 0;0 0 0 -1  1 0;0 0 0  1  1 0;0 0 0  1 -1 0];
            F = [1 2;2 3;3 4;4 1];
            p = PLQVC(V,E,f,F); 
        end
        function p = oneNormConjugate()
            V = [-1 -1; -1 1; 1 1; 1 -1];
            E = [1 2 1; 2 3 1; 3 4 1; 4 1 1];
            f = [0 0 0 0 0 0];
            F = [0 1; 0 1; 0 1; 0 1];
            p = PLQVC(V,E,f,F);
        end
        function p = energy()%half the norm squared
            p = PLQVC([1 0 1 0 0 0]);%x^2/2 + y^2/2
        end
        function p = cubic1()
            %cubic function defined on the 4 quadrants
            V = [0 0;-1 0;0 1;1 0;0 -1];
            E = [1 2 0;1 3 0;1 4 0;1 5 0];
            f = [-1 0 0 -1 1 0 1 0 0 0;
                 -1 0 0  2 1 0 2 0 0 0;
                  3 0 0  2 2 0 2 0 0 0;
                  3 0 0 -1 2 0 1 0 0 0];        
            Cw= ones(size(f,1),1)*(1./[1/6, 1/2, 1/2, 1/6, 1/2, 1, 1/2, 1, 1, 1]);
            f = f.* Cw;%fix after modifying weights internally
            F = [1 2;2 3;3 4;4 1]; 
            p = PLQVC(V,E,f,F);
        end
        function P = examples()
           %P=PLQVC.examples; for i=1:length(P),plot(P{i});title(i+": "+P{i}.nf+" faces");pause;end
           % 1 One norm
            P{1}=PLQVC.oneNorm();
            P{2}=PLQVC.oneNormConjugate();
            P{3}=PLQVC.energy();
            P{4}=PLQVC.cubic1();
            %all examples below were fixed after a weight change in the eval functions
            invW = 1./[1/6, 1/2, 1/2, 1/6, 1/2, 1, 1/2, 1, 1, 1];
            % 5 simple convex function (containing 4 pieces)
            %      x^2+ y^2     x<0    y<0
            %      x^2+2y^2     x<0    y>=0
            %      2x^2+2y^2    x>=0   y>=0
            %      2x^2+ y^2    x>=0   y<0
            V = [0 0;-1 0; 0 1;1 0;0 -1];
            E = [1 2 0;1 3 0;1 4 0;1 5 0];
            f = [1 0 1 0 0 0;1 0 2 0 0 0;2 0 2 0 0 0;2 0 1 0 0 0];
            %fix after modifying weights internally
            f = [zeros(10-size(f,2),size(f,1)),f];%padd with zeros to force being of same size as Cw
            Cw= ones(size(f,1),1)* invW;
            f = f .* Cw;
            F = [1 2;2 3;3 4;4 1];
            P{end+1} = PLQVC(V,E,f,F);
            % 6 Concave function at multiple points
            %     -x^2- y^2     x<0    y<0
            %     -x^2-2y^2     x<0    y>=0
            %     -2x^2-2y^2    x>=0   y>=0
            %     -2x^2- y^2    x>=0   y<0
            V = [0 0;-1 0;0 1;1 0;0 -1];
            E = [1 2 0;1 3 0;1 4 0;1 5 0];
            f = -[1 0 1 0 0 0;1 0 2 0 0 0;2 0 2 0 0 0;2 0 1 0 0 0];
            f = [zeros(size(f,1),10-size(f,2)),f];
            Cw= ones(size(f,1),1)* invW;
            f = f .* Cw;
            F = [1 2;2 3;3 4;4 1];
            P{end+1} = PLQVC(V,E,f,F);
            % 7 Finite polytope with infinity on the outside
            %       x^2+2xy+y^2     -1<=x<=1,-1<=y<=1
            %       infinity        otherwise
            V = [-1 -1;1 -1;1  1;-1  1];
            E = [1 2 1;2 3 1;3 4 1;4 1 1];
            F = [1 0;1 0;1 0;1 0];
            f = [1 2 1 0 0 0];
            f = [zeros(size(f,1),10-size(f,2)),f];
            Cw= ones(size(f,1),1)* invW;
            f = f .* Cw;
            P{end+1} = PLQVC(V,E,f,F);
            % 8 Difference indefinite
            %       x^2+y^2     x<0    y>=0
            %       x^2+xy+y^2  x>=0   y>=0
            %       infinity            y<0
            V = [0 0;0 1;  1 0;-1 0];
            E = [1 2 0;1 3 0;1 4 0];
            f = [1 0 1 0 0 0;1 1 1 0 0 0];
            f = [zeros(size(f,1),10-size(f,2)),f];
            Cw= ones(size(f,1),1)* invW;
            f = f .* Cw;
            F = [1 2; 2 0;0 1];
            P{end+1} = PLQVC(V,E,f,F);
            % 9 Single triangular Piece (Non box domain)
            %       x^2+2xy+y^2     (x,y) in Triangle with vertexes (0,0),(1,0),(0.5,1)
            %       infinity        otherwise
            V = [0 0;1 0;0.5 1;];
            E = [1 2 1;2 3 1;1 3 1];
            F = [1 0; 1 0; 0 1];
            f = [1 2 1 0 0 0];
            f = [zeros(size(f,1),10-size(f,2)),f];
            Cw= ones(size(f,1),1)* invW;
            f = f .* Cw;
            P{end+1} = PLQVC(V,E,f,F);
            % 10 Unbounded triangle piece (Non box domain)
            %       x^2+2xy+y^2     5y/6 <= x <= 10y    y>=0
            %       infinity        otherwise
            V = [0 0;1 0.1;.5 .6];
            E = [1 2 0;1 3 0];
            F = [1 0; 0 1];
            f = [1 2 1 0 0 0];
            f = [zeros(size(f,1),10-size(f,2)),f];
            Cw= ones(size(f,1),1)* invW;
            f = f .* Cw;
            P{end+1} = PLQVC(V,E,f,F);
            % 11 Fenchel Conjugate of a C1 2 piece function
            %       1/4(x^2-2xy+2y^2)   x>=y
            %       1/8(x^2-2xy+3y^2)   x<y
            V = [0 0; 1 1; -1 -1];
            E = [1 2 0; 1 3 0];
            f = [1/4 -1/2 1/2 0 0 0;
                 1/8 -1/4 3/8 0 0 0];
            f = [zeros(size(f,1),10-size(f,2)),f];
            Cw= ones(size(f,1),1)* invW;
            f = f .* Cw;
            F = [2 1;1 2];
            P{end+1} = PLQVC(V,E,f,F);
            % 12 Difference Indefinite Extended
            %       3x^2+xy+y^2     x<0    
            %       x^2+2xy+y^2     x>=0   
            V = [0 0; 0 1; 0 -1];
            E = [1 2 0; 1 3 0];
            f = [3 1 1 0 0 0;1 2 1 0 0 0];
            f = [zeros(size(f,1),10-size(f,2)),f];
            Cw= ones(size(f,1),1)* invW;
            f = f .* Cw;
            F = [1 2;2 1];
            P{end+1} = PLQVC(V,E,f,F);
            % 13 Simple cubic function
            %       -x^3-y^3+ x^2+ y^2      x<0    y<0
            %       -x^2+y^3+ x^2+2y^2      x<0    y>=0
            %       x^3+y^3+2x^2+2y^2       x>=0   y>=0
            %       x^3-y^3+2x^2+ y^2       x>=0   y<0
            V = [0 0; -1 0; 0 1; 1 0; 0 -1];
            E = [1 2 0; 1 3 0; 1 4 0; 1 5 0];
            f = [-1 0 0 -1 1 0 1 0 0 0;
                 -1 0 0  1 1 0 2 0 0 0;
                  1 0 0  1 2 0 2 0 0 0;
                  1 0 0 -1 2 0 1 0 0 0];
            Cw= ones(size(f,1),1)* invW;
            f = f .* Cw;
            F = [1 2;2 3;3 4;4 1];
            P{end+1} = PLQVC(V,E,f,F);
            % 14 Semi bounded cubic
            %       x^3+y^3  +y^2       0<x<1       0<=y<=1
            %       x^3+xy^2+y^3        x>=1        0<=y<=1
            %       infinity            otherwise
            V = [0 0; 1 0; 1 1; 0 1; 2 0; 2 1];
            E = [1 2 1; 3 4 1; 4 1 1; 2 3 1; 2 5 0; 3 6 0];
            f = [1 0 0 1 0 0 1 0 0 0;
                 1 0 1 1 0 0 0 0 0 0];      
            Cw= ones(size(f,1),1)* invW;
            f = f .* Cw;
            F = [1 0;1 0;1 0;1 2;2 0;0 2];
            P{end+1} = PLQVC(V,E,f,F);
            % 15 Semi bounded cubic (exactly same as above but slightly changed matrices)
            %       x^3+y^3  +y^2       0<x<1   0<=y<=1
            %       x^3+xy^2+y^3        x>=1    0<=y<=1
            %       infinity            otherwise
            V = [1 0; 0 0; 0 1; 1 1; 2 0; 2 1];
            E = [1 2 1; 2 3 1; 3 4 1; 4 1 1; 1 5 0; 4 6 0];
            f = [1 0 0 1 0 0 1 0 0 0;1 0 1 1 0 0 0 0 0 0];
            Cw= ones(size(f,1),1)* invW;
            f = f .* Cw;
            F = [0 1; 0 1; 0 1; 2 1; 2 0; 0 2];
            P{end+1} = PLQVC(V,E,f,F);
            % 16 half-star domain with 5 rays
            V = [0,0; -1,0; -1,-1; 0,1; 1,1; 1,0];
            E = [1,2,0; 1,3,0; 1,4,0; 1,5,0; 1,6,0];
            F = [0,1; 1,2; 2,3; 3,4; 4,0]; 
            f = zeros(4,1);
            P{end+1} = PLQVC(V,E,f,F);
            % 17 bounded domain with vertexes equispaced on the unit circle
            nv=10;N=(0:nv-1)';
            V = [cos(N*2*pi/nv), sin(N*2*pi/nv)];
            E = [(1:nv)', [(2:nv)';1], ones(nv,1)];
            F = [ones(nv,1), zeros(nv,1)];
            f = 0;
            P{end+1} = PLQVC(V,E,f,F);
            % 18 unbounded domain with 4 parallel rays
            V = [[(1:5)', zeros(5,1)];(1:4)',ones(4,1)];
            E = [1,2,1; 2,3,1; 3,4,1; 4,5,0;1,6,0;2,7,0;3,8,0;4,9,0];
            F = [1,0  ; 2,0  ;  3,0 ;   4,0;  0,1; 1,2 ; 2,3 ; 3,4 ];
            f = zeros(4,1);
            P{end+1} = PLQVC(V,E,f,F);
            % 19 unbounded grid of Tic-Tac-Toe
            V = [0,0; 1,0; -1,1; 0,1; 1,1; 2,1; -1,2; 0,2; 1,2; 2,2; 0,3; 1,3];
            E = [4,1,0; 5,2,0; 4,3,0; 4,5,1; 5, 6,0; 4,8,1; 5,9,1; 8,7,0; 8,9,1; 9,10,0; 8,11,0; 9,12,0];
            F = [2,1  ; 3,2  ; 1,4  ; 5,2  ; 6,3   ; 4,5  ; 5,6  ; 4,7  ; 8,5  ; 9,6   ; 7,8   ; 8,9   ];
            f = zeros(9,1);
            P{end+1} = PLQVC(V,E,f,F);
        end 
        function P = examples2()%degenerate cases; domain is always connected
           %P=PLQVC.examples; for i=1:length(P),plot(P{i});title(i+": "+P{i}.nf+" faces");pause;end
           %needle
           %1 needle function
            V = [0 0];
            E = [];
            f = [1 2 3];%x+2y+3
            F = [];           
            P{1}=PLQVC(V,E,f,F);
            %NEED TO FIX REPRESENTATION PROBLEM: NEEDLE, SEGMENT, RAY OK BUT LINE HAS 2 EDGES BUT 1 ROW IN f. SAME FOR
            %CHAIN OF EDGES. NEED TO REPRESENT CORRECTLY FUNCTIONS WITH DIMENSION 1 DOMAIN
            % 2 segment
            V = [0 0;1 0];
            E = [1 2 1];
            f = [1 0 0 0 0 0];%0.5 x^2
            F = [0 0];
            P{end+1} = PLQVC(V,E,f,F);
            % 3 ray
            V = [0 0;1 0];
            E = [1 2 0];
            f = [1 0 0 0 0 0];%x^2
            F = [0 0];
            P{end+1} = PLQVC(V,E,f,F);    
           % 4 line
            V = [0 0;1 0; -1 0];
            E = [1 2 0; 1 3 0];
            f = [1 0 0 0 0 0;1 0 0 0 0 0];%x^2
            F = [0 0;0 0];
            P{end+1} = PLQVC(V,E,f,F);               
            % 5 bounded slice boundary
            V = [0 0;1 0;0 1];
            E = [1 2 1;1 3 1];
            f = [1 0 0 0 0 0;0 0 1 0 0 0];
            F = [0 0;0 0];
            P{end+1} = PLQVC(V,E,f,F);
            % 6 half-bounded slice boundary
            V = [0 0;1 0;0 1];
            E = [1 2 0;1 3 1];
            f = [1 0 0 0 0 0;0 0 1 0 0 0];
            F = [0 0;0 0];
            P{end+1} = PLQVC(V,E,f,F);
            % 7 unbounded slice boundary
            V = [0 0;1 0;0 1];
            E = [1 2 0;1 3 0];
            f = [1 0 0 0 0 0;0 0 1 0 0 0];
            F = [0 0;0 0];
            P{end+1} = PLQVC(V,E,f,F);
            % 8 half-space
            V = [0 0;1 0;-1 0];
            E = [1 2 0;1 3 0];
            f = [1 0 0 0 0 0];%0.5 x^2
            F = [1 0;0 1];
            P{end+1} = PLQVC(V,E,f,F);
            % 9 3 parallel rays with dimension two domain
            V = [0 2; -0.5 2.5; 3 1.25; 2 1; 5 1.75; 5.5 3.25; 6 2.75; 2.5 1.75];
            E = [4 3 1; 5 3 1; 7 6 0; 1 2 0; 7 5 1; 1 4 1; 3 8 0];
            F = [1 0;0 2; 2 0; 0 1; 0 2; 1 0; 1 2];
            f = [0;0];        
            P{end+1} = PLQVC(V,E,f,F);
            % 10 chain of edges, domain is a staircase
            V = [0 0; 1 0; 1 1; 2 1; 2 2; 3 2];
            E = [1 2 1; 2 3 1; 3 4 1; 4 5 1; 5 6 1];
            F = zeros(5,2);
            f = [0 0 0 0 2 0 0 0 0 0; %x^2
                 0 0 0 0 0 0 2 0 0 1; %y^2 + 1
                 0 0 0 0 2 0 0 0 0 2; %x^2 + 2
                 0 0 0 0 0 0 2 0 0 6; %y^2 + 1
                 0 0 0 0 2 0 0 0 0 10];%x^2 + 2
            P{end+1} = PLQVC(V,E,f,F);
        end
        function P = examplesNonconvex()
            % 1 energy on nonconvex domain
            V = [0,0;1,1;2,1;0,2;2,2];
            E = [1,2,1; 2,3,1; 1,4,1; 3,5,1; 4,5,1];
            F = [1,0  ; 1,0  ; 0,1  ; 1,0  ; 0,1  ];
            f = [1,0,1,0,0,0];
            P{1} = PLQVC(V,E,f,F);   
            % 2 W epigraph of W function
            V = [-2,1;-1,0;0,1;1,0;2,1];
            E = [2,1,0;2,3,1;3,4,1;4,5,0];
            F = [0,1  ; 1,0 ; 1,0 ; 1,0 ];
            f = [1,0,1,0,0,0];
            P{end+1} = PLQVC(V,E,f,F);
            
        end
        function P = examplesDiscontinuous()
            % 1: constant function on each of 2 half-planes
            V = [0,0;1,0;-1,0];
            E = [1,2,0;1,3,0];
            F = [1,2;2,1];
            f = [0;1];
            P{1} = PLQVC(V,E,f,F);
            % 2: linear function on each of 2 half-planes
            V = [0,0;1,0;-1,0];
            E = [1,2,0;1,3,0];
            F = [1,2;2,1];
            f = [0 0 0;1 0 0];%0 on Face 1 (y>=0); x on Face 2 (y<=0)
            P{end+1} = PLQVC(V,E,f,F);
            % 3: quadratic function on each of 2 half-planes; continuous at 2 vertexes
            V = [0,0;1,0;-1,0];
            E = [1,2,0;1,3,0];
            F = [1,2;2,1];%x^2 on Face 1; y^2+x on Face 2
            f = [2 0 0 0 0 0; 0 0 2 1 0 0];
            P{end+1} = PLQVC(V,E,f,F);
            % 4: cubic function  on each of 2 half-planes; continuous at 3 vertexes
            V = [0,0;1,0;-1,0];
            E = [1,2,0;1,3,0];
            F = [1,2;2,1];%x^2 on Face 1; y^2+x on Face 2
            f = [6 0 0 0 0 0 0 -1 1 0; 6 0 0 0 0 0 0 -1 0 0];
            P{end+1} = PLQVC(V,E,f,F);
            
        end
  end
end % classdef

%     function [bool, flag] = checkDD() 
%         [bool, ~, ~,flag] = isDD(plq);
%     end
%     function [bool, flag] = checkConvexity()
%         [bool, ~,flag] = plqvc_isQuad2DConvex(plq);
%     end
%         end
