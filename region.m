classdef region
    % ineqs always stored as <= 0
    properties
        ineqs=symbolicFunction;
        nv;
        vx=sym.empty();
        vy=sym.empty();
        vars;
    end


     methods (Static)
         % ---- vertex-index guard ---------------------------------------------------------
         % HISTORY: this block also held probeAlong/probePerp, which turned a constraint SLOPE
         % into actual probe points for maxArray's tie-breaking. maxArray no longer probes at
         % all -- a sample cannot prove that one function dominates another, and taking it as
         % proof produced two silent wrong answers (read maxArray's header) -- so both helpers,
         % and maxFromPt/maxFromPts with them, have gone. The lesson they recorded is worth
         % keeping even though the code is not: a slope is a LOSSY encoding of a direction (a
         % vertical direction has no slope, the perpendicular of a horizontal one is vertical),
         % and a symbolic slope hides its own degeneracy from a syntactic `d == 0` test -- the
         % slope that broke f = x*y over conv{(0,0),(3,3),(1,2)} was
         %   2^(1/2)/4 + (4*2^(1/2) - 4)/(2*(4*2^(1/2) - 8)),
         % which IS zero but is not the symbol 0, so `-1/d` stayed unevaluated and only blew up
         % much later inside subsF's simplifyFraction. Test the geometry, never the syntax.

         function kk = probeVertexIndex (k, j, nv)
         % objective: the vertex index getNormalConeVertexQ's isZero fallback may safely index.
         %
         % Both halves of getNormalConeVertexQ pair the current vertex j with a NEIGHBOUR k
         % (j-1 for the first normal-cone half-plane, j+1 for the second). Each has an explicit
         % guard for k falling off the ends -- `if k < 1` / `if k > obj.nv` -- which probes from
         % vertex j instead. The isZero fallback that follows, however, re-probes from vertex k
         % UNGUARDED, so it indexed obj.vx(0) or obj.vx(nv+1) whenever that fallback fired at an
         % end vertex: the same "index used outside the guard that established it" shape already
         % corrected in getEdges, splitmax3 and poly2orderUnbounded. It is not hypothetical --
         % the second conjugation of the two-face unit square (f = x*y) hits it on the half-lens
         % {(s1+s2)^2 <= 4*s1, (s1+s2)^2 <= 4*s2, s2 <= s1}, nv = 2, with j = 2 and k = 3,
         % raising MATLAB:badsubscript from inside functionNDomain.conjugateOfPiecePoly.
         %
         % Out of range, fall back to j -- the very vertex the enclosing guard already chose.
         % The fallback still probes along a DIFFERENT constraint than the guard did (ineqs(j)
         % rather than ineqs(j+1)), which is the point of running it at all: it needs a feasible
         % point where eq does not vanish.
             kk = k;
             if kk < 1 || kk > nv
                 kk = j;
             end
         end

         function eq = coneNormalAt (obj, cI, at, opp, s1, s2)
         % One half-plane of an edge's normal cone: the line through vertex `at`, perpendicular to
         % constraint cI's tangent THERE, oriented so that vertex `opp` satisfies it.
         %
         % Factored out of getNormalConeEdgeQ/Q3, which each repeat it twice, so that
         % getNormalConeEdgeQE can build the same object from an explicit edge list. The tangent
         % is taken at the vertex rather than from the chord between the vertices -- that is the
         % whole difference between these routines and the polyhedral getNormalConeEdge, and it is
         % what a curved edge needs.
             slope = obj.slopeIneq(cI, [obj.vx(at), obj.vy(at)]);
             pslope = -1/slope;
             if pslope == -inf
                 pslope = inf;
             end
             if pslope ~= inf
                 q = obj.yIntercept(at, pslope);
                 eq = s2 - pslope*s1 - q;
             else
                 eq = s1 - obj.vx(at);
             end
             if isAlways(subs(eq, [s1,s2], [obj.vx(opp), obj.vy(opp)]) > 0)
                 eq = -eq;
             end
         end

         % ---- LP certificates ------------------------------------------------------------
         % Two of this class's set operations used to decide a GLOBAL question -- "is this
         % constraint implied by the others?" (simplifyUnboundedRegion) and "is the union of
         % these two regions exactly the intersection of what is left after deleting their
         % shared facet?" (merge) -- with a LOCAL syntactic proxy: does the constraint pass
         % through a finite vertex, do the facet's endpoints look convex. Both proxies delete
         % constraints that carry weight, so a region ends up claiming territory that was never
         % its own, with the wrong function value on it (SUPPORT_MATRIX.md section 1.2).
         %
         % Both questions are the same primitive: maximize a linear form over a polyhedron.
         % That is an LP, it is exact, and it decides unboundedness and infeasibility as
         % first-class answers rather than as failures -- which matters here, since these
         % regions are routinely unbounded. Every caller is written so that an UNDECIDED answer
         % means "keep the constraint" / "refuse the merge": over-describing a region is
         % harmless, under-describing it is the bug.

         function [val, st] = maxLinear (A, b, c)
         % max c*[x;y] over {A*[x;y] <= b}.
         %   st =  0 : optimal, val is the maximum
         %   st =  1 : unbounded above (val = inf)
         %   st = -1 : the constraint set is infeasible (val = -inf)
         %   st =  2 : undecided -- the solver did not converge; callers must treat this as
         %             "no certificate" and fall back to the conservative branch.
             c = double(c(:));
             if isempty(A)
                 if all(c == 0), val = 0; st = 0; else, val = inf; st = 1; end
                 return
             end
             persistent opts
             if isempty(opts)
                 opts = optimoptions('linprog', 'Display', 'none');
             end
             ws = warning('off', 'all');
             try
                 [~, fval, ef] = linprog(-c, double(A), double(b), [], [], [], [], opts);
             catch
                 ef = 0; fval = [];
             end
             warning(ws);
             switch ef
                 case 1,  val = -fval; st = 0;
                 case -3, val = inf;   st = 1;
                 case -2, val = -inf;  st = -1;
                 otherwise, val = NaN; st = 2;
             end
         end

         function l = impliedBy (Ac, bc, A, b)
         % Is every constraint Ac(i,:)*[x;y] <= bc(i) satisfied at every point of the
         % polyhedron {A*[x;y] <= b}? An infeasible polyhedron satisfies all of them
         % vacuously. Answers false whenever it cannot prove true.
             l = false;
             for i = 1:size(Ac,1)
                 [val, st] = region.maxLinear(A, b, Ac(i,:));
                 if st == -1
                     l = true;    % {A x <= b} is empty: nothing left to violate anything
                     return
                 end
                 if st ~= 0
                     return       % unbounded above, or undecided: no certificate
                 end
                 if val > bc(i) + 1.0d-9 * max(1, abs(bc(i)))
                     return
                 end
             end
             l = true;
         end

         function out = mergeTally (reason)
         % INSTRUMENTATION for Step 3's cell blow-up. `region.merge` is what is supposed to
         % collapse adjacent same-valued cells; when it refuses, the cell count grows and
         % functionNDomain.maxOfList's cost grows with it. Counting the refusals BY REASON is
         % what turns "merge never merges" into a specific gate to fix.
         %
         %   region.mergeTally('reset')      clear the counters
         %   s = region.mergeTally('get')    the counts, as a struct
         %   region.mergeTally('<reason>')   record one occurrence
         %
         % Always on: one struct-field increment per merge attempt is nothing beside the
         % symbolic comparisons merge already does, and a counter that has to be switched on
         % is a counter that is off when the interesting run happens.
             persistent T
             if isempty(T), T = struct(); end
             out = [];
             switch reason
                 case 'reset', T = struct();
                 case 'get',   out = T;
                 otherwise
                     if isfield(T, reason), T.(reason) = T.(reason) + 1;
                     else,                  T.(reason) = 1;
                     end
             end
         end
     end

%  57 methods
     methods

         % check when multiple regions are created by ineqs
         function obj = region(fs, vars) %, not)
        
            if nargin == 0
                return
            end
            m = size(fs,1);
            n = size(fs,2);
           
            nineq = 0;
            for i = 1:m
              for j = 1:n
                  f = symbolicFunction(fs(i,j));
             
                if (~f.isZero)

                  nineq = nineq + 1;
                  obj.ineqs(nineq) = f;
                end  
              end
            end
            
            % INSTRUMENTATION (CCA2_TRACE_BIGNUM): name whoever builds a region out of
            % DOUBLE-DERIVED numbers. A 15-digit integer literal in a constraint means some
            % exact quantity went through a double on the way here -- `sym` of a double is an
            % exact binary rational with a 2^53 denominator, and a few multiplications later
            % the coefficients are 145 digits long, every symbolic call is slow, and `isAlways`
            % can no longer prove that two forms of the same number are equal.
            %
            % That is not hypothetical: this is what found it. domain.mE / domain.cE were
            % DOUBLE arrays, so an exact edge slope became 0.6 and an exact zero y-intercept
            % became -9.06e-72; plq_1p.conjugateFunction's nCE == 1 branch builds a whole
            % conjugate out of those two numbers, and Step 3's merge then stopped working
            % because nothing could be compared. See domain.m's property comment.
            if ~isempty(getenv('CCA2_TRACE_BIGNUM'))
                for iBN = 1:nineq
                    if ~isempty(regexp(char(obj.ineqs(iBN).f), '\d{15,}', 'once'))
                        stBN = dbstack;
                        fprintf('BIGNUM region ctor ineq %d:', iBN);
                        for kBN = 2:min(9, numel(stBN))
                            fprintf(' %s:%d <-', stBN(kBN).name, stBN(kBN).line);
                        end
                        fprintf('\n');
                        break
                    end
                end
            end
            obj.vars = vars;
            obj = obj.normalize1;
            obj = obj.unique;
            %obj.ineqs.printL
            
            obj = obj.getVertices  ;
            
            
            if obj.nv == 0
                obj = region.empty;
            end
            % put simplify in here
         end

         function [l, pos] = isVertexIrrational (obj)
             l = false;
             n = 0;
             pos = [];
             for i = 1:obj.nv
                 if symFunType(obj.vx(i)) == 'plus' | symFunType(obj.vy(i)) == 'plus'
                     l = true;
                     n = n + 1;
                     pos(n) = i;
                 end
             end
         end
         function obj = regionWPts(obj, vx, vy, x, y)
             n = 0;
             for i = 1: size(vx,2)-1
                 if vx(i) == vx(i+1)
                   n = n + 1;
                   ineq(n) = x-vx(i)  ;
                 else
                   m = (vy(i+1)-vy(i))/(vx(i+1)-vx(i));
                   c = vy(i)-m*vx(i);
                   n = n + 1;
                   ineq(n) = y-m*x-c
                 end 
             end
             if vx(1) == vx(end)
               n = n + 1;
               ineq(n) = x-vx(1);  
             else
               m = (vy(end)-vy(1))/(vx(end)-vx(1));
               c = vy(1)-m*vx(1);
               n = n + 1;
               ineq(n) = y-m*x-c
             end
             meanx = sum(vx)/3;
             meany = sum(vy)/3;
             for i = 1:n
                 if double(subs(ineq(i),[x,y],[meanx,meany])) > 0
                     ineq(i) = -ineq(i);
                 end

             end
             obj = region(ineq,[x,y]);
         end

         function f = plus(obj1,obj2)
            % plus is region INTERSECTION (it unions the inequality lists), so an empty operand
            % gives an empty result. Without this guard `obj1.ineqs` on a region.empty() array
            % returns a 0-element comma-separated list and the indexing below errors instead.
            if isempty(obj1) || isempty(obj2)
                f = region.empty();
                return
            end
            l = [];
            for i = 1:size(obj1.ineqs,2)
              l = [l,obj1.ineqs(i).f];
            end
            for i = 1:size(obj2.ineqs,2)
               lunique = true;
               for j = 1:size(obj1.ineqs,2)
                   %obj1.ineqs(j)
                   %obj2.ineqs(i)
                 if (obj1.ineqs(j) == obj2.ineqs(i))
                     lunique = false;
                     break;
                 end
                 if (obj1.ineqs(j).f == -obj2.ineqs(i).f)
                     f = region.empty;
                     return;
                 end
               end
               if lunique
                   l = [l,obj2.ineqs(i).f];
               end
            end 
            %l
            f = region(l,obj1.vars);
         end

         function res =eq2(obj1,obj2)
             res = false;
             if ~eqVertices(obj1,obj2)
                 return
             end
             if size(obj1.ineqs,2) ~= size(obj2.ineqs,2)
                 return
             end
             for i = 1:size(obj1.ineqs,2)
                 marked(i) = false;
             end
             for i = 1:size(obj1.ineqs,2)
                f1 = obj1.ineqs(i).normalize(obj1.vars);
                 
                for j = 1:size(obj2.ineqs,2)
                    if marked(j)
                     continue
                    end
                    f2 = obj2.ineqs(j).normalize(obj2.vars);
                  if isAlways(f1 == f2)
                    marked(i) = true;
                  end
                end
             end
             
             for i = 1:size(obj1.ineqs,2)
                if ~ marked(i)
                    return
                end
             end
             res = true;
         end

        
         function obj = poly2orderUnbounded(obj)
            % disp('in poly2orderUnbounded')
             if obj.nv == 1
                 return
             end
           % rad = cart2pol(obj.vx(1:obj.nv), obj.vy(1:obj.nv))
           % radWrapped = mod(rad,2*pi)
           % radWrapped(radWrapped==0 & rad>0) = 2*pi
           % [~, sortIdx] = sort(radWrapped, 'descend') 
           % obj.vx(sortIdx) = obj.vx(1:obj.nv);
           % obj.vy(sortIdx) = obj.vy(1:obj.nv);

           x = obj.vx(1:obj.nv);
           y = obj.vy(1:obj.nv);
           cx = mean(double(x));
           cy = mean(double(y));
           a = atan2(y - cy, x - cx);
           [~, order] = sort(a);
           obj.vx(1:obj.nv) = x(order);
           obj.vy(1:obj.nv) = y(order);
           
           if obj.nv == size(obj.ineqs,2)
               return
           end
           %obj.vx(1:obj.nv)
           %obj.vy(1:obj.nv)
           % Search only up to obj.nv (not size(obj.ineqs,2), which indexes obj.vx/vy here and can
           % legitimately exceed obj.nv: obj.nv ~= size(obj.ineqs,2) -- the caller's own signal for
           % "call poly2orderUnbounded, not poly2order" -- can ALSO happen on an ORDINARY BOUNDED
           % region with redundant/non-essential inequality rows (tangent at an existing vertex,
           % not cutting a new one), not only on a genuinely unbounded region with ray edges. See
           % HISTORY below for the crash this produced and why "not found" is handled as a no-op
           % rather than an error.
           found = false;
           for i = 1:obj.nv
             edges = obj.getEdges(obj.vx(i),obj.vy(i));

             % HISTORY: edges(1)/edges(2) were indexed unconditionally, which THREW
             % ('Index exceeds the number of array elements') on the simplest unbounded region
             % there is -- a 2-constraint wedge such as the first quadrant {-x<=0, -y<=0}.
             % getEdges reports the constraints of obj.ineqs ACTIVE at a point, and the vertex
             % list of an unbounded region also contains the box-clip vertices this file's own
             % getVertices appends: the corners (+-intmax,+-intmax) and the per-vertex
             % projections (vx(i),+-intmax), (+-intmax,vy(i)). Those lie on the IMPLICIT
             % +-intmax box, which is not among obj.ineqs, so they have fewer than two active
             % constraints -- (intmax,intmax) on the first quadrant has none at all.
             %
             % Such a vertex simply cannot be the "seam" this loop looks for (a vertex adjacent
             % to a ray edge of obj.ineqs), so skipping it is not a workaround, it is the correct
             % reading. Note the ANGLE SORT above has already put the whole list in boundary
             % cyclic order, which is what getNormalConeVertex actually consumes -- it walks
             % consecutive pairs j, j+1 and wraps at nv. This loop only chooses where the cycle
             % STARTS, so failing to find a seam costs nothing.
             if numel(edges) < 2
                 continue
             end

             [nv1, vx, vy] = obj.vertexOfEdge(edges(1));
             [nv2, vx, vy] = obj.vertexOfEdge(edges(2));

             if nv1 == 1 | nv2 == 1
                 found = true;
                 break
             end

           end
           if ~found
               % No vertex here is adjacent to a ray edge: this region is not actually unbounded
               % (see comment above) -- there is no "seam" to rotate the vertex cycle around, so
               % leave the angle-sorted order from above as-is instead of reading past the end of
               % obj.vx/vy (obj.vx(i+1) below, with i==obj.nv, previously threw
               % 'MATLAB:badsubscript' -- found via cplqAdapterTest/testMaxMultiRegion's
               % testMax, a plq.biconjugateF call on a genuinely multi-triangle input).
               return
           end
           %obj.print
           %i
           %obj.vx(i+1),obj.vy(i+1)
           iNext = mod(i, obj.nv) + 1;   % wrap: the found vertex can legitimately be the LAST one
                                         % in angular order, same wraparound fix already applied
                                         % elsewhere in this codebase (see maxQuaPar.m HISTORY).
           edges = obj.getEdges(obj.vx(iNext),obj.vy(iNext));

           if numel(edges) >= 2           % same box-clip guard as the search loop above
               [nv1, vx, vy] = obj.vertexOfEdge(edges(1));
               [nv2, vx, vy] = obj.vertexOfEdge(edges(2));
               if nv1 == 1 | nv2 == 1
                 i = iNext;
               end
           end

           vx = obj.vx;
           vy = obj.vy;
           obj.vx(1:obj.nv-i+1) = vx(i:obj.nv); 
           obj.vy(1:obj.nv-i+1) = vy(i:obj.nv);
           obj.vx(obj.nv-i+2:obj.nv) = vx(1:i-1); 
           obj.vy(obj.nv-i+2:obj.nv) = vy(1:i-1); 
           
         end

         function obj = poly2order(obj)
            % obj.print
             vx(1) = obj.vx(1);
             vy(1) = obj.vy(1);
             for i = 1:obj.nv
               lineqs(i) = false ;
             end
             
             for i = 1:obj.nv-1
                 
                 
                 ineqs = obj.ineqThroughVertex1 (vx(i), vy(i));
                 for j = 1:2
                     if lineqs(ineqs(j))
                         continue
                     end
                     break;
                 end
                 
                 lineqs(ineqs(j)) = true;
                 v = obj.getEndpoints (ineqs(j));
                 if (v(1,1) == vx(i) & v(1,2) == vy(i))
                     vx(i+1) = v(2,1);
                     vy(i+1) = v(2,2);
                 else
                     vx(i+1) = v(1,1);
                     vy(i+1) = v(1,2);
                 end


             end
             obj.vx = vx;
             obj.vy = vy;
             cx = sum(obj.vx) / obj.nv;
             cy = sum(obj.vy) / obj.nv;
%obj.print
             theta1 = atan((vy(1)-cy)/(vx(1)-cx));
             if theta1 < 0
                 theta1 = theta1 + pi;
             end
             theta2 = atan((vy(2)-cy)/(vx(2)-cx));
             if theta2 < 0
                 theta2 = theta2 + pi;
             end
             if theta1 > theta2
                 return;
             end
             for i = 1:obj.nv
                 obj.vx(i) = vx(obj.nv-i+1);
                 obj.vy(i) = vy(obj.nv-i+1);
             end
         end

         function ineqs = ineqThroughVertex (obj,j)
           ineqs = [];
           for i = 1:size(obj.ineqs,2)
               % change this to exact
               if abs(double(subs(obj.ineqs(i).f,obj.vars,[obj.vx(j),obj.vy(j)]))) <= 1.0e-6
                   ineqs = [ineqs,i];
               end
           end
         end

         function ineqs = ineqThroughVertex1 (obj,vx, vy)
           ineqs = [];
           for i = 1:size(obj.ineqs,2)
               % change this to exact
               if abs(double(subs(obj.ineqs(i).f,obj.vars,[vx,vy]))) <= 1.0e-6
                   ineqs = [ineqs,i];
               end
           end
         end

         function v = getEndpoints (obj, i)
             n = 0;
             for j = 1:obj.nv
               if abs(double(subs(obj.ineqs(i).f,obj.vars,[obj.vx(j),obj.vy(j)]))) <= 1.0e-6
                   n = n + 1;
                   v(n,1) = obj.vx(j);
                   v(n,2) = obj.vy(j);
               end
             end
            
             v = sortrows(v);
         end

         function n = getN(obj)
             n = 0;
             for i = 1: obj.nv
                 if abs(obj.vx(i)) == intmax
                     continue
                 end
                 if abs(obj.vy(i)) == intmax
                     continue
                 end
                 n = n + 1;
             end
         end

         function l = in(obj1,obj2)
             l = false;
             for i = 1:obj2.nv
                 if ~obj1.ptFeasible (obj1.vars,[obj2.vx(i),obj2.vy(i)])
                     return
                 end
             end
             l = true;
         end

         function f = minus(obj1,obj2)
             
             if (obj1 == obj2)
                 f = [region.empty];
                 return
             end     
            l = [];
            % HISTORY: v1/v2 used to be fixed-shape 3D arrays, one 1x2 "slab"
            % per edge, on the assumption every edge has exactly 2 endpoints.
            % That's only true for a bounded polygon's edges (segments
            % between 2 vertices) -- an unbounded region can have a "ray"
            % edge with just 1 finite endpoint (the placeholder "vertex at
            % infinity" scheme in region.getVertices doesn't guarantee a
            % matching point on every arbitrary-slope unbounded edge), which
            % crashed the fixed-shape assignment below (testOpenconvex).
            % Cell arrays hold the variable-length endpoint lists instead;
            % the comparison further down now compares sizes too, treating a
            % genuine shape mismatch as "not the same edge" rather than
            % crashing.
            v1 = {};
            for i = 1:size(obj1.ineqs,2)
                v1{i} = obj1.getEndpoints(i);
                l = [l,obj1.ineqs(i).f];
            end
            l2 = [];
            rm = [];
            v2 = {};
            for i = 1:size(obj2.ineqs,2)
              v2{i} = obj2.getEndpoints(i);
              ladd = true;
              lsub = false;
              for j = 1:size(obj1.ineqs,2)
                %  i,j
                  %obj1.ineqs(j).f
                if (obj2.ineqs(i) == obj1.ineqs(j))
                    lv = isequal(size(v1{j}),size(v2{i}));
                    if lv
                        for i1 = 1:size(v1{j},1)
                            for i2 = 1:size(v1{j},2)
                                lv = lv & (v1{j}(i1,i2) == v2{i}(i1,i2));
                            end
                        end
                    end
                    if lv
                      rm = [rm,j];
                      ladd = false;
                    else
                      ladd = false;
                      lsub = true;
                    end
                %f = [f0.simplify];

                    break;
                end
              end
              %ladd
              %obj2.ineqs(i) == obj1.ineqs(j)
              %obj2.ineqs(i).f
              if ladd
                  l2 = [l2,obj2.ineqs(i).f];
                 
              end
              % if lsub
              %     l2 = [l2,-obj2.ineqs(i).f];
              % 
              % end
            end
            
            %l2
            l(rm) = [];
            % l
            % l2
            for i = 1: size(l2,2)
                l = [l,-l2(i)];
            end
            % l
            % size(l)
            if size(l,2) <= 2
                f = [region.empty];
                return
            end
            % disp('here')
            f0 = region(l,obj1.vars);
            % f0.print
            %f0.getN
            if ~isempty (f0) & f0.getN > 2 
                f = [f0];
                return
            end
            % disp('here2')
            rL = region.empty();
            
            % implement this with ripple 1/0 from 3 to n
            for i = 1:size(l,2)
                
                for j = i+1:size(l,2)
                    for k = j+1:size(l,2)
                        l0 = [l(i),l(j),l(k)];
                        f0 = region(l0,obj1.vars);
                        % i,j,k
                         % f0.print
                         
                        if isempty(f0)
                            continue;
                        end
                        if f0.getN < 3
                             continue;
                        end
                        if ~ obj1.in(f0)
                            continue
                        end
                        f0 = f0.simplify;
                        if f0.nv ~= size(f0.ineqs,2)
                            continue
                        end
                        rL = [rL,f0]; 
                    end

                end

            end
            % disp('triple')
            % for i = 1:size(rL,2)
            %     rL(i).print
            % end
            for i = 1:size(l,2)
                for j = i+1:size(l,2)
                    for k = j+1:size(l,2)
                       for il = k+1:size(l,2)
                        l0 = [l(i),l(j),l(k),l(il)];
                        f0 = region(l0,obj1.vars);
                        if isempty(f0)
                            continue;
                        end
                                   if f0.getN < 3
                             continue;
                        end
                        if ~ obj1.in(f0)
                            continue
                        end
                        f0 = f0.simplify;
                        if f0.nv ~= size(f0.ineqs,2)
                            continue
                        end
                        rL = [rL,f0]; 
                       end 
                    end

                end

            end

            for i = 1:size(l,2)
                for j = i+1:size(l,2)
                    for k = j+1:size(l,2)
                       for il = k+1:size(l,2)
                         for im = il+1:size(l,2)
               
                           l0 = [l(i),l(j),l(k),l(il),l(im)];
                           f0 = region(l0,obj1.vars);
                           if isempty(f0)
                             continue;
                           end
                           if f0.getN < 3
                             continue;
                           end
                           if ~ obj1.in(f0)
                            continue
                           end
                        
                           f0 = f0.simplify;
                        if f0.nv ~= size(f0.ineqs,2)
                            continue
                        end
                        rL = [rL,f0]; 
                         end
                       end 
                    end

                end

            end
            if (size(l,2) > 6)
                disp('minus to implement ')
                size(l,2)
                l
                obj1.print
                obj2.print
            end
            %rL

            if size(rL,2) > 0
                f = rL.uniqueL;
            else
                f = [];
                return
            end
            rm = [];
            for i = 1:size(f,2)
             if f(i) == obj1
                 rm = [rm,i];
             end
            end
            
            
            f(rm) = [];
             % disp('finak')
             % for i = 1:size(f,2)
             % f(i).print
             % end
            return

            f0 = region(l,obj1.vars);
            disp('in minus')
            l
            f0.print
            if isempty(f0)
                f = [f0];
                return
            end
            if f0.nv >= 3
              f = [f0.simplify];
              return
            else  
              f1 = f0.divideRegions(obj1);
              f = [];
              for i = 1: size(f1,2)
                  if f1(i).nv <= 2
                      continue
                  end
                  
                  f0 = f1(i).simplify;
                  if f0 == obj1
                      continue
                  end
                  f = [f,f0];
              end
            end

            
            
            
         end
         
         
         
         function l = eqVertices(obj1,obj2)
             l = false;
             %obj1.nv
             if obj1.nv ~= obj2.nv
                return
             end
             for i = 1:obj1.nv
                 l1(i) = false;
                 l2(i) = false;
             end

             for i = 1:obj1.nv
               for j = 1:obj2.nv
                   if l2(j)
                       continue
                   end
                   if (abs(obj1.vx(i) - obj2.vx(j))<1.0d-15) & (abs(obj1.vy(i) - obj2.vy(j)) < 1.0d-15)
                       l1(i) = true;
                       l2(j) = true;
                       break;
                   end
               end
             end
             
             if (all(l1)==true )
                 l = true;
             else
                 l = false;
             end
         end

         function l = eq(obj1,obj2)
             l = false;
             
             if (size(obj1.ineqs,2)~=size(obj2.ineqs,2))
                 return;
             end
             
             if ~ obj1.eqVertices(obj2)
                 return
             end
             l = true;
         end

         function obj = unique(obj)
             n = 0;
             
            duplicates = [];
            for i = 1:size(obj.ineqs,2)
                for j = i+1:size(obj.ineqs,2)
                    if (obj.ineqs(i) == obj.ineqs(j))
                        n = n + 1;
                        duplicates(n) = j;
                        lMark(j) = 1;
                    end
                end

            end
           
            obj.ineqs(duplicates) = [];
         end


         function obj = uniqueL(obj)
             n = 0;
             
            duplicates = [];
            for i = 1:size(obj,2)
                for j = i+1:size(obj,2)
                    if (obj(i) == obj(j))
                        n = n + 1;
                        duplicates(n) = j;
                        
                    end
                end

            end
           
            obj(duplicates) = [];
         end
         function m = slope (obj,i,j)
          m = (obj.vy(i)-obj.vy(j))/(obj.vx(i)-obj.vx(j));
         end

         function m = slopeIneq(obj,i,pt)
             if obj.ineqs(i).isLinear
                 c = obj.ineqs(i).getLinearCoeffs (obj.vars);
                 
                 if c(2) == 0
                     m = inf;
                 else
                     m = -c(1)/c(2);
                 end
             else
                  vars =  obj.vars;
                  drx1 = obj.ineqs(i).dfdx(vars(1));
                  drx2 = obj.ineqs(i).dfdx(vars(2));
             %     subs(drx2.f,vars,pt)
                  if abs(subs(drx2.f,vars,pt)) < 1.0d-6
                      
                      m = intmax;
                  else
                      m0 = -drx1.f/drx2.f;
                      
                      m = subs(m0,vars,pt);
                  end
             
                 
             end
         end

          function c = yIntercept (obj,i,m)
          c = obj.vy(i)-m*obj.vx(i);   
          end

         function set = ithSet (obj, logical_indices)
           set = true;
           if all(logical_indices == false)
             set = false;
           end
                   
           for j = 1:size(logical_indices,2)
              if (logical_indices(j))
                set = set & obj.ineqs(j).f<=0;
              end
           end
            
         end
         
       

         function print(obj)
             if isempty(obj)
                 disp("Empty region")
                 return
             end
             disp("Variables")
             obj.vars
             disp(["nVertices = ", num2str(obj.nv)]);
             fprintf("vx =  ")
             fprintf("%f  ", obj.vx);
             fprintf("\n")
             fprintf("vy =  ")
             fprintf("%f  ", obj.vy);
             fprintf("\n\n")
               disp("Intersection of following ineqs")
               obj.ineqs.printLIneq;
         end

         function printLatex(obj)
             disp("Variables")
             fprintf("\\[")
             for i = 1: size(obj.vars,2)
               fprintf(char(obj.vars(i)) );
               if i == size(obj.vars,2)
                 break;
               end
               fprintf(",");
             end
             n = 0;
             for i = 1: obj.nv
                 if (abs(obj.vx(i)) == intmax)
                     continue
                 end
                 if (abs(obj.vy(i)) == intmax)
                     continue
                 end
                 n = n + 1;
             end
             fprintf("\\]\n\\[\\text{Number of Vertices = }")
             fprintf(num2str(n));
             
             fprintf("\\]\n\\[v =  ")
             for i = 1: obj.nv
                 if (abs(obj.vx(i)) == intmax)
                     continue
                 end
                 if (abs(obj.vy(i)) == intmax)
                     continue
                 end
             fprintf("(" + obj.vx(i) + ","+obj.vy(i)+")");
             if i == obj.nv
                 break;
             end
             fprintf(",");
             end
             fprintf("\\]\n")
             
               disp("Intersection of following ineqs")
               obj.ineqs.printLIneqLatex;
         end

         
         function printMaple(obj)
             obj = obj.subsF;
             
             obj.ineqs.printLIneqM;
             
         end

         function fprint(obj, uNo)
             fprintf(uNo, num2str(obj.nv)+"\n");
             for i = 1:obj.nv
               fprintf(uNo, num2str(obj.vx(i)) + "  " + num2str(obj.vy(i)) + "\n")  
             end
             obj.ineqs.fprintLIneq(uNo);
         end

          function plot (obj)
             
             l1 = min(min(obj.vx),min(obj.vy));
             if l1 < -25
                 l1 = -25;
             end
             l2 = max(max(obj.vx),max(obj.vy));
             if l2 > 20
                 l2 = 20;
             end
             l1 = -25;
             l2 = 20;
           obj.ineqs.plotLIneq (obj.vars, [l1,l2])   ;

         end


         function [vx,vy] = plotRegionC (obj, textR, c)
            limitsx = [-25,17];
            limitsy = [-15,20];
            
            pts = 75;
            colors = ['b', 'r', 'g', 'm', 'c', 'y'];
            stepx = (limitsx(2)-limitsx(1))/pts;
            stepy = (limitsy(2)-limitsy(1))/pts;
            n = 0;
            vx = [];
            vy = [];
            ci = limitsx(1);
            for i = 1:pts
                cj = limitsy(1);
                for j = 1:pts
                  if obj.ptFeasible (obj.vars,[ci,cj]);
                      n = n+1;
                      vx(n) = ci;
                      vy(n) = cj;
                  end
                  cj = cj+stepy;    
                end
                ci = ci+stepx;    
            end
            %c = colors(1+mod(obj.getGlobalParameter,6))
            if n == 0
                disp ('region not displayed')
            end
            avx = sum(vx)/n;
            avy = sum(vy)/n;
            m = 0;
            for i = 1:n
                if (abs(vx(i) -avx)<0.1 & abs(vy(i) -avy)<0.1)
                    continue
                end    
                m = m+1;
                vx1(m) = vx(i);
                vy1(m) = vy(i);
            end
            %text = "R";
            text(avx,avy,textR,'FontSize',12, 'FontWeight', 'bold','Color', 'k')
            if m==0
                fill(vx, vy, c, 'FaceAlpha', 0.9, 'EdgeColor', c);
            else
                fill(vx1, vy1, c, 'FaceAlpha', 0.9, 'EdgeColor', c);
            end
            
         end

         function [vx,vy] = plotRegion (obj, textR)
            limitsx = [-6,6];
            limitsy = [-15,20];
            limitsx = [-25,17];
            limitsy = [-6,11];
            
            pts = 75;
            colors = ['b', 'r', 'g', 'm', 'c', 'y'];
            stepx = (limitsx(2)-limitsx(1))/pts;
            stepy = (limitsy(2)-limitsy(1))/pts;
            n = 0;
            vx = [];
            vy = [];
            ci = limitsx(1);
            for i = 1:pts
                cj = limitsy(1);
                for j = 1:pts
                  
                  if obj.ptFeasible (obj.vars,[ci,cj]);
                      n = n+1;
                      vx(n) = ci;
                      vy(n) = cj;
                  end
                  cj = cj+stepy;    
                end
                ci = ci+stepx;    
            end
            c = colors(1+mod(obj.getGlobalParameter,6))
            if n == 0
                obj.print 
                obj.ptFeasible (obj.vars,[2,-9])
                disp ('region not displayed')
            end
            avx = sum(vx)/n;
            avy = sum(vy)/n;
            m = 0;
            for i = 1:n
                if (abs(vx(i) -avx)<0.1 & abs(vy(i) -avy)<0.1)
                    continue
                end    
                m = m+1;
                vx1(m) = vx(i);
                vy1(m) = vy(i);
            end
            %text = "R";
            text(avx,avy,textR,'FontSize',12, 'FontWeight', 'bold','Color', 'k')
            if m==0
                fill(vx, vy, c, 'FaceAlpha', 0.9, 'EdgeColor', c);
            else
                fill(vx1, vy1, c, 'FaceAlpha', 0.9, 'EdgeColor', c);
            end
            
         end

        
         function obj = subsF(obj)
             x= sym('x');
             y= sym('y');
             for i = 1:size(obj.ineqs,2)
                 obj.ineqs(i) = obj.ineqs(i).subsF(obj.vars,[x,y]);
             end
         end
         function m = slopes (obj)
           n = 0 ;  
           for i = 1:size(obj.ineqs,2)
              if obj.ineqs(i).isLinear
                  c = obj.ineqs(i).getLinearCoeffs (obj.vars);
                  n = n + 1;
                  m(n) = -c(1)/c(2);
              end
           end
         end

         
         function m = slopeAtVertex (obj, pi, pt)
            % pt
           n = 0 ;  
           %pi
           %pt
           %obj.ineqs.printL
           for j = 1:size(pi,2)
              i = pi(j); 
              if obj.ineqs(i).isLinear
                  c = obj.ineqs(i).getLinearCoeffs (obj.vars);
                  n = n + 1;
                  m(n) = -c(1)/c(2);
              else
                  %check if point lies on curve
                  %
            %      obj.ineqs(i).print
                  vars =  obj.vars;
                  drx1 = obj.ineqs(i).dfdx(vars(1));
                  drx2 = obj.ineqs(i).dfdx(vars(2));
             %     subs(drx2.f,vars,pt)
                  if abs(subs(drx2.f,vars,pt)) < 1.0d-6
                      n = n + 1;
                      m(n) = intmax;
                  else
                      m0 = -drx1.f/drx2.f;
                      n = n + 1;
                      m(n) = subs(m0,vars,pt);
                  end
              end
           end
         end

         
         function m = slopes2 (obj)
         % One representative boundary slope per constraint. NOTE: its only caller was
         % maxArray's tie-break probing, which no longer exists (see maxArray's header), so
         % nothing in the toolbox calls this today. It is kept because the note below records a
         % measured performance trap that any future caller needs to know about.
         %
         % A CURVED constraint has no single slope, so its tangent is taken at a
         % vertex of this region that lies on it -- and a curved constraint need not have one.
         % That case was unreachable while simplifyUnboundedRegion deleted every constraint
         % missing a finite vertex; now that deletion requires a real redundancy certificate
         % (see redundantSubset), such a constraint survives, `pt` was left unassigned, and the
         % substitution below failed with 'Unrecognized function or variable pt'.
         %
         % Do NOT simply drop such a constraint. m is read pairwise by maxArray to build
         % bisector probe DIRECTIONS, and a shorter list means fewer probes, so maxArray
         % returns undecided more often, maxEqDom falls through to splitmax3, and every
         % undecided region splits in two -- which compounds round over round. Take the
         % tangent at this region's own finite-vertex centroid instead: an interior point, so
         % the direction it yields is a local one, and it always exists.
           n = 0 ;
           m = [];        % a region can consist entirely of skipped constraints
           % l = false;

           for i = 1:size(obj.ineqs,2)
              if obj.ineqs(i).isLinear
                  c = obj.ineqs(i).getLinearCoeffs (obj.vars);
                  n = n + 1;
                  m(n) = -c(1)/c(2);
              else
                  vars =  obj.vars;

                  pt = [];
                  for j = 1:obj.nv
                      if abs(obj.ineqs(i).subsF(vars,[obj.vx(j),obj.vy(j)]).f) < 1.0d-8
                          pt = [obj.vx(j),obj.vy(j)];
                          break;
                      end
                  end
                  if isempty(pt)
                      [nPc, pxc, pyc] = obj.finiteVertices;
                      if nPc == 0
                          continue          % nothing finite to anchor a direction to
                      end
                      pt = [sum(pxc)/nPc, sum(pyc)/nPc];
                  end
                  %obj.ineqs(i)
                  drx1 = obj.ineqs(i).dfdx(vars(1));
                  drx2 = obj.ineqs(i).dfdx(vars(2));
                  %drx1
                  %drx2
                  % check sign
                  m0 = drx1.f/drx2.f;
                  n = n + 1;
                  %disp("Slope of tangent")
                  %pt
                  % l = true;
                  if abs(subs(drx2.f,vars,pt)) < 1.0d-6
                      m(n) = intmax;
                  else
                      m(n) = subs(m0,vars,pt);
                  end
                  

              end
              % if l
              %      disp("Slope of tangent")
              %      pt
              %      m
              % end

           end
           
         end

         % max over a region
         function [l, fmax, index, lsing] = maxArray (obj, f1, f2)
        % objective: which of f1, f2 is the larger ON obj -- l true with fmax/index when one of
        %   them dominates, l false when neither does (the caller then splits on f1 = f2).
        %
        % THE TIE CASE, AND WHY IT NO LONGER SAMPLES. When f1 and f2 agree at every vertex, the
        % vertices say nothing about which is larger inside, and this routine used to fall back
        % to PROBING: step 0.1 off a vertex along bisectors of the constraint slopes, and take
        % the first feasible probe's verdict as the answer. A probe is evidence, never proof,
        % and that shortcut returned wrong numbers on two independent inputs:
        %
        %   * f = x*y on [0,1]^2 as TWO triangles. The triangles' conjugates overlap in the lens
        %     {(s1+s2)^2 <= 4*s1, (s1+s2)^2 <= 4*s2}, whose only vertices (0,0) and (1,1) both
        %     lie on s1 = s2. One probe put s2 on the whole lens, so f*(0.66,0.18) came out 0.18
        %     for a truth of 0.66 -- wrong at 800 of 40000 grid points, and silent.
        %   * The 4-cone fan of conjCPLQTest/step3UnboundedAssemblyAgreesWithItsOwnPieces. On the
        %     quadrant {s1<=0, s2>=0} the candidates are s2^2/2 and s1^2/2 + s2^2/4, whose
        %     difference s2^2/4 - s1^2/2 changes sign inside; the cone's single vertex is the
        %     origin, where both vanish, and the only feasible probe (-0.1, 0) reported the
        %     wrong one, giving f*(-0.5,2) = 1.125 for a truth of 2.
        %
        % What replaces it is a SOUND test that is also cheaper: ask the symbolic engine whether
        % f1 - f2 has one sign EVERYWHERE. That is sufficient, not necessary -- a difference
        % that is sign-definite only on obj is reported as undecided -- and being undecided is
        % safe, because splitting on f1 = f2 is always VALUE-correct: the two functions are
        % equal on the split boundary, so even a half that survives as a degenerate sliver
        % carries the right value, and maximumP/maxEqDom already drop a half that comes out
        % empty. Splitting unnecessarily costs a cell; deciding wrongly costs a wrong answer.
        %
        % Refusing to guess also removed the whole probe apparatus -- slopes2, probeAlong,
        % probePerp, maxFromPt -- and with it the slopes2 call this routine made on EVERY
        % comparison, which was pure overhead outside the tie case.
        %
        % TWO LIMITS, both deliberate and both measured rather than assumed:
        %
        %   * The refusal applies only to a POLYNOMIAL pair. A region's constraints must be
        %     polynomial, so splitting on a rational f1 = f2 cannot be represented at all -- see
        %     the comment at that branch. Refusing there produced symbolic:coeffs:NotAPolynomial
        %     out of region.normalize1 on the biconjugate path, where every operand is rational.
        %   * The non-tie path below still concludes from vertex values alone, which is exact
        %     only when obj is the convex hull of its vertices -- a curved-edge or unbounded
        %     region is not. Same class of hole, not observed to fire wrongly, and closing it
        %     would split far more aggressively.
        %
        % Both limits are the same trade: this routine now refuses to guess exactly where the
        % refusal is actionable, and no further. An earlier, unrestricted version of this test
        % broke conjCPLQTest/multiFaceUnboundedConvexFacesConjugateExactly, which is how the
        % first limit came to be measured rather than guessed at.
          lsing = false;
          fv1 = obj.funcVertices (f1);
          fv2 = obj.funcVertices (f2);
          for i = 1:size(fv1,2)
              sv1(i) = fv1(i).f;
              sv2(i) = fv2(i).f;
          end

          l = true;
          if all(abs(double(sv1 - sv2))< 1.0d-14)
              g = simplifyFraction(f1.f - f2.f);
              % The variables MUST be re-declared REAL before asking. A bare sym('s_1') is
              % COMPLEX to the symbolic engine, and over the complex numbers -s_1^2/2 <= 0 is
              % simply not true -- so isAlways cannot decide the sign of any quadratic
              % difference, which is precisely the case that matters here, and every comparison
              % would fall through to an unnecessary split. (Observed: the 4-cone fan then
              % produced the vacuous constraint -s_1^2 <= 0 as a split boundary and crashed
              % downstream in removeTangent.) Substituting fresh real symbols asks the question
              % that was actually meant without touching global assumptions.
              gv = symvar(g);
              if ~isempty(gv)
                  rv = sym('maxArrayReal_%d', [1 numel(gv)], 'real');
                  g = subs(g, gv, rv);
              end
              if isAlways(g >= 0, 'Unknown', 'false')
                  fmax = f1;
                  index = 1;
                  return
              elseif isAlways(g <= 0, 'Unknown', 'false')
                  fmax = f2;
                  index = 2;
                  return
              end
              % REFUSE ONLY WHERE THE REFUSAL CAN BE ACTED ON. Saying "neither dominates" makes
              % the caller split on f1 = f2, and a region's constraints must be POLYNOMIAL:
              % splitmax3 hands f1 - f2 straight to region(), whose normalize1 raises
              % symbolic:coeffs:NotAPolynomial on a rational one. Every second-pass conjugate is
              % rational, so on the biconjugate path that error is the common case, not a corner
              % -- which is why this branch checks first.
              %
              % KNOWN HOLE, deliberately left: for a RATIONAL pair whose sign this test cannot
              % settle, the vertex comparison below still decides, and on an all-vertices-tied
              % cell that means "f2 wins" for no better reason than the order of the operands.
              % That is the same unsound shape the probing had, now confined to the one case
              % where the sound answer cannot be represented. Closing it needs splitmax3 to clear
              % denominators (sound only where both are provably nonzero on the cell) -- see
              % biconjugateTest's failure-site comment.
              if f1.isPolynomial && f2.isPolynomial
                  l = false;
                  fmax = 0;
                  index = 0;
                  return
              end
          end


          if all(double(sv1) <= double(sv2))
              fmax = f2;
              index=2;
          elseif all(double(sv2) <= double(sv1))
              fmax = f1;
              index=1;
          else
              l = false;
              fmax = 0;
              index=0;
          end
         
        end
        
        function [r] = splitmax2 (obj, f1, f2) 
          fv1 = obj.funcVertices (f1);
          fv2 = obj.funcVertices (f2);
          f = f1-f2;
          vars = f.getVars;
          x = sym('x');
          y = sym('y');
          f = subs(f.f, vars ,[x,y]);
          
          fx = [];
          fy = [];
          n = 0;
          fxy = [];
          for i = 1:size(obj.ineqs,2)
            g = subs(obj.ineqs(i).f, vars,[x,y]);
            s = solve ([f==0,g==0],[x,y]);
            if isempty(s)
              continue;
            end
            
            if isempty(s.x)
              continue;
            end
            sx = double(s.x);
            sy = double(s.y);
            if ~obj.ptFeasible (vars,[sx,sy])
                continue
            end
            fx = [fx,sx];
            fy = [fy,sy];
            n = n+1;
            fxy(n,1)= sx;
            fxy(n,2)= sy;
          end
          fxy
          fxy = unique(fxy,"rows");
          fx = fxy(:,1);
          fy = fxy(:,2);
          if size(fx,1) == 2
              m = (fy(2)-fy(1))/(fx(2)-fx(1));
              c = fy(1) - m*fx(1);
              ineq = vars(2)-m*vars(1)-c;
          end
          for i = 1:size(fv1,2)
              if (double(fv1(i).f) > double(fv2(i).f))
                  if subs(ineq,vars,[obj.vx(i),obj.vy(i)]) <= 0
                    r = [ineq,-ineq];
                    return
                  else
                    r = [-ineq,ineq];
                    return  
                  end
              end

          end
          
        end

        function [r] = splitmax3 (obj, f1, f2) 

          fv1 = obj.funcVertices (f1);
          fv2 = obj.funcVertices (f2);
          %disp('splitmax3.1')
          %fv1.printL
          %disp('splitmax3.2')
          %fv2.printL
          f = f1-f2;
           vars = f.getVars;
          ineq = f;
          for i = 1:size(fv1,2)
              if isAlways(fv1(i).f >= fv2(i).f)   % correct but slow
              %if (double(fv1(i).f) >= double(fv2(i).f))
              %if (abs(double(fv1(i).f)) >= abs(double(fv2(i).f)))
                  %ineq.f
                  %subs(ineq.f,vars,[obj.vx(i),obj.vy(i)])
                  %
                  % ineq = f1 - f2 is RATIONAL whenever the operands are (every second-pass
                  % conjugate is), and its denominator can vanish at precisely this vertex --
                  % `subs` then raises symbolic:kernel:DivisionByZero rather than returning a
                  % value. That is not a reason to give up on the whole split: it says only that
                  % THIS vertex cannot decide the sign, so move to the next one, and if none can,
                  % fall through to the well-defined default below. The same singular-at-a-vertex
                  % situation is handled by directional limits in vertexOfEdge and
                  % simplifyUnboundedRegion; here there is nothing to take a limit OF, since the
                  % vertex is only being used to orient an inequality.
                  %
                  % Reachable since maxArray stopped guessing: it now reports "neither dominates"
                  % far more often, so splitmax3 runs on operand pairs it never used to see --
                  % including rational ones sharing a pole with the cell's own vertex.
                  try
                      neg = isAlways(subs(ineq.f,vars,[obj.vx(i),obj.vy(i)]) < 0);
                  catch
                      continue
                  end
                  if neg
                    r = [ineq,-ineq];
                    return
                  else
                    r = [-ineq,ineq];
                    return
                  end
              end

          end

          % HISTORY: falling off the end here left `r` UNASSIGNED, so the caller got
          % MATLAB:unassignedOutputs ("Output argument r ... not assigned") from
          % functionNDomain.maximumP rather than a split -- the same "output only ever assigned
          % inside the loop" shape as region.getEdges had. The case is not degenerate: the loop
          % exits without returning exactly when NO vertex has f1 >= f2, i.e. f2 is the maximum
          % on the whole region and there is nothing to split.
          %
          % That still has a well-defined answer in this routine's own convention, which the
          % caller relies on: ineqs(1) delimits where f(1) wins and ineqs(2) where f(2) does
          % (maximumP assigns objL(i).f(1) to the first half and f(2) to the second). With
          % f1 < f2 throughout, the first half must come out EMPTY and the second must be the
          % whole region, which is exactly [-ineq, ineq] for ineq = f1 - f2: the first gives
          % {f2 - f1 <= 0} = empty, the second {f1 - f2 <= 0} = everything. maximumP already
          % drops an empty half (see its own isempty(d1) guard), so this needs nothing there.
          r = [-ineq, ineq];
        end

        
        function [nl, v1, v2] = splitmaxArray (obj, f1, f2) 
          fv1 = obj.funcVertices (f1);
          fv2 = obj.funcVertices (f2);
          v1 = [];
          v2 = [];
          for i = 1:size(fv1,2)
              sv1(i) = fv1(i).f;
              sv2(i) = fv2(i).f;
              % replace with limit
              if (isnan(sv1(i)))
                  sv1(i) = 0;
                  sv2(i) = 0;
              elseif (isnan(sv2(i)))
                  sv1(i) = 0;
                  sv2(i) = 0;
              end
              if sv1(i) == sv2(i)
                  v2 = [v2,i];
                  v1 = [v1,i];
              elseif sv1(i) > sv2(i)
                  v1 = [v1,i];
              else    
                  v2 = [v2,i];
              end

          end
          
          if size(v1,2) == 3
              nl(1) = 1;
          else
              nl(1) = 0;
          end
          if size(v2) == 3
              nl(2) = 1;
          else
              nl(2) = 0;
          end
          
        end
        
        function [lelim, obj] = deleteIneq (obj, vars)
          lelim = false;

          for i = 1:size(obj.ineqs,2)
              l(i) = obj.ineqs(i).f;
              mark(i) = false;
              nPts = 0;

              for j = 1:obj.nv
                  if obj.ineqs(i).subsF(obj.vars,[obj.vx(j),obj.vy(j)]).isZero
                      nPts = nPts+1;
                  end
              end
              if nPts <= 1
                  mark(i)=true;
              end
          end
          for i = 1:size(obj.ineqs,2)
              if ~mark(i)
                  continue
              end
              if obj.ineqs(i).isQuad
                  continue;
              end
             l1 = l;
             l1(i) = [];
             r = region (l1,vars);
             if isempty(r)
                 return
             end
             if obj.eqVertices(r)
                 lelim = true;
                 obj = r;
                 return;
             end
          end


        end

        function obj = simplifyOpenRegion1 (obj, nP, px, py)
            % remove ineqs that dont go through a vertex
            n = 0;
            mark = [];
            for i = 1:size(obj.ineqs,2)
                l = false;
                for j = 1:nP
                    val = obj.ineqs(i).subsF(obj.vars,[px(j),py(j)]).f;
                    % HISTORY: for a genuinely non-polynomial piece (e.g. a
                    % sqrt/fractional-power term, as in testFractional), the
                    % substituted expression can retain symbolic structure
                    % that doesn't auto-cancel (e.g. sqrt(3) terms across
                    % different constraints) -- double() then errors
                    % ("Unable to convert expression containing symbolic
                    % variables") instead of evaluating. simplify() first
                    % resolves the cases that are genuinely just an
                    % unreduced constant; if it's still not a evaluable
                    % constant after that, this ineq truly isn't (near) zero
                    % at this vertex, so treat it as such rather than crash.
                    try
                        v = abs(double(simplify(val)));
                    catch
                        v = Inf;
                    end
                    if (v < 1.0d-8)
                        l = true;
                        break
                    end
                end
                if l
                    continue;
                end
                n = n + 1;
                mark(n) = i;
            end
            obj = obj.deleteIfRedundant(mark);

        end

        % ---- the instance side of the LP certificates (see the static block at the top) ----

        function [A, b, lin] = linearForm (obj)
        % Numeric matrix form of obj's LINEAR constraints: constraint j, stored as
        % ineqs(j) <= 0, becomes A(j,:)*[vars(1);vars(2)] <= b(j). lin(j) records whether
        % constraint j is affine at all -- a quadratic or rational facet has no such row, its
        % row is left NaN, and every caller must consult lin before using it.
        %
        % Coefficients are read by EVALUATION, not by coeffs/getLinearCoeffs: an affine g
        % satisfies g = g(0,0) + (g(1,0)-g(0,0)) x + (g(0,1)-g(0,0)) y identically, so three
        % substitutions recover it exactly regardless of how the expression happens to be
        % written. That independence is the point -- it is the same lesson the note at the top of
        % this file records: test the geometry, never the syntax.
        %
        % PERFORMANCE: the six probe points are substituted into the WHOLE constraint vector at
        % once, six round trips to the symbolic engine per region instead of six PER CONSTRAINT.
        % `subs` is elementwise, so the numbers are identical; only the call count changes.
            n = size(obj.ineqs,2);
            A = nan(n,2); b = nan(n,1); lin = false(1,n);
            if n == 0
                return
            end
            pts = [0 0; 1 0; 0 1; 1 1; 2 1; 1 2];
            G = [obj.ineqs.f];
            V = nan(6, n);
            for k = 1:6
                try
                    V(k,:) = double(subs(G, obj.vars, pts(k,:)));
                catch
                    % One unevaluable constraint must not cost the others their rows, so fall
                    % back to per-constraint evaluation for this probe point only.
                    for j = 1:n
                        try
                            V(k,j) = double(subs(obj.ineqs(j).f, obj.vars, pts(k,:)));
                        catch
                            V(k,j) = NaN;
                        end
                    end
                end
            end
            for j = 1:n
                g = obj.ineqs(j);
                if ~g.isPolynomial
                    continue                       % a rational facet is not affine
                end
                c0 = V(1,j);
                c1 = V(2,j) - c0;
                c2 = V(3,j) - c0;
                if any(isnan(V(:,j)))
                    continue                       % not numerically evaluable: no row
                end
                % Confirm affineness rather than trusting the degree bookkeeping. For a
                % quadratic part alpha x^2 + beta xy + gamma y^2 these three points force
                % beta=0, then alpha=0, then gamma=0 in turn, so passing all three means
                % the quadratic part is identically zero.
                chk = [1 1; 2 1; 1 2];
                aff = true;
                for k = 1:3
                    want = c0 + c1*chk(k,1) + c2*chk(k,2);
                    got  = V(3+k, j);
                    if abs(got - want) > 1.0d-9 * max(1, abs(want))
                        aff = false; break
                    end
                end
                if ~aff || ~all(isfinite([c0 c1 c2]))
                    continue
                end
                A(j,:) = [c1, c2];
                b(j)   = -c0;
                lin(j) = true;
            end
        end

        function [tf, why] = quadUnboundedBelow (obj, Q, L)
        % Is the quadratic q(x) = 1/2 x'Qx + L'x + c unbounded BELOW on this region?
        %
        % This is the gate Step 1 needs before it can take a convex envelope over an UNBOUNDED
        % face: conv q is -inf exactly when q is unbounded below, and only otherwise is there a
        % finite envelope to compute (the convex case of which is [GARDINER-11]/[GARDINER-13]).
        % The constant term is irrelevant and is not taken.
        %
        % THE TEST. Along a ray x0 + t*d the quadratic is
        %     q(x0) + t*(Q*x0 + L)'d + (t^2/2)*d'Q*d,
        % so q is unbounded below on a region P with recession cone R exactly when either
        %   (i)  some d in R has d'Q*d < 0 -- the quadratic itself decays along a direction the
        %        region never leaves; or
        %   (ii) some d in R has d'Q*d == 0 while the linear slope (Q*x+L)'d can be made
        %        negative somewhere in P, which is an LP (region.maxLinear).
        % Everything else is bounded below. Equivalently, and this is the usual statement:
        % diagonalize Q and compare its NEGATIVE-eigenvalue directions against R -- if the cone
        % {d : d'Q*d < 0} meets R, the value is -inf.
        %
        % WHY IT IS CLOSED FORM IN 2D, with no cone enumeration. Writing d = (cos t, sin t),
        %     d'Q*d = (Q11+Q22)/2 + ((Q11-Q22)/2)*cos(2t) + Q12*sin(2t),
        % a pure sinusoid in 2t; and each affine constraint a'x <= b contributes a'd <= 0 to R,
        % i.e. a half-circle of admissible t. R is therefore an arc, so the minimum of d'Q*d
        % over it is attained either at a constraint boundary (a'd == 0) or at the sinusoid's own
        % critical angle when that lies inside R. Both are finite, explicit candidate sets, so no
        % search and no cone-generator enumeration is needed. A bounded region contributes no
        % admissible direction at all and is reported bounded below, as it must be.
        %
        % NON-AFFINE FACETS ARE REFUSED, not dropped. Dropping them would ENLARGE the region and
        % so enlarge R, which is unsound in the one direction that matters here: it could report
        % -inf for a direction the true region never actually recedes along. Step 1's faces are
        % polyhedral, so this does not bite; anything else must be handled explicitly.
            tf = false; why = 'bounded below';
            if isempty(obj)
                return
            end
            [A, b, lin] = obj.linearForm;
            if ~all(lin)
                error('region:quadUnboundedBelow:nonAffineFacet', ...
                    ['quadUnboundedBelow needs every facet affine: a curved facet''s recession ' ...
                     'behaviour is not a half-plane, and dropping it would enlarge the recession ' ...
                     'cone and could wrongly certify -inf.']);
            end
            Q = double(Q); L = double(L(:));
            scaleA = max(1, max(abs(A(:))));
            tolA = 1e-9 * scaleA;

            % Candidate directions: every constraint boundary, plus the sinusoid's own critical
            % angles (4 of them mod 2*pi, from the two roots of tan(2t) = Q12/((Q11-Q22)/2)).
            cand = [];
            for j = 1:size(A,1)
                a = A(j,:);
                if norm(a) <= tolA, continue, end
                p = [-a(2), a(1)] / norm(a);       % both senses along the boundary line a'd = 0
                cand = [cand; p; -p]; %#ok<AGROW>
            end
            c2 = (Q(1,1) - Q(2,2))/2; s2 = Q(1,2);
            if abs(c2) + abs(s2) > 0
                t0 = 0.5 * atan2(s2, c2);
                for k = 0:3
                    tk = t0 + k*pi/2;
                    cand = [cand; cos(tk), sin(tk)]; %#ok<AGROW>
                end
            elseif isempty(cand)
                cand = [1 0; 0 1; -1 0; 0 -1];      % Q is a multiple of I: any direction will do
            end
            if isempty(cand)
                return                              % no constraints and no curvature: q is affine
            end

            % Keep only the directions the region actually recedes along.
            keep = true(size(cand,1),1);
            for i = 1:size(cand,1)
                if any(A * cand(i,:)' > tolA), keep(i) = false; end
            end
            cand = cand(keep,:);
            if isempty(cand)
                return                              % bounded region: nothing to recede along
            end

            vals = zeros(size(cand,1),1);
            for i = 1:size(cand,1)
                vals(i) = cand(i,:) * Q * cand(i,:)';
            end
            tolQ = 1e-9 * max(1, max(abs(Q(:))));
            if min(vals) < -tolQ
                tf = true;
                why = 'a recession direction has d''Qd < 0';
                return
            end
            % (ii) d'Qd == 0 along some recession direction: the ray is a straight line for q, so
            % it is the LINEAR slope that decides, and that slope varies over the region.
            for i = 1:size(cand,1)
                if vals(i) > tolQ, continue, end
                d = cand(i,:)';
                % min over P of (Q*x + L)'d  =  -max over P of (-Q*d)'x  -  ... (constant L'd)
                [mx, st] = region.maxLinear(A, b, (-Q*d)');
                if st == 1
                    tf = true; why = 'a d''Qd == 0 recession direction has unbounded negative slope';
                    return
                end
                if st ~= 0
                    error('region:quadUnboundedBelow:undecided', ...
                        ['the LP deciding the linear slope along a d''''Qd == 0 recession ' ...
                         'direction returned no certificate (status %d); refusing to guess.'], st);
                end
                if -mx + L'*d < -1e-9 * max(1, abs(L'*d))
                    tf = true; why = 'a d''Qd == 0 recession direction has a negative slope';
                    return
                end
            end
        end

        function [D, kind] = recessionRays (obj)
        % The EXTREME RAYS of this region's recession cone, as unit-length direction rows of D,
        % read from the region's INEQUALITIES -- never from obj.vx/obj.vy.
        %
        % That last point is the whole reason this exists. An unbounded region stores its
        % directions as vertices flagged +/-intmax (see finiteVertices, and getVertices' own
        % "putting intmax for inf to avoid Nans" comment). Those markers are fine for the
        % ordering and feasibility work region does with them, but they are NOT coordinates,
        % and Step 1 reads coordinates: plq_1p's envelope formulas are built from vertex
        % values, so a face carrying an intmax vertex produces an envelope with 2147483647 and
        % intmax^2 = 4611686014132420609 in it (measured -- see plq_1p's own header). The
        % recession directions have to come from the half-planes instead, and this returns them.
        %
        % kind reports the SHAPE of the cone R = {d : A*d <= 0}, since callers must branch on it:
        %   'bounded'    R = {0}: the region is bounded, D is empty.
        %   'ray'        R is a single ray (one extreme ray).
        %   'wedge'      R is a pointed 2-dimensional cone (two extreme rays).
        %   'nonpointed' R contains a line -- a half-plane, a full line, or all of R^2. D then
        %                holds whatever boundary directions were found, but the region has no
        %                apex to fan from and callers should refuse it rather than guess.
        %
        % WHY IT IS EXACT AND NEEDS NO CONE ENUMERATION. Writing d = (cos t, sin t), each
        % constraint a'x <= b contributes a'd <= 0, i.e. the half-circle of angles t at least
        % pi/2 away from a's own direction. R is the intersection of those half-circles, hence
        % an ARC, and an arc's endpoints are boundary points of some contributing half-circle.
        % So every extreme ray is one of the 2m angles phi_j +/- pi/2, and testing that finite
        % candidate set against A is exhaustive. This is the same closed-form observation
        % quadUnboundedBelow is built on, and the two are meant to be read together.
        %
        % Non-affine facets are REFUSED for quadUnboundedBelow's reason: dropping one enlarges
        % the region, and so enlarges the recession cone, which is the unsafe direction here.
            D = zeros(0,2); kind = 'bounded';
            if isempty(obj)
                return
            end
            [A, ~, lin] = obj.linearForm;
            if ~all(lin)
                error('region:recessionRays:nonAffineFacet', ...
                    ['recessionRays needs every facet affine: a curved facet''s recession ' ...
                     'behaviour is not a half-plane, and dropping it would enlarge the ' ...
                     'recession cone.']);
            end
            tolA = 1e-9 * max(1, max(abs(A(:))));
            rows = zeros(0,2);
            for j = 1:size(A,1)
                if norm(A(j,:)) > tolA
                    rows(end+1,:) = A(j,:) / norm(A(j,:)); %#ok<AGROW>
                end
            end
            if isempty(rows)
                kind = 'nonpointed';           % no constraints at all: R is the whole plane
                return
            end
            A = rows;

            % Candidate extreme directions: the boundary of each constraint's half-circle. The
            % VECTOR [-a2,a1] is kept alongside its angle and is what gets returned -- a
            % direction of a rational half-plane is itself rational, and rebuilding it as
            % (cos t, sin t) would not be. That is not cosmetic: the returned direction goes on
            % to BUILD the sub-face inequalities in fanUnboundedFace, where a 6.1e-17 in place
            % of a 0 turned `x <= 0` into `x - 4967757600021511/8.1e31*y <= 0` -- a genuinely
            % different, and no longer pointed, half-plane.
            candV = zeros(0,2); candTh = [];
            for j = 1:size(A,1)
                p = [-A(j,2), A(j,1)];
                candV = [candV; p; -p];                                        %#ok<AGROW>
                candTh = [candTh; atan2(p(2), p(1)); atan2(-p(2), -p(1))];     %#ok<AGROW>
            end
            tol  = 1e-9;
            isFeas = @(th) all(A * [cos(th); sin(th)] <= tol);

            keepF = arrayfun(isFeas, candTh);
            if ~any(keepF)
                kind = 'bounded';              % no admissible direction: nothing to recede along
                return
            end
            candV = candV(keepF,:); candTh = candTh(keepF);

            % A direction is an EXTREME ray when the arc stops there, i.e. at least one of its
            % two immediate neighbours leaves R. delta is well below the spacing of the
            % candidate set, which is what the classification actually has to resolve.
            delta = 1e-7;
            D = zeros(0,2);
            for i = 1:numel(candTh)
                th = candTh(i);
                if isFeas(th + delta) && isFeas(th - delta)
                    continue                   % interior of the arc, not an extreme ray
                end
                d = candV(i,:) / max(abs(candV(i,:)));   % rational, and O(1) in size
                dup = false;
                for k = 1:size(D,1)
                    if norm(D(k,:) - d) < 1e-9, dup = true; break, end
                end
                if ~dup, D(end+1,:) = d; end %#ok<AGROW>
            end

            % R contains a line exactly when A has a nontrivial null space, so rank decides
            % pointedness outright -- no need to test antipodal pairs.
            if rank(A, tolA) < 2
                kind = 'nonpointed';
                return
            end
            switch size(D,1)
                case 0,    kind = 'bounded';
                case 1,    kind = 'ray';
                otherwise, kind = 'wedge';
            end
        end

        function del = redundantSubset (obj, cand)
        % Of the candidate indices cand (into obj.ineqs), the subset that is PROVABLY
        % redundant and may therefore be deleted without changing the feasible set:
        %
        %     ineq i is redundant  <=>  max{ g_i(x) : g_j(x) <= 0 for all j ~= i }  <=  0.
        %
        % This replaces the old proxy "delete any constraint that does not pass through a
        % finite vertex", which is not a redundancy test at all for an unbounded or curved
        % region and was the dominant half of the wrong-partition defect: a region carrying
        % (s1+s2)^2/4 lost -s1-s2 <= 0 and then reported f*(-3,-3) = 9 where the truth is 0.
        %
        % Two deliberate asymmetries, both erring toward keeping a constraint:
        %   * a non-affine g_i is never reported redundant (no LP certificate for it);
        %   * non-affine g_j are DROPPED from the constraint set the test runs against. That
        %     relaxation only ENLARGES the feasible set, so redundancy in the relaxation
        %     implies redundancy here -- sound in the one direction that matters.
        % Candidates are processed one at a time against the constraints still standing, so
        % two copies of the same constraint cannot certify each other and both be deleted.
            del = [];
            if isempty(obj) || isempty(cand)
                return
            end
            [A, b, lin] = obj.linearForm;
            standing = true(1, size(obj.ineqs,2));
            for t = 1:numel(cand)
                i = cand(t);
                if i < 1 || i > numel(standing) || ~lin(i) || ~standing(i)
                    continue
                end
                use = lin & standing;
                use(i) = false;
                [val, st] = region.maxLinear(A(use,:), b(use), A(i,:));
                if st == 0 && val <= b(i) + 1.0d-9 * max(1, abs(b(i)))
                    del(end+1) = i;             %#ok<AGROW>
                    standing(i) = false;
                end
            end
        end

        function obj = deleteIfRedundant (obj, cand)
        % Delete exactly those of the candidate constraints cand that redundantSubset can
        % certify as redundant, and keep the rest. Every deletion site in this class goes
        % through here, so a heuristic upstream can only ever PROPOSE a deletion.
            del = obj.redundantSubset(cand);
            if ~isempty(del)
                obj.ineqs(del) = [];
            end
        end

        function [tf, why] = certifiesNonPositive (objP, h)
        % A SOUND yes/unknown answer to "does h <= 0 hold everywhere on objP?", for an h that is
        % NOT affine. `false` means "not certified" and never "false" -- every uncertain answer
        % is a refusal, exactly as the LP certificates at the top of this file are.
        %
        % WHY IT EXISTS. unionIsExact needs "every constraint of B' holds on A". For an affine
        % constraint that is an LP. For a curved one there was no certificate at all, so it
        % refused outright -- and `merge` grew two HEURISTICS to stand in for the missing
        % answer (an all-quadratics-must-match cross product, and a "does this conic cut the
        % other region away from a vertex" probe). MEASURED 2026-08-17 on the all-rational
        % control case in `.claude/step3cost.m`: those two refused 322 of 612 merge attempts at
        % one fold, against 4 that succeeded, and Step 3's cell count ran to 57 for 10 distinct
        % functions. The heuristics were the cost.
        %
        % THE ARGUMENT, which is why the test below is exactly these two checks.
        % Take P = the LINEAR RELAXATION of objP -- its affine facets only. P contains objP, so
        % `h <= 0 on P` implies `h <= 0 on objP`; dropping the curved facets can only make the
        % certificate harder to obtain, never wrong. That is the same relaxation unionIsExact
        % already applies to the region it tests OVER.
        %
        % Let h be a CONVEX quadratic, h(x) = 1/2 x'Qx + L'x + c with Q PSD, and let P have at
        % least one vertex, so rec(P) = {d : Ad <= 0} is pointed and generated by at most two
        % extreme rays, each lying on the boundary of one of P's own constraints. Write any
        % x in P as x0 + d with x0 in conv(vertices) and d in rec(P). Then
        %   * h <= max_i h(v_i) on conv(vertices), by convexity;
        %   * for an extreme ray e we require Qe = 0, so Q(sum beta_k e_k) = 0 and
        %     grad h(x0).d = x0'Qd + L'd = L'd, INDEPENDENT of x0;
        %   * hence h(x0 + d) = h(x0) + L'd, and requiring L'e <= 0 on each extreme ray gives
        %     h(x) <= max_i h(v_i).
        % So `h(v) <= 0 at every vertex` together with `Qe = 0 and L'e <= 0 on every extreme
        % ray` is SUFFICIENT. Both are finite checks in closed form.
        %
        % Everything outside that hypothesis is refused rather than approximated: a rational h,
        % a cubic or higher h, an indefinite or concave Q (whose max over P sits nowhere the
        % vertices can see), and a P with no vertex at all (a half-plane, a slab, the plane --
        % where this argument has no base point to stand on).
        %
        % The arithmetic here is NUMERIC on purpose. This is a DECISION, not a value that flows
        % into the geometry, and it is the same standing this file's LP certificates have; the
        % tolerances match region.impliedBy's.
            persistent qopts
            tf = false;
            why = 'notPolynomial';
            if isempty(objP) || ~h.isPolynomial
                return
            end
            vars = objP.vars;
            % Q, L, c by differentiation -- exact, and independent of how h is written.
            try
                Qs = hessian(h.f, vars);
                why = 'degreeAboveTwo';
                if any(has(Qs(:), vars))
                    return                      % degree 3 or higher: not a quadratic
                end
                why = 'unreadableParts';
                Q = double(Qs);
                L = double(subs(gradient(h.f, vars), vars, [0 0]));
                L = L(:);
                c = double(subs(h.f, vars, [0 0]));
            catch
                return
            end
            if any(~isfinite([Q(:); L(:); c]))
                return
            end
            Q = (Q + Q')/2;
            scaleQ = max(1, max(abs(Q(:))));
            evQ = eig(Q);
            isConvexH = min(evQ) >= -1.0d-9 * scaleQ;
            isConcaveH = max(evQ) <=  1.0d-9 * scaleQ;
            why = 'indefinite';
            if ~isConvexH && ~isConcaveH
                return          % neither argument below applies to a saddle
            end

            [A, b, lin] = objP.linearForm;
            A = A(lin,:); b = b(lin);           % the LINEAR RELAXATION: a superset of objP
            why = 'noConstraints';
            if isempty(b)
                return                          % the whole plane: no vertex to stand on
            end
            scaleA = max(1, max(abs([A(:); b(:)])));
            tolA = 1.0d-9 * scaleA;
            hAt = @(z) 0.5*(z(:).'*Q*z(:)) + L.'*z(:) + c;
            scaleH = max(1, max(abs([Q(:); L(:); c])));
            tolH = 1.0d-9 * scaleH;

            % VERTICES of the relaxation: every pair of constraints that meets in one point,
            % kept when that point is feasible for all of them.
            nA = size(A,1);
            V = zeros(0,2);
            for i = 1:nA
                for j = i+1:nA
                    M = [A(i,:); A(j,:)];
                    if abs(det(M)) <= tolA^2 + eps
                        continue                % parallel: no vertex here
                    end
                    z = M \ [b(i); b(j)];
                    if all(A*z <= b + tolA)
                        V(end+1,:) = z(:).';    %#ok<AGROW>
                    end
                end
            end
            % THE CONCAVE CASE, and it is the one that actually fires here. A parabolic facet
            % like -(s1+2*s2)^2 + 16*s1 has a NEGATIVE semidefinite Hessian, so the vertices
            % bound nothing -- a concave function is LARGEST in the middle. Measured
            % 2026-08-17 on the A.4/A.5 quadrilateral: every refusal this routine made came
            % back `curved_notConvex`, all of them this shape.
            %
            % Maximising a concave quadratic over a polyhedron is a CONVEX program -- minimise
            % -h, whose Hessian -Q is PSD -- so quadprog decides it outright, in the same
            % standing as the LP certificates this file already relies on (same toolbox, same
            % tolerance convention, and an unbounded or undecided answer is a refusal).
            if isConcaveH && ~isConvexH
                why = 'concaveQpUndecided';
                if isempty(qopts)
                    qopts = optimoptions('quadprog', 'Display', 'none');
                end
                ws = warning('off', 'all');
                try
                    [~, negMax, ef] = quadprog(-Q, -L, A, b, [], [], [], [], [], qopts);
                catch
                    ef = 0; negMax = [];
                end
                warning(ws);
                if ef == -2
                    tf = true; why = 'ok';      % P is empty: h <= 0 on it vacuously
                    return
                end
                if ef ~= 1
                    return                      % unbounded above, or no certificate
                end
                why = 'positiveSomewhereOnP';
                if -negMax + c > tolH
                    return                      % h really is positive somewhere on P
                end
                tf = true; why = 'ok';
                return
            end

            why = 'noVertex';
            if isempty(V)
                return                          % no vertex: outside the argument's hypothesis
            end
            for k = 1:size(V,1)
                why = 'positiveAtAVertex';
                if hAt(V(k,:)) > tolH
                    return                      % h is positive somewhere on P
                end
            end

            % EXTREME RAYS of rec(P) = {d : Ad <= 0}. In the plane each lies along the boundary
            % of one of P's own constraints, so both senses of each constraint direction are
            % the complete candidate list.
            for i = 1:nA
                nrm = norm(A(i,:));
                if nrm <= tolA, continue, end
                p = [-A(i,2), A(i,1)] / nrm;
                for sgn = [1 -1]
                    d = sgn * p(:);
                    if any(A*d > tolA)
                        continue                % not a recession direction
                    end
                    why = 'growsAlongARay';
                    if abs(d.'*Q*d) > 1.0d-9 * scaleQ
                        return                  % h grows quadratically along a ray of P
                    end
                    if L.'*d > tolH
                        return                  % h grows linearly along a ray of P
                    end
                end
            end
            tf = true;
            why = 'ok';
        end

        function [l, why] = unionIsExact (objA, objB, ia, ib)
        % Precondition: objA.ineqs(ia) == -objB.ineqs(ib), i.e. the two regions meet on the
        % facet {g = 0} of that shared constraint, with objA on {g <= 0} and objB on {g >= 0}.
        %
        % merge unions them by deleting that facet from both and intersecting what is left:
        % writing A = A' n {g<=0} = objA and B = B' n {g>=0} = objB, it returns M = A' n B'.
        %
        % M contains A u B in one direction ALWAYS: a point of A' n B' either has g <= 0, and
        % then lies in A' n {g<=0} = A, or has g >= 0, and then lies in B' n {g>=0} = B. So
        % merge can never LOSE a point. The other direction is the one that can fail, and it
        % needs exactly
        %       A subset B'   and   B subset A',
        % i.e. every constraint of B' holds everywhere on A and every constraint of A' holds
        % everywhere on B -- which is also precisely the statement that A u B is convex.
        %
        % Each of those is "is this linear form <= 0 over that polyhedron", an LP.
        %
        % Two roles, two different requirements, and only one of them is strict:
        %   * the constraints being TESTED (those of A' and B') must be affine -- there is no
        %     LP certificate for a curved one, so refuse outright if any is not;
        %   * the region being tested OVER enters only as its LINEAR RELAXATION, dropping any
        %     curved facet. That relaxation is a SUPERSET, so "valid on the relaxation"
        %     implies "valid on the region" -- sound, and it lets a pair meeting along a
        %     shared PARABOLIC facet still be decided, which is the shape this codebase
        %     actually produces (a parabolic arc between two rays).
        % Refusing only leaves the two regions separate, which is always correct, just less
        % compact -- so every uncertain answer is a refusal.
            l = false;
            [AA, bA, linA] = objA.linearForm;
            [AB, bB, linB] = objB.linearForm;
            keepA = true(1, numel(bA)); keepA(ia) = false;   % A' = A minus the shared facet
            keepB = true(1, numel(bB)); keepB(ib) = false;   % B' = B minus the shared facet

            % A subset B'? Affine constraints of B' by LP, curved ones by the closed-form
            % certificate. A curved constraint used to refuse outright here; that refusal is
            % what merge's two quadratic heuristics were standing in for.
            why = 'exactAnotInB';
            if ~region.impliedBy(AB(keepB & linB,:), bB(keepB & linB), AA(linA,:), bA(linA))
                return
            end
            for j = find(keepB & ~linB)
                [okC, whyC] = objA.certifiesNonPositive(objB.ineqs(j));
                if ~okC
                    why = ['curved_' whyC];
                    return
                end
            end
            why = 'exactBnotInA';
            if ~region.impliedBy(AA(keepA & linA,:), bA(keepA & linA), AB(linB,:), bB(linB))
                return
            end
            for i = find(keepA & ~linA)
                [okC, whyC] = objB.certifiesNonPositive(objA.ineqs(i));
                if ~okC
                    why = ['curved_' whyC];
                    return
                end
            end
            l = true;
            why = 'ok';
        end

        function [nP, px, py] = finiteVertices (obj)
        % The region's vertices that are genuine POINTS: an entry flagged intmax stands for a
        % direction of an unbounded region, not a vertex. This is the (nP,px,py) triple
        % removeTangent wants -- it uses them as candidate tangency points, so they must be
        % this region's own vertices, and must be read BEFORE simplifyUnboundedRegion drops
        % any of them.
            nP = 0;
            px = sym.empty();
            py = sym.empty();
            if isempty(obj)
                % An empty region still reaches here through mergeL's accumulation loops
                % (region.merge's own empty-operand guard leaves the accumulator empty);
                % obj.nv on a 0x0 region array is a comma-separated list with no values.
                return
            end
            for j = 1:obj.nv
                if abs(obj.vx(j)) == intmax || abs(obj.vy(j)) == intmax
                    continue
                end
                nP = nP + 1;
                px(nP) = obj.vx(j);
                py(nP) = obj.vy(j);
            end
        end

        function obj = removeTangent (obj, nP, px, py)
            n = 0;
            mark = [];
            vars = obj.vars;
            for i = 1:size(obj.ineqs,2)
              if ~ obj.ineqs(i).isQuad 
                continue
              end
              % get feasible point
              for j = 1:nP
                if ~obj.ineqs(i).subsF(vars,[px(j),py(j)]).isZero  
                  continue
                end
                % A TANGENT NEEDS A GRADIENT. At a vertex where this quadratic's gradient
                % VANISHES there is no tangent line -- every direction is tangent -- and whatever
                % `tangent` returns there is meaningless. Concluding from it deletes constraints
                % the region actually needs.
                %
                % That is not a corner case: it is the apex of a cone, which is exactly where the
                % Step 3 split conics of an unbounded fan meet. SUPPORT_MATRIX.md 8.2(e) records
                % the same trap in simplifyUnboundedRegion -- "the split conic's gradient vanishes
                % at exactly that vertex, so those directions are meaningless" -- fixed there by
                % region.witnessAwayFrom. This is its sibling.
                %
                % MEASURED, the 4-cone fan of conjCPLQTest/step3UnboundedAssemblyMatchesTheTruth:
                % the Step 3 cell {s2 <= 0, s2^2/2 - s1^2 <= 0, -s1 <= 0, s1^2 - 2*s2^2 <= 0} has
                % both quadratics' gradients vanishing at its own apex, the origin. removeTangent
                % deleted -s1 <= 0 there. The two quadratics are BLIND TO THE SIGN of s1 -- they
                % constrain |s1| -- so the region became symmetric under s1 -> -s1 and claimed the
                % mirror wedge: f*(-3,-2.4) came back 5.130 for a truth of 4.500.
                gI = [obj.ineqs(i).dfdx(vars(1)).subsF(vars,[px(j),py(j)]).f, ...
                      obj.ineqs(i).dfdx(vars(2)).subsF(vars,[px(j),py(j)]).f];
                if isAlways(gI(1) == 0) && isAlways(gI(2) == 0)
                    continue
                end
                tangent = obj.ineqs(i).tangent(px(j),py(j));
                %tangent
                tangent = tangent.normalize1;
                %disp('tan')
                %tangent
                % Get point in parabola.
                %
                % s0 is the solution set of a QUADRATIC in vars(2), so it can be EMPTY (the
                % parabola has no point at sx) or COMPLEX, and this used to take s0(1)
                % unconditionally: empty threw "Index exceeds array bounds", and a complex root
                % made the isAlways tests below undecidable, so `lin` and `tin` both came back
                % false and the routine went on to mark constraints from a nonsense probe point.
                % Same "root chosen before it is checked" shape as the py(1)-before-isempty
                % inversions corrected three times over in getNormalConeVertexQ. There is no
                % sensible answer from a point that does not exist: skip this vertex instead.
                sx = px(j) + 0.1;
                p = obj.ineqs(i).subsF([vars(1)],[sx]);
                s0 = solve(p.f,vars(2));
                sy = [];
                for r0 = 1:numel(s0)
                    try
                        ok0 = isreal(double(s0(r0)));
                    catch
                        ok0 = false;   % not numerically evaluable: not a usable probe point
                    end
                    if ok0
                        sy = s0(r0);
                        break
                    end
                end
                if isempty(sy)
                    continue
                end
                mx = (px(j) + sx)/2;
                my = (py(j) + sy)/2;
                lin = isAlways(obj.ineqs(i).subsF(vars,[mx,my]).f <= 0);
                for k = 1:size(obj.ineqs,2)
                  l1 = isAlways(obj.ineqs(k).f==tangent.f);
                  if isAlways(obj.ineqs(k).f==-tangent.f);
                      tangent = -tangent;
                      l1 = true;
                  end
                  %obj.ineqs(k).f
                  if ~l1
                      continue
                  end
                   [nvi, vxi, vyi] = obj.vertexOfEdge(k);
                  if nvi > 1
                      continue
                  end
                  tin = isAlways(tangent.subsF(vars,[mx,my]).f<= 0);
                  if lin & tin
                    n = n + 1;
                    mark(n) = k;
                  end  
                  if lin & ~tin
                    obj = region.empty;
                    return
                  end  
                  if ~lin & tin
                    % do nothing
                  end  
                  if ~lin & ~tin
                    n = n + 1;
                    mark(n) = i;
                  end  
                  %mark
                  break;
                end 
              end
            end
            %    disp('removeTangent')
                %obj.print
                %mark
            obj.ineqs(mark) = [];
            % check
            nv = obj.nv;
            obj = obj.getVertices;
   % obj.print
            %obj.print   
            if nv ~= obj.nv
                disp('Error in removeTangent')
                disp('This is okay when 3 vertices are on parabola - its removing middle one')
            end
        end
        
        function obj = simplifyOpenRegion (obj)
            nP = 0;
            for j = 1:obj.nv
              if abs(obj.vx(j)) == intmax
                continue
              end
              if abs(obj.vy(j)) == intmax
                continue
              end
              nP = nP+1;
              px(nP) = obj.vx(j);
              py(nP) = obj.vy(j);
            end
            obj = obj.simplifyOpenRegion1 (nP, px, py);
            for j = 1:nP
                nPoint(j) = 0;
            end
            keep = [];
            for i = 1:size(obj.ineqs,2)
                keep(i) = false;
            end
            for i = 1:size(obj.ineqs,2)
                for j = 1:nP
                    
                    if (abs(obj.ineqs(i).subsF(obj.vars,[px(j),py(j)]).f) > 1.0d-8)
                        continue;
                    end
                    nPoint(j)=nPoint(j)+1;
                    point(j,nPoint(j)) = i;
                    if obj.ineqs(i).isQuad
                        continue
                    end
                    
                end
            end
            if all(nPoint == 2)
                return;
            end
            for i = 1:size(obj.ineqs,2)
               if obj.ineqs(i).isQuad
                   return
               end
            end
            m0 = obj.slopes;
            n0 = 0;
            n1 = 0;
            markF0 = [];
            for i = 1:size(obj.ineqs,2)
                markF0(i) = i;
               if obj.ineqs(i).isQuad
                n0 = n0+1;
                % replace by slope of tangent and keep it
                m(n0) = intmax;
               else
                n0 = n0+1;
                n1 = n1+1;
                m(n0) = m0(n1);
               end
            end
           
            
            markF = [];
            for ip = 1:nP
                sx = px(ip);
                sy = py(ip);
                pi0 = point(ip,1:nPoint(ip));
                
                mark = [];
                for j = 1:nPoint(ip)
                    if obj.ineqs(pi0(j)).isQuad
                        mark = [mark,j];
                    end
                end
                
                pi0(mark) = [];
                mp = m(pi0);
                [sorted_m, indices] = sort(mp);
              for i = 1:size(indices,2)
                m1 = sorted_m(i);
                if i == size(indices,2)
                    m2 = sorted_m(1);
                else
                    m2 = sorted_m(i+1);
                end
                if (abs(m1)~= inf) & (abs(m2)~= inf)
                  d =  (m1+m2 )/2;
                else
                    % put temp on 14 th May
                    return
                  if (abs(m1)==inf)
                    d = tan((pi/2 + atan(m2))/2);
                  else
                      
                    d = tan((pi/2 + atan(m1))/2);
                  end
                end
                if i == size(indices,2)
                    d = -d;
                end
                c = sy - d * sx;
                tx = sx + 0.1;
                ty = d*tx+c   ;    
                if obj.ptFeasible (obj.vars,[tx,ty])
                   lm = false;
                   continue
                end
                tx = sx - 0.1;
                ty = d*tx+c   ;    
                if obj.ptFeasible (obj.vars,[tx,ty])
                  lm = false;
                  continue
                end


                if i == 1
                    m2 = sorted_m(size(indices,2));
                else
                    m2 = sorted_m(i-1);
                end
                if (abs(m1)~= inf) & (abs(m2)~= inf)
                  d =  (m1+m2 )/2;
                else
                  if (abs(m1)==inf)
                    d = tan((pi/2 + atan(m2))/2);
                  else
                      
                    d = tan((pi/2 + atan(m1))/2);
                  end
                end
                if i == 1
                    d = -d;
                end
                c = sy - d * sx;
                tx = sx + 0.1;
                ty = d*tx+c   ;    
                if obj.ptFeasible (obj.vars,[tx,ty])
                   lm = false;
                   continue
                end
                tx = sx - 0.1;
                ty = d*tx+c   ;    
                if obj.ptFeasible (obj.vars,[tx,ty])
                  lm = false;
                  continue
                end
                markF = [markF,pi0(indices(i))];
                
              end
              for i = 1:size(obj.ineqs,2)
                  keep0(i) = false;
              end
              for i = 1:size(pi0,2)
                  keep0(pi0(i)) = true;
              end
              for i = 1:size(markF,2)
                  keep0(markF(i)) = false;
              end
              keep = keep | keep0; 
            end
            markF0 = [];
            for i = 1:size(pi0,2)
                if keep(pi0(i))
                    continue
                end
                markF0 = [markF0,pi0(i)];
            end
            
            if size(markF0,2) == size(obj.ineqs,2)
                obj = region.empty;
                %disp("Singleton")
                return
                
            end
            obj.ineqs(markF0) = [];
            % check vertices


        end
            
        function tf = witnessAwayFrom (obj, px, py)
        % objective: exhibit a feasible point of obj that is NOT one of the listed vertices, so
        %   that "this region has collapsed to a point" can be REFUTED rather than assumed.
        %
        % [output] tf : true iff such a point was found. A false return means "no witness
        %               found", NOT "empty" -- this routine can only ever veto a claim of
        %               emptiness, never make one.
        %
        % WHY. simplifyUnboundedRegion decides a region has degenerated by stepping 0.1 off a
        % vertex along bisectors of the CONSTRAINT SLOPES there and testing feasibility; if no
        % such step is feasible it returns region.empty. A quadratic constraint whose gradient
        % VANISHES at that vertex has no slope there, and the probe directions built from it are
        % meaningless. That is not a corner case: cPLQ's Step 3 splits a cell on f1 = f2, and
        % when f1 and f2 are quadratics agreeing at the cell's apex -- e.g. s2^2/2 against
        % s1^2/2 + s2^2/4 on the quadrant {s1<=0, s2>=0} of conjCPLQTest's 4-cone fan -- the
        % split conic is singular at exactly that apex. The half where s2^2/2 wins,
        % {s1<=0, s2>=0, s1^2/2 - s2^2/4 <= 0}, is a genuine 2-D cone containing (-0.5,2),
        % (-0.1,3), (-1,4)... and was being declared empty, which is how f*(-0.5,2) went missing.
        %
        % Soundness: a genuinely empty region has no feasible point at all, so no probe here can
        % succeed and no true emptiness verdict can be overturned. The converse does not hold --
        % failing to find a witness proves nothing -- which is why this only ever vetoes.
            tf = false;
            if isempty(obj) || isempty(obj.ineqs)
                return
            end
            % px/py arrive as syms; the probe geometry is numeric, so pin them down here.
            try
                base = double([px(:), py(:)]);
            catch
                base = zeros(0,2);
            end
            if ~isempty(base)
                base = base(all(isfinite(base),2) & all(abs(base) ~= intmax, 2), :);
            end
            if isempty(base)
                base = [0 0];
            end
            ang = (0:11)' * (pi/6);
            C = [cos(ang), sin(ang)];
            for b = 1:size(base,1)
                for r = [0.05 0.5 5 50]
                    P = base(b,:) + r*C;
                    for k = 1:size(P,1)
                        if min(vecnorm(base - P(k,:), 2, 2)) < 1.0d-9
                            continue        % that IS a listed vertex; it witnesses nothing
                        end
                        if obj.ptFeasible(obj.vars, P(k,:))
                            tf = true;
                            return
                        end
                    end
                end
            end
        end

        function obj = simplifyUnboundedRegion (obj)
            
            if isempty(obj)
                return
            end
            nP = 0;
            for j = 1:obj.nv
              if abs(obj.vx(j)) == intmax
                continue
              end
              if abs(obj.vy(j)) == intmax
                continue
              end
              nP = nP+1;
              px(nP) = obj.vx(j);
              py(nP) = obj.vy(j);
            end
            
            for i = 1:size(obj.ineqs,2)
              for j = i+1:size(obj.ineqs,2)
                  if isAlways(obj.ineqs(i) == -obj.ineqs(j))
                      obj = region.empty();
                      return
                  end
              end
            
            end
            if nP == 0
                % NO FINITE VERTEX IS NOT THE SAME AS EMPTY. A half-plane, a slab and the whole
                % plane all have no corner, and this routine declared every one of them empty --
                % the same "a degenerate geometric object is not a geometric object" trap as
                % removeTangent's apex and getEdgeNosInf's lens, in the one direction that
                % silently deletes an answer.
                %
                % A HALF-PLANE IS EXACTLY WHAT A TANGENT VERTEX PRODUCES, so this is not a corner
                % case. getNormalConeVertexQ builds the cone at a vertex from the two edges meeting
                % there; when those edges are TANGENT -- an arc and its chord touching, which is
                % how a curvilinear piece ends -- both half-planes are the SAME one and the cone is
                % a half-plane with no vertex. Measured on piece 9 of f* for x*y over the
                % parallelogram conv{(0,0),(2,0),(2.5,1),(0.5,1)}: the cone at (1/4,7/8) is
                % {2x/3 + y >= 4/3}, its cell carries x/4 + 7y/8 - 1/2, and dropping it left the
                % conjugate uncovered on exactly that half-plane -- wrong or uncovered at 6 of 10
                % probe points against a brute-force sup, and the biconjugate over the whole
                % parallelogram then collapsed to nothing.
                %
                % Refute the verdict with a WITNESS rather than reversing it: a feasible point is a
                % certificate of non-emptiness (ptFeasible tests every constraint, conics
                % included), while failing to find one proves nothing -- so the old verdict stands
                % whenever no witness turns up, and this can only ever recover a region that was
                % being deleted. Nothing is simplified here: with no vertex there is no
                % vertex-based simplification to do, and over-describing a region is harmless
                % (this class's LP-certificate header says so).
                if ~obj.witnessAwayFrom(sym.empty(), sym.empty())
                    obj = region.empty();
                end
                return
            end
            obj = obj.simplifyOpenRegion1 (nP, px, py);
            
            if obj.nv == size(obj.ineqs,2)
                return
            end 
            if isempty(obj)
                return
            end
            
            for j = 1:nP
                nPoint(j) = 0;
            end
            keep = [];
            for i = 1:size(obj.ineqs,2)
                keep(i) = false;
            end

            % get Ineqs going through a vertex
            point = zeros(nP,0);   % nP rows up front: a vertex with NO active constraint would
                                   % otherwise leave `point` short and the ip-loop below would
                                   % index a row that does not exist (MATLAB:badsubscript).
            for i = 1:size(obj.ineqs,2)
                for j = 1:nP
                    % subsF returns NaN for a 0/0 substitution -- both numerator and denominator
                    % vanishing AT this vertex, which is a REMOVABLE singularity, not a missing
                    % value. Resolve it by the directional limit, exactly as vertexOfEdge does;
                    % isAlways(NaN == 0) is false, so without this the constraint is silently
                    % dropped and the vertex ends up with no active constraint at all.
                    fij = obj.ineqs(i).subsF(obj.vars,[px(j),py(j)]);
                    if isnan(fij.f)
                        fij = obj.ineqs(i).limitDirectional(obj.vars,[px(j),py(j)]);
                    end
                    if ~isnan(fij.f) && isAlways(fij.f == 0)
                    nPoint(j)=nPoint(j)+1;
                    point(j,nPoint(j)) = i;
                    end

                end
            end
            if all(nPoint == 0)
                obj.print
                nP
                px
                py
            end
            
             
            
            % make a function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if nP == 1 & nPoint(1)==2
                mp = obj.slopeAtVertex(point(1,1:2), [px,py]);
                theta = atan2(mp,1);
                t1 = theta(1);
                t2 = theta(2);
                ind = [1,2];
                d = tan((t1+t2)/2);
                
                if isnan(d) | isinf(d)
                    tx = px(1);
                    ty = py(1)+0.1;
                else
                    c = py(1) - d * px(1);
                    
                    tx = px(1) + 0.1;
                    ty = d*tx+c;
                end
                
                lm = obj.ptFeasible (obj.vars,[tx,ty]);
                if isnan(d)| isinf(d)
                    tx = px(1);
                    ty = py(1)-0.1;
                else
                
                    tx = px(1) - 0.1;
                    ty = d*tx+c   ;
                end
                
                lm = lm | obj.ptFeasible (obj.vars,[tx,ty]);
                d = -1/d;
                
                if isnan(d) | isinf(d)
                    tx = px(1);
                    ty = py(1)+0.1;
                else
                    c = py(1) - d * px(1);
                    
                    tx = px(1) + 0.1;
                    ty = d*tx+c;
                end
                
                lm = lm | obj.ptFeasible (obj.vars,[tx,ty]);
                if isnan(d)| isinf(d)
                    tx = px(1);
                    ty = py(1)-0.1;
                else
                
                    tx = px(1) - 0.1;
                    ty = d*tx+c   ;
                end
                
                lm = lm | obj.ptFeasible (obj.vars,[tx,ty]);
                 


                
               if ~lm && ~obj.witnessAwayFrom(px, py)
                   obj = region.empty();
               end

               return

            end


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if all(nPoint == 2)
                return;
            end
            
            % point
            % nPoint
            for ip = 1:nP
             %   ip
                
                sx = px(ip);
                sy = py(ip);
                pi0 = point(ip,1:nPoint(ip));
                markF = [];
                if size(pi0,2) == 2
                    keep(pi0) = true;
                    continue
                end
                if nPoint(ip) == 0
                    continue
                end
                

                mp = obj.slopeAtVertex(pi0, [sx,sy]);
                
                for i = 1:size(mp,2)
                    if mp(i) == -inf
                        mp(i) = inf;
                    end
                end

                %mp

                [sorted_mp, indices] = sort(double(mp));
                theta = atan2(double(mp),1);
                sorted_theta = atan2(sorted_mp,1);
                
              for i = 1:size(indices,2)
                 
                t1 = sorted_theta(i);
                 
                
                if i == size(indices,2)
                    t2 = sorted_theta(1);
                    ind = [pi0(indices(i)),pi0(indices(1))];
                    l90 = true;
                else
                    t2 = sorted_theta(i+1);
                    ind = [pi0(indices(i)),pi0(indices(i+1))];
                    l90 = false;
                end
                if t1 == t2
                    return
                    
                end
                
                if l90
                    d = tan(pi/2+ (t1+t2)/2);
                else
                    d = tan((t1+t2)/2);
                end
               
                if isnan(d) | isinf(d)
                    tx = sx;
                    ty = sy+0.1;
                else
                    c = sy - d * sx;
                    tx = sx + 0.1;
                    ty = d*tx+c;
                end
                %double(tx),double(ty) 
                %l = obj.ptFeasible (obj.vars,[tx,ty]);

                %if obj.ptFeasibleSubset (obj.vars,[tx,ty], ind)
                if obj.ptFeasible (obj.vars,[tx,ty])
                   lm = false;
                   %disp('c1')     
                   continue
                end
                if isnan(d)| isinf(d)
                    tx = sx;
                    ty = sy-0.1;
                else
                
                    tx = sx - 0.1;
                    ty = d*tx+c   ;
                end
                %double(tx),double(ty) 
                %l = obj.ptFeasible (obj.vars,[tx,ty])
                %if obj.ptFeasibleSubset (obj.vars,[tx,ty], ind)
                if obj.ptFeasible (obj.vars,[tx,ty])
                  lm = false;
                  %disp('c2')     
                  continue
                end
                


                
                if i == 1
                    %m2 = sorted_m(size(indices,2));
                    t2 = sorted_theta(size(indices,2));
                    ind = [pi0(indices(i)),pi0(indices(size(indices,2)))];
                    l90 = true;
                else
                    %m2 = sorted_m(i-1);
                    t2 = sorted_theta(i-1);
                    ind = [pi0(indices(i)),pi0(indices(i-1))];
                    l90 = false;
                end
                if t1 == t2
                    return
                end
                
                
                if l90
                    d = tan(pi/2+ (t1+t2)/2);
                else
                    d = tan((t1+t2)/2);
                end
 %               m1, m2, double(d)
         %       atan(m1)
         %       m1,m2,d
                % if i == 1
                %     d = -d;
                % end
               % d
                %isnan(d)
                if isnan(d)| isinf(d)
                    tx = sx;
                    ty = sy+0.1;
                else
                    c = sy - d * sx;
                    tx = sx + 0.1;
                    ty = d*tx+c;
                end
                
                %d,tx,ty
                %double(tx),double(ty) 
                %if obj.ptFeasibleSubset (obj.vars,[tx,ty], ind)
                if obj.ptFeasible (obj.vars,[tx,ty])
                   lm = false;
                   %disp('c3')     
                   continue
                end
                if isnan(d)| isinf(d)
                    tx = sx;
                    ty = sy-0.1;
                else
                
                    tx = sx - 0.1;
                    ty = d*tx+c   ;    
                end
                %double(tx),double(ty) 
                %if obj.ptFeasibleSubset (obj.vars,[tx,ty], ind)
                if obj.ptFeasible (obj.vars,[tx,ty])
                  lm = false;
                  %disp('c4')     
                  continue
                end
                %disp('end');
                %disp('marked')
                %i
             %   indices(i)
              %  pi0(indices(i))
                markF = [markF,pi0(indices(i))];
                
              end
              %markF
              for i = 1:size(obj.ineqs,2)
                  keep0(i) = false;
              end
              %keep0
              for i = 1:size(pi0,2)
                  keep0(pi0(i)) = true;
              end
              %keep0
              for i = 1:size(markF,2)
                  keep0(markF(i)) = false;
              end
              %keep0
              %keep
              keep = keep | keep0 ;
              %keep
              
            end
            %keep
            %markF
            %if ~isempty (markF)
                if size(markF,2) == size(obj.ineqs,2)
                %disp("Singleton1")
                %obj.print
                    % Every constraint's sector probes came back infeasible, which is taken to
                    % mean the region has collapsed to the vertex. Those probe DIRECTIONS come
                    % from constraint slopes, so they are meaningless for a constraint whose
                    % gradient vanishes at the vertex -- refute the verdict with an actual
                    % feasible point before acting on it. See witnessAwayFrom.
                    if obj.witnessAwayFrom(px, py)
                        return
                    end
                    obj = region.empty;
                   % disp("Singleton")
                    return

                end
            %end
            %obj.ineqs(markF) = [];
            %obj.ineqs = obj.ineqs(markF);
            % check vertices
            markF0 = [];
            for i = 1:size(obj.ineqs,2)
                if keep(i)
                    continue
                end
                markF0 = [markF0,i];
            end
            %markF0
            if size(markF0,2) == size(obj.ineqs,2)
                %disp("Singleton0")
                %obj.print
                if obj.witnessAwayFrom(px, py)   % same veto as above
                    return
                end
                obj = region.empty;
                %disp("Singleton")
                return

            end
            objTemp = obj;
            obj = obj.deleteIfRedundant(markF0);
            nv = obj.nv;
            obj = obj.getVertices;
            %obj = obj.removeTangent (nP, px, py);
            obj =  obj.removeInfV ;
            if nv ~= obj.nv
                disp ("Error in simplify");
                %sorted_m
                objTemp.print
                obj.print
            end


        end
        


        function obj = simplify (obj) 
          if isempty(obj)
              return
          end
          for i = 1:size(obj.ineqs,2)
            for j = i+1:size(obj.ineqs,2)
                if isAlways(obj.ineqs(i) == -obj.ineqs(j))
                    obj = region.empty();
                    return
                end
            end
          end
          for i = 1:size(obj.ineqs,2)
            [lelim, obj] = obj.deleteIneq (obj.vars);
            if ~lelim
                return
            end
          end

        end


        function [tx,ty] = getFeasiblePtNearV (obj, i)
            edges = obj.getEdges(obj.vx(i),obj.vy(i));
            mp = obj.slopeAtVertex(edges, [obj.vx(i),obj.vy(i)]);
            theta = atan2(mp,1);
            t1 = theta(1);
            t2 = theta(2);
            d = tan((t1+t2)/2);
                
            if isnan(d) | isinf(d)
              tx = obj.vx(i);
              ty = obj.vy(i)+0.1;
            else
              c = obj.vy(i) - d * obj.vx(i);
              tx = obj.vx(i) + 0.1;
              ty = d*tx+c;
            end
                
            
            if obj.ptFeasible (obj.vars,[tx,ty]);
              return
            end

            if isnan(d)| isinf(d)
              tx = obj.vx(i);
              ty = obj.vy(i)-0.1;
            else
              tx = obj.vx(i) - 0.1;
              ty = d*tx+c   ;
            end

            
            if obj.ptFeasible (obj.vars,[tx,ty]);
              return
            end

                
            d = -1/d;
                
            if isnan(d) | isinf(d)
              tx = obj.vx(i);
              ty = obj.vy(i)+0.1;
            else
              c = obj.vy(i) - d * obj.vx(i);
              tx = obj.vx(i) + 0.1;
              ty = d*tx+c;
            end
                
           
            if obj.ptFeasible (obj.vars,[tx,ty]);
              return
            end

            if isnan(d)| isinf(d)
              tx = obj.vx(i);
              ty = obj.vy(i)-0.1;
            else
              tx = obj.vx(i) - 0.1;
              ty = d*tx+c   ;
            end

           
            if obj.ptFeasible (obj.vars,[tx,ty]);
              return
            end

            
        end        

        %function simplifyClosedRegion
        %end
         function l = ptFeasible(obj, vars, point)
         % objective: true iff every point is feasible for every constraint (ineqs(i) <= 0).
         %
         % PERFORMANCE, and nothing else: this used to substitute and test ONE CONSTRAINT AT A
         % TIME, so a region with n constraints cost 2n round trips to the symbolic engine per
         % point. It is the single hottest routine in the conjugate -- profiling the two-face
         % unit square's f* put it at 35 s of 148, with 4022 of the run's 10809 `subs` calls --
         % and every one of those round trips is ~7 ms of interprocess overhead around a
         % substitution that is itself trivial.
         %
         % `subs` and `isAlways` are both ELEMENTWISE on a sym array, so substituting into the
         % whole constraint vector at once and testing it at once gives bit-identical results in
         % 2 round trips instead of 2n. Nothing here is approximated: the arithmetic is still
         % exact symbolic arithmetic on the same expressions, and `isAlways(v > 0)` on an array
         % decides each entry exactly as it did one at a time (including returning false, with
         % its usual warning, for an entry it cannot decide).
         %
         % The only behavioural difference is that the short-circuit is gone -- all constraints
         % are now evaluated even when the first already fails. That is a strictly better trade
         % here: n is small (a handful of constraints) and one batched round trip beats even a
         % single unbatched one.
           l = true;
           if isempty(obj.ineqs)
               return
           end
           g = [obj.ineqs.f];
           for j = 1:size(point,1)
               if any(isAlways(subs(g, vars, point(j,:)) > 0))
                   l = false;
                   return
               end
           end
         end

         function l = ptFeasibleSubset(obj, vars, point, index)
           l = true;
          % point
           
           
           for i = 1:size(index,2)
               %obj.ineqs(index(i)).f    
               for j = 1:size(point,1)
               %    subs ([obj.ineqs(index(i)).f],vars,point(j,:))
                if    (subs ([obj.ineqs(index(i)).f],vars,point(j,:)) > 0)
                %if double(subs ([obj.ineqs(i).f],vars,double(point(j,:)))) > 1.0e-12
                   l = false;
                   return
               end
               end
           end  
         end

         function obj = removeDenominator(obj)
           for i = 1:size(obj.ineqs,2)
              
              obj.ineqs(i) = obj.ineqs(i).removeDenominator2;
           end 

         end

         function l = isFeasibleWBPts(obj)
             l = false;
             for i = 1:obj.nv
               if ~ obj.ptFeasible(obj.vars, [obj.vx(i),obj.vy(i)])
                   return
               end
             end
             l = true;
         end

         function l = isFeasible(obj)
             l = false;
             for i = 1:size(obj.ineqs,2)
               for j = i+1:size(obj.ineqs,2)
                 if ( obj.ineqs(i) == unaryminus(obj.ineqs(j)))
                     return
                 end
                 if ( obj.ineqs(i) == obj.ineqs(j))
                     continue
                 end
                 
                 s = solve(obj.ineqs(i).f<=0,obj.ineqs(j).f<=0);
                 if isempty(s)
                     return;
                 end
               end 
 
             end 
             
             l = true;
         end

         % stupid way of doing this
         function obj = intersection (obj1, obj2)
         l = [];
           for i = 1:size(obj1.ineqs,2)
              l = [l,obj1.ineqs(i).f];
           end 
         for i = 1:size(obj2.ineqs,2)
             lf = true;
             for j = 1:size(obj1.ineqs,2)
               if obj1.ineqs(j)==obj2.ineqs(i)
                   lf = false;
                   break;
               end
             end
             if lf
               l = [l,obj2.ineqs(i).f];
             end
         end 
         obj = region(l, obj1.vars);
         if isempty(obj)
             return
         end
         %obj.print
         if (isFeasible(obj))
             %disp('feasible')
             obj = obj.unique;
             if obj.nv <= 2
                 %disp('degenerate ');
                 obj = region.empty;
             else
                 if obj.nv == 0
           disp('in intersection nv 0')
         end

             end
         else
             obj = region.empty;
         end

     end
     
     function obj = normalizeEdge(obj)
         for i = 1:size(obj.ineqs,2)  
             c = getLinearCoeffs (obj.ineqs(i),obj.vars);
             
             if (c(2) == 0)
                 % double * f not overloaded
                 obj.ineqs(i).f = (1/c(1)) * obj.ineqs(i).f;
             else
                 obj.ineqs(i).f = (1/c(2)) * obj.ineqs(i).f;
             end
         end
     end

     function obj = normalize1(obj)
         for i = 1:size(obj.ineqs,2)
             obj.ineqs(i) = obj.ineqs(i).normalize1;
         end
     end

     function l = isVertex(obj,V)
         l = true;
         for i = 1:size(V,1)
           for j = 1:obj.nv
             l1 = false;
             if (isAlways(obj.vx(j) == V(i,1))) & (isAlways(obj.vy(j) == V(i,2)))
                 l1 = true;
                 break
             end
           end
           if ~l1
               l = false;
               return
           end
         end
     end
     
     function V = getIntersectingFeasiblePts(obj, f0)
         
       nv=0;
       V = sym.empty();
       %obj.vx=[];
       %obj.vy=[];
       t1 = sym('t1');
       t2 = sym('t2');
       varsTemp = [t1,t2];
       f2 = f0.subsF (obj.vars,varsTemp);
       for i = 1:size(obj.ineqs,2)  
           f1 = obj.ineqs(i);
           f1 = f1.subsF (obj.vars,varsTemp);
               
               s = solve ([f1.f==0,f2.f==0],varsTemp);
               
               if isempty(s)
                   continue;
               elseif isempty(s.t1)
                   continue;
               elseif isempty(s.t2)
                   continue;
               end
               
               
               for k = 1:size(s.t1,1)
                    if ~ isreal(double(s.t1(k)))  % weird error - not detecting complex in symbolic z^4
                        continue
                    end
                    if ~ isreal(double(s.t2(k)))
                        continue
                    end
                    
                    if obj.ptFeasible(obj.vars, [s.t1(k),s.t2(k)])
                       
                          nv=nv+1;
                          V(nv,1) = s.t1(k);
                          V(nv,2) = s.t2(k);
                        
                    end
               end 
           
       end
     end
     
     % wont work for degree > 2
     % removed complex vertices
     
     function obj = linear3pt(obj)
       mark = [];  
      
       for i = 1:size(obj.ineqs,2) 
           
         if obj.ineqs(i).isQuad  
             continue
         end
        
         [nvi, vxi, vyi] = obj.vertexOfEdge(i);
         
         if nvi <= 2
             continue
         end
         % assuming nvi = 3
         [sortvx, ind] = sort(double(vxi));  % rational sorting goofed up
         x = (vxi(ind(1)) + vxi(ind(2)))/2;
         y = (vyi(ind(1)) + vyi(ind(2)))/2;
         %obj.ptFeasible (x,y) 
        
         if ~obj.ptFeasible (obj.vars,[x,y])
             mark(1) = ind(1);
         else
             mark(1) = ind(3);
         end
        
         break;

       end  
       
       if isempty(mark)
           return
       end
       %obj.print
       %mark
       % HISTORY: mark(1) comes out of the loop above as an index into vxi/vyi -- THIS EDGE's
       % own vertex list, at most 3 long -- and is then translated into an index into
       % obj.vx/obj.vy. That translation used to write its answer back into mark(1) while the
       % loop was still reading vxi(mark(1)) as the search key, so the first match turned the
       % key into an obj.vx index and the next iteration overran vxi ("Index exceeds the
       % number of array elements. Index must not exceed 3." -- f=xy over
       % conv{(0,0),(3,3),(1,2)}, reached from region.plus via functionNDomain.mtimes).
       % Capture the point first, stop at the first match, and delete nothing if the middle
       % collinear vertex is not one of this region's vertices after all (deleting an
       % arbitrary vertex was the old fallthrough).
       tx = vxi(mark(1));
       ty = vyi(mark(1));
       iv = 0;
       for i = 1:obj.nv
           if (obj.vx(i) == tx) & (obj.vy(i) == ty) %#ok<AND2> (sym comparison, matches the
               iv = i;                              % elementwise form used throughout region.m)
               break
           end
       end
       if iv == 0
           return
       end
       obj.vx(iv) = [];
       obj.vy(iv) = [];
       obj.nv = obj.nv-1;
     end 

     function obj = removeInfV (obj)
         n = 0;
         for i = 1:obj.nv
             
             if isAlways (abs(obj.vx(i)) == intmax)
                 continue
             end
             
             if isAlways (abs(obj.vy(i)) == intmax)
                 continue
             end
             n = n + 1;
             obj.vx(n) = obj.vx(i);
             obj.vy(n) = obj.vy(i);
             
         end
         obj.nv = n;
     end 

     function obj = getVertices(obj)

       obj.nv=0;
       % HISTORY: obj.vx/obj.vy must be reset here, not just obj.nv -- this
       % method is called more than once on the same (already-populated)
       % object (e.g. region.removeTangent calls it again after deleting
       % some ineqs). Without clearing them, the "vertex at infinity"
       % placeholder-point phase below (which always appends via
       % obj.nv=obj.nv+1; obj.vx(obj.nv)=..., never replaces) piles new
       % placeholder points on top of the ones a previous call already
       % added, instead of recomputing from scratch -- producing duplicate
       % vertices at the same coordinates. Those duplicates then make
       % region.getEndpoints return more than 2 points for an edge that
       % genuinely only has 2, which crashes region.minus's fixed-size
       % v1(i,:,:) assignment (testOpenconvex).
       obj.vx=sym.empty();
       obj.vy=sym.empty();
       t1 = sym('t1');
       t2 = sym('t2');
       varsTemp = [t1,t2];

       % PERFORMANCE: rename each constraint's variables to t1,t2 ONCE, not once per PAIR.
       % subsF is ~19 ms (two substitutions, two simplifyFractions and an isAlways each), and
       % the inner loop used to redo the j-th constraint's rename for every i < j -- n(n-1)/2
       % renames where n suffice. Identical results: the rename does not depend on i.
       nIn = size(obj.ineqs,2);
       fT = symbolicFunction.empty();
       for i = 1:nIn
           fT(i) = obj.ineqs(i).subsF (obj.vars,varsTemp);
       end
       % PERFORMANCE, part 2: two AFFINE constraints meet in one point given by a determinant,
       % and calling solve() to discover that costs ~30 ms of round trip per pair -- 322 of the
       % 438 solve() calls in a two-face conjugate come from right here. The coefficients are
       % read by EVALUATION at (0,0), (1,0), (0,1), exactly as region.linearForm does and for
       % the same reason (an affine g satisfies g = g(0,0) + (g(1,0)-g(0,0))t1 + (g(0,1)-g(0,0))t2
       % identically, however it happens to be written), and all three substitutions are batched
       % across every constraint at once: 3 round trips for the whole region.
       %
       % Everything stays EXACT -- these are sym values, not doubles, and the closed form is the
       % same rational arithmetic solve() would have done. The fast path is taken only for a pair
       % that is affine/affine with a NONZERO determinant; a curved constraint, parallel lines or
       % coincident lines all fall through to solve() and keep their existing behaviour, whatever
       % it is, rather than having a new answer invented for them here.
       affT = false(1,nIn); cf = sym(zeros(nIn,3));
       if nIn > 0
           GT = [fT.f];
           try
               e00 = subs(GT, varsTemp, [0 0]);
               e10 = subs(GT, varsTemp, [1 0]);
               e01 = subs(GT, varsTemp, [0 1]);
               for i = 1:nIn
                   % isLinear is degreeNum == 1 AND degreeDen == 0, so it already rules out a
                   % rational facet; no separate isPolynomial round trip needed.
                   if ~fT(i).isLinear
                       continue
                   end
                   cf(i,:) = [e10(i)-e00(i), e01(i)-e00(i), e00(i)];
                   affT(i) = true;
               end
           catch
               affT(:) = false;      % anything unevaluable: every pair uses solve(), as before
           end
       end
       for i = 1:nIn
           f1 = fT(i);
           for j = i+1:nIn
               f2 = fT(j);
               s = [];
               if affT(i) && affT(j)
                   dt = cf(i,1)*cf(j,2) - cf(j,1)*cf(i,2);
                   if ~isAlways(dt == 0)
                       s = struct('t1', (cf(i,2)*cf(j,3) - cf(j,2)*cf(i,3))/dt, ...
                                  't2', (cf(j,1)*cf(i,3) - cf(i,1)*cf(j,3))/dt);
                   end
               end
               if isempty(s)
                   s = solve ([f1.f==0,f2.f==0],varsTemp);
               end

               if isempty(s)
                   continue;
               elseif isempty(s.t1)
                   continue;
               elseif isempty(s.t2)
                   continue;
               end
               
               
               for k = 1:size(s.t1,1)
                    if ~ isreal(double(s.t1(k)))  % weird error - not detecting complex in symbolic z^4
                        continue
                    end
                    if ~ isreal(double(s.t2(k)))
                        continue
                    end
                    
                    %if obj.ptFeasible(obj.vars, [s.t1(k),s.t2(k)])
                    if obj.ptFeasible(obj.vars, [s.t1(k),s.t2(k)])
                       %i,j,obj.nv
                       %isAlways(subs ([obj.ineqs(i).f],obj.vars,[s.t1(k),s.t2(k)])<=0)
                       %double(subs ([obj.ineqs(i).f],obj.vars,[s.t1(k),s.t2(k)]))
                  
                          obj.nv=obj.nv+1;
                          obj.vx(obj.nv) = simplify(s.t1(k));
                          obj.vy(obj.nv) = simplify(s.t2(k));
                        
                    end
               end 
           end
               
       
       end
       
       if obj.nv ~= 0
           %disp('b4 sort')
           %obj.nv
           %obj.vx
           %double(obj.vx)
       for i = 1:obj.nv
           V(i,1) = obj.vx(i);
           V(i,2) = obj.vy(i);
       end 
       % double(obj.vx)
       % double(obj.vy)
       % obj.nv
       % V
       V = unique(V,"rows");
       if obj.nv ~= size(V,1)
         obj.nv = size(V,1);
         obj.vx = V(:,1)';
         obj.vy = V(:,2)';
       end
     
           obj = obj.linear3pt; 
          % obj.nv
       end  
       %[vx,vy] = poly2cw(obj.vx,obj.vy);
       %obj = obj.poly2orderUnbounded;

       % putting intmax for inf to avoid Nans 
       % intmax + intmax = intmax
       %disp('ptFeasile')
       %obj.print
       n = obj.nv;
       % if n== 0
       %     return
       % end
       %return
       % to remove inf - fix + in convex envelope solve
     
       if obj.ptFeasible(obj.vars, [intmax,intmax])
          obj.nv=obj.nv+1;
          obj.vx(obj.nv) = intmax;
          obj.vy(obj.nv) = intmax;
          %disp("++")
       end
       if obj.ptFeasible(obj.vars, [intmax,-intmax])
          obj.nv=obj.nv+1;
          obj.vx(obj.nv) = intmax;
          obj.vy(obj.nv) = -intmax;
          %disp("+-")
       end
       if obj.ptFeasible(obj.vars, [-intmax,intmax])
          obj.nv=obj.nv+1;
          obj.vx(obj.nv) = -intmax;
          obj.vy(obj.nv) = intmax;
          %disp("-+")
       end
       if obj.ptFeasible(obj.vars, [-intmax,-intmax])
          obj.nv=obj.nv+1;
          obj.vx(obj.nv) = -intmax;
          obj.vy(obj.nv) = -intmax;
          %disp("--")
       end

       % A SLAB -- bounded in one direction, infinite in the other -- reaches here with NOTHING.
       % Its constraints are parallel, so the pairwise phase above finds no intersection, and no
       % CORNER of the +-intmax box satisfies the two bounds at once, so none of the four tests
       % above fires either. nv stays 0, and the constructor reads that as "no vertices, hence
       % infeasible" and returns region.empty. That is simply wrong: {0 <= y <= 1} is a perfectly
       % good nonempty region, and it made every slab FACE unrepresentable -- measured, the
       % symptom surfaced two layers away as fanUnboundedFace refusing a slab because
       % recessionRays reports 'bounded' for an empty region.
       %
       % The box clip is the right idea; it was only ever applied at the box's CORNERS. Apply it
       % along the box's EDGES too: intersect each affine constraint with x = +-intmax and
       % y = +-intmax and keep the feasible hits. For the slab that yields (+-intmax, 0) and
       % (+-intmax, 1), exactly the four markers the rest of region already expects for an
       % unbounded shape.
       %
       % GUARDED BY nv == 0 deliberately, so that nothing which already produces a vertex can
       % change: this only ever runs where the alternative is a wrong `empty`.
       %
       % And only for a certified 2-DIMENSIONAL region. nv == 0 also covers a genuinely infeasible
       % set and a DEGENERATE one (a line or a point, e.g. {x <= 0, -x <= 0}), and both must keep
       % coming out empty: a measure-zero cell contributes nothing to a sup, but if it survives
       % into a max assembly it can carry a function on territory that is not its own -- the
       % over-claim failure mode this file's LP certificates exist to prevent.
       %
       % A polyhedron is 2-dimensional exactly when no constraint is an IMPLICIT EQUALITY, i.e.
       % when min a_i'z over the region is strictly below b_i for every i. That is the same LP
       % primitive maxLinear already provides. Anything the LP cannot certify -- a curved facet,
       % an infeasible set, a non-converged solve -- falls through to the previous behaviour
       % (leave nv at 0, hence empty) rather than guessing.
       if obj.nv == 0
           [Alp, blp, linlp] = obj.linearForm;
           full2D = ~isempty(Alp) && all(linlp);
           if full2D
               for i = 1:size(Alp,1)
                   [v, st] = region.maxLinear(Alp, blp, -Alp(i,:));
                   if st == 1
                       continue                     % min a_i'z = -inf: not an implicit equality
                   end
                   if st ~= 0 || abs(-v - blp(i)) <= 1e-9*max(1, abs(blp(i)))
                       full2D = false; break        % infeasible, undecided, or tight everywhere
                   end
               end
           end
       else
           full2D = false;
       end
       if full2D
           M = sym(intmax);
           for i = 1:nIn
               if ~affT(i)
                   continue                 % a curved facet has no single box-edge crossing
               end
               a = cf(i,1); b = cf(i,2); c = cf(i,3);
               cand = sym.empty(0,2);
               if ~isAlways(b == 0)
                   for sgn = [1 -1]
                       cand(end+1,:) = [sgn*M, -(a*sgn*M + c)/b]; %#ok<AGROW>
                   end
               end
               if ~isAlways(a == 0)
                   for sgn = [1 -1]
                       cand(end+1,:) = [-(b*sgn*M + c)/a, sgn*M]; %#ok<AGROW>
                   end
               end
               for k = 1:size(cand,1)
                   if obj.ptFeasible(obj.vars, cand(k,:))
                       obj.nv = obj.nv + 1;
                       obj.vx(obj.nv) = cand(k,1);
                       obj.vy(obj.nv) = cand(k,2);
                   end
               end
           end
       end

       for i = 1:n
           if nargin == 2
               %intmax,obj.vy(i)
               obj.ptFeasible(obj.vars, [intmax,obj.vy(i)]);
           end
           if obj.ptFeasible(obj.vars, [obj.vx(i),intmax])
             obj.nv=obj.nv+1;
             obj.vx(obj.nv) = obj.vx(i);
             obj.vy(obj.nv) = intmax;
           end
           if obj.ptFeasible(obj.vars, [obj.vx(i),-intmax])
             obj.nv=obj.nv+1;
             obj.vx(obj.nv) = obj.vx(i);
             obj.vy(obj.nv) = -intmax;
           end
           if obj.ptFeasible(obj.vars, [intmax,obj.vy(i)])
             obj.nv=obj.nv+1;
             obj.vx(obj.nv) = intmax;
             obj.vy(obj.nv) = obj.vy(i);
           end
           if obj.ptFeasible(obj.vars, [-intmax,obj.vy(i)])
             obj.nv=obj.nv+1;
             obj.vx(obj.nv) = -intmax;
             obj.vy(obj.nv) = obj.vy(i);
           end

       end

       % HISTORY: the "vertex at infinity" placeholder loop above can
       % legitimately re-derive a point that a separate mechanism already
       % put in obj.vx/obj.vy -- e.g. once some other operation has added an
       % explicit +-intmax boundary ineq (turning what was an unbounded
       % direction into a literal, very-large-but-finite line), the ordinary
       % pairwise-intersection phase earlier in this function can solve that
       % boundary ineq against another edge and land on the exact same
       % (x,+-intmax) point this per-vertex loop also generates as an
       % "at infinity" projection -- and only the finite-vertex subphase
       % above (before this loop) deduplicates via unique(rows); nothing
       % deduplicates the combined finite+infinite list once this loop has
       % run. The resulting duplicate vertices make region.getEndpoints
       % return more than 2 points for an edge that only has 2 distinct
       % ones, which crashes region.minus's fixed-size v1(i,:,:) assignment
       % (testOpenconvex). Dedup the full list here, mirroring the earlier
       % finite-only dedup.
       if obj.nv ~= 0
         V = [obj.vx(:), obj.vy(:)];
         V = unique(V,"rows");
         if obj.nv ~= size(V,1)
           obj.nv = size(V,1);
           obj.vx = V(:,1)';
           obj.vy = V(:,2)';
         end
       end

     end

     function [edgeNo] = getEdgeNosInf2(obj, vars)
         edgeNo = zeros(size(obj.ineqs,2),1)
         for i = 1:obj.nv
           e1 = obj.getEdges(obj.vx(i),obj.vy(i))
         end
         [nv1, vx, vy] = obj.vertexOfEdge(e1(1))
         [nv2, vx, vy] = obj.vertexOfEdge(e1(2))
         if nv1 == 1
           edgeNo(1) = e1(1) ;  
           edgeNo(2) = e1(2) ;  
           cEdge = e1(2);
         else
           edgeNo(2) = e1(1) ;  
           edgeNo(1) = e1(2) ;  
           cEdge = e1(1);
         end
         for i = 2:obj.nv
             e1 = obj.getEdges(obj.vx(i),obj.vy(i))
             if e1(1) == cEdge
                 edgeNo(i+1) = e1(2);
                 cEdge = e1(2);
             else
                 edgeNo(i+1) = e1(1);
                 cEdge = e1(1);
             end

         end
     end

     function [edgeNo] = getEdgeNosInf(obj, vars)
         edgeNo = zeros(size(obj.ineqs,2),1);
         e1 = obj.getEdges(obj.vx(1),obj.vy(1));

         
         % return
         % if obj.nv == 1 & size(obj.ineqs,2) ==2
         %     edgeNo = [1,2];
         %     return
         % end
         [nv1, vx, vy] = obj.vertexOfEdge(e1(1));
         [nv2, vx, vy] = obj.vertexOfEdge(e1(2));
         if nv1 == 1 | nv2 == 1
             add = 1;
             edgeNo(1) = 1 ;
         else
             add = 0;
         end
         %add
         for i = 1:size(obj.ineqs,2)
            %obj.ineqs(i).print
           [nv, vx, vy] = obj.vertexOfEdge(i);
           % nv == 0 -- a constraint with NO vertex on this region -- bounds no EDGE of it, so it
           % has no edge number, and 0 reports that. Without this, vx is empty and vx(start)
           % below overruns ('Index exceeds array bounds'). Unreachable while
           % simplifyUnboundedRegion deleted every constraint missing a finite vertex; reachable
           % now that deletion needs a redundancy certificate (see region.redundantSubset).
           %
           % Callers must handle the 0 rather than feed it to an index -- see
           % functionNDomain.conjugateOfPiecePoly, which parks such constraints above every real
           % edge slot.
           if nv == 0
               edgeNo(i) = 0;
               continue
           end
           start = 1;
           %% change condition
           if nv > 1
               if isAlways(vx(2) == obj.vx(obj.nv)) & isAlways(vy(2) == obj.vy(obj.nv)) & ~ (isAlways(vx(1) == obj.vx(obj.nv-1)) & isAlways(vy(1) == obj.vy(obj.nv-1)))
                  start = 2;
               end
           end
            for j = 1:obj.nv
                  if obj.vx(j)==vx(start) & obj.vy(j)==vy(start)
                      break;
                  end    
            end
            edgeNo(i) = j +add;
           
         end
         if nv1 == 1
             edgeNo(e1(1)) = 1;
         elseif nv2 == 1
             edgeNo(e1(2)) = 1;
         end
         
     end


     function [edgeNo] = getEdgeNos(obj, vars)
         
       edgeNo = zeros(obj.nv,1);
       for j = 1: obj.nv-1
         slope = obj.slope(j,j+1);
         if slope == -inf
           slope = inf;
         end
         if slope == inf
           edge = vars(1) -obj.vx(j) ; 
         else
           q = obj.yIntercept (j,slope);
           edge = vars(2)-slope*vars(1)-q;
         end
         edge = symbolicFunction(edge);
         edge = edge.normalize(vars);
         %edge
         for k = 1: size(obj.ineqs,2)
           e0 = obj.ineqs(k);
           e0 = e0.normalize(vars);
           %e0 = e0.normalize (vars);
           if (isAlways(e0.f == edge.f) | isAlways(e0.f == -edge.f))
             break;
           end
         end
         edgeNo(j)=k;
       end
       j = obj.nv;
       slope = obj.slope(j,1);
       if slope == inf
         edge = vars(1) -obj.vx(j)  ;
       else
         q = obj.yIntercept (j,slope);
         edge = vars(2)-slope*vars(1)-q;
       end
       edge = symbolicFunction(edge);
       edge = edge.normalize(vars);
       %edge
       for k = 1: size(obj.ineqs,2)
         e0 = obj.ineqs(k);
         e0 = e0.normalize (vars);
         if (isAlways(e0.f == edge.f) | isAlways(e0.f == -edge.f))
           break;
         end
       end
       edgeNo(j)=k;
     end

      function [nv, vx, vy] = vertexOfEdge(obj,ind)
              nv = 0;
              vx = sym.empty();
              vy = sym.empty();
              %disp('vertex of edge')
              %obj.nv
              %obj.ineqs(ind)
               % HISTORY: subsF signals a 0/0 substitution by returning NaN (see its den==0
               % branch), and isZero -> getDen -> numden(NaN) then errors out of MuPAD with
               % "Invalid return value 'undefined'". That is not an edge case here: a rational
               % constraint whose numerator and denominator both vanish AT one of the region's own
               % vertices is the normal shape of these envelope/conjugate pieces -- the
               % singularity is removable, and the value is the limit. funcVertices already uses
               % exactly this NaN-then-limit idiom; this site simply did not.
               for i = 1:obj.nv
                 f1 = obj.ineqs(ind).subsF (obj.vars,[obj.vx(i),obj.vy(i)]);
                 if isnan(f1.f)
                     % 0/0 at this vertex. subsF flags that with NaN; the value is the LIMIT.
                     % limitDirectional, not limit: the latter iterates univariate limits and
                     % returns NaN again for a bivariate 0/0 (see symbolicFunction.limit).
                     f1 = obj.ineqs(ind).limitDirectional (obj.vars,[obj.vx(i),obj.vy(i)]);
                 end
                 if isnan(f1.f)
                     continue   % genuinely no limit: the vertex is not on this edge
                 end
                 if f1.isZero
                     nv = nv + 1;
                     vx(nv) = obj.vx(i);
                     vy(nv) = obj.vy(i);
                 end
               end
      end

      function edges = getEdges(obj,vx,vy)
          % HISTORY: `edges` was only ever assigned INSIDE the loop, so a point with no active
          % constraint did not return an empty list -- it raised MATLAB:unassignedOutputs,
          % "Output argument edges ... not assigned". That is not a hypothetical point: the
          % vertex list of an unbounded region contains the box-clip vertices getVertices
          % appends, and the corner (intmax,intmax) of the first quadrant {-x<=0,-y<=0} lies on
          % the implicit +-intmax box rather than on either constraint, so nothing is active
          % there. "No constraint is active at this point" is a perfectly good answer and callers
          % can test numel(edges); erroring instead made poly2orderUnbounded unusable on the
          % simplest unbounded region there is.
          edges = [];
          n = 0;
          for i = 1:size(obj.ineqs,2)
              %obj.ineqs(i)
              %obj.ineqs(i).subsF(obj.vars,[vx,vy])
              %isZero(obj.ineqs(i).subsF(obj.vars,[vx,vy]))
              if isZero(obj.ineqs(i).subsF(obj.vars,[vx,vy]))
                  n = n + 1;
                  edges(n) = i;
              end
          end
          
      end

      function f = divideRegions(obj,obj2)

       obj.nv=0;
         t1 = sym('t1');
         t2 = sym('t2');
         varsTemp = [t1,t2];
         ir = 0;
         % in this loop we are getting all points of intersection of pair
         % of ineqs.
         % remove duplicate points 
         for i = 1:size(obj.ineqs,2)  
           f1 = obj.ineqs(i);
           f1 = f1.subsF (obj.vars,varsTemp);
           for j = i+1:size(obj.ineqs,2)  
               f2 = obj.ineqs(j);
               f2 = f2.subsF (obj.vars,varsTemp);
               s = solve ([f1.f==0,f2.f==0],varsTemp);

               if isempty(s)
                   continue;
               elseif isempty(s.t1)
                   continue;
               elseif isempty(s.t2)
                   continue;
               end
               lpt = false;
               for k  = 1:ir
                   if all(points(k,:) == [s.t1,s.t2])
                       lpt = true;
                       break;
                   end
               end
               if lpt
                   continue
               end
               if ~obj2.ptFeasible (obj2.vars,[s.t1,s.t2])
                 continue
               end

               ir = ir + 1;
               index(ir) = 0;
               points(ir,1:2) = [s.t1,s.t2]; 
               for k = 1:size(obj.ineqs,2)  
                 if double(subs ([obj.ineqs(k).f],obj.vars,[s.t1,s.t2])) < 1.0e-12
                     index(ir) = index(ir)+1;
                     r(ir,index(ir)) = obj.ineqs(k).f;
                 end
               end

           end

         end
 
          is = 0;
         s = r;
         for i=1:ir
                 leq = false;
                 for j = 1:is
                     if indexS(j) ~= index(i)
                         continue;
                     end
                     leq = true;
                     for k = 1:index(i)
                         if ~ isequal(r(i,k),s(j,k))
                             leq = false;
                             break;
                         end
                     end
                     if leq
                       break;
                     end
                 end
                 if leq
                    continue;
                 end
                 is = is+1;

                 indexS(is) = index(i);
                 for j = 1:indexS(is)
                     s(is,j) = r(i,j);
                 end

         end

         f = [];
         for i = 1:is
           f = [f, region(s(i,:),obj.vars)]; 
         end

       end

  
 
     function fv = funcVertices (obj, f)
         n = 0;
         for i =1:obj.nv
             if (abs(obj.vx(i))==intmax)
                 continue
             end
             if (abs(obj.vy(i))==intmax)
                 continue
             end
             n = n + 1;
             fv(n) = f.subsF(obj.vars,[obj.vx(i),obj.vy(i)]);
             if (isnan(fv(n).f))
                 fv(n) = f.limit(obj.vars,[obj.vx(i),obj.vy(i)]);
             end
         end
     end


     function [l, maxf, index, lSing] = maximum(obj, f)
         
         f(1).f = simplifyFraction(f(1).f);
         f(2).f = simplifyFraction(f(2).f);
         % f(1) == f(2)
         % isAlways(f(1)==f(2))
          if isAlways(f(1)==f(2))
              l = true;
              lSing = false;
              maxf = f(1);
              index = 1;
              return
          end
          [l,  maxf, index, lSing] = obj.maxArray (f(1), f(2)) ;
     end

     function [f2, fe, d] = splitMax(obj, f, expr)
       [nl, v1, v2] = obj.splitmaxArray (f(1), f(2)) 
       f2 = [];
       fe = [];
       d = [];
       v{1}=v1;
       v{2}=v2;
       for i = 1:2
           if nl(i) == 0
               continue
           end
           f2 = [f2,f(i)];
           fe = [fe,expr(i)];
           d0 = obj.regionWPts(obj.vx(v{i}), obj.vy(v{i}), obj.vars(1), obj.vars(2));
           d = [d,d0];
           % disp('Regiopn')
           % d0.print
       end
       
     end
     
     function e = getOtherEdgeAtVertex (obj, ind, vertex)
         e = 0;
         %disp('vertex')
         %vertex
         for i = 1:size(obj.ineqs,2)
            if i == ind
                continue
            end
            f1 = obj.ineqs(i).subsF (obj.vars,vertex);
            %double(obj.ineqs(i).subsF (obj.vars,vertex).f)
            if f1.isZero
                
                e = i;
                break
            end

         end
     end


     function l = isconvex (obj, obj2, i, j, vx, vy)
       l = false;
      
       edgeiNo = obj.getOtherEdgeAtVertex (i,[vx,vy]);
       if edgeiNo == 0
          return
       end
       % if ~obj.ineqs(edgeiNo).isLinear 
       %   return
       % end
       edgejNo = obj2.getOtherEdgeAtVertex (j,[vx,vy]);
       if edgejNo == 0
         return
       end
       % if ~obj2.ineqs(edgejNo).isLinear 
       %   return
       % end
       
       if isAlways (obj.ineqs(edgeiNo) == obj2.ineqs(edgejNo))
           l = true;
           return
       end

       if obj.slopeIneq(edgeiNo,[vx,vy]) == inf
           x1 = vx;
           y1 = vy + 0.1;
           if ~ obj.ptFeasible (obj.vars,[x1,y1])
             y1 = vy - 0.1;
         
           end
       
       else

       x1 = vx+0.1;
       f1 = obj.ineqs(edgeiNo).subsF([obj.vars(1)],[x1]);
       y0 = solve(f1, obj.vars(2));
       %if (size(y0,1) == 1)
           y1 = y0(1);
       %end
       
       if ~ obj.ptFeasible (obj.vars,[x1,y1])
         x1 = vx-0.1;
         f1 = obj.ineqs(edgeiNo).subsF([obj.vars(1)],[x1]);
         y0 = solve(f1, obj.vars(2));
       %if (size(y0,1) == 1)
           y1 = y0(1);
       %end
       
       end
       end
       %x1,y1
       if obj2.slopeIneq(edgejNo,[vx,vy]) == inf
           x2 = vx;
           y2 = vy + 0.1;
           if ~ obj2.ptFeasible (obj2.vars,[x2,y2])
             y2 = vy - 0.1;
         
           end
       
       else

       
       x2 = vx+0.1;
       f1 = obj2.ineqs(edgejNo).subsF([obj2.vars(1)],[x2]);
       y0 = solve(f1, obj2.vars(2));
       %if (size(y0,1) == 1)
           y2 = y0(1);
       %end
       if ~ obj2.ptFeasible (obj2.vars,[x2,y2])
         x2 = vx-0.1;
         f1 = obj2.ineqs(edgejNo).subsF([obj2.vars(1)],[x2]);
         y0 = solve(f1, obj2.vars(2));
         %if (size(y0,1) == 1)
           y2 = y0(1);
         %end
       
       end
       end
       %x2,y2
       xm = (x1 + x2)/2;
       ym = (y1 + y2)/2;
       if obj.ptFeasible (obj.vars,[xm,ym]) | obj2.ptFeasible (obj2.vars,[xm,ym])
          l = true;
       end
       %disp('is convex')
       %l
      end


     function [l,obj] = merge (obj, obj2)
     % Union two regions that share a facet, by deleting the facet from both and intersecting
     % what is left. Writing obj = A = A' n {g<=0} and obj2 = B = B' n {g>=0} for the shared
     % constraint g, the recipe returns M = A' n B'.
     %
     % That recipe is not unconditionally the union, and used to be applied unconditionally.
     % M never LOSES a point -- any x in A' n B' has g(x) <= 0, putting it in A, or g(x) >= 0,
     % putting it in B -- but it can GAIN points that belonged to neither operand, and the
     % merged region then carries its function value on territory that was never its own. For
     % f=xy over conv{(0,0),(3,3),(1,2)}, three same-valued Step 3 regions merged into one
     % covering s=(1,1), where none of them does, and the assembled partition returned 1.0
     % instead of 1.125.
     %
     % M = A u B exactly when A subset B' and B subset A' -- equivalently, when A u B is
     % convex -- and unionIsExact now decides that by LP before any facet is deleted. See the
     % LP-certificate block at the top of this file for why the test is an LP and why refusing
     % (leaving the two regions separate, which is always correct) is the safe answer.
     %
     % HISTORY, worth keeping: guarding this alone -- with the far stronger "the two constraint
     % sets must be EQUAL after the facet deletion" -- was tried, is provably sound, and made
     % the measured result WORSE, 36 -> 125 wrong points of 289. Refusing merges leaves more
     % regions standing, and simplifyUnboundedRegion was at the time deleting non-redundant
     % constraints from every one of them. That is why the two were fixed together, and it is
     % why the weaker, exact condition here is the right one: it refuses only what it must.
         l = false;
         % HISTORY: an empty (0x0) region can reach here from an upstream
         % copy-through that didn't check isempty after a simplification
         % step revealed a domain to be infeasible (see functionNDomain.
         % mergeL/maximumP's own HISTORY comments -- this was unreachable
         % before functionNDomain.maximumP stopped discarding lSing pairs).
         % obj.ineqs/obj2.ineqs on a 0x0 region array below has no elements
         % to index, which errors rather than behaving like an empty set;
         % merging with an empty region is a no-op, so just say so.
         if isempty(obj) || isempty(obj2)
             region.mergeTally('emptyOperand');
             return
         end
         n = 0;
         lQuad1 = false;
         nmq1 = 0;
         for i =1:size(obj.ineqs,2)
             if obj.ineqs(i).isQuad
                 lQuad1 = true;
                 nmq1 = nmq1 + 1;
                 mq1(nmq1) = i;
             end
         end
         lQuad2 = false;
         nmq2 = 0;
         %obj2
         %size(obj2.ineqs,2)
         for i =1:size(obj2.ineqs,2)
             %obj2.ineqs(i)
             if obj2.ineqs(i).isQuad
                 lQuad2 = true;
                 nmq2 = nmq2 + 1;
                 mq2(nmq2) = i;
             end
         end
         if (lQuad1 & lQuad2)
             marki = [];
             markj = [];
             for i = 1:nmq1
               for j = 1:nmq2
                   if obj.ineqs(mq1(i)) == -obj2.ineqs(mq2(j))
                       n = n + 1;
                       marki(n) = mq1(i);
                       markj(n) = mq2(j);
                   end
               end
             end
             % Exactly ONE shared facet, and the union must actually be the intersection of
             % what is left. Both conditions are needed and neither is cosmetic. With two
             % shared facets g1,g2 the "M never loses a point" argument breaks outright: a
             % point with g1<=0 and g2>=0 lies in neither operand yet survives into M. And
             % without unionIsExact the single-facet case is the over-claiming defect this
             % function's header describes.
             if n == 1
               [okU, whyU] = obj.unionIsExact(obj2, marki(1), markj(1));
               if okU
                 l = true;
                 obj3 = obj;
                 obj.ineqs(marki) = [];
                 obj2.ineqs(markj) = [];
                 obj = obj+obj2;
                 if isempty(obj)
                     l = false;
                     obj=obj3;
                     region.mergeTally('quadFacetEmptyUnion');
                 else
                     region.mergeTally('okQuadFacet');
                 end
                 return
               end
               region.mergeTally(['quadFacet_' whyU]);
               return
             end
             if n > 1
               region.mergeTally('quadFacetMultiShared');
               return          % shares a quadratic facet, but not exactly/convexly: no merge
             end
          end
         % HISTORY 2026-08-17 -- TWO HEURISTIC REFUSALS REMOVED, and they were the blow-up.
         % `quadMismatch` demanded, as a CROSS PRODUCT, that every quadratic of one region equal
         % every quadratic of the other; `quadCutsOther` refused when one region's conic met the
         % other anywhere but at a vertex. Neither proves the union is not convex -- two adjacent
         % cells each carrying a different parabolic arc ELSEWHERE have a perfectly convex union
         % and were refused by the first outright. Both existed only because unionIsExact had no
         % certificate for a CURVED constraint and refused those outright too.
         %
         % MEASURED on the all-rational control case in `.claude/step3cost.m` (x*y over
         % conv{(0,0),(3,3),(1,2)} via convEnvCPLQ + ratPolToPlq, no surds and no doubles
         % anywhere): at fold 3 the two refused 322 of 612 attempts while 4 succeeded, and
         % Step 3's cell count reached 57 for 10 DISTINCT FUNCTIONS.
         %
         % region.certifiesNonPositive supplies the missing certificate, and unionIsExact is the
         % EXACT criterion -- M = A' n B' equals A u B precisely when A subset B' and B subset A'
         % -- so dropping a conservative pre-check cannot make a wrong merge happen: every
         % candidate still has to pass it. What survives here is only the computation of lQuad,
         % which the single-vertex shared-facet case below still consults; a mismatch now leaves
         % lQuad TRUE, so that case stays as closed as it was.
         lQuad = lQuad1 | lQuad2;
         if lQuad1 & lQuad2
             lSameQuads = true;
             for i = 1:nmq1
               for j = 1:nmq2
                   if obj.ineqs(mq1(i)).f ~= obj2.ineqs(mq2(j)).f
                       lSameQuads = false;
                   end
               end
             end
             if lSameQuads
                 lQuad = false;
             end
         end
         marki = [];
         markj = [];
         for i =1:size(obj.ineqs,2)

           for j =1:size(obj2.ineqs,2)
             if obj.ineqs(i) == -obj2.ineqs(j)
                 [nvi, vxi, vyi] = obj.vertexOfEdge(i); % move to outer loop
                 [nvj, vxj, vyj] = obj2.vertexOfEdge(j);
                 if nvi ~= nvj
                     continue
                 end
                 if nvi == 1 & (~lQuad)
                     if obj.isconvex (obj2, i, j, vxi(1), vyi(1))
                     l = true;
                     n = n + 1;
                     marki(n) = i;
                     markj(n) = j;
                     end
                     break;
                 end
                 if nvi == 2
                 [tvxi, sortedIndices] = sort(abs(vxi));
                 vxi = vxi(sortedIndices);
                 vyi = vyi(sortedIndices);
                 [tvxj, sortedIndices] = sort(abs(vxj));
                 vxj = vxj(sortedIndices);
                 vyj = vyj(sortedIndices);
                 if all(vxi == vxj) & all(vyi == vyj)
                     if ~obj.isconvex (obj2, i, j, vxi(1), vyi(1))
                       continue
                     end
                    if nvi == 2   
                        if ~obj.isconvex (obj2, i, j, vxi(2), vyi(2))
                            continue
                        end
                    end
                    l = true;
                       n = n + 1;
                       marki(n) = i;
                       markj(n) = j;
                       break;
                     end
                 end  
             end
           end
         end
         % isconvex above is a LOCAL probe: it midpoint-tests a step either side of the shared
         % facet's endpoints, which is necessary for the union to be convex but nowhere near
         % sufficient -- it says nothing about the two regions' other constraints, and those
         % are what decide whether A' n B' over-claims. Same single-facet requirement as the
         % quadratic branch above, for the same reason.
         if ~l
           region.mergeTally('noSharedFacet');
         elseif n > 1
           l = false;
           region.mergeTally('multiSharedFacet');
         else
           [okU, whyU] = obj.unionIsExact(obj2, marki(1), markj(1));
           if ~okU
             l = false;
             region.mergeTally(['lin_' whyU]);
           else
             obj3 = obj;
             obj.ineqs(marki) = [];
             obj2.ineqs(markj) = [];
             obj = obj+obj2;
             obj = obj.simplifyUnboundedRegion;
             if isempty(obj)
                 disp('empty')
                 l = false;
                 obj = obj3;
                 region.mergeTally('linEmptyUnion');
             else
                 region.mergeTally('okLinear');
             end
           end
         end
     end



     function l = hasNegativeIneqs(obj)
         l = true;
         for i = 1:size(obj.ineqs,2)
           for j = i+1:size(obj.ineqs,2)
               if (obj.ineqs(i) == -obj.ineqs(j))
                   return;
               end
           end
         end
         l = false;
     end

     end

     methods % limits
       function [l,limf] = limitOfFAtVertices (obj, f)
       % PREALLOCATED AS SYM, and that line is the whole point of this comment. `limf` used to
       % be built by assignment alone, so its CLASS was decided by whichever branch ran first --
       % and the `limf(j) = 0` branch below writes a DOUBLE. One vertex where the two iterated
       % limits disagree therefore turned the array double, and every exact value written after
       % it was silently rounded on the way in.
       %
       % MEASURED 2026-08-17 on the A.4 sub-triangle of the general convex quadrilateral: the
       % envelope's gradient is 0/0 at the first vertex, so `l` came back [0 1 1] and the exact
       % gradient limits at the other two arrived as 0.70778 and 0.81131 --
       % 7307585874000779/9007199254740992 and its neighbour, 2^53 denominators. Those go
       % straight into getSubdiffVertexT1/T2, hence into every conjugate cell of the piece.
       %
       % Why that is a defect and not an inefficiency: two cells that SHARE A FACET can then
       % carry two different doubles of the same exact number, and region.merge finds a facet by
       % asking whether one constraint is the negation of another. Measured, one ULP apart, on
       % `4 - 2*sqrt(2)`. The facet becomes invisible and Step 3's cell count grows without
       % bound. Same defect as domain.mE/cE; DECISIONS.md 2026-08-17 has the chain.
         vars = obj.vars;
         vars2 = [vars(2),vars(1)];
         l = false(1, obj.nv);
         limf = sym(zeros(1, obj.nv));
         for j = 1: obj.nv

           l1 = f.limit(vars,[obj.vx(j),obj.vy(j)]);
           l2 = f.limit(vars2,[obj.vy(j),obj.vx(j)]);

           if (isAlways(l1.f == l2.f))
             l(j) = true;
             limf(j) = l1.f;
           else
             l(j) = false;
             limf(j) = sym(0);
           end
         end
       end
       %limf
     end

     methods % normal cone
         function NC = getNormalConeVertexQ(obj, s1, s2, eIdx)
             % obj = obj.envelope(i).d
             %
             % EXPLICIT EDGE LIST (optional 4th argument). Without it this routine reads the
             % constraint bounding an edge off a SLOT: the cone at vertex j is built from
             % ineqs(j) and ineqs(j+1), wrapped modulo the number of constraints. That works
             % only while "constraint j bounds edge j" holds, which is the correspondence
             % conjugateOfPiecePoly's scatter maintains -- and which it cannot maintain for a
             % LENS, where two genuine edges join the SAME pair of vertices and therefore
             % collide on one edge number (see functionNDomain.edgeIndexList).
             %
             % Given eIdx, edge j is the edge from vertex j to vertex j+1 (cyclically) and its
             % constraint is ineqs(eIdx(j)) -- so every index below that stood for "the
             % constraint of edge j" is looked up in the list instead. Nothing else changes:
             % hand this eIdx = 1:nv and it is the routine it was, index for index.
            if nargin < 4
                eIdx = [];
            end
            NC = sym(zeros(obj.nv,2));
            meanx = sum(obj.vx)/obj.nv;
            meany = sum(obj.vy)/obj.nv;

            % Assuming edges are ordered
            % V1-V2 is edge 1
            % obj.print
            %disp('normalConeQ')
            vars = obj.vars;
             for j = 1: obj.nv
                % WHICH EDGE IS "THE FIRST CONSTRAINT AT VERTEX j". This routine's own probe
                % points settle it: the first half probes at vertex j-1 and the second at vertex
                % j+1, so the first half is the edge ARRIVING at vertex j and the second the edge
                % LEAVING it. In eIdx's convention edge e runs from vertex e to vertex e+1, so the
                % arriving edge is e = j-1 and the leaving edge is e = j.
                cj = j;
                if ~isempty(eIdx)
                    cj = eIdx(mod(j-2, obj.nv) + 1);
                end
                slope = obj.slopeIneq(cj,[obj.vx(j),obj.vy(j)]);
                pslope = -1/slope;
                if pslope == -inf
                    pslope = inf;
                end
                if pslope ~= inf
                    q = obj.yIntercept (j,pslope);
                    eq = s2 - pslope*s1 - q;
                else
                    eq = s1 - obj.vx(j);
                end
                eq = symbolicFunction(eq);
                if obj.nv > 1
                    k = j-1;
                    % With an explicit edge list the region is CLOSED -- edge nv runs from the
                    % last vertex back to the first -- so vertex 1's predecessor is vertex nv,
                    % not "off the end". Wrapping here is what makes the k<1 probe construction
                    % below unnecessary in that case; without a list nothing changes.
                    if k < 1 && ~isempty(eIdx)
                        k = obj.nv;
                    end
                    if k < 1

                      vs = obj.ineqs(cj).getVars();
                  if size(vs,2) == 1 & isAlways(vs(1) == vars(2))
                      py = obj.vy(j);
                      px = obj.vx(j)+0.1 ;
                      if ~obj.ptFeasible(obj.vars,[double(px),double(py)])
                          px = obj.vx(j)-0.1 ;
                      end
                      if eq.subsF([s1,s2],[px,py]).isZero
                          py = obj.vy(j)+0.1;
                          if ~obj.ptFeasible(obj.vars,[double(px),double(py)])
                            py = obj.vy(j)-0.1 ;
                          end
                      end
                  elseif size(vs,2) == 1 & isAlways(vs(1) == vars(1))
                      py = obj.vy(j)+0.1 ;
                      px = obj.vx(j);
                      
                      if ~obj.ptFeasible(obj.vars,[double(px),double(py)])
                          py = obj.vy(j)-0.1 ;
                      end
                      % added as next eq can be 0 eg. x-4, y+5
                      %disp('zero check')
                      %eq.subsF([s1,s2],[px,py]).isZero
                      if eq.subsF([s1,s2],[px,py]).isZero
                          px = obj.vx(j)+0.1;
                       %   obj.ptFeasible(obj.vars,[double(px),double(py)])
                          if ~obj.ptFeasible(obj.vars,[double(px),double(py)])
                            px = obj.vx(j)-0.1 ;
                          end
                        %  obj.ptFeasible(obj.vars,[double(px),double(py)])
                        %  disp('done')
                      end
                  else
                      px = obj.vx(j) - 0.1;
                      ey = subs(obj.ineqs(cj).f,obj.vars(1),px);
                      py = solve(ey,obj.vars(2));
                      % py can have >1 root when obj.ineqs(j) is quadratic in y (this is
                      % getNormalConeVertexQ, the curved/quadratic-edge case) -- reduce to a
                      % single candidate BEFORE the isempty check (indexing an empty sym array
                      % errors), matching the same fix already applied a few lines below in this
                      % same function (see HISTORY: found via testMaxMultiRegion/testMax, a
                      % plq.biconjugateF call, throwing horzcat:dimensionMismatch when py had 2
                      % elements here but px stayed a scalar).
                      if isempty(py)
                          py = obj.vy(j);
                      else
                          py = py(1);
                      end
                      %obj.ptFeasible(obj.vars,[double(px),double(py)])
                      if ~obj.ptFeasible(obj.vars,[double(px),double(py)])
                          px = obj.vx(j) + 0.1;
                          ey = subs(obj.ineqs(cj).f,obj.vars(1),px);
                          py = solve(ey,obj.vars(2));
                          if isempty(py)
                            py = obj.vy(j);
                          else
                            py = py(1);
                          end
                      end
                      %obj.ptFeasible(obj.vars,[double(px),double(py)])
                  end
                  %[tx,ty] = getFeasiblePtNearV (obj, j);
                else
                  px = obj.vx(k);
                  py = obj.vy(k);
                  %obj.ptFeasible(obj.vars,[double(px),double(py)])
                  %obj.ptFeasible(obj.vars,[px,py])
                end
                if eq.subsF([s1,s2],[px,py]).isZero
                    kk = region.probeVertexIndex(k, j, obj.nv);
                    ckk = kk;
                    if ~isempty(eIdx)
                        ckk = eIdx(kk);
                    end
                    px = obj.vx(kk) - 0.1;
                      ey = subs(obj.ineqs(ckk).f,obj.vars(1),px);
                      py = solve(ey,obj.vars(2));
                      % see the matching HISTORY comment above: reduce a possibly-multi-root py
                      % (quadratic-in-y ineqs(kk)) to a single candidate before use.
                      if isempty(py)
                          py = obj.vy(kk);
                      else
                          py = py(1);
                      end
                      if ~obj.ptFeasible(obj.vars,[double(px),double(py)])
                          px = obj.vx(kk) + 0.1;
                          ey = subs(obj.ineqs(ckk).f,obj.vars(1),px);
                          py = solve(ey,obj.vars(2));
                          if isempty(py)
                            py = obj.vy(kk);
                          else
                            py = py(1);
                          end
                      end
                end
                qd = 0;
                  if isAlways(eq.subsF([s1,s2],[px,py]).f < 0)
                    eq = -eq;
                  end
                end
                % end
                NC(j,1) = eq.f;

                % THE SECOND CONSTRAINT AT THIS VERTEX, WRAPPED. This routine's convention is
                % that vertex j lies on constraints j and j+1, which holds for the UNBOUNDED
                % layout conjugateOfPiecePoly builds (slot 1 reserved for the ray, so there are
                % nv+1 slots for nv vertices). On a BOUNDED region there are exactly nv slots and
                % the last vertex's second constraint is slot 1, not slot nv+1 -- indexing it
                % unwrapped raised MATLAB:badsubscript, which is why the only caller used to send
                % every bounded region to the POLYHEDRAL getNormalConeVertex instead. That is
                % wrong for a bounded region with a curved edge: the polyhedral routine builds the
                % cone from the CHORD between vertices, and the chord is not the boundary there.
                % Wrapping is a strict generalisation -- for nIn = nv+1 and j <= nv it is exactly
                % j+1, so no previously-working input changes.
                nIn = size(obj.ineqs,2);
                jNext = mod(j, nIn) + 1;
                cNext = jNext;
                if ~isempty(eIdx)
                    % The edge LEAVING vertex j, which in eIdx's convention is edge j itself;
                    % jNext stays what the probe below wants, the next VERTEX.
                    jNext = mod(j, obj.nv) + 1;
                    cNext = eIdx(j);
                end
                slope = obj.slopeIneq(cNext,[obj.vx(j),obj.vy(j)]);
                pslope = -1/slope;
                if pslope == -inf
                    pslope = inf;
                end
                if pslope ~= inf
                    q = obj.yIntercept (j,pslope);
                    eq = s2 - pslope*s1 - q;
                else
                    eq = s1 - obj.vx(j);
                end
                eq = symbolicFunction(eq);
                if obj.nv > 1
                k = jNext;
                if k >obj.nv
                    % if size(vs,2) == 1 & isAlways(vs(1) == 's_2')
                    %   py = obj.vy(j) - 0.1;
                    %   px = obj.vx(j);
                    %   if ~obj.ptFeasible(obj.vars,[px,py])
                    %       py = obj.vy(j) + 0.1;
                    % 
                    %   end
                    % else
                    %   px = obj.vx(j) - 0.1;
                    %   ey = subs(obj.ineqs(jNext).f,obj.vars(1),px);
                    %   py = solve(ey,obj.vars(2));
                    %   if ~obj.ptFeasible(obj.vars,[px,py])
                    %     px = obj.vx(j) + 0.1;
                    %     ey = subs(obj.ineqs(j).f,obj.vars(1),px);
                    %     py = solve(ey,obj.vars(2));
                    %     size(py)
                    %     if size(py,1) > 1
                    %       py = py(1);
                    %     end
                    %     if isempty(py)
                    %       py = obj.vy(j);
                    %     end
                    %   end
                    % end


                   vs = obj.ineqs(cNext).getVars();

                   if size(vs,2) == 1 & isAlways(vs(1) == vars(2))
                      py = obj.vy(j);
                      px = obj.vx(j)+0.1 ;
                      if ~obj.ptFeasible(obj.vars,[px,py])
                          px = obj.vx(j)-0.1 ;
                      end
                   elseif size(vs,2) == 1 & isAlways(vs(1) == vars(1))
                      
                      py = obj.vy(j)+0.1 ;
                      px = obj.vx(j);
                      obj.ptFeasible(obj.vars,[px,py]);
                      if ~obj.ptFeasible(obj.vars,[px,py])
                          py = obj.vy(j)-0.1 ;
                      end
                   else
                      px = obj.vx(j) - 0.1;
                      ey = subs(obj.ineqs(cNext).f,obj.vars(1),px);
                      py = solve(ey,obj.vars(2));
                      % see the matching HISTORY comment above: reduce a possibly-multi-root py
                      % (quadratic-in-y ineqs(j+1)) to a single candidate before use.
                      if isempty(py)
                          py = obj.vy(j);
                      else
                          py = py(1);
                      end
                      if ~obj.ptFeasible(obj.vars,[double(px),double(py)])
                          px = obj.vx(j) + 0.1;
                          ey = subs(obj.ineqs(cj).f,obj.vars(1),px);
                          py = solve(ey,obj.vars(2));
                          % HISTORY: this read `py = py(1); if isempty(py) ... end`, indexing
                          % BEFORE the emptiness guard, so an empty solve threw
                          % "Index exceeds array bounds" and the guard below it was dead code.
                          % The identical block a few lines above orders these correctly; this
                          % copy did not. Reachable from functionNDomain.conjugateOfPiecePoly
                          % once a Case C conjugate is rich enough to get this far.
                          if isempty(py)
                              py = obj.vy(j);
                          else
                              py = py(1);
                          end
                      end
                    end

                else
                  px = obj.vx(k);
                  py = obj.vy(k);
                end
                
                %isAlways(subs(eq,[s1,s2],[px,py]) == 0)
                if eq.subsF([s1,s2],[px,py]).isZero
                    kk = region.probeVertexIndex(k, j, obj.nv);
                    ckk = kk;
                    if ~isempty(eIdx)
                        ckk = eIdx(kk);
                    end
                    px = obj.vx(kk) - 0.1;
                      ey = subs(obj.ineqs(ckk).f,obj.vars(1),px);
                      py = solve(ey,obj.vars(2));
                      % see the matching HISTORY comment above: reduce a possibly-multi-root py
                      % (quadratic-in-y ineqs(kk)) to a single candidate before use.
                      if isempty(py)
                          py = obj.vy(kk);
                      else
                          py = py(1);
                      end
                      if ~obj.ptFeasible(obj.vars,[double(px),double(py)])
                          px = obj.vx(kk) + 0.1;
                          ey = subs(obj.ineqs(ckk).f,obj.vars(1),px);
                          py = solve(ey,obj.vars(2));
                          % HISTORY: this read `py = py(1); if isempty(py) ... end`, indexing
                          % BEFORE the emptiness guard -- the same inversion already corrected
                          % twice elsewhere in this function; the guard below it was dead code
                          % and an empty solve threw "Index exceeds array bounds" instead.
                          if isempty(py)
                            py = obj.vy(kk);
                          else
                            py = py(1);
                          end
                      end
                end
                
                
                  if isAlways(eq.subsF([s1,s2],[px,py]).f < 0)
                    eq = -eq;
                  end
                
                end

                
                NC(j,2) = eq.f;
                
                

             end
             
             if obj.nv == 1
                 [tx,ty] = getFeasiblePtNearV (obj, 1);
                 j = 1;
                 if isAlways(subs(NC(j,1),[s1,s2],[tx,ty]) < 0)
                         NC(j,1) = - NC(j,1);
                     end
                     if isAlways(subs(NC(j,2),[s1,s2],[tx,ty]) < 0)
                         NC(j,2) = - NC(j,2);
                     end
             end
                
        end

         function NC = getNormalConeEdgeQE(obj, s1, s2, eIdx)
         % EDGE normal cones from an EXPLICIT edge list. Edge j joins vertex j to vertex j+1
         % (cyclically) and is bounded by ineqs(eIdx(j)).
         %
         % This is getNormalConeEdgeQ3 and getNormalConeEdgeQ with their one difference removed.
         % Both build the same two half-planes per edge -- the perpendicular to the EDGE'S OWN
         % constraint at each of the edge's two endpoints, oriented by the other endpoint -- and
         % differ only in WHICH slot they believe that constraint sits in: ineqs(j) for a region
         % with nv constraints, ineqs(j+1) for one with nv+1 (slot 1 reserved for a ray). Neither
         % convention can express a LENS, where two edges join the same pair of vertices; see
         % functionNDomain.edgeIndexList. Hand this eIdx = 1:nv and it is Q3; hand it 2:nv and it
         % is Q, edge for edge.
            NC = sym(zeros(max(numel(eIdx), obj.nv), 3));
            for j = 1:numel(eIdx)
                a = j;
                b = mod(j, obj.nv) + 1;
                NC(j,1) = region.coneNormalAt(obj, eIdx(j), a, b, s1, s2);
                NC(j,2) = region.coneNormalAt(obj, eIdx(j), b, a, s1, s2);
            end
         end

         function NC = getNormalConeEdgeQ(obj, s1, s2)
             % obj = obj.envelope(i).d
            NC = sym(zeros(obj.nv,3));
            meanx = sum(obj.vx)/obj.nv;
            meany = sum(obj.vy)/obj.nv;

            % Assuming edges are ordered
            % V1-V2 is edge 1
            %disp('in edge')
             for j = 1: obj.nv-1
                 k = j+1;
                 %disp('f')
                 %obj.vx(j),obj.vy(j)
                slope = obj.slopeIneq(k,[obj.vx(j),obj.vy(j)]);
                pslope = -1/slope;
                if pslope == -inf
                    pslope = inf;
                end
                if pslope ~= inf
                    q = obj.yIntercept (j,pslope);
                    eq = s2 - pslope*s1 - q;
                else
                    eq = s1 - obj.vx(j);
                end
                if isAlways(subs(eq,[s1,s2],[obj.vx(k),obj.vy(k)]) > 0)
                    eq = -eq;
                end
                NC(j,1) = eq;
                %disp('s')
                %obj.vx(k),obj.vy(k)
                slope = obj.slopeIneq(k,[obj.vx(k),obj.vy(k)]);
                pslope = -1/slope;
                if pslope == -inf
                    pslope = inf;
                end
                if pslope ~= inf
                    q = obj.yIntercept (k,pslope);
                    eq = s2 - pslope*s1 - q;
                else
                    eq = s1 - obj.vx(k);
                end
                if isAlways(subs(eq,[s1,s2],[obj.vx(j),obj.vy(j)]) > 0)
                    eq = -eq;
                end
                NC(j,2) = eq;

             end
             
                
         end

 function NC = getNormalConeEdgeQ3(obj, s1, s2)
             % obj = obj.envelope(i).d
            NC = sym(zeros(obj.nv,3));
            meanx = sum(obj.vx)/obj.nv;
            meany = sum(obj.vy)/obj.nv;

            % Assuming edges are ordered
            % V1-V2 is edge 1
            %disp('in edge')
             for j = 1: obj.nv-1
                 
                 k = j+1;
                 %disp('f')
                 %obj.vx(j),obj.vy(j)
                slope = obj.slopeIneq(j,[obj.vx(j),obj.vy(j)]);
                pslope = -1/slope;
                if pslope == -inf
                    pslope = inf;
                end
                if pslope ~= inf
                    q = obj.yIntercept (j,pslope);
                    eq = s2 - pslope*s1 - q;
                else
                    eq = s1 - obj.vx(j);
                end
                if isAlways(subs(eq,[s1,s2],[obj.vx(k),obj.vy(k)]) > 0)
                    eq = -eq;
                end
                NC(j,1) = eq;
                %disp('s')
                %obj.vx(k),obj.vy(k)
                slope = obj.slopeIneq(j,[obj.vx(k),obj.vy(k)]);
                pslope = -1/slope;
                if pslope == -inf
                    pslope = inf;
                end
                if pslope ~= inf
                    q = obj.yIntercept (k,pslope);
                    eq = s2 - pslope*s1 - q;
                else
                    eq = s1 - obj.vx(k);
                end
                if isAlways(subs(eq,[s1,s2],[obj.vx(j),obj.vy(j)]) > 0)
                    eq = -eq;
                end
                NC(j,2) = eq;

             end
             j = obj.nv;
             k = 1;
                 %disp('f')
                 %obj.vx(j),obj.vy(j)
                slope = obj.slopeIneq(j,[obj.vx(j),obj.vy(j)]);
                pslope = -1/slope;
                if pslope == -inf
                    pslope = inf;
                end
                if pslope ~= inf
                    q = obj.yIntercept (j,pslope);
                    eq = s2 - pslope*s1 - q;
                else
                    eq = s1 - obj.vx(j);
                end
                if isAlways(subs(eq,[s1,s2],[obj.vx(k),obj.vy(k)]) > 0)
                    eq = -eq;
                end
                NC(j,1) = eq;
               % disp('s')
               % obj.vx(k),obj.vy(k)
                slope = obj.slopeIneq(j,[obj.vx(k),obj.vy(k)]);
                pslope = -1/slope;
                if pslope == -inf
                    pslope = inf;
                end
                if pslope ~= inf
                    q = obj.yIntercept (k,pslope);
                    eq = s2 - pslope*s1 - q;
                else
                    eq = s1 - obj.vx(k);
                end
                if isAlways(subs(eq,[s1,s2],[obj.vx(j),obj.vy(j)]) > 0)
                    eq = -eq;
                end
                NC(j,2) = eq;
                
        end

        function NC = getNormalConeVertex(obj, s1, s2)
             % obj = obj.envelope(i).d
            
            NC = sym(zeros(obj.nv,2));
            meanx = sum(obj.vx)/obj.nv;
            meany = sum(obj.vy)/obj.nv;
            %obj.vx
            %obj.vy
             for j = 1: obj.nv-1
                slope = obj.slope(j,j+1);
                pslope = -1/slope;
                if pslope == -inf
                    pslope = inf;
                end
                if pslope ~= inf
                    q = obj.yIntercept (j,pslope);
                    eq = s2 - pslope*s1 - q;
                else
                    eq = s1 - obj.vx(j);
                end
                eq = simplifyFraction(eq);
                if isAlways(subs(eq,[s1,s2],[obj.vx(j+1),obj.vy(j+1)]) < 0)
                    eq = -eq;
                end
                NC(j,1) = eq;
                if pslope ~= inf
                    q = obj.yIntercept (j+1,pslope);
                    eq = s2 - pslope*s1 - q;
                else
                    eq = s1 - obj.vx(j+1);
                end
                eq = simplifyFraction(eq);
                if isAlways(subs(eq,[s1,s2],[obj.vx(j),obj.vy(j)]) < 0)
                    eq = -eq;
                end
                %disp('here')

                j+1;
                NC(j+1,2) = eq;
             end
             j = obj.nv;
             % if nv = 1, slope calculation is wrong - temp fix put in
             % place 

             % if j == 1
             % else
                 % We fill NC (nv,1) and NC (1,2)
                slope = obj.slope(j,1);
                pslope = -1/slope;
                if pslope == -inf
                    pslope = inf;
                end
                
                if pslope ~= inf
                   q = obj.yIntercept (j,pslope);
                   eq = s2 - pslope*s1 - q;
                else
                   eq = s1 - obj.vx(j);
                end
                if isAlways(subs(eq,[s1,s2],[obj.vx(1),obj.vy(1)]) < 0)
                    eq = -eq;
                end
              
                NC(j,1) = eq;
                if pslope ~= inf
                    q = obj.yIntercept (1,pslope);
                    eq = s2 - pslope*s1 - q;
                else
                    eq = s1 - obj.vx(1);
                end
                if isAlways(subs(eq,[s1,s2],[obj.vx(j),obj.vy(j)]) < 0)
                    eq = -eq;
                end
                
                NC(1,2) = eq;
             % end
            
                
        end

        function [NC] = getNormalConeEdge(obj, s1, s2)
          NC = sym(zeros(obj.nv,2));
          for j = 1: obj.nv-1
            slope = obj.slope(j,j+1);
            pslope = -1/slope;
            if pslope == -inf
              pslope = inf;
            end
            if pslope ~= inf
              q = obj.yIntercept (j,pslope);
              eq = s2 - pslope*s1 - q;
            else
              eq = s1 - obj.vx(j);
            end
            if isAlways(subs(eq,[s1,s2],[obj.vx(j+1),obj.vy(j+1)]) > 0    )
              eq = -eq;
            end
            NC(j,1) = eq;
            if pslope ~= inf
              q = obj.yIntercept (j+1,pslope);
              eq = s2 - pslope*s1 - q;
            else
              eq = s1 - obj.vx(j+1);
            end
            if isAlways(subs(eq,[s1,s2],[obj.vx(j),obj.vy(j)]) > 0   )
              eq = -eq;
            end
            NC(j,2) = eq;
          end
          % added for unbounded
          if obj.nv ~= size(obj.ineqs,2)
              return
          end
          j = obj.nv;
          slope = obj.slope(j,1);
          pslope = -1/slope;
          if pslope == -inf
            pslope = inf;
          end
          if pslope ~= inf
            q = obj.yIntercept (j,pslope);
            eq = s2 - pslope*s1 - q;
          else
            eq = s1 - obj.vx(j);
          end
          if isAlways(subs(eq,[s1,s2],[obj.vx(1),obj.vy(1)]) > 0   )
            eq = -eq;
          end
          NC(j,1) = eq;
          if pslope ~= inf
            q = obj.yIntercept (1,pslope);
            eq  = s2 - pslope*s1 - q;
          else
            eq = s1 - obj.vx(1);
          end
          if isAlways(subs(eq,[s1,s2],[obj.vx(j),obj.vy(j)]) > 0   )
            eq = -eq;
          end
          NC(j,2) = eq;
        end

        function [NC] = getNormalConeEdgeQ0(obj, s1, s2)
          NC = sym(zeros(obj.nv,3));
          for j = 1: obj.nv-1
            slope = obj.slope(j,j+1)
            pslope = -1/slope;
            if pslope == -inf
              pslope = inf;
            end
            if pslope ~= inf
              q = obj.yIntercept (j,pslope);
              eq = s2 - pslope*s1 - q;
            else
              eq = s1 - obj.vx(j);
            end
            if isAlways(subs(eq,[s1,s2],[obj.vx(j+1),obj.vy(j+1)]) > 0    )
              eq = -eq;
            end
            NC(j,1) = eq;
            if pslope ~= inf
              q = obj.yIntercept (j+1,pslope);
              eq = s2 - pslope*s1 - q;
            else
              eq = s1 - obj.vx(j+1);
            end
            if isAlways(subs(eq,[s1,s2],[obj.vx(j),obj.vy(j)]) > 0   )
              eq = -eq;
            end
            NC(j,2) = eq;
            obj.ineqs(j+1)
            if obj.ineqs(j+1).isQuad
                mx = (obj.vx(j)+obj.vx(j+1))/2
                my = (obj.vy(j)+obj.vy(j+1))/2
                obj.ineqs(j+1).subsF(obj.vars,[mx,my])
                obj.ineqs(j+1).subsF(obj.vars,[mx,my])<0
                if isAlways(obj.ineqs(j+1).subsF(obj.vars,[mx,my])<0)
                   %  tang = obj.ineqs(j+1).tangentOfSlope (slope)
                   c = sym('c')
                   
                   f = obj.ineqs(j+1).subsF(obj.vars(2),slope*obj.vars(1)+c);
                   ans = solve(f.f,obj.vars(1));
                
                   f = ans(1)-ans(2);
                   c1 = solve(f,c);
                   %tang 
                   eq2 = s2-slope*s1-c1;
                else
                   eq2 = s2 - obj.vy(j) - slope*(s1 - obj.vx(j));
                end
                NC(j,3) = eq2;
                
            end
            
          end
          % added for unbounded
          if obj.nv ~= size(obj.ineqs,2)
              return
          end
          j = obj.nv;
          slope = obj.slope(j,1);
          pslope = -1/slope;
          if pslope == -inf
            pslope = inf;
          end
          if pslope ~= inf
            q = obj.yIntercept (j,pslope);
            eq = s2 - pslope*s1 - q;
          else
            eq = s1 - obj.vx(j);
          end
          if isAlways(subs(eq,[s1,s2],[obj.vx(1),obj.vy(1)]) > 0   )
            eq = -eq;
          end
          NC(j,1) = eq;
          if pslope ~= inf
            q = obj.yIntercept (1,pslope);
            eq  = s2 - pslope*s1 - q;
          else
            eq = s1 - obj.vx(1);
          end
          if isAlways(subs(eq,[s1,s2],[obj.vx(j),obj.vy(j)]) > 0   )
            eq = -eq;
          end
          NC(j,2) = eq;
        end


     function NC = adjustNormalConeUnB(obj, NC, edgeNo, s1, s2)

       if size(obj.ineqs,2) == obj.nv 
           return
       end
       %NCOP = sym(zeros(size(obj.ineqs,2),2));
       %NC(2:obj.nv+1,:) = NC(1:obj.nv,:)     
       % linear ineq so passing arbitrary point
       slope = obj.slopeIneq(edgeNo(1),[0,0]);
       pslope = -1/slope;
       if pslope == -inf
         pslope = inf;
       end
       if pslope ~= inf
         q = obj.yIntercept (1,pslope);
         eq = s2 - pslope*s1 - q;
       else
         eq = s1 - obj.vx(1);
       end
       eq
       [tx,ty] = getFeasiblePtNearV (obj, 1)
       isAlways(subs(eq,[s1,s2],[tx,ty]) < 0   )
       if isAlways(subs(eq,[s1,s2],[tx,ty]) < 0   )
            eq = -eq;
       end
       NC (1,1) = eq;
       %NCOP (1,2) = 0;
       %size(edgeNo,1)
       slope = obj.slopeIneq(edgeNo(size(edgeNo,1)),[0,0]);
       pslope = -1/slope;
       if pslope == -inf
         pslope = inf;
       end
       if pslope ~= inf
         q = obj.yIntercept (1,pslope);
         eq = s2 - pslope*s1 - q;
       else
         eq = s1 - obj.vx(1);
       end
       eq
       isAlways(subs(eq,[s1,s2],[tx,ty]) < 0   )
       if isAlways(subs(eq,[s1,s2],[tx,ty]) < 0   )
            eq = -eq;
       end
       %NCOP (obj.nv+1,1) = 0;
       NC (obj.nv,2) = eq;

       %NC = NCOP
       
         return
       if size(obj.ineqs,2) > obj.nv
           NC(2:obj.nv+1,:) = NC(1:obj.nv,:)     
             

           edge = edgeNo(1)
           j = 1;
           % linear ineq so passing arbitrary point
           slope = obj.slopeIneq(edge,[0,0])
           pslope = -1/slope;
           if pslope == -inf
             pslope = inf;
           end
           if pslope ~= inf
            q = obj.yIntercept (1,pslope);
            eq = s2 - pslope*s1 - q;
          else
            eq = s1 - obj.vx(j);
          end
          if subs(eq,[s1,s2],[obj.vx(1),obj.vy(1)]) > 0   
            eq = -eq;
          end
          NC(1,1) = sym(0);
          NC(1,2) = eq;

% make routine 
          edge = edgeNo(obj.nv)
          j = obj.nv;
           % linear ineq so passing arbitrary point
           slope = obj.slopeIneq(edge,[0,0])
           pslope = -1/slope;
           if pslope == -inf
             pslope = inf;
           end
           if pslope ~= inf
            q = obj.yIntercept (1,pslope);
            eq = s2 - pslope*s1 - q;
          else
            eq = s1 - obj.vx(j);
          end
          if subs(eq,[s1,s2],[obj.vx(1),obj.vy(1)]) > 0   
            eq = -eq;
          end
          NC(obj.nv+1,1) = eq;
          NC(obj.nv+1,2) = 0;
          
       end
     end


     % function expr = conjugateExprEdgesT1Poly2 (obj, f, vars)
     %        expr = []
     %        for i = 1:obj.nv
     %           expr = [expr,conjugateExpr(obj.ineqs(i).f,f.f,vars(1),vars(2))]
     %        end
     %    end
     end
     
        
end