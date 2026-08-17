% rename to piecewise function5

% piece
% piecewise function
% list of piecewise functions

classdef functionNDomain
    properties
      f = symbolicFunction.empty();
      d = region.empty();
    end
% 15 methods
     methods
         function obj = functionNDomain(f, d)
              % disp('in')
              % f
              % 
             obj.f = f;
              % d
             obj.d = d;
              % disp('out')
         end


         


         function print (obj)
             disp("Function")
             obj.f.printL;
             disp("Domain")
             obj.d.print;
         end

         function printLatex1 (obj)
             disp("Function")
             obj.f.printLatexWB;
             disp("Domain")
             obj.d.printLatex;
         end

         function printL(objL)
             for i = 1:size(objL,2)
                 i
                 objL(i).print;
             end
         end

         function printLatex (objL)
             for i = 1:size(objL,2)
                 i
                 objL(i).printLatex1;
             end
         end

         function printLLatex(objL)
             for i = 1:size(objL,2)
                 objL(i).printLatex;
             end
         end
         function printM(objL)
             colorList = ["red","blue","yellow","green","purple","cyan","orange","brown","crimson", "pink","tan","aquamarine","navy", "palegreen"];
             fprintf("display(inequal({");
             f = objL(1).f.f;
             j = 1;
             for i = 1:size(objL,2)
                 if objL(i).f.f == f
                   objL(i).d.printMaple;
                   fprintf(",");
                   
                 else
                   fprintf("},x=-15..15,y=-15..15,color=[");
                   fprintf(colorList(j)) ;
                   fprintf(",");
                   fprintf(colorList(j)); 
                   fprintf(",");
                   fprintf(colorList(j)); 
                   fprintf("],nolines),inequal({")  ;
                   f = objL(i).f.f;
                   j = j+1;
                   objL(i).d.printMaple;
                   fprintf(","); 
                 end
                 
             end
             fprintf("},x=-15..15,y=-15..15,color=[");
                   fprintf(colorList(j)) ;
                   fprintf(",");
                   fprintf(colorList(j)); 
                   fprintf(",");
                   fprintf(colorList(j)); 
                   fprintf("],nolines))")  ;
             fprintf("\n");
             
         end

         function printM2(objL)
             %1. AliceBlue. 2. AntiqueWhite. 3. Aquamarine. 4. Azure. 5. Beige ; 7. Black. 8. BlanchedAlmond. 9. Blue. 10. BlueViolet.
             colorList = ["red","blue","yellow","green","purple","cyan","orange","brown","crimson", "pink","tan","aquamarine","navy", "palegreen"];
             fprintf("display(inequal({");
             %f = objL(1).f.f;
             j = 1;
             
             for i = 1:size(objL,2)
                 
                 
                 
              %   if objL(i).f.f == f
                   objL(i).d.printMaple;
                   fprintf(",");
                   
              %   else
                   fprintf("},x=-15..15,y=-15..15,color=[");
                   fprintf(colorList(j)) ;
                   if i == size(objL,2)
                       fprintf("],nolines))")  ;
                       break
                   end
                   fprintf("],nolines),inequal({")  ;
                  % f = objL(i).f.f;
                   j = j+1;
                   if j > 14 
                       j = 1;
                   end
                %   objL(i).d.printMaple;
                 %  fprintf(","); 
               %  end
                 
             end
            
             fprintf("\n");
             
         end

          function plotDomain(obj)
             figure;
             colors = ['b', 'r', 'g', 'm', 'c', 'y'];
             n = 0
             f = obj (1).f
             c = colors(mod(n,6)+1)
             for i =1:size(obj,2)
                i
                if (f.f ~= obj (i).f.f)
                  n = n + 1
                  c = colors(mod(n,6)+1)
                  f = obj (i).f
                end
                obj (i).d.plot;
                textR = "R"+num2str(i);
                textR="";
                obj (i).d.plotRegionC(textR,c);
             end
          end

          % for closed regions and objL1 = objL2
          function [objL,index] = times (objL1, objL2)
             n = 0;
             objL=functionNDomain.empty();
             for i = 1:size(objL1,2)
                 markedi(i) = false;
             end
             for i = 1:size(objL1,2)
               for j = i+1:size(objL1,2)
                %   i,j
                 rf = objL1(i).d + objL2(j).d;
                 if size(rf.ineqs,2) < 3
                     rf = region.empty;
                 end
                  if isempty(rf)
                   continue
                 end
                 
                
                 % disp('b4 simplify')
                 % rf.print
                 %rf = rf.simplify;
                 if rf.nv <= 2
                     rf = region.empty;
                 end
                  if isempty(rf)
                   continue
                  end
                  rf = rf.simplify;
                 markedi(i) = true;
                 markedi(j) = true;
                 % i,j
                 % rf.print
                 n = n + 1;
                 objL(n) = functionNDomain([objL1(i).f(1), objL2(j).f(1)],rf);
                 index(n,1:2) = [i,j];

                   r = objL1(i).d - rf;
                 if ~ isempty(r)
                   for k = 1:size(r,2) 
                     if size(r(k).ineqs,2) < 3
                         continue
                     end
                     if r(k).nv < 3
                         continue
                     end
                     
                     n = n + 1;  
                     %  for iv = 1:r(k).nv
                     % 
                     % 
                     %      if (r(k).nv==3 & size(r(k).ineqs,2)==4)
                     %         disp('error1')
                     %         objL1(i).d.print
                     %         rf.print
                     %         r(k).print
                     %     end
                     % 
                     %      if (abs(r(k).vx(iv))==intmax)
                     %         disp('open in minus')
                     %         r(k).print
                     %         objL1(i).d.print
                     %         rf.print
                     % 
                     %     end
                     %     if (abs(r(k).vy(iv))==intmax)
                     %         disp('open in minus')
                     %         r(k).print
                     %         objL1(i).d.print
                     %         rf.print
                     % 
                     %     end
                     % end
                     
                     objL(n) = functionNDomain([objL1(i).f(1)],r(k));
                     index(n,1) = [i];
                   end
                 end
                  
                  r = objL2(j).d - rf;
                 
                 if ~ isempty(r)
                   for k = 1:size(r,2)  
                           if size(r(k).ineqs,2) < 3
                         continue
                     end
                     if r(k).nv < 3
                         continue
                     end
                 
                     n = n + 1  ;
                     
                         % disp("t2")
                         % r(k).print
                     % for iv = 1:r(k).nv
                     %     if (r(k).nv==3 & size(r(k).ineqs,2)==4)
                     %         disp('error2')
                     %         objL2(i).d.print
                     %         rf.print
                     %         r(k).print
                     %     end
                     %     if (abs(r(k).vx(iv))==intmax)
                     %         disp('open in minus')
                     %         r(k).print
                     %     end
                     %     if (abs(r(k).vy(iv))==intmax)
                     %         disp('open in minus')
                     %         r(k).print
                     %     end
                     % end
                     objL(n) = functionNDomain([objL2(j).f(1)],r(k));
                     index(n,1) = [j];
                   end
                 end
               end
             end
             for i = 1:size(objL1,2)
                 if markedi(i) 
                     continue
                 end
                 % disp('in times')
                 % i
                 n = n + 1;
                 % objL1(i).print
                 objL(n) = objL1(i);
                 index(n,1) = [i];

             end
             
                      
          end


          function [objR,index2] = maximumPC(objL, index) %, f, r2)
           n = 0;
           for i = 1:size(objL,2)
          %     i,size(objL(i).f,2)
          
             if size(objL(i).f,2) == 1
               n = n + 1;
               objR(n) = objL(i);
               %disp("orig")
               %objR(n).print
               index2(n) = index(i);
               continue;
             end
             [l, fmax, ind, lSing] = objL(i).d.maximum(objL(i).f);
             if lSing
                continue
             end
             %l
             if l
                % disp("l1")
               n = n + 1;
               objR(n) = functionNDomain([fmax],objL(i).d);
               %objR(n).print
               index2(n) = index(i,ind);
               continue
             end
             %disp('maxsplit')
             %objL(i).f(1).print
             %objL(i).f(2).print
             ineqs = objL(i).d.splitmax3 (objL(i).f(1),objL(i).f(2));
             ineqs1 = sym.empty ;          
             for k = 1: size(objL(i).d.ineqs,2)
               ineqs1(k) = objL(i).d.ineqs(k).f;
             end
             
             ineqs1(size(objL(i).d.ineqs,2)+1) = ineqs(1).f;
             d1 = region(ineqs1,objL(i).d.vars);
             d1 = d1.simplify;
             %d1.print
             if ~ isempty(d1)
             n = n + 1;
             %disp("l2")
             objR(n) = functionNDomain([objL(i).f(1)],d1);
             %objR(n).print
             
             index2(n) = index(i,1);
             end
             ineqs1(size(objL(i).d.ineqs,2)+1) = ineqs(2).f;
             d1 = region(ineqs1,objL(i).d.vars);
             d1 = d1.simplify;
             if ~ isempty(d1)
              %   d1.print
             n = n + 1;
             
             objR(n) = functionNDomain([objL(i).f(2)],d1);
             index2(n) = index(i,2);
             %disp("l3")
             %objR(n).print
             
             end
           end
           if n == 0
               objR = functionNDomain.empty();
              return
           end
           return
           
         end
          % will work only when entire R2 is covered
         function objL = mtimes (objL1, objL2)
             n = 0;
             objL=functionNDomain.empty();
             for i = 1:size(objL1,2)
               for j = 1:size(objL2,2)
                  % i,j
                 
                % if i == 1 & j == 4
                %disp('mtimes')
                 %    objL1(i).d.print
                 %    objL2(j).d.print
                %     rf.print
                % end
                rf = objL1(i).d + objL2(j).d;
                 if isempty(rf)
                   continue
                 end
              %   if n == 2
              % disp("out")
              % objL1(i).f(1).print
              % objL2(j).f(1).print
                  
              %   end
             % rf.print
                 rf = rf.simplifyUnboundedRegion ;
                 %rf.print
                 %rf = rf.simplifyOpenRegion ();
                 if isempty(rf)
                     disp("empty")
                   continue
                 end
               %  if n == 2
                  %   disp('times')
               %  rf.print
                 %end
                 n = n + 1;
                 objL(n) = functionNDomain([objL1(i).f(1), objL2(j).f(1)],rf);


               end
             end

         end

         function objR3 = maximumP(objL, lmerge) %, f, r2)
             % disp("in maximumP")
             % objL.printL
           % HISTORY: objR3 must start initialized. The `if n == 0 ... return` below assigns
           % objR (a different variable) and returns, leaving the OUTPUT unassigned, so MATLAB
           % raises "Output argument objR3 ... not assigned" instead of returning a well-defined
           % empty result -- the same shape as region.getEdges' and region.splitmax3' own bugs.
           % n == 0 simply means no piece survived, which is a legitimate answer, and it is
           % reachable from the biconjugate's max (plq.biconjugateF -> maxOfList).
           objR3 = functionNDomain.empty();
           n = 0;
           for i = 1:size(objL,2)
             %  i
             if size(objL(i).f,2) == 1
               n = n + 1;
               objR(n) = objL(i);
               continue;
             end
             %objL(i).f.printL
             [l, fmax, ind, lSing] = objL(i).d.maximum(objL(i).f);

             % HISTORY (tie-point fix): region.maxArray sets lSing when it
             % cannot find an interior probe point where f1~=f2 to decide
             % which of f1,f2 dominates -- but for a region whose only
             % vertices are curve (e.g. parabola) intersections, an affine
             % f1-f2 can be tied at every vertex while still genuinely
             % changing sign in the interior (e.g. two independent triangles'
             % conjugates overlapping in a lens bounded by two parabolas
             % meeting exactly where f1==f2). Discarding the pair here (the
             % old behaviour) silently dropped that whole overlap region from
             % the assembled partition, producing a real coverage gap at the
             % tie point (not a tolerance artifact) -- see cplqAdapterTest.m /
             % conjCPLQTest.m's exact-tie-point cases and
             % .claude/SESSION_HANDOFF.md. splitmax3 makes its own sign
             % decision directly from a region vertex and handles the tied
             % case correctly (falls back to a well-defined split rather than
             % guessing), so when lSing fires we must still fall through to
             % it instead of skipping -- mirrors what maxEqDom already does
             % (it captures lSing but never checks it, always falling
             % through to splitmax3 when l is false).
             if ~lSing && l
               n = n + 1;
               objR(n) = functionNDomain([fmax],objL(i).d);
               continue
             end
             ineqs = objL(i).d.splitmax3 (objL(i).f(1),objL(i).f(2));
             ineqs1 = sym.empty ;
             for k = 1: size(objL(i).d.ineqs,2)
               ineqs1(k) = objL(i).d.ineqs(k).f;
             end
             ineqs1(size(objL(i).d.ineqs,2)+1) = ineqs(1).f;
             d1 = region(ineqs1,objL(i).d.vars);

              d1 = d1.simplifyUnboundedRegion ;
             % HISTORY (tie-point fix follow-up): one of the two split halves
             % can simplify down to empty (the split boundary can coincide
             % with an existing edge of a degenerate/tangent domain) -- must
             % not add an empty-domain piece here, same class of gap as
             % mergeL's (see its own HISTORY comment); an empty-domain
             % functionNDomain surviving into mergeL crashes region.merge
             % (obj2.ineqs(i) on a 0x0 region has no elements to index).
             if ~isempty(d1)
               n = n + 1;
               objR(n) = functionNDomain([objL(i).f(1)],d1);
             end
             ineqs1(size(objL(i).d.ineqs,2)+1) = ineqs(2).f;
             d1 = region(ineqs1,objL(i).d.vars);
             d1 = d1.simplifyUnboundedRegion ;
             if ~isempty(d1)
               n = n + 1;
               objR(n) = functionNDomain([objL(i).f(2)],d1);
             end
 
             
           end
           
           if n == 0
              return
           end
           
           if ~lmerge
               objR3 = jSort(objR);
               return
            end

           % INSTRUMENTATION (CCA2_TRACE_MAXP). Step 3's cost is driven by the CELL COUNT, and
           % the count is set here: the split loop above makes cells and mergeL below is the
           % only thing that ever unmakes them. Printing in/afterSplit/afterMerge per fold is
           % what says whether a fold's growth is cells being created or merges being refused;
           % region.mergeTally says WHY they were refused.
           lTrace = ~isempty(getenv('CCA2_TRACE_MAXP'));
           if lTrace, tMerge = tic; end
           objR2 = mergeL(objR);
            objR3 = mergeL(objR2);
           if lTrace
             fprintf('[maxP] in=%d afterSplit=%d merge1=%d merge2=%d (%.0f s)\n', ...
                     size(objL,2), n, size(objR2,2), size(objR3,2), toc(tMerge));
           end
           % objR2 = mergeL(objR3);
            % disp("aft merge")
            % objR.printL
         end

         function [objL2,index,lCh] = maxEqDom(objL)  
           lCh = false;  
           ia(1) = 1;
           n = 0;
           for i = 1:size(objL,2)
              marked(i) = false;
           end
           ja = [];
           % ja has indices of all equal functions , ia by col no
           for i = 1:size(objL,2)
              if (marked(i))
                  ia(i+1) = n+1;
                  continue
              end
              
              for j = i+1:size(objL,2)
                  
                  if isAlways(objL(i).d == objL(j).d)
                      n = n+1;
                      ja(n) = j;
                      marked(j) =true;
                  end
              end
              ia(i+1) = n+1;
            end
            %ia
            %ja
            for i = 1:size(objL,2)
              marked(i) = false;
            end
            m = 0;
            for i = 1:size(objL,2)
               % i, marked(i), m
              if marked(i)
                continue;
              end
              if ia(i) == ia(i+1)
                m = m + 1;
                objL2(m) = objL(i);
                index(m) = i;
                marked(i) = true;
                continue;    
              end

              
              for j=ia(i):ia(i+1)-1
                [l, fmax, ind, lSing] = objL(i).d.maximum([objL(i).f, objL(ja(j)).f]);
                %l
                lCh = true;
                marked(i) = true;
                marked(ja(j)) = true;
               if l
                 m = m + 1;
                 objL2(m) = functionNDomain(fmax,objL(i).d);
                 if ind == 1
                   index(m) = i;
                 else
                   index(m) = ja(j);  
                 end
                 
                 
               else
                   
                 %ineqs = objL(i).d.splitmax2 (objL(i).f, objL(ja(j)).f);
                 ineqs = objL(i).d.splitmax3 (objL(i).f,objL(ja(j)).f);
                 ineqs1 = sym.empty;
                 for k = 1: size(objL(i).d.ineqs,2)
                    ineqs1(k) = objL(i).d.ineqs(k).f;
                 end
                            
                 ineqs1(size(objL(i).d.ineqs,2)+1) = ineqs(1).f;
                 d1 = region(ineqs1,objL(i).d.vars);
                 d1 = d1.simplify; 
                 
                 m = m + 1;
                 objL2(m) = functionNDomain(objL(i).f,d1);
                 index(m) = i;
                 
                 ineqs1(size(objL(i).d.ineqs,2)+1) = ineqs(2).f;
                 d1 = region(ineqs1,objL(i).d.vars);
                 d1 = d1.simplify; 
                 m = m + 1;
                 objL2(m) = functionNDomain(objL(ja(j)).f,d1);
                 index(m) = ja(j);
                 
               end  
                  
                 
              end
            end
            %disp('list size')
            %m 
          end
          
          function [objL2,index,lCh] = maxEqFun(objL)  
           lCh = false;  
           ia(1) = 1;
           n = 0;
           for i = 1:size(objL,2)
              marked(i) = false;
           end
           ja = [];
           % ja has indices of all equal functions , ia by col no
           for i = 1:size(objL,2)
              if (marked(i))
                  ia(i+1) = n+1;
                  continue
              end
              
              for j = i+1:size(objL,2)
                  
                  if isAlways(objL(i).f == objL(j).f)
                      n = n+1;
                      ja(n) = j;
                      marked(j) =true;
                  end
              end
              ia(i+1) = n+1;
            end
            for i = 1:size(objL,2)
              marked(i) = false;
            end
            m = 0;
            for i = 1:size(objL,2)
               % i, marked(i), m
              if marked(i)
                continue;
              end
              if ia(i) == ia(i+1)
                m = m + 1;
                objL2(m) = objL(i);
                index(m) = i;
                marked(i) = true;
                continue;    
              end

              
              eqFun = functionNDomain.empty;
              nF = 1;
              eqFun(nF) = objL(ia(i));
              for j=ia(i):ia(i+1)-1
                nF = nF+1;
                eqFun(nF) = objL(ja(j));
              end
              %disp('eqFun')
              %eqFun.printL
              [eqDom,indexD] = eqFun.*eqFun;
              [eqDom] = eqDom.unique ;
              
              %disp('eqDom1')
              %eqDom.printL
              eqFun = eqDom;
              [eqDom,indexD] = eqFun.*eqFun
              eqDom = eqDom.unique
              %disp('eqDom2')
              %eqDom.printL
              for j = 1:size(eqDom,2)
                m = m + 1;
                objL2(m) = eqDom(j);
                index(m) = i;
              end
                  
                 
              
            end
            for i = 1:size(objL2,2)
                size(objL2(i).f)
                if size(objL2(i).f,2) == 2
                    %disp('checking')
                    %objL2(i).f.printL
                    objL2(i).f([2]) = [];
                    %disp('checking2')
                    %objL2(i).f.printL
                end
                %disp('checking3')
                %objL2(i).f

            end
            %objL2.printL;
            %disp('list size')
            %m 
          end
         
          
          function objL2 = addEq(objL)
          % HISTORY: objL2 must start initialized -- if objL is empty, or no
          % entry ever satisfies both isempty(r) and r.nv==size(r.ineqs,2)
          % below, the loop body never executes/never assigns objL2, and
          % MATLAB errors with "output argument not assigned" rather than
          % returning a well-defined empty result (hit via
          % plq.biconjugateF -> bc.addEq when bc, the merged
          % conjugateOfPiecePoly result, comes out empty for a given
          % domain -- testMaxMultiRegion/testPCE2's domain is one such case,
          % a separate, still-open gap in the biconjugate pipeline upstream
          % of this function -- see testPCE2's own comment).
          objL2 = functionNDomain.empty();
          ia(1) = 1;
          n = 0;
          for i = 1:size(objL,2)
              marked(i) = false;
          end
          ja = [];
          % ja has indices of all equal functions , ia by col no
          for i = 1:size(objL,2)
              if (marked(i))
                  ia(i+1) = n+1;
                  continue
              end
              
              for j = i+1:size(objL,2)
                  % i,j
                  % simplifyFraction(objL(i).f.f - objL(j).f.f)
                  % isAlways(objL(i).f.f == objL(j).f.f)
                  % isAlways(simplifyFraction(objL(i).f.f - objL(j).f.f) == 0)
                  %if isAlways(objL(i).f.f == objL(j).f.f)
                  if isAlways(simplifyFraction(objL(i).f.f - objL(j).f.f) == 0)
                      n = n+1;
                      ja(n) = j;
                      marked(j) =true;
                  end
              end
              ia(i+1) = n+1;
          end
      %     ia
      %     ja
          for i = 1:size(objL,2)
              marked(i) = false;
          end
          m = 0;
          for i = 1:size(objL,2)
            if marked(i)
                continue;
            end
            r = objL(i).d;
            % m = m + 1;
            % objL2(m) = objL(i);
            marked(i) = true;
            for j=ia(i):ia(i+1)-1
                %m = m + 1;
                %objL2(m) = objL(ja(j));
                if ~isempty(r)
                r = r + objL(ja(j)).d;
           %     r.print
                r = r.simplifyUnboundedRegion;
           %     ja(j)
            %    r.print
                end
                marked(ja(j)) = true;
            end
            if ~isempty(r)
                r =  r.removeInfV;
                if r.nv == size(r.ineqs,2)
            m = m + 1;
            objL2(m) = functionNDomain(objL(i).f, r);
                end
            end
           end
             
          end   
         


         function objL2 = jSort(objL)  
          ia(1) = 1;
          n = 0;
          for i = 1:size(objL,2)
              marked(i) = false;
          end
          ja = [];
          % ja has indices of all equal functions , ia by col no
          for i = 1:size(objL,2)
              if (marked(i))
                  ia(i+1) = n+1;
                  continue
              end
              
              for j = i+1:size(objL,2)
                  
                  if isAlways(objL(i).f.f == objL(j).f.f)
                      n = n+1;
                      ja(n) = j;
                      marked(j) =true;
                  end
              end
              ia(i+1) = n+1;
          end
          %ia
          %ja
          for i = 1:size(objL,2)
              marked(i) = false;
          end
          m = 0;
          for i = 1:size(objL,2)
            if marked(i)
                continue;
            end
            m = m + 1;
            objL2(m) = objL(i);
            marked(i) = true;
            for j=ia(i):ia(i+1)-1
                m = m + 1;
                objL2(m) = objL(ja(j));
                marked(ja(j)) = true;
            end
           end
             
          end
          
          function [objL2,index] = unique(objL)  
          ia(1) = 1;
          n = 0;
          for i = 1:size(objL,2)
              marked(i) = false;
          end
          ja = [];
          % ja has indices of all equal functions , ia by col no
          for i = 1:size(objL,2)
              if (marked(i))
                  ia(i+1) = n+1;
                  continue
              end
              
              for j = i+1:size(objL,2)
                  
                  if isAlways(objL(i).f(1).f == objL(j).f(1).f) & (objL(i).d == objL(j).d)
                      n = n+1;
                      ja(n) = j;
                      marked(j) =true;
                  end
              end
              ia(i+1) = n+1;
          end
          for i = 1:size(objL,2)
              marked(i) = false;
          end
          m = 0;
          for i = 1:size(objL,2)
            %i = ia(j)  
            if marked(i)
                continue;
            end
            m = m + 1;
            objL2(m) = objL(i);
            
            marked(i) = true;
             for j=ia(i):ia(i+1)-1
            %     m = m + 1;
            %     objL2(m) = objL(ja(j));
                 marked(ja(j)) = true;
             end
           end
             
          end
         
          function [objL2,index] = mergeL(objL)
          ia(1) = 1;
          n = 0;
          for i = 1:size(objL,2)
              marked(i) = false;
          end
          ja = [];
          % ja has indices of all equal functions , ia by col no
          for i = 1:size(objL,2)
              if (marked(i))
                  ia(i+1) = n+1;
                  continue
              end
              
              for j = i+1:size(objL,2)
                  
                  if isAlways(simplifyFraction(objL(i).f.f - objL(j).f.f)==0)
                      n = n+1;
                      ja(n) = j;
                      marked(j) =true;
                  end
              end
              ia(i+1) = n+1;
          end
           
          m = 0;
          for i = 1:size(objL,2)
              marked(i) = false;
          end
          for i = 1:size(objL,2)
              
            if  marked(i)
                continue
            end
            %disp('first')
            if (ia(i) == ia(i+1))
                if marked(i)
                   continue
                end

                marked(i) = true;
                % HISTORY (tie-point fix follow-up): objL(i).d can already be
                % empty here (e.g. a maximumP split half from an earlier,
                % unrelated merge round that simplified to empty) even when
                % there's no same-valued partner to merge with -- don't copy
                % it forward uncheck (crashes region.merge/removeTangent
                % later when reused as a same-valued partner elsewhere; see
                % maximumP/mergeL's other HISTORY comments for the class of
                % bug this belongs to).
                if ~isempty(objL(i).d)
                  m = m + 1;
                  objL2(m) = objL(i);
                  index(m) = i;
                end
            else
              
               r = objL(i).d;
               lmerge = true;
               while lmerge
                 lmerge = false;
                 for j=ia(i):ia(i+1)-1
                   if marked(ja(j))
                       continue
                   end
                   
                   [l,r] = r.merge (objL(ja(j)).d);
                   if l
                     marked(ja(j)) = true;
                     lmerge = true;
                   end
                 end
               end
               marked(i) = true;
               % HISTORY (tie-point fix follow-up): r (seeded from objL(i).d)
               % can already be empty here -- e.g. objL(i).d itself was an
               % empty domain surviving an earlier unguarded copy-through, in
               % which case region.merge's own empty-operand guard leaves r
               % unchanged (still empty) throughout the accumulation loop
               % above. r.nv below (and removeTangent further down) both
               % assume a scalar, non-empty region; skip this piece entirely
               % rather than erroring, matching the isempty guards already
               % used by region.merge and functionNDomain.mtimes elsewhere.
               % This path was unreachable before maximumP stopped discarding
               % lSing pairs (see maximumP's own HISTORY comment).
               if ~isempty(r)
               % Removing inf from vertices - to be removed later
               [nP, px, py] = r.finiteVertices;
             %  disp('b4 simp')
               r = r.simplifyUnboundedRegion;
             %  disp('aft simp')
               if ~isempty(r)
                 m = m + 1;
                 r = r.removeTangent (nP, px,py);
                 objL2(m) = functionNDomain([objL(i).f],r);
                 index(m) = i;
               end
               end
               for j=ia(i):ia(i+1)-1
                 if marked(ja(j))
                   continue
                 end
                 r = objL(ja(j)).d ;
                 lmerge = true;
                 while lmerge
                   lmerge = false;
                   for k=j+1:ia(i+1)-1
                     if marked(ja(k))
                       continue
                     end
                     [l,r] = r.merge (objL(ja(k)).d);
                     if l
                       marked(ja(k)) = true;
                       lmerge = true;
                     end
                   end
                 end
                % (nP,px,py) must come from THIS r, read before simplifyUnboundedRegion --
                % they are removeTangent's candidate tangency points. This loop used to reuse
                % whatever the FIRST accumulated region left in those variables above, which is
                % a different region, and is undefined outright whenever that first region was
                % empty and skipped ('Unrecognized function or variable nP' -- f=xy over
                % conv{(0,0),(3,3),(1,2)}).
                [nP, px, py] = r.finiteVertices;
                r = r.simplifyUnboundedRegion;
                marked(i) = true;
                marked(ja(j)) = true;
                % HISTORY (tie-point fix follow-up): see the matching guard
                % above -- an accumulated merge can still simplify to empty.
                if ~isempty(r)
                  m = m + 1;
                  r = r.removeTangent (nP, px,py);
                  objL2(m) = functionNDomain([objL(i).f],r);
                  index(m) = i;
                end

                 
               end

            end
        end
      end
     end

     methods % conjugate

        function [pc,ia] = conjugateOfPiecePoly (obj)
            disp("in conjugateOfPiecePoly")

            % EXPLICIT EDGE LIST PER PIECE, empty where the slot conventions still hold. Built in
            % the loop below and consumed in the second one -- both are in this function, so no
            % class or cross-file contract changes to carry it. See functionNDomain.edgeIndexList.
            eIdxAll = cell(1, size(obj,2));

            for i=1:size(obj,2)
                d = obj(i).d;

                % (0) DROP PROVABLY REDUNDANT CONSTRAINTS FIRST, from this local copy.
                % This routine is edge-indexed throughout (see the ADAPTER note below), and a
                % redundant constraint bounds no edge -- so by the adapter's own principle it has
                % no place here. Doing it up front, rather than only when the edge numbering
                % happens to collide, is what the biconjugate needs: the pieces of f* are affine
                % functions on cones, and such a cone routinely arrives with a redundant facet
                % (on the unit box with f = x*y, f*'s piece {s2>=1, s1>=1, s1<=s2} has s2>=1
                % implied). With it present the edge indexing mis-assigns and the conjugate comes
                % out on {x>=1, x+y<=2} instead of the true {y<=1, x+y<=2}.
                %
                % LOCAL COPY ONLY. The assembly path keeps these constraints -- preserving them
                % is what took Step 3 from 57 wrong grid points to 0 -- and nothing here is
                % pushed back into region.redundantSubset's own policy. redundantSubset is used
                % purely as an LP-certified test, and it never deletes every copy of a duplicated
                % constraint (regionTest pins that).
                delR = d.redundantSubset(1:size(d.ineqs,2));
                if ~isempty(delR)
                    keepR = true(1, size(d.ineqs,2));
                    keepR(delR) = false;
                    if any(keepR)
                        d.ineqs = d.ineqs(keepR);
                        d = d.getVertices;
                    end
                end


                d = d.removeTangent(d.nv,d.vx,d.vy);
                % removeTangent can delete every constraint and hand back region.empty. An empty
                % region has no nv/ineqs to index -- obj.nv on a 0x0 region array is a
                % comma-separated list with 0 values, which is exactly what removeInfV errored on
                % ('a dot indexing expression produced ... 0 values') -- and a piece whose domain
                % is empty contributes nothing to the conjugate anyway, so carry it through
                % untouched. Same class of guard as region.plus/merge/finiteVertices and this
                % file's own mergeL/maximumP; reachable here now that region.redundantSubset
                % keeps constraints the old delete-anything-without-a-finite-vertex rule removed.
                if isempty(d)
                    obj(i).d = d;
                    continue
                end
                % IS THIS REGION BOUNDED? Decide it HERE, because removeInfV is about to delete
                % the only evidence. getVertices marks an unbounded region by appending box-clip
                % vertices at +-intmax; removeInfV strips them, after which every remaining vertex
                % is finite and a bounded and an unbounded region are indistinguishable by their
                % vertex lists. resolveLensEdges needs the answer, and getting it wrong is not
                % harmless: it drops constraints that bound no EDGE, which on a genuinely
                % unbounded region are load-bearing. Measured -- reading it after removeInfV made
                % f** of x*y over the two-face square exact at (0.9,0.6) and (0.6,0.6) and +inf at
                % (0.25,0.25) and (0.1,0.1), trading one hole in the domain for another.
                wasBounded = true;
                for vI = 1:d.nv
                    if isAlways(abs(d.vx(vI)) == intmax) || isAlways(abs(d.vy(vI)) == intmax)
                        wasBounded = false; break
                    end
                end
                d = d.removeInfV;
                %d = d.poly2orderUnbounded;
               
                if d.nv == size(d.ineqs,2)
                    d = d.poly2order;
                else
                    d = d.poly2orderUnbounded;
                end
                edgeNo = d.getEdgeNosInf(d.vars);
                % This scatter re-indexes the constraints BY EDGE NUMBER, which is the contract
                % the endNv loop further down relies on: for a bounded domain edge j is ineqs(j),
                % for an unbounded one edge j is ineqs(j+1) (getEdgeNosInf's own `add`, which is
                % why slot 1 is deliberately left empty there). Only slots 1..nv are ever read as
                % edges, so the scatter growing d.ineqs past that is harmless by design.
                %
                % Do NOT replace the scatter with a sort. Tried, and it broke a third,
                % previously-passing test (testMaxMultiRegion/testFractional): edgeNo values can
                % exceed numel(ineqs) (hence the deliberate growth above) and can repeat (hence
                % last-write-wins), and sorting reproduces neither.
                %
                % ADAPTER, in two parts, for constraints that bound no edge of d. THIS ROUTINE IS
                % EDGE-INDEXED THROUGHOUT -- this scatter, the isQuad chord rewrite below,
                % getNormalConeVertex, getSubdiffVertexT1/T2/T2Q all address d.ineqs by edge --
                % so a constraint that bounds no edge has no place here at all. Both parts drop
                % from this function's LOCAL copy of d only: the conjugate/assembly path, where
                % these constraints are load-bearing (SUPPORT_MATRIX.md section 1.2), keeps them.
                % Do NOT push either back into region.redundantSubset -- the constraints it now
                % preserves are what took the Step 3 assembly from 57 wrong grid points to 0.
                %
                % (1) NO vertex on d: getEdgeNosInf reports 0, and feeding a 0 to the scatter
                % errors ('Array indices must be positive integers') -- it broke
                % testMaxMultiRegion/testMax and testcPLQ/testRectBiconj. Parking such a
                % constraint ABOVE the real slots was tried and is worse than useless: the isQuad
                % branch then derives a chord for it from d0.vx(1),d0.vx(2), vertices that have
                % nothing to do with it, and the comparison is undecidable.
                %
                % (2) EXACTLY ONE vertex, COLLIDING with a constraint that has two. Here
                % getEdgeNosInf reports not 0 but that vertex's own index -- precisely the slot
                % the real edge leaving that vertex claims -- and because the scatter is
                % last-write-wins the intruder EVICTS a genuine edge, while the evicted
                % constraint's original slot is left holding a stale copy of its old occupant.
                % testcPLQ/testRectBiconj piece 23 arrived as the triangle
                %     (9*s2)/5 - s1 + 5,   -s1 - 7*s2 - 4,   s1 + 2*s2 - 4
                % plus a quadratic active only at the single vertex (139/44,-45/44). edgeNo came
                % out [3 1 1 2] -- the quadratic and -s1-7*s2-4 both claiming slot 1 -- and the
                % scatter returned [quad, s1+2*s2-4, (9*s2)/5-s1+5, s1+2*s2-4]: -s1-7*s2-4 GONE,
                % s1+2*s2-4 present twice. The isQuad branch then chorded d0.vx(1) to d0.vx(2),
                % two vertices the quadratic does not join, so solve() returned a COMPLEX
                % conjugate pair for the third point and symbolicFunction.gtd -- a bare
                % `if (obj1.f>obj2)` -- cannot take an undecidable sym ('Conversion to logical
                % from sym is not possible').
                %
                % The test is the COLLISION, not the vertex count on its own, and that is not a
                % detail: a RAY edge of an unbounded region is also active at exactly one finite
                % vertex and is load-bearing, but getEdgeNosInf reserves slot 1 for it (its `add`)
                % so it never collides. Nor can boundedness be used to tell the two apart -- piece
                % 24 of this same test is unbounded (recession direction (2,-1)) yet carries no
                % vertex at infinity for removeInfV to find. Keying on the collision also bounds
                % the blast radius: a piece whose edgeNo has no repeats is untouched.
                nOn = zeros(size(edgeNo));
                for k = 1:size(d.ineqs,2)
                    nOn(k) = d.vertexOfEdge(k);
                end
                hasEdge = edgeNo > 0;
                keepE = hasEdge;
                for k = 1:numel(edgeNo)
                    if ~hasEdge(k)
                        continue
                    end
                    % Strictly-better rival only, so no two constraints can evict each other and
                    % a tie leaves the slot exactly as it is today.
                    if any(hasEdge & edgeNo == edgeNo(k) & nOn > nOn(k))
                        keepE(k) = false;
                    end
                end
                % (3) A TIE. The rule above drops only a STRICTLY worse rival, so two colliding
                % constraints with the SAME vertex count both survive and the scatter -- being
                % last-write-wins -- silently destroys one of them anyway. That is not a corner
                % case: it is what breaks the BICONJUGATE. The pieces of f* are affine functions
                % on cones, and such a cone routinely arrives carrying a redundant constraint.
                % On the unit box with f = x*y, f*'s piece {s1<=s2, s2<=0, s1<=0} (in which
                % s1<=0 is implied by the other two) gives edgeNo = [1 2 2], and the scatter
                % turns [s1-s2, s2, s1] into [s1-s2, s1, s1]: s2 GONE, s1 twice. The conjugate
                % of that piece then comes out on {x+y>=0, y<=0} instead of the true polar cone
                % {x>=0, x+y>=0}, and the max over the pieces collapses to nothing.
                %
                % Resolve the tie by LP-CERTIFIED REDUNDANCY: among the colliding constraints,
                % drop one that region.redundantSubset proves adds nothing to the feasible set.
                % This uses redundantSubset as a TEST, on this function's own local copy of d --
                % it does not change what redundantSubset preserves for the assembly path, which
                % is what the note above forbids. At most one constraint is dropped per slot, so
                % a slot can never be emptied, and a piece whose edgeNo has no repeats is still
                % untouched.
                still = find(keepE);
                slots = edgeNo(still);
                for sIdx = unique(slots(:))'
                    grp = still(slots == sIdx);
                    if numel(grp) < 2
                        continue
                    end
                    del = d.redundantSubset(grp);
                    if ~isempty(del)
                        keepE(del(1)) = false;   % one only: dropping it may un-redundant the rest
                    end
                end

                % (4) TWO GENUINE EDGES BETWEEN THE SAME TWO VERTICES. The three rules above all
                % resolve a collision by DROPPING a constraint, which is right when the collision
                % comes from a redundant or edge-less one. It is wrong when both constraints bound
                % a real edge -- and that is exactly a LENS: one arc and its chord, or two arcs,
                % joining the SAME pair of vertices. getEdgeNosInf numbers an edge by one of its
                % endpoint VERTICES, so two edges on the same pair are indistinguishable to it and
                % get the same number; the scatter below is last-write-wins, so one of them is then
                % destroyed and the other stored twice.
                %
                % MEASURED on the half-lens {(s1+s2)^2 <= 4*s2, s2 <= s1} carrying s1 -- which is
                % what f* of x*y over the unit square is made of: nv = 2, both constraints have
                % both vertices on them, both get edgeNo 1, and the conjugate comes out as TWO
                % IDENTICAL cells built from whichever constraint came last. Fed the arc first it
                % is the conjugate of the CHORD, finite only on a strip.
                %
                % So spread them instead of dropping one: both are edges, they simply need
                % distinct numbers.
                [edgeNo, keepE] = functionNDomain.spreadCollidingEdges(d, edgeNo, keepE, nOn, wasBounded);

                % (5) STILL COLLIDING, and both are genuine edges. spreadCollidingEdges resolves
                % a collision when some edge number is free; on the pipeline's own lens every
                % number it could move to is held by another constraint, so it leaves things
                % alone and the last-write-wins scatter below then DESTROYS one of the two edges.
                %
                % There is nothing left to rearrange at that point -- the slot conventions cannot
                % represent two edges on one vertex pair at all -- so stop using slots for this
                % piece and derive the edge list from the geometry instead. The scatter is what
                % maintains the slot correspondence, so it is skipped exactly when the list
                % replaces it; the constraints themselves are all kept, which is what the vertex
                % cones and every feasibility test below still need.
                %
                % SCOPE. edgesStillCollide is the lens signature and nothing else reaches it, and
                % edgeIndexList prefers today's slot wherever it is geometrically valid -- so a
                % piece that works today keeps the same indices whether it takes this branch or
                % not.
                eList = [];
                if functionNDomain.edgesStillCollide(edgeNo, keepE, nOn)
                    dTry = d;
                    dTry.ineqs = dTry.ineqs(keepE);
                    [eTry, okL] = functionNDomain.edgeIndexList(dTry, wasBounded);
                    if okL
                        d = dTry;
                        eList = eTry;
                    end
                end
                if isempty(eList)
                    d.ineqs = d.ineqs(keepE);
                    edgeNo  = edgeNo(keepE);
                    d.ineqs(edgeNo) = d.ineqs;

                    % (6) A BOUNDED region whose constraint COUNT does not match its vertex count.
                    % Both slot conventions are chosen by that count, so such a region is read with
                    % the UNBOUNDED one -- endNv = nv-1, so it is built one edge cell SHORT, and
                    % every edge is read one slot along. A bounded piece's conjugate is finite
                    % everywhere, so a missing edge cell is a hole where the answer must be finite.
                    %
                    % MEASURED on piece 9 of f* for x*y over the parallelogram
                    % conv{(0,0),(2,0),(2.5,1),(0.5,1)}: a bounded region with 3 vertices and 3
                    % genuine edges, carrying a fourth constraint -- a conic touching it at ONE
                    % vertex. 4 constraints, 3 vertices, so the count says "unbounded"; the
                    % conjugate came back with 4 cells instead of 6, wrong or uncovered at 6 of 10
                    % probe points against a brute-force sup, and the biconjugate over the whole
                    % parallelogram then collapsed to nothing (QuaParCPLQ:conj:emptyResult).
                    %
                    % SCOPE, and why this cannot disturb what works. It runs only when the count
                    % DISAGREES with the vertex count, so every piece the conventions already
                    % describe correctly is untouched; edgeIndexList refuses an unbounded region
                    % outright and returns ok = false whenever the geometry does not settle the
                    % whole list, in which case the old reading stands; and the readers below give
                    % the count-matching branches precedence over the list.
                    if size(d.ineqs,2) ~= d.nv
                        [eTry, okL] = functionNDomain.edgeIndexList(d, wasBounded);
                        if okL
                            eList = eTry;
                        end
                    end
                end
                eIdxAll{i} = eList;

                obj(i).d = d;

            end
         
             % redesign so u can have a list of piecewise functions
             ia(1) = 1;
             pc = functionNDomain.empty();
             x=sym('x');
             y=sym('y');
             for i = 1:size(obj,2)
               f = obj(i).f;
               vars = f.getVars;
               d = obj(i).d;
               d0 = d;
               % A piece whose domain came out empty (see the removeTangent guard in the first
               % loop above) has nothing here to read: d0.ineqs on a 0x0 region array is a
               % comma-separated list with 0 values, which is what the isQuad test below errored
               % on. It contributes no conjugate piece -- but ia must stay a running index into
               % pc with one entry per INPUT piece, so record that this one added none
               % (ia(i+1)==ia(i), the same convention maxEqDom uses for a skipped entry) rather
               % than leaving a hole for callers to trip over.
               if isempty(d)
                   ia(i+1) = size(pc,2)+1;
                   continue
               end
               if f.isQuad

                   for j = 1:size(d0.ineqs,2)
                       if (d0.ineqs(j).isQuad)
                           % CHORDING THE FIRST TWO VERTICES, whatever the conic actually joins,
                           % is deliberate here for want of anything better -- and it is NOT the
                           % cause of the two remaining wrong values on the parallelogram's piece
                           % 9. Both alternatives were measured 2026-08-16 and are recorded in
                           % DECISIONS.md: chording the vertices the conic actually touches
                           % (region.vertexOfEdge) makes that piece WORSE, 2 wrong of 10 -> 3, and
                           % skipping the rewrite altogether changes nothing there at all. The
                           % remaining defect is functionNDomain.getInterior on a SINGULAR
                           % quadratic; do not attack this line again for it.
                           m0 = (d0.vy(1)-d0.vy(2))/(d0.vx(1)-d0.vx(2));
                           c = d0.vy(1) - m0*d0.vx(1);
                           nineq = symbolicFunction(vars(2)-m0*vars(1)-c);
                           mx = (d0.vx(1)+d0.vx(2))/2;
                           d1 = d0.ineqs(j).subsF(vars(1),mx);
                           my2 = solve(d1.f,vars(2));
                           %d0.ptFeasible(vars,[mx,my2(1)])
                           if (d0.ptFeasible(vars,[mx,my2(1)]))
                               my = my2(1);
                           else
                               my = my2(2);
                           end
                           if isAlways(nineq.subsF(vars,[mx,my])>0)
                               nineq = -nineq;
                           end
                           d0.ineqs(j) = nineq;
                       end
                   end
               end
               % WHICH NORMAL-CONE ROUTINE. getNormalConeVertex builds the cone at a vertex from
               % the SLOPE OF THE CHORD to the next vertex (region.slope); getNormalConeVertexQ
               % builds it from the slope of the CONSTRAINT ITSELF at that vertex
               % (region.slopeIneq), i.e. the curve's own tangent. For a region with a parabolic
               % edge only the second is right -- the chord is not the boundary there.
               %
               % The test used to be `size(d.ineqs,2) == d.nv`, a constraint COUNT standing in for
               % "is this region unbounded", which is not the question. It fails exactly on a
               % bounded region with a curved edge: the half-lens has nv = 2 and, once its two
               % edges are given distinct numbers, exactly 2 constraints -- so the count sends a
               % CURVED region to the polyhedral routine, its vertex cones come out wrong, and
               % every vertex cell is dropped as empty. Ask about the constraints instead.
               hasQuadCon = false;
               for kq = 1:size(d.ineqs,2)
                   if d.ineqs(kq).isQuad, hasQuadCon = true; break, end
               end
               % PRECEDENCE. The count-matching branches come FIRST, so a region the slot
               % conventions already describe correctly keeps exactly the reading it has; the
               % explicit list is consulted only where the count disagrees with the vertex count,
               % which is the case those conventions cannot express.
               eIdx = eIdxAll{i};
               if size(d.ineqs,2) == d.nv && ~hasQuadCon
                   NCV = d.getNormalConeVertex(x, y);
               elseif ~isempty(eIdx)
                   NCV = d0.getNormalConeVertexQ(x, y, eIdx);
               else
                   NCV = d0.getNormalConeVertexQ(x, y);
               end
               [subdV,undV] =  obj(i).getSubdiffVertexT1 (NCV, [x,y]);
               exprs = obj(i).conjugateExprVerticesT1 ([x,y], undV );
               
               for j = 1:size(exprs,2)
                   r = region(subdV(j,:),[x,y]);
                   r = r.simplifyUnboundedRegion;
                   if ~isempty(r)
                     pc = [pc,functionNDomain([symbolicFunction(exprs(j))],r)];
                   end
               end
               if obj(i).d.nv > 1
                   % ONE STATEMENT OF THE CONVENTION, not two. edgeIneq(j) is the constraint
                   % bounding edge j, and the loop below reads it from there instead of
                   % re-deriving the same rule from endNv -- which is how the two could disagree
                   % once a third branch existed.
                   if size(d0.ineqs,2) == d.nv
                       NCE = d0.getNormalConeEdgeQ3(x, y);
                       [subdE,unR] = obj(i).getSubdiffVertexT2 (NCE, [x,y]);
                       edgeIneq = 1:obj(i).d.nv;
                   elseif ~isempty(eIdx)
                       NCE = d0.getNormalConeEdgeQE(x, y, eIdx);
                       [subdE,unR] = obj(i).getSubdiffVertexT2Q (NCE, [x,y]);
                       edgeIneq = eIdx;
                   else
                       NCE = d0.getNormalConeEdgeQ(x, y);
                       [subdE,unR] = obj(i).getSubdiffVertexT2Q (NCE, [x,y]);
                       % d0 is a copy of d = obj(i).d, so reaching this branch already means
                       % size(obj(i).d.ineqs,2) ~= nv -- the old test here was dead.
                       edgeIneq = 2:obj(i).d.nv;
                   end
                   endNv = numel(edgeIneq);
                   for j = 1:endNv % fix this  obj(i).d.nv-1 ?
                       ineq = subs(obj(i).d.ineqs(edgeIneq(j)).f,obj(i).d.vars,[x,y]);

                       f0 = subs(obj(i).f.f,obj(i).d.vars,[x,y]);
                       [expr] = simplifyFraction(conjugateExpr(ineq,f0,x,y));
                       ineq1 = subdE(j,:);
                   if obj(i).f.isQuad
                       edgeInt = obj(i).getInterior(x,y,obj(i).d.ineqs(edgeIneq(j))) ;
                       s = solve(ineq1,[x,y]);
                       px = s.x;
                       py = s.y;
                       % temp fix
                       if isAlways (edgeInt == y/4 - x/2 + 1) 
                           edgeInt = -edgeInt;
                       end
                       if isAlways(subs(edgeInt,[x,y],[px,py])>0   )
                           edgeInt = -edgeInt;
                       end
                       edgeInt = simplifyFraction(edgeInt);
                       ineq1 = [ineq1,edgeInt];
                       r = region(ineq1, [x,y]);
                    else
                       edgeInt = obj(i).getInterior(x,y,obj(i).d.ineqs(edgeIneq(j))) ;
                       ineq1 = [ineq1,edgeInt];
                       r = region(ineq1, [x,y]);
                   end
                   r = r.simplifyUnboundedRegion;
                   if ~isempty(r)
                     pc = [pc,functionNDomain([symbolicFunction(expr)],r)];
                   end

                end
               end
              ia(i+1) = size(pc,2)+1;
             end
             
             return
             

        end

        function ineq = singularEdgeCut(obj, x, y, edgeIneq)
        % The multiplier cut described in getInterior's header, or [] when this is not the case
        % it is for. Deliberately narrow: it fires ONLY for a quadratic f whose Hessian is
        % singular and an AFFINE edge, which is exactly the case the elimination cannot do.
        % Everything else keeps the behaviour it had.
            ineq = [];
            vars = obj.d.vars;
            try
                Q = double(hessian(obj.f.f, vars));
            catch
                return
            end
            if any(~isfinite(Q(:))), return, end
            Q = (Q + Q')/2;
            if abs(det(Q)) > 1.0d-9 * max(1, max(abs(Q(:))))^2
                return                      % nonsingular: the elimination is already right
            end
            e = edgeIneq;
            if ~e.isLinear, return, end
            cf = e.getLinearCoeffs(vars);   % [a1 a2 a0] for a1*v1 + a2*v2 + a0 <= 0
            a1 = cf(1); a2 = cf(2); a0 = cf(3);
            if isAlways(a1 == 0) && isAlways(a2 == 0), return, end
            n = [a1; a2];
            d = [-a2; a1];
            Qs = sym(Q);
            dQd = simplify(d.' * Qs * d);
            if isAlways(dQd == 0)
                return                      % f is affine along this edge: no interior cut
            end
            % A point on the edge line, picked from whichever coefficient is nonzero.
            if isAlways(a1 == 0)
                x0 = [sym(0); -a0/a2];
            else
                x0 = [-a0/a1; sym(0)];
            end
            L = sym(zeros(2,1));
            L(1) = subs(obj.f.dfdx(vars(1)).f, vars, [0 0]);
            L(2) = subs(obj.f.dfdx(vars(2)).f, vars, [0 0]);
            sv = [x; y];
            tstar = ((sv - (Qs*x0 + L)).' * d) / dQd;
            xstar = x0 + d*tstar;
            mu    = n.' * (sv - (Qs*xstar + L));
            ineq  = simplifyFraction(expand(-mu));   % mu >= 0
        end

        function ineq = getInterior(obj,x,y,edgeIneq)
        % The cut that separates an EDGE cell of the conjugate from the interior.
        %
        % SINGULAR QUADRATIC (4th argument, the edge's own primal constraint). The elimination
        % below reads the interior side off `x = d1f, y = d2f`. When f is a NONSINGULAR
        % quadratic the gradient map is a bijection, the elimination is vacuous, and the answer
        % is the identically-zero expression `region` then drops -- which is right, because the
        % subdifferential cone already bounds the cell.
        %
        % When f is SINGULAR it is neither. The gradient map collapses the plane onto a LINE,
        % `solve` returns that line, and the elimination hands back the map's IMAGE -- a curve
        % that separates nothing, added to the cell as if it did. Measured on piece 9 of f* for
        % x*y over conv{(0,0),(2,0),(2.5,1),(0.5,1)}, whose f = s1*s2/2 + s1^2/8 + s2^2/2 has
        % det(hessian) = 0: the chord's edge cell and the arc's then claim the SAME region, the
        % chord's is checked first, and it is wrong there -- 4 of 10 probe points wrong or
        % double-covered against a brute-force sup.
        %
        % WHAT THE CUT ACTUALLY IS, and it needs no inversion. The edge cell holds the s whose
        % maximiser sits in the relative interior of that edge, and the boundary with the
        % interior cell is where the KKT multiplier vanishes:
        %       x* = x0 + t* d,   t* = <s - grad f(x0), d> / (d'Qd),   mu = <n, s - grad f(x*)>,
        % with n the edge's outward normal and d its direction. mu is AFFINE in s, so the cell
        % stays representable; `d'Qd = 0` means f is affine along that edge, its max sits at an
        % endpoint, and the vertex cones already cover it -- so no cut, which is what the
        % elimination happens to give too. This is the same construction conjConvexOverPiece
        % uses for its edge cells, and it never inverts Q.
                if nargin >= 4 && obj.f.isQuad
                    ineqS = obj.singularEdgeCut(x, y, edgeIneq);
                    if ~isempty(ineqS)
                        ineq = ineqS;
                        return
                    end
                end
                g(1) = obj.f.dfdx(obj.d.vars(1));
                g(2) = obj.f.dfdx(obj.d.vars(2));
                eq1 = x - g(1).f;
                eq2 = y - g(2).f;
                %if obj.f.isQuad
                    s12 = solve([eq1,eq2],obj.d.vars);
                    % HISTORY: s12 is the struct solve() returns for a
                    % multi-variable system (fields named after each dual
                    % variable). useFallback covers two cases the same way:
                    % (1) no solution at all (the original isempty check),
                    % and (2) a solution exists but sits exactly on a
                    % singularity of eq1/eq2's own denominator (e.g. a
                    % fractional-power/sqrt term whose argument vanishes
                    % there -- testFractional's domain has such a point),
                    % making the direct substitution below divide by zero.
                    % Both fall back to the same coeffs-based linear
                    % elimination, which needs no point evaluation at all.
                    useFallback = isempty(s12) || isempty(s12.s_1);
                    if ~useFallback
                        try
                            % subs needs the field VALUES as an array
                            % matching obj.d.vars's order, not the struct
                            % itself (subs errored: "Substitution expression
                            % X must be a symbolic, cell, or numeric array").
                            ineq = subs(eq2,obj.d.vars,[s12.s_1,s12.s_2]);
                        catch
                            useFallback = true;
                        end
                    end
                    if useFallback
                       c1 = coeffs(eq1,obj.d.vars(1));
                       c2 = coeffs(eq2,obj.d.vars(1));

                       % ELIMINATE vars(1) BETWEEN eq1 AND eq2. coeffs returns [const, slope]
                       % when the expression depends on vars(1) and just [const] when it does
                       % not, so the two sizes must be tested SEPARATELY. This tested only
                       % size(c1,2) and then indexed c2(2), raising MATLAB:badsubscript for any
                       % f with an s1^2 term but no s1*s2 term -- eq1 = x - s1 has two
                       % coefficients, eq2 = y - <no s1> has one. Same "index used outside the
                       % guard that established it" shape as region.probeVertexIndex's case.
                       %
                       % When only ONE of the two depends on vars(1) there is nothing to
                       % eliminate: the other equation is already free of it and IS the
                       % relation, which is what the elimination would have produced. When
                       % NEITHER does, the expression below is identically zero -- no extra cut,
                       % which is the right answer for an affine f (constant gradient, so the
                       % edge cell needs no interior side) and is what every case that reaches
                       % here today relies on.
                       n1 = size(c1,2);
                       n2 = size(c2,2);
                       if n1 == 2 && n2 == 2
                           ineq = c2(2)*eq1 - c1(2)*eq2;
                       elseif n1 == 2
                           ineq = eq2;
                       elseif n2 == 2
                           ineq = eq1;
                       else
                           ineq = c2*eq1 - c1*eq2;
                       end

                    end
                %else
                %    ineq = eq1+eq2;
                %end
        end
            
     end
     

     
     methods
       function [lg,limg] = limitOfGradientAtVertices (obj)
       % limg is PREALLOCATED AS SYM for the same reason region.limitOfFAtVertices is: an array
       % whose class is decided by its first assignment will silently round every exact value
       % written into it afterwards, and these values become the conjugate's cell boundaries.
           g = obj.f.gradient(obj.d.vars);
           lg = false(2, obj.d.nv);
           limg = sym(zeros(2, obj.d.nv));
           for i = 1:2   % Size of variables - change it
               [lg(i,:),limg(i,:)] = obj.d.limitOfFAtVertices (g(i));
           end

       end
     end
     methods % subdifferentials
      
       function [subdV,undV] = getSubdiffVertexT1 (obj, NCV, dualVars)
            subdV = sym(zeros(size(NCV,1),3));
            undV = zeros(obj.d.nv,1);
            [lg,limg] = obj.limitOfGradientAtVertices ;
           
            for j = 1:obj.d.nv %size(obj.d.vars,2)
              if ~lg(1,j)
                undV(j)=true;
                continue;
              end
              if ~lg(2,j)
                  undV(j)=true;
                  continue;
              end
            
              undV(j)=false;
              f = symbolicFunction(NCV(j,1));
              coef = f.getLinearCoeffs (dualVars);
              if (coef(1) == 0)
                 if coef(2) > 0
                 subdV(j,1) = dualVars(2)-limg(2,j);
                 else
                     subdV(j,1) = -(dualVars(2)-limg(2,j));
                 end
              elseif (coef(2) == 0) 
                subdV(j,1) = coef(1)*(dualVars(1) - limg(1,j));  
              elseif (coef(2) < 0)
                m = diff(NCV(j,1),dualVars(1));
                c = yIntercept(m, [limg(1,j),limg(2,j)]);
                subdV(j,1) = -1 * (dualVars(2) - m*dualVars(1) - c);
              else
                m = -diff(NCV(j,1),dualVars(1));
                c = yIntercept(m, [limg(1,j),limg(2,j)]);
                subdV(j,1) = dualVars(2) - m*dualVars(1) - c;
              end 
              f = symbolicFunction(NCV(j,2));
              coef = f.getLinearCoeffs (dualVars);
              if (coef(1) == 0)
                subdV(j,2) = dualVars(2)-limg(2,j);
                 if coef(2) > 0
                 subdV(j,2) = dualVars(2)-limg(2,j);
                 else
                     subdV(j,2) = -(dualVars(2)-limg(2,j));
                 end
              elseif (coef(2) == 0) 
                subdV(j,2) = coef(1) * (dualVars(1) - limg(1,j));  
              elseif (coef(2) < 0)
                m = diff(NCV(j,2),dualVars(1));
                c = yIntercept(m, [limg(1,j),limg(2,j)]);
                subdV(j,2) = -1 * (dualVars(2) - m*dualVars(1) - c);
              else
                m = -diff(NCV(j,2),dualVars(1));
                c = yIntercept(m, [limg(1,j),limg(2,j)]);
                subdV(j,2) = dualVars(2) - m*dualVars(1) - c;
              end 
            end
        end
  
        function [subdV,undV] = getSubdiffVertexT2 (obj, NCV, dualVars)
            subdV = sym(zeros(obj.d.nv,3));
            undV = zeros(obj.d.nv,1);
            [lg,limg] = obj.limitOfGradientAtVertices;
             
            for j = 1:obj.d.nv
              %  j
              if ~lg(1,j)
                  undV(j)=true;
                  continue;
              end
              if ~lg(2,j)
                  undV(j)=true;
                  continue;
              end
              undV(j)=false;
              f = symbolicFunction(NCV(j,1));
              coef = f.getLinearCoeffs (dualVars);
              if (coef(1) == 0)
                 if coef(2) > 0
                 %subdV(j,1) = dualVars(2)-limg(j,2);
                 subdV(j,1) = dualVars(2)-limg(2,j);
                 else
                     %subdV(j,1) = -(dualVars(2)-limg(j,2));
                     subdV(j,1) = -(dualVars(2)-limg(2,j));
                 end
              elseif (coef(2) == 0) 
                %subdV(j,1) = coef(1)*(dualVars(1) - limg(j,1));  
                subdV(j,1) = coef(1)*(dualVars(1) - limg(1,j));  
              elseif (coef(2) < 0)
                m = diff(NCV(j,1),dualVars(1));
                %c = yIntercept(m, [limg(j,1),limg(j,2)]);
                c = yIntercept(m, [limg(1,j),limg(2,j)]);
                subdV(j,1) = -1 * (dualVars(2) - m*dualVars(1) - c);
              else
                m = -diff(NCV(j,1),dualVars(1));
                %c = yIntercept(m, [limg(j,1),limg(j,2)]);
                c = yIntercept(m, [limg(1,j),limg(2,j)]);
                subdV(j,1) = dualVars(2) - m*dualVars(1) - c;
              end 
              %%%%%%%%%%%%%%%%%%%%%%%%%


              f = symbolicFunction(NCV(j,2));
              k = j+1;
              if k > obj.d.nv
                  k = 1;
              end
              coef = f.getLinearCoeffs (dualVars);
              %%%%%%%%%%%%%%%%%%%%%%%%%
              % if (coef(1) == 0)
              %    subdV(j,2) = dualVars(2)-limg(k,2);
              % elseif (coef(1) < 0)
              %   m = diff(NCV(j,2),dualVars(1));
              %   c = yIntercept(m, [limg(k,1),limg(k,2)]);
              %   subdV(j,2) = -1 * (dualVars(2) - m*dualVars(1) - c);
              % else
              %   m = -diff(NCV(j,2),dualVars(1));
              %   c = yIntercept(m, [limg(k,1),limg(k,2)]);
              %   subdV(j,2) = dualVars(2) - m*dualVars(1) - c;
              % end 
              %%%%%%%%%%%%%%%%%%%%%%%%%

              if (coef(1) == 0)
                subdV(j,2) = dualVars(2)-limg(2,k);
                 if coef(2) > 0
                 subdV(j,2) = dualVars(2)-limg(2,k);
                 else
                     subdV(j,2) = -(dualVars(2)-limg(2,k));
                 end
              elseif (coef(2) == 0) 
                subdV(j,2) = coef(1) * (dualVars(1) - limg(1,k));  
              elseif (coef(2) < 0)
                m = diff(NCV(j,2),dualVars(1));
                c = yIntercept(m, [limg(1,k),limg(2,k)]);
                subdV(j,2) = -1 * (dualVars(2) - m*dualVars(1) - c);
              else
                m = -diff(NCV(j,2),dualVars(1));
                c = yIntercept(m, [limg(1,k),limg(2,k)]);
                subdV(j,2) = dualVars(2) - m*dualVars(1) - c;
              end 
            end
        end

        function [subdV,undV] = getSubdiffVertexT2Q (obj, NCV, dualVars)
            subdV = sym(zeros(obj.d.nv,3));
            undV = zeros(obj.d.nv,1);
            %g = obj.envelope(i).f.gradient;
            [lg,limg] = obj.limitOfGradientAtVertices;

            for j = 1:obj.d.nv
              if ~lg(1,j)
                  undV(j)=true;
                  continue;
              end
              if ~lg(2,j)
                  undV(j)=true;
                  continue;
              end
              undV(j)=false;
              f = symbolicFunction(NCV(j,1));
              coef = f.getLinearCoeffs (dualVars);
              if (coef(1) == 0)
                 if coef(2) > 0
                 %subdV(j,1) = dualVars(2)-limg(j,2);
                 subdV(j,1) = dualVars(2)-limg(2,j);
                 else
                     %subdV(j,1) = -(dualVars(2)-limg(j,2));
                     subdV(j,1) = -(dualVars(2)-limg(2,j));
                 end
              elseif (coef(2) == 0) 
                %subdV(j,1) = coef(1)*(dualVars(1) - limg(j,1));  
                subdV(j,1) = coef(1)*(dualVars(1) - limg(1,j));  
              elseif (coef(2) < 0)
                m = diff(NCV(j,1),dualVars(1));
                %c = yIntercept(m, [limg(j,1),limg(j,2)]);
                c = yIntercept(m, [limg(1,j),limg(2,j)]);
                subdV(j,1) = -1 * (dualVars(2) - m*dualVars(1) - c);
              else
                m = -diff(NCV(j,1),dualVars(1));
                %c = yIntercept(m, [limg(j,1),limg(j,2)]);
                c = yIntercept(m, [limg(1,j),limg(2,j)]);
                subdV(j,1) = dualVars(2) - m*dualVars(1) - c;
              end 
              %%%%%%%%%%%%%%%%%%%%%%%%%


              f = symbolicFunction(NCV(j,2));
              k = j+1;
              if k > obj.d.nv
                  k = 1;
              end
              coef = f.getLinearCoeffs (dualVars);
              %%%%%%%%%%%%%%%%%%%%%%%%%
              % if (coef(1) == 0)
              %    subdV(j,2) = dualVars(2)-limg(k,2);
              % elseif (coef(1) < 0)
              %   m = diff(NCV(j,2),dualVars(1));
              %   c = yIntercept(m, [limg(k,1),limg(k,2)]);
              %   subdV(j,2) = -1 * (dualVars(2) - m*dualVars(1) - c);
              % else
              %   m = -diff(NCV(j,2),dualVars(1));
              %   c = yIntercept(m, [limg(k,1),limg(k,2)]);
              %   subdV(j,2) = dualVars(2) - m*dualVars(1) - c;
              % end 
              %%%%%%%%%%%%%%%%%%%%%%%%%

              if (coef(1) == 0)
                subdV(j,2) = dualVars(2)-limg(2,k);
                 if coef(2) > 0
                 subdV(j,2) = dualVars(2)-limg(2,k);
                 else
                     subdV(j,2) = -(dualVars(2)-limg(2,k));
                 end
              elseif (coef(2) == 0) 
                subdV(j,2) = coef(1) * (dualVars(1) - limg(1,k));  
              elseif (coef(2) < 0)
                m = diff(NCV(j,2),dualVars(1));
                c = yIntercept(m, [limg(1,k),limg(2,k)]);
                subdV(j,2) = -1 * (dualVars(2) - m*dualVars(1) - c);
              else
                m = -diff(NCV(j,2),dualVars(1));
                c = yIntercept(m, [limg(1,k),limg(2,k)]);
                subdV(j,2) = dualVars(2) - m*dualVars(1) - c;
              end 

              f = symbolicFunction(NCV(j,3));
              if f.isZero
                  continue
              end
              coef = f.getLinearCoeffs (dualVars);
              if (coef(1) == 0)
                 if coef(2) > 0
                 %subdV(j,1) = dualVars(2)-limg(j,2);
                 subdV(j,3) = dualVars(2)-limg(2,j);
                 else
                     %subdV(j,1) = -(dualVars(2)-limg(j,2));
                     subdV(j,3) = -(dualVars(2)-limg(2,j));
                 end
              elseif (coef(2) == 0) 
                %subdV(j,1) = coef(1)*(dualVars(1) - limg(j,1));  
                subdV(j,3) = coef(1)*(dualVars(1) - limg(1,j));  
              elseif (coef(2) < 0)
                m = diff(NCV(j,3),dualVars(1));
                %c = yIntercept(m, [limg(j,1),limg(j,2)]);
                c = yIntercept(m, [limg(1,j),limg(2,j)]);
                subdV(j,3) = -1 * (dualVars(2) - m*dualVars(1) - c);
              else
                m = -diff(NCV(j,3),dualVars(1));
                %c = yIntercept(m, [limg(j,1),limg(j,2)]);
                c = yIntercept(m, [limg(1,j),limg(2,j)]);
                subdV(j,3) = dualVars(2) - m*dualVars(1) - c;
              end 
            end
        end

        
        function subdV = getSubDiffVertexSpT1(obj, subdV, undV, crs)
          for j = 1:obj.d.nv
              if (~undV(j))
                  continue
              end
              em = j-1;
              if em == 0
                em = obj.d.nv;
              end
              ep = j+1;
              if ep == obj.d.nv+1
                ep = 1;
              end
              
              
              subdV(j,1) = subdV(em,1);
              subdV(j,2) = subdV(ep,2);
              subdV(j,3) = crs;  
          end    
        end
       
        function [subdE, unR, crs] = getSubDiffEdgeT1(obj, subdE, edgeNo, unDV, crs, dualvars)
          unR = zeros(obj.d.nv,1);
          for j = 1:obj.d.nv-1
             if unDV(j)
               unR(j) = false;
               continue
             end
             if unDV(j+1)
               continue
             end
              
             unR(j) = true;
             subdE(j,3) = crs;

          end    
          j = obj.d.nv;
          if unDV(obj.d.nv)
            unR(j) = false;
            return
          end
          if unDV(1) 
            return
          end
          unR(j) = true;

            
           subdE(j,3) = crs;
        end




     end

     methods % conjugate exprs quad
         function expr = conjugateExprVerticesT1 (obj, dualVars, unV )
            vars = obj.d.vars;
            for j = 1:obj.d.nv
                if unV(j)
                    expr(j) = obj.d.vx(j)*dualVars(1) + obj.d.vy(j)*dualVars(2) - obj.f.limit ( vars,[obj.d.vx(j),obj.d.vy(j)]).f;
                else
                    expr(j) = obj.d.vx(j)*dualVars(1) + obj.d.vy(j)*dualVars(2) - obj.f.subsF(vars,[obj.d.vx(j),obj.d.vy(j)]).f;
                end
                
            end
         end     

         % polyhedral
         function expr = conjugateExprEdgesT1Poly (obj, dualVars, edgeNo, psi0, psi1, psi2, expr )
            vars = obj.f.getVars;
            s1 = dualVars(1);
            s2 = dualVars(2);
            for j = 1:obj.d.nv
                no = edgeNo(j);
                mq = obj.d.ineqs(no).getLinearCoeffs (vars);
                if mq(2) == 0 
                    edgeCoef = obj.d.ineqs(no).getLinearCoeffs (vars);
                    c = -edgeCoef(3);
                    c2 = psi2(3)/(2*psi1(2));
                    c3 = -psi0(2)*c2;
                    c7 = psi1(2)*c2;
                    %c4 = -c3;
                    d2 = -psi1(2)*psi0(2)*c2+psi1(3)+c*psi1(1);
                    d3 = -psi0(2)^2*c2 + c*psi0(1) + psi0(3);
                    c5 = c7^2/psi2(3);
                    % document has c4 
                    %c6 = 2*c7*d2/psi2(0) + c4;
                    c6 = 2*c7*d2/psi2(3) - c3;
                    d4 = d2^2/psi2(3)+d3;
                    a = c2-c5;
                    b = c3-c6; % can be simplified = c3 - (2*c7*d2/psi2(0) - c3)
                    %d = -d4;
                    expr(obj.d.nv+j) = a*s2^2 + c*s1+b*s2-d4;
                else
                m = -mq(1)/mq(2);
                q = -mq(3)/mq(2);
                %psi2(1) + m*psi2(2)
                if psi2(1) + m*psi2(2) == 0
                
                  t0 = (-psi0(1)-m*psi0(2))/(2*(psi1(1)+m*psi1(2)));
                  t1 = 1/(2*(psi1(1)+m*psi1(2)));
                  t2 = m/(2*(psi1(1)+m*psi1(2)));
                  gamma10 = t1*(psi2(3)+q*psi2(2))/(psi1(1)+m*psi1(2));
                  gamma01 = t2*(psi2(3)+q*psi2(2))/(psi1(1)+m*psi1(2));
                  gamma00 = (t0*(psi2(3)+q*psi2(2))-psi1(3)-q*(psi1(2)))/(psi1(1)+m*psi1(2));
                  zeta11 = -(psi1(1)*gamma10+m*psi1(2)*gamma10)^2/(psi2(3)+q*psi2(2)) + gamma10;
                  zeta12 = -(2*(psi1(1)*gamma01+m*psi1(2)*gamma01)*(psi1(1)*gamma10+m*psi1(2)*gamma10))/(psi2(3)+q*psi2(2)) + gamma01 + m *gamma10;
                  zeta22 = -(psi1(1)*gamma01+m*psi1(2)*gamma01)^2/(psi2(3)+q*psi2(2)) + gamma01 *m;
                  zeta10 = -2*(psi1(1)*gamma01+m*psi1(2)*gamma10)*(psi1(3)+psi1(1)*gamma00+psi1(2)*(q+m*gamma00))/(psi2(3)+q*psi2(2)) - m*psi0(2)*gamma10 + gamma00 - psi0(1)*gamma10;
                  zeta01 = -(2*(psi1(1)*gamma01+m*psi1(2)*gamma01)*(psi1(3)+psi1(1)*gamma00+psi1(2)*(q+m*gamma00)))/((psi2(3)+q*psi2(2))) - m*psi0(2)*gamma01 - psi0(1)*gamma01 + m*gamma00+q;
                  zeta00 = -(psi1(3)+psi1(1)*gamma00 +psi1(2)*(q+m*gamma00))^2/(psi2(3)+q*psi2(2)) -psi0(3) - psi0(1)*gamma00 - psi0(2)*(q+m*gamma00);
                  expr(obj.d.nv+j) = simplify(zeta11*s1^2 + zeta12*s1*s2 + zeta22*s2^2 + zeta10*s1 + zeta01*s2 + zeta00);
                else
                  zeta00 = (psi2(1) + m*psi2(2))^2  ;
                  delta1 = -2*(psi1(3)*psi2(1) - psi1(1)*psi2(3) + m*psi1(3)*psi2(2) - m*psi1(2)*psi2(3) - q*psi1(1)*psi2(2) + q*psi1(2)*psi2(1))*(psi0(1)*psi2(1) + psi1(1)^2 + m*(psi1(2)^2*m + psi0(1)*psi2(2) + psi0(2)*psi2(1) + 2*psi1(1)*psi1(2) + psi0(2)*psi2(2)*m))  ;
                  delta0 = 2*psi1(1)*psi2(3)*(psi1(1)+m*psi1(2)) -2*psi1(3)*psi2(1)*(psi1(1)+m*psi1(2))-psi0(3)*psi2(1)*(psi2(1)+m*psi2(2)) ...
                  + psi0(1)*psi2(3)*(psi2(1)+m*psi2(2)) -2*m*psi1(3)*psi2(2)*(psi1(1)+m*psi1(2)) + 2*m*psi1(2)*psi2(3)*(psi1(1)+m*psi1(2)) ...
                  - m * psi0(3)*psi2(2) * (psi2(1)+m*psi2(2)) + m * psi0(2)*psi2(3)*(psi2(1)+m*psi2(2) ) + 2*q*psi1(1)*psi2(2)*(psi1(1)+m*psi1(2)) ...
                  - 2*q*psi1(2)*psi2(1)* (psi1(1)+m*psi1(2)) + q*psi0(1)*psi2(2)*(psi2(1)+m*psi2(2))-q*psi0(2)*psi2(1)*(psi2(1)+m*psi2(2));
                  si1 = 2*(psi2(1) + m*psi2(2)) * (psi1(3)*psi2(1) - psi1(1)*psi2(3) + m*psi1(3)*psi2(2) - m*psi1(2)*psi2(3) - q*psi1(1)*psi2(2) + q*psi1(2)*psi2(1))*s1 + 2*m*(m*psi2(2) + psi2(1)) * (psi1(3)*psi2(1) - psi1(1)*psi2(3) + m*psi1(3)*psi2(2) - m* psi1(2)*psi2(3) - q*psi1(1)*psi2(2) + q*psi1(2)*psi2(1))*s2 + delta1;
                  si1_2 = -(psi2(1) + m*psi2(2))*s1 -m *(psi2(1)+m*psi2(2))*s2+psi0(1)*psi2(1) +psi1(1)^2 + m*(m*psi1(2)^2 + psi0(1)*psi2(2)+psi0(2)*psi2(1)+2*psi1(1)*psi1(2)+m*psi0(2)*psi2(2));
                  si0 = (-psi2(3)*(psi2(1)+m*psi2(2))-q*psi2(2)*(psi2(1)+m*psi2(2)))*s1 + (q*psi2(1)*(psi2(1)+m*psi2(2))-m*psi2(3)*(psi2(1)+m*psi2(2)))*s2 + delta0;  
                  expr(obj.d.nv+j) = simplify((si1 / (zeta00 * sqrt(si1_2))) + si0);
                end
                end
            end
        end

        

     end

    
   
     

     methods (Static)
         function [edgeNo, keepE] = spreadCollidingEdges(d, edgeNo, keepE, nOn, wasBounded) %#ok<INUSD>
         % Give two GENUINE edges that collided on one edge number distinct numbers, instead of
         % letting conjugateOfPiecePoly's scatter destroy one. See the call site for the shape and
         % the measurement.
         %
         % ONLY constraints with at least TWO of the region's vertices on them count as edges
         % here. A constraint active at one vertex only is either a RAY of an unbounded region --
         % which getEdgeNosInf already reserves slot 1 for, so it never collides -- or the
         % edge-less kind the three rules above exist to drop. Leaving those alone is what keeps
         % this from touching any piece the old code handled.
         %
         % THE NUMBER EACH GETS. getEdgeNosInf's convention is "edge j starts at vertex j", so two
         % edges sharing vertices j1 and j2 want those two indices -- one traversed out of each.
         % The group keeps the number it already has for its first member and the others take the
         % index of another vertex they are actually on, so this can never invent an edge where
         % the region has none.
         %
         % SCOPE. An extra is moved only to a number no surviving constraint already claims. If
         % there is none, nothing changes and the old behaviour stands, so this can only ever add
         % an edge that was being destroyed -- never remove or renumber one that was fine.
             if isempty(edgeNo), return, end

             % SCOPE. Everything below runs only when some edge number is claimed by two or more
             % constraints that each have at least two of the region's vertices on them. That is
             % the lens signature; nothing else reaches it, which is what keeps a change to
             % vendored, index-fragile code away from pieces that were already right.
             still = find(keepE(:))';
             isLens = false;
             for sIdx = unique(edgeNo(still))'
                 g = still(edgeNo(still) == sIdx);
                 if numel(g(nOn(g) >= 2)) >= 2, isLens = true; break, end
             end
             if ~isLens, return, end

             % (i) Within a colliding group, drop any constraint whose curve between the shared
             % vertices LEAVES the region: it meets the region only at those vertices and bounds
             % no edge. On the half-lens that is the SECOND conic, which passes through both
             % vertices and is redundant -- and which LP redundancy cannot see, because it is not
             % linear.
             for sIdx = unique(edgeNo(still))'
                 g = still(edgeNo(still) == sIdx);
                 g = g(nOn(g) >= 2);
                 if numel(g) < 2, continue, end
                 for k = g
                     if keepE(k) && ~functionNDomain.boundsAnEdge(d, k, keepE)
                         keepE(k) = false;
                     end
                 end
             end

             % (ii) NOT DONE: freeing a slot by dropping the constraint that holds it.
             % On the pipeline's own lens the two numbers the real edges need are held by
             % constraints with a single vertex on them -- and dropping those was tried and is
             % UNSOUND. A constraint active at exactly one vertex of a convex region can still be
             % essential: removing it enlarges the region, and an enlarged piece of f* has a
             % SMALLER conjugate domain. Measured on f = x*y over the two-face square: with the
             % drop, f** is exact at (0.25,0.25) and (0.1,0.1) and +inf at (0.9,0.6) and (0.6,0.6);
             % without it, exactly the other way round. Both are 5 of 7, and only one of them is
             % sound, so the drop is out.
             % Deciding boundedness correctly does NOT rescue it -- the harmed piece is genuinely
             % bounded too. And it is not needed: what a lens actually needs is an explicit EDGE
             % LIST rather than a slot convention, which conjugateOfPiecePoly now derives when
             % this routine leaves a collision unresolved (functionNDomain.edgeIndexList). So no
             % constraint is dropped, and the slot the two edges were fighting over stops
             % mattering.

             still = find(keepE(:))';
             for sIdx = unique(edgeNo(still))'
                 grp = still(edgeNo(still) == sIdx);
                 grp = grp(nOn(grp) >= 2);
                 if numel(grp) < 2, continue, end
                 for gi = 2:numel(grp)
                     k = grp(gi);
                     [~, vxk, vyk] = d.vertexOfEdge(k);
                     for t = 1:numel(vxk)
                         cand = 0;
                         for j = 1:d.nv
                             if isAlways(d.vx(j) == vxk(t)) && isAlways(d.vy(j) == vyk(t))
                                 cand = j; break
                             end
                         end
                         if cand == 0 || cand == sIdx, continue, end
                         if any(keepE(:) & edgeNo(:) == cand), continue, end   % already claimed
                         edgeNo(k) = cand;
                         break
                     end
                 end
             end
         end

         function [eIdx, ok] = edgeIndexList(d, wasBounded)
         % WHICH CONSTRAINT BOUNDS WHICH EDGE, read off the region's own geometry rather than off
         % a slot. eIdx(j) is the index into d.ineqs of the constraint bounding the edge from
         % vertex j to vertex j+1 (cyclically, when the region is bounded). ok is false when the
         % geometry does not settle the whole list, and the caller then leaves the slot
         % conventions exactly as they are.
         %
         % WHY THIS EXISTS. Everything downstream in conjugateOfPiecePoly is edge-indexed --
         % getNormalConeVertexQ reads the constraints of the edges meeting at a vertex,
         % getNormalConeEdgeQ/Q3 the constraint of each edge, and the edge loop that same
         % constraint again -- under one of two slot conventions chosen by the COUNT
         % `size(ineqs,2) == nv`. That count stands in for "is this region unbounded", and on a
         % LENS it is simply the wrong question: two vertices, two genuine edges (an arc and its
         % chord), five constraints. The count says unbounded, so edge 1 is read from ineqs(2),
         % the arc is never read as an edge at all, and what comes out is the conjugate of the
         % CHORD -- finite on a strip where the truth is finite everywhere.
         %
         % HOW THE LIST IS BUILT. A constraint bounds the edge between two consecutive vertices
         % when both vertices lie on it AND its own curve between them stays in the region --
         % region.vertexOfEdge answers the first, boundsAnEdge the second (the second matters:
         % the lens's OTHER conic passes through both vertices too, and meets the region nowhere
         % else). Where several qualify, today's slot is preferred whenever it is among them, so
         % a piece whose indexing already works keeps exactly the list it has.
             eIdx = []; ok = false;
             if isempty(d) || d.nv < 2
                 return
             end
             nv = d.nv;
             m  = size(d.ineqs,2);
             if m < 1
                 return
             end
             % BOUNDED ONLY. The list is consumed by getNormalConeVertexQ, which walks all nv
             % vertices and asks for the edge leaving each one -- so it needs a CLOSED cycle,
             % nE = nv. An unbounded region has nv-1 edges between its finite vertices and two
             % rays, the first vertex has no predecessor edge, and there is no entry for the
             % caller to look up; indexing eIdx(nv) there would run off the end. Unbounded lenses
             % are not a shape the pipeline has produced -- a lens is bounded by construction --
             % so refuse rather than invent a convention for one.
             if ~wasBounded
                 return
             end
             nE = nv;
             % today's convention, kept as the preferred answer wherever it is still valid
             if m == nv
                 prefer = 1:nv;
             else
                 prefer = 2:(nv+1);
             end

             % which vertices lie on which constraint
             onV = false(m, nv);
             for k = 1:m
                 [~, vxk, vyk] = d.vertexOfEdge(k);
                 for t = 1:numel(vxk)
                     for vI = 1:nv
                         if isAlways(d.vx(vI) == vxk(t)) && isAlways(d.vy(vI) == vyk(t))
                             onV(k, vI) = true;
                         end
                     end
                 end
             end

             eIdx = zeros(1, nE);
             for j = 1:nE
                 a = j;
                 b = mod(j, nv) + 1;
                 cand = find(onV(:,a)' & onV(:,b)');
                 cand = cand(~ismember(cand, eIdx(1:j-1)));
                 % PASSING THROUGH BOTH VERTICES IS NOT ENOUGH. On the lens BOTH conics do -- the
                 % second one meets the region at those two points and nowhere else. Keep only
                 % constraints whose own curve between the vertices stays inside; boundsAnEdge
                 % reports anything undecidable as an edge, so this only ever drops a candidate it
                 % has positive evidence against. Filtering BEFORE the preference below matters:
                 % preferring today's slot unconditionally would hand edge 2 of the lens that
                 % second conic, which bounds nothing.
                 cand = cand(arrayfun(@(k) functionNDomain.boundsAnEdge(d, k, true(1, m)), cand));
                 if isempty(cand)
                     eIdx = []; return
                 end
                 p = prefer(j);
                 if p <= m && any(cand == p)
                     eIdx(j) = p;
                 else
                     eIdx(j) = cand(1);
                 end
             end
             ok = true;
         end

         function tf = edgesStillCollide(edgeNo, keepE, nOn)
         % Do two constraints that each bound a genuine EDGE still share one edge number, after
         % spreadCollidingEdges has moved everything it can?
         %
         % This is the LENS signature and the only condition under which conjugateOfPiecePoly
         % departs from its slot conventions. It is deliberately narrow: when it is false the
         % scatter stores every edge in a distinct slot, the two conventions are consistent, and
         % nothing about the piece changes.
             tf = false;
             still = find(keepE(:))';
             if isempty(still)
                 return
             end
             for sIdx = unique(edgeNo(still))'
                 g = still(edgeNo(still) == sIdx);
                 if numel(g(nOn(g) >= 2)) >= 2
                     tf = true; return
                 end
             end
         end

         function tf = boundsAnEdge(d, k, keepE)
         % Does constraint k bound a genuine EDGE of d, or does it only touch the region at
         % vertices? Decided from a point strictly BETWEEN two of its vertices, ON its own curve --
         % the rule this repository keeps re-learning (QuaPar.chordCuts, pieceRecessionRays), and
         % for the same reason: the shared vertices lie on every candidate at once and so say
         % nothing. Anything undecidable is reported as a genuine edge, so this can only ever drop
         % a constraint it has positive evidence against.
             tf = true;
             M = functionNDomain.pointBetweenOnCurve(d, k);
             if isempty(M), return, end
             for j = 1:numel(keepE)
                 if j == k || ~keepE(j), continue, end
                 try
                     v = d.ineqs(j).subsF(d.vars, M);
                     if isAlways(v.f > 0), tf = false; return, end
                 catch
                     return                      % undecidable: keep k
                 end
             end
         end

         function M = pointBetweenOnCurve(d, k)
         % A point ON {ineqs(k) = 0} strictly between two of the region's vertices that lie on it.
         % Empty when there is no such point this can produce, which the caller reads as
         % "undecidable".
             M = [];
             [nvk, vxk, vyk] = d.vertexOfEdge(k);
             if nvk < 2, return, end
             A = [vxk(1), vyk(1)];
             B = [vxk(2), vyk(2)];
             mid = (A + B)/2;
             c = d.ineqs(k);
             try
                 if isAlways(c.subsF(d.vars, mid).f == 0)
                     M = mid; return                 % straight edge: the chord IS the curve
                 end
             catch
             end
             try
                 t = sym('t_pbc');
                 dir = [-(B(2)-A(2)), B(1)-A(1)];    % perpendicular bisector of AB
                 P = mid + t*dir;
                 sol = solve(c.subsF(d.vars, P).f == 0, t, 'Real', true);
                 if isempty(sol), return, end
                 sol = double(sol);
                 sol = sol(isfinite(sol));
                 if isempty(sol), return, end
                 [~, i0] = min(abs(sol));            % the root nearest the chord
                 M = sym(double(mid) + sol(i0)*double(dir));
             catch
                 M = [];
             end
         end

         function objR = maxOfList(groups)
         % THE pointwise maximum of several piecewise functions. `groups` is a cell array; each
         % entry is a functionNDomain array whose cells are disjoint (one piecewise function).
         % Returns their pointwise max as a single functionNDomain array.
         %
         % This is the max step of the [JOGO]/[COAP] algorithm, and it is used TWICE -- once to
         % assemble f* from the per-triangle conjugates, and once to assemble f** from the
         % per-piece conjugates of f*. It lives here, in one place, so the two callers
         % (plq.maximumConjugate and plq.biconjugateF) cannot drift apart.
         %
         % The conjugate's max is the harder of the two: its operands can carry parabolic edges,
         % so the result is a general quadratic-on-parabolic subdivision. The biconjugate's max
         % only ever sees polyhedral operands and returns a polyhedral subdivision -- a strict
         % special case. So the conjugate's implementation covers both, and there is nothing to
         % specialize for the biconjugate.
         %
         % `*` here is functionNDomain.mtimes, which does NOT multiply: it intersects each pair
         % of regions and stores BOTH functions on the common region, so that maximumP has two
         % values to compare on one domain. That pairing is what makes the max computable.
             objR = functionNDomain.empty();
             if isempty(groups), return, end
             first = true;
             for j = 1:numel(groups)
                 g = groups{j};
                 if isempty(g), continue, end
                 if first
                     objR = g;
                     first = false;
                 else
                     objR = objR * g;
                     objR = objR.maximumP(true);
                 end
             end
         end
     end

end
