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

           objR2 = mergeL(objR);
            objR3 = mergeL(objR2);
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
               %%%%%%%%%%%%%%%
            nP = 0;
            for j = 1:r.nv
              if abs(r.vx(j)) == intmax
                continue
              end
              if abs(r.vy(j)) == intmax
                continue
              end
              nP = nP+1;
              px(nP) = r.vx(j);
              py(nP) = r.vy(j);
            end
            %%%%%%%%%%%%%%%%%%%%%%
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
            

            for i=1:size(obj,2)
                d = obj(i).d;
                
                d = d.removeTangent(d.nv,d.vx,d.vy);
                d = d.removeInfV;
                %d = d.poly2orderUnbounded;
               
                if d.nv == size(d.ineqs,2)
                    d = d.poly2order;
                else
                    d = d.poly2orderUnbounded;
                end
                edgeNo = d.getEdgeNosInf(d.vars);
                d.ineqs(edgeNo) = d.ineqs;

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
               if f.isQuad
                   
                   for j = 1:size(d0.ineqs,2)
                       if (d0.ineqs(j).isQuad)
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
               % fix this for quad
               if size(d.ineqs,2) == d.nv
                   NCV = d.getNormalConeVertex(x, y);
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
                   if size(d0.ineqs,2) == d.nv
                       NCE = d0.getNormalConeEdgeQ3(x, y);
                       [subdE,unR] = obj(i).getSubdiffVertexT2 (NCE, [x,y]);
                   else
                       NCE = d0.getNormalConeEdgeQ(x, y);
                       [subdE,unR] = obj(i).getSubdiffVertexT2Q (NCE, [x,y]);
                   end
                   endNv = obj(i).d.nv-1;
                   if size(obj(i).d.ineqs,2) == obj(i).d.nv
                       endNv = obj(i).d.nv;
                   end
                   for j = 1:endNv % fix this  obj(i).d.nv-1 ?
                       if endNv ==  obj(i).d.nv
                           ineq = subs(obj(i).d.ineqs(j).f,obj(i).d.vars,[x,y]);
                       else
                           ineq = subs(obj(i).d.ineqs(j+1).f,obj(i).d.vars,[x,y]);
                       end
                   
                       f0 = subs(obj(i).f.f,obj(i).d.vars,[x,y]);
                       [expr] = simplifyFraction(conjugateExpr(ineq,f0,x,y));
                       ineq1 = subdE(j,:);
                   if obj(i).f.isQuad
                       edgeInt = obj(i).getInterior(x,y) ;
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
                       edgeInt = obj(i).getInterior(x,y) ;
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

        function ineq = getInterior(obj,x,y)
                g(1) = obj.f.dfdx(obj.d.vars(1));
                g(2) = obj.f.dfdx(obj.d.vars(2));
                eq1 = x - g(1).f;
                eq2 = y - g(2).f;
                %if obj.f.isQuad
                    s12 = solve([eq1,eq2],obj.d.vars);
                    if isempty (s12) | isempty(s12.s_1)
                       c1 = coeffs(eq1,obj.d.vars(1));
                       c2 = coeffs(eq2,obj.d.vars(1));
                       
                       if size(c1,2) == 2
                           ineq = c2(2)*eq1 - c1(2)*eq2;
                       else
                           ineq = c2*eq1 - c1*eq2;
                       end    
                       %ineq = c2(2)*eq1 - c1(2)*eq2;
                       
                    else
                    ineq = subs(eq2,obj.d.vars,s12);
                    end
                %else
                %    ineq = eq1+eq2;
                %end
        end
            
     end
     

     
     methods
       function [lg,limg] = limitOfGradientAtVertices (obj)
           g = obj.f.gradient(obj.d.vars);
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

    
   
     
end