classdef plq
  properties
      nPieces = 0;
      % Untyped (not plq_1p.empty()/plq_1piece.empty()): plq_1piece is an older, independently-
      % tested parallel implementation of the same per-piece API (testMaxMultiRegion.m constructs
      % plq objects from plq_1piece pieces directly), so this property must accept either class.
      pieces = [];
      maxConjugate = functionNDomain.empty();
      biconjugate = functionNDomain.empty();
      biconjfia = [];
      lMax = false
      lBiconj = false

  end
  % 11 methods
  methods
      function obj = plq(ps)
          % Assign the whole array in one shot (not element-by-element indexed assignment into
          % obj.pieces): obj.pieces starts as literal [] (class double), so obj.pieces(i)=ps(i)
          % tries to convert ps(i) (a plq_1p or plq_1piece object) to double and fails. A single
          % whole-array assignment lets obj.pieces take on ps's own class instead.
          obj.nPieces = size(ps,2);
          obj.pieces = ps;
      end

      % io
      function print(obj)
          disp("")
          disp("")
          disp("")
          disp("")
          disp("")
          for i = 1:obj.nPieces
             disp(["Piece ", num2str(i)])
              obj.pieces(i).print;
          end
          if ~obj.lMax
              return
          end
          obj.maxConjugate.printL
          if ~obj.lBiconj
              return
          end
          
          obj.biconjugate.printL
      end

      function printLatex(obj)
          disp("")
          disp("")
          disp("")
          disp("")
          disp("")
          for i = 1:obj.nPieces
             fprintf("Piece " + num2str(i) + "\n")
             disp(" ")
              obj.pieces(i).printLatex;
          end
          %return
          if ~obj.lMax
              return
          end
          
          obj.maxConjugate.printLLatex
          if ~obj.lBiconj
              return
          end
          
          obj.biconjugate.printL
      end

      function plotMaxd(obj)
          figure;
          for i =1:size(obj.maxf,1)
            obj.maxd(i,1).plotRegion;
          end
      end

      
      function plotDomain(obj)
           %figure;
           for i = 1:obj.nPieces
               obj.pieces(i).plotDomain
           end
      end

      
      function printDomainMaple(obj)
          disp('printDomainMaple')
          for i=1: obj.nPieces
            obj.pieces(i).Mprint;
          end
          %return
          if ~obj.lMax
              return
          end
          disp('maxConj')
          % HISTORY: lMax/lBiconj only track whether .maximum/.biconjugateF
          % were called, not whether they produced a non-empty result --
          % functionNDomain.addEq (used by biconjugateF) can legitimately
          % return an empty array for a domain where the merged biconjugate
          % comes out empty (see testMaxMultiRegion/testMaxThesis), and
          % printM unconditionally indexes objL(1), erroring on an empty
          % array rather than just having nothing to print.
          if ~isempty(obj.maxConjugate)
            obj.maxConjugate.printM;
          end
          %obj.maxConjugate.printM2;
          if ~obj.lBiconj
              return
          end
          if ~isempty(obj.biconjugate)
            obj.biconjugate.printM;
          end
          %obj.maxConjugate.printM2;

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


      
      function obj = triangulate(obj)
        pieces = plq_1p.empty();
        for i = 1:obj.nPieces
            pieces = obj.pieces(i).triangulate (pieces);
        end
        obj = plq(pieces);
      end

      function obj = convexEnvelope(obj)

        for i=1: obj.nPieces
           
          obj.pieces(i)=obj.pieces(i).convexEnvelope;
          
        end
        return
      end


      function obj = conjugate(obj)

        for i=1: obj.nPieces
          obj.pieces(i)=obj.pieces(i).convexEnvelope;

          obj.pieces(i)=obj.pieces(i).conjugate;
        end
        return
      end



      function obj = maximum(obj)

        for i=1: obj.nPieces
          
          obj.pieces(i)=obj.pieces(i).convexEnvelope;
          obj.pieces(i)=obj.pieces(i).conjugate;
          obj.pieces(i) = obj.pieces(i).maximumConjugate;
        end
        
        obj = obj.maximumConjugate;
        obj.lMax = true;
      end

      
      function obj = maximumConjugate(obj)
      % Step 3 of the conjugate: the max of the per-piece conjugates. The loop that used to be
      % written out here now lives in functionNDomain.maxOfList, because biconjugateF needs the
      % SAME max and must not carry a second copy of it.
          groups = cell(1, obj.nPieces);
          for j = 1:obj.nPieces
              groups{j} = obj.pieces(j).maxConjugate;
          end
          obj.maxConjugate = functionNDomain.maxOfList(groups);
      end

      function obj = biconjugateF(obj)
         % for i = 1:obj.nPieces
         %    obj.pieces(i).biconjugateP
         % end
        [bc,ia] = obj.maxConjugate.conjugateOfPiecePoly;

        % HISTORY: this was `bc.mergeL` then `bc.addEq`, and it performed NO MAX -- so the last
        % step of the algorithm ("max of all those conjugates to obtain f**") was simply absent.
        % conjugateOfPiecePoly returns the conjugates of the pieces of f* CONCATENATED, with
        % ia(k):ia(k+1)-1 delimiting piece k's own cells; cells within one block are disjoint,
        % cells in different blocks overlap because they belong to different functions. addEq
        % groups by equal function ACROSS blocks and adds their equations, which intersects
        % regions belonging to different conjugates: on the unit box with f = x*y the four cells
        % carrying x+y-1 come from four different blocks and their intersection degenerates to
        % x = 1, so the group is dropped and the biconjugate comes back EMPTY. (It is worse in
        % pristine cPLQ, where the same input errors outright in addEq.)
        %
        % The fix is to run the max the algorithm calls for, and to run the SAME max Step 3 uses
        % -- functionNDomain.maxOfList, shared with maximumConjugate above. The blocks are
        % exactly the groups it wants, which is what ia is for.
        groups = cell(1, numel(ia)-1);
        for k = 1:numel(ia)-1
            groups{k} = bc(ia(k):ia(k+1)-1);
        end
        obj.biconjugate = functionNDomain.maxOfList(groups);
        obj.lBiconj = true;
        return
% 
     end
  end

   
end
