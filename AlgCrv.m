classdef AlgCrv < AlgCon & Crv
   % AlgCrv: an ALGEBRAIC function on an ALGEBRAIC subdivision -- the type `biconj` returns when the
   % input's pieces meet along a CURVED tie.
   %
   % It is the (Alg, Crv) cell of the lattice, one step above AlgCon on the SUBDIVISION axis:
   %
   %       AlgCrv  (algebraic on an algebraic subdivision)   < AlgCon & Crv   <- this class
   %          |
   %       AlgCon  (algebraic on conic)                      < Alg & Conic
   %          |
   %       RatCon  (rational on conic)                       < AlgCon & Rat
   %
   % WHY IT EXISTS, in one line: the convex envelope's BRIDGE face carries the value -- measured
   % strictly above the best quadratic face at 60 of 60 sampled points -- and the boundary between
   % it and a quadratic face has degree 4, which a conic subdivision cannot hold as an edge.
   % AlgCon.m records why the FUNCTION axis had to grow; Crv.m records why this one did.
   %
   % THE ALGEBRAIC EDGE TABLE IS `AlgCon.EcA`, declared in the one place the mesh lives and pinned
   % empty unless the object is a Crv -- exactly the rule `den`, `Ec`, `FP` and `FB` already follow.
   % This class is the concrete cell where that pin is released, plus the evaluator that reads such
   % an edge. The EXACT mesh it also needs is declared below rather than inherited, for a reason
   % that is measured and recorded there and in TODO.md; that is a wart, and it is a wart with a
   % price tag rather than an oversight.
   %
   % AN FC ROW MAY NOW NAME EITHER TABLE. A row is [edgeIndex, side, branch, table], where table 0
   % (or absent, which is every row that predates this class) means `Ec`/`EcQ`'s CONIC and table 1
   % means `EcA`'s algebraic curve. For a conic row the condition is the one QuaCon documents; for
   % an algebraic row it is
   %
   %        side * ( q(x) - value_FB(x) ) >= 0 ,
   %
   % with q and FP/FB taken from the EcA entry. `algCurveSide` decides that exactly at a rational
   % point, by counting how many roots of the fibre lie at or below q -- no root is ever computed
   % and no irrational value is ever formed. `eval` below reads it in DOUBLES, deliberately and for
   % the same reason QuaCon.eval does: the operations that must leave Q are the PREDICATES that
   % decide a mesh, and reading a finished function at a point is not one of them.

   properties
        % THE EXACT MESH, DECLARED HERE AND NOT INHERITED -- which is worth saying plainly, because
        % AlgCon's header says the mesh is declared once and this looks like a breach of it.
        %
        % It is not where it should be, and the cost of moving it was MEASURED rather than assumed.
        % `fN`, `fD`, `EcQ`, `Vname` and `FC` live on QuaCon, not on AlgCon, and this class cannot
        % extend QuaCon because a QuaCon is a `Rat` and a Rat may not carry an implicit face
        % (`AlgCon:FPMustBeEmpty`). Moving the five declarations up to AlgCon -- where the header
        % says they belong -- was tried on 2026-09-20 and took the fast bucket from 692 passing to
        % 363 passing, 331 failed and 326 INCOMPLETE, i.e. classes that no longer load. So they are
        % declared here, duplicated, until that move is done properly as its own piece of work with
        % its own measurement. TODO.md carries it.
        fN (:,10){mustBeNumeric}   % nf x 10 integer numerators, ratQ.canon form against fD
        fD (:,1) {mustBeNumeric}   % nf x 1  positive integer denominators
        EcQ (:,6){mustBeNumeric}   % ne x 6  primitive integer conics, ratQ.conic form
        Vname (:,3){mustBeInteger,mustBeNonnegative} % nv x 3 [edgeA edgeB rootIdx]
        FC                         % cell(nf,1); FC{k} is m x 4 of [edgeIndex, sign, branch, table]
                                   % -- table 0 selects the conic EcQ, table 1 the algebraic EcA.
                                   % Fewer columns are accepted and read with the rest zero.
   end

   methods
      function obj = AlgCrv(Vname, E, EcQ, fN, fD, F, FC, FP, FB, EcA)
      % objective: the mesh, exactly as AlgCon takes it, plus the algebraic edge table.
      % [input]  as RatCon's constructor, then FP/FB (the implicit faces) and EcA (the new edges)
      %
      % THE CONSTRUCTOR PROTOCOL of RatCon applies unchanged: the LEAF writes the mesh and then
      % calls setPinnedDefaults, because MATLAB re-runs a shared base constructor once per
      % inheritance path and the second run would otherwise clobber the first.
         if nargin == 0, return, end
         if nargin < 8,  FP = {};  end
         if nargin < 9,  FB = [];  end
         if nargin < 10, EcA = struct('FP', {}, 'FB', {}, 'num', {}, 'den', {}); end

         obj.nv = size(Vname,1);  obj.nf = numel(FC);  obj.ne = size(EcQ,1);
         obj.Vname = Vname;  obj.E = E;  obj.EcQ = EcQ;
         obj.fN = fN;  obj.fD = fD;  obj.F = F;  obj.FC = FC;
         obj.FP = FP;  obj.FB = FB;
         obj.EcA = EcA;
         obj = obj.setPinnedDefaults();
      end

      function [val, idx] = eval(obj, X)
      % objective: the function's value at points X, +Inf where no face owns the point.
      % [input]  X : n x 2
      % [output] val : n x 1;  idx : n x 1 face index, 0 where none owns it
      %
      % A HOLE IS +INFINITY, not a gap: outside `dom f**` the supremum diverges, and leaving the
      % cell out of the mesh is how every producer here says so.
         if nargin ~= 2
            error('AlgCrv:nargin', 'expects (obj, X).');
         end
         if ~isnumeric(X) || size(X,2) ~= 2
            error('AlgCrv:points', 'X must be n x 2; got %s.', mat2str(size(X)));
         end
         n = size(X,1);
         val = inf(n,1);  idx = zeros(n,1);
         for i = 1:n
            for k = 1:obj.nf
               if ~obj.faceOwns(k, X(i,:)), continue, end
               v = obj.faceValue(k, X(i,:));
               if isnan(v), continue, end
               if isinf(val(i)) || v > val(i), val(i) = v;  idx(i) = k; end
            end
         end
      end
   end

   methods (Access = private)
      function tf = faceOwns(obj, k, x)
      % Every row of face k, conic or algebraic. The tolerance is QuaCon.eval's: a point within it
      % of a boundary is accepted, and on a boundary two faces agree anyway, the envelope being
      % continuous there.
         tf = true;
         rows = obj.FC{k};
         for t = 1:size(rows,1)
            tbl = 0;
            if size(rows,2) >= 4, tbl = rows(t,4); end
            if tbl == 1
               a = obj.EcA(rows(t,1));
               v = algBranch(a.FP, x, a.FB);
               if ~isfinite(v), tf = false; return, end
               q = ratQ.evalFace(a.num, a.den, x);
               if rows(t,2) * (q - v) < -1e-9, tf = false; return, end
            else
               c = obj.EcQ(rows(t,1), :);
               br = 0;
               if size(rows,2) >= 3, br = rows(t,3); end
               if br ~= 0
                  v = ratQ.pairLineEval(c, br, x);
               else
                  v = c(1)*x(1)^2 + c(2)*x(1)*x(2) + c(3)*x(2)^2 + ...
                      c(4)*x(1) + c(5)*x(2) + c(6);
               end
               if rows(t,2) == 0
                  if abs(v) > 1e-9, tf = false; return, end      % an EQUALITY row: a THIN face
               elseif rows(t,2) * v < -1e-9
                  tf = false; return
               end
            end
         end
      end

      function v = faceValue(obj, k, x)
      % FP/FB when face k carries no formula, fN/fD when it does -- AlgCon's own rule.
         if ~isempty(obj.FP) && numel(obj.FP) >= k && ~isempty(obj.FP{k})
            v = algBranch(obj.FP{k}, x, obj.FB(k));
         else
            v = ratQ.evalFace(obj.fN(k,:), obj.fD(k), x);
         end
      end
   end
end
