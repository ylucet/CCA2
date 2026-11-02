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
   % THE ALGEBRAIC EDGE TABLE IS `AlgCon.EcA`, declared in the one place the subdivision lives and pinned
   % empty unless the object is a Crv -- exactly the rule `den`, `Ec`, `FP` and `FB` already follow.
   % This class is the concrete cell where that pin is released, plus the evaluator that reads such
   % an edge. The EXACT subdivision it also needs (`fN`, `fD`, `EcQ`, `Vname`, `FC`) is inherited from
   % AlgCon, where it is declared once; until 2026-09-22 it was duplicated here (DECISIONS).
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
   % decide a subdivision, and reading a finished function at a point is not one of them.
   %
   % TABLE 2 IS `EcI`, A TIE BETWEEN TWO IMPLICIT VALUES (2026-10-31): value_1(x) >= value_2(x), with
   % FP1/FB1 and FP2/FB2 taken from the entry. A minimum over pairs orders an irrational face of one
   % pair against an irrational face of another, and the curve where the two tie is never formed --
   % the row names the two branches. It lives here and not on AlgCon because only this type has it;
   % `algCrvVertices` and `algCrvIncidence` refuse it by name, since a tie row has no polynomial
   % they could intersect.
   %
   % TABLE 3 IS `Tri`, THE BICONJUGATE'S TRIANGLE FACES (2026-11-01). A row [t e side 3] says
   % side * s_e(x) >= 0, s_e the e-th side of triangle Tri(t) oriented +1 toward the opposite corner
   % -- triangleFace's struct, whose corners are affine in a NAMED dual point s0 (triangleFaceSide
   % decides a side exactly at a rational point; eval reads the corner cache in doubles). A face
   % with FT(k) = t > 0 IS triangle t: its value is <s0, x> - g(s0), read from the cache. Erratum
   % Example 2's degree-4 triangle needs both; a triangle with a rational s0 is an ordinary affine
   % face on three lines and uses neither.

   properties
      EcI (:,1) struct = struct('FP1', {}, 'FB1', {}, 'FP2', {}, 'FB2', {})   % table 2: value_1 >= value_2
      Tri (:,1) struct = struct('kind', {}, 's', {}, 'c', {}, 'T', {})        % table 3: triangleFace structs
      FT  (:,1) double = zeros(0,1)                                            % face k is triangle FT(k), 0 if not
   end

   methods
      function obj = AlgCrv(Vname, E, EcQ, fN, fD, F, FC, FP, FB, EcA, EcI, Tri, FT)
      % objective: the subdivision, exactly as AlgCon takes it, plus the algebraic edge table.
      % [input]  as RatCon's constructor, then FP/FB (the implicit faces), EcA (the new edges) and EcI
      %          (the tie rows, table 2)
      %
      % THE CONSTRUCTOR PROTOCOL of RatCon applies unchanged: the LEAF writes the subdivision and then
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
         if nargin >= 11, obj.EcI = EcI; end
         if nargin >= 12 && ~isempty(Tri), obj.Tri = Tri(:); end
         if nargin >= 13, obj.FT = FT(:); else, obj.FT = zeros(numel(FC), 1); end
         obj = obj.setPinnedDefaults();
      end

      function [val, idx] = eval(obj, X)
      % objective: the function's value at points X, +Inf where no face owns the point.
      % [input]  X : n x 2
      % [output] val : n x 1;  idx : n x 1 face index, 0 where none owns it
      %
      % A HOLE IS +INFINITY, not a gap: outside `dom f**` the supremum diverges, and leaving the
      % cell out of the subdivision is how every producer here says so.
         if nargin ~= 2
            error('AlgCrv:nargin', 'expects (obj, X).');
         end
         if ~isnumeric(X) || size(X,2) ~= 2
            error('AlgCrv:points', 'X must be n x 2; got %s.', mat2str(size(X)));
         end
         n = size(X,1);
         val = inf(n,1);  idx = zeros(n,1);
         % FACE BY FACE over all points at once (2026-11-01): erratum Example 2's subdivision has
         % thousands of faces, and the point-by-face loop called algBranch once per pair. The rows
         % are read with the same tolerances and the max keeps the same (first) face on a tie.
         for k = 1:obj.nf
            M = obj.facePoints(k, X);
            if ~any(M), continue, end
            v = nan(n,1);  v(M) = obj.faceValues(k, X(M,:));
            up = M & ~isnan(v) & (isinf(val) | v > val);
            val(up) = v(up);  idx(up) = k;
         end
      end
   end

   methods (Access = private)
      function M = facePoints(obj, k, X)
      % faceOwns for every row of X at once.
         M = true(size(X,1), 1);
         rows = obj.FC{k};
         for t = 1:size(rows,1)
            if ~any(M), return, end
            P = X(M,:);  ok = true(size(P,1), 1);
            tbl = 0;
            if size(rows,2) >= 4, tbl = rows(t,4); end
            if tbl == 3
               T = obj.Tri(rows(t,1)).T;  i = rows(t,2);  j = mod(i,3) + 1;  o = 6 - i - j;
               e = T(j,:) - T(i,:);
               so = sign(e(1)*(T(o,2) - T(i,2)) - e(2)*(T(o,1) - T(i,1)));
               cr = (e(1)*(P(:,2) - T(i,2)) - e(2)*(P(:,1) - T(i,1))) * so;
               ok = ~(rows(t,3) * cr < -1e-9 * max(1, norm(e)) * max(1, sqrt(sum(P.^2, 2))));
            elseif tbl == 2
               r = obj.EcI(rows(t,1));
               v1 = algBranch(r.FP1, P, r.FB1);  v2 = algBranch(r.FP2, P, r.FB2);
               ok = isfinite(v1) & isfinite(v2) & ~(rows(t,2) * (v1 - v2) < -1e-9);
            elseif tbl == 1
               a = obj.EcA(rows(t,1));
               v = algBranch(a.FP, P, a.FB);
               q = ratQ.evalFace(a.num, a.den, P);
               ok = isfinite(v) & ~(rows(t,2) * (q - v) < -1e-9);
            else
               c = obj.EcQ(rows(t,1), :);
               br = 0;
               if size(rows,2) >= 3, br = rows(t,3); end
               if br ~= 0
                  v = ratQ.pairLineEval(c, br, P);
               else
                  v = c(1)*P(:,1).^2 + c(2)*P(:,1).*P(:,2) + c(3)*P(:,2).^2 + ...
                      c(4)*P(:,1) + c(5)*P(:,2) + c(6);
               end
               if rows(t,2) == 0
                  ok = ~(abs(v) > 1e-9);                     % an EQUALITY row: a THIN face
               else
                  ok = ~(rows(t,2) * v < -1e-9);
               end
            end
            M(M) = ok;
         end
      end

      function v = faceValues(obj, k, X)
      % faceValue for every row of X at once.
         if numel(obj.FT) >= k && obj.FT(k) > 0
            v = X * obj.Tri(obj.FT(k)).s(:) - obj.Tri(obj.FT(k)).c;
         elseif ~isempty(obj.FP) && numel(obj.FP) >= k && ~isempty(obj.FP{k})
            v = algBranch(obj.FP{k}, X, obj.FB(k));
         else
            v = ratQ.evalFace(obj.fN(k,:), obj.fD(k), X);
         end
      end

   end
end
