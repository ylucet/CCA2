classdef (Abstract) PAlgCon < PAlg & PConic
   properties
      f = []
      den = []
      FK = []          % the new algebraic-face data, pinned to empty on a Rat
   end
   methods
      function obj = set.FK(obj, val)
         if ~isempty(val) && isa(obj,'PRat')
            error('PAlgCon:FKMustBeEmpty', 'a %s is rational, so FK must be empty.', class(obj));
         end
         obj.FK = val;
      end
      function obj = set.den(obj, val)
         if ~isempty(val) && isa(obj,'PQua')
            trivial = repmat([0 0 1], size(val,1), 1);
            if any(val ~= trivial,'all')
               error('PAlgCon:denMustBeTrivial','%s: den must be [0 0 1].', class(obj));
            end
         end
         obj.den = val;
      end
      function k = kind(obj), k = class(obj); end
   end
end
