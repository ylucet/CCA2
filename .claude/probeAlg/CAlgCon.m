classdef CAlgCon < PAlg & PConic
   properties
      f = []
      FP = []
   end
   methods
      function obj = CAlgCon(f, FP)
         if nargin == 0, return, end          % leaf-only protocol: no-arg writes NOTHING
         obj.f = f;
         if nargin > 1, obj.FP = FP; end
      end
      function obj = set.FP(obj, val)
         if ~isempty(val) && isa(obj,'PRat')
            error('CAlgCon:FPMustBeEmpty','a %s is rational, so FP must be empty.', class(obj));
         end
         obj.FP = val;
      end
   end
end
