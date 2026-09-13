classdef PBiconj < PAlgCon
   methods
      function obj = PBiconj(f, FK)
         if nargin == 0, return, end
         obj.f = f;  obj.FK = FK;
      end
   end
end
