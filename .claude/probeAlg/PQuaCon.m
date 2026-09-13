classdef PQuaCon < PRatCon & PQua
   methods
      function obj = PQuaCon(f)
         if nargin == 0, return, end
         obj.f = f;  obj.den = repmat([0 0 1], size(f,1), 1);
      end
   end
end
