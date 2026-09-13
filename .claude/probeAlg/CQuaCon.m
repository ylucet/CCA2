classdef CQuaCon < CRatCon & PQua
   methods
      function obj = CQuaCon(f)
         if nargin == 0, return, end
         obj.f = f;
      end
   end
end
