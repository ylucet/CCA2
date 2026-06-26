classdef PLQVC < QuaPoly
   % PLQVC is retained as a backward-compatible ALIAS of QuaPoly (quadratic on a polyhedral
   % subdivision). The implementation now lives in QuaPoly.m; new code should use QuaPoly.
   % PLQVC(V,E,f,F) and PLQVC(f) keep working, and the released static factories
   % (PLQVC.oneNorm, PLQVC.energy, PLQVC.examples, ...) are inherited from QuaPoly.
   % See QuaPoly.m and DESIGN.md.
   methods
       function obj = PLQVC(varargin)
           obj@QuaPoly(varargin{:});
       end
   end
end % classdef
