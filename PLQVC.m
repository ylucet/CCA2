classdef PLQVC < QuaPol
   % PLQVC is retained as a backward-compatible ALIAS of QuaPol (quadratic on a polyhedral
   % subdivision). The implementation now lives in QuaPol.m; new code should use QuaPol.
   % PLQVC(V,E,f,F) and PLQVC(f) keep working, and the released static factories
   % (PLQVC.oneNorm, PLQVC.energy, PLQVC.examples, ...) are inherited from QuaPol.
   % See QuaPol.m and DESIGN.md.
   methods
       function obj = PLQVC(varargin)
           obj@QuaPol(varargin{:});
       end
   end
end % classdef
