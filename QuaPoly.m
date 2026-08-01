classdef QuaPoly < QuaPol
   % QuaPoly is retained as a backward-compatible ALIAS of QuaPol (QUAdratic on a POLyhedral
   % subdivision). The implementation lives in QuaPol.m; new code should use QuaPol.
   % QuaPoly(V,E,f,F) and QuaPoly(f) keep working, and the static factories
   % (QuaPoly.energy, QuaPoly.oneNorm, QuaPoly.examples, ...) are inherited from QuaPol.
   %
   % WHY THIS EXISTS. QuaPol.m's rename note (2026-07-27) says no compatibility shim was left
   % "because CCA2 has no tagged release, so nothing external can depend on it". That premise was
   % wrong at the time it was written: the SCIP feasibility spike (AI/spike/SCIP) bridges to this
   % toolbox through MATLAB Engine, and its glue `src/cca2ConvexEnvelope.m` constructs
   % `QuaPoly(Vin, Ein, fin, Fin)`. The rename broke that bridge silently -- an
   % "Unrecognized function or variable 'QuaPoly'" at the first call, from a project that cannot
   % be edited from here and had every reason to expect the name it was written against.
   %
   % So the alias goes in rather than the caller being changed, for the same reason PLQVC's does:
   % a downstream consumer should not have to move because an internal naming scheme was
   % regularized. This one is cheap -- QuaPol adds no properties of its own, so the subclass is a
   % constructor forward and nothing else.
   %
   % Keep it through 0.1 at least. Removing it is a breaking change for a real, identified caller,
   % and should be done (if ever) by deprecation with that caller updated first.
   methods
       function obj = QuaPoly(varargin)
           obj@QuaPol(varargin{:});
       end
   end
end % classdef
