classdef (Abstract) Pol < Par
   % Pol: the POLyhedral specialization of the subdivision-type axis -- every edge is a straight
   % line, ray or segment, i.e. a Par whose per-edge conics are all pinned to zero.
   %
   % Pol < Par because a polyhedral subdivision IS a parabolic one with no curvature. That
   % subsumption propagates for free: any object marked Pol also answers true to isa(.,'Par').
   %
   % A property-less ABSTRACT MARKER -- see Rat.m for why the axes are markers, and Par.m for where
   % subdivision-axis behaviour should live.
   %
   % The pinning is ENFORCED, not merely documented: RatPar's `set.Ec` rejects any assignment of a
   % nonzero conic to an object that isa(.,'Pol'), from the single place `Ec` is defined. So a
   % RatPol or QuaPol cannot be made to lie about its own type.
   %
   % NAMING CAVEAT: "Pol" here abbreviates POLYHEDRAL (the subdivision), NOT polynomial. The
   % adjacent axis name Rat/Qua invites exactly that misreading, since the natural opposite of
   % "rational" is "polynomial" -- but that distinction is the Rat/Qua axis, and this one is purely
   % about whether edges are straight.
end
