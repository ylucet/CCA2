classdef (Abstract) Qua < Rat
   % Qua: the QUAdratic specialization of the function-type axis -- a polynomial, i.e. a Rat whose
   % denominator is pinned to 1 (`den` is [0 0 1] on every face).
   %
   % Qua < Rat because a polynomial IS a rational function with unit denominator. That subsumption
   % propagates for free: any object marked Qua also answers true to isa(.,'Rat').
   %
   % A property-less ABSTRACT MARKER -- see Rat.m for why the axes are markers rather than
   % data-carrying classes, and for where function-axis behaviour should live.
   %
   % The pinning is ENFORCED, not merely documented: RatPar's `set.den` rejects any assignment of a
   % non-trivial denominator to an object that isa(.,'Qua'), from the single place `den` is
   % defined. So a QuaPar or QuaPol cannot be made to lie about its own type.
   %
   % Note "Qua" names the DENOMINATOR being trivial, not the numerator's degree -- the numerator is
   % stored in the cubic basis and degree is a runtime check, not a type (see Rat.m's closing note).
end
