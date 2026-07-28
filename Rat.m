classdef (Abstract) Rat
   % Rat: the FUNCTION-type axis of a piecewise function -- rational (a polynomial numerator over a
   % linear denominator). One of the two axes every type in this family is built from; the other is
   % the SUBDIVISION-type axis, Par/Pol. See RatPar.m for the full picture.
   %
   % Rat is the general case: `den` may be any linear form. Its specialization Qua pins the
   % denominator to 1 (see Qua.m).
   %
   % A property-less ABSTRACT MARKER, deliberately. It carries no data and, for now, no methods:
   % all state lives once on RatPar, so that the concrete classes can inherit along BOTH axes
   % without MATLAB's multiple-inheritance member conflicts (a property defined in two superclasses
   % is fatal and unresolvable -- MATLAB:class:conflictingSuperClassProperty). What a marker buys
   % is the ability to ask about ONE axis independently of the other:
   %
   %       isa(f, 'Rat')   -- is the function rational?      (true for every type here)
   %       isa(f, 'Qua')   -- is the function a polynomial?  (denominator pinned to 1)
   %       isa(f, 'Par')   -- parabolic subdivision?         (true for every type here)
   %       isa(f, 'Pol')   -- polyhedral subdivision?        (conics pinned to 0)
   %
   % rather than having to enumerate combinations (`isa(f,'QuaPol') || isa(f,'QuaPar')`). RatPar's
   % `set.den`/`set.Ec` validators use exactly these tests to enforce each type's pinned values
   % from a single definition site.
   %
   % WHERE BEHAVIOUR SHOULD GO. When function-axis behaviour is factored out of the concrete
   % classes, it belongs here (on Rat) or on Qua: evaluating a face's formula, matrixForm,
   % evalHessian, scalarMul, negate, addQuadratic, degree -- everything that touches only the
   % coefficients `f`/`den` and never the mesh. Subdivision-axis behaviour (point location,
   % orderEdges, createDom, clipping) belongs on Par/Pol instead.
   %
   % NOT MODELLED AS A TYPE: the numerator's DEGREE. `f` is stored in the 10-wide cubic basis, so
   % a cubic numerator is representable, but nothing in this toolbox DISPATCHES on degree -- every
   % consumer either rejects a cubic outright (conjPieceCPLQ, quaPolToPlq, assertOperable) or
   % ignores the distinction (isConvex is the sole acceptor). A type should track what changes the
   % ALGORITHM, not what changes the data's shape; degree changes only a guard, so it stays a
   % runtime check rather than becoming a `Cub` class. Making it a type would double the lattice
   % (degree x denominator x subdivision = 8 combinations) for no dispatch benefit. Revisit if an
   % operator ever PRODUCES cubics rather than merely tolerating them.
end
