classdef (Abstract) Alg
   % Alg: the most general FUNCTION-type axis -- ALGebraic, i.e. a face need not be given by a
   % formula at all. It is cut out IMPLICITLY, by a polynomial relation p(x,y,z) = 0 with integer
   % coefficients, together with a BRANCH INDEX saying which of that polynomial's real roots in z
   % is the function's value. One of the two axes every type in this family is built from; the
   % other is the SUBDIVISION-type axis, Pol/Par/Conic. See AlgCon.m for the full picture.
   %
   % The function axis is a three-level chain, most special first:
   %
   %       Qua   <   Rat   <   Alg
   %       polynomial  numerator over    z is a named ROOT of an integer
   %       (den = 1)   a linear den      polynomial p(x,y,z), no formula
   %
   % WHY THIS LEVEL EXISTS. `biconj` of a QuaPol is the closed convex envelope, and where two
   % CONVEX pieces bridge a gradient jump the envelope's face is not rational. Derived in closed
   % form and checked against the primal definition to 1.5e-11 (DECISIONS.md 2026-09-05,
   % `.claude/bridgeFaceProbe.m`): for q1 = |v|^2/2 and q2 = |v-a|^2 + c the two pieces tie on a
   % CIRCLE and the envelope over it is
   %
   %       co f(x) = <2a, x-2a> + R*|x-2a| + 2c ,        R = 2*sqrt(|a|^2 - c) ,
   %
   % affine PLUS A NORM. The best affine fit to it is off by 2.9 and the best quadratic fit by
   % 0.36, so neither Qua nor Rat can hold it -- that is a measurement, not a modelling choice.
   %
   % AND WHY IT IS NOT A NEW NUMBER FIELD. `R` is irrational but `R^2 = 4(|a|^2-c)` is not, so
   % written IMPLICITLY the same face is a quadric cone with RATIONAL coefficients:
   %
   %       (z - <2a,x-2a> - 2c)^2 - R^2*|x-2a|^2 = 0 .
   %
   % "Not a polynomial in x" is not "not algebraic", and the difference is what makes this level
   % cheap: every stored coefficient stays an integer, so equality stays bitwise and every mesh
   % predicate stays an exact integer sign -- the properties the whole exact pipeline rests on.
   % The general (different-Hessian) case is bigger and still implicit and still rational:
   % measured degrees (x,y,z) from (2,2,2) isotropic to (6,6,4) anisotropic, integer coefficients
   % under 1e9, vanishing on the face to 1.4e-16 with a ROOT GAP OF EXACTLY ZERO -- the face value
   % IS a root of the stored polynomial (DECISIONS.md 2026-09-06).
   %
   % ONE FACE KIND, NOT FOUR. The 2026-08-23 plan listed `fkind` in {QUAD, RAT, RULED, AFFINE},
   % the AFFINE case being `<p,x> - f*(p)` with `p` a dual vertex of degree up to 4 -- irrational
   % coefficients unless the cell NAMES that vertex, which is what "a decoration of the QuaCon"
   % was for. The implicit form SUBSUMES it: `p` satisfies its own minimal polynomial over Q, so
   % eliminating `p` from `z = <p,x> - f*(p)` leaves a relation in (x,y,z) with rational
   % coefficients, of degree at most 4 in z. So an AFFINE face is an implicit face like any other
   % and needs no special case -- and the S4 vertex of doc/QuaConExample.md, where three pieces
   % tie at one dual point of degree 4, is stored the same way as the two-piece bridge that sweeps
   % a curve. `src` survives on AlgCon as PROVENANCE (which dual vertex a face came from, for the
   % classification table and for debugging), never as the definition of a value.
   %
   % WHAT IS *NOT* RELAXED HERE. Alg widens the storable FACE and nothing else. The subdivision
   % stays conic: the bridge REGION was measured to be the intersection of exactly two sign
   % conditions on two conics (`.claude/bridgeRegionSides.m`, `.claude/bridgeRegionCross.m`) and
   % the non-bridge regions are polyhedral, so no edge needs a shape `Conic` cannot already hold.
   % That is why the type is AlgCon and not `AlgAlg`: only the FUNCTION axis grows.
   %
   % A property-less ABSTRACT MARKER -- see Rat.m for why the axes are markers rather than
   % data-carrying classes. Rat < Alg, so `isa(f,'Alg')` is true for every existing type, which is
   % the honest reading: a polynomial IS an algebraic function. The test that separates them is
   % `isa(f,'Rat')`, and AlgCon's `set.FP` validator uses exactly that to pin the algebraic payload
   % to empty on every rational type, from one definition site.
end
