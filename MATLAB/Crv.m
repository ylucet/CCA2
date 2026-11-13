classdef (Abstract) Crv < Conic
   % Crv: the top of the SUBDIVISION-type axis -- an edge may be an arbitrary ALGEBRAIC CURVE, not
   % just a conic. One of the two axes every type in this family is built from; the other is the
   % FUNCTION-type axis, Qua/Rat/Alg. See RatCon.m for the full picture and AlgCon.m for where the
   % data lives.
   %
   % The subdivision axis is now a four-level chain, most special first:
   %
   %       Pol   <   Par    <   Con    <   Crv
   %       every     every      an edge     an edge may be a curve of ANY degree, named
   %       edge      edge is    may be      as (rational polynomial, BRANCH) because the
   %       straight  a          ANY conic   curve that is rational names a PAIR of branches
   %                 parabola               and a face needs one of them
   %
   % WHY THIS LEVEL EXISTS, and it is not a limitation being relaxed casually -- Conic.m makes the
   % same argument one level down, and this one has the same shape.
   %
   % The CONVEX ENVELOPE of a piecewise quadratic has faces a conjugate does not: where the two
   % pieces' tie curve is CURVED, the supremum is attained ON that tie, and the resulting BRIDGE
   % face carries no formula at all -- it is the branch of an implicit polynomial, which is what the
   % `Alg` level of the FUNCTION axis exists for. That much AlgCon already held.
   %
   % What it could not hold is the EDGE between such a face and an ordinary quadratic one. Measured
   % (DECISIONS.md 2026-09-19): a bridge face is strictly above the best quadratic face owning the
   % same point at 60 of 60 sampled points, so the two really do compete, and on the c = -2 fixture
   % they compete inside a cell whose other winner belongs to a DIFFERENT cell pair. The boundary is
   %
   %       { x : FP(x, y, q(x,y)) = 0 } ,   q the quadratic face's value,
   %
   % which has degree 4 in x and y on the real fixtures -- a curve, not a conic. So the subdivision
   % axis had to grow, exactly as the function axis did for the value.
   %
   % THE EDGE IS NAMED, NOT STORED AS COORDINATES, which is the same move every irrational object in
   % this toolbox gets. `FP(x, y, q)` vanishes wherever q meets EITHER root of FP, so the rational
   % polynomial names the PAIR of branches: measured, its sign agrees with the product over both
   % branches everywhere the fibre is real, and with one branch only 63% to 97% of the time. The
   % edge therefore carries the polynomial AND a branch index, and the side of it at a rational point
   % is decided by COUNTING how many roots of the fibre lie at or below q -- Sturm, over the
   % integers, with no root ever computed (`algCurveSide`). `Vname` names a vertex, `FB` names a
   % value, `QuaCon.FC`'s third column names one line of a degenerate conic; this is that same
   % convention one degree up.
   %
   % WHAT A Crv PINS. Nothing. It is the general level, so `EcA` is unconstrained here and pinned
   % EMPTY on every type that is not a Crv -- enforced by AlgCon's `set.EcA`, which reads the TRAIT
   % rather than the class name so that it covers types that do not exist yet.
end
