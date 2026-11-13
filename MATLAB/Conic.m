classdef (Abstract) Conic
   % Conic: the most general SUBDIVISION-type axis -- CONic, i.e. each edge may be an arbitrary
   % conic (per-edge `Ec` with no constraint on its discriminant). One of the two axes every type
   % in this family is built from; the other is the FUNCTION-type axis, Rat/Qua. See RatCon.m for
   % the full picture.
   %
   % The subdivision axis is a three-level chain, most special first:
   %
   %       Pol  <  Par  <  Con
   %       every edge     every edge is a      an edge may be ANY conic
   %       straight       parabola or a line   (ellipse, hyperbola included)
   %
   % WHY THIS LEVEL EXISTS, AND WHY IT IS NOT A LIMITATION BEING RELAXED CASUALLY. Par.m records
   % -- correctly, for everything [COAP]/[JOGO] Theorem 6 actually covers -- that a conjugate's
   % edges are parabolas, "and that is a property of the mathematics, not a limitation of the
   % storage". That claim is TRUE for every comparison the theorem's proof covers (two ADJACENT
   % pieces of f, whose difference conic is degenerate: CONJ_FIELD_PROOF.md Theorem 3) and FALSE
   % in general. Step 3 of the conjugate algorithm maximises over ALL pairs of pieces, so once f
   % has three or more pieces some compared pair is non-adjacent, and then both compared functions
   % can be elliptic quadratics with a non-degenerate difference.
   %
   % That is not speculative: CONJ_FIELD_PROOF.md 7.4b exhibits a continuous five-piece PLQ whose
   % f* has an arc of POSITIVE LENGTH lying on the ellipse {g1 = g3}, traced against the exact
   % per-piece definition of the conjugate rather than against any CCA2 code, every traced point
   % on the conic to <= 2e-15. doc/QuaConExample.md 2 reduces it to THREE pieces. Such an f* is
   % not storable as a QuaPar, so `conj` needs this level to hold its own answers.
   %
   % A property-less ABSTRACT MARKER -- see Rat.m for why the axes are markers rather than
   % data-carrying classes, and Par.m for where subdivision-axis behaviour should live.
   %
   % WHY THE NAME IS SPELLED OUT WHERE THE OTHER AXIS MARKERS ARE ABBREVIATED (Pol, Par, Rat, Qua).
   % `CON` is a RESERVED DEVICE NAME on Windows, for every extension, so `Con.m` is the console and
   % not a file: it can be created and listed by a POSIX layer, but any Windows API that opens it
   % BLOCKS reading the console instead of returning the source. MATLAB reports this as
   %     Previously accessible file "...\Con.m" is now inaccessible.
   % Measured, not guessed -- a plain read of the path hung until it was killed. The compound class
   % names RatCon and QuaCon are ordinary filenames and keep the three-letter abbreviation.
   %
   % WHAT IS *NOT* RELAXED HERE. Con widens the storable EDGE SHAPE and nothing else. It does not
   % license a caller to hand an ellipse to QuaPar (whose constructor still asserts b^2-4ac=0),
   % and it does not change any existing type's behaviour: Par < Con, so every object that was a
   % Par before is a Con now, and none of them can carry a conic it could not carry before.
end
