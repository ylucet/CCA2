classdef (Abstract) RatCon < AlgCon & Rat
   % RatCon: rational function on a CONIC subdivision -- the level of the lattice at which every
   % face still has a FORMULA.
   %
   % IT NO LONGER DECLARES THE DATA. Until 2026-09-13 this class was the lattice's root and the one
   % place the mesh properties lived. The function axis then gained a top level (Qua < Rat < Alg,
   % see Alg.m), so the root moved UP to AlgCon and the properties, their `set.` validators, the
   % edge-flag statics and `kind`/`isMeshed`/`setPinnedDefaults` moved with it -- unchanged, and
   % still declared exactly once. Nothing else about this class changed, and no type below it
   % changed behaviour: everything that was a RatCon still is one, and is additionally an Alg,
   % which is the honest reading (a rational function IS algebraic).
   %
   % WHAT THIS LEVEL NOW MEANS, and it is the useful half of the split. `isa(f,'Rat')` is the test
   % that separates "every face has a formula" from "a face may be a named ROOT of a polynomial",
   % and AlgCon's `set.FP` validator uses exactly that to pin the algebraic payload to empty here
   % and everywhere below. So this class is what makes the pinning expressible.
   %
   % The rest of this header is the original and still applies to the lattice as a whole.
   %
   % PURPOSE -- one static return type. Before this class existed, `conj` returned a DIFFERENT
   % class depending on the shape of its input (QuaPol for a full-domain strictly convex
   % quadratic, QuaPar for a single bounded triangle, QuaParCPLQ for a general bounded domain), so
   % no caller could write `g = f.conj()` without dispatching on class first. Every one of those
   % is now a RatPar, so `conj`/`biconj`/`convEnv`/... are statically typed: they return a RatPar,
   % and `kind(g)` says what it actually is. See RETURN_TYPE.md.
   %
   % TWO AXES, NOT ONE LIST. A piecewise function is built from two independent choices: the
   % FUNCTION type (rational, or its polynomial specialization) and the SUBDIVISION type (parabolic,
   % or its polyhedral specialization). Those axes are modelled as property-less abstract markers --
   % Rat with Qua < Rat, and Par with Pol < Par -- so either axis can be queried on its own:
   %
   %       isa(f,'Qua')   the denominator is pinned to 1
   %       isa(f,'Pol')   every edge conic is pinned to 0
   %
   % The four named types are the four cells of that 2x2 grid, and the inheritance is the grid's own
   % partial order -- a genuine diamond, because that is what a product lattice IS:
   %
   %                 RatCon  (rational on conic)        < Rat & Conic   <- this class
   %                    |
   %                 RatPar  (rational on parabolic)     < RatCon & Par
   %                 /    \
   %      (rational on    (quadratic on
   %       polyhedral)     parabolic)
   %          RatPol        QuaPar                       < RatPar & Pol / RatPar & Qua
   %                 \    /
   %                 QuaPol  (quadratic on polyhedral)   < RatPol & QuaPar
   %
   %       QuaCon      quadratic on a CONIC subdivision              < RatCon & Qua
   %       QuaParCPLQ  a QuaPar still in cPLQ's own symbolic form    < RatPar & Qua
   %
   % The subdivision axis has THREE levels, Pol < Par < Conic, because a conjugate's edges are only
   % parabolic when the two pieces compared are adjacent (CONJ_FIELD_PROOF.md Theorem 3); Step 3
   % compares non-adjacent pieces too, and those give genuine ellipses. See Conic.m. Con was added
   % ABOVE the previous top, so no existing type's behaviour changed: RatPar and everything under
   % it still is-a Par and still cannot hold a non-parabolic conic.
   %
   % QuaPol IS-A RatPol (a polynomial is a rational with unit denominator) and IS-A QuaPar (a
   % polyhedral subdivision is a parabolic one with no curvature), so both edges of the diamond are
   % mathematically honest and `isa` reports the full lattice.
   %
   % DATA LIVES HERE, ONCE. `den` and `Ec` are declared in this class alone -- NOT overridden in the
   % children -- because MATLAB makes a property defined in two superclasses fatal AND unresolvable
   % (MATLAB:class:conflictingSuperClassProperty; a child cannot redefine it either:
   % MATLAB:class:RedefinedProperty). Had RatPol declared `den` and QuaPar declared its own pinned
   % `den`, QuaPol could not exist. Instead each type's pinned values are enforced by the
   % `set.den`/`set.Ec` validators below, which read the object's own TRAITS -- so one definition
   % site serves the whole lattice and no type can be made to lie about itself.
   %
   % QuaParCPLQ deliberately does NOT inherit from QuaPar: it carries the same mathematical type
   % (marked Qua, and Par via this class) but not the mesh, so making `isa(x,'QuaPar')` true for it
   % would invite callers to read V/E/F that it does not have. `isMeshed()` is the intended test.
   %
   % NO STORED TYPE FLAG, BY DESIGN. `kind()` is a METHOD returning the object's actual class name,
   % not a property set at construction: a stored flag can drift out of sync with the real class
   % (after a copy, a subclass change, or a constructor path that forgets to set it), and a flag
   % that disagrees with the object is worse than no flag. Deriving it from `class(obj)` makes that
   % impossible. It also keeps MATLAB's own dispatch as the single source of truth -- which the
   % operator code already depends on heavily (e.g. maxQuaPar/infConv/moreau/proxAverage compose
   % with QuaParCPLQ only because MATLAB dispatches on the operand's actual class).
   %
   % VALUE SEMANTICS. RatPar is a VALUE class (not a handle), matching all four children: every
   % operator returns a new object and never mutates its argument. Do not change this without
   % auditing every operator, which all rely on copy-on-assign.
   %
   % SCOPE. Introducing RatPar does NOT mean implementing every possible RatPar. CCA2's goal is
   % QuaPol conjugate/biconjugate; the RatPols and QuaPars that actually arise from it are very
   % special cases (e.g. a parabolic edge only ever occurs surrounded by two parallel rays, and
   % hyperbolic edges never arise at all -- see [COAP]/[JOGO] and SUPPORT_MATRIX.md section 0).
   % RatPar exists here as a COMMON TYPE, not as a general rational-cubic-on-parabolic engine.
end
