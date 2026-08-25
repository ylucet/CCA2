classdef (Abstract) Par < Conic
   % Par: the SUBDIVISION-type axis of a piecewise function -- PARabolic, i.e. each edge may be a
   % parabola arc (per-edge conic `Ec`). One of the two axes every type in this family is built
   % from; the other is the FUNCTION-type axis, Rat/Qua. See RatPar.m for the full picture.
   %
   % Par is the MIDDLE of a three-level chain, Pol < Par < Conic: `Ec` may carry a parabola on any
   % edge but nothing of higher discriminant. Its specialization Pol pins every conic to zero,
   % making every edge straight (see Pol.m); its generalization Con drops the discriminant
   % condition entirely (see Conic.m).
   %
   % Only PARABOLAS at this level, never ellipses or hyperbolas -- enforced by
   % QuaPar.assertParabolicEdges, which checks b^2-4ac=0.
   %
   % CORRECTED 2026-08-24. This header used to say that was "a property of the mathematics, not a
   % limitation of the storage", on the grounds that a parabolic edge only ever occurs surrounded
   % by two parallel rays and hyperbolic edges never arise. That claim holds for every comparison
   % [COAP]/[JOGO] Theorem 6's proof actually covers -- two ADJACENT pieces of f, whose difference
   % conic is degenerate (CONJ_FIELD_PROOF.md Theorem 3) -- and fails in general, because Step 3
   % maximises over ALL pairs of pieces. From three pieces up, some compared pair is non-adjacent
   % and {g_i = g_j} is a genuine ellipse; CONJ_FIELD_PROOF.md 7.4b traces such an arc of positive
   % length against the exact definition of the conjugate. Those functions are stored as QuaCon,
   % which is a RatCon and NOT a RatPar. Par remains exactly as strict as it always was.
   %
   % A property-less ABSTRACT MARKER -- see Rat.m for why the axes are markers rather than
   % data-carrying classes.
   %
   % WHERE BEHAVIOUR SHOULD GO. When subdivision-axis behaviour is factored out of the concrete
   % classes, it belongs here (on Par) or on Pol: point location, createP/orderEdges, edgeConics,
   % createDom, isDomBounded, plotDomain, and the clipping primitives (facePoly, polyConstraints,
   % clipByFace, clipPolyHalfPlane/clipArcByHalfPlane, subdivision overlay) -- everything that
   % touches only V/E/Ec/F/P and never the face coefficients. Function-axis behaviour belongs on
   % Rat/Qua instead.
end
