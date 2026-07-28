classdef (Abstract) Par
   % Par: the SUBDIVISION-type axis of a piecewise function -- PARabolic, i.e. each edge may be a
   % parabola arc (per-edge conic `Ec`). One of the two axes every type in this family is built
   % from; the other is the FUNCTION-type axis, Rat/Qua. See RatPar.m for the full picture.
   %
   % Par is the general case: `Ec` may carry a genuine conic on any edge. Its specialization Pol
   % pins every conic to zero, making every edge straight (see Pol.m).
   %
   % Only PARABOLAS, never ellipses or hyperbolas -- and that is a property of the mathematics, not
   % a limitation of the storage. CCA2's goal is the conjugate/biconjugate of a QuaPol, and the
   % subdivisions arising from it are very special: a parabolic edge only ever occurs surrounded by
   % two parallel rays, and hyperbolic edges never arise at all (see [COAP]/[JOGO],
   % SUPPORT_MATRIX.md section 0, and QuaPar.assertParabolicEdges which checks b^2-4ac=0).
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
