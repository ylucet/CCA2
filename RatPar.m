classdef (Abstract) RatPar < RatCon & Par
   % RatPar: rational function on a PARABOLIC subdivision -- the level of the lattice at which
   % every edge conic is a parabola or a line (b^2 - 4ac = 0).
   %
   % IT DECLARES NOTHING. All mesh data (V/E/f/F/P/dom/nv/ne/nf/fIsConvex) and both axis-varying
   % properties (`den`, `Ec`) live once on RatCon -- see RatCon.m's "DATA LIVES HERE, ONCE" note
   % for why a property defined in two superclasses is fatal in MATLAB and why that forces a
   % single declaration site. This class exists to carry the `Par` marker, nothing else, and it is
   % what every pre-existing concrete type (RatPol, QuaPar, QuaPol, QuaParCPLQ) inherits.
   %
   % WHY IT IS NO LONGER THE TOP. The subdivision axis grew a level ABOVE it on 2026-08-24:
   %
   %       Pol  <  Par  <  Conic
   %
   % because `conj` provably produces edges that are not parabolas. An edge of f* between two
   % faces g_i, g_j lies on {g_i = g_j}, whose quadratic part is (H_i - H_j)/2 and whose
   % discriminant is therefore -det(H_i - H_j). [COAP]/[JOGO] Theorem 6 makes that vanish by
   % comparing only ADJACENT pieces of f (CONJ_FIELD_PROOF.md Theorem 3), but Step 3 maximises
   % over ALL pairs, so from three pieces up some compared pair is non-adjacent and the difference
   % conic is a genuine ellipse. CONJ_FIELD_PROOF.md 7.4b traces such an arc, of positive length,
   % against the exact definition of the conjugate. Conic.m and RatCon.m carry the detail.
   %
   % NOTHING BELOW THIS LINE OF THE LATTICE CHANGED. Adding Conic above RatPar is a pure widening:
   % every object that was a RatPar is still a RatPar and still a Par, `QuaPar.assertParabolicEdges`
   % still runs in QuaPar's constructor, and no existing type can hold a conic it could not hold
   % before. The new capability is reachable only through the new type QuaCon, which is a RatCon
   % and deliberately NOT a RatPar.
   %
   % SEE ALSO: RatCon.m (the data, the validators, the constructor protocol, and the full lattice
   % picture), Conic.m / Par.m / Pol.m (the subdivision axis), Rat.m / Qua.m (the function axis),
   % RETURN_TYPE.md (why one static return type was needed at all).
end
