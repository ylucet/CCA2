function p = ratPolToPlq(obj)
% ratPolToPlq  Convert a CCA2 RatPol (Step 1's convex-envelope output, faces possibly RATIONAL)
%   into a cPLQ `plq` object, so cPLQ's own symbolic Step 2/3 (plq.maximum -> conjugate +
%   maximumConjugate) can conjugate it. The rational-face counterpart of quaPolToPlq.m.
%
% [input]  obj : RatPol, bounded faces, quadratic numerator (cubic rejected). Each face k carries
%                the numerator obj.f(k,5:10) in the stored [x^2 xy y^2 x y 1] convention (so the
%                polynomial is 1/2(c1 x^2 + 2 c2 xy + c3 y^2) + c4 x + c5 y + c6, matching
%                QuaPol.matrixForm) over the linear denominator obj.den(k,:) = [a b c] meaning
%                a x + b y + c.
% [output] p   : plq, one plq_1p piece per face (domain = the face's own polygon, function = the
%                face's rational expression as a sym in x,y).
%
% WHY THIS EXISTS: CCA2's Step 2 (conjPieceCPLQ) has no rational-piece branch, but cPLQ's does --
% that is the whole point of starting from cPLQ. Step 1's envelope of an indefinite triangle
% contains rational faces whenever it has to split (2- or 3-convex-edge), so handing those faces
% to cPLQ is what lets conjCPLQ finish the single-triangle conjugate. See conjCPLQ.m's
% conjEnvelopeViaCPLQ and SUPPORT_MATRIX.md section 1.2.
%
% COEFFICIENT CLEANING: CCA2's Step 1 is double precision (the design's exact rational arithmetic
% is a Phase 2 follow-up -- see conjCPLQ.m's NOTE on arithmetic), while cPLQ downstream is exact
% symbolic. Feeding raw doubles across that boundary turns rounding dust into exact nonzero
% symbolic terms: a coefficient of 2.6e-14 sitting next to one of 2.0 becomes a genuine xy term
% that MuPAD then carries through numden/solve, and the 0/0 structure these envelope formulas have
% AT a face vertex (the denominator of COAP eq.16 vanishes there, and so must the numerator) stops
% cancelling. Coefficients below RELTOL times the largest magnitude in the same row are therefore
% dropped before the hand-off. This is a numeric->symbolic hand-off convention, not a modelling
% choice: it only ever removes terms that double precision cannot distinguish from zero. It is not
% a substitute for doing Step 1 exactly.
    obj.assertOperable();
    RELTOL = 1e-12;
    x = sym('x'); y = sym('y');
    pieces = plq_1p.empty();
    for k = 1:obj.nf
        iVs = faceVertexIndices(obj, k);
        V = obj.V(iVs, :);
        if any(abs(obj.f(k,1:4)) > sqrt(eps))
            error('ratPolToPlq:cubic', 'Face %d is cubic; ratPolToPlq requires degree<=2.', k);
        end
        c  = dropDust(obj.f(k,5:10), RELTOL);
        dn = dropDust(obj.den(k,:),  RELTOL);
        num = 0.5*(c(1)*x^2 + 2*c(2)*x*y + c(3)*y^2) + c(4)*x + c(5)*y + c(6);
        den = dn(1)*x + dn(2)*y + dn(3);
        if den == 0
            error('ratPolToPlq:zeroDenominator', 'Face %d has an identically zero denominator.', k);
        end
        pieces(k) = plq_1p(domain(V, x, y), symbolicFunction(simplifyFraction(num/den)));
    end
    p = plq(pieces);
end

function v = dropDust(v, reltol)
% Zero out entries that double precision cannot distinguish from zero relative to the row's own
% scale. A row that is entirely zero is left alone.
    m = max(abs(v));
    if m > 0
        v(abs(v) < reltol*m) = 0;
    end
end

function iVs = faceVertexIndices(obj, k)
% Vertex indices around face k, in the order of its ordered edge list obj.P{k} (same convention as
% quaPolToPlq.m's/convEnvCPLQ.m's/conjCPLQ.m's own file-local helper of the same name, duplicated
% here since it is file-local in each of those too).
    face = obj.P{k}; iVs = zeros(1, numel(face));
    for i = 1:numel(face)
        j = face(i);
        if j > 0, iVs(i) = obj.E(j,1); else, iVs(i) = obj.E(-j,2); end
    end
end
