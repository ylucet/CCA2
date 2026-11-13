function d = transformDomain(d0, Minv, vars)
% transformDomain  Push a `domain` through the linear map z = Minv*x, i.e. x = M*z.
%
% [input]  d0    : domain (its polygon's facets must be affine).
%          Minv  : 2x2, the map taking primal points to the new frame.
%          vars  : 1x2 sym, the variable names to build the result with.
% [output] d     : domain over {z : x = M z lies in d0}.
%
% Facets, not vertices. A constraint g(x) <= 0 becomes g(Mz) <= 0, which is a substitution and
% is exactly as valid for a RAY as for a segment -- so an unbounded face survives the change of
% frame with no special handling, and the +-intmax direction vertices are re-derived by region
% in the new frame rather than being transformed (multiplying intmax by M would be meaningless).
% This is the same reason quaPolToPlq builds faces from half-planes.
    M = inv(Minv);
    x = vars(1); y = vars(2);
    sub = [M(1,1)*x + M(1,2)*y, M(2,1)*x + M(2,2)*y];

    r0 = d0.polygon;
    g = sym.empty(1,0);
    for j = 1:size(r0.ineqs,2)
        g(j) = expand(subs(r0.ineqs(j).f, r0.vars, sub));
    end
    d = domain();
    d = d.domainEdge(g, vars);
end
