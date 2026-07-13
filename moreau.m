function h = moreau(f, mu, engine)
% moreau  Moreau envelope e_mu f(x) = inf_z f(z) + (1/(2*mu))||x-z||^2.
%
% objective: compute e_mu*f as a SINGLE conjugate via the "expand the square" identity
%   [HIRIART-URRUTY-07] (DESIGN.md II.6):
%     e_mu*f(x) = inf_z [f(z) + (1/(2*mu))||x-z||^2]
%               = (1/(2*mu))||x||^2 - (1/mu) * sup_z [x.z - (mu*f(z) + 1/2||z||^2)]
%               = (1/(2*mu))||x||^2 - (1/mu) * (mu*f + 1/2||.||^2)*(x)
%   Deliberately NOT infConv(f, (1/(2*mu))||.||^2): that route is a BICONJUGATE (two conj calls)
%   and per infConv.m's own caveat is only valid for f convex. This formula is pure algebraic
%   regrouping of the definition of conjugate -- no biconjugation, no convexity assumption on f
%   at all -- which is why it calls conj only ONCE.
%
% [input]  f      : QuaPoly | QuaPar, operable (degree<=2); need NOT be convex
%          mu     : positive real scalar (the envelope parameter)
%          engine : conjugate engine, 'cplq' (default) | 'pqp' | 'graph'
% [output] h      : QuaPoly | QuaPar, the Moreau envelope e_mu*f
    if nargin < 3, engine = 'cplq'; end
    g = addQuadratic(scalarMul(f, mu), eye(2), [0;0], 0);       % T_mu f = mu*f + 1/2||.||^2
    h = addScaledEnergy(scalarMul(conj(g, engine), -1/mu), 1/(2*mu));
end
