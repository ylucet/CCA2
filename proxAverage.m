function h = proxAverage(f, g, lambda, mu, engine)
% proxAverage  Proximal average P = P_mu(f,g;lambda,1-lambda), lambda in [0,1].
%
% objective: compute P via the derivation in DESIGN.md II.6 -- P is characterized by its own
%   Moreau envelope being the convex combination of f's and g's: e_mu*P = lambda*e_mu*f +
%   (1-lambda)*e_mu*g. Writing T_mu*h := mu*h + 1/2||.||^2 (moreau's own inner term) and
%   substituting moreau's expand-the-square identity into both sides, the (1/(2mu))||.||^2 terms
%   cancel (since lambda+(1-lambda)=1), giving (T_mu*P)* = lambda*(T_mu*f)* + (1-lambda)*(T_mu*g)*.
%   The right-hand side is always convex/proper/lsc (a conic combination of conjugates), so
%   PROVIDED T_mu*P is itself convex (guaranteed when f,g are convex and mu>0), Fenchel-Moreau
%   lets us recover T_mu*P by conjugating both sides; solving for P gives the code below: two
%   conjugations feed a weighted add, a third conjugation recovers T_mu*P, then P falls out after
%   subtracting back 1/2||.||^2 and rescaling by 1/mu.
%
% [input]  f, g   : QuaPoly | QuaPar, operable (degree<=2), BOTH CONVEX
%          lambda : real scalar in [0,1] (weight on f; g gets 1-lambda)
%          mu     : positive real scalar (the shared Moreau-envelope parameter)
%          engine : conjugate engine, 'cplq' (default) | 'pqp' | 'graph'
% [output] h      : QuaPoly | QuaPar, the proximal average P
%
% CAVEAT (DESIGN.md II.6): like infConv, the last step is a genuine biconjugation (not pure
%   algebra, unlike moreau/lasryLions), so this is valid (returns the actual proximal average)
%   only for f,g both convex; NOT checked here.
%
% METHOD: conj(gf,engine)/conj(gg,engine) may each come back as either QuaPoly or QuaPar (see
%   infConv.m's own header for why), and add only accepts same-class operands, so both are
%   promoted via toQuaPar.m before the weighted add -- same plumbing as infConv.
    if nargin < 5, engine = 'cplq'; end
    gf = addQuadratic(scalarMul(f, mu), eye(2), [0;0], 0);   % T_mu f := mu*f + 1/2||.||^2
    gg = addQuadratic(scalarMul(g, mu), eye(2), [0;0], 0);   % T_mu g := mu*g + 1/2||.||^2
    cf = toQuaPar(conj(gf, engine));
    cg = toQuaPar(conj(gg, engine));
    s  = scalarMul(cf, lambda).add(scalarMul(cg, 1-lambda));
    h  = addScaledEnergy(scalarMul(conj(s, engine), 1/mu), -1/(2*mu));
end
