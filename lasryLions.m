function h = lasryLions(f, lambda, mu, engine)
% lasryLions  Lasry-Lions double envelope h = -e_mu( -e_lambda(f) ), a pure composition of moreau.
%
% objective: regularize f by an inner Moreau envelope (parameter lambda), negate, take an outer
%   Moreau envelope (parameter mu), negate again (DESIGN.md II.6). No new geometry: moreau itself
%   needs no convexity assumption on its input ([HIRIART-URRUTY-07]), so lasryLions inherits that
%   -- f need not be convex, only lambda,mu small enough that each moreau call's own
%   positive-definiteness condition holds (see moreau.m).
%
% [input]  f      : QuaPoly | QuaPar, operable (degree<=2); need NOT be convex
%          lambda : positive real scalar (inner envelope parameter)
%          mu     : positive real scalar (outer envelope parameter)
%          engine : conjugate engine, 'cplq' (default) | 'pqp' | 'graph'
% [output] h      : QuaPoly | QuaPar, the Lasry-Lions double envelope
    if nargin < 4, engine = 'cplq'; end
    h = negate(moreau(negate(moreau(f, lambda, engine)), mu, engine));
end
