function h = infConv(f, g, engine)
% infConv  Infimal convolution (f box g)(x) = inf_z f(z) + g(x-z), via conj(add(conj f,conj g)).
%
% objective: compute f box g using the identity (f box g)* = f* + g*, hence
%   f box g = (f* + g*)* whenever f box g is itself convex ([JOGO]/[COAP]; DESIGN.md II.6). This
%   holds when f and g are BOTH CONVEX (a sum of conjugates is convex, so the final conj recovers
%   f box g exactly, not merely its convex envelope); NOT checked here -- see the caveat below.
%
% [input]  f, g   : QuaPoly | QuaPar, operable (degree<=2), BOTH CONVEX
%          engine : conjugate engine passed through to conj/add's downstream conj calls,
%                   'cplq' (default) | 'pqp' | 'graph'
% [output] h      : QuaPoly | QuaPar, the inf-convolution f box g
%
% CAVEAT (DESIGN.md II.6): valid (returns the TRUE inf-convolution) only for f,g both convex. For
%   nonconvex f or g this formula still runs but computes conj(conj f + conj g), which is the
%   CONVEX ENVELOPE of the actual inf-convolution, not the inf-convolution itself.
%
% METHOD: conj(f,engine) and conj(g,engine) may each come back as either QuaPoly (conjCPLQ's
%   full-domain-quadratic shortcut) or QuaPar (its general single-piece case) -- see conjCPLQ.m's
%   own header. QuaPoly.add and QuaPar.add only accept same-class operands, so both conjugates
%   are first promoted to QuaPar (toQuaPar below): QuaPoly and QuaPar share the same V/E/f/F
%   layout (QuaPar is a strict superset, all-zero Ec -- DESIGN.md II.2), so this promotion is a
%   lossless relabeling, not a geometric operation. addQuaPar then combines them and the result is
%   conjugated back through the same engine.
    if nargin < 3, engine = 'cplq'; end
    cf = toQuaPar(conj(f, engine));
    cg = toQuaPar(conj(g, engine));
    h = conj(cf.add(cg), engine);
end

% ================================================================================================
function q = toQuaPar(obj)
% Promote a QuaPoly to the equivalent QuaPar (all-zero Ec); a QuaPar is returned unchanged.
    if isa(obj, 'QuaPar'), q = obj; return; end
    if ~isa(obj, 'QuaPoly')
        error('infConv:unsupportedType', ...
            'infConv: conj(f,engine) returned a %s; expected QuaPoly or QuaPar.', class(obj));
    end
    if obj.nv == 0
        q = QuaPar(obj.f);
    else
        q = QuaPar(obj.V, obj.E, zeros(obj.ne, 6), obj.f, obj.F);
    end
end
