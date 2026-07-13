function q = toQuaPar(obj)
% toQuaPar  Promote a QuaPoly to the equivalent QuaPar (all-zero Ec); a QuaPar is returned
% unchanged.
%
% objective: QuaPoly.add/QuaPar.add (and other same-class-only operators) can't combine a QuaPoly
%   with a QuaPar directly, even though QuaPoly and QuaPar share the same V/E/f/F layout (QuaPar
%   is a strict superset, DESIGN.md II.2) -- this arises whenever conj(f,engine) may come back as
%   either type (conjCPLQ's full-domain-quadratic shortcut returns QuaPoly; its general
%   single-piece case returns QuaPar -- see conjCPLQ.m's own header) and the caller needs to add
%   two such conjugates together (infConv.m, proxAverage.m). This promotion is a lossless
%   relabeling, not a geometric operation.
%
% [input]  obj : QuaPoly | QuaPar
% [output] q   : QuaPar (obj unchanged if it already is one)
    if isa(obj, 'QuaPar'), q = obj; return; end
    if ~isa(obj, 'QuaPoly')
        error('toQuaPar:unsupportedType', 'toQuaPar: expected a QuaPoly or QuaPar, got %s.', class(obj));
    end
    if obj.nv == 0
        q = QuaPar(obj.f);
    else
        q = QuaPar(obj.V, obj.E, zeros(obj.ne, 6), obj.f, obj.F);
    end
end
