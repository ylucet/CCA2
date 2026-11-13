function q = toQuaPar(obj)
% toQuaPar  Promote a QuaPol to the equivalent QuaPar (all-zero Ec); a QuaPar or QuaParCPLQ is
% returned unchanged.
%
% objective: QuaPol.add/QuaPar.add (and other same-class-only operators) can't combine a QuaPol
%   with a QuaPar directly, even though QuaPol and QuaPar share the same V/E/f/F layout (QuaPar
%   is a strict superset, DESIGN.md II.2) -- this arises whenever conj(f,engine) may come back as
%   either type (conjCPLQ's full-domain-quadratic shortcut returns QuaPol; its general
%   single-piece case returns QuaPar; its general multi-face/non-triangular Case C returns
%   QuaParCPLQ -- see conjCPLQ.m's own header) and the caller needs to add two such conjugates
%   together (infConv.m, proxAverage.m). This promotion is a lossless relabeling, not a geometric
%   operation; QuaParCPLQ is passed through as-is since it already implements its own add/
%   scalarMul/conj (see QuaParCPLQ.m) rather than being convertible to a true QuaPar's V/E/Ec/F/P.
%
% [input]  obj : QuaPol | QuaPar | QuaParCPLQ
% [output] q   : QuaPar (obj unchanged if it already is one), or QuaParCPLQ unchanged
    if isa(obj, 'QuaPar'), q = obj; return; end
    % QuaParCPLQ (Case C's own conjugate wrapper, see QuaParCPLQ.m) already provides add/
    % scalarMul/conj/etc. itself -- pass it through unchanged rather than erroring, so
    % infConv/proxAverage's toQuaPar(conj(f,engine)) call works when conj(f,engine) is Case C.
    if isa(obj, 'QuaParCPLQ'), q = obj; return; end
    if ~isa(obj, 'QuaPol')
        error('toQuaPar:unsupportedType', 'toQuaPar: expected a QuaPol or QuaPar, got %s.', class(obj));
    end
    if obj.nv == 0
        q = QuaPar(obj.f);
    else
        q = QuaPar(obj.V, obj.E, zeros(obj.ne, 6), obj.f, obj.F);
    end
end
