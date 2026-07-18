function [val, idx] = evalFunctionNDomain(fnd, s)
% evalFunctionNDomain  Numerically evaluate a cPLQ functionNDomain array (a list of (region,
%   symbolicFunction) pairs -- e.g. a plq's .maxConjugate after .maximum, or one piece's own
%   .conjugates after .conjugate) at a numeric dual point s, without converting the symbolic
%   result back into a CCA2 QuaPar. Part of the Phase 1 cPLQ integration (DESIGN.md II.5.1,
%   .claude/SESSION_HANDOFF.md): lets the existing "numeric sup-sampling ground truth" test
%   convention this codebase already uses everywhere apply directly to cPLQ's own output.
%
% [input]  fnd : functionNDomain array (nx1 or 1xn).
%          s   : 1x2 numeric point [s1, s2] (matching the dual variable order the piece's own
%                region/function were built with -- normally [s_1, s_2] as plq_1p.conjugateFunction
%                names them).
% [output] val : the piecewise value at s (NaN if s is not in any region's closure).
%          idx : index into fnd of the (first) region found containing s (0 if none).
%
% Membership: evaluate every inequality of a region's own ineqs at s (each ineqs(j) represents
%   ineqs(j) <= 0) with a small numeric tolerance, so points on a shared boundary between two
%   adjacent regions are accepted by whichever is checked first (consistent, since both regions
%   agree on the value there for a continuous piecewise function).
    tol = 1e-6;
    val = NaN; idx = 0;
    for i = 1:numel(fnd)
        r = fnd(i).d;
        if isempty(r) || isempty(r.ineqs), continue; end
        vars = r.vars;
        inside = true;
        for j = 1:numel(r.ineqs)
            iv = double(subs(r.ineqs(j).f, vars, s));
            if iv > tol
                inside = false;
                break;
            end
        end
        if inside
            val = double(subs(fnd(i).f.f, vars, s));
            idx = i;
            return
        end
    end
end
