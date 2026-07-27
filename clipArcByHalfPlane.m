function [status, Xnew] = clipArcByHalfPlane(Ec, X0, X1, nrm, c)
% clipArcByHalfPlane  Clip a parabola ARC against a half-plane -- the missing geometric primitive
%   maxQuaPar.m's own TODO calls out (a curved/parabolic input edge, e.g. conjBilinearXYoneCE's or
%   conjIndefiniteQuadTriangle-with-1-convex-edge's output), needed before maxQuaPar can combine
%   QuaPar conjugates that have curved edges. Part of the Phase 2 scoping (DESIGN.md II.5.1 /
%   .claude/SESSION_HANDOFF.md "Phase 2" -- go one small step at a time, validate each in
%   isolation before composing).
%
% STATUS: implemented and validated STANDALONE (clipArcByHalfPlaneTest.m) against hand-derived and
%   independently-constructed (rotated/shifted) parabola cases, AND now WIRED into maxQuaPar.m --
%   see its clipPolyHalfPlaneCurved, which calls this per curved cell edge while threading the
%   arc's position through each of clipPolyHalfPlane's keep/cut branches, and maxQuaPar.m's own
%   header for the end-to-end validation of that path. One caveat carried there: within
%   maxQuaPar's supported regime the 'cut' status below is never actually reached (the other
%   operand's face boundaries are TANGENT to the parabola there, so a clip keeps the arc whole or
%   reduces it to a point), so that branch stays covered by this file's own direct tests rather
%   than end to end.
%
%   The (u,v)-frame construction itself now lives in parabolaArcFrame.m, shared with maxQuaPar.m's
%   arc midpoint/tangent and arc-bulge tests so the two cannot drift apart.
%
% [input]  Ec  : 1x6 conic [a b c d e f] for a*x^2+b*xy+c*y^2+d*x+e*y+f=0, a genuine PARABOLA
%                (b^2-4ac==0, checked). X0, X1 : 1x2 points, both already known to lie on Ec (a
%                QuaPar curved edge's own two given endpoints).
%          nrm, c : half-plane {nrm*x'<=c} to clip by, nrm a 1x2 row vector.
% [output] status : 'inside'  -> the whole arc is retained; Xnew = [X0;X1] unchanged.
%                    'outside' -> the whole arc is removed; Xnew = [].
%                    'cut'     -> exactly one endpoint was outside and is replaced by the single
%                                 new crossing point; Xnew = [X0new; X1new] (only one row changed).
%
% METHOD: rewrite the parabola in its own (u,v) frame, v = q(u) = -(lam*u^2+du*u+f)/dv -- see
%   parabolaArcFrame.m, which owns that construction (it is shared with maxQuaPar.m's curved-edge
%   clipping) and documents why u is a GLOBAL, monotonic parameter for every point of the conic, so
%   ordering two points by u correctly orders them along the finite ("short") arc between them
%   whenever that arc doesn't need to detour through infinity -- always true here, since a QuaPar
%   curved edge is by construction never an unbounded curved ray (see QuaPar.m's own header).
%   nrm*X(u)'-c is then an explicit quadratic in u (fr.lineCoeffs): solved directly for u, then
%   converted back to a point via fr.point -- no line-conic-intersection-then-reproject detour
%   needed, and no ambiguity about which of the (up to 2) roots of the FULL conic-vs-line
%   intersection is the one actually on this arc (picked by u-range membership instead).
    tol = 1e-8;
    fr = parabolaArcFrame(Ec, 'clipArcByHalfPlane');
    Xu = fr.point;

    u0 = fr.uOf(X0); u1 = fr.uOf(X1);
    v0 = nrm*X0' - c; v1 = nrm*X1' - c;
    in0 = v0 <= tol; in1 = v1 <= tol;

    if in0 && in1
        status = 'inside'; Xnew = [X0; X1]; return
    end
    if ~in0 && ~in1
        status = 'outside'; Xnew = []; return
    end

    % Exactly one endpoint outside: solve the quadratic-in-u restriction of nrm*X(u)'-c directly
    % (substitute q(u)=-(lam u^2+du u+f)/dv into nrm*X(u)'-c), then keep whichever root lies
    % within this arc's own [u0,u1] span (up to a small tolerance).
    A = fr.lineCoeffs(nrm, c);
    roots_u = solveQuadLocal(A(1), A(2), A(3));
    uLo = min(u0,u1); uHi = max(u0,u1);
    marg = 1e-6*(1+abs(uLo)+abs(uHi));
    cand = roots_u(roots_u >= uLo-marg & roots_u <= uHi+marg);
    if isempty(cand)
        error('clipArcByHalfPlane:internal', ...
            'No crossing of the clip line found within the arc''s own u-span.');
    end
    ustar = cand(1);
    if numel(cand) > 1
        % Both roots landed in range (only possible right at a near-tangential boundary case);
        % pick whichever is closer to the OUTSIDE endpoint's own u (the genuine crossing nearest
        % where the sign actually flips).
        if in0, uOut = u1; else, uOut = u0; end
        [~,ix] = min(abs(cand - uOut));
        ustar = cand(ix);
    end
    Xstar = Xu(ustar);
    if in0
        Xnew = [X0; Xstar];
    else
        Xnew = [Xstar; X1];
    end
    status = 'cut';
end

function r = solveQuadLocal(A,B,C)
    tol = 1e-9*(1+abs(A)+abs(B)+abs(C));
    if abs(A) <= tol
        if abs(B) <= tol, r = []; return; end
        r = -C/B; return
    end
    disc = B^2 - 4*A*C;
    if disc < -tol, r = []; return; end
    disc = max(disc,0);
    r = [(-B+sqrt(disc))/(2*A), (-B-sqrt(disc))/(2*A)];
end
