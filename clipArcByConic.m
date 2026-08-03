function [status, Xnew, uCross] = clipArcByConic(Ec, X0, X1, Ecut)
% clipArcByConic  Clip a parabola ARC against the region bounded by ANOTHER conic -- the primitive
%   arc-vs-arc face clipping needs, and the exact analogue of clipArcByHalfPlane for a curved cut.
%
% [input]  Ec        : 1x6 conic [a b c d e f] of the arc's own parabola (b^2-4ac == 0, checked by
%                      parabolaArcFrame). X0, X1 : 1x2 endpoints, both already on Ec.
%          Ecut      : 1x6 conic of the CUTTING curve. The region KEPT is
%                      { x : evalConic(Ecut,x) >= 0 }, which is facePoly's own sign convention --
%                      it normalizes a face's curveEc so the conic is > 0 on that face's interior,
%                      so passing polyL.curveEc here keeps polyL's inside.
% [output] status    : 'inside'  -> the whole arc is kept; Xnew = [X0;X1].
%                      'outside' -> the whole arc is removed; Xnew = [].
%                      'cut'     -> exactly one crossing strictly inside the arc; Xnew = [A;B]
%                                   with the outside endpoint replaced by it.
%                      'twice'   -> TWO crossings strictly inside, i.e. both endpoints on the same
%                                   side and the arc leaving and re-entering. Xnew holds the two
%                                   crossing points in arc order. Reported, never guessed: the
%                                   surviving set is then either two components or one bounded by
%                                   two sub-arcs, and the caller must decide.
%          uCross    : the crossing parameters in the arc's own u-frame, for callers that need to
%                      re-derive positions without a second root solve.
%
% METHOD. parabolaArcFrame(Ec) supplies the intrinsic (u,v) frame in which every point of the
% parabola is point(u), with u a GLOBAL monotone parameter -- so ordering by u orders along the
% arc, and "is this root on THIS arc" is a range test rather than a reprojection. Its conicCoeffs
% restricts a second conic to that parametrization, giving an explicit QUARTIC in u (point(u) is
% quadratic, so a conic in it has degree 4). Crossings are the real roots of that quartic inside
% the arc's u-range. This is exactly the pattern lineCoeffs already uses for a straight cut; the
% only change is degree 2 -> 4.
%
% WHY A QUARTIC IS NOT OVERKILL HERE. Two parabolas can meet in up to 4 points in general. In this
% pipeline both curves are PARABOLAS (Stage C: no ellipse or hyperbola can arise as a dual
% boundary), and what matters is only how many roots land INSIDE the arc's own u-range; that is
% what the status distinguishes. Roots outside the range belong to other arcs of the same conic
% and are correctly ignored.
%
% SIGN CONVENTION AND TANGENCY. A root of even multiplicity is a TANGENCY, not a crossing: the arc
% touches the cutting curve without changing side. Those are excluded by requiring an actual sign
% change across the root, tested on the quartic itself rather than by inspecting multiplicities,
% which is numerically the more robust of the two.

    fr = parabolaArcFrame(Ec, 'clipArcByConic');
    u0 = fr.uOf(X0);
    u1 = fr.uOf(X1);
    lo = min(u0, u1); hi = max(u0, u1);
    if hi - lo <= 0
        error('clipArcByConic:degenerateArc', ...
            'the arc''s two endpoints have the same frame parameter u; it has no extent.');
    end

    q = fr.conicCoeffs(Ecut);                 % [A4 A3 A2 A1 A0], value of Ecut along the arc
    scale = max(1, max(abs(q)));
    tol   = 1e-9 * scale;
    span  = hi - lo;

    r = roots(q);
    r = r(abs(imag(r)) <= 1e-8*max(1,abs(r)));
    r = sort(real(r));

    % Keep only roots strictly inside the arc, and only those the value actually CHANGES SIGN
    % across -- a tangency touches without crossing and must not split the arc.
    inside = [];
    for k = 1:numel(r)
        uk = r(k);
        if uk <= lo + 1e-9*span || uk >= hi - 1e-9*span, continue, end
        h = 1e-6*span;
        vBefore = polyval(q, max(lo, uk - h));
        vAfter  = polyval(q, min(hi, uk + h));
        if sign2(vBefore, tol) * sign2(vAfter, tol) < 0
            inside(end+1) = uk; %#ok<AGROW>
        end
    end
    inside = sort(inside);

    if isempty(inside)
        % No crossing: the whole arc is on one side. Decide it at an INTERIOR sample, not at an
        % endpoint -- an endpoint can sit exactly on the cutting curve (the two faces share it),
        % where the sign is 0 and says nothing.
        vMid = polyval(q, 0.5*(lo+hi));
        if vMid >= -tol
            status = 'inside'; Xnew = [X0; X1]; uCross = [];
        else
            status = 'outside'; Xnew = []; uCross = [];
        end
        return
    end

    if numel(inside) == 1
        uc = inside(1);
        Xc = fr.point(uc);
        % Replace whichever endpoint is on the removed side.
        if polyval(q, lo + 1e-6*span) >= -tol
            keepU = lo;                        % the lo end survives
        else
            keepU = hi;
        end
        if keepU == lo
            A = fr.point(lo); B = Xc;
        else
            A = Xc;          B = fr.point(hi);
        end
        % Return in the SAME orientation the caller gave, so X0 stays X0's end of the arc.
        if u0 <= u1, Xnew = [A; B]; else, Xnew = [B; A]; end
        status = 'cut'; uCross = uc;
        return
    end

    if numel(inside) == 2
        status = 'twice';
        if u0 <= u1, uCross = inside; else, uCross = flip(inside); end
        Xnew = [fr.point(uCross(1)); fr.point(uCross(2))];
        return
    end

    error('clipArcByConic:tooManyCrossings', ...
        ['the cutting conic crosses this arc %d times; a QuaPar edge cannot be split into that ' ...
         'many sub-arcs by one clip.'], numel(inside));
end

% ------------------------------------------------------------------------------------------
function s = sign2(v, tol)
    if v >  tol, s =  1; elseif v < -tol, s = -1; else, s = 0; end
end
