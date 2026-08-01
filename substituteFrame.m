function cj = substituteFrame(cj0, frame)
% substituteFrame  Read a conjugate computed in the x*y frame back in the original dual variables.
%
% [input]  cj0   : functionNDomain array, h* in the z-frame, in variables s_1, s_2.
%          frame : struct('M',M,'a',a,'c0',c0) from xyFrame.
% [output] cj    : functionNDomain array for f*(s) = h*(M's - a) - c0, same variables.
%
% Both the FUNCTION and every REGION inequality are substituted -- the dual cells move with the
% frame just as the values do, and substituting only the functions would leave each piece
% carrying the right value on the wrong cell. The substitution is simultaneous, via temporaries,
% because s_1 and s_2 each appear on both sides: writing s_1 first and then s_2 would feed the
% already-rewritten s_1 back into the second replacement.
    M = frame.M; a = frame.a(:); c0 = frame.c0;
    s1 = sym('s_1'); s2 = sym('s_2');
    t1 = sym('cca2_frame_t1'); t2 = sym('cca2_frame_t2');

    Mt = M.';                       % transpose, NOT ctranspose: these are real matrices but the
                                    % expressions they build are symbolic (see conjConvexOverPiece)
    e1 = Mt(1,1)*t1 + Mt(1,2)*t2 - a(1);
    e2 = Mt(2,1)*t1 + Mt(2,2)*t2 - a(2);

    cj = functionNDomain.empty();
    for i = 1:numel(cj0)
        fi = subs(cj0(i).f.f, [s1 s2], [e1 e2]);
        fi = subs(fi, [t1 t2], [s1 s2]) - c0;

        r0 = cj0(i).d;
        g = sym.empty(1,0);
        for j = 1:size(r0.ineqs,2)
            gj = subs(r0.ineqs(j).f, [s1 s2], [e1 e2]);
            g(j) = expand(subs(gj, [t1 t2], [s1 s2]));
        end
        rk = region(g, [s1 s2]);
        if isempty(rk)
            continue                % the cell degenerated under the map; nothing to carry over
        end
        cj = [cj, functionNDomain(symbolicFunction(simplify(expand(fi))), rk)]; %#ok<AGROW>
    end
end
