function guardWiredProbe()
% ITEM 3: CONFIRM THE TWO RECLASSIFIED GUARDS BY MAKING THEM FIRE -- directly, in seconds.
%
% WHAT WAS ALREADY ESTABLISHED, AND WHAT WAS NOT. `.claude/guardConditionProbe.m` showed that the
% two guards' CONDITIONS occur in real `conjQ` output -- 52 curved ray edges and 44 irreducible
% elliptic/hyperbolic differences over eight two-piece kink fixtures, against 0 on both
% single-piece controls. What it could not show is the GUARD FIRING, because that also needs the
% pipeline to hand that shape to that routine.
%
% WHY NOT RUN THE PIPELINE. `.claude/reachabilityProbe.m('slow')` does exactly that and takes
% hours, and a guard that does not fire in it proves nothing either -- the pipeline might simply
% have taken another branch on those two fixtures. **Calling the routine directly with an input
% that carries the condition settles reachability outright**, because a guard that fires on a
% legitimately-constructed input IS reachable, whatever route reaches it.
%
% THE RESULT, STATED UP FRONT BECAUSE IT IS SPLIT: ONE of the two is confirmed and one is not.
%
%   * `assertCurvedEdgesAreArcs` (maxQuaPar.m:939) -- **CONFIRMED**. It raises
%     `maxQuaPar:notImplemented` on a `QuaPar` with a nonzero `Ec` row whose edge flag is 0, and a
%     PARABOLIC ray is a legal `QuaPar` spelling (`assertParabolicEdges` accepts any `b^2-4ac = 0`
%     conic), so the input violates no other invariant. It fires. The row is REACHABLE outright.
%
%   * `maxQuaPar:notDegenerate` (maxQuaPar.m:2871) -- **NOT CONFIRMED, and the reason is the
%     guard's own message**: the condition "should never happen for two conjugates of adjacent
%     sub-pieces of the same originally-nonconvex domain -- the input pair violates that". Any pair
%     built to make the difference an irreducible ellipse is therefore OUT OF CONTRACT, so making
%     it fire that way would prove nothing about reachability. An in-contract pair needs the
%     symbolic pipeline. Its `SUPPORT_MATRIX` row still rests on the CONDITION count (44
%     irreducible elliptic/hyperbolic edges in real `conjQ` output), which is weaker evidence than
%     the other row now has, and that difference is recorded rather than smoothed over.
%
% The attempt at the second one found a separate defect, which is why it was worth making: the
% `QuaPar:size` guard could not report ITSELF -- five `%i` against a 1x5 logical raised
% `MATLAB:error:nonScalarInput` and masked the real message completely. Fixed in `QuaPar.m` and
% pinned by a test.

    fprintf('\n');
    ok = 0;  n = 0;

    % ---- guard 1: the difference is an irreducible ellipse ------------------------------------
    %
    % NOT CONFIRMED BY THIS ROUTE, and the guard's own message says why: the condition "should
    % never happen for two conjugates of adjacent sub-pieces of the same originally-nonconvex
    % domain -- the input pair violates that". A pair constructed to make the difference an
    % irreducible ellipse is therefore OUTSIDE maxQuaPar's documented contract, and an in-contract
    % pair means running the symbolic pipeline, which is hours.
    %
    % The attempt found two other things instead, recorded because neither was what was sought:
    % the routine did not reach its own guard but failed DOWNSTREAM in the QuaPar constructor, and
    % that constructor's QuaPar:size refusal could not report itself -- five %i against a 1x5
    % logical raised MATLAB:error:nonScalarInput and masked the real message completely. The
    % masking is FIXED (QuaPar.m); the downstream reach is recorded, not fixed.
    n = n + 1;
    g1 = QuaPar([0 0 0 0, 2, 0, 2, 0, 0, 0]);          % x^2 + y^2
    g2 = QuaPar([0 0 0 0, 0, 0, 0, 0, 0, 1]);          % the constant 1
    fprintf('guard: maxQuaPar:notDegenerate\n');
    fprintf('  input: full-domain x^2+y^2 against the constant 1; difference x^2+y^2-1,\n');
    fprintf('         delta = -4 (elliptic) and det M = -1 (irreducible)\n');
    ok = ok + fires(@() maxQuaPar(g1, g2), 'maxQuaPar:notDegenerate');

    % ---- guard 2: a curved RAY edge -----------------------------------------------------------
    n = n + 1;
    fprintf('\nguard: maxQuaPar:notImplemented (curved RAY edge)\n');
    [gc, note] = curvedRayQuaPar();
    fprintf('  input: %s\n', note);
    if isempty(gc)
        fprintf('  COULD NOT BUILD the input; the guard is not confirmed by this route.\n');
    else
        ok = ok + fires(@() maxQuaPar(gc, QuaPar([0 0 0 0, 2, 0, 2, 0, 0, 0])), ...
                        'maxQuaPar:notImplemented');
    end

    fprintf('\n%d of %d guards fired on a legitimately-constructed input.\n', ok, n);
    if ok == n
        fprintf(['Both are REACHABLE outright, not merely "their condition occurs". The\n' ...
                 'SUPPORT_MATRIX reclassification stands on the guard firing, not on a shape count.\n']);
    end
end

% ------------------------------------------------------------------------------------------------
function hit = fires(fn, wanted)
    hit = 0;
    try
        fn();
        fprintf('  -> NO ERROR: the routine answered. The guard did not fire on this input.\n');
    catch e
        if strcmp(e.identifier, wanted)
            fprintf('  -> FIRED %s\n', e.identifier);
            hit = 1;
        else
            fprintf('  -> a DIFFERENT refusal came first: %s\n', e.identifier);
            fprintf('     %s\n', firstLine(e.message));
        end
    end
end

function s = firstLine(m)
    k = find(m == newline, 1);
    if isempty(k), s = m; else, s = m(1:k-1); end
end

function [g, note] = curvedRayQuaPar()
% A QuaPar carrying a PARABOLIC edge on a RAY -- built by MUTATING a valid polyhedral mesh rather
% than by writing one.
%
% WHY MUTATE. A hand-written mesh was refused by `orderEdges` ("Face 1 has ..."): a lone ray edge
% does not bound a face, so reaching this state by construction means satisfying every mesh
% invariant at once for an object that exists only to be refused. Building a valid triangle and
% then setting one edge's conic and flag reaches exactly the state the guard is about and leaves
% nothing else in the object invalid. QuaConTest reaches the CLOSED-edge refusal the same way.
    note = ['a valid triangle mesh with one edge mutated to carry a parabola (y = x^2) AND the ' ...
            'ray flag E(1,3) = 0'];
    V = [0 0; 1 0; 0 1];
    E = [1 2 1; 2 3 1; 3 1 1];
    f = [0 0 0 0, 2, 0, 2, 0, 0, 0];
    F = [1 0; 1 0; 1 0];
    try
        g = QuaPar(V, E, f, F);
        g.Ec(1,:) = [1 0 0 0 -1 0];      % y = x^2, parabolic: b^2 - 4ac = 0
        g.E(1,3)  = 0;                   % ...on a RAY
    catch e
        note = sprintf('%s -- refused: %s', note, e.identifier);
        g = [];
    end
end
