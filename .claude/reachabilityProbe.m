function reachabilityProbe(varargin)
% ITEM 1: ARE THE `N/R` GUARDS ACTUALLY UNREACHABLE? Re-tested against MULTI-PIECE input.
%
% WHY THIS HAD TO BE RE-RUN. `SUPPORT_MATRIX.md` section 0 classified guards as **N/R** ("not
% reachable from a QuaPol conjugate/biconjugate") on the strength of a scope claim that a QuaPol
% conjugate never produces hyperbolic edges. That claim was REFUTED on 2026-09-06: it holds for a
% SINGLE piece and fails for a two-piece KINK, where `conjQ` returns 2 to 4 hyperbolic edges plus
% elliptic ones. Every N/R classification resting on it is therefore unverified, and the one that
% names an ellipse or hyperbola outright -- `maxQuaPar:notDegenerate`, "the difference of the two
% candidate quadratics is an irreducible ellipse/hyperbola" -- is the obvious suspect.
%
% THE METHOD, and it is deliberately blunt. Run the WIRED pipelines on fixtures that a user could
% legitimately write, and record every named refusal that comes back. A guard that fires here is
% reachable, full stop -- no argument about scope can survive an input that reaches it. A guard
% that does NOT fire is not thereby proved unreachable, and this probe never claims that; it
% partitions the N/R list into "refuted, with the input that does it" and "still unrefuted".
%
% THE FIXTURES ARE ORDINARY. A square and a quadrilateral, each split by a diagonal into two
% positive definite quadratic pieces that agree on the shared edge, with the kink swept convex /
% flat / non-convex by one scalar. That construction is DECISIONS 2026-09-05's, and its whole point
% is that it produces a legal continuous multi-piece QuaPol for every value of the knob.

    fixtures = {};
    for c = [-2 -0.5 0.5 2]
        fixtures{end+1} = {sprintf('square2 kink c=%+g', c), kink('square2', c)}; %#ok<AGROW>
        fixtures{end+1} = {sprintf('poly2   kink c=%+g', c), kink('poly2',   c)}; %#ok<AGROW>
    end
    fixtures{end+1} = {'square1 single PD',         single([2 0; 0  2])};
    fixtures{end+1} = {'square1 single indefinite', single([2 0; 0 -2])};

    % THE SYMBOLIC PAIR IS OPTED INTO, NOT DEFAULTED. conjCPLQ runs for the better part of an hour
    % on one fixture (see .claude/suite.sh's bucket notes), so ten fixtures times two operators is
    % a night on its own. The exact pair answers in seconds and reaches the conjQ/biconjQ guards;
    % the symbolic pair is what reaches the maxQuaPar ones, and is run separately on two fixtures.
    %
    %   reachabilityProbe()          exact operators, every fixture
    %   reachabilityProbe('slow')    add conjCPLQ/biconjCPLQ, first two fixtures only
    ops = {'conjQ', @conjQ; 'biconjQ', @biconjQ};
    slow = nargin >= 1 && strcmp(varargin{1}, 'slow');
    if slow
        ops = [ops; {'conjCPLQ', @conjCPLQ; 'biconjCPLQ', @biconjCPLQ}];
        fixtures = fixtures(1:2);
    end

    seen = containers.Map('KeyType','char','ValueType','any');
    fprintf('\n%-26s %-12s %s\n', 'fixture', 'operator', 'outcome');
    fprintf('%s\n', repmat('-', 1, 92));
    for i = 1:numel(fixtures)
        nm = fixtures{i}{1};  f = fixtures{i}{2};
        for k = 1:size(ops,1)
            try
                g = ops{k,2}(f);
                fprintf('%-26s %-12s ok (%s)\n', nm, ops{k,1}, describe(g));
            catch e
                fprintf('%-26s %-12s %s\n', nm, ops{k,1}, e.identifier);
                if ~isKey(seen, e.identifier), seen(e.identifier) = {}; end
                v = seen(e.identifier);  v{end+1} = sprintf('%s / %s', nm, ops{k,1});
                seen(e.identifier) = v;
            end
        end
    end

    fprintf('\n\n=== EVERY IDENTIFIER THAT FIRED, with a witness ===\n');
    ks = keys(seen);
    for i = 1:numel(ks)
        v = seen(ks{i});
        fprintf('  %-46s  %d witness(es), e.g. %s\n', ks{i}, numel(v), v{1});
    end
    fprintf(['\nAn identifier listed here is REACHABLE from a legitimate QuaPol, whatever\n' ...
             'SUPPORT_MATRIX.md says about it. Absence from this list is NOT a proof of\n' ...
             'unreachability -- it means this probe did not reach it.\n']);
end

% ------------------------------------------------------------------------------------------------
function s = describe(g)
    try
        s = sprintf('%s nf=%d', class(g), g.nf);
    catch
        s = class(g);
    end
end

function f = single(H)
    V = [0 0; 1 0; 1 1; 0 1];
    f = QuaPol(V, [1 2 1; 2 3 1; 3 4 1; 4 1 1], ...
               [0 0 0 0, H(1,1), H(1,2), H(2,2), 0, 0, 0], [1 0; 1 0; 1 0; 1 0]);
end

function f = kink(shape, c)
% Two positive definite pieces that agree on the shared edge FOR EVERY c, because the perturbation
% is a multiple of that edge's own linear form, which vanishes there. One knob sweeps convex, flat
% and non-convex without ever producing a mesh that is not a function. DECISIONS 2026-09-05.
    q0 = [0 0 0 0, 2, 0, 2, 0, 0, 0];
    if shape == "square2"
        V = [0 0; 1 0; 1 1; 0 1];  lin = [-1 1 0];
    else
        V = [-2 0; 2 0; 1 2; -1 2];  lin = [0 1 0];
    end
    E = [1 2 1; 2 3 1; 3 4 1; 4 1 1; 1 3 1];
    F = [1 0; 1 0; 2 0; 2 0; 2 1];
    q1 = q0;  q1(8) = q1(8) + c*lin(1);  q1(9) = q1(9) + c*lin(2);  q1(10) = q1(10) + c*lin(3);
    f = QuaPol(V, E, [q0; q1], F);
end
