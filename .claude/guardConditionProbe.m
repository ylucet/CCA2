function guardConditionProbe()
% ITEM 1, PART 2: test each N/R guard's CONDITION directly, on real conjQ output.
%
% WHY NOT JUST RUN THE PIPELINE. Running conjCPLQ to see whether a guard fires takes the better
% part of an hour per fixture, and a guard that does not fire proves nothing. But most of the N/R
% rows in `SUPPORT_MATRIX.md` are conditions on the SHAPE of a conjugate's subdivision -- "a curved
% ray edge", "a face with more than one curved edge", "the difference of the two candidate
% quadratics is an irreducible ellipse/hyperbola" -- and `conjQ` produces that subdivision exactly,
% in seconds. So the condition is tested where it lives.
%
% WHAT THIS ESTABLISHES AND WHAT IT DOES NOT. Finding the condition in a conjugate that `conjQ`
% actually returns from a legitimate `QuaPol` shows the SHAPE is reachable. It does not by itself
% prove the specific guard fires, because that also depends on the pipeline handing that shape to
% that routine. Stated per row below, and `.claude/reachabilityProbe.m('slow')` settles the ones
% that need the wired path.
%
% AN EDGE CONIC IS THE DIFFERENCE OF TWO CANDIDATE QUADRATICS, by construction: edge j of a QuaCon
% is the locus where the faces on its two sides agree. So `maxQuaPar:notDegenerate`'s condition and
% "an irreducible non-degenerate edge conic" are the same object, and irreducibility is the exact
% integer test det(M) ~= 0 on the conic's 3x3 matrix -- a singular M is a degenerate conic (a pair
% of lines, a single point), which is what "notDegenerate" is named after.

    fprintf('\n%-20s %6s %8s %8s %9s %8s %9s\n', 'fixture', 'ne', 'curved', 'irred', ...
            'curvedRay', 'multi', 'hyperbol');
    fprintf('%s\n', repmat('-', 1, 80));
    tot = zeros(1,5);
    for c = [-2 -0.5 0.5 2]
        for shape = ["square2" "poly2"]
            g = conjQ(kink(shape, c));
            r = scan(g);
            fprintf('%-20s %6d %8d %8d %9d %8d %9d\n', ...
                    sprintf('%s c=%+g', shape, c), g.ne, r(1), r(2), r(3), r(4), r(5));
            tot = tot + r;
        end
    end
    for nm = ["PD" "indefinite"]
        g = conjQ(single(nm));
        r = scan(g);
        fprintf('%-20s %6d %8d %8d %9d %8d %9d\n', ...
                sprintf('1-piece %s', nm), g.ne, r(1), r(2), r(3), r(4), r(5));
    end

    fprintf('\nTOTALS over the eight KINK fixtures: curved %d, irreducible %d, curved RAY %d,\n', ...
            tot(1), tot(2), tot(3));
    fprintf('faces with MORE THAN ONE curved edge %d, hyperbolic %d\n\n', tot(4), tot(5));
    fprintf(['Each nonzero total is a SUPPORT_MATRIX N/R row whose condition occurs in a\n' ...
             'conjugate conjQ really returns from a legitimate QuaPol. The single-piece rows\n' ...
             'are printed alongside because they are what the original scope claim was true of.\n']);
end

% ------------------------------------------------------------------------------------------------
function r = scan(g)
% r = [curved, irreducible curved, curved RAY, faces with >1 curved edge, hyperbolic]
    curved = false(g.ne,1);  irred = false(g.ne,1);  hyp = false(g.ne,1);
    for j = 1:g.ne
        k = g.edgeKind(j);
        curved(j) = ~strcmp(k, 'line');
        hyp(j)    = strcmp(k, 'hyperbolic');
        if curved(j)
            q = g.EcQ(j,:);
            M = [2*q(1), q(2), q(4); q(2), 2*q(3), q(5); q(4), q(5), 2*q(6)];
            irred(j) = ratQ.detExact(M) ~= 0;
        end
    end

    % A CURVED RAY: an edge whose conic is not a line and whose segment/ray flag says ray.
    % `SUPPORT_MATRIX` row "Curved ray edge -- QuaPar has no unbounded curved edge".
    curvedRay = 0;
    if ~isempty(g.E)
        curvedRay = sum(curved & (g.E(:,3) == 0));
    end

    % A FACE WITH MORE THAN ONE CURVED EDGE: row "one Ec slot per piece". FC{k} lists the edges
    % bounding face k, so this is a count over the faces, not over the edges.
    multi = 0;
    for k = 1:g.nf
        e = g.FC{k}(:,1);
        if sum(curved(e)) > 1, multi = multi + 1; end
    end

    r = [sum(curved), sum(irred), curvedRay, multi, sum(hyp)];
end

function f = single(which)
    if which == "PD", H = [2 0; 0 2]; else, H = [2 0; 0 -2]; end
    V = [0 0; 1 0; 1 1; 0 1];
    f = QuaPol(V, [1 2 1; 2 3 1; 3 4 1; 4 1 1], ...
               [0 0 0 0, H(1,1), H(1,2), H(2,2), 0, 0, 0], [1 0; 1 0; 1 0; 1 0]);
end

function f = kink(shape, c)
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
