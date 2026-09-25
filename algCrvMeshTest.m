classdef algCrvMeshTest < matlab.unittest.TestCase
% [slow] The EDGE LIST of the real curved-tie AlgCrv (c = -2, the diagonal-kink fixture that
% `biconjQ` returns as an AlgCrv), checked against the DEFINITION of a face boundary rather than
% against stored numbers: at every vertex point of a face's boundary exactly two of its arcs end.
%
% BUCKET: slow, although nothing here is symbolic. The fixture costs about 7 minutes to assemble and
% its vertices about 4 more (`algCrvVertices`, degree-16 resultants), which does not fit the normal
% bucket's 10 minutes on its own. The stages are separate tests so a break names its stage.
%
% WHY THIS FIXTURE: before arcs were traced, 31 of its 56 edges carried more than two vertices, its
% conic table spells y = 1 as y - 1, 1 - y, y^2 - 1 and (y - 1)^2, and it has two faces of zero area.
% Every one of those broke the degree condition on the way (DECISIONS 2026-09-22).

    methods (TestClassSetup)
        function build(t)
            fA = tri([0 0; 1 0; 1 1], 0);  fB = tri([0 0; 1 1; 0 1], -2);
            env = envelopePrune(envelopeAssemble(conjQ(fA), conjQ(fB)));
            t.G = assembleAlgCrv(envelopeMeshAlg(env));
            [t.V, t.info] = algCrvVertices(t.G);
            t.X = t.info.X;
        end
    end

    properties
        G
        V
        X
        info
    end

    methods (Test)
        function theVERTICESAreNamedWithoutRefusal(t)
            t.verifyEqual(size(t.X), [size(t.V,1), 2], 'a coordinate for every name');
            t.verifyEqual(size(unique(t.V, 'rows'), 1), size(t.V, 1), 'no name twice');
            t.verifyEqual(t.info.refused, 0, 'no edge pair refused');
        end

        function fStarStarIsFINITEOnAllOfDomF(t)
            % dom f** = conv dom f = the unit square, so every point of it, BOUNDARY INCLUDED, must
            % be owned by a face. On the domain edges x = 0 and y = 1 the bridge's two branches
            % meet, the fibre is a double root, and algBranch used to drop it as complex: f** was
            % +Inf at (0, 0.5), (0, 0.7) and (0.41, 1) (2026-09-25).
            s = linspace(0, 1, 401).';
            B = [s 0*s; 1+0*s s; s 1+0*s; 0*s s];
            [u, w] = meshgrid(linspace(0.013, 0.987, 23));
            X = [B; u(:) w(:)];
            v = t.G.eval(X);
            t.verifyEqual(sum(~isfinite(v)), 0, ...
                sprintf('points of dom f where f** is not finite: %s', mat2str(X(~isfinite(v),:), 4)));
        end

        function everyFaceWithAREAHasAClosedBoundary(t)
            [E, F, info] = algCrvIncidence(t.G, t.V, t.X);
            t.verifyEqual(info.untraced, 0, 'every edge with a vertex on it is traced');
            area = sampledArea(t.G);
            t.verifyGreaterThan(sum(area > 0), 1, 'the fixture has at least two faces with area');
            for k = find(area > 0)
                mine = find(F(:,1) == k | F(:,2) == k);
                t.verifyNotEmpty(mine, sprintf('face %d has area, so it has a boundary', k));
                t.verifyTrue(all(E(mine,5) == 1), ...
                             sprintf('face %d is bounded (sampled), so none of its arcs is a ray', k));
                P = zeros(0,2);
                for r = mine.'
                    for c = [3 4]
                        if E(r,c) > 0, P(end+1,:) = round(t.X(E(r,c),:), 6); end %#ok<AGROW>
                    end
                end
                [~, ~, j] = unique(P, 'rows');
                t.verifyEqual(unique(accumarray(j, 1)).', 2, ...
                              sprintf('face %d: every vertex ends exactly two of its arcs', k));
            end
        end
    end
end

% ==================================================================================================
function a = sampledArea(g)
% Grid points each face's evaluator assigns to it -- an independent read of where the faces are.
% The grid is shifted by an IRRATIONAL offset so no point lies on a rational line: a face of zero
% area (the fixture has two, on y = 1 and x = 0) must get no points, and an aligned grid gives it some.
    [xx, yy] = meshgrid(linspace(-3, 4, 71) + 0.01 * sqrt(2));
    [~, idx] = g.eval([xx(:) yy(:)]);
    a = accumarray(idx(idx > 0), 1, [g.nf 1]).';
end

function f = tri(V, c)
    q = [0 0 0 0, 2, 0, 2, 0, 0, 0];  q(8) = q(8) - c;  q(9) = q(9) + c;
    f = QuaPol(V, [1 2 1; 2 3 1; 3 1 1], q, [1 0; 1 0; 1 0]);
end
