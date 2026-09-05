classdef domainIsConvexDTest < matlab.unittest.TestCase
% Unit tests for domainIsConvexD, the doubles domain-convexity test createDom now uses.
%
% BUCKET: fast (a few 2x2 determinants; no symbolic call anywhere).
%
% It is a predicate on the MESH ARRAYS, so it is tested on mesh arrays directly. That is what makes
% the recession-direction branch reachable at all: for a cone whose ray direction points sit on the
% rays, a violating direction always comes with a violating vertex, and the vertex loop returns
% first. Placing the direction point closer in separates the two, which is exactly the case the
% branch exists for.

    methods (Test)

        function theHalfStarIsRejectedAndTheWitnessIsInTheFixture(testCase)
        % QuaPol.examples{16}: five rays from the origin sweeping more than a half-turn. The
        % boundary walk this replaced called it convex. (-1,-1) and (1,0) are in it, (0,-0.5) is
        % not -- see PLQVCTest/testDomIsConvexIsFalseForAnUNBOUNDEDNonConvexDomain.
            V = [0,0; -1,0; -1,-1; 0,1; 1,1; 1,0];
            E = [1,2,0; 1,3,0; 1,4,0; 1,5,0; 1,6,0];
            F = [0,1; 1,2; 2,3; 3,4; 4,0];
            testCase.verifyFalse(domainIsConvexD(V, E, F));
        end

        function aRECESSIONDirectionAloneCanRefuteIt(testCase)
        % The branch the vertex loop cannot reach. The domain is bounded below by the x-axis, so
        % every vertex satisfies y >= 0 -- including (1, 0.5), the direction point of the ray from
        % (0,1). But that ray's DIRECTION is (1,-0.5), which leaves the half-plane, so the domain
        % runs out of its own boundary constraint and is not convex.
            V = [-1 0; 1 0; 0 1; 1 0.5];
            E = [1 2 1;        % the x-axis segment, domain above it
                 3 4 0];       % a ray from (0,1) through (1,0.5): direction (1,-0.5)
            % Edge 1 runs (-1,0) -> (1,0), so its direction is +x and its LEFT is +y. The finite
            % face therefore goes in column 1, which puts the domain in y >= 0.
            F = [1 0;          % boundary edge
                 1 2];         % interior edge
            testCase.verifyFalse(domainIsConvexD(V, E, F), ...
                'the ray leaves the half-plane its own boundary edge defines');
        end

        function aCONVEXConeIsAccepted(testCase)
        % The companion, so the previous test is about the direction and not about the shape being
        % unbounded: the same boundary edge with a ray that stays inside.
            V = [-1 0; 1 0; 0 1; 1 2];
            E = [1 2 1; 3 4 0];      % ray from (0,1) through (1,2): direction (1,1), stays in y>=0
            F = [1 0; 1 2];
            testCase.verifyTrue(domainIsConvexD(V, E, F));
        end

        function aBOUNDEDNonConvexPolygonIsRejectedAndItsConvexHullIsNot(testCase)
        % The case the boundary walk already got right, kept so a regression there is visible.
            square = [0 0; 1 0; 1 1; 0 1];
            E = [1 2 1; 2 3 1; 3 4 1; 4 1 1];
            F = [1 0; 1 0; 1 0; 1 0];
            testCase.verifyTrue(domainIsConvexD(square, E, F));

            % the same square with one corner pushed INSIDE: a reflex vertex
            dart = [0 0; 1 0; 0.4 0.4; 0 1];
            testCase.verifyFalse(domainIsConvexD(dart, E, F));
        end

        function aTHINMeshIsREFUSEDRatherThanCalledConvex(testCase)
        % Every edge with +infinity on both sides means the domain has dimension < 2, and the
        % containment test says nothing about it -- the loop would find no boundary edge and return
        % a vacuous true. Measured before wiring this in (.claude/domConvexDiff.m): four fixtures
        % of the corpus, examples2{5,6,7,10}, are exactly that shape and orderEdges rules three of
        % them non-convex. createDom decides them before reaching here; this raises so that a
        % future caller cannot get the vacuous answer by accident.
            V = [0 0; 1 0; 1 1];
            E = [1 2 1; 2 3 1];
            F = [0 0; 0 0];
            testCase.verifyError(@() domainIsConvexD(V, E, F), 'PLQ:domainIsConvexD:thinMesh');
        end

        function itAgreesWithTheEXACTTestOnEveryTwoDimensionalFixture(testCase)
        % The differential check that matters: isConvexQ runs the same algorithm over the integers
        % with no tolerance, so where both can run they must agree about the DOMAIN. isConvexQ
        % reports on the whole function, so only its domain verdict is compared -- read off the
        % reason string it returns.
            L = [QuaPol.examples(), QuaPol.examples2(), QuaPol.examplesNonconvex(), ...
                 QuaPol.examplesDiscontinuous()];
            nCompared = 0;
            for k = 1:numel(L)
                f = L{k};
                if f.dom.dim < 2 || f.nv == 0, continue, end
                try, [tf, why] = isConvexQ(f); catch, continue, end
                exactSaysDomainBad = ~tf && contains(why, 'domain');
                d = domainIsConvexD(f.V, f.E, f.F);
                nCompared = nCompared + 1;
                if exactSaysDomainBad
                    testCase.verifyFalse(d, sprintf('fixture %d: exact rules the domain out', k));
                else
                    testCase.verifyTrue(d, sprintf('fixture %d: exact accepts the domain', k));
                end
            end
            testCase.verifyGreaterThan(nCompared, 10, 'too few fixtures compared');
        end

    end
end
