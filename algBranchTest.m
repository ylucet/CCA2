classdef algBranchTest < matlab.unittest.TestCase
% Unit tests for algBranch, the evaluator of ONE BRANCH of an implicit algebraic face.
%
% BUCKET: fast (a quartic root-find per point; no symbolic call anywhere).
%
% THE ORACLE IS A CLOSED FORM IN EVERY CASE. Each fixture is a polynomial whose branches are known
% by hand, so agreement is evidence about this routine and not about a shared expression. The last
% fixture is the real one: the bridge face DECISIONS.md 2026-09-05 derived and checked against the
% primal definition of the envelope, written implicitly as Alg.m records it.

    methods (Test)

        function theTwoBranchesOfACONEAreTheTwoSIGNSOfANorm(testCase)
        % p = z^2 - (x^2 + y^2), whose branches are -|x| and +|x| -- the shape every bridge face
        % has (affine plus a norm), with the affine part zero. Branch 1 is the LOWER root and
        % branch 2 the upper, which is the ascending convention algBranch documents.
            P = zeros(3,3,3);
            P(1,1,3) =  1;      % z^2
            P(3,1,1) = -1;      % -x^2
            P(1,3,1) = -1;      % -y^2
            X = [3 4; 1 0; 0 1; -5 12; 0.5 0.5; -2 -2];
            r = sqrt(sum(X.^2, 2));
            testCase.verifyEqual(algBranch(P, X, 1), -r, 'AbsTol', 1e-12, 'lower branch');
            testCase.verifyEqual(algBranch(P, X, 2),  r, 'AbsTol', 1e-12, 'upper branch');
        end

        function aBRANCHTHATDOESNOTEXISTThereIsNaNAndNotAWrongNumber(testCase)
        % p = z^2 - (x - 2) has NO real root where x < 2, and exactly one (a double) at x = 2. The
        % answer there is "this face does not reach this point", and the one thing it must never be
        % is a plausible number. A third branch never exists at all.
            P = zeros(2,1,3);
            P(1,1,3) =  1;      % z^2
            P(2,1,1) = -1;      % -x
            P(1,1,1) =  2;      % +2
            X = [1 0; 0 0; -5 0; 3 0; 6 0];
            z1 = algBranch(P, X, 1);
            testCase.verifyTrue(all(isnan(z1(1:3))), 'no real root where x < 2');
            testCase.verifyEqual(z1(4:5), [-1; -2], 'AbsTol', 1e-12, 'lower branch where it exists');
            testCase.verifyEqual(algBranch(P, X, 2), [NaN;NaN;NaN;1;2], 'AbsTol', 1e-12);
            testCase.verifyTrue(all(isnan(algBranch(P, X, 3))), 'there is no third branch');
        end

        function theBRANCHCOUNTIsReportedWITHMULTIPLICITY(testCase)
        % The second output exists for the producer, not the reader: a face is only well-formed
        % where the number of real roots is constant, so the assembly has to be able to SEE the
        % count rather than infer it from a value that came back NaN.
        %
        % WITH MULTIPLICITY, and that is a decision rather than an accident. At x = 0 the two
        % branches of z^2 = x MEET, and `roots` returns the double root as two entries. Collapsing
        % them to one would need a tolerance -- `roots` separates a true double root by about
        % sqrt(eps), never exactly -- and that tolerance would be deciding a STRUCTURAL question
        % (do two branches merge here), which is the one kind of question this project never lets a
        % tolerance answer. So the count is the multiset's, branches 1 and 2 agree at the merge, and
        % the merge locus is the discriminant curve, which a well-formed face has on its BOUNDARY.
            P = zeros(2,1,3);
            P(1,1,3) = 1;  P(2,1,1) = -1;          % z^2 - x
            [~, info] = algBranch(P, [-1 0; 0 0; 4 0; 9 0], 1);
            testCase.verifyEqual(info.nreal, [0; 2; 2; 2], 'multiplicity counted');
            testCase.verifyEqual(info.allReal{4}, [-3; 3], 'AbsTol', 1e-12, 'sorted ascending');
            % the two branches AGREE at the merge, which is what makes the index continuous there
            X0 = [0 0];
            testCase.verifyEqual(algBranch(P, X0, 1), algBranch(P, X0, 2), 'AbsTol', 1e-7);
        end

        function aQUARTICBranchIsSelectedByIndexAndNotByMagnitude(testCase)
        % Four distinct real roots at every point, so every index is exercised and an implementation
        % that sorted by |z| rather than by z would get two of them wrong.
        % p(z) = (z+3)(z+1)(z-1)(z-2) shifted by x: roots are {-3,-1,1,2} + x.
            base = poly([-3 -1 1 2]);              % highest power first
            base = fliplr(base);                   % algBranch's layout: ascending in z
            X = [0 0; 0.25 0; -0.5 0];
            for t = 1:size(X,1)
                s = X(t,1);
                % substitute z -> z - s by expanding the shifted product directly
                c = fliplr(poly([-3 -1 1 2] + s));
                P = zeros(1,1,5);  P(1,1,:) = c;
                want = sort([-3 -1 1 2] + s).';
                for bIdx = 1:4
                    testCase.verifyEqual(algBranch(P, X(t,:), bIdx), want(bIdx), 'AbsTol', 1e-10, ...
                        sprintf('shift %g, branch %d', s, bIdx));
                end
            end
            testCase.verifyEqual(numel(base), 5, 'the quartic has five coefficients');
        end

        function theMEASUREDBRIDGEFACEIsReproducedFromItsImplicitForm(testCase)
        % THE FIXTURE THAT MATTERS. DECISIONS.md 2026-09-05 derived, for q1 = |v|^2/2 and
        % q2 = |v-a|^2 + c, the envelope over the tie circle
        %       co f(x) = <2a, x-2a> + R*|x-2a| + 2c ,      R = 2*sqrt(|a|^2 - c) ,
        % and checked it against the primal definition to 1.5e-11. Alg.m stores the same face
        % implicitly as the quadric cone
        %       (z - <2a,x-2a> - 2c)^2 - R^2*|x-2a|^2 = 0 ,   R^2 = 4(|a|^2 - c)  RATIONAL.
        % This asserts that the implicit form with the UPPER branch reproduces the closed form --
        % i.e. that the representation this type is being built around actually holds the face.
            a = [3; 1];  c = -2;                    % |a|^2 - c = 12 > 0, so R^2 = 48, R irrational
            R2 = 4*(a.'*a - c);
            L  = 2*a;                               % the affine part's gradient
            k  = -(L.'*(2*a)) + 2*c;                % so the affine part is <L,x> + k

            % (z - <L,x> - k)^2 - R2*|x - 2a|^2 = 0, expanded into the (x,y,z) coefficient array.
            % Rather than expand by hand -- which is how a fixture silently encodes the same slip
            % as the code -- build it from two elementary pieces and combine them numerically.
            P = zeros(3,3,3);
            % (z - L1 x - L2 y - k)^2
            v = [-k, -L(1), -L(2), 1];              % [1, x, y, z] coefficients of the inner form
            idx = {[1 1 1], [2 1 1], [1 2 1], [1 1 2]};
            for i = 1:4
                for j = 1:4
                    e = idx{i} + idx{j} - 1;
                    P(e(1), e(2), e(3)) = P(e(1), e(2), e(3)) + v(i)*v(j);
                end
            end
            % - R2 * ((x - 2a1)^2 + (y - 2a2)^2)
            P(3,1,1) = P(3,1,1) - R2;
            P(2,1,1) = P(2,1,1) + R2*2*(2*a(1));
            P(1,1,1) = P(1,1,1) - R2*(2*a(1))^2;
            P(1,3,1) = P(1,3,1) - R2;
            P(1,2,1) = P(1,2,1) + R2*2*(2*a(2));
            P(1,1,1) = P(1,1,1) - R2*(2*a(2))^2;

            X = [7 4; 6 1; 6.5 3; 10 -2; 2 5; 6 1.001];
            d = X - (2*a).';
            want = (X*L) + k + sqrt(R2)*sqrt(sum(d.^2, 2));     % the closed form, upper branch
            testCase.verifyEqual(algBranch(P, X, 2), want, 'RelTol', 1e-10, 'AbsTol', 1e-10);

            % and the LOWER branch is the other sign of the norm -- the face that is NOT the
            % envelope. Storing the polynomial without the branch index would not distinguish them,
            % which is the whole reason the index is part of the representation.
            wantLo = (X*L) + k - sqrt(R2)*sqrt(sum(d.^2, 2));
            testCase.verifyEqual(algBranch(P, X, 1), wantLo, 'RelTol', 1e-10, 'AbsTol', 1e-10);
            testCase.verifyGreaterThan(max(want - wantLo), 1, 'the two branches are far apart');
        end

        function anILLFORMEDFACEIsRefusedByNameRatherThanEvaluated(testCase)
        % Every refusal here names the PRODUCER's defect, because none of these can arise from a
        % well-formed AlgCon -- which is exactly why they must not be allowed to return a number.
            Pz = zeros(3,3,3);  Pz(1,1,3) = 1;
            testCase.verifyError(@() algBranch(Pz, [0 0 0], 1), 'algBranch:shape');
            testCase.verifyError(@() algBranch(Pz, [0 0], 0), 'algBranch:branch');
            testCase.verifyError(@() algBranch(Pz, [0 0], 1.5), 'algBranch:branch');
            testCase.verifyError(@() algBranch([], [0 0], 1), 'algBranch:polynomial');
            testCase.verifyError(@() algBranch(zeros(2,2,2,2), [0 0], 1), 'algBranch:polynomial');
            % degree 0 in z: {p = 0} is a condition on (x,y), not a value
            testCase.verifyError(@() algBranch(ones(2,2), [0 0], 1), 'algBranch:noZ');
        end

        function aVANISHINGFIBREIsNaNAndNotAnError(testCase)
        % Where p(x,y,.) is identically zero every z satisfies it, so the face names no value --
        % a legitimate question with no answer. It is reported as NaN with a zero root count,
        % NOT raised: the assembly needs to see it to locate the defect, and stopping hides it.
            P = zeros(2,1,3);
            P(2,1,3) = 1;  P(2,1,1) = -1;          % x*z^2 - x = x*(z^2-1): identically 0 at x = 0
            [z, info] = algBranch(P, [0 0; 1 0], 1);
            testCase.verifyTrue(isnan(z(1)), 'the whole fibre satisfies p there');
            testCase.verifyEqual(info.nreal(1), 0);
            testCase.verifyEqual(z(2), -1, 'AbsTol', 1e-12, 'and the ordinary point still answers');
        end

    end
end
