classdef QuaParCPLQUnitTest < matlab.unittest.TestCase
% Unit tests for `QuaParCPLQ`'s OPERATOR SURFACE -- eval, scalarMul, negate, addQuadratic,
% addScaledEnergy, add.
%
% WHY. QuaParCPLQ sat at 3.9% -- two lines of fifty-one -- after every bucket, fast, normal and
% slow alike. It exists so that a Case C conjugate composes with the generic operators
% infConv/moreau/proxAverage call by plain function name, which means MATLAB dispatches to it
% without any of those files knowing it exists. A class reached only by dispatch, and never
% directly by a test, is the easiest place in the toolbox for a wrong sign to live.
%
% WHAT IS ASSERTED. These operators are an ALGEBRA, so the assertions are the algebra's own
% identities, evaluated pointwise:
%       eval(scalarMul(f,a))     = a * eval(f)
%       eval(negate(f))          = -eval(f)
%       eval(addQuadratic(f,Q))  = eval(f) + q
%       eval(addScaledEnergy(f)) = eval(f) + alpha|x|^2
%       eval(add(f,g))           = eval(f) + eval(g)
% Nothing pins a cell list or an expression: every check is a value at a point, which is what a
% caller of infConv actually consumes.
%
% NOT COVERED HERE: `conj`. It drives the full symbolic conjugateOfPiecePoly/mergeL/addEq recipe
% and the class header records a known open defect in it (region.getNormalConeVertexQ). A test
% there would be a research item on an unfixed bug, not a unit test of this class.
%
% BUCKET: fast. The fixture is two affine pieces on two half-planes; every operator above is
% symbolic-expression arithmetic with no solve.

    properties (Constant)
        X = sym('x')
        Y = sym('y')
    end

    methods (Static)
        function v = vars()
            v = [QuaParCPLQUnitTest.X, QuaParCPLQUnitTest.Y];
        end

        function f = twoHalfPlanes(exprLeft, exprRight)
        % A full-plane piecewise function: exprLeft on {x <= 0}, exprRight on {x >= 0}. Two
        % pieces that TILE R^2 is the shape `add` is documented to be correct for (both operands'
        % pieces jointly covering the same region), and the shape a conjugate of a
        % compact-domain primal actually has.
            v = QuaParCPLQUnitTest.vars();
            rL = region(v(1), v);                      % x <= 0
            rR = region(-v(1), v);                     % -x <= 0
            fnd = [functionNDomain(symbolicFunction(exprLeft),  rL), ...
                   functionNDomain(symbolicFunction(exprRight), rR)];
            f = QuaParCPLQ(fnd);
        end

        function P = probes()
        % Points away from x = 0, where both pieces claim the boundary and which piece answers is
        % not part of any contract.
            P = [-2 1; -0.5 -3; -1 0; 3 2; 0.75 -1.5; 2 0];
        end
    end

    methods (Test)

        function evalReadsThePieceThatOwnsThePoint(testCase)
        % The base case everything below is asserted through: eval must return the value of the
        % piece whose region contains the point, and report a nonzero piece index for it.
            v = QuaParCPLQUnitTest.vars();
            f = QuaParCPLQUnitTest.twoHalfPlanes(-v(1) + 2*v(2), 3*v(1) - v(2) + 1);
            P = QuaParCPLQUnitTest.probes();
            for i = 1:size(P,1)
                x = P(i,1); y = P(i,2);
                if x < 0, want = -x + 2*y; else, want = 3*x - y + 1; end
                [got, idx] = f.eval(P(i,:));
                testCase.verifyEqual(got, want, 'AbsTol', 1e-12, sprintf( ...
                    'eval at (%g,%g)', x, y));
                testCase.verifyGreaterThan(idx, 0, sprintf( ...
                    'eval at (%g,%g) found no piece, but the two pieces tile the plane', x, y));
            end
        end

        function scalarMulAndNegateScaleTheValue(testCase)
            v = QuaParCPLQUnitTest.vars();
            f = QuaParCPLQUnitTest.twoHalfPlanes(-v(1) + 2*v(2), 3*v(1) - v(2) + 1);
            P = QuaParCPLQUnitTest.probes();
            for a = [2, 0.5, -3]
                g = scalarMul(f, a);
                for i = 1:size(P,1)
                    testCase.verifyEqual(g.eval(P(i,:)), a * f.eval(P(i,:)), 'AbsTol', 1e-12, ...
                        sprintf('scalarMul by %g at (%g,%g)', a, P(i,1), P(i,2)));
                end
            end
            h = negate(f);
            for i = 1:size(P,1)
                testCase.verifyEqual(h.eval(P(i,:)), -f.eval(P(i,:)), 'AbsTol', 1e-12, sprintf( ...
                    'negate at (%g,%g)', P(i,1), P(i,2)));
            end
        end

        function addQuadraticAddsTheQuadraticEverywhere(testCase)
        % A full-domain quadratic is finite everywhere, so it must be added to EVERY piece --
        % including the unbounded one. A version that touched only the first piece would pass any
        % check made at a single point, so both half-planes are probed.
            v = QuaParCPLQUnitTest.vars();
            f = QuaParCPLQUnitTest.twoHalfPlanes(-v(1) + 2*v(2), 3*v(1) - v(2) + 1);
            A = [2 1; 1 4]; b = [-1; 0.5]; c = 3;
            g = addQuadratic(f, A, b, c);
            P = QuaParCPLQUnitTest.probes();
            for i = 1:size(P,1)
                z = P(i,:).';
                q = 0.5*(z.'*A*z) + b.'*z + c;
                testCase.verifyEqual(g.eval(P(i,:)), f.eval(P(i,:)) + q, 'AbsTol', 1e-10, ...
                    sprintf('addQuadratic at (%g,%g)', z(1), z(2)));
            end
        end

        function addScaledEnergyIsTheQuadraticItClaims(testCase)
        % addScaledEnergy(alpha) is defined as addQuadratic(2*alpha*I, 0, 0), i.e. + alpha|x|^2.
        % Asserted against that closed form rather than against addQuadratic, so a change to the
        % factor of 2 -- the easiest slip in the file -- is caught.
            v = QuaParCPLQUnitTest.vars();
            f = QuaParCPLQUnitTest.twoHalfPlanes(v(2), v(1));
            P = QuaParCPLQUnitTest.probes();
            for alpha = [1, 0.25, 3]
                g = addScaledEnergy(f, alpha);
                for i = 1:size(P,1)
                    want = f.eval(P(i,:)) + alpha*(P(i,1)^2 + P(i,2)^2);
                    testCase.verifyEqual(g.eval(P(i,:)), want, 'AbsTol', 1e-10, sprintf( ...
                        'addScaledEnergy(%g) at (%g,%g)', alpha, P(i,1), P(i,2)));
                end
            end
        end

        function addIsPointwiseAddition(testCase)
        % Two full-plane tilings, which is the case `add` documents itself as correct for. The
        % operands deliberately have DIFFERENT breakpoints in their values so that a result which
        % kept only one operand's pieces would be wrong somewhere.
            v = QuaParCPLQUnitTest.vars();
            f = QuaParCPLQUnitTest.twoHalfPlanes(-v(1) + 2*v(2), 3*v(1) - v(2) + 1);
            g = QuaParCPLQUnitTest.twoHalfPlanes(v(1) + v(2) - 4, -2*v(2) + 5);
            h = add(f, g);
            P = QuaParCPLQUnitTest.probes();
            for i = 1:size(P,1)
                [want1] = f.eval(P(i,:));
                [want2] = g.eval(P(i,:));
                [got, idx] = h.eval(P(i,:));
                testCase.verifyGreaterThan(idx, 0, sprintf( ...
                    'add left (%g,%g) uncovered, but both operands cover the plane', P(i,1), P(i,2)));
                testCase.verifyEqual(got, want1 + want2, 'AbsTol', 1e-10, sprintf( ...
                    'add at (%g,%g)', P(i,1), P(i,2)));
            end
        end

        function addRefusesAForeignClassByName(testCase)
        % The header is explicit that mixed composition is NOT supported and must error rather
        % than silently do the wrong thing. Asserted by identifier, so the refusal is a contract
        % and not an accident of whatever MATLAB happens to raise.
            v = QuaParCPLQUnitTest.vars();
            f = QuaParCPLQUnitTest.twoHalfPlanes(v(1), v(2));
            other = QuaPol([1 0 1 0 0 0]);             % a full-domain quadratic, a different class
            testCase.verifyError(@() add(f, other), 'QuaParCPLQ:add:unsupportedType', ...
                'mixing a QuaParCPLQ with another mesh class must be refused by name');
        end
    end
end
