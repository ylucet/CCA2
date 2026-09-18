classdef bridgeFaceAlgTest < matlab.unittest.TestCase
% Tests for bridgeFaceAlg -- the producer that fills AlgCon.FP/.FB.
%
% BUCKET: fast. Everything here is closed-form numerics on 2x2 data; no symbolic call, no solver.
%
% THE ORACLE IS THE DEFINITION, AND IT SHARES NO ARITHMETIC WITH THE SUBJECT. The face is
%
%       z(x) = max { <s,x> - g1(s) : g1(s) = g2(s) } ,
%
% so the oracle maximises exactly that, over a parameterisation of the tie curve built from the
% quadratics themselves -- a LINE solved in closed form when the two Hessians agree, an ELLIPSE
% swept in theta and refined with fminbnd otherwise. bridgeFaceAlg reaches the same number through
% a discriminant in the Lagrange multiplier and a Sturm count; the two meet only at the answer.
%
% WHAT IS PINNED HERE, beyond "it runs":
%   - p VANISHES on constructed face points and does NOT vanish just off them. A polynomial that
%     vanishes everywhere would pass the first half alone, which is why the second half exists.
%   - the REGION equals the set where the maximiser is on the tie, against an independent
%     classifier that maximises each piece separately.
%   - the value is the TOP branch, and -- the finding that forced the negative index -- the
%     ASCENDING index of that same value is NOT constant on the face. That last test is the one
%     that goes red if anyone "tidies" FB back to a positive number.

    methods (Test)

        % ---- the polynomial ----------------------------------------------------------------

        function theImplicitPolynomialVanishesOnTheFaceAndNotBesideIt(testCase)
        % THE TWO HALVES ARE ASKED IN DIFFERENT CURRENCIES, and that is deliberate.
        %   ON the face, p is EVALUATED: the residual relative to the sum of its own term
        %   magnitudes is the sharpest statement available, and it is about p alone.
        %   OFF the face, a residual proves nothing -- the leading coefficient in z can be small
        %   and then |p| is small far from the surface too (measured: 1.1e-9 a full 0.1 above it).
        %   So that half asks about ROOTS: is z + 0.1 clear of every root of p(x,y,.)? It is a
        %   question about the surface, so it is asked of the surface.
        % The root DISTANCE is not used for the first half because it would then be measuring
        % `roots`, whose accuracy on a quartic with 1e6 coefficients is about 1e-8 -- the
        % root-finder's error, not the polynomial's.
            for f = fixtures()
                fx = f{1}; %#ok<*PFBNS>
                face = bridgeFaceAlg(fx.n1, fx.d1, fx.n2, fx.d2);
                [onFace, offFace] = faceResiduals(fx, face.FP, 0.1);
                testCase.verifyLessThan(onFace, 1e-12, ...
                    sprintf('%s: p must vanish at every constructed face point', fx.name));
                testCase.verifyGreaterThan(offFace, 0.05, ...
                    sprintf('%s: 0.1 above the face must be clear of every root', fx.name));
            end
        end

        function theDegenerateTieTakesTheQUADRATICPathInTheMultiplier(testCase)
        % Two pieces with the SAME Hessian tie along a LINE, not a conic, and the mu^3 coefficient
        % is half the determinant of the tie's own 3x3 matrix -- zero exactly then. This is the
        % common case, not an exotic one: any two translates of one quadratic are such a pair, so
        % the quadratic branch of the discriminant has to work.
            fx = fixtureByName('sameHessian');
            [~, info] = bridgeFaceAlg(fx.n1, fx.d1, fx.n2, fx.d2);
            testCase.verifyEqual(info.muDegree, 2, ...
                'a degenerate tie must drop Psi to degree 2 in mu');
        end

        % ---- the region --------------------------------------------------------------------

        function theRegionIsWhereTheMaximiserSitsOnTheTie(testCase)
            for f = fixtures()
                fx = f{1};
                face = bridgeFaceAlg(fx.n1, fx.d1, fx.n2, fx.d2);
                [agree, disagree] = regionAgreement(fx, face, 1500);
                testCase.verifyGreaterThan(agree, 100, ...
                    sprintf('%s: too few points to conclude anything', fx.name));
                testCase.verifyEqual(disagree, 0, ...
                    sprintf('%s: the H-form and the maximiser disagree somewhere', fx.name));
            end
        end

        % ---- the value, and the branch ------------------------------------------------------

        function theValueIsTheTopBranchOfTheImplicitPolynomial(testCase)
            for f = fixtures()
                fx = f{1};
                face = bridgeFaceAlg(fx.n1, fx.d1, fx.n2, fx.d2);
                testCase.verifyEqual(face.FB, -1, ...
                    sprintf('%s: a supremum is the top branch', fx.name));
                [worst, n] = valueAgreement(fx, face);
                testCase.verifyGreaterThan(n, 50, ...
                    sprintf('%s: too few region points to conclude', fx.name));
                testCase.verifyLessThan(worst, 1e-9, ...
                    sprintf('%s: algBranch(FP,-1) must reproduce the definition', fx.name));
            end
        end

        function theASCENDINGIndexCannotNameThisFace(testCase)
        % THE MEASUREMENT THAT FORCED THE NEGATIVE INDEX, pinned as a test rather than left in a
        % document. Two real roots can collide and LEAVE the real line inside the face, which drops
        % the count by two and shifts every index above the departed pair. So the value's ascending
        % index takes more than one value on ONE face, while its distance from the TOP never
        % changes. Anything that re-spells FB as a positive index makes this red.
            fx = fixtureByName('anisotropic');
            face = bridgeFaceAlg(fx.n1, fx.d1, fx.n2, fx.d2);
            idx = ascendingIndexOfValue(fx, face);
            testCase.verifyGreaterThan(numel(unique(idx)), 1, ...
                'the ascending index must be seen to VARY -- that is the whole reason FB is -1');
            testCase.verifyEqual(unique(idx(:)) > 0, true(numel(unique(idx)),1), ...
                'and every one of them is a real index, not a miss');
        end

        % ---- the deliverable: it fills an AlgCon -----------------------------------------------

        function theFaceFillsAnAlgConAndEvaluatesFromItAlone(testCase)
        % What the item was FOR. AlgCon has carried FP/FB since 2026-09-13 with nothing to put in
        % them; this is the first object built from a produced face, and its value is read back
        % from the payload alone -- no formula, no provenance.
            fx = fixtureByName('isotropic');
            face = bridgeFaceAlg(fx.n1, fx.d1, fx.n2, fx.d2);
            a = AlgCon();
            a.FP = {face.FP};
            a.FB = face.FB;
            testCase.verifyEqual(a.FB, -1);
            x = face.wit(1:2) / face.wit(3);
            got  = algBranch(a.FP{1}, x(:).', a.FB);
            want = maxOverTie(fx, x(:));
            testCase.verifyEqual(got, want, 'RelTol', 1e-9, ...
                'the value must come back from FP/FB alone');
        end


        % ---- the CUT: the face restricted to the arc inside its own cell pair ------------------

        function theCUTRestrictsTheFaceToTheArcInsideThePair(testCase)
        % THE ORACLE CHANGES WITH THE PAIR, which is the whole point. Without P the face is the max
        % over the WHOLE tie; with P it is the max over the tie points IN P, and for a line tie both
        % are closed form -- the objective is a concave quadratic along the line, so the restricted
        % maximiser is the stationary one CLIPPED to the interval. Nothing is sampled.
            fx = fixtureByName('sameHessian');
            cases = { 'quadrant', [1 0 0; 0 1 0], 1
                      'square',   [1 0 1; -1 0 1; 0 1 1; 0 -1 1], 2 };
            for c = 1:size(cases,1)
                P = cases{c,2};
                [face, info] = bridgeFaceAlg(fx.n1, fx.d1, fx.n2, fx.d2, P);
                testCase.verifyEqual(info.cut, 'line', sprintf('%s: a line tie', cases{c,1}));
                testCase.verifyEqual(info.nCut, cases{c,3}, ...
                    sprintf('%s: one half-plane per FINITE arc endpoint', cases{c,1}));
                [worst, n] = restrictedAgreement(fx, face, P);
                testCase.verifyGreaterThan(n, 50, sprintf('%s: too few region points', cases{c,1}));
                testCase.verifyLessThan(worst, 1e-9, ...
                    sprintf('%s: the cut face must equal the RESTRICTED definition', cases{c,1}));
            end
        end

        function theUNCUTFaceCLAIMSGroundTheCutOneDoesNot(testCase)
        % The cut has to be seen to REMOVE something, or the test above would pass on a cut that
        % does nothing. It also says why the cut is not optional: on the ground it gives up, the
        % uncut face is claiming a value the restricted problem cannot attain.
            fx = fixtureByName('sameHessian');
            P  = [1 0 1; -1 0 1; 0 1 1; 0 -1 1];
            full = bridgeFaceAlg(fx.n1, fx.d1, fx.n2, fx.d2);
            cut  = bridgeFaceAlg(fx.n1, fx.d1, fx.n2, fx.d2, P);
            ctr = full.wit(1:2)/full.wit(3);
            rng(23);
            dropped = 0;  overclaim = 0;
            for i = 1:2000
                x = ctr(:) + (rand(2,1)*2-1)*5;
                if ~inHForm(full, x) || inHForm(cut, x), continue, end
                dropped = dropped + 1;
                want = restrictedMax(fx, P, x);
                got  = algBranch(full.FP, x(:).', full.FB);
                if isfinite(want) && got > want + 1e-9*max(1,abs(want)), overclaim = overclaim + 1; end
            end
            testCase.verifyGreaterThan(dropped, 20, 'the cut must actually remove region');
            testCase.verifyEqual(overclaim, dropped, ...
                'and on every point it removes, the UNCUT face overstates the restricted sup');
        end

        function EVERYTieKindIsCUTCorrectlyByThePredicate(testCase)
        % THE ITEM, FINISHED. For every kind of tie the pipeline produces -- line, ellipse,
        % hyperbola -- and for two shapes of cell pair, the face must agree with the RESTRICTED
        % definition (the max over the tie points INSIDE the pair) and must have no HOLES and keep
        % nothing whose base lies outside. Values and membership are different questions and both
        % are asked; a cut that is too small passes the first alone.
        %
        % This is what two constructions over sign conditions could not do. The difference is that
        % the bounded conditions live in the DUAL plane -- mu in [-1,0] and s in a polygon -- where
        % everything is straight, instead of being forced into curves in x.
        %
        % THE LAST TWO FIXTURES ARE THE DEGENERATE CELLS, which are the MAJORITY on real input: a
        % conjugate has one full-rank cell, one rank-1 cell per edge and one affine cell per vertex.
        % Both used to be refused outright, because their gradient map is singular and its end of
        % the sweep is a LINE or a POINT rather than a conic. With membership decided from the
        % geometry that only costs a row of the outer bound, so they now answer.
            for nm = {'sameHessian', 'isotropic', 'anisotropic', 'hyperbolic', ...
                      'vertexCell', 'edgeCell'}
                fx = fixtureByName(nm{1});
                for Pc = {[1 0 0; 0 1 0], [1 0 1; -1 0 3; 0 1 2; 0 -1 3]}
                    P = Pc{1};
                    face  = bridgeFaceAlg(fx.n1, fx.d1, fx.n2, fx.d2, P);
                    full  = bridgeFaceAlg(fx.n1, fx.d1, fx.n2, fx.d2);
                    [worst, nv] = restrictedAgreementAny(fx, face, P);
                    testCase.verifyGreaterThan(nv, 15, sprintf('%s: too few region points', nm{1}));
                    testCase.verifyLessThan(worst, 1e-9, ...
                        sprintf('%s: the cut face must equal the RESTRICTED definition', nm{1}));
                    [inside, holes, over] = baseMembershipAny(fx, face, full, P);
                    testCase.verifyGreaterThan(inside, 15, sprintf('%s: too few base-inside', nm{1}));
                    testCase.verifyEqual(holes, 0, sprintf('%s: a point whose base IS in the pair was dropped', nm{1}));
                    testCase.verifyEqual(over, 0, sprintf('%s: a point whose base is OUTSIDE was kept', nm{1}));
                end
            end
        end

        function theROWBasedCutStillREFUSESAConicTieAndTheProducerCatchesIt(testCase)
        % The division of labour, pinned. `bridgeArcCut` produces EXACT cut rows and can only do so
        % for a LINE tie -- a conic tie's cut boundary is a segment, which no sign condition on a
        % whole curve expresses, and the two constructions that tried are refuted in its header. It
        % keeps refusing; the producer catches that and leaves the sign conditions as an OUTER
        % BOUND, because membership does not read them.
            fx = fixtureByName('isotropic');
            d  = 2;  N1 = [2 0; 0 2];  m1 = [0;0];  N2 = [4 0; 0 4];  m2 = [-4; 2];
            tie = ratQ.diffConic(fx.n1, fx.d1, fx.n2, fx.d2);
            testCase.verifyError(@() bridgeArcCut(tie, [1 0 0], N1, m1, N2, m2, d, [1 1 1], ...
                                                  [0 0 0 1 0 0], [0 0 0 1 0 -1], 0), ...
                'bridgeArcCut:conicTie');
            % ...and the producer still answers, because the face is not its bound
            face = bridgeFaceAlg(fx.n1, fx.d1, fx.n2, fx.d2, [1 0 0; 0 1 0]);
            testCase.verifyTrue(bridgeMember(face.pair, face.wit(1:2)/face.wit(3)), ...
                'the witness must be IN the face the predicate describes');
        end

        function theSIGNConditionsAreAnOUTERBoundAndNotTheFace(testCase)
        % `cells` is kept, and what it is worth stating: it CONTAINS the face, so it is a legitimate
        % cheap pre-filter, and for a conic tie it is strictly larger -- which is exactly why
        % membership does not read it.
            fx = fixtureByName('isotropic');
            P  = [1 0 0; 0 1 0];
            face = bridgeFaceAlg(fx.n1, fx.d1, fx.n2, fx.d2, P);
            rng(71);
            ctr = face.wit(1:2).'/face.wit(3);
            inBoth = 0;  outside = 0;  strictlyLarger = 0;
            for i = 1:1500
                x = ctr + (rand(2,1)*2-1)*5;
                inF = inHForm(face, x);
                inC = inCells(face, x);
                if inF && inC, inBoth = inBoth + 1; end
                if inF && ~inC, outside = outside + 1; end
                if inC && ~inF, strictlyLarger = strictlyLarger + 1; end
            end
            testCase.verifyGreaterThan(inBoth, 20, 'the test must actually reach the face');
            testCase.verifyEqual(outside, 0, 'no point of the face may fall outside the bound');
            testCase.verifyGreaterThan(strictlyLarger, 0, ...
                'and for a conic tie the bound is STRICTLY larger -- the reason it is not the face');
        end

        function aPAIRThatExcludesTheWholeTieIsRefusedByName(testCase)
        % A facet with the entire tie on its far side means the pair holds no tie point at all, so
        % it carries no bridge and should never have been classified as one. That is the caller's
        % defect and it is named as such.
            fx = fixtureByName('sameHessian');
        % The row is PARALLEL to the tie and on its far side -- which is what "excludes the whole
        % tie" has to mean here, since any row that crosses the line merely moves the arc, however
        % far away it puts it.
            testCase.verifyError(@() bridgeFaceAlg(fx.n1, fx.d1, fx.n2, fx.d2, [-3 1 -5]), ...
                'bridgeFaceAlg:noWitness');
            testCase.verifyError(@() bridgeArcCut(ratQ.conic([0 0 0 3 -1 -2]), [-3 1 -5], ...
                                                  eye(2), [0;0], eye(2), [1;1], 1, [0 0 1], ...
                                                  [0 0 0 1 0 0], [0 0 0 1 0 -1], 0), ...
                'bridgeArcCut:emptyArc');
        end

        function bridgeArcCutCHECKSItsOwnArguments(testCase)
            testCase.verifyError(@() bridgeArcCut([0 0 0 1 1 0], [1 0 0]), 'bridgeArcCut:nargin');
            testCase.verifyError(@() bridgeArcCut(ratQ.conic([0 0 0 0 0 1]), [1 0 0], ...
                                                  eye(2), [0;0], eye(2), [1;1], 1, [0 0 1], ...
                                                  [0 0 0 1 0 0], [0 0 0 1 0 -1], 0), ...
                'bridgeArcCut:degenerateTie');
        end


        function theCUTHasNOHOLESAndKeepsNothingItShouldNot(testCase)
        % THE DIRECTION THE VALUE CHECKS CANNOT SEE. A cut that is too SMALL agrees with the
        % definition at every point it keeps -- it simply answers on less ground -- so agreement is
        % only half the statement. The other half is membership, and it is asked of the BASE POINT:
        %
        %     x belongs to this face  <=>  x is in the uncut sweep AND its base s(x) lies in P.
        %
        % The base of a point of the uncut sweep is the maximiser over the WHOLE tie (measured over
        % 4758 points when the branch index was settled), so both directions are checkable without
        % touching the cut's own arithmetic. This is the test that caught the conic-tie construction
        % being wrong; the line tie passes it in both directions.
            fx = fixtureByName('sameHessian');
            cases = { [1 0 1; -1 0 1; 0 1 1; 0 -1 1], 'square'
                      [1 0 0; 0 1 0],                 'quadrant'
                      [0 1 1; 0 -1 2],                'strip' };
            for c = 1:size(cases,1)
                P = cases{c,1};
                face = bridgeFaceAlg(fx.n1, fx.d1, fx.n2, fx.d2, P);
                full = bridgeFaceAlg(fx.n1, fx.d1, fx.n2, fx.d2);
                [inside, holes, overKept] = baseMembership(fx, face, full, P);
                testCase.verifyGreaterThan(inside, 100, ...
                    sprintf('%s: too few points with a base inside the pair', cases{c,2}));
                testCase.verifyEqual(holes, 0, ...
                    sprintf('%s: the cut dropped a point whose base IS in the pair', cases{c,2}));
                testCase.verifyEqual(overKept, 0, ...
                    sprintf('%s: the cut kept a point whose base is OUTSIDE the pair', cases{c,2}));
            end
        end

        % ---- what it refuses, by name ------------------------------------------------------

        function aSINGULARGradientMapOMITSABoundInsteadOfRefusing(testCase)
        % WHAT USED TO BE TWO REFUSALS. A VERTEX cell's dual face is AFFINE, so grad g is constant
        % and its end of the sweep is a POINT; an EDGE cell's Hessian has rank 1, so grad g collapses
        % the plane onto a LINE. Neither end is a conic, so neither contributes a row to the outer
        % bound -- and since membership is decided from the pair's own data, that costs a row of the
        % BOUND and nothing of the FACE. These two cell kinds are the majority on real input.
            for nm = {'vertexCell', 'edgeCell'}
                fx = fixtureByName(nm{1});
                face = bridgeFaceAlg(fx.n1, fx.d1, fx.n2, fx.d2);
                testCase.verifyEqual(size(face.cells{1},1), 1, ...
                    sprintf('%s: one end of the sweep is not a conic, so one row is missing', nm{1}));
                testCase.verifyEqual(face.FB, -1, sprintf('%s: still the top branch', nm{1}));
                testCase.verifyTrue(bridgeMember(face.pair, face.wit(1:2)/face.wit(3)), ...
                    sprintf('%s: and the face is still a face', nm{1}));
            end
        end

        function aTANGENTIALPairCarriesNOTwoDimensionalFaceAndSaysSo(testCase)
        % NOT A GAP, AN ANSWER. The pair walk counts a TANGENTIAL tie as a bridge -- and tangential
        % ties are the MAJORITY (DECISIONS.md 2026-09-13) -- but where the tie only TOUCHES the
        % pair the arc is a single point, so the sweep is a single segment with no area. Measured
        % on the standard fixtures: 10 of 29 bridge pairs are of this kind. An assembly must skip
        % them and take those pairs' faces from the two sides alone, which is why the refusal says
        % what it means instead of reporting a failed search.
        %
        % Here: g1 = |s|^2 against the constant 1, so the tie is the unit circle, and the pair asks
        % for s1 >= 1 -- which the circle reaches at exactly one point.
            testCase.verifyError(@() bridgeFaceAlg([0 0 0 0 2 0 2 0 0 0], 1, ...
                                                   [0 0 0 0 0 0 0 0 0 1], 1, [1 0 -1]), ...
                'bridgeFaceAlg:tangentialArc');
            % ...while the same pair with room to cross carries a perfectly good face
            face = bridgeFaceAlg([0 0 0 0 2 0 2 0 0 0], 1, [0 0 0 0 0 0 0 0 0 1], 1, [1 0 0]);
            testCase.verifyTrue(bridgeMember(face.pair, face.wit(1:2)/face.wit(3)));
        end

        function TWOFacesThatCannotTIEAreRefusedByName(testCase)
        % Same Hessian AND same linear part: g1 - g2 is a CONSTANT, so { g1 = g2 } is empty or the
        % whole plane. There is no tie curve and so no bridge face, and Psi degenerates in mu.
            testCase.verifyError(@() bridgeFaceAlg([0 0 0 0 1 0 1 2 3 0], 1, ...
                                                   [0 0 0 0 1 0 1 2 3 5], 1), ...
                'bridgeFaceAlg:straightTie');
        end

        function aTIEWithNoREALPointsIsRefusedByName(testCase)
        % Two concentric quadratics that never take the same value: the tie conic is imaginary, so
        % there is no curve to sweep and no bridge. The pair filter would not have sent this here,
        % which is why the refusal names the caller's mistake rather than guessing a region.
            testCase.verifyError(@() bridgeFaceAlg([0 0 0 0 1 0 1 0 0 0], 1, ...
                                                   [0 0 0 0 2 0 2 0 0 5], 1), ...
                'bridgeFaceAlg:noWitness');
        end

        function AN2POWER53CeilingIsANamedREFUSALAndNotARawOverflow(testCase)
        % The one input class this producer cannot answer, and it must say so in its own words. The
        % discriminant multiplies four coefficients that already carry both faces' denominators
        % cubed, so a pair with large data puts the FACE'S OWN POLYNOMIAL past 2^53 -- the answer
        % being large, not an intermediate being clumsy. A raw `ratQ:overflow` escaping here would
        % look exactly like the routine being wrong.
            testCase.verifyError(@() bridgeFaceAlg([0 0 0 0 8 -2 4 -4 2 2], 2, ...
                                                   [0 0 0 0 6 -2 3 -5 -3 -7], 7), ...
                'bridgeFaceAlg:polynomialOverflow');
        end

        function aCUBICFaceIsRefusedByName(testCase)
        % A conjugate face is quadratic. A cubic here means the caller passed something that is not
        % a conjugate cell, and that is worth saying at the boundary.
            testCase.verifyError(@() bridgeFaceAlg([1 0 0 0 1 0 1 0 0 0], 1, ...
                                                   [0 0 0 0 2 0 2 1 1 1], 1), ...
                'bridgeFaceAlg:cubicFace');
        end

        function theARGUMENTSAreChecked(testCase)
            testCase.verifyError(@() bridgeFaceAlg([0 0 0 0 1 0 1 0 0 0], 1), ...
                'bridgeFaceAlg:nargin');
            testCase.verifyError(@() bridgeFaceAlg([1 2 3], 1, [0 0 0 0 1 0 1 0 0 0], 1), ...
                'bridgeFaceAlg:shape');
            testCase.verifyError(@() bridgeFaceAlg([0 0 0 0 1 0 1 0 0 0], 0, ...
                                                   [0 0 0 0 2 0 2 1 1 1], 1), ...
                'bridgeFaceAlg:denominator');
            testCase.verifyError(@() bridgeFaceAlg([0 0 0 0 1 0 1 0 0 0], 1, ...
                                                   [0 0 0 0 2 0 2 1 1 1], 1, [1 0]), ...
                'bridgeFaceAlg:pairShape');
        end
    end
end

% ==================================================================================================
% fixtures
% ==================================================================================================
function F = fixtures()
    F = {fixtureByName('isotropic'), fixtureByName('anisotropic'), ...
         fixtureByName('bothShifted'), fixtureByName('sameHessian')};
end

function fx = fixtureByName(name)
    switch name
        case 'isotropic'                        % A2 = 2*A1: the concentric-circle instance
            fx = mk(name, [0 0 0 0 1 0 1 0 0 0], 1, [0 0 0 0 4 0 4 -4 2 3], 2);
        case 'anisotropic'                      % a genuinely different Hessian: degree 4 in z
            fx = mk(name, [0 0 0 0 1 0 1 0 0 0], 1, [0 0 0 0 3 1 2 -4 2 3], 1);
        case 'bothShifted'                      % neither piece centred, both denominators > 1
            fx = mk(name, [0 0 0 0 4 1 3 2 5 -2], 3, [0 0 0 0 2 0 2 -1 1 1], 1);
        case 'vertexCell'                       % g2 AFFINE: its dual cell is a VERTEX cell
            fx = mk(name, [0 0 0 0 2 0 2 0 0 0], 1, [0 0 0 0 0 0 0 -3 1 2], 1);
        case 'edgeCell'                         % g2's Hessian has RANK 1: an EDGE cell
            fx = mk(name, [0 0 0 0 2 0 2 0 0 0], 1, [0 0 0 0 1 0 0 -3 1 2], 1);
        case 'hyperbolic'                       % Delta indefinite: the tie is a hyperbola
            fx = mk(name, [0 0 0 0 3 0 1 0 0 0], 1, [0 0 0 0 1 0 3 -4 2 3], 1);
        case 'sameHessian'                      % the tie is a LINE: the degenerate-conic path
            fx = mk(name, [0 0 0 0 2 0 2 0 0 0], 1, [0 0 0 0 2 0 2 -3 1 2], 1);
        otherwise
            error('bridgeFaceAlgTest:fixture', 'no fixture named %s', name);
    end
end

function fx = mk(name, n1, d1, n2, d2)
    fx = struct('name', name, 'n1', n1, 'd1', d1, 'n2', n2, 'd2', d2, ...
                'A1', [n1(5) n1(6); n1(6) n1(7)]/d1, 'b1', [n1(8); n1(9)]/d1, 'k1', n1(10)/d1, ...
                'A2', [n2(5) n2(6); n2(6) n2(7)]/d2, 'b2', [n2(8); n2(9)]/d2, 'k2', n2(10)/d2);
end

% ==================================================================================================
% the oracle: the tie curve, parameterised from the two quadratics
% ==================================================================================================
function [s, tangent] = tiePoint(fx, t)
% A point of { g1 = g2 } at parameter t, and the parameterisation's own derivative is not needed --
% the maximisations below are over t directly.
    [Dl, dl, kp] = tieData(fx);
    if all(Dl(:) == 0)
        % a LINE: dl's + kp = 0
        s0  = -kp * dl / (dl.'*dl);
        dir = [-dl(2); dl(1)];
        s = s0 + t*dir;  tangent = dir;
        return
    end
    [V, L, R] = ellipseFrame(fx);
    s = tieCentre(fx) + V*[sqrt(R/L(1))*cos(t); sqrt(R/L(2))*sin(t)];
    tangent = V*[-sqrt(R/L(1))*sin(t); sqrt(R/L(2))*cos(t)];
end

function [Dl, dl, kp] = tieData(fx)
    Dl = fx.A1 - fx.A2;  dl = fx.b1 - fx.b2;  kp = fx.k1 - fx.k2;
end

function c0 = tieCentre(fx)
    [Dl, dl, ~] = tieData(fx);
    c0 = -Dl \ dl;
end

function [V, L, R] = ellipseFrame(fx)
% The tie in centred form: with u = s - c0, the tie is u'Dl u = R. When Dl is definite and R has
% its sign, that is an ellipse and (V, L) is the eigenframe of the matrix that makes it one.
    [Dl, dl, kp] = tieData(fx);
    c0 = tieCentre(fx);
    R  = -2*((c0.'*Dl*c0)/2 + dl.'*c0 + kp);
    M  = Dl;
    if R < 0, M = -Dl;  R = -R; end
    [V, D] = eig(M);
    L = diag(D);
    if any(L <= 0)
        error('bridgeFaceAlgTest:notAnEllipse', ...
            ['the oracle parameterises a LINE or an ELLIPSE, and this tie is neither. A ' ...
             'hyperbolic tie needs its own parameterisation rather than a silent sample.']);
    end
end

function z = maxOverTie(fx, x)
% max over the tie of <s,x> - g1(s): the definition, maximised over the parameterisation above.
    g1 = @(s) (s.'*fx.A1*s)/2 + fx.b1.'*s + fx.k1;
    obj = @(t) -objAt(fx, g1, t, x);                       % fminbnd MINIMISES
    [Dl, ~, ~] = tieData(fx);
    if all(Dl(:) == 0)
        % a line: the objective is a concave quadratic in t and its maximiser is closed form
        [s0, dir] = deal(tiePoint(fx, 0), lineDir(fx));
        a = -(dir.'*fx.A1*dir)/2;
        b = dir.'*x - (s0.'*fx.A1*dir + fx.b1.'*dir);
        if a >= 0, z = inf; return, end                    % unbounded along the tie
        t = -b/(2*a);
        z = objAt(fx, g1, t, x);
        return
    end
    ts = linspace(0, 2*pi, 721);
    vals = arrayfun(@(t) objAt(fx, g1, t, x), ts);
    [~, i] = max(vals);
    lo = ts(max(i-1,1));  hi = ts(min(i+1, numel(ts)));
    tb = fminbnd(obj, lo, hi, optimset('TolX', 1e-14));
    z = max(objAt(fx, g1, tb, x), vals(i));
end

function d = lineDir(fx)
    [~, dl, ~] = tieData(fx);
    d = [-dl(2); dl(1)];
end

function v = objAt(fx, g1, t, x)
    s = tiePoint(fx, t);
    v = s.'*x - g1(s);
end

% ==================================================================================================
% the checks the tests read
% ==================================================================================================
function [onFace, offFace] = faceResiduals(fx, FP, h)
% CONSTRUCTED face points: at a tie point s the whole segment between the two gradients lies on the
% face, with value <s,x> - g1(s). Nothing here solves anything, so it cannot agree with the subject
% by accident.
% [output] onFace  : the worst RELATIVE residual |p| at a constructed face point
%          offFace : the smallest distance from (that value + h) to ANY root of p(x,y,.)
    g1 = @(s) (s.'*fx.A1*s)/2 + fx.b1.'*s + fx.k1;
    [Dl, ~, ~] = tieData(fx);
    if all(Dl(:) == 0), ts = linspace(-3, 3, 9); else, ts = linspace(0, 2*pi, 9); end
    onFace = 0;  offFace = inf;
    for t = ts
        s = tiePoint(fx, t);
        for lam = [0 0.25 0.5 0.75 1]
            x = (1-lam)*(fx.A1*s + fx.b1) + lam*(fx.A2*s + fx.b2);
            z = s.'*x - g1(s);
            [~, info] = algBranch(FP, x(:).', -1);
            rr = info.allReal{1};
            if isempty(rr)
                onFace = inf;  offFace = 0;  return         % no root at all: a definite failure
            end
            onFace  = max(onFace,  relEval(FP, x, z));
            offFace = min(offFace, min(abs(rr - (z + h))));
        end
    end
end

function [agree, disagree] = regionAgreement(fx, face, n)
    ctr = face.wit(1:2)/face.wit(3);
    rng(11);
    agree = 0;  disagree = 0;
    for i = 1:n
        x = ctr(:) + (rand(2,1)*2-1)*5;
        if inHForm(face, x) == maximiserIsOnTheTie(fx, x)
            if inHForm(face, x), agree = agree + 1; end
        else
            disagree = disagree + 1;
        end
    end
end

function [worst, n] = valueAgreement(fx, face)
    ctr = face.wit(1:2)/face.wit(3);
    rng(13);
    worst = 0;  n = 0;
    for i = 1:600
        x = ctr(:) + (rand(2,1)*2-1)*5;
        if ~inHForm(face, x), continue, end
        want = maxOverTie(fx, x);
        if ~isfinite(want), continue, end
        got = algBranch(face.FP, x(:).', face.FB);
        n = n + 1;
        worst = max(worst, abs(got - want)/max(1, abs(want)));
    end
end

function idx = ascendingIndexOfValue(fx, face)
    ctr = face.wit(1:2)/face.wit(3);
    rng(17);
    idx = [];
    for i = 1:800
        x = ctr(:) + (rand(2,1)*2-1)*5;
        if ~inHForm(face, x), continue, end
        want = maxOverTie(fx, x);
        if ~isfinite(want), continue, end
        [~, info] = algBranch(face.FP, x(:).', -1);
        rr = info.allReal{1};
        if isempty(rr), continue, end
        [gap, j] = min(abs(rr - want));
        if gap > 1e-8*max(1, abs(want)), continue, end
        idx(end+1) = j; %#ok<AGROW>
    end
end

function [worst, n] = restrictedAgreement(fx, face, P)
    ctr = face.wit(1:2)/face.wit(3);
    rng(29);
    worst = 0;  n = 0;
    for i = 1:1200
        x = ctr(:) + (rand(2,1)*2-1)*4;
        if ~inHForm(face, x), continue, end
        want = restrictedMax(fx, P, x);
        if ~isfinite(want), continue, end
        got = algBranch(face.FP, x(:).', face.FB);
        n = n + 1;
        worst = max(worst, abs(got - want)/max(1, abs(want)));
    end
end

function v = restrictedMax(fx, P, x)
% max over { s on the tie AND s in P } of <s,x> - g1(s), in CLOSED FORM for a line tie: the
% objective is a concave quadratic along the line, so the answer is its stationary point clipped to
% the interval P leaves. No sampling, and no arithmetic shared with the subject.
    [Dl, dl, kp] = tieData(fx);
    if any(Dl(:) ~= 0)
        error('bridgeFaceAlgTest:notALineTie', 'this oracle is for a line tie');
    end
    w  = [dl(2); -dl(1)];
    s0 = -kp * dl / (dl.'*dl);
    lo = -inf;  hi = inf;
    for i = 1:size(P,1)
        a = P(i,1:2).';  c = P(i,3);
        A = a.'*w;  B = a.'*s0 + c;
        if A == 0
            if B < 0, v = -inf; return, end
            continue
        end
        if A > 0, lo = max(lo, -B/A); else, hi = min(hi, -B/A); end
    end
    if lo > hi, v = -inf; return, end
    alpha = w.'*fx.A1*w;
    if alpha <= 0, v = inf; return, end
    beta  = w.'*x - (w.'*fx.A1*s0 + fx.b1.'*w);
    gamma = s0.'*x - ((s0.'*fx.A1*s0)/2 + fx.b1.'*s0 + fx.k1);
    tau = min(max(beta/alpha, lo), hi);
    if ~isfinite(tau), v = inf; return, end
    v = -alpha/2*tau^2 + beta*tau + gamma;
end

function [inside, holes, overKept] = baseMembership(fx, face, full, P)
% Walk points of the UNCUT sweep, find each one's base on the tie, and compare "the base is in P"
% with "the H-form keeps it". Both disagreements are counted, and they are opposite defects.
    [Dl, dl, kp] = tieData(fx);
    w  = [dl(2); -dl(1)];
    s0 = -kp * dl / (dl.'*dl);
    T  = s0 + w*linspace(-200, 200, 40001);
    vq = 0.5*sum(T.*(fx.A1*T), 1) + fx.b1.'*T + fx.k1;
    ctr = face.wit(1:2).'/face.wit(3);
    rng(43);
    inside = 0;  holes = 0;  overKept = 0;
    for i = 1:3000
        x = ctr + (rand(2,1)*2-1)*4;
        if ~inHForm(full, x), continue, end
        [~, j] = max(sum(T.*x, 1) - vq);
        sB = T(:,j);
        slack = P(:,1:2)*sB + P(:,3);
        if all(slack > 1e-6)
            inside = inside + 1;
            if ~inHForm(face, x), holes = holes + 1; end
        elseif all(slack < -1e-6) && inHForm(face, x)
            overKept = overKept + 1;
        end
    end
end

function [worst, n] = restrictedAgreementAny(fx, face, P)
% Against the RESTRICTED definition, for ANY tie kind: sweep the tie, keep the points inside P, and
% refine the bracket in the curve's OWN parameter -- bracketing on a coordinate instead is what made
% an earlier version of this oracle report the subject's error as 1e-3.
    ctr = face.wit(1:2).'/face.wit(3);
    rng(73);
    worst = 0;  n = 0;
    for i = 1:500
        x = ctr + (rand(2,1)*2-1)*4;
        if ~inHForm(face, x), continue, end
        want = tieMax(fx, P, x);
        if ~isfinite(want), continue, end
        got = algBranch(face.FP, x(:).', face.FB);
        if ~isfinite(got), continue, end
        n = n + 1;
        worst = max(worst, abs(got - want)/max(1, abs(want)));
    end
end

function [inside, holes, overKept] = baseMembershipAny(fx, face, full, P)
    ctr = face.wit(1:2).'/face.wit(3);
    rng(79);
    inside = 0;  holes = 0;  overKept = 0;
    for i = 1:500
        x = ctr + (rand(2,1)*2-1)*4;
        if ~inHForm(full, x), continue, end
        [~, sB] = tieMax(fx, zeros(0,3), x);
        if isempty(sB), continue, end
        slack = P(:,1:2)*sB + P(:,3);
        if all(slack > 1e-6)
            inside = inside + 1;
            if ~inHForm(face, x), holes = holes + 1; end
        elseif all(slack < -1e-6) && inHForm(face, x)
            overKept = overKept + 1;
        end
    end
end

function [v, sBest] = tieMax(fx, P, x)
% max over { the tie, inside P } of <s,x> - g1(s), and where it is attained. Coarse sweep then four
% rounds of refinement, bracketed in the PARAMETER the curve is generated by.
    [Dl, dl, kp] = tieData(fx);
    isLine = all(Dl(:) == 0);
    if isLine, lo = -60; hi = 60; else, lo = -25; hi = 25; end
    v = -inf;  sBest = [];  best = NaN;
    % THE FIRST SWEEP MUST BE FINE ENOUGH TO FIND THE RIGHT BASIN; the refinements after it need
    % almost nothing. Making every round coarse was measured wrong: the search then converges to a
    % LOCAL maximum of a curve that has several, and reports the subject as being in error.
    % Rounds after the first cost ~200 evaluations each and contract the bracket 100-fold, so six
    % of them are free. They are NEEDED because the restricted maximiser often sits AT the arc's
    % endpoint, where the value's error is LINEAR in the bracket width rather than quadratic.
    res = [4001 201 201 201 201 201 201];
    for round = 1:numel(res)
        ss = linspace(lo, hi, res(round));
        [S, tt] = tiePoints(Dl, dl, kp, ss);
        if isempty(S), return, end
        val = sum(S.*x, 1) - (0.5*sum(S.*(fx.A1*S), 1) + fx.b1.'*S + fx.k1);
        if ~isempty(P), val(~all(P(:,1:2)*S + P(:,3) >= 0, 1)) = -inf; end
        [w, j] = max(val);
        if ~isfinite(w), return, end
        if w > v, v = w;  sBest = S(:,j);  best = tt(j); end
        h = (hi - lo)/(res(round) - 1);
        lo = best - 2*h;  hi = best + 2*h;
    end
end

function [S, tt] = tiePoints(Dl, dl, kp, ss)
    if all(Dl(:) == 0)
        w = [dl(2); -dl(1)];  s0 = -kp*dl/(dl.'*dl);
        S = s0 + w*ss;  tt = ss;  return
    end
    S = zeros(2,0);  tt = [];
    for s1 = ss
        a = Dl(2,2);  b = 2*Dl(1,2)*s1 + 2*dl(2);  c = Dl(1,1)*s1^2 + 2*dl(1)*s1 + 2*kp;
        if abs(a) < 1e-14
            if abs(b) > 1e-14, S(:,end+1) = [s1; -c/b];  tt(end+1) = s1; end %#ok<AGROW>
            continue
        end
        disc = b^2 - 4*a*c;
        if disc < 0, continue, end
        r = sqrt(disc);
        S(:,end+1) = [s1; (-b + r)/(2*a)];  tt(end+1) = s1; %#ok<AGROW>
        if r > 0, S(:,end+1) = [s1; (-b - r)/(2*a)];  tt(end+1) = s1; end %#ok<AGROW>
    end
    keep = all(abs(S) < 1e4, 1);
    S = S(:, keep);  tt = tt(keep);
end

function tf = inHForm(face, x)
% MEMBERSHIP IS DECIDED FROM THE GEOMETRY, not from the sign conditions. The face's boundary
% includes a SEGMENT -- a bounded piece of a line -- which no condition on a whole curve expresses,
% and two constructions that tried were measured wrong. `bridgeMember` asks the defining question
% instead: is there a sweeping segment through x whose base tie point lies in the pair. The name is
% kept because every test above reads the same thing: "is x in this face".
    tf = bridgeMember(face.pair, x(:).');
end

function tf = inCells(face, x)
% The OUTER BOUND, which is what `cells` is: an intersection of sign conditions that CONTAINS the
% face. Used only by the test that pins that relationship.
    for b = 1:numel(face.cells)
        rows = face.cells{b};
        ok = true;
        for r = 1:size(rows, 1)
            c = rows(r,1:6);
            if size(rows,2) >= 8 && rows(r,8) ~= 0
                v = ratQ.pairLineEval(c, rows(r,8), x(:).');
            else
                v = c(1)*x(1)^2 + c(2)*x(1)*x(2) + c(3)*x(2)^2 + c(4)*x(1) + c(5)*x(2) + c(6);
            end
            if rows(r,7)*v < 0, ok = false; break, end
        end
        if ok, tf = true; return, end
    end
    tf = false;
end

function tf = maximiserIsOnTheTie(fx, x)
% The independent classifier: maximise <s,x> - g_k(s) for each piece SEPARATELY and ask whether
% either unconstrained maximiser is feasible for its own piece. If neither is, the maximiser of the
% max is on the tie. No conic, no polynomial, no multiplier.
    s1 = fx.A1 \ (x - fx.b1);
    if gval(fx.A1, fx.b1, fx.k1, s1) >= gval(fx.A2, fx.b2, fx.k2, s1) - 1e-12
        tf = false;  return
    end
    s2 = fx.A2 \ (x - fx.b2);
    if gval(fx.A2, fx.b2, fx.k2, s2) >= gval(fx.A1, fx.b1, fx.k1, s2) - 1e-12
        tf = false;  return
    end
    tf = true;
end

function v = gval(A, b, k, s), v = (s.'*A*s)/2 + b.'*s + k; end

function r = relEval(P, x, z)
% |p| at (x,z), RELATIVE to the sum of the term magnitudes: a residual means nothing beside
% coefficients of size 1e6 unless it is read against the terms that produced it.
    sz = size(P);  sz(end+1:3) = 1;
    v = 0;  s = 0;
    for i = 1:sz(1)
        for j = 1:sz(2)
            for k = 1:sz(3)
                c = P(i,j,k);
                if c == 0, continue, end
                t = c * x(1)^(i-1) * x(2)^(j-1) * z^(k-1);
                v = v + t;  s = s + abs(t);
            end
        end
    end
    if s == 0, s = 1; end
    r = abs(v)/s;
end
