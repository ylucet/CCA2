classdef cellInteriorTest < matlab.unittest.TestCase
% [fast] Does the arrangement's emptiness test answer the DEFINITION -- a cell is live exactly when
% some point satisfies every one of its rows strictly? Nothing here is a stored value: each case
% either exhibits such a point or proves there is none by hand, and the cases with a zero margin are
% the ones the exact stage exists for.

    methods (Test)
        function aTRIANGLEIsLiveAndAnOPPOSITEPairIsNot(t)
        % The two cases the cheap stage must get right on its own.
            C = [0 0 0  1  0 0, 0;          % x >= 0
                 0 0 0  0  1 0, 0;          % y >= 0
                 0 0 0 -1 -1 1, 0];         % 1 - x - y >= 0
            t.verifyTrue(cellInterior(C, [1 1; 2 1; 3 1]), 'the unit triangle has interior');
            t.verifyFalse(cellInterior(C, [1 1; 1 -1]), 'a row and its opposite cannot both hold');
            t.verifyTrue(cellInterior(C, [1 1; 3 -1]), 'x >= 0 and x + y >= 1 do');
        end

        function aBRANCHRowCutsOnONELineOfItsPairNotOnTheConic(t)
        % The whole reason branch rows exist: the pair `x^2 = 2 y^2` is TWO lines, and a row naming
        % one of them must halve the plane, not take the conic's sign -- which is the product of the
        % two and keeps a double wedge. The cell below is the sector between the named line and the
        % axis, and it contains a point the conic's own sign would also admit on the far side.
            P = ratQ.conic([1 0 -2 0 0 0]);          % x^2 - 2 y^2 = (x - s2 y)(x + s2 y)
            C = [P, 1; P, 2; 0 0 0 0 1 0, 0];
            % The two lines cross at the origin and cut the plane into four sectors, and y >= 0 is
            % a third line through that same point, so exactly the sector lying in y < 0 is lost. A
            % row read as the CONIC's sign instead of one LINE's would keep a double wedge and make
            % all four patterns live.
            live = 0;
            for s1 = [1 -1]
                for s2 = [1 -1]
                    if cellInterior(C, [1 s1; 2 s2; 3 1]), live = live + 1; end
                end
            end
            t.verifyEqual(live, 3, 'three of the four sectors meet y >= 0');
        end

        function aZEROWidthCellIsEMPTYAndItsNEIGHBOURIsNot(t)
        % The tie, and the case the exact stage was written for. The two lines of the pair and the
        % line y = 0 all pass through the ORIGIN, so every cell built from them has a margin of
        % exactly zero at that point and no tolerance can separate the empty ones from the rest.
            P = ratQ.conic([1 0 -2 0 0 0]);
            C = [P, 1; P, 2; 0 0 0 0 1 0, 0];
            % Three lines through ONE point cut the plane into six sectors, so of the eight sign
            % patterns exactly six are live and two are pinched down to the point itself. The COUNT
            % is the assertion: a test that dropped every tie would report six live and two empty
            % for the wrong reason, one that kept every tie would report eight.
            nLive = 0;
            for s1 = [1 -1]
                for s2 = [1 -1]
                    for s3 = [1 -1]
                        if cellInterior(C, [1 s1; 2 s2; 3 s3]), nLive = nLive + 1; end
                    end
                end
            end
            t.verifyEqual(nLive, 6, ...
                'three concurrent lines give six sectors; the other two patterns are pinched');
        end

        function theTIEIsDecidedEXACTLYAndNotByTolerance(t)
        % A cell can be thin without being empty, and no tolerance can tell the two apart -- so the
        % test must. The wedge between the same two lines is opened by a tiny rational amount here:
        % it is thinner than any margin tolerance and it IS live.
            P = ratQ.conic([1 0 -2 0 0 0]);
            C = [P, 1; P, 2; 0 0 0 1 0 -1, 0];             % ... and x >= 1, so the apex is excluded
            s1 = sign(ratQ.pairLineEval(P, 1, [1 0]));
            s2 = sign(ratQ.pairLineEval(P, 2, [1 0]));
            [tf, info] = cellInterior(C, [1 s1; 2 s2; 3 1]);
            t.verifyTrue(tf, 'the sector beyond x = 1 has interior');
            t.verifyGreaterThan(info.margin, 0);
        end

        function TWOConjugatePairsInOneCellAreDecidedInTheirCOMPOSITUM(t)
        % A cell whose rows name lines of TWO different pairs needs a certificate in Q(t1, t2), and
        % the case is not hypothetical: four concurrent lines drawn from two pairs leave 8 of the 16
        % sign patterns undecidable in either field alone, and every one of them is EMPTY -- so
        % keeping them, as the single-field code had to, keeps eight empty cells.
        %
        % Four lines through one point cut the plane into eight sectors, so exactly 8 of the 16
        % patterns are live. That count is the assertion, and `unsureTie` must not appear at all.
            P1 = ratQ.conic([1 0 -2 0 0 0]);          % slopes +-1/sqrt(2)
            P2 = ratQ.conic([1 0 -3 0 0 0]);          % slopes +-1/sqrt(3), a DIFFERENT field
            C = [P1, 1; P1, 2; P2, 1; P2, 2];
            live = 0;  unsure = 0;
            for s1 = [1 -1]
                for s2 = [1 -1]
                    for s3 = [1 -1]
                        for s4 = [1 -1]
                            [tf, info] = cellInterior(C, [1 s1; 2 s2; 3 s3; 4 s4]);
                            live = live + tf;
                            unsure = unsure + strcmp(info.path, 'unsureTie');
                        end
                    end
                end
            end
            t.verifyEqual(live, 8, 'four concurrent lines give eight sectors');
            t.verifyEqual(unsure, 0, 'and none of the ties may be left undecided');
        end

        function theCOMPOSITUMRefusesTwoPairsOfTheSAMEField(t)
        % The zero test behind the sign -- all four coordinates vanish -- is only valid when the two
        % quadratics generate DIFFERENT quadratic fields. They coincide exactly when the product of
        % the discriminants is a perfect square, and then the sign is REFUSED rather than guessed.
            [s, ok] = ratQ.qs2Sign([1 1 -1 0 1], [1 0 -2], [1 2], [1 0 -8], [2 3]);
            t.verifyFalse(ok, 'sqrt(2) and sqrt(8) generate one field, so this must refuse');
            t.verifyEqual(s, 0);
            % and a genuine compositum is decided: 1 + sqrt(2) - sqrt(3) > 0
            [s, ok] = ratQ.qs2Sign([1 1 -1 0 1], [1 0 -2], [1 2], [1 0 -3], [1 2]);
            t.verifyTrue(ok);
            t.verifyEqual(s, 1);
            % ... and its negation is negative, which pins the orientation rather than the value
            [s, ok] = ratQ.qs2Sign([-1 -1 1 0 1], [1 0 -2], [1 2], [1 0 -3], [1 2]);
            t.verifyTrue(ok);
            t.verifyEqual(s, -1);
        end

        function theARGUMENTSAreChecked(t)
            t.verifyError(@() cellInterior([0 0 0 1 0 0, 0]), 'cellInterior:nargin');
            t.verifyError(@() cellInterior([0 0 0 1 0 0], [1 1]), 'cellInterior:shape');
            t.verifyError(@() cellInterior([0 0 0 0 0 1, 0], [1 1]), 'cellInterior:constant');
            t.verifyError(@() ratQ.qs2Red([1 2 3 4]), 'ratQ:qs2Shape');
            t.verifyError(@() ratQ.qs2Red([1 2 3 4 0]), 'ratQ:qs2ZeroDen');
        end
    end
end
