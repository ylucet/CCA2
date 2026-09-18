classdef cornerWedgeRowsTest < matlab.unittest.TestCase
% [fast] Does a corner's ROW SPELLING describe the same wedge as the predicate it was built from?
% The predicate `mu = M \ (p - x) >= 0` is the oracle: it is what `cornerFaceAlg` derived the face
% from and what `envelopeEval` answers with, so the rows are correct exactly when they agree with it
% everywhere. Nothing here is a stored value.

    methods (Test)
        function theROWSAreTheSameWEDGEAsThePredicate(t)
            [faces, tie] = fixture();
            t.verifyNotEmpty(faces, 'the fixture must produce corners at all');
            rng(11);  X = [rand(300,1)*4-2, rand(300,1)*4-2];
            total = 0;
            for k = 1:numel(faces)
                f = faces(k);
                rows = cornerWedgeRows(f, f.src.n, f.src.d, tie, f.src.side);
                t.verifyEqual(size(rows), [2 8], 'a wedge has exactly two edges');
                for i = 1:size(X,1)
                    mu = f.M \ (f.p - X(i,:).');
                    byRows = true;
                    for r = 1:2
                        if rows(r,7)*rowValue(rows(r,:), X(i,:)) < 0, byRows = false; break, end
                    end
                    t.verifyEqual(byRows, all(mu >= 0), ...
                                  sprintf('face %d, point %d', k, i));
                    total = total + 1;
                end
            end
            t.verifyGreaterThan(total, 0);
        end

        function eachEDGEPassesThroughTheCornerItBounds(t)
        % The sharper statement, and the one that caught the defect this file was written for: edge i
        % is the line { mu_i = 0 }, which CONTAINS the corner's gradient p. A row whose line misses p
        % is the wrong line -- the conjugate one, or an arithmetic slip -- however plausible its
        % signs look on a sample.
            [faces, tie] = fixture();
            for k = 1:numel(faces)
                f = faces(k);
                rows = cornerWedgeRows(f, f.src.n, f.src.d, tie, f.src.side);
                for i = 1:2
                    dir = f.M(:, 3-i);                   % along the edge where mu_i = 0
                    P = [f.p.' + dir.'/norm(dir); f.p.' - dir.'/norm(dir)];
                    for j = 1:2
                        v = rowValue(rows(i,:), P(j,:));
                        t.verifyLessThan(abs(v), 1e-9 * max(1, max(abs(rows(i,1:6)))), ...
                                         sprintf('face %d edge %d point %d', k, i, j));
                    end
                end
            end
        end

        function aNEGATIVELeadingCoefficientIsNotASIGNFLIP(t)
        % The defect this file exists for: every product in Q(sigma) carries the corner quadratic's
        % leading coefficient in its DENOMINATOR, so `al < 0` made every element negative-denominated
        % -- and a normalisation that flipped the denominator twice negated them all. The two sides of
        % a tie have opposite `al`, so asking for both is what exercises it.
        %
        % `al` is the TIE restricted to the facet direction, so it is the tie's own orientation that
        % sets its sign and not `side` -- the same curve written with the opposite sign is the
        % fixture that exercises the other branch.
            [fp, tp] = fixture(+1, +1);
            [fm, tm] = fixture(+1, -1);
            t.verifyNotEmpty(fp);  t.verifyNotEmpty(fm);
            t.verifyGreaterThan(fp(1).sigma.quad(1), 0);
            t.verifyLessThan(   fm(1).sigma.quad(1), 0, ...
                         'the second fixture must give al < 0, or this test proves nothing');
            for kk = 1:(numel(fp) + numel(fm))
                if kk <= numel(fp), f = fp(kk);  tie = tp;
                else,               f = fm(kk - numel(fp));  tie = tm;  end
                rows = cornerWedgeRows(f, f.src.n, f.src.d, tie, f.src.side);
                xw = f.p - f.M * [1; 1];                 % mu = [1; 1]: strictly inside the wedge
                for r = 1:2
                    t.verifyGreaterThan(rows(r,7)*rowValue(rows(r,:), xw.'), 0, ...
                                        'the witness must be strictly inside every row');
                end
            end
        end

        function theARGUMENTSAreChecked(t)
            [faces, tie] = fixture();
            t.verifyError(@() cornerWedgeRows(faces(1), 1, 1), 'cornerWedgeRows:nargin');
            t.verifyError(@() cornerWedgeRows(struct('a',1), faces(1).src.n, 1, tie, 1), ...
                          'cornerWedgeRows:shape');
        end
    end
end

% ==================================================================================================
function [faces, tie] = fixture(side, orient)
% g(s) = s1^2 + s2^2 on the square |s_i| <= 1, tied against 2 s1^2 + 2 s2^2 = 3 -- whose meets with
% the facets are at sigma = +-1/sqrt(2), IRRATIONAL, which is the whole point. `orient` writes the
% same curve with the opposite sign, which is what flips the sign of `al`.
    if nargin < 1, side = +1; end
    if nargin < 2, orient = +1; end
    n = [0 0 0 0, 2, 0, 2, 0, 0, 0];
    P = [1 0 1; -1 0 1; 0 1 1; 0 -1 1];
    tie = orient * [2 0 2 0 0 -3];
    faces = cornerFaceAlg(n, 1, P, tie, side);
end

function v = rowValue(row, x)
    if row(8) == 0
        v = row(4)*x(1) + row(5)*x(2) + row(6);
    else
        v = ratQ.pairLineEval(row(1:6), row(8), x);
    end
end
