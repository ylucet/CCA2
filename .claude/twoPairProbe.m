function twoPairProbe()
% Can a cell whose tight subsystem mixes TWO conjugate pairs actually be built? `cellInterior` skips
% such a certificate and keeps the cell; the case has never occurred in the pipeline, so before
% building the compositum arithmetic the question is whether it is reachable at all -- and what the
% cost of keeping it is.
%
% Two pairs through a common point: x^2 - 2y^2 (slopes +-1/sqrt2) and x^2 - 3y^2 (+-1/sqrt3). Three
% concurrent lines drawn from BOTH pairs pinch a cell exactly as one pair's three did, but now the
% certificate needs sqrt2 and sqrt3 together.
    P1 = ratQ.conic([1 0 -2 0 0 0]);
    P2 = ratQ.conic([1 0 -3 0 0 0]);
    C = [P1, 1; P1, 2; P2, 1; P2, 2];
    n = 0;  kept = 0;  unsure = 0;
    for s1 = [1 -1]
        for s2 = [1 -1]
            for s3 = [1 -1]
                for s4 = [1 -1]
                    rows = [1 s1; 2 s2; 3 s3; 4 s4];
                    [tf, info] = cellInterior(C, rows);
                    n = n + 1;
                    kept = kept + tf;
                    unsure = unsure + strcmp(info.path, 'unsureTie');
                    % the truth, by sampling: four lines through the origin cut six... eight sectors
                    th = linspace(0, 2*pi, 20001).';
                    X = [cos(th), sin(th)];
                    v = [s1*ratQ.pairLineEval(P1,1,X), s2*ratQ.pairLineEval(P1,2,X), ...
                         s3*ratQ.pairLineEval(P2,1,X), s4*ratQ.pairLineEval(P2,2,X)];
                    truth = any(all(v > 1e-12, 2));
                    if tf ~= truth
                        fprintf('  %+d %+d %+d %+d -> live %d (%s) but sampled %d\n', ...
                                s1, s2, s3, s4, tf, info.path, truth);
                    end
                end
            end
        end
    end
    fprintf('%d patterns | kept %d | reached the unsure path %d\n', n, kept, unsure);
end
