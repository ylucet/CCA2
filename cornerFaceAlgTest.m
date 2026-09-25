classdef cornerFaceAlgTest < matlab.unittest.TestCase
% [fast] cornerFaceAlg, checked against the DEFINITION of a corner face rather than stored output:
% at the corner s -- a point of the tie on one of the cell's facets -- the face's value is the
% AFFINE function <s,x> - g(s), so the implicit polynomial must vanish on it and the chosen branch
% must equal it.
%
% The inputs are the corner calls biconjQ makes on the poly2 kink at c = -2 (captured 2026-09-25).
% Five of them used to raise ratQ:overflow: the square in the implicit polynomial passed 2^53 while
% the polynomial, content removed, does not.

    methods (Test)
        function theCORNERSOfThePOLY2KinkBuildAndCarryTheirValue(t)
            calls = {{[0 0 0 0 1 0 0 0 0 0],     2, [1 0 4;-1 0 4;0 -1 0;3 2 16;-3 -2 10;2 -3 2], [8 -24 -8 16 -128 -96],  1}, ...
                     {[0 0 0 0 9 6 4 -4 32 24], 26, [1 0 4;-1 0 4;0 -1 0;3 2 16;-3 -2 10;2 -3 2], [8 -24 -8 16 -128 -96], -1}, ...
                     {[0 0 0 0 1 0 1 0 0 0],     2, [0 1 0;-2 -1 8;2 -3 8;1 0 -2;3 2 -10],         [1 0 1 -4 -8 4],          1}, ...
                     {[0 0 0 0 0 0 0 1 2 -1],    1, [0 1 0;-2 -1 8;2 -3 8;1 0 -2;3 2 -10],         [1 0 1 -4 -8 4],         -1}, ...
                     {[0 0 0 0 1 0 1 0 0 0],     2, [0 1 0;-2 -1 8;2 -3 8;3 2 16;-3 -2 10;2 -3 2], [8 -24 18 16 -128 -96],  1}, ...
                     {[0 0 0 0 9 6 4 -4 32 24], 26, [0 1 0;-2 -1 8;2 -3 8;3 2 16;-3 -2 10;2 -3 2], [8 -24 18 16 -128 -96], -1}};
            nFaces = 0;
            for k = 1:numel(calls)
                a = calls{k};
                faces = cornerFaceAlg(a{:});
                n = a{1};  d = a{2};
                N = [n(5) n(6); n(6) n(7)];  m = [n(8); n(9)];  q = n(10);
                g = @(s) ((s.'*N*s)/2 + m.'*s + q) / d;
                for f = faces
                    nFaces = nFaces + 1;
                    s = f.s;
                    for off = [0 0; 0.1 0; 0 0.1; -0.07 0.05].'
                        x = f.p + off;                               % around grad g at the corner
                        z = s.'*x - g(s);
                        t.verifyLessThan(relResidual(f.FP, x, z), 1e-10, ...
                            sprintf('call %d: the polynomial vanishes on the corner''s value', k));
                    end
                    % THE BRANCH only at the corner's own point: the index names this corner inside
                    % its region, and an offset can leave it -- for an affine g the region at p is a
                    % cone, and beside it the top root is the CONJUGATE corner's value. The
                    % polynomial vanishes on both, which is why the residual above is asked everywhere.
                    t.verifyEqual(algBranch(f.FP, f.p.', f.FB), s.'*f.p - g(s), 'AbsTol', 1e-7, ...
                        sprintf('call %d: the chosen branch IS the corner''s value at grad g(s)', k));
                end
            end
            t.verifyGreaterThan(nFaces, 0, 'the fixture must produce corners, or this proves nothing');
        end
    end
end

function r = relResidual(P, x, z)
% |P(x,y,z)| relative to the sum of its terms' magnitudes.
    v = 0;  a = 0;
    for i = 1:size(P,1)
        for j = 1:size(P,2)
            for k = 1:size(P,3)
                term = P(i,j,k) * x(1)^(i-1) * x(2)^(j-1) * z^(k-1);
                v = v + term;  a = a + abs(term);
            end
        end
    end
    r = abs(v) / max(a, realmin);
end
