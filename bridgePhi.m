%BRIDGEPHI  The tie residual along the sweeping family, as an integer polynomial in the multiplier.
%
% objective: the one polynomial both halves of the bridge rest on, so it is defined once.
%   `bridgeWitness` asks whether it has a root in [-1,0] -- that is membership in the swept region --
%   and `bridgeFaceAlg` reads the stationary multipliers off the same roots to get the value. A
%   second spelling of it would be a second definition of what the face IS.
%
% [input]  N1,m1,N2,m2 : the two faces over the common denominator d
%          Dl, dl, kp  : the tie's own data, over the same d
%          xn, xd      : the point, as xn/xd
% [output] Phi : the quartic, ASCENDING in mu, integer -- divided by its POSITIVE content when the raw
%                coefficients pass 2^53 (same roots, same signs; 2026-10-11)
%          Dt  : det(M_mu), ascending, integer
%          U1, U2 : the components of U = xd det(M_mu) s(mu), ascending -- what lets a caller ask
%                   WHERE the base point is, since a's + c has the sign of (a'U + c xd det M) times
%                   the sign of det M, both polynomials in mu
%
% [input]  exact : (optional, default true) false lets xn be a DOUBLE and skips the integrality
%                  check, so the same formula serves `bridgeMember`'s cache. One definition of the
%                  polynomial, two precisions -- a second spelling of it would be a second
%                  definition of what the face IS.
%
% At the multiplier mu the base point is s = adj(M_mu)(d x - c_mu)/det(M_mu), and the tie equation
% s'Dl s + 2 dl's + 2 kp = 0 clears to
%
%       Phi(mu) = U' Dl U + 2 xd det(M_mu) dl'U + 2 kp xd^2 det(M_mu)^2 ,
%
% with no division performed anywhere.

function [Phi, Dt, U1, U2] = bridgePhi(N1, m1, N2, m2, Dl, dl, kp, d, xn, xd, exact)
% Phi as above, with det(M_mu) and the two components of U = xd det(M_mu) s(mu) alongside it -- all
% ASCENDING integer vectors in mu. U is what lets a CALLER ask where the base point is: the row
% a's + c has the sign of (a'U + c xd det M) times the sign of det M, both polynomials in mu, so
% "is the base point inside the pair" is a sign at a root and `ratQ.signAt` decides it.
    if nargin < 11 || exact
        % DOUBLES FIRST, BIG ONLY WHERE THEY RUN OUT (2026-10-11). Phi carries d, xd and the pair's
        % data to the fourth power, so its raw coefficients pass 2^53 on erratum Example 2 (69-76
        % bits) while the reduced quartic needs 18-22. ratQ.exactFormula re-runs this same formula
        % in BigInteger (bigZ) and divides out the positive content -- no root and no sign changes,
        % and those are all an exact caller reads. U and det(M_mu) are used unreduced: they must fit.
        [~, Dt, U1, U2] = phiRaw(N1, m1, N2, m2, Dl, dl, kp, d, xn, xd);
        ratQ.chk([U1, U2, Dt], 'U and det(M_mu)');
        Phi = ratQ.exactFormula(@phiRaw, 'Phi', N1, m1, N2, m2, Dl, dl, kp, d, xn, xd);
    else
        [Phi, Dt, U1, U2] = phiRaw(N1, m1, N2, m2, Dl, dl, kp, d, xn, xd);
    end
end

function [Phi, Dt, U1, U2] = phiRaw(N1, m1, N2, m2, Dl, dl, kp, d, xn, xd)
% the formula itself; integer inputs as doubles or as bigZ
    e11 = [N1(1,1), N1(1,1)-N2(1,1)];  e12 = [N1(1,2), N1(1,2)-N2(1,2)];
    e22 = [N1(2,2), N1(2,2)-N2(2,2)];
    Dt  = vadd(conv(e11, e22), -conv(e12, e12));
    a11 = e22;  a12 = -e12;  a22 = e11;

    w1 = [d*xn(1) - xd*m1(1), -xd*(m1(1)-m2(1))];     % xd (d x - c_mu), component by component
    w2 = [d*xn(2) - xd*m1(2), -xd*(m1(2)-m2(2))];
    U1 = vadd(conv(a11, w1), conv(a12, w2));
    U2 = vadd(conv(a12, w1), conv(a22, w2));

    Phi = vadd(vadd(Dl(1,1)*conv(U1,U1), 2*Dl(1,2)*conv(U1,U2)), Dl(2,2)*conv(U2,U2));
    Phi = vadd(Phi, 2*xd*conv(Dt, vadd(dl(1)*U1, dl(2)*U2)));
    Phi = vadd(Phi, 2*kp*xd^2*conv(Dt, Dt));
end

% --------------------------------------------------------------------------------------------------
function c = vadd(a, b)
% Add two ASCENDING integer vectors of possibly different degree.
    n = max(numel(a), numel(b));
    c = [a, zeros(1, n-numel(a))] + [b, zeros(1, n-numel(b))];
end
