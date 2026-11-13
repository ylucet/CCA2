function [M, a, c0] = xyFrame(Q, L, c)
% xyFrame  The affine change of variables that turns an INDEFINITE quadratic into exactly x*y.
%
% [input]  Q, L, c : q(x) = 1/2 x'Qx + L'x + c with Q indefinite (det Q < 0 in 2D).
% [output] M       : 2x2, with q(M z) = z1*z2 + a'z + c0. Equivalently M'QM = [0 1; 1 0].
%          a, c0   : the leftover affine part, a = M'L and c0 = c.
%
% WHY. cPLQ's Step 1 closed forms (plq_1p.convexEnvelope1's nCE == 0/1/2 branches) and its
% "convex edge" test (an edge of positive finite SLOPE) are written for f = x*y and for nothing
% else -- they never reference the piece's own function. That is exact for cPLQ, every caller of
% which passes x*y, and silently wrong in CCA2, where quaPolToPlq builds a general quadratic per
% face; measured, Case C returned f*(0.3,0.4) = 0.4 where the truth is 0.125. Rather than
% re-derive those closed forms for a general q -- which is the whole [COAP] appendix again --
% move the PROBLEM into the frame the formulas are already correct in. Convex and concave q are
% handled directly elsewhere, so indefinite is the only case that needs this.
%
% THE CONSTRUCTION, in two congruences.
%   1. Q = P D P' with D = diag(l1,l2), l1 > 0 > l2 (indefinite, so the signs are opposite).
%      S = P * diag(1/sqrt(l1), 1/sqrt(-l2))  gives  S'QS = diag(1,-1).
%   2. diag(1,-1) is u^2 - v^2, and u = w1 + w2, v = w1 - w2 turns that into 4*w1*w2. So with
%      T = [1 1; 1 -1]/sqrt(2),  T' diag(1,-1) T = [0 1; 1 0],  whose quadratic form
%      1/2 z'[0 1;1 0]z is exactly z1*z2.
%   M = S*T then satisfies M'QM = [0 1; 1 0], i.e. 1/2 (Mz)'Q(Mz) = z1*z2.
% The linear part comes along unchanged: L'(Mz) = (M'L)'z, and the constant is untouched.
%
% HOW THE CONJUGATE COMES BACK. With g(z) = f(Mz),
%       f*(s) = sup_x <s,x> - f(x) = sup_z <s,Mz> - g(z) = sup_z <M's,z> - g(z) = g*(M's),
% so a conjugate computed in the z-frame is read at M's. Stripping the affine part as well,
% h(z) = z1*z2 + I_D'(z) gives g = h + <a,.> + c0 and hence
%       f*(s) = h*(M's - a) - c0,
% which is the substitution plq_1p.conjugate applies. Note M is NOT orthogonal in general (step 1
% rescales by the eigenvalues), so this is a shear-and-scale, not merely a rotation; that is
% required -- a pure rotation only reaches l1*u^2 + l2*v^2, not u*v.

    evd = eig(double((double(Q) + double(Q)')/2));
    if ~(min(evd) < 0 && max(evd) > 0)
        error('xyFrame:notIndefinite', ...
            ['xyFrame needs an INDEFINITE Q (one positive and one negative eigenvalue); got ' ...
             'eigenvalues [%g %g]. Convex and concave faces are handled without a frame change.'], ...
            evd(1), evd(2));
    end

    % Q ALREADY x*y: only the affine part differs, so no change of variables is needed at all --
    % co(x*y + l) = co(x*y) + l, and stripping l is exactly what the caller does with (a,c0).
    % Worth taking explicitly: the general construction below would return a perfectly valid but
    % IRRATIONAL M for this Q (any M with M'JM = J), and irrational vertices cost the symbolic
    % pipeline dearly for no gain.
    if norm(double(Q) - [0 1; 1 0]) <= 1e-12
        M  = sym(eye(2));
        a  = sym(double(L(:)), 'r');
        c0 = sym(double(c), 'r');
        return
    end

    % EXACT, not floating point. M's entries carry 1/sqrt(lambda) factors, which are irrational
    % whenever the eigenvalues are -- that is inherent (reducing an indefinite form to the
    % hyperbolic form u*v needs a square root; it is rational only when the ratio of the
    % eigenvalues happens to be a rational square). The choice that matters is how the
    % irrational number is CARRIED. As a double it becomes, under sym, an exact binary rational
    % with a 2^52 denominator, and every downstream solve/simplify in cPLQ then works with
    % 17-digit numerators -- measured to take the pipeline from seconds to over 20 minutes on
    % q = x^2+3xy-2y^2+x. Built symbolically it stays a clean surd like 3^(1/2)/2, which is both
    % faster and consistent with this codebase's exact-arithmetic target (the nCE==2 envelope
    % already carries sqrt terms of exactly this kind).
    Qs = sym((double(Q) + double(Q)')/2, 'r');   % 'r': recover the intended rationals, not 2^-52
    Ls = sym(double(L(:)), 'r');
    cs = sym(double(c), 'r');

    [P, D] = eig(Qs);
    d = simplify(diag(D));
    if simplify(d(1)) > 0, ip = 1; im = 2; else, ip = 2; im = 1; end

    % Normalize each eigenvector so that S'*Q*S = diag(1,-1): for a symmetric Q the columns are
    % already orthogonal, and p'Qp = lambda*(p'p), so dividing by sqrt(|lambda|*(p'p)) is exactly
    % what puts +-1 on the diagonal.
    S = sym(zeros(2,2));
    for k = [1 2]
        if k == 1, j = ip; else, j = im; end
        p = P(:,j);
        S(:,k) = p / sqrt(abs(d(j)) * (p.' * p));
    end
    T = sym([1 1; 1 -1]) / sqrt(sym(2));
    M = simplify(S * T);

    a  = simplify(M.' * Ls);
    c0 = cs;
end
