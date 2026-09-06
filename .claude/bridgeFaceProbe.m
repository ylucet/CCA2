function bridgeFaceProbe()
% WHAT FUNCTION CLASS does a face of `co f` belong to, when two CONVEX quadratic pieces bridge?
%
% biconjQ.m's header asserts that the faces which cannot be written rationally are AFFINE cells
% `<p,x> - f*(p)` with `p` a dual VERTEX of degree up to 4, and that naming `p` (the planned
% AlgAlg type) is therefore enough. This checks that assertion on the simplest possible bridge,
% and it is the whole reason the script exists: if the bridge face is not affine, AlgAlg as
% specified cannot hold it.
%
% THE INSTANCE, chosen so every quantity is available in closed form and no domain clipping can
% muddy the question -- both pieces on the whole plane, so `min` alone defines f:
%
%     q1(v) = |v|^2 / 2                    (Hessian I)
%     q2(v) = |v - a|^2 + c                (Hessian 2I -- DIFFERENT, which is the point)
%
% THE DERIVATION being tested. q1*(s) = |s|^2/2 and q2*(s) = <s,a> + |s|^2/4 - c, so the two
% pieces tie where |s|^2/4 - <s,a> + c = 0 -- a CIRCLE of centre 2a and radius R = 2*sqrt(|a|^2-c),
% not a line. On that tie, co f(x) = <s,x> - |s|^2/2 = <s, x - 2a> + 2c, and the sup of a LINEAR
% function over a circle is attained at the far point, giving
%
%     co f(x) = <2a, x - 2a> + R*|x - 2a| + 2c .
%
% A NORM, i.e. the square root of a quadratic. Not affine, not quadratic, not rational.
%
% THE ORACLE IS INDEPENDENT of that derivation: the primal definition
%   co f(x) = inf over t in [0,1] and (1-t)u + t v = x of (1-t)q1(u) + t q2(v),
% whose inner inf-convolution is available in closed form for quadratics, leaving a ONE-dimensional
% minimisation that is done here by brute force on a fine grid of t. It shares no algebra with the
% dual derivation above -- it never forms a conjugate at all.
%
%   inf-conv of |z|^2/(2*lam1) and |z - lam2*a|^2/lam2  =  |x - lam2*a|^2 / (2*(lam1 + lam2/2)),
%   so   co f(x) = min over t in [0,1] of  |x - t*a|^2/(2 - t) + t*c .

    a = [1.3; -0.7];
    c = -0.5;
    R = 2*sqrt(a.'*a - c);

    rng(20260905);
    X = 6*(rand(4000,2) - 0.5);

    ts = linspace(0, 1, 200001).';          % the oracle's only parameter
    nBridge = 0;  worst = 0;  worstAff = 0;
    got = zeros(size(X,1),1);  want = zeros(size(X,1),1);  tstar = zeros(size(X,1),1);
    for i = 1:size(X,1)
        x = X(i,:).';
        vals = sum((x.' - ts.*a.').^2, 2)./(2 - ts) + ts*c;
        [got(i), j] = min(vals);
        tstar(i) = ts(j);
        want(i) = 2*(a.'*(x - 2*a)) + R*norm(x - 2*a) + 2*c;
    end

    % The closed form is the BRIDGE's formula, so it is only claimed where the bridge is active --
    % t strictly inside (0,1). Outside, one piece is exposed on its own and co f = that piece.
    inB = tstar > 1e-6 & tstar < 1 - 1e-6;
    nBridge = sum(inB);
    worst = max(abs(got(inB) - want(inB)));

    fprintf('a = (%g, %g), c = %g, R = %g\n', a(1), a(2), c, R);
    fprintf('points %d, of which BRIDGE (0 < t < 1): %d\n', size(X,1), nBridge);
    fprintf('max |oracle - closed form| on the bridge: %.3e\n', worst);
    fprintf('   (the oracle grid step in t is %.1e, so ~1e-5 is the grid floor)\n', ts(2)-ts(1));

    % IS IT AFFINE? Fit the best affine function on the bridge points and report the residual.
    % An affine face -- what AlgAlg is built to hold -- would fit to the grid floor.
    Ab = [X(inB,:), ones(nBridge,1)];
    resAff = Ab*(Ab\got(inB)) - got(inB);
    worstAff = max(abs(resAff));

    % IS IT QUADRATIC? Same question for the full quadratic basis, which is what a QuaCon face is.
    Xb = X(inB,:);
    Aq = [Xb(:,1).^2, Xb(:,1).*Xb(:,2), Xb(:,2).^2, Xb, ones(nBridge,1)];
    resQ = Aq*(Aq\got(inB)) - got(inB);

    fprintf('\nbest AFFINE fit to the bridge face : max residual %.3e\n', worstAff);
    fprintf('best QUADRATIC fit to the bridge face: max residual %.3e\n', max(abs(resQ)));
    fprintf(['\nA face that neither an affine nor a quadratic form can represent. The closed form\n' ...
             'above says why: it carries |x - 2a|, the square root of a quadratic.\n']);

    % THE SAME-HESSIAN CONTROL. With q2 = |v-a|^2/2 + c the two conjugates differ by an AFFINE
    % function, so the tie is a LINE rather than a circle, and the sup of a concave quadratic over
    % a line IS quadratic. This is the case that stays rational, and it is worth pinning next to
    % the one that does not -- the difference is the Hessian and nothing else.
    fprintf('\n--- control: the SAME Hessian on both pieces ---\n');
    got2 = zeros(size(X,1),1);  t2 = zeros(size(X,1),1);
    for i = 1:size(X,1)
        x = X(i,:).';
        % perspective inf-conv with both widths lam: |x - t*a|^2/2 + t*c
        vals = sum((x.' - ts.*a.').^2, 2)/2 + ts*c;
        [got2(i), j] = min(vals);
        t2(i) = ts(j);
    end
    inB2 = t2 > 1e-6 & t2 < 1 - 1e-6;
    if any(inB2)
        Xb2 = X(inB2,:);
        Aq2 = [Xb2(:,1).^2, Xb2(:,1).*Xb2(:,2), Xb2(:,2).^2, Xb2, ones(sum(inB2),1)];
        r2 = Aq2*(Aq2\got2(inB2)) - got2(inB2);
        fprintf('bridge points %d, best QUADRATIC fit max residual %.3e\n', sum(inB2), max(abs(r2)));
    else
        fprintf('no interior-t bridge points in the control\n');
    end
end
