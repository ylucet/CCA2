function bridgeEliminationSpike()
% SPIKE for the implicit algebraic face: can the bridge's data be produced EXACTLY, over the
% rationals, from the two pieces' own integer coefficients? If not the type is moot -- there would
% be nothing exact to store in it. Nothing is built here; this answers one question on one instance
% so the next session starts from a working kernel rather than from the derivation again.
%
% THE ELIMINATION, for two pieces strictly convex on the whole plane (no clipping, which is the
% case whose answer is already known in closed form and can therefore be checked):
%
%   q_k(v) = (v - a_k)'H_k(v - a_k)/2 + c_k,      q_k*(s) = <s,a_k> + s'inv(H_k)s/2 - c_k.
%
% The pieces TIE where q_1*(s) = q_2*(s) -- a conic C(s) = 0. At a tie the primal point lies on the
% segment joining the two tangency points y_k(s) = a_k + inv(H_k)s, so
%
%   D(s; x) = det[ x - y_1(s) , y_2(s) - y_1(s) ] = 0
%
% is a second conic carrying x as a parameter, and the face value is z = <s,x> - q_1*(s).
%
% WHAT IS ESTABLISHED HERE
%   1. C has EXACT INTEGER coefficients from integer inputs -- clearing inv(H_k) costs a factor
%      det(H_k), an integer, so nothing leaves Q. Checked against the derived circle.
%   2. D likewise, and the pair (C,D) does select the right s: solving them and reading off z
%      reproduces the closed form.
%   3. The BRANCH SELECTOR is t in [0,1] -- x must lie BETWEEN the two tangency points. Two conics
%      meet in up to four points and the others are not the bridge.
%
% WHAT IS NOT HERE, and is the next thing to build: the RESULTANT eliminating s symbolically to give
% p(x,y,z) = 0. ratQ has univariate machinery (polyMul, polyPseudoRem, polyGcd, Sturm) but no
% bivariate resultant.

    a2n = [13; -7];   a2d = 10;   c2n = -1;   c2d = 2;
    a = a2n/a2d;      c = c2n/c2d;

    fprintf('q1(v) = |v|^2/2 ,  q2(v) = |v-a|^2 + c   with a = (%d,%d)/%d, c = %d/%d\n\n', ...
            a2n(1), a2n(2), a2d, c2n, c2d);

    % ---- the tie conic C(s), cleared to integers -------------------------------------------------
    % C(s) = q1*(s) - q2*(s) = |s|^2/4 - <s,a> + c;  multiplied by the POSITIVE 4*a2d*c2d.
    k = a2d * c2d;
    C = ratQ.conic([ratQ.chk(k, 'C s1^2'), 0, ratQ.chk(k, 'C s2^2'), ...
                    ratQ.chk(-4*c2d*a2n(1), 'C s1'), ratQ.chk(-4*c2d*a2n(2), 'C s2'), ...
                    ratQ.chk(4*a2d*c2n, 'C 1')]);
    fprintf('tie conic C(s), primitive integers: %s   (%s)\n', mat2str(C), ratQ.conicKind(C));

    ctr = 2*a;   R = 2*sqrt(a.'*a - c);
    th  = linspace(0, 2*pi, 2001).';
    S   = ctr.' + R*[cos(th), sin(th)];
    val = C(1)*S(:,1).^2 + C(2)*S(:,1).*S(:,2) + C(3)*S(:,2).^2 + ...
          C(4)*S(:,1) + C(5)*S(:,2) + C(6);
    fprintf('   max |C| on the derived circle: %.3e   (coefficient scale %g)\n\n', ...
            max(abs(val)), max(abs(C)));

    % ---- the second conic D(s; x), and the branch selector ---------------------------------------
    % x is CONSTRUCTED on the bridge rather than picked, so the test means something: take an s on
    % the tie circle, form its two tangency points, and sit x at their midpoint. A first version
    % used x = (0.4, 1.1) and both roots came back with t outside [0,1] -- the branch test working
    % correctly on a point that is simply not a bridge point, where co f = q1(x) instead.
    s0 = ctr + R*[cos(0.7); sin(0.7)];
    x  = 0.5*s0 + 0.5*(a + s0/2);
    D = @(s) (x(1)-s(1))*(2*a2n(2) - a2d*s(2)) - (x(2)-s(2))*(2*a2n(1) - a2d*s(1));
    fprintf('sample x = (%g, %g)\n', x(1), x(2));

    % Roots by walking the circle for sign changes of D and bisecting -- adequate for a spike whose
    % point is the algebra, not the solver. A first attempt took the sample MINIMISING |D|, which is
    % not a root at all: it returned the opposite point of the circle and z = -21.1 against 0.683.
    gD = arrayfun(@(t) D(ctr + R*[cos(t); sin(t)]), th);
    zBest = -inf;
    for i = 1:numel(th)-1
        if gD(i) == 0 || sign(gD(i)) == sign(gD(i+1)), continue, end
        lo = th(i);  hi = th(i+1);
        for it = 1:80
            mid = (lo+hi)/2;
            if sign(D(ctr + R*[cos(mid); sin(mid)])) == sign(gD(i)), lo = mid; else, hi = mid; end
        end
        sr  = ctr + R*[cos((lo+hi)/2); sin((lo+hi)/2)];
        y1  = sr;   y2 = a + sr/2;   seg = y2 - y1;
        tt  = (seg.' * (x - y1)) / (seg.' * seg);
        zr  = sr.'*x - (sr.'*sr)/2;
        tag = '';
        if tt >= -1e-9 && tt <= 1+1e-9
            tag = '   <- on the bridge';
            if zr > zBest, zBest = zr; end
        end
        fprintf('   root s = (%9.6f, %9.6f)   t = %8.4f   z = %11.6f%s\n', ...
                sr(1), sr(2), tt, zr, tag);
    end

    zCf = 2*(a.'*(x - 2*a)) + R*norm(x - 2*a) + 2*c;
    fprintf('\n   z from the elimination : %.10f\n', zBest);
    fprintf('   z from the closed form : %.10f\n', zCf);
    fprintf('   difference             : %.3e\n\n', abs(zBest - zCf));

    fprintf(['Both conics are exact over Z from exact inputs, and the pair selects the right s once\n' ...
             'the t in [0,1] branch test is applied. The missing piece is the bivariate RESULTANT\n' ...
             'that eliminates s to give p(x,y,z) = 0 -- ratQ has no such routine today.\n']);
end
