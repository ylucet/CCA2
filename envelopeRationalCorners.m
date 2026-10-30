function env = envelopeRationalCorners(env)
% envelopeRationalCorners  Move every corner face whose corner point is RATIONAL from `env.bridge`
%   to `env.quad`, as the rational affine face it is -- so a subdivider compares its VALUE exactly.
%
% [input]  env : as `envelopeAssemble` returns it
% [output] env : the same envelope; such corners now sit in `.quad` (num, den, con, pair) and
%                no longer in `.bridge`
%
% WHY (2026-10-17). Both subdividers carried corners for COVERAGE only, on the strength of a measurement
% that a corner never exceeds the best quadratic face owning the same point (DECISIONS.md
% 2026-09-19). That is false in general: on every one-pair atom envelope a corner is the maximum in
% slivers near a domain vertex, by up to 8e-6, and the subdivisions came out that much too LOW there. A
% corner's value is <s,x> - g(s) at its corner point s, and whenever s is rational -- the facet
% crosses the tie at a rational parameter, including when the discriminant is a perfect square --
% that is an affine face with rational coefficients: a `.quad` face in every respect, whose
% comparisons with the quadratic faces are conics and with the bridges are the same algebraic rows
% a quadratic face gets. Only an irrational corner is left in `.bridge`.
    keep = true(1, numel(env.bridge));
    for k = 1:numel(env.bridge)
        c = env.bridge{k};
        if ~strcmp(c.kind, 'corner') || isempty(c.con) || ~isfield(c, 'sigma'), continue, end
        [sn, sd, ok] = rationalCorner(c.sigma);
        if ~ok, continue, end
        n = c.src.n;  d = c.src.d;
        N = [n(5) n(6); n(6) n(7)];  m = [n(8); n(9)];  qq = n(10);
        A = @abs;
        % <s,x> - g(s), g = (s'Ns/2 + m's + qq)/d, over 2 d sd^2
        g0 = ratQ.chkB(sn.'*N*sn + 2*sd*(m.'*sn) + 2*sd^2*qq, ...
                       A(sn).'*A(N)*A(sn) + 2*sd*(A(m).'*A(sn)) + 2*sd^2*A(qq), 'corner value');
        num = ratQ.chk([0 0 0 0, 0 0 0, 2*d*sd*sn(1), 2*d*sd*sn(2), -g0], 'corner value');
        den = ratQ.chk(2*d*sd^2, 'corner value');
        gg = gcd(den, gcd(gcd(num(8), num(9)), num(10)));
        env.quad(end+1) = struct('num', num / gg, 'den', den / gg, 'con', c.con, ...
                                 'pair', c.pairIndex);
        keep(k) = false;
    end
    env.bridge = env.bridge(keep);
end

function [sn, sd, ok] = rationalCorner(sg)
% s = (s0n + sigma w)/na as integers sn/sd, when sigma is rational.
    sn = [0; 0];  sd = 1;  ok = false;
    al = sg.quad(1);  be = sg.quad(2);  ga = sg.quad(3);
    if sg.root == 0 || al == 0
        if be == 0, return, end
        p = -ga;  q = be;                               % sigma = -ga/be
    else
        D = ratQ.chkB(be^2 - 4*al*ga, be^2 + 4*abs(al)*abs(ga), 'corner discriminant');
        r = round(sqrt(D));
        if D < 0 || r^2 ~= D, return, end               % an irrational corner stays a bridge
        p = -be + sg.root*r;  q = 2*al;
    end
    sn = ratQ.chkB(sg.s0n(:)*q + p*sg.w(:), abs(sg.s0n(:))*abs(q) + abs(p)*abs(sg.w(:)), 'corner point');
    sd = ratQ.chk(sg.na*q, 'corner point');
    if sd < 0, sn = -sn;  sd = -sd; end
    gs = gcd(gcd(sn(1), sn(2)), sd);
    sn = sn / gs;  sd = sd / gs;
    ok = true;
end
